/*  File: modset.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2018
 *-------------------------------------------------------------------
 * Description: package to handle sets of "mod" sequence hashes
 * Exported functions:
 * HISTORY:
 * Last edited: Jan 26 22:33 2019 (rd109)
 * Created: Tue Nov  6 17:31:14 2018 (rd109)
 *-------------------------------------------------------------------
 */

#include "modset.h"
  
Modset *modsetCreate (Seqhash *sh, int bits, U32 size)
{
  if (bits < 20 || bits > 34) die ("table bits %d must be between 20 and 34", bits) ;
  Modset *ms = new0 (1, Modset) ;
  ms->hasher = sh ;
  ms->tableBits = bits ;
  ms->tableSize = (U64)1 << ms->tableBits ;
  ms->tableMask = ms->tableSize - 1 ;
  ms->index = new0 (ms->tableSize, U32) ;
  if (size >= (ms->tableSize >> 2)) die ("Modset size %u is too big for %d bits", size, bits) ;
  else if (size) ms->size = size ;
  else ms->size = (ms->tableSize >> 2) - 1 ;
  ms->value = new (ms->size, U64) ;
  ms->depth = new0 (ms->size, U16) ;
  ms->info = new0 (ms->size, U8) ;
  return ms ;
}

void modsetDestroy (Modset *ms)
{ free (ms->index) ; free (ms->value) ; free (ms->info) ; free (ms) ; }

BOOL modsetPack (Modset *ms)	/* compress per-item arrays */
{ if (ms->size == ms->max+1) return FALSE ;
  resize (ms->value, ms->size, ms->max+1, U64) ;
  resize (ms->depth, ms->size, ms->max+1, U16) ;
  resize (ms->info, ms->size, ms->max+1, U8) ;
  ms->size = ms->max+1 ;
  return TRUE ;
}

U32 modsetIndexFind (Modset *ms, U64 hash, int isAdd)
{
  U64 diff = 0 ;
  U64 offset = hash & ms->tableMask ;
  U32 index = ms->index[offset] ;
  while (index && (ms->value[index] != hash))
    { if (!diff) diff = ((hash >> ms->tableBits) & ms->tableMask) | 1 ; /* odd so comprime */
      offset = (offset + diff) & ms->tableMask ;
      index = ms->index[offset] ;
    }
  if (!index && isAdd)
    { index = ms->index[offset] = ++ms->max ;
      if (ms->max >= ms->size) die ("hashTableSize %u is too small for %u", ms->size, ms->max) ;
      ms->value[index] = hash ;
    }
  return index ;
}

void modsetDepthPrune (Modset *ms, int min, int max)
{
  U32 i ;
  U32 N = ms->max ; ms->max = 0 ;
  memset (ms->index, 0, ms->tableSize*sizeof(U32)) ;
  for (i = 1 ; i <= N ; ++i)	/* NB index runs from 1..max */
    if (ms->depth[i] >= min && (!max || ms->depth[i] < max))
      { modsetIndexFind (ms, ms->value[i], TRUE) ;
	ms->info[ms->max] = ms->info[i] ;
	ms->depth[ms->max] = ms->depth[i] ;
      }
  fprintf (stderr, "  pruned Modset from %d to %d with min %d <= depth < max %d\n",
	   N, ms->max, min, max) ;
}

void modsetWrite (Modset *ms, FILE *f)
{ if (fwrite ("MSHSTv1",8,1,f) != 1) die ("failed to write modset header") ;
  if (fwrite (&ms->tableBits,sizeof(int),1,f) != 1) die ("failed to write bits") ;
  U32 size = ms->max+1 ; if (fwrite (&size,sizeof(U32),1,f) != 1) die ("failed to write size") ;
  seqhashWrite (ms->hasher, f) ;
  if (fwrite (ms->index,sizeof(U32),ms->tableSize,f) != ms->tableSize) die ("fail write index") ;
  if (fwrite (ms->value,sizeof(U64),ms->max+1,f) != ms->max+1) die ("failed to write value") ;
  if (fwrite (ms->depth,sizeof(U16),ms->max+1,f) != ms->max+1) die ("failed to write depth") ;
  if (fwrite (ms->info,sizeof(U8),ms->max+1,f) != ms->max+1) die ("failed to write info") ;
}

Modset *modsetRead (FILE *f)
{ char name[8] ;
  if (fread (name,8,1,f) != 1) die ("failed to read modset header") ;
  if (strcmp (name, "MSHSTv1")) die ("bad reference header") ;
  int bits ; if (fread (&bits,sizeof(int),1,f) != 1) die ("failed to read bits") ;
  U32 size ; if (fread (&size,sizeof(U32),1,f) != 1) die ("failed to read size") ;
  Seqhash *sh = seqhashRead (f) ;
  Modset *ms = modsetCreate (sh, bits, size) ;
  if (fread (ms->index,sizeof(U32),ms->tableSize,f) != ms->tableSize) die ("failed read index") ;
  if (fread (ms->value,sizeof(U64),size,f) != size) die ("failed to read value") ;
  if (fread (ms->depth,sizeof(U16),size,f) != size) die ("failed to read depth") ;
  if (fread (ms->info,sizeof(U8),size,f) != size) die ("failed to read info") ;
  ms->max = size - 1 ;
  return ms ;
}

BOOL modsetMerge (Modset *ms1, Modset *ms2)
{
  U32 i, index ;
  /* first check that the hashers are identical */
  Seqhash *sh1 = ms1->hasher, *sh2 = ms2->hasher ;
  if (sh1->w != sh2->w || sh1->k != sh2->k || sh1->factor1 != sh2->factor1) return FALSE ;
  /* need to expand size of ms1 to make space */
  U64 newSize = ms1->max + ms2->max + 1 ;
  if (newSize >= (ms1->tableSize >> 2)) newSize = (ms1->tableSize >> 2) - 1 ;
  resize (ms1->value, ms1->size, newSize, U64) ; 
  resize (ms1->depth, ms1->size, newSize, U16) ; 
  resize (ms1->info, ms1->size, newSize, U8) ;
  ms1->size = newSize ;
  /* then pass through ms2 adding into ms1 */
  for (i = 1 ; i <= ms2->max ; ++i)
    { index = modsetIndexFind (ms1, ms2->value[i], TRUE) ;
      U32 d = ms1->depth[index] + (U32) ms2->depth[i] ; if (d > U16MAX) d = U16MAX ;
      ms1->depth[index] = d ;
      int c = msCopy(ms1,index) + msCopy(ms2,i) ; if (c > 3) c = 3 ;
      ms1->info[index] &= 0x3 ; ms1->info[index] |= c ;
    }
  return TRUE ;
}

void modsetSummary (Modset *ms, FILE *f)
{
  seqhashReport (ms->hasher, f) ;
  fprintf (f, "MS table bits %d size %llu number of entries %u",
	   ms->tableBits, ms->tableSize, ms->max) ;
  if (!ms->max) { fputc ('\n', f) ; return ; }
  U32 i, copy[4] ; copy[0] = copy[1] = copy[2] = copy[3] = 0 ;
  Array h = arrayCreate (256, U32) ;
  for (i = 1 ; i <= ms->max ; ++i)
    { if (ms->depth[i] < arrayMax(h)) ++arr(h,ms->depth[i],U32) ; /* more efficient to check */
      else ++array(h,ms->depth[i],U32) ;
      ++copy[msCopy(ms,i)] ;
    }
  U64 sum = 0, tot = 0 ;
  for (i = 0 ; i < arrayMax(h) ; ++i) { sum += arr(h,i,U32) ; tot += i * arr(h,i,U32) ; }
  I64 htot = tot / 2 ;
  for (i = 0 ; i < arrayMax(h) ; ++i) { htot -= i*arr(h,i,U32) ; if (htot < 0) break ; }
  fprintf (f, " total count %llu\nMS average depth %.1f N50 depth %u",
	   tot, tot / (double)sum , i) ;
  if (copy[0] < ms->max)
    fprintf (f, " copy0 %u copy1 %u copy2 %u copyM %u", copy[0], copy[1], copy[2], copy[3]) ;
  fputc ('\n', f) ;
  arrayDestroy (h) ;
}

/***************************************************/
