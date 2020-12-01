/*  File: modrep.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2020
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Jul  6 14:01 2020 (rd109)
 * Created: Fri Jun 19 07:38:24 2020 (rd109)
 *-------------------------------------------------------------------
 */

#include <math.h>
#include "modset.h"
#include "seqio.h"

FILE *outFile ;
bool isVerbose = false ;

typedef struct {
  int     len ;
  Modset *ms ;
  int    *pos ;
  bool   *isF ;
} Ref ;

static Ref* refCreate (char *seqFileName, char *modFileName)
{
  int i, n = 0 ;
  FILE *f ;
  Ref *ref = new0 (1, Ref) ;

  if (!(f = fopen (modFileName, "r")))
    die ("failed to open mod file %s", modFileName) ;
  if (!(ref->ms = modsetRead (f)))
    die ("failed to read reference modset from file %s", modFileName) ;
  fclose (f) ;

  ref->pos = new0 (ref->ms->max+1, int) ;
  ref->isF = new0 (ref->ms->max+1, bool) ;
  SeqIO *si = seqIOopenRead (seqFileName, dna2indexConv, false) ; /* false for no qualities */
  if (!si) die ("can't open reference sequence file %s", seqFileName) ;
  if (seqIOread (si))
    { U64 hash, index ; int loc ; bool isF ;
      SeqhashRCiterator *mi = modRCiterator (ref->ms->hasher, sqioSeq(si), si->seqLen) ;
      while (modRCnext (mi, &hash, &loc, (bool*)&isF))
	if ((index = modsetIndexFind (ref->ms, hash, false)))
	  { if (ref->pos[index]) die ("duplicate mod entry at position %d in ref", loc) ;
	    ref->pos[index] = loc ;
	    ref->isF[index] = isF ;
	    if (loc >= ref->len) ref->len = loc+1 ;
	    ++n ;
	  }
    }
  else die ("can't read reference sequence") ;
  if (seqIOread (si)) die ("multiple sequences in ref file - only one allowed") ;

  fprintf (stderr, "found %d of %d locations in ref length %llu\n",
	   n, ref->ms->max, si->seqLen) ;

  seqIOclose (si) ;
  return ref ;
}

typedef struct {
  int k ;
  int x ;
  int n ;
} Hit ;

typedef struct {
  int   i ;	// original read index
  int   len ;	// length
  bool  isF ;	// if it was in same direction as reference - flipped if not
  int   loc0 ;	// location of mod index 1 in this read
  Array hits ;	// of Hit
} Read ;

int readOrder (const void *a, const void *b)
{ Read *ra = (Read*)a, *rb = (Read*)b ;
  int i ;
  int n = arrayMax (ra->hits) ; if (arrayMax (rb->hits) < n) n = arrayMax (rb->hits) ;
  Hit *ha = arrp(ra->hits,0,Hit), *hb = arrp(rb->hits,0,Hit) ;
  for (i = 0 ; i < n ; ++i, ++ha, ++hb)
    if (ha->k > hb->k) return 1 ;
    else if (ha->k < hb->k) return -1 ;
  if (arrayMax (ra->hits) > n) return 1 ;
  else if (arrayMax (rb->hits) > n) return -1 ;
  return 0 ;
}

typedef struct {
  int    n, nPost, nPre ;
  Array  pre, post ;
} Mod ;

static void cleanMods (Mod *mods, Array reads, Modset *ms)
{
  int i, j, k ;
  Read *r ;
  int thresh =  arrayMax(reads) / 2 ;
  
  // remove mods present under 5 times or more than half the possible times
  int nMod0, nMod1, nMod2, nMod3 ;
  nMod0 = 0 ; nMod1 = 0 ; nMod2 = 0 ; nMod3 = 0 ;
  for (i = 0 ; i < ms->max ; ++i)
    if (!mods[i].n) nMod0++ ;
    else if (mods[i].n < 5)  { mods[i].n = 0 ; nMod1++ ; }
    else if (mods[i].n > thresh) { mods[i].n = 0 ; nMod2++ ; }
    else
      { if (!mods[i].pre)
	  { mods[i].pre = arrayCreate (8, Hit) ;
	    mods[i].post = arrayCreate (8, Hit) ;
	  }
	nMod3++ ;
      }
  printf ("NMOD mod0 %d modSmall %d modBig %d modGood %d\n", nMod0, nMod1, nMod2, nMod3) ;
  
  // pack r->mods to remove 0s
  for (i = 0, r = arrp(reads, 0, Read) ; i < arrayMax(reads) ; ++i, ++r) 
    { Hit *h = arrp(r->hits,0,Hit) ;
      for (j = 0, k = 0 ; j < arrayMax(r->hits) ; ++j, ++h)
	if (mods[h->k].n)
	  arr(r->hits, k++, Hit) = arr(r->hits, j, Hit) ;
      arrayMax (r->hits) = k ;
    }
}

static inline void addHit (Array a, int k, int dx) // keeps the most common first in the list
{
  int i ;
  Hit *h0 = arrp(a,0,Hit), *h = h0 ;
  for (i = 0 ; i < arrayMax(a) ; ++i, ++h)
    if (h->k == k)
      { ++h->n ;
	h->x += dx ;
	if (h != h0 && h->n > h0->n) // move to front of the list
	  { Hit hTemp = *h ;
	    while (h > h0) { *h = h[-1] ; --h ; }
	    *h0 = hTemp ;
	  }
	return ;
      }
  h = arrayp(a,i,Hit) ; // may reallocate and so be in a new location
  h->k = k ;
  h->n = 1 ;
  h->x = dx ;
}

static void buildPrePost (Mod *mods, Array reads, Modset *ms)
{
  int i, j ;
  Read *r ;

  for (i = 0 ; i < ms->max ; ++i) // first clean arrays
    if (mods[i].pre)
      { arrayMax (mods[i].pre) = 0 ; mods[i].nPre = 0 ;
	arrayMax (mods[i].post) = 0 ; mods[i].nPost = 0 ;
      }
  for (i = 0, r = arrp(reads, 0, Read) ; i < arrayMax(reads) ; ++i, ++r) 
    { for (j = 1 ; j < arrayMax(r->hits) ; ++j)
	{ int k0 = arrp(r->hits,j-1,Hit)->k, k1 = arrp(r->hits,j,Hit)->k ;
	  int dx = arrp(r->hits,j,Hit)->x - arrp(r->hits,j-1,Hit)->x ;
	  addHit (mods[k0].post, k1, dx) ; ++mods[k0].nPost ;
	  addHit (mods[k1].pre, k0, dx) ; ++mods[k1].nPre ;
	}
    }
}

void analyzeSequences3 (char *seqFileName, char *modFileName, Ref *ref)
{
  int     i, j, k ;
  Read   *r ;
  Mod    *m ;
  Hit    *h ;
  Modset *ms ;
  Array   reads = arrayCreate (12000, Read) ;

  { FILE *f = fopen (modFileName, "r") ;
    if (!f) die ("failed to open mod file %s", modFileName) ;
    if (!(ms = modsetRead (f)))
      die ("failed to read modset from file %s", modFileName) ;
    fclose (f) ;
  }

  Mod *mods = new0 (ms->max, Mod) ;
  
  SeqIO *si = seqIOopenRead (seqFileName, dna2indexConv, false) ; /* false for no qualities */
  if (!si) die ("can't open sequence file %s", seqFileName) ;
  int nRead = 0, nBad = 0 ;
  bool *isDup = new (ms->max, bool) ;
  while (seqIOread (si))
    { ++nRead ;
      U64 hash, index ; int loc ; bool isF ;
      SeqhashRCiterator *mi = modRCiterator (ref->ms->hasher, sqioSeq(si), si->seqLen) ;
      int seqF = 0, seqR = 0, n = 0 ;
      while (modRCnext (mi, &hash, 0, (bool*)&isF) && n < 100)
	if ((index = modsetIndexFind (ref->ms, hash, false))) // in reference
	  { if ((isF && ref->isF[index]) || (!isF && !ref->isF[index])) ++seqF ;
	    else ++seqR ;
	    ++n ;
	  }
      seqhashRCiteratorDestroy (mi) ;
      if (n < 100 || (seqF > 10 && seqR > 10))
	{ ++nBad ;
	  printf ("BADREAD %5d len %5d n %d F %4d R %4d\n",
		  nRead, (int)si->seqLen, n, seqF, seqR) ;
	  continue ;
	}

      r = arrayp (reads, arrayMax(reads), Read) ;
      r->i = nRead-1 ;
      r->loc0 = 0 ;
      r->len = si->seqLen ;
      if (seqF < seqR)		// reverse complement it, in place
	{ char *s = sqioSeq(si) ;
	  for (i = 0, j = si->seqLen-1 ; i < j ; ++i, --j)
	    { char t = 3-s[i] ; s[i] = 3-s[j] ; s[j] = t ; }
	  if (i == j) s[i] = 3-s[i] ;
 	}
      else r->isF = true ;

      mi = modRCiterator (ref->ms->hasher, sqioSeq(si), si->seqLen) ;
      r->hits = arrayCreate (500, Hit) ;
      bzero (isDup, ms->max) ; // reset array
      while (modRCnext (mi, &hash, &loc, 0))
	if ((index = modsetIndexFind (ms, hash, false)))
	  { ++mods[index].n ;
	    if (isDup[index]) ++mods[index].nPre ; else isDup[index] = true ;
	    h = arrayp(r->hits, arrayMax(r->hits), Hit) ;
	    h->k = index ; h->x = loc ;
	  }
      seqhashRCiteratorDestroy (mi) ;
    }
  free (isDup) ;
  fprintf (stderr, "read %d reads, %d bad, %d good: ", nRead, nBad, arrayMax(reads)) ;
  int nMod = 0, nDup = 0, tDup = 0 ;
  for (i = 0 ; i < ms->max ; ++i)
    if (mods[i].nPre)
      { ++nDup ;
	tDup += mods[i].nPre ;
	mods[i].n = 0 ;
      }
    else ++nMod ;
  fprintf (stderr, "mods total %d good %d dup %d avdup %.1f\n",
	   (int)ms->max, nMod, nDup, nDup?(tDup/(double)nDup):0.) ;
  timeUpdate (stderr) ;

  // next make a covering of maximally deep mods
  int minMax = 0 ;
  for (i = 0 ; i < arrayMax(reads) ; ++i)
    { r = arrp(reads,i,Read) ;
      int max = 0 ;
      for (j = 0 ; j < arrayMax(r->hits) ; ++j)
	if (mods[arrp(r->hits,j,Hit)->k].n > max) max = mods[arrp(r->hits,j,Hit)->k].n ;
      if (!minMax || (max < minMax)) minMax = max ;
    }
  fprintf (stderr, "minimum max for a read is %d\n", minMax) ;
  

  for (i = 0 ; i < arrayMax(reads) ; ++i) arrayDestroy (arr(reads, i, Read).hits) ;
  arrayDestroy (reads) ;
  for (i = 0 ; i < ms->max ; ++i)
    if (mods[i].pre) { arrayDestroy (mods[i].pre) ; arrayDestroy (mods[i].post) ; }
  free (mods) ;
  seqIOclose (si) ;
  modsetDestroy (ms) ;
}

/************ OLD VERSIONS ************/

void analyzeSequences1 (char *seqFileName, char *modFileName, Ref *ref)
{
  int     i, j, k ;
  Read   *r ;
  Mod    *m ;
  Hit    *h ;
  Modset *ms ;
  Array   reads = arrayCreate (12000, Read) ;

  { FILE *f = fopen (modFileName, "r") ;
    if (!f) die ("failed to open mod file %s", modFileName) ;
    if (!(ms = modsetRead (f)))
      die ("failed to read modset from file %s", modFileName) ;
    fclose (f) ;
  }

  Mod *mods = new0 (ms->max, Mod) ;
  
  SeqIO *si = seqIOopenRead (seqFileName, dna2indexConv, false) ; /* false for no qualities */
  if (!si) die ("can't open sequence file %s", seqFileName) ;
  int nRead = 0, nBad = 0 ;
  while (seqIOread (si))
    { ++nRead ;
      U64 hash, index ; int loc ; bool isF ;
      SeqhashRCiterator *mi = modRCiterator (ref->ms->hasher, sqioSeq(si), si->seqLen) ;
      int seqF = 0, seqR = 0, n = 0 ;
      while (modRCnext (mi, &hash, 0, (bool*)&isF) && n < 100)
	if ((index = modsetIndexFind (ref->ms, hash, false))) // in reference
	  { if ((isF && ref->isF[index]) || (!isF && !ref->isF[index])) ++seqF ;
	    else ++seqR ;
	    ++n ;
	  }
      seqhashRCiteratorDestroy (mi) ;
      if (n < 100 || (seqF > 10 && seqR > 10))
	{ ++nBad ;
//	  printf ("BAD %5d len %5d start %d end %d F %4d R %4d\n",
//		 nRead, (int)si->seqLen, start, end, seqF, seqR) ;
	  continue ;
	}
      if (seqF < seqR)		// reverse complement it, in place
	{ char *s = sqioSeq(si) ;
	  for (i = 0, j = si->seqLen-1 ; i < j ; ++i, --j)
	    { char t = 3-s[i] ; s[i] = 3-s[j] ; s[j] = t ; }
	  if (i == j) s[i] = 3-s[i] ;
	}

      mi = modRCiterator (ref->ms->hasher, sqioSeq(si), si->seqLen) ;
#ifdef SEGMENT      
      while (modRCnext (mi, &hash, &loc, 0) && (loc + 5000 < si->seqLen))
	if ((index = modsetIndexFind (ref->ms, hash, false)) && index == 1)
	  { int locEnd = loc + 5000 ;
	    r = arrayp (reads, arrayMax(reads), Read) ;
	    r->i = nRead-1 ;
	    r->loc0 = loc ;
	    r->len = si->seqLen ;
	    r->hits = arrayCreate (500, Hit) ;
	    while (modRCnext (mi, &hash, &loc, 0) && (loc < locEnd))
	      if ((index = modsetIndexFind (ms, hash, false)))
		{ ++mods[index].n ;
		  h = arrayp(r->hits, arrayMax(r->hits), Hit) ;
		  h->k = index ; h->x = loc - r->loc0 ;
		}
	  }
#else 
      r = arrayp (reads, arrayMax(reads), Read) ;
      r->i = nRead-1 ;
      r->loc0 = 0 ;
      r->len = si->seqLen ;
      r->hits = arrayCreate (500, Hit) ;
      while (modRCnext (mi, &hash, &loc, 0))
	if ((index = modsetIndexFind (ms, hash, false)))
	  { ++mods[index].n ;
	    h = arrayp(r->hits, arrayMax(r->hits), Hit) ;
	    h->k = index ; h->x = loc ;
	  }
#endif      
      seqhashRCiteratorDestroy (mi) ;
    }

  fprintf (stderr, "read %d reads, %d bad, %d good: ", nRead, nBad, arrayMax(reads)) ;
  timeUpdate (stderr) ;

  cleanMods (mods, reads, ms) ; // remove mods with count < 5 or > half the readcount
  
  // pack r->mods to remove runs
  int K = ms->hasher->k ;
  for (i = 0, r = arrp(reads, 0, Read) ; i < arrayMax(reads) ; ++i, ++r) 
    { int xNext = 0 ;
      h = arrp(r->hits,0,Hit) ;
      for (j = 0, k = 0 ; j < arrayMax(r->hits) ; ++j, ++h)
	if (h->x >= xNext)
	  { arr(r->hits, k++, Hit) = arr(r->hits, j, Hit) ;
	    xNext = h->x + K ;
	  }
	else
	  --mods[h->k].n ;
      arrayMax (r->hits) = k ;
    }
  cleanMods (mods, reads, ms) ;
  
  // next remove redundant mods, and bad mods
  buildPrePost (mods, reads, ms) ;
  for (i = 0 ; i < ms->max ; ++i)
    if (mods[i].n)
      { h = arrp(mods[i].pre, 0 , Hit) ;
	if (h->n == mods[i].n && h->n == mods[h->k].nPost) // no new info in this mod
	  mods[i].n = 0 ;
	else
	  { bool isBad = true ;
	    int nThresh = mods[i].n / 2 ;
	    for (j = 0 ; j < arrayMax(mods[i].pre) ; ++j, ++h)
	      if (isBad && h->n >= 5 && (h->n > nThresh || h->n > mods[h->k].nPost / 2))
		isBad = false ;
	    h = arrp(mods[i].post, 0 , Hit) ;
	    for (j = 0 ; j < arrayMax(mods[i].post) ; ++j, ++h)
	      if (isBad && h->n >= 5 && (h->n > nThresh || h->n > mods[h->k].nPre / 2))
		isBad = false ;
	    if (isBad) mods[i].n = 0 ;
	  }
      }
  cleanMods (mods, reads, ms) ;
  
  // now remove reads that contain links with support < 5 - a bit radical
  buildPrePost (mods, reads, ms) ;
  int kp ;
  r = arrp(reads, 0, Read) ;
  for (i = 0, k = 0 ; i < arrayMax (reads) ; ++i, ++r)
    { h = arrp(r->hits, 0, Hit) ;
      for (j = 1 ; j < arrayMax(r->hits) ; ++j)
	{ Hit *hp = arrp(mods[h->k].post, 0, Hit) ;
	  int n = arrayMax(mods[h->k].post) ;
	  ++h ; // increment here in the loop, so that now it points to r->hits[j]
	  for (kp = 0 ; kp < n ; ++kp, ++hp)
	    if (hp->k == h->k) break ; // found it!
	  assert (kp < n) ; // else some screwup
	  if (hp->n < 5) break ; // a weak link
	}
      if (j < arrayMax(r->hits)) // read contained a weak link
	arrayDestroy (r->hits) ;
      else
	{ arr(reads,k,Read) = *r ; ++k ; }
    }
  fprintf (stderr, "reduced %d reads to %d reads\n", arrayMax(reads), k) ;
  arrayMax(reads) = k ;
  // rebuild Mods->n
  for (i = 0 ; i < ms->max ; ++i) mods[i].n = 0 ;
  r = arrp(reads, 0, Read) ;
  for (i = 0 ; i < arrayMax (reads) ; ++i, ++r)
    { h = arrp(r->hits, 0, Hit) ;
      for (j = 1 ; j < arrayMax(r->hits) ; ++j, ++h)
	++mods[h->k].n ;
    }
  cleanMods (mods, reads, ms) ;
    
  // next remove redundant mods, and bad mods
  buildPrePost (mods, reads, ms) ;
  for (i = 0 ; i < ms->max ; ++i)
    if (mods[i].n)
      { h = arrp(mods[i].pre, 0 , Hit) ;
	if (h->n == mods[i].n && h->n == mods[h->k].nPost) // no new info in this mod
	  mods[i].n = 0 ;
	else
	  { bool isBad = true ;
	    int nThresh = mods[i].n / 2 ;
	    for (j = 0 ; j < arrayMax(mods[i].pre) ; ++j, ++h)
	      if (isBad && h->n >= 5 && (h->n > nThresh || h->n > mods[h->k].nPost / 2))
		isBad = false ;
	    h = arrp(mods[i].post, 0 , Hit) ;
	    for (j = 0 ; j < arrayMax(mods[i].post) ; ++j, ++h)
	      if (isBad && h->n >= 5 && (h->n > nThresh || h->n > mods[h->k].nPre / 2))
		isBad = false ;
	    if (isBad) mods[i].n = 0 ;
	  }
      }
  cleanMods (mods, reads, ms) ;

  // report
  buildPrePost (mods, reads, ms) ;
  for (i = 0 ; i < ms->max ; ++i)
    if (mods[i].n)
      { printf ("MOD %d n %d pre %d (", i, mods[i].n, mods[i].nPre) ;
	h = arrp(mods[i].pre, 0 , Hit) ;
	for (j = 0 ; j < arrayMax(mods[i].pre) ; ++j, ++h)
	  printf (" %d:%d|%d:%d", h->k, h->n, mods[h->k].nPost, h->x/h->n) ;
	printf (") post %d (", mods[i].nPost) ;
	h = arrp(mods[i].post, 0 , Hit) ;
	for (j = 0 ; j < arrayMax(mods[i].post) ; ++j, ++h)
	  printf (" %d:%d|%d:%d", h->k, h->n, mods[h->k].nPre, h->x/h->n) ;
	printf (")\n") ;
      }

  arraySort (reads, readOrder) ;

  int block = 0 ;
  for (i = 0, r = arrp(reads, 0, Read) ; i < arrayMax(reads) ; ++i, ++r)
    { if (i && readOrder (r, r-1))
	{ printf ("BLOCK %3d", block) ; block = 0 ;
	  --r ;
	  for (j = 0 ; j < arrayMax(r->hits) ; ++j)
	    printf ("\t%5d", arrp(r->hits,j,Hit)->k) ;
	  ++r ;
	  putchar ('\n') ;
	}
      ++block ;
      printf ("READ %5d n %3d mods", r->i, arrayMax(r->hits)) ;
      for (j = 0 ; j < arrayMax(r->hits) ; ++j)
	printf ("\t%5d", arrp(r->hits,j,Hit)->k) ;
      putchar ('\n') ;
    }

  for (i = 0 ; i < arrayMax (reads) ; ++i) arrayDestroy (arr(reads, i, Read).hits) ;
  arrayDestroy (reads) ;
  for (i = 0 ; i < ms->max ; ++i)
    if (mods[i].pre) { arrayDestroy (mods[i].pre) ; arrayDestroy (mods[i].post) ; }
  free (mods) ;
  seqIOclose (si) ;
  modsetDestroy (ms) ;
}

/****************************/

int boundary[] = { 1, 		// pos     8 depth 4416
		   961,		// pos  6166 depth 3555
		   1951,	// pos 13919 depth 3289
		   2961 } ;	// pos 25031 depth 4535

void analyzeSequences2 (char *seqFileName, char *modFileName, Ref *ref)
{
  int i, j ;
  Modset *ms ;
  Array   reads = arrayCreate (12000, Read) ;

  { FILE *f = fopen (modFileName, "r") ;
    if (!f) die ("failed to open mod file %s", modFileName) ;
    if (!(ms = modsetRead (f)))
      die ("failed to read modset from file %s", modFileName) ;
    fclose (f) ;
  }
  Mod *mods = new0 (ms->max, Mod) ;
  U64 *hits = new (ms->max, U64) ;

  //  for (i = 0 ; i < ref->ms->max ; ++i)
  //    if ((j = modsetIndexFind (ms, ref->ms->value[i], false)))
  //      printf ("REF %4d pos %5d depth %5d\n", i, ref->pos[i], ms->depth[j]) ;

  SeqIO *si = seqIOopenRead (seqFileName, dna2indexConv, false) ; /* false for no qualities */
  if (!si) die ("can't open sequence file %s", seqFileName) ;
  int n1 = 0, n2 = 0, n3 = 0, n4 = 0 ;
  while (seqIOread (si))
    { SeqhashRCiterator *mi = modRCiterator (ref->ms->hasher, sqioSeq(si), si->seqLen) ;
      U64 hash, index ; int loc ; bool isF ;
      bool is0 = false, is1 = false, is2 = false, is3 = false ;
      while (modRCnext (mi, &hash, &loc, (bool*)&isF))
	if ((index = modsetIndexFind (ref->ms, hash, false))) // in reference
	  { if (index == boundary[0]) is0 = true ;
	    else if (index == boundary[1]) is1 = true ;
	    else if (index == boundary[2]) is2 = true ;
	    else if (index == boundary[3]) is3 = true ;
	  }
      if (is0 && is1) ++n1 ;
      if (is1 && is2) ++n2 ;
      if (is2 && is3) ++n3 ;
      if (is3 && is0) ++n4 ;
    }
  printf ("n1 %d n2 %d n3 %d n4 %d\n", n1, n2, n3, n4) ;
  
  modsetDestroy (ms) ;
}

void usage (void)
{ fprintf (stderr, "Usage: modrep <commands>\n") ;
  fprintf (stderr, "Commands are executed in order - set parameters before using them!\n") ;
  fprintf (stderr, "  -v | --verbose : toggle verbose mode\n") ;
  fprintf (stderr, "  -o | --output <output_filename> : '-' for stdout\n") ;
  fprintf (stderr, "  -R | --ref <seq_file> <mod_file>\n") ;
  fprintf (stderr, "  -s1 | --seq1 <seq_file> <mod_file>: analyse reads\n") ;
  fprintf (stderr, "  -s2 | --seq2 <seq_file> <mod_file>: analyse reads\n") ;
  fprintf (stderr, "  -s3 | --seq3 <seq_file> <mod_file>: analyse reads\n") ;
}

int main (int argc, char *argv[])
{
  --argc ; ++argv ;		/* eat program name */
  if (!argc) usage () ;

  outFile = stdout ;
  timeUpdate (stdout) ;		/* initialise timer */
  dna2indexConv['N'] = dna2indexConv['n'] = 0 ; /* to get 2-bit encoding */

  int i ;			/* generically useful variables */
  FILE *f ;
  Ref *ref = 0 ;

  while (argc)
    { if (**argv != '-')
	die ("option/command %s does not start with '-': run without arguments for usage", *argv) ;
      fprintf (stderr, "COMMAND %s", *argv) ;
      for (i = 1 ; i < argc && *argv[i] != '-' ; ++i) fprintf (stderr, " %s", argv[i]) ;
      fputc ('\n', stderr) ;
      
#define ARGMATCH(x,y,n)	((!strcmp (*argv, x) || (!strcmp (*argv,y))) && argc >= n && (argc -= n, argv += n))
      if (ARGMATCH("-v","--verbose",1)) isVerbose = !isVerbose ;
      else if (ARGMATCH("-o","--output",2))
	{ if (!strcmp (argv[-1], "-"))
	    outFile = stdout ;
	  else if (!(outFile = fopen (argv[-1], "w")))
	    { fprintf (stderr, "can't open output file %s - resetting to stdout\n", argv[-1]) ;
	      outFile = stdout ;
	    }
	}
      else if (ARGMATCH("-R","--ref",3))
	ref = refCreate (argv[-2], argv[-1]) ;
      else if (ARGMATCH("-s1","--seq1",3))
	if (ref)
	  analyzeSequences1 (argv[-2], argv[-1], ref) ;
        else
	  die ("you must read reference data with -R before command -s") ;
      else if (ARGMATCH("-s2","--seq2",3))
	if (ref)
	  analyzeSequences2 (argv[-2], argv[-1], ref) ;
        else
	  die ("you must read reference data with -R before command -s") ;
      else if (ARGMATCH("-s3","--seq3",3))
	if (ref)
	  analyzeSequences3 (argv[-2], argv[-1], ref) ;
        else
	  die ("you must read reference data with -R before command -s") ;
      else
	die ("unknown option %s", *argv) ;
    }

  fprintf (stderr, "total resources used: ") ; timeTotal (stderr) ;
}

/************* end of file ************/
