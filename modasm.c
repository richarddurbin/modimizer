/*  File: modasm.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2018
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Sep 20 15:17 2020 (rd109)
 * * Jan  4 00:14 2019 (rd109): completed findOverlaps and markBadReads and cleaned Read object
 * Created: Tue Nov  6 18:30:49 2018 (rd109)
 *-------------------------------------------------------------------
 */

#include "modset.h"
#include "seqio.h"
#include <math.h>

#ifdef OMP
#include <omp.h>
#endif

#define TOPBIT  0x80000000	/* set for FORWARD orientation, unset for reverse orientation */
#define TOPMASK 0x7fffffff

int numThreads = 1 ;		/* default to serial - reset if multi-threaded */
FILE *outFile ;			/* initialise to stdout at start of main() */
bool isVerbose = false ;

/* drop the read names, and don't store the sequences: mod indexes and spacings only */

typedef struct {
  int len ;			/* sequence length in base pairs */
  int nHit ;			/* number of ms mod hit */
  U32 *hit ;			/* list of hits (ms indexes, using top bit for orientation) */
  U16 *dx ;	     		/* lists of distances to previous hit (increasing order) */
  union {
    U8 bad ;			/* allows single test for any type of bad read */
    struct {
      U8 badRepeat : 1 ;	/* contains repeated single copy mod(es) */
      U8 badOrder10 : 1 ;	/* matches 10 or more other reads in incorrect order */
      U8 badOrder1 : 1 ;	/* matches 1 or more other good reads in incorrect order */
      U8 badNoMatch : 1 ;	/* doesn't match any other reads on UNIQUE  */
      U8 badLowHit : 1 ;	/* hit:miss ratio too low - below 5% */
      U8 badLowCopy1 : 1 ;	/* low number of copy 1 hits */
    } ;
  } ;
  union {
    U8 otherFlags ;
    struct {
      U8 isRDNA : 1 ;
      /* space here for more... */
    } ;
  } ;
  U16 pad1 ;			/* for future use - to pad out flags to 32 bits */
  int nMiss ;			/* number of mods in sequence that miss ms */
  int contained ;		/* index in rs->reads of read that contains this one */
  int nCopy[4] ;		/* number of copy 0, 1, 2, M hits - could create on demand */
  U32 pad2[4] ;			/* for future use */
} Read ;

typedef struct {
  union  {
    U8 isRDNA ;            // all the below indicate in the rDNA repeat units
    struct  {
      U8 isRefRDNA : 1 ;   // in the rDNA reference
      U8 isCoreRDNA : 1 ;  // defined as 2.75-4.75k depth
      U8 isVarRDNA : 1 ;   // under 2.75k and real
      U8 isMultiRDNA : 1 ; // over 4.75k and real
    } ;
  } ;
  int rDNApos ; // (inferred approximate) position in canonical rDNA reference: -1 if inconsistent
  int nGood ;   // number of good partners
  int nMod2 ;   // number of poor partners
  int nBadLD ;  // number of times it appears as a bad partner in an LD test of another mod
  int nSplit ;  // number of mods found both sides of this mod
  int nSplitLD ; // number of times it appears as a split partner in an LD test of another mod
} ModInfo ;

typedef struct {		/* data structure for long read set */
  Modset *ms ;
  ModInfo *modInfo ;
  Array reads ;			/* of Read */
  U64   totHit ;	        /* total number of hits */
  U32   **inv ;			/* inverse array: for each mod, a pointer to a list of reads */
  U32   *invSpace ;		/* inv[i] points into this */
} Readset ;

static void invBuild (Readset *rs) ; /* forward declaration */

Readset *readsetCreate (Modset *ms, int n) /* create empty - n is used to initialise arrays */
{
  Readset *rs = new0 (1, Readset) ;
  rs->ms = ms ;
  rs->reads = arrayCreate (n, Read) ;
  arrayp(rs->reads,0,Read)->nHit = 0 ; /* burn first read, so can use 0 as null */
  /* create inv, invSpace in invBuild() */
  return rs ;
}

void readsetDestroy (Readset *rs)
{
  int i, n = arrayMax(rs->reads) ; Read *r = arrp(rs->reads,0,Read) ;
  for (i = 1 ; i < n ; ++i, ++r) if (r->nHit) { free(r->hit) ; free(r->dx) ; }
  arrayDestroy (rs->reads) ;
  if (rs->inv) { free (rs->inv) ; free (rs->invSpace) ; }
  if (rs->modInfo) free (rs->modInfo) ;
  free (rs) ;
}

void readsetWrite (Readset *rs, char *root)
{
  FILE *f ;
  if (!(f = fopenTag (root, "mod", "w"))) die ("can't open file %s.mod", root) ;
  modsetWrite (rs->ms, f) ; fclose (f) ;
  if (!(f = fopenTag (root, "readset", "w"))) die ("can't open file %s.readset", root) ;
  if (fwrite("RSMSHv2",8,1,f) != 1) die ("failed to write readset header") ;
  if (fwrite (&rs->totHit,sizeof(U64),1,f) != 1) die ("failed to write totHit") ;
  arrayWrite (rs->reads, f) ;
  int i, n ; Read *r = arrp(rs->reads,1,Read) ;
  for (i = 1 ; i < arrayMax(rs->reads) ; ++i, ++r)
    if ((n = r->nHit))
      { if (fwrite(r->hit,sizeof(U32),n,f) != n) die ("failed write hits %d", n) ;
	if (fwrite(r->dx,sizeof(U16),n,f) != n) die ("failed write dx %d", n) ;
      }
  fclose (f) ;
}

Readset *readsetRead (char *root)
{
  FILE *f ;
  if (!(f = fopenTag (root, "mod", "r"))) die ("can't open file %s.mod", root) ;
  Modset *ms = modsetRead (f) ; fclose (f) ;
  if (!(f = fopenTag (root, "readset", "r"))) die ("can't open file %s.readset", root) ;
  char name[8] ;
  if (fread (name,8,1,f) != 1) die ("failed to read readset header") ;
  if (strcmp (name, "RSMSHv2")) die ("bad readset header %s != RSMSHv2", name) ;
  Readset *rs = readsetCreate (ms, 16) ;
  if (fread (&rs->totHit,sizeof(U64),1,f) != 1) die ("failed to read totHit") ;
  arrayDestroy (rs->reads) ; rs->reads = arrayRead (f) ;
  int i, n ; Read *r = arrp(rs->reads,1,Read) ;
  for (i = 1 ; i < arrayMax(rs->reads) ; ++i, ++r)
    if ((n = r->nHit))
      { r->hit = new (n, U32) ; if (fread(r->hit,sizeof(U32),n,f) != n) die ("failed read hits") ;
	r->dx = new (n, U16) ;	if (fread(r->dx,sizeof(U16),n,f) != n) die ("failed read dx") ;
      }
  fclose (f) ;
  invBuild (rs) ;
  return rs ;
}

void readsetFileRead (Readset *rs, char *filename)
{
  char *seq ;			/* ignore the name for now */
  int len ;
  Array hitsA = arrayCreate (1024, U32) ; /* reuse these to build the lists of hits and dx */
  Array dxA = arrayCreate (1024, U16) ;

  memset (rs->ms->depth, 0, (rs->ms->max+1)*sizeof(U16)) ; /* rebuild depth from this file */
  dna2indexConv['N'] = dna2indexConv['n'] = 0 ; /* to get 2-bit encoding */
  SeqIO *si = seqIOopenRead (filename, dna2indexConv, false) ;
  while (seqIOread (si))
    { Read *read = arrayp(rs->reads, arrayMax(rs->reads), Read) ;
      read->len = si->seqLen ;
      SeqhashRCiterator *mi = modRCiterator (rs->ms->hasher, sqioSeq(si), si->seqLen) ;
      hitsA = arrayReCreate (hitsA, 1024, U32) ;
      dxA = arrayReCreate (dxA, 1024, U16) ;
      U64 kmer ; int lastPos = 0, pos ; bool isForward ;
      while (modRCnext (mi, &kmer, &pos, &isForward))
	{ U32 index = modsetIndexFind (rs->ms, kmer, false) ;
	  if (index)
	    { array(hitsA,read->nHit,U32) = isForward ? (index | TOPBIT) : index ;
	      array(dxA,read->nHit,U16) = pos - lastPos ; lastPos = pos ;
	      ++read->nHit ;
	      U16 *di = &rs->ms->depth[index] ; ++*di ; if (!*di) *di = U16MAX ;
	    }
	  else ++read->nMiss ;
	}
      seqhashRCiteratorDestroy (mi) ;
      if (read->nHit)
	{ read->hit = new (read->nHit, U32) ;
	  memcpy (read->hit, arrp(hitsA,0,U32), read->nHit*sizeof(U32)) ;
	  read->dx = new(read->nHit,U16) ;
	  memcpy (read->dx, arrp(dxA,0,U16), read->nHit*sizeof(U16)) ;
	  rs->totHit += read->nHit ;
	}
    }
  seqIOclose (si) ;

  invBuild (rs) ;
  arrayDestroy (hitsA) ; arrayDestroy (dxA) ;
}

void readsetStats (Readset *rs)
{
  U32 n = arrayMax(rs->reads)-1 ;
  if (!n) { fprintf (stderr, "stats called on empty readset\n") ; return ; }

  modsetSummary (rs->ms, outFile) ;
  
  U32 i, j ; 
  int nUnique0 = 0, nUnique1 = 0 ;
  U64 totLen = 0, totMiss = 0, lenUnique0 = 0, lenUnique1 = 0 ;
  U64 totCopy[4] ; totCopy[0] = totCopy[1] = totCopy[2] = totCopy[3] = 0 ;
  U32 nBad = 0, nBadRepeat = 0, nBadOrder10 = 0, nBadOrder1 = 0 ;
  U32 nBadNoMatch = 0, nBadLowHit = 0, nBadLowCopy1 = 0 ;
  
  for (i = 1 ; i <= n ; ++i)
    { Read *read = arrp(rs->reads, i, Read) ;
      totLen += read->len ;
      totMiss += read->nMiss ;
      for (j = 0 ; j < 4 ; ++j) totCopy[j] += read->nCopy[j] ;
      if (read->nCopy[1] == 0) { ++nUnique0 ; lenUnique0 += read->len ; }
      else if (read->nCopy[1] == 1) { ++nUnique1 ; lenUnique1 += read->len ; }
      if (read->bad)
	{ ++nBad ;
	  if (read->badRepeat) ++nBadRepeat ;
	  if (read->badOrder10) ++nBadOrder10 ;
	  if (read->badOrder1) ++nBadOrder1 ;
	  if (read->badNoMatch) ++nBadNoMatch ;
	  if (read->badLowHit) ++nBadLowHit ;
	  if (read->badLowCopy1) ++nBadLowCopy1 ;
	}
    }
  fprintf (outFile, "RS %d sequences, total length %llu (av %.1f)\n",
	   n, totLen, totLen/(double)n) ;
  fprintf (outFile, "RS %llu mod hits, %.1f bp/hit, frac hit %.2f, av hits/read %.1f\n",
	   rs->totHit, totLen/(double)rs->totHit, rs->totHit/(double)(totMiss+rs->totHit),
	   rs->totHit/(double)n) ;
  fprintf (outFile, "RS hit distribution %.2f copy0, %.2f copy1, %.2f copy2, %.2f copyM\n",
	   totCopy[0]/(double)rs->totHit, totCopy[1]/(double)rs->totHit,
	   totCopy[2]/(double)rs->totHit, totCopy[3]/(double)rs->totHit) ;
  U32 nUniqueMulti = n - nUnique0 - nUnique1 ;
  fprintf (outFile, "RS num reads and av_len with 0 copy1 hits %d %.1f with 1 copy1 hits %d %.1f"
	   " >1 copy1 hits %d %.1f av copy1 hits %.1f\n",
	   nUnique0, lenUnique0/(double)nUnique0, nUnique1, lenUnique1/(double)nUnique1,
	   nUniqueMulti, (totLen - lenUnique0 - lenUnique1)/(double)nUniqueMulti,
	   (totCopy[1]-nUnique1)/(double)nUniqueMulti) ;
  fprintf (outFile, "RS bad %u : %u repeat, %u order10, %u order1, ",
	   nBad, nBadRepeat, nBadOrder10, nBadOrder1) ;
  fprintf (outFile, "%u no_match, %u low_hit, %u low_copy1\n",
	   nBadNoMatch, nBadLowHit, nBadLowCopy1) ;
  U32 nCopy[4], hitCopy[4], hit2Copy[4] ;
  U64 depthCopy[4] ;
  for (j = 0 ; j < 4 ; ++j) nCopy[j] = hitCopy[j] = hit2Copy[j] = depthCopy[j] = 0 ;
  for (i = 1 ; i <= rs->ms->max ; ++i)
    { j = msCopy(rs->ms,i) ;
      ++nCopy[j] ;
      if (rs->ms->depth[i] > 0) ++hitCopy[j] ;
      if (rs->ms->depth[i] > 1) { ++hit2Copy[j] ; depthCopy[j] += rs->ms->depth[i] ; }
    }
  fprintf (outFile, "RS mod frac hit hit>1 av: copy0 %.3f %.3f %.1f copy1 %.3f %.3f %.1f copy2 %.3f %.3f %.1f copyM %.3f %.3f %.1f\n", 
   hitCopy[0]/(double)nCopy[0], hit2Copy[0]/(double)nCopy[0], depthCopy[0]/(double)hit2Copy[0], 
   hitCopy[1]/(double)nCopy[1], hit2Copy[1]/(double)nCopy[1], depthCopy[1]/(double)hit2Copy[1], 
   hitCopy[2]/(double)nCopy[2], hit2Copy[2]/(double)nCopy[2], depthCopy[2]/(double)hit2Copy[2], 
   hitCopy[3]/(double)nCopy[3], hit2Copy[3]/(double)nCopy[3], depthCopy[3]/(double)hit2Copy[3]) ;
}

static void invBuild (Readset *rs) /* run this after ms->depth built from reads->hit[] */
{
  Modset *ms = rs->ms ;
  U32 i, j ;
  U64 offset = 0 ;
  if (rs->inv) bzero (rs->inv, (ms->max+1)*sizeof(U32*)) ;
  else rs->inv = new0 (ms->max+1, U32*) ;
  if (!rs->invSpace) rs->invSpace = new (rs->totHit, U32) ;
  for (i = 1 ; i <= ms->max ; ++i)
    if (ms->depth[i] && ms->depth[i] < U16MAX)
      { rs->inv[i] = rs->invSpace + offset ;
	offset += ms->depth[i] ;
      }
  Read *read = arrp(rs->reads,1,Read) ;
  for (i = 1 ; i < arrayMax(rs->reads) ; ++i, ++read)
    { for (j = 0 ; j < 4 ; ++j) read->nCopy[j] = 0 ; /* build nCopy here in case msCopy changed */
      U32 *x = read->hit ;
      for (j = 0 ; j < read->nHit ; ++j, ++x)
	{ U32 y = *x & TOPMASK ;
	  ++read->nCopy[msCopy(rs->ms,y)] ;
	  if (ms->depth[y] < U16MAX) *rs->inv[y]++ = i ; /* this is the inverse map */
	}
    }
  offset = 0 ; /* now recreate rs->inv because we moved the pointers */
  for (i = 1 ; i <= ms->max ; ++i)
    if (ms->depth[i] && ms->depth[i] < U16MAX)
      { rs->inv[i] = rs->invSpace + offset ;
	offset += ms->depth[i] ;
      }
}

/************************************************************/

typedef struct {
  U32  iy ;			/* index of overlap read in rs->reads */
  U16  nHit ;			/* number of shared mod hits */
  bool isPlus ;			/* relative direction */
  bool isContained ;		/* x is contained in y */
  U16  nBadOrder ;
  U16  nBadFlip ;
} Overlap ;

static int compareOverlap (const void* a, const void* b)
{
  Overlap *oa = (Overlap*)a, *ob = (Overlap*)b ;
  return (ob->nHit - oa->nHit) ;
}

/* in ~58k reads from SK1 
   	around 5.5k have high nDeep - expected ~10% rDNA so OK
	around 5k have low nHits - presumably rubbish, they are not short
	around 2k have lots of bad reads (more than 5 or 10) and are corrupt - many have repeats
	if we remove all these, there are ~1700 with a bad read - probably unfortunate errors
   removing these we get ~42k good reads - lets try to assembly these
*/

Array findOverlaps (Readset *rs, Read *x, int reportLevel) /* returns array of Overlap */
/* reportLevel 0 for none, 1 for per read (RR lines), 2 for per overlap (RR and RH lines) */
{
  static Array omapA = 0, hmapA = 0, xPosA = 0 ; /* avoid many-fold allocation and free */
  int *omap ;  /* map from read2 id to index in olap */
  omapA = arrayReCreate (omapA, arrayMax(rs->reads), int) ; omap = arrp(omapA, 0, int) ;
  U16 *hmap ; /* map from hashes to index in read1 hits */ 
  hmapA = arrayReCreate (hmapA, rs->ms->max+1, U16) ; hmap = arrp(hmapA, 0, U16) ;
  U32 *xPos ; /* convert dx into absolute position NB +1 */
  xPosA = arrayReCreate (xPosA, x->nHit+1, U32) ; xPos = arrp(xPosA, 0, U32) ;

  int j, k ;	/* j is index in list of read hits, k is index over reads containing hxx */
  int nRepeat = 0 ; /* number of repeat mods */
  Array olap = arrayCreate (256, Overlap) ;
  arrayp(olap,0,Overlap)->nHit = 0 ; /* burn the 0 position */
  Overlap *o ;
  U32 *hx = x->hit ; U16 *dx = x->dx ;

  for (j = 0 ; j < x->nHit ; ++j, ++hx, ++dx)
    { U32 hxx = *hx & TOPMASK ;
      xPos[j+1] = xPos[j] + *dx ; /* j+1 not j because hmap[] below is j+1 */
      if (msIsCopy1 (rs->ms, hxx))
	{ if (hmap[hxx]) { ++nRepeat ; x->badRepeat = 1 ; continue ; }
	  hmap[hxx] = j+1 ; /* note the +1 : needed to distinguish from missing */
	  U32 *r2 = rs->inv[hxx] ;
	  for (k = 0 ; k < rs->ms->depth[hxx] ; ++k, ++r2)
	    { if (!omap[*r2])
		{ o = arrayp(olap, omap[*r2] = arrayMax(olap), Overlap) ;
		  o->iy = *r2 ;
		}
	      else o = arrp(olap, omap[*r2], Overlap) ;
	      ++o->nHit ;
	    }
	}
    }

  /* now consider overlap for each read with at least 3 hits */
  int nGood = 0, nBad = 0 ;
  hx = x->hit ;
  arraySort (olap, compareOverlap) ;
  o = arrp(olap, 0, Overlap) ;
  for (k = 1 ; k < arrayMax(olap) ; ++k, ++o) /* NB drop out of loop before reach end */
    { if (o->nHit < 3) break ;
      Read *y = arrp(rs->reads,o->iy,Read) ;
      if (y->bad) continue ;
      int nPlus = 0, nMinus = 0 ;
      U16 ihx ;
      U32 *hy = y->hit ;
      for (j = 0 ; j < y->nHit ; ++j, ++hy) /* run through once to find direction */
	if ((ihx = hmap[*hy & TOPMASK]))
	  { if ((*hy & TOPBIT) == (hx[--ihx] & TOPBIT)) ++nPlus ; else ++nMinus ; }
      hy = y->hit ; /* reset hy for loop to come */
      U16 *dy = y->dx ; double yPos = *dy ;
      if (nPlus > nMinus)
	{ o->isPlus = true ;
	  o->nBadFlip = nMinus ;
	  int last = 0, lastDiff ;
	  for (j = 0 ; j < y->nHit ; ++j, ++hy, yPos += *++dy)
	    if ((ihx = hmap[*hy & TOPMASK]))
	      { lastDiff = xPos[ihx] - yPos ;
		if (!last && lastDiff < 0) o->isContained = true ; /* x starts in y */
		if (ihx < last) { ++o->nBadOrder ; --nPlus ; } last = ihx ;
	      }
	  if (o->isContained && x->len - lastDiff > y->len) o->isContained = false ;
	}
      else if (nMinus && !nPlus)
	{ o->isPlus = false ;
	  o->nBadFlip = nPlus ;
	  int last = x->nHit, lastDiff ;
	  for (j = 0 ; j < y->nHit ; ++j, ++hy, yPos += *++dy)
	    if ((ihx = hmap[*hy & TOPMASK]))
	      { lastDiff = x->len - xPos[ihx] - yPos ;
		if (!last && lastDiff < 0) o->isContained = true ; /* x starts in y */
		if (ihx > last) { ++o->nBadOrder ; --nMinus ; } last = ihx ;
	      }
	  if (o->isContained && x->len - lastDiff > y->len) o->isContained = false ;
	}
      if (o->nBadOrder || o->nBadFlip) ++nBad ; else ++nGood ;
      
      if (reportLevel > 1)
	{ fprintf (outFile, "RH\t%u\tlen %d\t%s\t+ %d\t- %d\tbadOrder %d",
		   o->iy, y->len, (o->nBadOrder + o->nBadFlip) ? "BAD" : "GOOD",
		   nPlus, nMinus, o->nBadOrder) ;
	  fprintf (outFile, "\t%s\n", o->isContained ? "CONTAINED" : "OVERLAP") ;
	}
    }
  arrayMax(olap) = k ;

  if (!nGood && !nBad)
    { x->badNoMatch = 1 ;
      if (x->nHit < 10) x->badLowHit = 1 ;
      else if (x->nCopy[1] < 10) x->badLowCopy1 = 1 ;
    }

  if (reportLevel > 0)
    { I64 iix = x - arrp(rs->reads,0,Read) ;
      if (iix < 0 || iix > arrayMax(rs->reads)) iix = 0 ; /* read not in readset */
      U32 ix = iix ;
      fprintf (outFile, "RR %6u\tlen %d\tnHit %3d\tnMiss %3d\t", ix, x->len, x->nHit, x->nMiss) ;
      fprintf (outFile, "nCpy %d %d %d %d\t", x->nCopy[0], x->nCopy[1], x->nCopy[2], x->nCopy[3]);
      fprintf (outFile, "nRepeatMod %d\tnGood %4d\tnBad %4d\n", nRepeat, nGood, nBad) ;
    }

  return olap ;
}

void printOverlap (Readset *rs, U32 ix, U32 iy)
{
  Read *x = arrp(rs->reads, ix, Read), *y = arrp(rs->reads, iy, Read) ;
  fprintf (outFile, "RR overlaps_for %u\tlen %d\tnHit %d\tnMiss %d\tnCopy %d %d %d %d\n",
	   ix, x->len, x->nHit, x->nMiss, x->nCopy[0], x->nCopy[1], x->nCopy[2], x->nCopy[3]) ;
  fprintf (outFile, "RR overlaps_for %u\tlen %d\tnHit %d\tnMiss %d\tnCopy %d %d %d %d\n",
	   iy, y->len, y->nHit, y->nMiss, y->nCopy[0], y->nCopy[1], y->nCopy[2], y->nCopy[3]) ;
  int j, k ;	/* j is index in list of read 1 hits, k is index in list of read 2 hits */
  int xPos = 0 ;
  int xLast = -1, yLast = -1 ;
  U32 *hx = x->hit ;
  for (j = 0 ; j < x->nHit ; ++j, ++hx)
    { U32 hxx = *hx & TOPMASK ;
      xPos += x->dx[j] ;
      if (msIsCopy1 (rs->ms, hxx))
	{ int yPos = 0 ;
	  U32 *hy = y->hit ;
	  for (k = 0 ; k < y->nHit ; ++k, ++hy)
	    { U32 hyy = *hy & TOPMASK ;
	      yPos += y->dx[k] ;
	      if (hxx == hyy)
		{ bool isSame = ((*hx&TOPBIT) == (*hy&TOPBIT)) ;
		  fprintf (outFile, "RO\t%8x %5d %c\t",
			   hxx, rs->ms->depth[hxx], isSame ? '+' : '-') ;
		  fprintf (outFile, "%u %u %c\t",
			   ix, xPos, (*hx & TOPBIT) ? 'F' : 'R') ;
		  fprintf (outFile, "%u %u %c",
			   iy, yPos, (*hy & TOPBIT) ? 'F' : 'R') ;
		  if (xLast >= 0)
		    { I64 dirn = (xPos-xLast)*(yPos-yLast) ;
		      if ((isSame && dirn < 0) || (!isSame && dirn > 0))
			printf ("\tX xLast %d yLast %d yLen %d", xLast, yLast, y->len) ;
		    }
		  xLast = xPos ; yLast = yPos ;
		  fputc ('\n', outFile) ;
		}
	    }
	}
    }
}

void cluster (Readset *rs)
{
  int i, j ;
  int *link = new0 (arrayMax(rs->reads), int) ; // link to lower read in same cluster
  int nOverlapMade = 0, nNonEmpty = 0 ;
  for (i = 1 ; i < arrayMax(rs->reads) ; ++i)
    if (!link[i])
      { Array olap = findOverlaps (rs, arrp(rs->reads,i,Read), 0) ;
	Overlap *o = arrp(olap,0,Overlap) ;
	int iLink = i ; // the minimal location i links to
	for (j = 1 ; j < arrayMax(olap) ; ++j, ++o) // NB 0'th position is burned
	  { if (o->iy == i) continue ;
	    U32 z = o->iy ;
	    while (link[z])
	      { if (link[z] == iLink) break ;
		z = link[z] ;
	      }
	    if (!link[z])
	      { if (z+1 > iLink) link[z] = iLink ;
		else link[iLink-1] = z ;
	      }
	  }
	++nOverlapMade ;
	if (arrayMax(olap) > 1) ++nNonEmpty ;
	arrayDestroy (olap) ;
      }
  printf ("made %d overlap arrays, of which %d nonEmpty\n", nOverlapMade, nNonEmpty) ;
  int nClus = 0 ;
  int *clus = new0 (arrayMax(rs->reads), int) ; // could reuse link, but this helps readability
  for (i = 1 ; i < arrayMax(rs->reads) ; ++i)
    if (link[i]) clus[i] = clus[link[i]] ;
    else clus[i] = ++nClus ;
  free (link) ;
  int *clusSize = new0 (nClus, int) ;
  for (i = 1 ; i < arrayMax(rs->reads) ; ++i) ++clusSize[clus[i]] ;
  int nProperCluster = 0 ;
  int *properClus = new0 (nClus, int) ;
  for (i = 0 ; i < nClus ; ++i)
    if (clusSize[i] > 1)
      { properClus[i] = ++nProperCluster ;
	printf ("proper cluster %d size %d\n", nProperCluster, clusSize[i]) ;
	clusSize[nProperCluster] = clusSize[i] ;
      }
  printf ("found %d clusters of which %d are proper\n", nClus, nProperCluster) ;
  clusSize[0] = 0 ;
  for (i = 1 ; i < arrayMax(rs->reads) ; ++i) clus[i] = properClus[clus[i]] ;
  free (properClus) ;
  free (clusSize) ;
  free (clus) ;
}

/************************************************************/

void cleanMods (Readset *rs)
{
  Modset *ms = rs->ms ;
  int i, j ;
  Read *r ;
  bool *isInRead = new (ms->max+1, bool) ;
  int w = ms->hasher->w ;

  r = arrp(rs->reads, 0, Read) ;
  for (i = 1 ; i < arrayMax(rs->reads) ; ++i, ++r)
    { bzero (isInRead, ms->max+1) ;
      U32 *h = r->hit ;
      U16 *dx = r->dx ;
      int lastDepth, thisDepth ;
      U32 hh, hhLast ;
      for (j = 0 ; j < r->nHit ; ++j, ++h, ++dx)
	{ U32 hh = *h & TOPMASK ;
	  if (isInRead[hh])  msSetRepeat(ms,hh) ;
	  isInRead[hh] = true ;
	  if (j && *dx < w && j+1 < r->nHit && dx[1] < w) msSetInternal(ms,hh) ;
	  thisDepth = ms->depth[hh] ;
	  if (j)
	    { if (lastDepth > 2*thisDepth) msSetMinor(ms,hh) ;
	      if (thisDepth > 2*lastDepth) msSetMinor(ms,hhLast) ;
	    }
	  lastDepth = thisDepth ;
	  hhLast = hh ;
	}
    }
  free (isInRead) ;

  int nRep = 0, nInt = 0, nMinor = 0 ;
  for (i = 0 ; i < rs->ms->max+1 ; ++i)
    { if (msIsRepeat(ms,i)) ++nRep ;
      if (msIsInternal(ms,i)) ++nInt ;
      if (msIsMinor(ms,i)) ++nMinor ;
    }
  
  invBuild (rs) ;	// need this to rebuild read[i]->nCopy[]

  printf ("set %d repeated, %d internal, %d minor_variant mods\n", nRep, nInt, nMinor) ;
}

/************************************************************/

typedef struct {
  U32 mod ; // use (-2) and (-4) for the two ends
  int dx ;
} Test ;

static bool inline checkMod (Modset *ms, U32 h)
{ static U32 MASK = MS_REPEAT | MS_RDNA ;
  static U32 CHECK = MS_RDNA ;
  return !msIsCopy0(ms,h) && (ms->info[h] & MASK) == CHECK ;
}

static void inline addTest (Array test, Modset *ms, U32 h, int x)
{ h &= TOPMASK ;
  if (checkMod (ms,h))
    { Test *t = arrayp (test, arrayMax(test), Test) ;
      t->mod = h ;
      t->dx = x ;
    }
}

static int testSort (const void *a, const void *b)
{ Test *ta = (Test*)a, *tb = (Test*)b ;
  if (ta->mod < tb->mod) return -1 ;
  if (ta->mod > tb->mod) return 1 ;
  if (ta->dx < tb->dx) return -1 ;
  if (ta->dx > tb->dx) return 1 ;
  return 0 ;
}

static int intSort (const void *a, const void *b)
{ int ia = *(int*)a, ib = *(int*)b ;
  if (ia < ib) return -1 ;
  if (ia > ib) return 1 ;
  return 0 ;
}

void testMods (Readset *rs, int minDepth, int maxDepth)
{
  int i, j, k, kk ;
  Modset *ms = rs->ms ;
  ModInfo *mi = rs->modInfo ;
  Array test = 0 ; // of Test
  Array start = 0, end = 0 ; // will get number of starts, ends of length <= i
  static int RUN = 0 ; ++RUN ;
  static char yName[12], zName[12] ;
  sprintf (yName, "YY-TEST%d", RUN) ;
  FILE *yFile = fopen (yName,"w") ;
  sprintf (zName, "ZZ-TEST%d", RUN) ;
  FILE *zFile = fopen (zName,"w") ;

  if (!mi) die ("need to run -R first") ;
  for (i = 0 ; i < ms->max+1 ; ++i, ++mi)
    mi->nGood = mi->nMod2 = mi->nBadLD = mi->nSplit = mi->nSplitLD = 0 ;

  int nTested = 0 ;
  for (i = 0 ; i < ms->max+1 ; ++i)
    if (ms->depth[i] >= minDepth && ms->depth[i] < maxDepth && checkMod(ms,i))
      { ++nTested ;
	U32 *rj = rs->inv[i] ;
	test = arrayReCreate (test, 4096, Test) ;
	start = arrayReCreate (start, 20000, int) ;
	end = arrayReCreate (end, 20000, int) ;
	for (j = 0 ; j < ms->depth[i] ; ++j, ++rj)
	  { Read *r = arrp(rs->reads, *rj, Read) ;
	    int x = 0 ;
	    int it = arrayMax(test) ;
	    for (k = 0 ; k < r->nHit ; ++k)
	      { x += r->dx[k] ;
		if ((r->hit[k] & TOPMASK) == i)
		  { if (r->hit[k] & TOPBIT) // forward
		      { ++array(start,x,int) ;
			++array(end, r->len - x - ms->hasher->w,int) ;
			while (it < arrayMax(test)) { arrp(test,it,Test)->dx -= x ; ++it ; }
			x = 0 ;
			while (++k < r->nHit) addTest (test, ms, r->hit[k], x += r->dx[k]) ;
		      }
		    else // reversed
		      { ++array(start, r->len - x - ms->hasher->w,int) ;
			++array(end,x,int) ;
			while (it < arrayMax(test))
			  { arrp(test,it,Test)->dx = x - arrp(test,it,Test)->dx ; ++it ; }
			x = 0 ;
			while (++k < r->nHit) addTest (test, ms, r->hit[k], x -= r->dx[k]) ;
		      }
		  }
		else
		  addTest (test, ms, r->hit[k], x) ;
	      }
	  }
	//	arr(end,arrayMax(end)-1, int) = 0 ; arr(start,arrayMax(start)-1, int) = 0 ;
	assert (arr(end,arrayMax(end)-1,int) > 0) ;
	assert (arr(start,arrayMax(start)-1,int) > 0) ;
	if (arrayMax(end) > 1)
	  for (kk = arrayMax(end)-1 ; kk-- ; ) arr(end,kk,int) += arr(end,kk+1,int) ;
	if (arrayMax(start) > 1)
	  for (kk = arrayMax(start)-1 ; kk-- ; ) arr(start,kk,int) += arr(start,kk+1,int) ;
	arraySort (test, testSort) ; 
	Test *t = arrp(test,0,Test) ;
	int nMod = 0, nMod2 = 0, nGood = 0, nSplit = 0 ;
	for (k = 0 ; k < arrayMax(test) ; )
	  { ++nMod ;
	    int n = k, xmin, xmax ;
	    U32 m = t->mod ;
	    if (t->dx > 0)
	      { xmin = t->dx ;
		assert (xmin < arrayMax(end)) ;
		while (t->mod == m && k < arrayMax(test)) { ++k ; ++t ; }
		n = k - n ;
		xmax = (t-1)->dx ;
		if (n < ms->depth[m] && n*2 < arr(end,xmin,int))
		  { ++nMod2 ;
		    if (RUN > 3) ++rs->modInfo[m].nBadLD ;
		  }
		if (n == ms->depth[m] || n >= 0.8*arr(end,xmin,int)) ++nGood ;
		if (n == 1 && arr(end,xmin,int) >= 10) ++rs->modInfo[i].nBadLD ;
		fprintf (zFile,
			 "i %d depth %d m %d depth %d + count %d min %d at %d max %d at %d\n",
			 i, ms->depth[i], m, ms->depth[m], n,
			 arr(end,xmin,int), xmin, arr(end,xmax,int), xmax) ;
	      }
	    else
	      { xmax = -t->dx ;
		while (t->mod == m && k < arrayMax(test)) { ++k ; ++t ; }
		n = k - n ;
		xmin = - (t-1)->dx ;
		if (xmin < 0)
		  { ++nSplit ; ++rs->modInfo[m].nSplitLD ;
		    xmin = xmax ;
		    //		    if (RUN > 1) { printf ("SPLIT %d d %d %d d %d x", i, ms->depth[i], m, ms->depth[m]) ; for (kk = k-n ; kk < k ; ++kk) printf (" %d", t[kk-k].dx) ; printf ("\n") ; }
		  } // shouldn't happen - repeat?
		assert (xmin < arrayMax(start)) ;
		if (xmin < 0) { n = 0 ; xmin = 0 ; }
		if (n < ms->depth[m] && n*2 < arr(start,xmin,int))
		  { ++nMod2 ;
		    if (RUN > 3) ++rs->modInfo[m].nBadLD ;
		  }
		else if (n == 1 && arr(start,xmin,int) >= 10) ++rs->modInfo[m].nBadLD ;
		if (n == ms->depth[m] || n >= 0.8*arr(start,xmin,int)) ++nGood ;
		fprintf (zFile,
			 "i %d depth %d m %d depth %d - count %d min %d at %d max %d at %d\n",
			 i, ms->depth[i], m, ms->depth[m], n,
			 arr(start,xmin,int), xmin, arr(start,xmax,int), xmax) ;
	      }
	  }
	rs->modInfo[i].nGood = nGood ;
	rs->modInfo[i].nMod2 = nMod2 ;
	rs->modInfo[i].nSplit = nSplit ;
      }
  
  mi = rs->modInfo ;
  int nZero1 = 0, nZero2 = 0, nZero3 = 0 ;
  for (i = 0 ; i < ms->max+1 ; ++i, ++mi)
    { if (mi->nGood || mi->nMod2)
	fprintf (yFile, "TEST %d depth %d nGood %d nMod2 %d nBadLD %d nSplit %d\n",
		 i, ms->depth[i], mi->nGood, mi->nMod2, mi->nBadLD, mi->nSplit) ;
      if (mi->nGood < mi->nMod2) { msSetCopy0(ms,i) ; ++nZero1 ; }
      if (mi->nSplit > 10) { msSetCopy0(ms,i) ; ++nZero2 ; }
      if (RUN == 2 || RUN == 6)
	{ if (mi->nBadLD > 20 || mi->nSplitLD > 10)
	    { fprintf (yFile, "BADLD %d depth %d nBadLD %d nSplitLD %d\n",
		       i, ms->depth[i], mi->nBadLD, mi->nSplitLD) ;
	      msSetCopy0(ms,i) ; ++nZero3 ;
	    }
	}
      if (RUN == 3 || RUN == 7)
	{ if (mi->nMod2 > 25) { msSetCopy0(ms,i) ; ++nZero1 ; }
	  if (mi->nSplit) { msSetCopy0(ms,i) ; ++nZero2 ; }
	  if (mi->nBadLD > 10)
	    { fprintf (yFile, "BADLD %d depth %d nBadLD %d nSplitLD %d\n",
		       i, ms->depth[i], mi->nBadLD, mi->nSplitLD) ;
	      msSetCopy0(ms,i) ; ++nZero3 ;
	    }
	}
      if (RUN == 4 || RUN == 8)
	{ if (mi->nBadLD > 6)
	  if (mi->nSplit) { msSetCopy0(ms,i) ; ++nZero2 ; }
	    { fprintf (yFile, "BADLD %d depth %d nBadLD %d nSplitLD %d\n",
		       i, ms->depth[i], mi->nBadLD, mi->nSplitLD) ;
	      msSetCopy0(ms,i) ; ++nZero3 ;
	    }
	}
    }
  printf ("RUN %d tested %d mods and zeroed %d bad>good %d split %d LD\n",
	  RUN, nTested, nZero1, nZero2, nZero3) ;

  invBuild (rs) ;	// need this to rebuild read[i]->nCopy[] after altering copy0
    
  arrayDestroy (test) ;
  fclose (yFile); fclose (zFile) ;
}

/************************************************************/

void refFlag (Readset *rs, char *filename)
{
  int i, j ;
  SeqIO *si = seqIOopenRead (filename, dna2indexConv, false) ; /* false for no qualities */
  if (!si) die ("failed to open ref seq file %s", filename) ;
  Modset *ms = rs->ms ;
  if (!rs->modInfo) rs->modInfo = new0 (ms->max+1, ModInfo) ;
  ModInfo *mi ;
  int *rCount = new0 (ms->max+1, int) ;

  while (seqIOread (si))
    { SeqhashRCiterator *mit = modRCiterator (ms->hasher, sqioSeq(si), si->seqLen) ;
      U64 kmer ; int pos ; U32 index ;
      while (modRCnext (mit, &kmer, &pos, 0))
	if ((index = modsetIndexFind (ms, kmer, false))) // false for do not add
	  { mi = &(rs->modInfo[index]) ;
	    msSetRDNA(ms,index) ; 
	    mi->isRefRDNA = 1 ; mi->rDNApos = pos ;
	    if (ms->depth[index] > 4750) mi->isMultiRDNA = 1 ;
	    else if (ms->depth[index] > 2750) mi->isCoreRDNA = 1 ;
	    else mi->isVarRDNA = 1 ;
	  }
      seqhashRCiteratorDestroy (mit) ;
    }
  seqIOclose (si) ;

  int nRDNAreads = 0 ;
  for (i = 1 ; i < arrayMax(rs->reads) ; ++i)
    { Read *r = arrp(rs->reads, i, Read) ;
      int n = 0, n1 = 0, n100 = 0, n200 = 0, m1 = 0, m100 = 0, m200 = 0 ;
      for (j = 0 ; j < r->nHit ; ++j)
	{ mi = &(rs->modInfo[r->hit[j] & TOPMASK]) ;
	  if (mi->isCoreRDNA && mi->isRefRDNA)
	    { if (!n) n1 = j ;
	      ++n ;
	      if (n == 100) n100 = j ;
	      if (n == 200) { n200 = j ; break ; }
	    }
	}
      if (n200)
	{ n = 0 ;
	  for (j = r->nHit ; --j ; )
	  { mi = &(rs->modInfo[r->hit[j] & TOPMASK]) ;
	    if (mi->isCoreRDNA && mi->isRefRDNA)
	      { if (!n) m1 = j ;
		++n ;
		if (n == 100) m100 = j ;
		if (n == 200) { m200 = j ; break ; }
	      }
	  }
	}
      if (m200 > n200)
	{ int lastPos = 0 ;
	  for (j = n200 ; j < m200 ; ++j)
	    { U32 h = r->hit[j] & TOPMASK ;
	      mi = &(rs->modInfo[h]) ;
	      if (mi->isRDNA)
		{ int p = mi->rDNApos ;
		  if (mi->isRefRDNA)
		    lastPos = p ;
		  else if (p > 0 && p < lastPos + 50 && p > lastPos - 50)
		    { mi->rDNApos = (rCount[h] * p + lastPos) / (rCount[h] + 1) ;
		      ++rCount[h] ;
		    }
		  else
		    mi->rDNApos = -1 ;
		}
	      else
		{ msSetRDNA(ms,h) ;
		  if (ms->depth[h] > 4750) mi->isMultiRDNA = 1 ;
		  else if (ms->depth[h] > 2750) mi->isCoreRDNA = 1 ;
		  else mi->isVarRDNA = 1 ;
		  mi->rDNApos = lastPos ;
		  rCount[h] = 1 ;
		}
	    }
	  r->isRDNA = 1 ;
	  ++nRDNAreads ;
	}
      //      if (n1 > 100 || m1 < r->nHit-100)	printf ("  read %d nmod %d nrefcore %d n1 %d n100 %d n200 %d m200 %d m100 %d m1 %d\n", i, r->nHit, n, n1, n100, n200, m200, m100, m1) ;
    }

  int nRDNA = 0, nRef = 0, nGoodPos = 0 ;
  int nRefC = 0, nRefV0 = 0, nRefV1 = 0, nRefM = 0, nOthC = 0, nOthV0 = 0, nOthV1 = 0, nOthM = 0 ;
  mi = rs->modInfo ;
  for (i = 0 ; i < ms->max + 1 ; ++i, ++mi)
    if (mi->isRDNA)
      { ++nRDNA ;
	if (mi->isRefRDNA)
	  { ++nRef ;
	    if (mi->isCoreRDNA) ++nRefC ;
	    else if (mi->isMultiRDNA) ++nRefM ;
	    else if (msIsCopy0(ms,i)) ++nRefV0 ;
	    else ++nRefV1 ;
	  }
	else
	  { if (mi->isCoreRDNA) ++nOthC ;
	    else if (mi->isMultiRDNA) ++nOthM ;
	    else if (msIsCopy0(ms,i)) ++nOthV0 ;
	    else ++nOthV1 ;
	    if (mi->rDNApos > 0) ++nGoodPos ;
	  }
      }
  printf ("total nRDNAreads %d other reads %d\n", nRDNAreads, arrayMax(rs->reads)-1-nRDNAreads) ;
  printf ("total nRDNAmods %d nRDNAref %d other mods %d\n", nRDNA, nRef, ms->max+1-nRDNA) ;
  printf ("  nRefC %d nRefM %d nRefVcopy>0 %d nRefVcopy0 %d\n", nRefC, nRefM, nRefV1, nRefV0) ;
  printf ("  nOthC %d nOthM %d nOthVcopy>0 %d nOthVcopy0 %d", nOthC, nOthM, nOthV1, nOthV0) ;
  printf (" nGoodPos %d\n", nGoodPos) ;
}

/* next function contains adhoc rules to reset bits to support further analyses */

void resetBits (Readset *rs, int op)
{
  int i, n = 0 ;
  Modset *ms = rs->ms ;
  ModInfo *mi = rs->modInfo ;
  
  switch (op)
    {
    case 1: printf ("resetting rDNA core kmers to copy1, rest to copy0:") ;
      for (i = 0 ; i < ms->max + 1 ; ++i, ++mi)
	if (mi->isCoreRDNA)
	  { msSetCopy1 (ms, i) ; ++n ; }
	else
	  msSetCopy0 (ms, i) ;
      printf (" %d kept\n", n) ;
      break ;
    case 2: printf ("resetting non-repetitive rDNA core kmers to copy1, rest to copy0:") ;
      for (i = 0 ; i < ms->max + 1 ; ++i, ++mi)
	if (mi->isCoreRDNA && !msIsRepeat(ms,i))
	  { msSetCopy1 (ms, i) ; ++n ; }
	else
	  msSetCopy0 (ms, i) ;
      printf (" %d kept\n", n) ;
      break ;
    case 3: printf ("resetting rDNA core kmers not repeated in read 1 to copy1: ") ;
      for (i = 0 ; i < ms->max + 1 ; ++i, ++mi)
	if (mi->isCoreRDNA)
	  { msSetCopy1 (ms, i) ; ++n ; }
	else
	  msSetCopy0 (ms, i) ;
      bool *z = new0(ms->max+1,bool) ;
      Read *r1 = arrp(rs->reads,1,Read) ;
      for (i = 0 ; i < r1->nHit ; ++i)
	{ U32 h = r1->hit[i] & TOPMASK ;
	  if (!msIsCopy1(ms,h)) continue ;
	  if (z[h]) { msSetCopy0(ms,h) ; --n ; }
	  else z[h] = true ;
	}
      printf (" %d kept\n", n) ;
      free (z) ;
      break ;
    }
  
  invBuild (rs) ;		/* need to rebuild with new copy numbers */
}

/************************************************************/

void readProperties (Readset *rs) // studies properties of reads in terms of copy 1 mods
{
  int i, j ;
  U32 h ;
  Read *read = arrp(rs->reads,1,Read) ;
  Modset *ms = rs->ms ;
  int *f = new (ms->max+1, int) ;
  int *r = new (ms->max+1, int) ;
  
  for (i = 1 ; i < arrayMax(rs->reads) ; ++i, ++read)
    { memset (f, 0, (ms->max+1)*sizeof(int)) ;
      memset (r, 0, (ms->max+1)*sizeof(int)) ;
      for (j = 0 ; j < read->nHit ; ++j)
	{ h = read->hit[j] & TOPMASK ;
	  if (!msIsCopy1(ms,h)) continue ;
	  if (read->hit[j] & TOPBIT) ++f[h] ;
	  else ++r[h] ;
	}
      int n = 0, n2Rev = 0, n2Tan = 0, nMoreTan = 0, nMoreRev = 0 ;
      for (h = 0 ; h < ms->max+1 ; ++h)
	if (f[h] + r[h])
	  { ++n ;
	    if (f[h] + r[h] == 1) continue ;
	    if (f[h] == 1 && r[h] == 1) ++n2Rev ;
	    else if ((f[h] == 2 && r[h] == 0) || (f[h] == 0 && r[h] == 2)) ++n2Tan ;
	    else if (f[h] > 0 && r[h] > 0) ++nMoreRev ;
	    else
	      { ++nMoreTan ;
		printf ("MT i %d h %d count %d\n", i, h, f[h]+r[h]) ;
	      }
	  }
      printf ("READ %d n %d n2Tan %d n2Rev %d nMoreTan %d nMoreRev %d\n",
	      i, n , n2Tan, n2Rev, nMoreTan, nMoreRev) ;
      if (nMoreTan > 5)
	{ printf ("RM %d nMoreTan %d", i, nMoreTan) ;
	  for (h = 0 ; h < ms->max+1 ; ++h) if (f[h] + r[h] > 2) printf (" %d", h) ;
	  putchar ('\n') ;
	}
    }
  free (f) ; free (r) ;
}

/************************************************************/

typedef struct {
  U32 from, to ; // hits, i.e. with TOPBIT set if forward, use 0 for end of read
  U32 i, x ;     // i is the read index, x is the position of ->to in it
} Link ;

int compareLink (const void *a, const void *b)
{
  Link *la = (Link*) a, *lb = (Link*) b ;

  if (la->from < lb->from) return -1 ;
  if (la->from > lb->from) return 1 ;
  if (la->to < lb->to) return -1 ;
  if (la->to > lb->to) return 1 ;
  if (la->i < lb->i) return -1 ;
  if (la->i > lb->i) return 1 ;
  if (la->x < lb->x) return -1 ;
  if (la->x > lb->x) return 1 ;
  die ("problem in compareLink") ; // shouldn't get here!
  return 0 ;
}

static inline char* modText (Readset *rs, U32 h, bool isReverse)
{
  static char buf[64] ;
  int m = (int)(h & TOPMASK) ;
  ModInfo *mi = &rs->modInfo[m] ;
  if (!(h & TOPBIT)) isReverse = !isReverse ;
  if (mi->isRefRDNA)
    sprintf (buf, "%d %c d %d C%d P %d",
	     m, isReverse ? 'R' : 'F', (int)rs->ms->depth[m], msCopy(rs->ms,m), mi->rDNApos) ;
  else
    sprintf (buf, "%d %c d %d C%d p %d",
	     m, isReverse ? 'R' : 'F', (int)rs->ms->depth[m], msCopy(rs->ms,m), mi->rDNApos) ;
    
  return buf ;
}

typedef struct {
  int  read ;
  int  start, end ;
  bool isForward ;
  int  nHit ;
} Layout ;

static int compareLayout (const void *a, const void *b)
{ return ((Layout*)a)->start - ((Layout*)b)->start ; }

typedef struct {
  int  iRead, iLayout ;
  int  x, dx ;
} Active ;

static int inline addActive (HASH hActive, Array active, Array layout, int i, int x, int offset)
{
  int n ;
  hashAdd (hActive, HASH_INT(i), &n) ;
  Active *a = arrayp (active, n, Active) ;
  a->iRead = i ; 
  a->iLayout = arrayMax(layout) ;
  a->x = x ;
  printf ("  added %d x %d\n", i, x) ;
  Layout *y = arrayp(layout, a->iLayout, Layout) ;
  y->read = i ;
  y->start = offset-x ;
  return n ;
}

static int compareInt (const void *a, const void *b)
{ return *(int*)a - *(int*)b ; }

static void assembleFrom (Readset *rs, Array links,
			  U32 from, int offset, bool isReverse,
			  int *iForward, int *iReverse)
{
  int     i, ia ;
  Link   *l ;
  HASHKEY hk ;
#define lStart(h) arrp(links, ((h) & TOPBIT) ? iForward[(h) & TOPMASK] : iReverse[h], Link)

  Array layout = arrayCreate (1024, Layout) ;
  Array active = arrayCreate (64, Active) ; // list of active reads
  HASH  hActive = hashCreate (4096) ; // plenty of space so few clashes
  Array dd = arrayCreate (64, int) ; // used to find median value of d when !isBestUniform

  // code to initialise xFrom assumes from is single copy, i.e. 1 or 0 copies per read
  hashStats () ;
  HASH hash = hashCreate (64) ;
  for (l = lStart(from) ; l->from == from ; ++l)
     if (l->to) // almost always
      hashAdd (hash, HASH_INT(l->to), &ia) ;
    else // have to look for from in read - ugly but rare
      { Read *r = arrp(rs->reads,l->i,Read) ;
	int x = 0 ;
	for (i = 0 ; i < r->nHit ; ++i)
	  { x += r->dx[i] ;
	    if ((r->hit[i] & TOPMASK) == (from & TOPMASK))
	      { if ((r->hit[i] & TOPBIT) != (from & TOPBIT)) x = r->len - x ;
		addActive (hActive, active, layout, l->i, x, offset) ;
		break ;
	      }
	  }
      }
  hashStats () ;
  hashInitIterator (hash) ;
  while (hashNextKeyValue (hash, &hk, 0)) // run through mods that follow 'from'
    { U32 to = (U32)(HASH_INT(hk.i).i) ^ TOPBIT ; // NB HASH_INT() is self-inverse
      for (l = lStart(to) ; l->from == to ; ++l)
	if (l->to == (from ^ TOPBIT))
	  addActive (hActive, active, layout, l->i, arrp(rs->reads,l->i,Read)->len-l->x, offset) ;
    }
  hashDestroy (hash) ;

  // main loop will update from, i.e. move assembly along by one mod
  while (true) 
    { U32  bestTo = 0, lastTo = 0 ;
      int  dBest = 0, nBest = 0 ;
      bool isBestUniform ;
      int  ia, d, dMin = 0, dSum = 0, nLast = 0, iLast = -1 ;

      printf ("FROM %s pos %d active %d",
	      modText (rs,from,isReverse), offset, hashCount (hActive)) ;

      hashInitIterator (hActive) ; // clear "active"->dx
      while (hashNextKeyValue (hActive, &hk, &ia)) arrp(active, ia, Active)->dx = 0 ;

      for (l = lStart(from) ; l->from == from ; ++l)
	if (hashFind (hActive, HASH_INT(l->i), &ia)) // only consider active reads
	  { Active *a = arrp (active, ia, Active) ;
	    d = l->x - a->x ;

	    if (isVerbose)
	      { printf ("\n  TO %s i %d x %d dx %d", modText (rs,l->to,isReverse), l->i, l->x, d) ;
		if (l->to == 0) printf (" end %d", l->i) ;
	      }
	      
	    if (l->to != lastTo)
	      { if (lastTo && 2*nLast > hashCount(hActive) && (!dBest || dMin < dBest))
		  { dBest = dMin ; bestTo = lastTo ; nBest = nLast ;
		    isBestUniform = (dSum == nBest*dBest) ;
		  }
		lastTo = l->to ; nLast = 0 ; iLast = -1 ; dMin = 0 ; dSum = 0 ;
	      }
	  
	    if (d > 0 && l->i != iLast)
	      { ++nLast ; iLast = l->i ; dSum += d ;
		if (dMin == 0 || d < dMin) dMin = d ;
		a->dx = d ;

		Layout *y = arrp(layout, a->iLayout, Layout) ;
		++y->nHit ; printf (" hit %d", y->nHit) ;
		y->end = offset - l->x ; // will add the read length at the end
	      }
	  }
      if (lastTo && 2*nLast > hashCount(hActive) && (!dBest || dMin < dBest)) // check final
	{ dBest = dMin ; bestTo = lastTo ; nBest = nLast ;
	  isBestUniform = (dSum == nBest*dBest) ;
	}
      if (isVerbose) putchar ('\n') ;

      if (!nBest) break ;  // loop exit - insufficient support
      
      // Update xFrom. Need to add new reads here. How do I avoid adding them on repeats?
      // I think I can only add at copy1 SNPs. Then do I fill in backwards? Not for now.
      if (isBestUniform) // standard case where all the deltas agree
	{ hashInitIterator (hActive) ;
	  while (hashNextKeyValue (hActive, &hk, &ia))
	    { Active *a = arrp(active,ia,Active) ;
	      a->x += dBest ;
	      if (a->x > arrp(rs->reads,a->iRead,Read)->len)
		{ hashRemove (hActive, hk) ;
		  printf ("\nEND %d pos %d end %d\n", a->iRead, offset,
			  arrp(rs->reads,a->iRead,Read)->len +
			  arrp(layout,a->iLayout,Layout)->end) ;
		}
	    }
	}
      else
	{ arrayMax(dd) = 0 ; hashInitIterator (hActive) ;  // set dBest to the median dx
	  while (hashNextKeyValue (hActive, &hk, &ia))
	    { Active *a = arrp(active,ia,Active) ;
	      if (a->dx) array(dd,arrayMax(dd),int) = a->dx ;
	    }
	  arraySort (dd, compareInt) ;
	  dBest = arr(dd, nBest/2, int) ;
	  
	  hashInitIterator (hActive) ;
	  while (hashNextKeyValue (hActive, &hk, &ia))
	    { Active *a = arrp(active,ia,Active) ;
	      if (!a->dx || a->dx == dBest)
		a->x += dBest ; // default
	      else if (a->dx > dBest-10 && a->dx < dBest+10) // things just a little off
		{ printf (" dx %d %d", i, a->dx - dBest) ; a->x += a->dx ; }
	      else
		{ printf (" xx %d %d", i, a->dx - dBest) ; a->x += a->dx ; --nBest ; }
	      if (a->x > arrp(rs->reads,a->iRead,Read)->len)
		{ hashRemove (hActive, hk) ;
		  printf ("\nEND %d pos %d end %d\n", a->iRead, offset,
			  arrp(rs->reads,a->iRead,Read)->len +
			  arrp(layout,a->iLayout,Layout)->end) ;
		}
	    }
	}
      if (msIsCopy1 (rs->ms, bestTo & TOPMASK)) // add any new reads to active set
	{ for (l = lStart(from) ; l->to < bestTo ; ++l) ;
	  for (iLast = -1 ; l->from == from && l->to == bestTo ; ++l)
	    if (!hashFind (hActive, HASH_INT(l->i), 0))
	      addActive (hActive, active, layout, l->i, l->x, offset) ;
	}

      printf (" BEST %s nBest %d dBest %d", modText(rs,bestTo,isReverse), nBest, dBest) ;

      // finally terminate the line and update major loop variables
      putchar ('\n') ;
      from = bestTo ;
      if (isReverse) offset -= dBest ;
      else offset += dBest ;
    }
  printf ("\nDONE\n") ;

  arraySort (layout, compareLayout) ;
  for (i = 0 ; i < arrayMax(layout) ; ++i)
    { Layout *y = arrp(layout,i,Layout) ;
      Read *r = arrp(rs->reads,y->read,Read) ;
      y->end += r->len ;
      printf ("LAYOUT %d start %d end %d n %d / %d\n",
	      y->read, y->start, y->end, y->nHit, r->nHit) ;
    }

  free (dd) ;
  arrayDestroy (active) ; arrayDestroy (layout) ;
  hashDestroy (hActive) ;
}

void assembleFromMod (Readset *rs, U32 seed, int offset)
{
  int i, j ;
  U32 *h ;
  Modset *ms = rs->ms ;
  Link *l ;

  printf ("assembling mod %d depth %d\n", seed, ms->depth[seed]) ;
  if (!msIsCopy1(ms,seed)) die ("seed copy number %d != 1", msCopy(ms,seed)) ;

  Array reads = arrayCreate (1024, U32) ;
  for (i = 0 ; i < ms->depth[seed] ; ++i) // simple set of reads that contain the seed
    array (reads, arrayMax(reads), U32) = rs->inv[seed][i] ;
  
  // next build an array of all the links seen in the reads
  Array links = arrayCreate (100000, Link) ;
  for (i = 0 ; i < arrayMax(reads) ; ++i)
    { U32 ir = arr(reads, i, U32) ;
      Read *read = arrp(rs->reads, ir, Read) ;
      int len = read->len ;
      int x = 0, xLast ;
      U32 last = 0 ;
      for (j = 0, h = read->hit ; j < read->nHit ; ++j, ++h)
	{ x += read->dx[j] ;
	  if (!msIsCopy0 (ms, *h & TOPMASK))
	    { l = arrayp(links, arrayMax(links), Link) ; l->i = ir ;
	      l->from = *h ^ TOPBIT ; l->to = 0 ; l->x = len ;
	      last = *h ; xLast = x ;
	      break ;
	    }
	}
      for ( ++j, ++h ; j < read->nHit ; ++j, ++h)
	{ x += read->dx[j] ;
	  if (!msIsCopy0 (ms, *h & TOPMASK))
	    { l = arrayp(links, arrayMax(links), Link) ; l->i = ir ;
	      l->from = last ; l->to = *h ; l->x = x ;
	      l = arrayp(links, arrayMax(links), Link) ; l->i = ir ;
	      l->from = *h ^ TOPBIT ; l->to = last ^ TOPBIT ; l->x = len - xLast ;
	      last = *h ; xLast = x ;
	    }
	}
      if (last)
	{ l = arrayp(links, arrayMax(links), Link) ; l->i = ir ;
	  l->from = last ; l->to = 0 ; l->x = len ;
	}
    }
  arraySort (links, compareLink) ;

  // now run through and find start points for each mod and its reverse
  int *iForward = new0 (ms->max+1, int), *iReverse = new0 (ms->max+1, int) ;
  U32 last = 0, to ;
  for (i = 0, l = arrp (links, 0, Link) ; i < arrayMax(links) ; ++i, ++l)
    if (l->from != last)
      { if (l->from & TOPBIT) iForward[l->from & TOPMASK] = i ;
	else iReverse[l->from] = i ;
	last = l->from ;
      }
  array(links,arrayMax(links),Link).from = U32MAX ; // terminator to simplify loop tests below

  // then start from seed and build forwards
  assembleFrom (rs, links, seed | TOPBIT, offset, false, iForward, iReverse) ;
  // and then do the same in reverse
  // assembleFrom (rs, links, seed, offset, true, iForward, iReverse) ;

  arrayDestroy (reads) ; arrayDestroy (links) ;
  free (iForward) ; free (iReverse) ;
}

/************************************************************/

/* mark bad reads based on one of a set of criteria
   most complex is lack of consistent alignment to another read
   complex because a priori it is not clear which of the two is bad
   sort out in three phases: first egregious ones, which have many bad overlaps
   then ones with multiple bad overlaps, then the remainder with just one
*/

void markBadReads (Readset *rs)	/* new version just makes overlaps once */
{
  int i, ix ;
  Read *x, *x1 = arrp(rs->reads,1,Read) ;

  for (ix = 1, x = x1 ; ix < arrayMax(rs->reads) ; ++ix, ++x)
    x->bad = 0 ; /* initialise by clearing all bad flags */

  int *badList = new0 (arrayMax(rs->reads)*10, int) ;
  int *nBad = new0 (arrayMax(rs->reads), int), *lBad = new0 (arrayMax(rs->reads), int) ;
  
  /* pass through findOverlaps and record how many bad overlaps in bady and 
     the bad overlaps themselves to reads possibly still good in badList
   */
  for (ix = 1, x = x1 ; ix < arrayMax(rs->reads) ; ++ix, ++x)
    { Array olap = findOverlaps (rs, x, 0) ;
      Overlap *o = arrp(olap, 0, Overlap) ;
      for (i = 0 ; i < arrayMax(olap) ; ++i, ++o)
	if (o->nBadFlip || o->nBadOrder)
	  { int iy = o->iy ;
	    ++nBad[iy] ;
	    if (nBad[iy] < 10 && lBad[ix] < 10)
	      badList[10*ix + lBad[ix]++] = iy ;
	  }
      arrayDestroy (olap) ;
    }

  /* first pass, set badOrder10 for all 10 or more - these are clearly bad */
  int N = 0 ;
  for (ix = 1, x = x1 ; ix < arrayMax(rs->reads) ; ++ix, ++x)
    if (nBad[ix] >= 10) { x->badOrder10 = 1 ; ++N ; lBad[ix] = 0 ; }
  printf ("MB  %d with >=10 bad overlaps\n", N) ;

  for (ix = 1, x = x1 ; ix < arrayMax(rs->reads) ; ++ix, ++x) /* remove bad from lists */
    for (i = lBad[ix] ; i-- ; )
      if (arrp(rs->reads, badList[10*ix+i], Read)->bad)
	badList[10*ix+i] = badList[10*ix + --lBad[ix]] ;
  
  /* second pass, set badOrder1 for all 2 or more - removes some more singleton errors */
  N = 0 ;
  for (ix = 1, x = x1 ; ix < arrayMax(rs->reads) ; ++ix, ++x)
    if (lBad[ix] >= 2) { x->badOrder1 = 1 ; ++N ; lBad[ix] = 0 ; }
  printf ("MB  %d with multiple bad overlaps\n", N) ;

  for (ix = 1, x = x1 ; ix < arrayMax(rs->reads) ; ++ix, ++x) /* remove bad from lists */
    for (i = lBad[ix] ; i-- ; )
      if (arrp(rs->reads, badList[10*ix+i], Read)->bad)
	badList[10*ix+i] = badList[10*ix + --lBad[ix]] ;

  /* final pass, mark anything left with any bad overlaps as bad */
  N = 0 ;
  for (ix = 1, x = x1 ; ix < arrayMax(rs->reads) ; ++ix, ++x)
    if (lBad[ix] > 0) { x->badOrder1 = 1 ; ++N ; lBad[ix] = 0 ; }
  printf ("MB  %d with single bad overlaps\n", N) ;

  free (badList) ; free (nBad) ; free (lBad) ;
}

int badOverlaps (Readset *rs, Read *x)
{
  int i, nBad = 0 ;
  Array olap = findOverlaps (rs, x, 0) ;
  Overlap *o = arrp(olap,0,Overlap) ;
  for (i = 0 ; i < arrayMax(olap) ; ++i, ++o)
    if (o->nBadOrder || o->nBadFlip) ++nBad ;
  arrayDestroy (olap) ;
  return nBad ;
}

void markBadReadsOld (Readset *rs)
{
  int ix ;
  Read *x, *x1 = arrp(rs->reads,1,Read) ;

  for (ix = 1, x = x1 ; ix < arrayMax(rs->reads) ; ++ix, ++x)
    x->bad = 0 ; /* initialise by clearing all bad flags */
  
  /* first pass through findOverlaps - set badOrder10 */
  int nBad = 0 ;
  for (ix = 1, x = x1 ; ix < arrayMax(rs->reads) ; ++ix, ++x)
    if (badOverlaps (rs, x) >= 10) { x->badOrder10 = 1 ; ++nBad ; }
  printf ("MB  %d with >=10 bad overlaps\n", nBad) ;


  nBad = 0 ;
  for (ix = 1, x = x1 ; ix < arrayMax(rs->reads) ; ++ix, ++x)
    if (badOverlaps (rs, x) > 1 && !x->badOrder10) { x->badOrder1 = 1 ; ++nBad ; }
  printf ("MB  %d with multiple bad overlaps\n", nBad) ;

  /* third pass, set badOrder1 for anything bad  - must be simple pairs - might be too harsh */
  nBad = 0 ;
  for (ix = 1, x = x1 ; ix < arrayMax(rs->reads) ; ++ix, ++x)
    if (badOverlaps (rs, x) > 0 && !x->badOrder10) { x->badOrder1 = 1 ; ++nBad ; }
  printf ("MB  %d with single bad overlaps\n", nBad) ;
}


/************************************************************/

/* identify containing read with maximal number of hits 
   could look for high hit:miss rate 
   and expected:observed overlap hit count (to avoid phase switching) 
*/

void markContained (Readset *rs) 
{
  int ix, io ;
  Read *x, *x1 = arrp(rs->reads,1,Read) ;
  int nContained = 0, nNotContained = 0 ;
  U64 totLen = 0 ;
  Array overlaps = 0 ;
  
  for (ix = 1, x = x1 ; ix < arrayMax(rs->reads) ; ++ix, ++x)
    { if (overlaps) { arrayDestroy (overlaps) ; overlaps = 0 ; }
      if (x->bad) continue ;
      overlaps = findOverlaps (rs, x, 0) ;
      int maxHit = 0 ;
      for (io = 0 ; io < arrayMax (overlaps) ; ++io)
	{ Overlap *o = arrp(overlaps,io,Overlap) ;
	  if (o->iy == ix) continue ; /* don't allow self-containment! */
	  Read *y = arrp(rs->reads,o->iy,Read) ;
	  if (!o->isContained || o->nHit <= maxHit) continue ;
	  x->contained = o->iy ; maxHit = o->nHit ;
	}
      if (x->contained) ++nContained ; else { ++nNotContained ; totLen += x->len ; }
    }
  printf ("MC  found %d contained reads, leaving %d not contained, av length %.1f\n",
	  nContained, nNotContained, nNotContained ? totLen/(double)nNotContained : 0.) ;
}

/************************************************************/

/* initial assembly strategy is to find all overlaps of a read, find all mods in them, 
   and put them in order, with estimated spacing, to make a new assemblyFragment, 
   which is a type of consensus read
*/

typedef struct {
  U32 hit ;
  U32 count ;
  int pos ;
  int upCount ;			/* count of upstream hits  */
  Array downHits ;		/* list of downstream hits */
} AssemblyHit ;

void assembleFromRead (Readset *rs, U32 ix)
{
  Read *x = arrp(rs->reads,ix,Read) ;
  Read *z = new0 (1,Read) ;	/* build the assembly here */
  int io, iy, ih ;

  /* first collect all the hits in all the overlapping reads */
  Array aHits = arrayCreate (1024, AssemblyHit) ;
  AssemblyHit *ah ;
  HASH  hitHash = hashCreate (1024) ;
  Array overlaps = findOverlaps (rs, x, 1) ;
  for (io = 0 ; io < arrayMax (overlaps) ; ++io)
    { Overlap *o = arrp(overlaps,io,Overlap) ;
      Read *y = arrp(rs->reads,o->iy,Read) ;
      int yPos = 0 ;
      Array lastDown = 0 ;
      if (o->isPlus)
	{
	  for (iy = 0 ; iy < y->nHit ; ++iy)
	    { U32 hit = y->hit[iy] & TOPMASK ; yPos += y->dx[iy] ;
	      hashAdd (hitHash, HASH_INT(hit), &ih) ;
	      ah = arrayp(aHits,ih,AssemblyHit) ;
	      if (!ah->count)
		{ ah->hit = hit ;
		  ah->downHits = arrayCreate (8, int) ;
		}
	      ++ah->count ;
	      if (iy) ++ah->upCount ;
	      if (lastDown) array(lastDown, arrayMax(lastDown), int) = ih ;
	      lastDown = ah->downHits ;
	    }
	}
      else
	{
	}
    }

  /* now want to select the longest clearly supported chain */
  /* need a full chain because we want to rely on the order for future matches */

  Array poHits = arrayCreate (arrayMax(aHits), int) ;
  for (ih = 0, ah = arrp(aHits,0,AssemblyHit) ; ih < arrayMax(aHits) ; ++ih, ++ah)
    if (!ah->upCount) array(poHits, arrayMax(poHits), int) = ih ;
  
  /* reporting */
  
  double totCount = 0. ;
  int i, j, countA[20][20], countB[20][20] ;
  for (i = 20 ; i-- ;) for (j = 20 ; j-- ;) { countA[i][j] = 0 ; countB[i][j] = 0 ; }
  for (ih = 0 ; ih < hashCount(hitHash) ; ++ih)
    { ah = arrp(aHits,ih,AssemblyHit) ;
      ah->pos /= ah->count ;	/* just an estimate for now */
      totCount += ah->count ;
      if (!msIsCopy1 (rs->ms, ah->hit)) continue ;
      i = ah->count ; if (i > 19) i = 19 ;
      j = rs->ms->depth[ah->hit] ; if (j > 19) j = 19 ; ++countA[i][j] ;
      j = (10*ah->count - 1) / rs->ms->depth[ah->hit] ; ++countB[i][j] ;
    }
  totCount /= hashCount(hitHash) ;
  printf ("AR  %d total hits - mean count %.1f\n", hashCount(hitHash), totCount) ;
  for (i = 0 ; i < 20 ; ++i)
    { printf ("AH  %2d\t", i) ;
      for (j = 0 ; j < 20 ; ++j) if (j < i) printf ("    ") ; else printf ("%4d", countA[i][j]) ;
      printf ("    ") ;
      for (j = 0 ; j < 10 ; ++j) printf ("%4d", countB[i][j]) ;
      printf ("\n") ;
    }

  hashDestroy (hitHash) ;
  arrayDestroy (aHits) ;
  arrayDestroy (overlaps) ;
}

/************************************************************/

void usage (void)
{ fprintf (stderr, "Usage: modasm <commands>\n") ;
  fprintf (stderr, "Commands are executed in order - set parameters before using them!\n") ;
  fprintf (stderr, "  -v | --verbose : toggle verbose mode\n") ;
  fprintf (stderr, "  -t | --threads <number of threads for parallel ops> [%d]\n", numThreads) ;
  fprintf (stderr, "  -o | --output <output filename> : '-' for stdout\n") ;
  fprintf (stderr, "  -m | --modset <mod file>\n") ;
  fprintf (stderr, "  -f | --seqfile <file of reads: fasta/q, can be gzipped, or binary>\n") ;
  fprintf (stderr, "  -w | --write <file stem> : writes assembly files\n") ;
  fprintf (stderr, "  -r | --read <file stem> : read assembly files\n") ;
  fprintf (stderr, "  -S | --stats : give readset stats\n") ;
  fprintf (stderr, "  -o1 | --overlap1 <read> : find overlaps for given read\n") ;
  fprintf (stderr, "  -o2 | --overlap2 <k> : give overlap stats for every k'th read\n") ;
  fprintf (stderr, "  -o3 | --overlap3 <read1> <read2> : print details of overlap\n") ;
  fprintf (stderr, "  -b | --markBadReads : identify and categorise bad reads\n") ;
  fprintf (stderr, "  -c | --markContained : identify contained reads\n") ;
  fprintf (stderr, "  -a1 | --assemble1 <read> : assemble starting from given read\n") ;
  fprintf (stderr, "  -a2 | --assemble2 <mod> : assemble starting from given mod\n") ;
  fprintf (stderr, "  -u | --cluster : single linkage cluster reads using good overlaps\n") ;
  fprintf (stderr, "  -C | --cleanmods : set repeat and minor allele flags\n") ;
  fprintf (stderr, "  -T | --testmods <minDepth> <maxDepth> : set copy0 if not read-LD consistent\n") ;
  fprintf (stderr, "  -R | --ref <ref seq file> : set rDNA info\n") ;
  fprintf (stderr, "  -rb | --resetbits <n> : various cookery operations - see code\n") ;
  fprintf (stderr, "  -P | --readProperties : info about reads\n") ;
}

int main (int argc, char *argv[])
{
  --argc ; ++argv ;		/* eat program name */

  outFile = stdout ;
  
  timeUpdate (stdout) ;		/* initialise timer */
#ifdef OMP
  numThreads = omp_get_max_threads () ;
  omp_set_num_threads (numThreads) ;
#endif

  if (!argc) usage () ;

  int i ;			/* generically useful variables */
  FILE *f ;
  Modset *ms = 0 ;
  Readset *rs = 0 ;

  while (argc) {
    if (**argv != '-')
      die ("option/command %s does not start with '-': run without arguments for usage", *argv) ;
    //fprintf (stderr, "COMMAND %s", *argv) ;
    //for (i = 1 ; i < argc && *argv[i] != '-' ; ++i) fprintf (stderr, " %s", argv[i]) ;
    //fputc ('\n', stderr) ;
    
#define ARGMATCH(x,y,n)	((!strcmp (*argv, x) || (!strcmp (*argv,y))) && argc >= n && (argc -= n, argv += n))
    if (ARGMATCH("-t","--threads",2))
      {
#ifdef OMP
	numThreads = atoi(argv[-1]) ;
	if (numThreads > omp_get_max_threads ()) numThreads = omp_get_max_threads () ;
	omp_set_num_threads (numThreads) ;
#else
	fprintf (stderr, "  can't set thread number - not compiled with OMP\n") ;
#endif
      }
    else if (ARGMATCH("-v","--verbose",1)) isVerbose = !isVerbose ;
    else if (ARGMATCH("-o","--output",2))
      { if (!strcmp (argv[-1], "-"))
	  outFile = stdout ;
	else if (!(outFile = fopen (argv[-1], "w")))
	  { fprintf (stderr, "can't open output file %s - resetting to stdout\n", argv[-1]) ;
	    outFile = stdout ;
	  }
      }
    else if (ARGMATCH("-m","--modset",2))
      { if (!(f = fzopen (argv[-1], "r"))) die ("failed to open mod file %s", argv[-1]) ;
	if (ms) modsetDestroy (ms) ;
	ms = modsetRead (f) ;
	if (ms->max >= TOPBIT) die ("too many entries in modset") ;
	modsetSummary (ms, outFile) ;
      }
    else if (ARGMATCH("-f","--seqfile",2))
      { if (ms)
	  { if (rs) readsetDestroy (rs) ;
	    rs = readsetCreate (ms, 1<<16) ;
	    readsetFileRead (rs, argv[-1]) ;
	    fclose (f) ;
	  }
	else fprintf (stderr, "** need to read a modset before a sequence file\n") ;
      }
    else if (ARGMATCH("-r","--read",2))
      { if (rs) readsetDestroy (rs) ;
	rs = readsetRead (argv[-1]) ;
      }
    else if (ARGMATCH("-w","--write",2))
      readsetWrite (rs, argv[-1]) ;
    else if (ARGMATCH("-S","--stats",1))
      readsetStats (rs) ;
    else if (ARGMATCH("-o1","--overlaps1",2))
      { Array a = findOverlaps (rs, arrp(rs->reads,(U32)atoi(argv[-1]), Read), 2) ;
	arrayDestroy (a) ;
      }
    else if (ARGMATCH("-o2","--overlaps2",2))
      { int d = atoi(argv[-1]) ;
	U32 ix ;
	Array a ;
	for (ix = d ; ix < arrayMax(rs->reads) ; ix += d)
	  { a = findOverlaps (rs, arrp(rs->reads,ix,Read), 1) ; arrayDestroy (a) ; }
      }
    else if (ARGMATCH("-o3","--overlap",3))
      printOverlap (rs, (U32)atoi(argv[-2]), (U32)atoi(argv[-1])) ;
    else if (ARGMATCH("-b","--markBadReads",1)) markBadReads (rs) ;
    else if (ARGMATCH("-c","--markContained",1)) markContained (rs) ;
    else if (ARGMATCH("-a1","--assemble1",2)) assembleFromRead (rs, (U32)atoi(argv[-1])) ;
    else if (ARGMATCH("-a2","--assemble2",3))
      assembleFromMod (rs, (U32)atoi(argv[-2]), atoi(argv[-1])) ;
    else if (ARGMATCH("-u","--cluster",1)) cluster (rs) ;
    else if (ARGMATCH("-C","--cleanmods",1)) cleanMods (rs) ;
    else if (ARGMATCH("-T","--testmods",3)) testMods (rs, atoi(argv[-2]), atoi(argv[-1])) ;
    else if (ARGMATCH("-R","--ref",2)) refFlag (rs, argv[-1]) ;
    else if (ARGMATCH("-rb","--resetbits",2)) resetBits (rs, atoi(argv[-1])) ;
    else if (ARGMATCH("-P","--readProperties",1)) readProperties (rs) ;
    else die ("unkown command %s - run without arguments for usage", *argv) ;

    timeUpdate (outFile) ;
  }

  fprintf (outFile, "total resources used: ") ; timeTotal (outFile) ;
  if (outFile != stdout) { printf ("total resources used: ") ; timeTotal (stdout) ; }
}

/************* end of file ************/
