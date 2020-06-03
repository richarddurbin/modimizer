/*  File: modasm.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2018
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Feb 17 10:13 2019 (rd109)
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
BOOL isVerbose = FALSE ;

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
  U8 otherFlags ;		/* for future use */
  U16 pad1 ;			/* for future use - to pad out flags to 32 bits */
  int nMiss ;			/* number of mods in sequence that miss ms */
  int contained ;		/* index in rs->reads of read that contains this one */
  int nCopy[4] ;		/* number of copy 0, 1, 2, M hits - could create on demand */
  U32 pad2[4] ;			/* for future use */
} Read ;

typedef struct {		/* data structure for long read set */
  Modset *ms ;
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
  for (i = 0 ; i < n ; ++i, ++r) if (r->nHit) { free(r->hit) ; free(r->dx) ; }
  arrayDestroy (rs->reads) ;
  if (rs->inv) { free (rs->inv) ; free (rs->invSpace) ; }
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
  int i, n ; Read *r = arrp(rs->reads,0,Read) ;
  for (i = 0 ; i < arrayMax(rs->reads) ; ++i, ++r)
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
  int i, n ; Read *r = arrp(rs->reads,0,Read) ;
  for (i = 0 ; i < arrayMax(rs->reads) ; ++i, ++r)
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
  SeqIO *si = seqIOopenRead (filename, dna2indexConv, FALSE) ;
  while (seqIOread (si))
    { Read *read = arrayp(rs->reads, arrayMax(rs->reads), Read) ;
      read->len = si->seqLen ;
      SeqhashRCiterator *mi = modRCiterator (rs->ms->hasher, sqioSeq(si), si->seqLen) ;
      hitsA = arrayReCreate (hitsA, 1024, U32) ;
      dxA = arrayReCreate (dxA, 1024, U16) ;
      U64 hash ; int lastPos = 0, pos ; BOOL isForward ;
      while (modRCnext (mi, &hash, &pos, &isForward))
	{ U32 index = modsetIndexFind (rs->ms, hash, FALSE) ;
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
  rs->inv = new0 (ms->max+1, U32*) ;
  rs->invSpace = new (rs->totHit, U32) ;
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
  U32 iy ;			/* index of overlap read in rs->reads */
  int nHit ;			/* number of shared mod hits */
  BOOL isPlus ;			/* relative direction */
  BOOL isBad ;
  BOOL isContained ;		/* x is contained in y */
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
      o->isBad = FALSE ;
      hy = y->hit ; /* reset hy for loop to come */
      U16 *dy = y->dx ; double yPos = *dy ;
      if (nPlus && !nMinus)
	{ o->isPlus = TRUE ;
	  int last = 0, lastDiff ;
	  for (j = 0 ; j < y->nHit ; ++j, ++hy, yPos += *++dy)
	    if ((ihx = hmap[*hy & TOPMASK]))
	      { lastDiff = xPos[ihx] - yPos ;
		if (!last && lastDiff < 0) o->isContained = TRUE ; /* x starts in y */
		if (ihx < last) { o->isBad = TRUE ; --nPlus ; } last = ihx ;
	      }
	  if (o->isContained && x->len - lastDiff > y->len) o->isContained = FALSE ;
	}
      else if (nMinus && !nPlus)
	{ o->isPlus = FALSE ;
	  int last = x->nHit, lastDiff ;
	  for (j = 0 ; j < y->nHit ; ++j, ++hy, yPos += *++dy)
	    if ((ihx = hmap[*hy & TOPMASK]))
	      { lastDiff = x->len - xPos[ihx] - yPos ;
		if (!last && lastDiff < 0) o->isContained = TRUE ; /* x starts in y */
		if (ihx > last) { o->isBad = TRUE ; --nMinus ; } last = ihx ;
	      }
	  if (o->isContained && x->len - lastDiff > y->len) o->isContained = FALSE ;
	}
      if (nPlus && nMinus) o->isBad = TRUE ;
      if (o->isBad) ++nBad ; else ++nGood ;
      
      if (reportLevel > 1)
	{ fprintf (outFile, "RH\t%u\tlen %d\t%s\tnPlus %d\tnMinus %d\t",
		   o->iy, y->len, o->isBad ? "BAD" : "GOOD", nPlus, nMinus) ;
	  fprintf (outFile, "%s\n", o->isContained ? "CONTAINED" : "OVERLAP") ;
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
  U32 *hx = x->hit ;
  for (j = 0 ; j < x->nHit ; ++j)
    { U32 hxx = *hx & TOPMASK ;
      xPos += x->dx[j] ;
      if (msIsCopy1 (rs->ms, hxx))
	{ U32 *hy = y->hit ;
	  int yPos = 0 ;
	  for (k = 0 ; k < y->nHit ; ++k, ++hy)
	    { U32 hyy = *hy & TOPMASK ;
	      yPos += y->dx[0] ;
	      if (hxx == hyy)
		{ fprintf (outFile, "RO\t%8x %5d %c\t",
			   hxx, rs->ms->depth[hxx], (*hx&TOPBIT) == (*hy&TOPBIT) ? '+' : '-') ;
		  fprintf (outFile, "%u %u %c\t",
			   ix, xPos, (*hx & TOPBIT) ? 'F' : 'R') ;
		  fprintf (outFile, "%u %u %c\n",
			   iy, yPos, (*hy & TOPBIT) ? 'F' : 'R') ;
		}
	    }
	}
    }
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
  Read *x, *x0 = arrp(rs->reads,0,Read) ;

  for (ix = 0, x = x0 ; ix < arrayMax(rs->reads) ; ++ix, ++x)
    x->bad = 0 ; /* initialise by clearing all bad flags */

  int *badList = new0 (arrayMax(rs->reads)*10, int) ;
  int *nBad = new0 (arrayMax(rs->reads), int), *lBad = new0 (arrayMax(rs->reads), int) ;
  
  /* pass through findOverlaps and record how many bad overlaps in bady and 
     the bad overlaps themselves to reads possibly still good in badList
   */
  for (ix = 0, x = x0 ; ix < arrayMax(rs->reads) ; ++ix, ++x)
    { Array olap = findOverlaps (rs, x, 0) ;
      for (i = 0 ; i < arrayMax(olap) ; ++i)
	if (arrp(olap,i,Overlap)->isBad)
	  { int iy = arrp(olap,i,Overlap)->iy ;
	    ++nBad[iy] ;
	    if (nBad[iy] < 10 && lBad[ix] < 10)
	      badList[10*ix + lBad[ix]++] = iy ;
	  }
      arrayDestroy (olap) ;
    }

  /* first pass, set badOrder10 for all 10 or more - these are clearly bad */
  int N = 0 ;
  for (ix = 0, x = x0 ; ix < arrayMax(rs->reads) ; ++ix, ++x)
    if (nBad[ix] >= 10) { x->badOrder10 = 1 ; ++N ; lBad[ix] = 0 ; }
  printf ("MB  %d with >=10 bad overlaps\n", N) ;

  for (ix = 0, x = x0 ; ix < arrayMax(rs->reads) ; ++ix, ++x) /* remove bad from lists */
    for (i = lBad[ix] ; i-- ; )
      if (arrp(rs->reads, badList[10*ix+i], Read)->bad)
	badList[10*ix+i] = badList[10*ix + --lBad[ix]] ;
  
  /* second pass, set badOrder1 for all 2 or more - removes some more singleton errors */
  N = 0 ;
  for (ix = 0, x = x0 ; ix < arrayMax(rs->reads) ; ++ix, ++x)
    if (lBad[ix] >= 2) { x->badOrder1 = 1 ; ++N ; lBad[ix] = 0 ; }
  printf ("MB  %d with multiple bad overlaps\n", N) ;

  for (ix = 0, x = x0 ; ix < arrayMax(rs->reads) ; ++ix, ++x) /* remove bad from lists */
    for (i = lBad[ix] ; i-- ; )
      if (arrp(rs->reads, badList[10*ix+i], Read)->bad)
	badList[10*ix+i] = badList[10*ix + --lBad[ix]] ;

  /* final pass, mark anything left with any bad overlaps as bad */
  N = 0 ;
  for (ix = 0, x = x0 ; ix < arrayMax(rs->reads) ; ++ix, ++x)
    if (lBad[ix] > 0) { x->badOrder1 = 1 ; ++N ; lBad[ix] = 0 ; }
  printf ("MB  %d with single bad overlaps\n", N) ;

  free (badList) ; free (nBad) ; free (lBad) ;
}

int badOverlaps (Readset *rs, Read *x)
{
  int i, nBad = 0 ;
  Array olap = findOverlaps (rs, x, 0) ;
  for (i = 0 ; i < arrayMax(olap) ; ++i) if (arrp(olap,i,Overlap)->isBad) ++nBad ;
  arrayDestroy (olap) ;
  return nBad ;
}

void markBadReadsOld (Readset *rs)
{
  int ix ;
  Read *x, *x0 = arrp(rs->reads,0,Read) ;

  for (ix = 0, x = x0 ; ix < arrayMax(rs->reads) ; ++ix, ++x)
    x->bad = 0 ; /* initialise by clearing all bad flags */
  
  /* first pass through findOverlaps - set badOrder10 */
  int nBad = 0 ;
  for (ix = 0, x = x0 ; ix < arrayMax(rs->reads) ; ++ix, ++x)
    if (badOverlaps (rs, x) >= 10) { x->badOrder10 = 1 ; ++nBad ; }
  printf ("MB  %d with >=10 bad overlaps\n", nBad) ;


  nBad = 0 ;
  for (ix = 0, x = x0 ; ix < arrayMax(rs->reads) ; ++ix, ++x)
    if (badOverlaps (rs, x) > 1 && !x->badOrder10) { x->badOrder1 = 1 ; ++nBad ; }
  printf ("MB  %d with multiple bad overlaps\n", nBad) ;

  /* third pass, set badOrder1 for anything bad  - must be simple pairs - might be too harsh */
  nBad = 0 ;
  for (ix = 0, x = x0 ; ix < arrayMax(rs->reads) ; ++ix, ++x)
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
  Read *x, *x0 = arrp(rs->reads,0,Read) ;
  int nContained = 0, nNotContained = 0 ;
  U64 totLen = 0 ;
  Array overlaps = 0 ;
  
  for (ix = 0, x = x0 ; ix < arrayMax(rs->reads) ; ++ix, ++x)
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
    fprintf (stderr, "COMMAND %s", *argv) ;
    for (i = 1 ; i < argc && *argv[i] != '-' ; ++i) fprintf (stderr, " %s", argv[i]) ;
    fputc ('\n', stderr) ;
    
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
      { if (!(f = fopen (argv[-1], "r"))) die ("failed to open mod file %s", argv[-1]) ;
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
    else die ("unkown command %s - run without arguments for usage", *argv) ;

    timeUpdate (outFile) ;
  }

  fprintf (outFile, "total resources used: ") ; timeTotal (outFile) ;
  if (outFile != stdout) { printf ("total resources used: ") ; timeTotal (stdout) ; }
}

/************* end of file ************/
