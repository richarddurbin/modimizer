/*  File: seqhash.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2018
 *-------------------------------------------------------------------
 * Description: seqhash package - uses random bit patterns for bases and rotate/XOR
 	compile with -DTEST to test, and -DDEBUG to debug, -DNDEBUG turns off asserts
	see test main() at end for standard usage pattern
 * Exported functions: see seqhash.h
 * HISTORY:
 * Last edited: Aug 13 19:21 2020 (rd109)
 * Created: Sat Feb 24 19:20:18 2018 (rd)
 *-------------------------------------------------------------------
 */

#include "seqhash.h"

static inline U64 rotateLeft (U64 x) { return (x << 2) | (x >> 62) ; }
static inline U64 rotateRight (U64 x) { return (x << 62) | (x >> 2) ; }

Seqhash *seqhashCreate (int k, int w, int seed)
{
  assert (sizeof (U64) == 8) ;
  Seqhash *sh = new0 (1, Seqhash) ;
  sh->k = k ; if (k < 1 || k >= 32) die ("seqhash k %d must be between 1 and 32\n", k) ;
  sh->w = w ; if (w < 1) die ("seqhash w %d must be positive\n", w) ;
  sh->seed = seed ;
  sh->mask = ((U64)1 << (2*k)) - 1 ;
  int i ;
  
  srandom (seed) ;
  sh->factor1 = (random() << 32) | random() | 0x01 ;
  sh->shift1 = 64 - 2*k ;
  sh->factor2 = (random() << 32) | random() | 0x01 ;
  sh->shift2 = 2*k ;
  for (i = 0 ; i < 4 ; ++i) { sh->patternRC[i] = (3-i) ; sh->patternRC[i] <<= 2*(k-1) ; }
  return sh ;
}

#include <stdio.h>

void seqhashWrite (Seqhash *sh, FILE *f)
{ if (fwrite ("SQHSHv2",8,1,f) != 1) die ("failed to write seqhash header") ;
  if (fwrite (sh,sizeof(Seqhash),1,f) != 1) die ("failed to write seqhash") ;
}

Seqhash *seqhashRead (FILE *f)
{ Seqhash *sh = new (1, Seqhash) ;
  char name[8] ;
  if (fread (name,8,1,f) != 1) die ("failed to read seqhash header") ;
  if (strcmp (name, "SQHSHv2")) die ("seqhash read mismatch") ;
  if (fread (sh,sizeof(Seqhash),1,f) != 1) die ("failed to read seqhash") ;
  return sh ;
}

void seqhashReport (Seqhash *sh, FILE *f)
{ fprintf (f, "SH k %d  w/m %d  s %d\n", sh->k, sh->w, sh->seed) ; }

/************** basic hash functions *************/

static inline U64 hashRC (SeqhashRCiterator *si, bool *isForward)
{ U64 hashF = seqhash (si->sh, si->h) ;
  U64 hashR = seqhash (si->sh, si->hRC) ;
#ifdef DEBUG
  printf ("hashRC: h %lx hRC %lx hashF %lx hashR %lx\n", si->h, si->hRC, hashF, hashR) ;
#endif
  if (hashF < hashR) { *isForward = true ; return hashF ; }
  else { *isForward = false ; return hashR ; }
}

static inline U64 advanceHashRC (SeqhashRCiterator *si, bool *isForward)
{ Seqhash *sh = si->sh ;
  if (si->s < si->sEnd)
    { si->h = ((si->h << 2) & sh->mask) | *(si->s) ;
      si->hRC = (si->hRC >> 2) | sh->patternRC[*(si->s)] ;
      return hashRC (si, isForward) ;
    }
  else
    return U64MAX ;
}

/************ iterator to run across a sequence ***********/

SeqhashRCiterator *minimizerRCiterator (Seqhash *sh, char *s, int len)
{
  assert (s && len >= 0) ;
  SeqhashRCiterator *si = (SeqhashRCiterator*) mycalloc (sizeof(SeqhashRCiterator), 1) ;
  si->sh = sh ;
  si->s = s ; si->sEnd = s + len ;
  si->hashBuf = new0 (sh->w, U64) ;
  si->fBuf = new0 (sh->w, bool) ;
  if (len < sh->k) { si->isDone = true ; return si ; } /* edge case */

  int i ;			/* preinitialise the hashes for the first kmer */
  for (i = 0 ; i < sh->k ; ++i, ++si->s)
    { si->h = (si->h << 2) | *si->s ;
      si->hRC = (si->hRC >> 2) | sh->patternRC[*(si->s)] ;
    }
  
  /* store first w hashes in hashBuf and set ->iMin */
  U64 min = hashRC (si, si->fBuf) ;
  si->iMin = 0 ;
  for (i = 1 ; i < sh->w ; ++i, ++si->s)
    { si->hashBuf[i] = advanceHashRC (si, &si->fBuf[i]) ;
      if (si->hashBuf[i] < min) { min = si->hashBuf[i] ; si->iMin = i ; }
    }

  return si ;
}

bool minimizerRCnext (SeqhashRCiterator *si, U64 *u, int *pos, bool *isF) /* returns u,pos,isF */
{
  if (si->isDone) return false ; /* we are done */

#ifdef DEBUG
  printf ("base %d, iStart %d, iMin %d\n", si->base, si->iStart, si->iMin) ;
  int j ; for (j = 0 ; j < si->sh->w ; ++j) printf ("  %x", si->hashBuf[j]) ;
  printf ("\n") ;
#endif

  assert (u) ;
  *u = si->hashBuf[si->iMin] ;
  if (pos) { *pos = si->base + si->iMin ; if (si->iMin < si->iStart) *pos += si->sh->w ; }
  if (isF) *isF = si->fBuf[si->iMin] ;
  if (si->s >= si->sEnd) { si->isDone = true ; return true ; }

  int i ;	    		/* next update hashBuf splitting into two cases */
  U64 min = *u ;    /* save this here for end case - see below */
  if (si->iMin >= si->iStart)
    for (i = si->iStart ; i <= si->iMin ; ++i, ++si->s)
      si->hashBuf[i] = advanceHashRC (si, &si->fBuf[i]) ;
  else
    { for (i = si->iStart ; i < si->sh->w ; ++i, ++si->s)
	si->hashBuf[i] = advanceHashRC (si, &si->fBuf[i]) ;
      si->base += si->sh->w ;
      for (i = 0 ; i <= si->iMin ; ++i, ++si->s)
	si->hashBuf[i] = advanceHashRC (si, &si->fBuf[i]) ;
    }
  si->iStart = si->iMin + 1 ;
  if (si->iStart == si->sh->w) { si->iStart = 0 ; si->base += si->sh->w ; }

  /* finally find new min to set up for next call */
  if (si->hashBuf[si->iMin] != U64MAX) /* there was a full new window */
    min = U64MAX ;
  else				/* otherwise, keep the last min */
    si->iMin = -1 ;
  for (i = 0 ; i < si->sh->w ; ++i)
    if (si->hashBuf[i] < min) { min = si->hashBuf[i] ; si->iMin = i ; }
  if (si->iMin == -1)		/* our old min was not beaten - we are done */
    si->isDone = true ;
  
  return true ;
}

SeqhashRCiterator *modRCiterator (Seqhash *sh, char *s, int len)
{
  assert (s && len >= 0) ;
  SeqhashRCiterator *si = (SeqhashRCiterator*) mycalloc (sizeof(SeqhashRCiterator), 1) ;
  si->sh = sh ;
  si->s = s ; si->sEnd = s + len ;
  si->hashBuf = new0 (sh->w, U64) ;
  si->fBuf = new0 (sh->w, bool) ;
  if (len < sh->k) { si->isDone = true ; return si ; } /* edge case */

  int i ;			/* preinitialise the hashes for the first kmer */
  for (i = 0 ; i < sh->k ; ++i, ++si->s)
    { si->h = (si->h << 2) | *si->s ;
      si->hRC = (si->hRC >> 2) | sh->patternRC[*(si->s)] ;
    }
  
  U64 hash = hashRC(si, si->fBuf) ;
  while ((hash % sh->w) && si->s < si->sEnd)
    { hash = advanceHashRC (si, si->fBuf) ; ++si->iMin ; ++si->s ; }
  if (!(hash % sh->w)) *si->hashBuf = hash ;
  else si->isDone = true ;

  return si ;
}

bool modRCnext (SeqhashRCiterator *si, U64 *kmer, int *pos, bool *isF)
{
  if (si->isDone) return false ; /* we are done */

  if (kmer) { if (*si->fBuf) *kmer = si->h ; else *kmer = si->hRC ; }
  if (pos) *pos = si->iMin ; if (isF) *isF = *si->fBuf ;

  if (si->s >= si->sEnd) { si->isDone = true ; return true ; }
    
  U64 u = advanceHashRC (si, si->fBuf) ; ++si->iMin ; ++si->s ;
  int w = si->sh->w ;
  while ((u % w) && si->s < si->sEnd)
    { u = advanceHashRC (si, si->fBuf) ; ++si->iMin ; ++si->s ; }
  if (!(u % w)) *si->hashBuf = u ;
  else si->isDone = true ;

  return true ;
}

char *seqString (U64 kmer, int len)
{
  static char trans[4] = { 'a', 'c', 'g', 't' } ;
  static char buf[33] ;
  assert (len <= 32) ;
  buf[len] = 0 ;
  while (len--) { buf[len] = trans[kmer & 0x3] ; kmer >>= 2 ; }
  return buf ;
}

/************** short test program, illustrating standard usage *************/

#ifdef TEST

#include "seqio.h"

int main (int argc, char *argv[])
{
  char *seq, *id ;
  int len ;
  U64 u, k ; int pos ; bool isF ;

  SeqIO *sio = seqIOopen ("-", dna2indexConv, 0) ;
  
  Seqhash *sh = seqhashCreate (16, 32, 0) ;
  while (seqIOread (sio))
    { printf ("\nread sequence %s length %d\n", sqioId(sio), sio->seqLen) ;
      SeqhashRCiterator *si = modRCiterator (sh, sqioSeq(sio), sio->seqLen) ;
      while (modRCnext (si, &u, &pos, &isF, &k))
	printf ("\t%08lx\t%s\t%d\t%c\n", u, seqString (k, 16), pos, isF?'F':'R') ;
      seqhashRCiteratorDestroy (si) ;
    }
  seqIOclose (sio) ;
}

#endif

/**************** end of file ****************/
