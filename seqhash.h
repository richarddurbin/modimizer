/*  File: seqhash.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2018
 *-------------------------------------------------------------------
 * Description: header file for seqhash package - minimizers and moshers
 * Exported functions: see below
 * HISTORY:
 * Last edited: Aug 13 22:13 2020 (rd109)
 * Created: Mon Mar  5 08:43:45 2018 (rd)
 *-------------------------------------------------------------------
 */

#include "utils.h"

typedef struct {
  int seed ;			/* seed */
  int k ;			/* kmer */
  int w ;			/* window */
  U64 mask ;			/* 2*k bits */
  int shift1, shift2 ;
  U64 factor1, factor2 ;
  U64 patternRC[4] ;		/* one per base */
} Seqhash ;

typedef struct {
  Seqhash *sh ;
  char *s, *sEnd ;     		/* sequence currently being hashed, end marker */
  U64 h, hRC ;			/* current hash values */
  U64 *hashBuf ;		/* buffer of length w holding hashes for current window */
  bool *fBuf ;			/* buffer of length w holding isForward for current window */
  int base ;			/* start of buf in sequence */
  int iStart, iMin ;		/* position in buf of start of current window, next min */
  bool isDone ;
} SeqhashRCiterator ;

Seqhash *seqhashCreate (int k, int w, int seed) ;
static void seqhashDestroy (Seqhash *sh) { free (sh) ; }

void seqhashWrite (Seqhash *sh, FILE *f) ;
Seqhash *seqhashRead (FILE *f) ;
void seqhashReport (Seqhash *sh, FILE *f) ;

// iterator to extract minimizers from a sequence
// NB sequence must continue to exist through the life of the iterator
SeqhashRCiterator *minimizerRCiterator (Seqhash *sh, char *s, int len) ;
bool minimizerRCnext (SeqhashRCiterator *si, U64 *u, int *pos, bool *isF) ; /* return u,pos,isF */

// modimizer extracts hashes that are divisible by m->w
// this is faster and more robust to errors - same mean density without evenness guarantees
SeqhashRCiterator *modRCiterator (Seqhash *sh, char *s, int len) ;
bool modRCnext (SeqhashRCiterator *si, U64 *kmer, int *pos, bool *isF) ;
/* returns any/all of kmer, pos, isF - get hash from seqhash(sh,kmer) */

static void seqhashRCiteratorDestroy (SeqhashRCiterator *si)
{ free (si->hashBuf) ; free (si->fBuf) ; free (si) ; }

// utilities
static inline U64 seqhash (Seqhash *sh, U64 k) { return ((k * sh->factor1) >> sh->shift1) ; }
char *seqString (U64 kmer, int len)  ;
static inline char* seqhashString (Seqhash *sh, U64 k) { return seqString (k, sh->k) ; }

/******* end of file ********/
