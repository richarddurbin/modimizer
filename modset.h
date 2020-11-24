/*  File: modset.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2018
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Aug 13 19:15 2020 (rd109)
 * Created: Tue Nov  6 17:31:35 2018 (rd109)
 *-------------------------------------------------------------------
 */

#include "utils.h"
#include "seqhash.h"

/* object to hold sets of modimizers */
typedef struct {
  Seqhash *hasher ;
  int tableBits ;		/* max 34 so size < 2^32 so index is 32bit */
  U32 size ;			/* size of value array, and other arrays over this set */
  U64 tableSize ; 		/* = 1 << tableBits */
  U64 tableMask ;		/* = tableSize - 1 */
  U32 *index ;			/* this is the primary table - size tableSize */
  U64 *value ;			/* the kmer values */
  U16 *depth ;			/* depth at each index */
  U8  *info ;			/* bits for various things */
  U32 max ;			/* number of entries in the set - must be less than size */
} Modset ;

Modset *modsetCreate (Seqhash *sh, int bits, U32 size) ;
void modsetDestroy (Modset *ms) ; 
void modsetWrite (Modset *ms, FILE *f) ;
Modset *modsetRead (FILE *f) ;

/* this is the key low level function, both to insert new hashes and find existing ones */
U32 modsetIndexFind (Modset *ms, U64 kmer, int isAdd) ;

/* the following act on the whole set */
void modsetSummary (Modset *ms, FILE *f) ;
bool modsetPack (Modset *ms)	; /* reduce size to max+1 and compress value; TRUE if changes */
void modsetDepthPrune (Modset *ms, int min, int max) ;
bool modsetMerge (Modset *ms1, Modset *ms2) ;

/* info fields */
/* bits 1 and 2 for copy number in {0,1,2,M} with 0 for errors */
/* bit 3 for minor variants, less than half depth of a neighbour in at least one read */
/* bit 4 for repeats within a read */
/* bit 5 for internal within a read - both neighbours within w */
#define MS_MINOR 4
#define MS_REPEAT 8
#define MS_INTERNAL 0x10
#define MS_RDNA 0x20
static inline void msSetCopy0 (Modset *ms, U32 i) { ms->info[i] &= 0xfc ; }
static inline void msSetCopy1 (Modset *ms, U32 i) { ms->info[i] = (ms->info[i] & 0xfc) | 1 ; }
static inline void msSetCopy2 (Modset *ms, U32 i) { ms->info[i] = (ms->info[i] & 0xfc) | 2 ; }
static inline void msSetCopyM (Modset *ms, U32 i) { ms->info[i] |= 3 ; }
static inline void msSetMinor (Modset *ms, U32 i) { ms->info[i] |= MS_MINOR ; }
static inline void msSetRepeat (Modset *ms, U32 i) { ms->info[i] |= MS_REPEAT ; }
static inline void msSetInternal (Modset *ms, U32 i) { ms->info[i] |= MS_INTERNAL ; }
static inline void msSetRDNA (Modset *ms, U32 i) { ms->info[i] |= MS_RDNA ; }
static inline bool msIsCopy0  (Modset *ms, U32 i) { return ((ms->info[i] & 3) == 0) ; }
static inline bool msIsCopy1  (Modset *ms, U32 i) { return ((ms->info[i] & 3) == 1) ; }
static inline bool msIsCopy2  (Modset *ms, U32 i) { return ((ms->info[i] & 3) == 2) ; }
static inline bool msIsCopyM  (Modset *ms, U32 i) { return ((ms->info[i] & 3) == 3) ; }
static inline int  msCopy (Modset *ms, U32 i) { return (ms->info[i] & 3) ; }
static inline bool msIsMinor  (Modset *ms, U32 i) { return (ms->info[i] & MS_MINOR) ; }
static inline bool msIsRepeat  (Modset *ms, U32 i) { return (ms->info[i] & MS_REPEAT) ; }
static inline bool msIsInternal  (Modset *ms, U32 i) { return (ms->info[i] & MS_INTERNAL) ; }
static inline bool msIsRDNA  (Modset *ms, U32 i) { return (ms->info[i] & MS_RDNA) ; }

/*************************/
