/*  File: seqio.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2018
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Jun 20 13:03 2020 (rd109)
 * Created: Sat Nov 10 08:51:49 2018 (rd109)
 *-------------------------------------------------------------------
 */

#ifndef SEQIO_DEFINED
#define SEQIO_DEFINED

#include "utils.h"
#include <zlib.h>

typedef enum { UNKNOWN, FASTA, FASTQ, BINARY, ONE, BAM } SeqIOtype ;
extern char* seqIOtypeName[] ; // = { "unknown", "fasta", "fastq", "binary", "onecode", "bam" } ;

typedef struct {
  SeqIOtype type ;
  bool isWrite ;
  U64 nSeq, totIdLen, totDescLen, totSeqLen, maxIdLen, maxDescLen, maxSeqLen ;
  U64 idLen, descLen, seqLen ;
  U64 idStart, descStart, seqStart, qualStart ;
  bool isQual ;			/* if set then convert qualities by subtracting 33 (FASTQ) */
  int qualThresh ;		/* used for binary representation of qualities */
  /* below here private */
  U64 bufSize ;
  U64 nb ;			/* nb is how many characters left to read in the buffer */
  U64 line, recStart ;		/* recStart is the offset for the current record */
  int fd ;			/* file descriptor, if gzf is not set */
  gzFile gzf ;
  char *buf, *b ;		/* b is current pointer in buf */
  int *convert ;
  char *seqBuf, *qualBuf ;	/* used in modes BINARY, VGP, BAM */
  char unpackConvert[4] ;	/* for unpacking */
  U32 seqExpand[256] ;		/* lookup for unpacking sequence */
  U64 qualExpand[256] ;		/* lookup for unpacking qual */
  void *handle;			/* used for VGP, BAM */
} SeqIO ;

/* Reads FASTA or FASTQ, gzipped or not. */
/* Philosophy here is to read blocks of 8Mb and provide direct access into the buffer. */
/* So the user does not own the pointers. */
/* Add 0 terminators to ids.  Convert sequences in place if convert != 0, and quals if isQual. */
/* Potential for future storage of (packed?) sequences with counts up front. */ 

SeqIO *seqIOopenRead (char *filename, int* convert, bool isQual) ; /* can use "-" for stdin */
bool seqIOread (SeqIO *si) ;
#define sqioId(si)   ((si)->buf+(si)->idStart)
#define sqioDesc(si) ((si)->buf+(si)->descStart)
#define sqioSeq(si)  ((si)->type >= BINARY ? (si)->seqBuf : (si)->buf+(si)->seqStart)
#define sqioQual(si) ((si)->type >= BINARY ? (si)->qualBuf : (si)->buf+(si)->qualStart)

SeqIO *seqIOopenWrite (char *filename, SeqIOtype type, int* convert, int qualThresh) ;
void seqIOwrite (SeqIO *si, char *id, char *desc, U64 seqLen, char *seq, char *qual) ;
void seqIOflush (SeqIO *si) ;	/* NB writes are buffered, so need this to ensure in file */

void seqIOclose (SeqIO *si) ;	/* will flush file opened for writing */

/* used for binary file writing/reading, but may be useful elsewhere */
U64  sqioSeqPack (char *s, U8 *u, U64 len, int *convert) ; /* compress s into (len+3)/4 u */
void sqioSeqUnpack (U8 *u, char *s, U64 len, SeqIO *si) ; /* uncompress (len+3)/4 u into s */
U64  sqioQualPack (char *q, U8 *u, U64 len, int thresh) ; /* compress q into (len+7)/8 u */
void sqioQualUnpack (U8 *u, char *q, U64 len, SeqIO *si) ; /* uncompress (len+7)/8 u into q */

extern int dna2textConv[] ;
extern int dna2textAmbigConv[] ;
extern int dna2textAmbig2NConv[] ;
extern int dna2indexConv[] ;
extern int dna2index4Conv[] ;
extern int dna2binaryConv[] ;
extern int dna2binaryAmbigConv[] ;
static const char index2char[] = "acgtn" ;
static const char binary2char[] = "-ACMGRSVTWYHKDBN" ;
extern int aa2textConv[] ;
extern int aa2indexConv[] ;
static const char index2aa[] = "ACDEFGHIKLMNPQRSTVWYX*" ;
extern int noConv[] ;

#endif	/* SEQIO_DEFINED */

/******************************************************************/
