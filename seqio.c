/*  File: seqio.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2018
 *-------------------------------------------------------------------
 * Description: buffered package to read arbitrary sequence files - much faster than readseq
 * Exported functions:
 * HISTORY:
 * Last edited: Dec  7 23:39 2020 (rd109)
 * Created: Fri Nov  9 00:21:21 2018 (rd109)
 *-------------------------------------------------------------------
 */

#include "seqio.h"
#include <fcntl.h>
#include <unistd.h>

#ifdef ONEIO
#include "ONElib.h"
#endif

#ifdef BAMIO
bool bamFileOpenRead (char* filename, SeqIO *si) ;
bool bamRead (SeqIO *si) ;
void bamFileClose (SeqIO *si) ;
#endif

// global
char* seqIOtypeName[] = { "unknown", "fasta", "fastq", "binary", "onecode", "bam" } ;

SeqIO *seqIOopenRead (char *filename, int* convert, bool isQual)
{
  SeqIO *si = new0 (1, SeqIO) ;
  if (!strcmp (filename, "-")) si->gzf = gzdopen (fileno (stdin), "r") ;
  else si->gzf = gzopen (filename, "r") ;
  if (!si->gzf) { free(si) ; return 0 ; }
  si->bufSize = 1<<24 ;
  si->b = si->buf = new (si->bufSize, char) ;
  si->convert = convert ;
  si->isQual = isQual ;
  si->nb = gzread (si->gzf, si->buf, si->bufSize) ;
  if (!si->nb)
    { fprintf (stderr, "sequence file %s unreadable or empty\n", filename) ;
      seqIOclose (si) ;
      return 0 ;
    }
  si->line = 1 ;
  if (*si->buf == '>')
    { si->type = FASTA ; si->isQual = false ;
      if (!si->convert) si->convert = dna2textAmbigConv ; /* default: need to remove whitespace */
    }
  else if (*si->buf == '@')
    {
#ifdef BAMIO // problem: SAM file headers can start with this...
      if (si->buf[3] == '\t' &&
	  ( (si->buf[1] == 'H' && si->buf[2] == 'D') ||
	    (si->buf[1] == 'S' && si->buf[2] == 'Q') ||
	    (si->buf[1] == 'R' && si->buf[2] == 'G') ||
	    (si->buf[1] == 'P' && si->buf[2] == 'G') ||
	    (si->buf[1] == 'C' && si->buf[2] == 'O'))) // then almost certainly a SAM file
	{ gzclose (si->gzf) ; si->gzf = 0 ;
	  if (!bamFileOpenRead (filename, si))
	    { fprintf (stderr, "failed to open file %s as SAM/BAM/CRAM\n", filename) ;
	      seqIOclose (si) ;
	      return 0 ;
	    }
	  si->type = BAM ; // important that this is after successful open
	}
      else
	si->type = FASTQ ;
#else
      si->type = FASTQ ;
#endif
    }
  else if (*si->buf == 'b')
    { si->type = BINARY ;
      if (!si->convert) si->convert = dna2textConv ;
      if (si->nb <= 64) die ("binary file too short\n") ;
      si->qualThresh = si->buf[1] ; si->b += 8 ;
      si->nSeq = *(U64*)si->b ; si->b += 8 ;
      si->totIdLen = *(U64*)si->b ; si->b += 8 ;
      si->totDescLen = *(U64*)si->b ; si->b += 8 ;
      si->totSeqLen = *(U64*)si->b ; si->b += 8 ;
      si->maxIdLen = *(U64*)si->b ; si->b += 8 ;
      si->maxDescLen = *(U64*)si->b ; si->b += 8 ;
      si->maxSeqLen = *(U64*)si->b ; si->b += 8 ;
      si->nb -= 64 ;
      si->unpackConvert[0] = si->convert['A'] ; 
      si->unpackConvert[1] = si->convert['C'] ; 
      si->unpackConvert[2] = si->convert['G'] ; 
      si->unpackConvert[3] = si->convert['T'] ; 
      { int i = 256 ; U8 u ;
	for (u = 0 ; i-- ; u++)
	  { sqioSeqUnpack (&u, (char*)&si->seqExpand[u], 4, si) ;
	    if (si->isQual)
	      sqioQualUnpack (&u, (char*)&si->qualExpand[u], 8, si) ;
	  }
      }
      si->seqBuf = new0 (si->maxSeqLen+1, char) ;
      if (si->isQual) si->qualBuf = new0 (si->maxSeqLen+1, char) ;
      { U64 maxBufSize = 3*sizeof(int) + 5 + si->maxIdLen + si->maxDescLen + si->maxSeqLen / 4 ;
	if (si->qualThresh) maxBufSize += (si->maxSeqLen / 8) ; /* 5 = 1 + 1 + 3 for pad */
	if (maxBufSize > si->nb)
	  { maxBufSize = ((maxBufSize >> 20) + 1) << 20 ; /* so a clean number of megabytes */
	    char *newBuf = new (maxBufSize, char) ; memcpy (newBuf, si->b, si->nb) ;
	    si->b = si->buf = newBuf ; si->bufSize = maxBufSize ;
	    si->nb += gzread (si->gzf, si->b + si->nb, si->bufSize - si->nb) ;
	  }
      }
    }
#ifdef ONEIO
  else if (*si->buf == '1')
    { gzclose (si->gzf) ; si->gzf = 0 ;
      OneFile *vf = oneFileOpenRead (filename, 0, "seq", 1) ;
      if (!vf)
	{ fprintf (stderr, "failed to open ONE seq file %s\n", filename) ;
	  seqIOclose (si) ;
	  return 0 ;
	}
      si->type = ONE ; // important that this is after successful open
      si->handle = vf ;
      if (vf->info['S']->given.count)
	{ si->nSeq = vf->info['S']->given.count ;
	  si->totSeqLen = vf->info['S']->given.total ;
	  si->maxSeqLen = vf->info['S']->given.max ;
	  si->seqBuf = new0 (si->maxSeqLen+1,char) ;
	  si->qualBuf = new0 (si->maxSeqLen+1,char) ;
	}
      while (oneReadLine (vf) && vf->lineType != 'S') ; // move up to first sequence line
      si->seqStart = 0 ;
    }
#endif
#ifdef BAMIO
  else
    { gzclose (si->gzf) ; si->gzf = 0 ;
      if (!bamFileOpenRead (filename, si))
	{ fprintf (stderr, "failed to open file %s as SAM/BAM/CRAM\n", filename) ;
	  seqIOclose (si) ;
	  return 0 ;
	}
      si->type = BAM ; // important that this is after successful open
    }
#else
  else
    { fprintf (stderr, "sequence file %s is unknown type\n", filename) ;
      seqIOclose (si) ;
      return 0 ;
    }
#endif
  return si ;
}

void seqIOclose (SeqIO *si)
{ if (si->isWrite)
    { if (si->type <= BINARY)
	seqIOflush (si) ;
      if (si->type == BINARY)	/* write header */
	{ if (lseek (si->fd, 0, SEEK_SET)) die ("failed to seek to start of binary file") ;
	  si->b = si->buf; 
	  *si->b++ = 'b' ; *si->b++ = si->qualThresh ; si->b += 6 ;
	  *(U64*)si->b = si->nSeq  ; si->b += 8 ;
	  *(U64*)si->b = si->totIdLen ; si->b += 8 ;
	  *(U64*)si->b = si->totDescLen ; si->b += 8 ;
	  *(U64*)si->b = si->totSeqLen ; si->b += 8 ;
	  *(U64*)si->b = si->maxIdLen ; si->b += 8 ;
	  *(U64*)si->b = si->maxDescLen ; si->b += 8 ;
	  *(U64*)si->b = si->maxSeqLen ; si->b += 8 ;
	  seqIOflush (si) ;
	}
#ifdef ONEIO
      if (si->type == ONE)
	oneFileClose ((OneFile*)si->handle) ;
#endif
#ifdef BAMIO
      if (si->type == BAM)
	bamFileClose (si->handle) ;
#endif
    }
  free (si->buf) ;
  if (si->seqBuf) free (si->seqBuf) ;
  if (si->qualBuf) free (si->qualBuf) ;
  if (si->gzf) gzclose (si->gzf) ;
  if (si->fd) close (si->fd) ;
  free (si) ;
}

/********** local routines for seqIOread() ***********/
 
static void bufRefill (SeqIO *si)
{ si->b -= si->recStart ;		/* will be position after move */
  memmove (si->buf, si->buf + si->recStart, si->b - si->buf) ;
  si->idStart -= si->recStart ; si->descStart -= si->recStart ; /* adjust all the offsets */
  si->seqStart -= si->recStart ; si->qualStart -= si->recStart ;
  si->recStart = 0 ;
  si->nb = gzread (si->gzf, si->b, si->buf + si->bufSize - si->b) ;
}

static void bufDouble (SeqIO *si)
{
  char *newbuf = new (si->bufSize*2, char) ;
  memcpy (newbuf, si->buf, si->bufSize) ;
  si->b = newbuf + si->bufSize ; si->nb = si->bufSize ; /* rely on being at end of old buf */
  free (si->buf) ; si->buf = newbuf ;
  si->nb = gzread (si->gzf, si->b, si->bufSize) ;
  si->bufSize *= 2 ;
}

#define bufAdvanceEndRecord(si) \
  { ++si->b ; \
    if (!--si->nb) \
      { if (si->recStart) bufRefill (si) ; \
        else if (si->b == si->buf + si->bufSize) bufDouble (si) ; \
      } \
  }

#define bufAdvanceInRecord(si) \
  { bufAdvanceEndRecord(si) ; \
    if (!si->nb) \
      { fprintf (stderr, "incomplete sequence record line %" PRIu64 "\n", si->line) ; return false ; } \
  } 

static void bufHardRefill (SeqIO *si, U64 n) /* like bufRefill() but for bufConfirmNbytes() */
{					     /* NB buf should be big enough because of header */
  si->b -= si->recStart ;		/* will be position after move */
  memmove (si->buf, si->buf + si->recStart, si->b - si->buf) ;
  si->recStart = 0 ; si->b = si->buf ;
  si->nb += gzread (si->gzf, si->b + si->nb, si->bufSize - si->nb) ;
  if (si->nb < n) die ("incomplete sequence record %" PRIu64 "", si->line) ;
}

#define bufConfirmNbytes(si, n) { if (si->nb < n) bufHardRefill (si, n) ; }

#include <ctype.h>

bool seqIOread (SeqIO *si)
{

#ifdef ONEIO
  if (si->type == ONE)
    { OneFile *vf = (OneFile*) si->handle ;
      if (vf->lineType != 'S') return false ; // at end of file
      si->seqLen = oneLen(vf) ; // otherwise we are at an 'S' line
      if (si->seqLen > si->maxSeqLen)
	{ if (si->maxSeqLen) { free (si->seqBuf) ; if (si->isQual) free (si->qualBuf) ; }
	  si->maxSeqLen = si->seqLen ;
	  si->seqBuf = new0 (si->maxSeqLen+1, char) ;
	  if (si->isQual) si->qualBuf = new0 (si->maxSeqLen+1, char) ;
	}
      if (si->convert)
	{ char *s = si->seqBuf, *e = s + si->seqLen, *sv = oneString(vf) ;
	  while (s < e) *s++ = si->convert[(int)*sv++] ;
	}
      else
	memcpy (si->seqBuf, oneString(vf), si->seqLen) ;
      if (!oneReadLine (vf)) return false ;
      if (si->isQual)
	{ while (vf->lineType != 'Q' && vf->lineType != 'S') if (!oneReadLine (vf)) return false;
	  if (vf->lineType == 'Q')
	    { char *q = si->qualBuf, *e = q + si->seqLen, *qv = oneString(vf) ;
	      while (q < e) *q++ = *qv++ - 33 ;
	    }
	}
      while (vf->lineType != 'S') if (!oneReadLine (vf)) return false ;
      return true ;
    }
#endif
#ifdef BAMIO
  if (si->type == BAM) return bamRead (si) ;
#endif

  if (!si->nb) return false ;
  si->recStart = si->b - si->buf ;
  
  if (si->type == BINARY)
    { if (si->line > si->nSeq) return false ; /* have already read all the sequences */
      bufConfirmNbytes (si, (U64)(3*sizeof(int))) ;
      si->idLen = *((int*)si->b) ; si->b += sizeof(int) ;
      si->descLen = *((int*)si->b) ; si->b += sizeof(int) ;
      si->seqLen = *((int*)si->b) ; si->b += sizeof(int) ;
      si->nb -= 3*sizeof(int) ;
      U64 nBytes = si->idLen + 1 + si->descLen + 1 + (si->seqLen+3)/4 ;
      if (si->qualThresh) nBytes += (si->seqLen+7)/8 ; /* NB qualThresh not isQual */
      nBytes = 4*((nBytes+3)/4) ;
      bufConfirmNbytes (si, nBytes) ;
      si->idStart = si->b - si->buf ;
      si->descStart = si->idStart + si->idLen + 1 ;
      si->seqStart = si->descStart + si->descLen + 1 ;
      sqioSeqUnpack ((U8*)(si->buf+si->seqStart), si->seqBuf, si->seqLen, si) ;
      if (si->isQual)
	{ si->qualStart = si->seqStart + (si->seqLen + 3) / 4 ;
	  sqioQualUnpack ((U8*)(si->buf+si->qualStart), si->qualBuf, si->seqLen, si) ;
	}
      ++si->line ;
      si->b += nBytes ; si->nb -= nBytes ;
      return true ;
    }
  
  /* if get to here then this is a text file, FASTA or FASTQ */
  
  if (si->type == FASTA)
    { if (*si->b != '>') die ("no initial > for FASTA record line %" PRIu64 "", si->line) ; }
  else if (si->type == FASTQ)
    { if (*si->b != '@') die ("no initial @ for FASTQ record line %" PRIu64 "", si->line) ; }
  bufAdvanceInRecord(si) ; si->idStart = si->b - si->buf ;
  while (!isspace(*si->b)) bufAdvanceInRecord(si) ;
  si->idLen = si->b - sqioId(si) ;
  if (*si->b != '\n') /* a space or tab - whatever follows on this line is description */
    { *si->b = 0 ; bufAdvanceInRecord(si) ;
      si->descStart = si->b - si->buf ;
      while (*si->b != '\n') bufAdvanceInRecord(si) ;
      si->descLen = si->b - sqioDesc(si) ;
    }
  else { si->descLen = si->descStart = 0 ; }
  *si->b = 0 ;
  ++si->line ; bufAdvanceInRecord(si) ;	              /* line 2 */
  si->seqStart = si->b - si->buf ;
  if (si->type == FASTA)
    { while (si->nb && *si->b != '>')
	{ while (*si->b != '\n') bufAdvanceInRecord(si) ;
	  ++si->line ; bufAdvanceEndRecord(si) ;
	}
      char *s = sqioSeq(si), *t = s ;
      while (s < si->b) if ((*t++ = si->convert[(int)*s++]) < 0) --t ;
      si->seqLen = t - sqioSeq(si) ;
    }
  else if (si->type == FASTQ)
    { while (*si->b != '\n') bufAdvanceInRecord(si) ;
      si->seqLen = si->b - sqioSeq(si) ;
      if (si->convert)
	{ char *s = sqioSeq(si) ;
	  while (s < si->b) { *s = si->convert[(int)*s] ; ++s ; }
	}
      ++si->line ; bufAdvanceInRecord(si) ; 	      /* line 3 */
      if (*si->b != '+') die ("missing + FASTQ line %" PRIu64 "", si->line) ;
      while (*si->b != '\n') bufAdvanceInRecord(si) ; /* ignore remainder of + line */
      ++si->line ; bufAdvanceInRecord(si) ;	      /* line 4 */
      si->qualStart = si->b - si->buf ;
      while (*si->b != '\n') bufAdvanceInRecord(si) ;
      if (si->b - si->buf - si->qualStart != si->seqLen)
	die ("qual not same length as seq line %" PRIu64 "", si->line) ;
      if (si->isQual) { char *q = sqioQual(si), *e = q + si->seqLen ; while (q < e) *q++ -= 33 ; }
      ++si->line ; bufAdvanceEndRecord(si) ;
    }

  ++si->nSeq ;
  return true ;
}

/*********************** open for writing ***********************/

#ifdef ONEIO
static char *schemaText =
  "1 3 def 1 0  schema for seqio\n"
  ".\n"
  "P 3 seq SEQUENCE\n"
  "S 3 irp   read pairs\n"
  "S 3 pbr   pacbio reads\n"
  "S 3 10x   10X Genomics data\n"
  "S 3 ctg   contigs from an assembly\n"
  "S 3 kmr   kmers\n"
  "D g 2 3 INT 6 STRING  group: count, name (e.g. use for flow cell/lane grouping)\n"
  "D S 1 3 DNA           sequence: the DNA string\n"
  "D I 1 6 STRING        id: (optional) sequence identifier\n"
  "D Q 1 6 STRING        quality: Q values (ascii string = q+33)\n" ;
#endif

SeqIO *seqIOopenWrite (char *filename, SeqIOtype type, int* convert, int qualThresh)
{
  SeqIO *si = new0 (1, SeqIO) ;
  int nameLen = strlen (filename) ;

  si->type = type ;
  si->isWrite = true ;
  si->convert = convert ;
  si->isQual = (qualThresh > 0) ;
  si->qualThresh = qualThresh ;
  if (si->type == FASTA && si->isQual)
    { fprintf (stderr, "warning : can't write qualities to FASTA file %s\n", filename) ;
      si->isQual = false ;
    }

  if (si->type == ONE ||
      (nameLen > 5 && *(filename+nameLen-5) == '.' && *(filename+nameLen-4) == '1'))
#ifdef ONEIO
    { OneSchema *schema = oneSchemaCreateFromText (schemaText) ;
      char *oneType = "seq" ;
      if (nameLen > 5 && *(filename+nameLen-5) == '.' && *(filename+nameLen-4) == '1')
	oneType = filename+nameLen-3 ;
      OneFile *vf = oneFileOpenWriteNew (filename, schema, oneType, true, 1) ;
      oneSchemaDestroy (schema) ;
      if (!vf) { free (si) ; return 0 ; }
      si->handle = vf ;
      char *commandLine = getCommandLine() ;
      if (!commandLine) commandLine = "-" ;
      oneAddProvenance (vf, "seqio", "1.0", commandLine, 0) ;
      oneWriteHeader (vf) ;
      return si ;
    }
#else
    { warn ("sorry, seqio not compiled with ONElib so can't write 1seq") ;
      free (si) ; return 0 ;
    }
#endif

  if (si->type == BAM)
    { warn ("sorry, seqio can't write BAM") ;
      free (si) ; return 0 ;
    }
  
  if (!strcmp (filename, "-"))
    { si->fd = fileno (stdout) ;
      if (si->fd == -1) { free (si) ; return 0 ; }
    }
  else if (!strcmp (filename, "-z"))
    { si->gzf = gzdopen (fileno(stdout), "w") ;
      if (!si->gzf) { free (si) ; return 0 ; }
    }
  else if (nameLen > 3 && !strcmp (filename+nameLen-3, ".gz"))
    { si->gzf = gzopen (filename, "w") ; nameLen -= 3 ; filename[nameLen] = 0 ;
      if (!si->gzf) { free (si) ; return 0 ; }
    }
  else
    { si->fd = open (filename, O_CREAT | O_TRUNC | O_WRONLY, 00644) ;
      if (si->fd == -1) { free (si) ; return 0 ; }
    }

  if (si->type == UNKNOWN)
    { if (nameLen > 3 && !strcmp (filename+nameLen-3, ".fa")) si->type = FASTA ;
      else if (nameLen > 3 && !strcmp (filename+nameLen-3, ".fq")) si->type = FASTQ ;
      else si->type = BINARY ;
    }
  if (si->type == BINARY && si->gzf)
    { fprintf (stderr, "can't write a gzipped binary file\n") ; seqIOclose (si) ; return 0 ; }
  
  si->nb = si->bufSize = 1<<24 ;
  si->b = si->buf = new (si->bufSize, char) ;

  if (si->type == BINARY)
    { si->b += 64 ; si->nb -= 64 ; /* make space for header - written on file close */
    }

  return si ;
}

void seqIOflush (SeqIO *si)	/* writes buffer to file and resets to 0 */
{
  if (!si->isWrite) return ;
  U64 retVal, nBytes = si->b - si->buf ;
  if (si->gzf) retVal = gzwrite (si->gzf, si->buf, nBytes) ;
  else retVal = write (si->fd, si->buf, nBytes) ;
  if (retVal != nBytes) die ("seqio write error %" PRIu64 " not %" PRIu64 " bytes written", retVal, nBytes) ;
  si->b = si->buf ;
  si->nb = si->bufSize ;
}

static void writeExtend (SeqIO *si, U64 len)
{
  while (si->bufSize < len) si->bufSize <<= 1 ;
  free (si->buf) ;
  si->b = si->buf = new (si->bufSize, char) ;
  si->nb = si->bufSize ;
}
  
void seqIOwrite (SeqIO *si, char *id, char *desc, U64 seqLen, char *seq, char *qual)
{
  U64 len, pad = 0 ;
  assert (si->isWrite) ;

  ++si->nSeq ;
  si->idLen = id ? strlen(id) : 0 ;
  si->totIdLen += si->idLen ; if (si->idLen > si->maxIdLen) si->maxIdLen = si->idLen ;
  si->descLen = desc ? strlen(desc) : 0 ;
  si->totDescLen += si->descLen ; if (si->descLen > si->maxDescLen) si->maxDescLen = si->descLen ;
  si->seqLen = seqLen ;
  si->totSeqLen += si->seqLen ; if (si->seqLen > si->maxSeqLen) si->maxSeqLen = si->seqLen ;

#ifdef ONEIO
  if (si->type == ONE)
    { OneFile *vf = (OneFile*)(si->handle) ;
      static U64 bufLen = 0 ;
      static char *buf = 0 ;
      if (seqLen > bufLen)
	{ if (buf) free (buf) ;
	  bufLen = seqLen ;
	  buf = new(bufLen+1,char) ;
	}
      if (si->convert)
	{ int i ;
	  for (i = 0 ; i < seqLen ; ++i) buf[i] = si->convert[(int)(seq[i])] ;
	  oneWriteLine (vf, 'S', seqLen, buf) ;
	}
      else
	oneWriteLine (vf, 'S', seqLen, seq) ;
      if (id)
	{ oneWriteLine (vf, 'I', strlen (id), id) ;
	  if (desc) oneWriteComment (vf, desc) ;
	}
      if (qual && si->isQual)
	{ int i ;
	  for (i = 0 ; i < seqLen ; ++i) buf[i] = qual[i] + 33 ;
	  oneWriteLine (vf, 'Q', seqLen, buf) ;
	}
      return ;
    }
#endif
  
  if (si->type == FASTA)
    len = 3 + si->idLen + (desc ? 1+si->descLen : 0) + seqLen ;   /* >id[ desc]\nseq\n */
  else if (si->type == FASTQ)
    len = 6 + si->idLen + (desc ? 1+si->descLen : 0) + 2*seqLen ; /* @id[ desc]\nseq\n+\nqual\n */
  else	/* binary */
    { U64 nBytes = si->idLen + si->descLen + 2 + (si->seqLen+3)/4 ;
      if (si->isQual) nBytes += (si->seqLen+7)/8 ;
      pad = 3 - ((nBytes+3) % 4) ; /* 0,1,2,3 maps to 0,3,2,1 */
      len = 3*sizeof(int) + nBytes + pad ;
    }
  if (len > si->nb) seqIOflush (si) ;
  if (len > si->nb) writeExtend (si, len) ;

  if (si->type == FASTA)
    { *si->b++ = '>' ;
      if (id) { strcpy (si->b, id) ; si->b += si->idLen ; }
      if (desc) { *si->b++ = ' ' ; strcpy (si->b, desc) ; si->b += si->descLen ; }
      *si->b++ = '\n' ;
      memcpy (si->b, seq, seqLen) ;
      if (si->convert) while (seqLen--) { *si->b = si->convert[(int)*si->b] ; ++si->b ; }
      else si->b += seqLen ;
      *si->b++ = '\n' ;
    }
  else if (si->type == FASTQ)
    { *si->b++ = '@' ;
      if (id) { strcpy (si->b, id) ; si->b += si->idLen ; }
      if (desc) { *si->b++ = ' ' ; strcpy (si->b, desc) ; si->b += si->descLen ; }
      *si->b++ = '\n' ;
      memcpy (si->b, seq, seqLen) ;
      if (si->convert) while (seqLen--) { *si->b = si->convert[(int)*si->b] ; ++si->b ; }
      else si->b += seqLen ;
      *si->b++ = '\n' ;
      *si->b++ = '+' ;
      *si->b++ = '\n' ;
      seqLen = si->seqLen ; while (seqLen--) *si->b++ = qual ? *qual++ + 33 : 33 ;
      *si->b++ = '\n' ;
    }
  else				/* binary */
    { int *ib = (int*)si->b ; si->b += 3 * sizeof(int) ;
      *ib++ = si->idLen ; *ib++ = si->descLen ; *ib++ = seqLen ;
      if (si->idLen) { strcpy (si->b, id) ; si->b += si->idLen ; } *si->b++ = 0 ;
      if (si->descLen) { strcpy (si->b, desc) ; si->b += si->descLen ; } *si->b++ = 0 ;
      si->b += sqioSeqPack (seq, (U8*)si->b, si->seqLen, si->convert) ;
      if (si->isQual) si->b += sqioQualPack (qual, (U8*)si->b, si->seqLen, si->qualThresh) ;
      si->b += pad ;
    }
  si->nb -= len ;
}

/********** some routines to pack sequence and qualities for binary representation ***********/

U64 sqioSeqPack (char *s, U8 *u, U64 len, int *convert) /* compress s into (len+3)/4 u */
{
  U8 *u0 = u ;
  int i ;
  if (!convert) convert = dna2index4Conv ;
  while (len > 4)
    { *u = 0 ; for (i = 0 ; i < 4 ; ++i) *u = (*u << 2) | convert[(int)*s++] ;
      len -= 4 ; ++u ;
    }
  if (len)
    { *u = 0 ; for (i = 0 ; i < len ; ++i) *u = (*u << 2) | convert[(int)*s++] ;
      ++u ;
    }
  return (u-u0) ;
}

void sqioSeqUnpack (U8 *u, char *s, U64 len, SeqIO *si) /* uncompress (len+3)/4 u into s */
{
  int i ;
  while (len > 4)		/* NB needs to be > here not >= so can prime seqExpand */
    { *(U32*)s = si->seqExpand[*u] ;
      ++u ; s += 4 ; len -= 4 ;
    }
  if (len) for (i = len ; --i ; ) { s[i] = si->unpackConvert[*u & 0x3] ; *u >>= 2 ; }
}

U64 sqioQualPack (char *q, U8 *u, U64 len, int thresh) /* compress q into (len+7)/8 u  */
{
  U8 *u0 = u ;
  int i ;
  while (len > 8)
    { *u = 0 ; for (i = 0 ; i < 8 ; ++i) { if (*q++ >= thresh) *u |= 1 ; *u <<= 1 ; }
      len -= 8 ; ++u ;
    }
  if (len)
    { *u = 0 ; for (i = 0 ; i < len ; ++i) { if (*q++ >= thresh) *u |= 1 ; *u <<= 1 ; }
      ++u ;
    }
  return (u-u0) ;
}

void sqioQualUnpack (U8 *u, char *q, U64 len, SeqIO *si) /* uncompress (len+7)/8 u into q */
{
  int i ;
  while (len > 8)		/* NB needs to be > here not >= so can prime qualExpand */
    { *(U64*)q = si->qualExpand[*u] ;
      ++u ; q += 8 ; len -= 8 ;
    }
  if (len) for (i = len ; --i ; ) { q[i] = (*u & 0x1) ? si->qualThresh : 0 ; *u >>= 1 ; }
}

/*********** standard conversion tables **************/

int dna2textConv[] = {
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2, 'A',  -2, 'C',  -2,  -2,  -2, 'G',  -2,  -2,  -2,  -2,  -2,  -2, 'N',  -2,
  -2,  -2,  -2,  -2, 'T',  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,
  -2, 'a',  -2, 'c',  -2,  -2,  -2, 'g',  -2,  -2,  -2,  -2,  -2,  -2, 'n',  -2,
  -2,  -2,  -2,  -2, 't',  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,
} ;

int dna2textAmbigConv[] = {
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, '-',  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2, 'A', 'B', 'C', 'D',  -2,  -2, 'G', 'H',  -2,  -2, 'K',  -2, 'M', 'N',  -2,
  -2,  -2, 'R', 'S', 'T',  -2, 'V', 'W',  -2, 'Y',  -2,  -2,  -2,  -2,  -2,  -2,
  -2, 'A', 'B', 'C', 'D',  -2,  -2, 'G', 'H',  -2,  -2, 'K',  -2, 'M', 'N',  -2,
  -2,  -2, 'R', 'S', 'T',  -2, 'V', 'W',  -2, 'Y',  -2,  -2,  -2,  -2,  -2,  -2,
} ;

int dna2textAmbig2NConv[] = {
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2, 'A', 'N', 'C', 'N',  -2,  -2, 'G', 'N',  -2,  -2, 'N',  -2, 'N', 'N',  -2,
  -2,  -2, 'N', 'N', 'T',  -2, 'N', 'N',  -2, 'N',  -2,  -2,  -2,  -2,  -2,  -2,
  -2, 'A', 'N', 'C', 'N',  -2,  -2, 'G', 'N',  -2,  -2, 'N',  -2, 'N', 'N',  -2,
  -2,  -2, 'N', 'N', 'T',  -2, 'N', 'N',  -2, 'N',  -2,  -2,  -2,  -2,  -2,  -2,
} ;

int dna2indexConv[] = {
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,   0,  -2,   1,  -2,  -2,  -2,   2,  -2,  -2,  -2,  -2,  -2,  -2,   4,  -2,
  -2,  -2,  -2,  -2,   3,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,
  -2,   0,  -2,   1,  -2,  -2,  -2,   2,  -2,  -2,  -2,  -2,  -2,  -2,   4,  -2,
  -2,  -2,  -2,  -2,   3,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,
} ;

int dna2index4Conv[] = {	/* sends N,n to A */
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,   0,  -2,   1,  -2,  -2,  -2,   2,  -2,  -2,  -2,  -2,  -2,  -2,   0,  -2,
  -2,  -2,  -2,  -2,   3,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,
  -2,   0,  -2,   1,  -2,  -2,  -2,   2,  -2,  -2,  -2,  -2,  -2,  -2,   0,  -2,
  -2,  -2,  -2,  -2,   3,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,
} ;

int dna2binaryConv[] = {
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,   1,  -2,   2,  -2,  -2,  -2,   4,  -2,  -2,  -2,  -2,  -2,  -2,  15,  -2,
  -2,  -2,  -2,  -2,   8,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,
  -2,   1,  -2,   2,  -2,  -2,  -2,   4,  -2,  -2,  -2,  -2,  -2,  -2,  15,  -2,
  -2,  -2,  -2,  -2,   8,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,
} ;

int dna2binaryAmbigConv[] = {
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,   0,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,   1,  14,   2,  13,  -2,  -2,   4,  11,  -2,  -2,  12,  -2,   3,  15,  -2,
  -2,  -2,   5,   6,   8,  -2,   7,   9,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,
  -2,   1,  14,   2,  13,  -2,  -2,   4,  11,  -2,  -2,  12,  -2,   3,  15,  -2,
  -2,  -2,   5,   6,   8,  -2,   7,   9,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,
} ;

int aa2textConv[] = {
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2, 'A', 'X', 'C', 'D', 'E', 'F', 'G', 'H', 'I',  -2, 'K', 'L', 'M', 'N',  -2,
 'P', 'Q', 'R', 'S', 'T',  -2, 'V', 'W', 'X', 'Y', 'X',  -2,  -2,  -2,  -2,  -2,
  -2, 'A', 'X', 'C', 'D', 'E', 'F', 'G', 'H', 'I',  -2, 'K', 'L', 'M', 'N',  -2,
 'P', 'Q', 'R', 'S', 'T',  -2, 'V', 'W', 'X', 'Y', 'X',  -2,  -2,  -2,  -2,  -2,
} ;

int aa2indexConv[] = {
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2, 
  -2,   0,  20,   1,   2,   3,   4,   5,   6,   7,  -2,   8,   9,  10,  11,  -2,
  12,  13,  14,  15,  16,  -2,  17,  18,  20,  19,  20,  -2,  -2,  -2,  -2,  -2,
  -2,   0,  20,   1,   2,   3,   4,   5,   6,   7,  -2,   8,   9,  10,  11,  -2,
  12,  13,  14,  15,  16,  -2,  17,  18,  20,  19,  20,  -2,  -2,  -2,  -2,  -2,
} ;

int noConv[] = {
   0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,
  48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,
  64,  65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,  76,  77,  78,  79,
  80,  81,  82,  83,  84,  85,  86,  87,  88,  89,  90,  91,  92,  93,  94,  95,
  96,  97,  98,  99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111,
 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127
} ;

/************ BAM access code ****************/

#ifdef BAMIO

#include "sam.h"

typedef struct {
  samFile *f ;
  sam_hdr_t *h ;
  bam1_t *b ;
} BamFile ;

bool bamFileOpenRead (char* filename, SeqIO *si)
{
  static char *referenceFileName = 0 ; // may need to set this somehow for CRAM
  BamFile *bf = new0 (1, BamFile) ;

  bf->f = sam_open (filename, "r") ;
  if (!bf->f) return false ;

  uint32_t rf = SAM_FLAG | SAM_SEQ ;
  if (si->isQual) rf |= SAM_QUAL ;
  if (hts_set_opt (bf->f, CRAM_OPT_REQUIRED_FIELDS, rf))
    { fprintf (stderr, "BamFileOpen failed to set CRAM_OPT_REQUIRED_FIELDS\n") ;
      bamFileClose (si) ;
      return false ;
    }
  if (hts_set_opt (bf->f, CRAM_OPT_DECODE_MD, 0))
    { fprintf (stderr, "BamFileOpen failed to set CRAM_OPT_DECODE_MD\n") ;
      bamFileClose (si) ;
      return false ;
    }
  if (referenceFileName && hts_set_fai_filename (bf->f, referenceFileName) != 0)
    { fprintf (stderr, "BamFileOpen failed to set reference genome from %s\n",
	       referenceFileName) ;
      bamFileClose (si) ;
      return false ;
    }

  if (!(bf->h = sam_hdr_read (bf->f)))
    { fprintf(stderr, "BamFileOpen failed to read header\n") ;
      bamFileClose (si) ;
      return false ;
    }

  if (!(bf->b = bam_init1 ()))
    { fprintf(stderr, "BamFileOpen failed to create record buffer\n") ;
      bamFileClose (si) ;
      return false ;
    }

  si->handle = (void*) bf ;
  return true ;
}

static const char binaryAmbigComplement[16] =
  { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 } ;
static const char binaryAmbig2text[] = "=ACMGRSVTWYHKDBN" ;

bool bamRead (SeqIO *si)
{
  int i ;
  BamFile *bf = (BamFile*) si->handle ;

  int res = sam_read1 (bf->f, bf->h, bf->b) ;
  if (res < -1)
    die ("BamFileRead failed to read bam record %d\n", si->nSeq) ;
  if (res == -1) // end of file
    return false ;
  
  si->seqLen = bf->b->core.l_qseq ;
  if (si->seqLen > si->maxSeqLen)
    { if (si->maxSeqLen) { free (si->seqBuf) ; if (si->isQual) free (si->qualBuf) ; }
      si->maxSeqLen = si->seqLen ;
      si->seqBuf = new0 (si->maxSeqLen+1, char) ;
      if (si->isQual) si->qualBuf = new0 (si->maxSeqLen+1, char) ;
    }

  char *bseq = (char*) bam_get_seq (bf->b) ;
  char *s = si->seqBuf ;
  if (bf->b->core.flag & BAM_FREVERSE)
    for (i = si->seqLen ; i-- ; )
      *s++ = binaryAmbig2text[(int)binaryAmbigComplement[(int)bam_seqi(bseq,i)]] ;
  else
    for (i = 0 ; i < si->seqLen ; ++i)
      *s++ = binaryAmbig2text[bam_seqi(bseq,i)] ;
  if (si->convert)
    while (s-- > si->seqBuf) 		// NB relies on s == seqBuf+seqLen at start
      *s = si->convert[(int)*s] ;
  
  if (si->isQual)
    { char *bq = (char*) bam_get_qual (bf->b) ;
      if (*bq == '\xff') bzero (si->qualBuf, si->seqLen) ;
      else if (bf->b->core.flag & BAM_FREVERSE)
	{ char *q = si->qualBuf + si->seqLen ;
	  while (q-- > si->qualBuf) *q = *bq ;
	}
      else
	memcpy (si->qualBuf, bq, si->seqLen) ;
    }

  // NB bam_get_qname(bf->b) returns the sequence name

  ++si->nSeq ;
  return true ;
}

void bamFileClose (SeqIO *si)
{ BamFile *bf = (BamFile*) si->handle ;
  if (bf->b) bam_destroy1 (bf->b) ;
  if (bf->h) sam_hdr_destroy (bf->h) ;
  if (bf->f) sam_close (bf->f) ;
  free (bf) ;
}

#endif

/******** end of file ********/
