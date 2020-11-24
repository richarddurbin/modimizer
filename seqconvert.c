/*  File: seqconvert.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2019
 *-------------------------------------------------------------------
 * Description: utility to convert between sequence formats
 * Exported functions:
 * HISTORY:
 * Last edited: Jun 20 13:10 2020 (rd109)
 * Created: Sun Feb 17 10:23:37 2019 (rd109)
 *-------------------------------------------------------------------
 */

#include "seqio.h"

int main (int argc, char *argv[])
{
  storeCommandLine (argc, argv) ;
  --argc ; ++argv ;

  timeUpdate (stderr) ;

  if (!argc || !strcmp(*argv,"-h") || !strcmp(*argv,"--help"))
    { fprintf (stderr, "Usage: seqconvert [-fa|fq|b|1] [-Q T] [-z] [-S] [-o outfile] [infile]\n") ;
      fprintf (stderr, "   .gz ending outfile name implies gzip compression\n") ;
      fprintf (stderr, "   -fa output as fasta, -fq as fastq, -b as binary, -1 as ONEcode\n") ;
      fprintf (stderr, "      else .fa or .fq in outfile name imply fasta, fastq else binary\n") ;
      fprintf (stderr, "   -Q sets the quality threshold for single bit quals in -b option [0]\n") ;
      fprintf (stderr, "   -S silent - else it reports to stderr on what it is doing\n") ;
      fprintf (stderr, "   NB gzip is not compatible with binary\n") ;
      fprintf (stderr, "   if no infile then use stdin\n") ;
      fprintf (stderr, "   if no -o option then use stdout and -z implies gzip\n");
      exit (0) ;
    }
  
  SeqIOtype type = UNKNOWN ;
  bool isVerbose = true ;
  bool isGzip = false ;
  char *inFileName = "-" ;
  char *outFileName = "-z" ;
  int qualThresh = 0 ;
  while (argc)
    { if (!strcmp (*argv, "-fa")) type = FASTA ;
      else if (!strcmp (*argv, "-fq")) type = FASTQ ;
      else if (!strcmp (*argv, "-b")) type = BINARY ;
      else if (!strcmp (*argv, "-1")) type = ONE ;
      else if (!strcmp (*argv, "-Q") && argc >1)
	{ --argc ; ++argv ; qualThresh = atoi (*argv) ; }
      else if (!strcmp (*argv, "-z")) isGzip = true ;
      else if (!strcmp (*argv, "-o") && argc > 1)
	{ --argc ; ++argv ; outFileName = *argv ; }
      else if (!strcmp (*argv, "-S")) isVerbose = false ;
      else if (argc == 1 && **argv != '-') inFileName = *argv ;
      else die ("unknown option %s - run without arguments for help\n", *argv) ;
      --argc ; ++argv ;
    }

  if (!strcmp(outFileName, "-z") && !isGzip) outFileName = "-" ; /* remove 'z' */
  SeqIO *siOut = seqIOopenWrite (outFileName, type, 0, qualThresh) ;
  if (!siOut) die ("failed to open output file %s", outFileName) ;
  bool isQual = (siOut->type == BINARY && qualThresh > 0) ||
    siOut->type == FASTQ || siOut->type == ONE ;
  SeqIO *siIn = seqIOopenRead (inFileName, 0, isQual) ;
  if (!siIn) die ("failed to open input file %s", inFileName) ;
  if (isVerbose)
    { fprintf (stderr, "reading from file type %s", seqIOtypeName[siIn->type]) ;
      if (siIn->type == BINARY)
	fprintf (stderr, "  with %llu sequences totLen %llu", siIn->nSeq, siIn->totSeqLen) ;
      fprintf (stderr, "\n") ;
    }

  while (seqIOread (siIn))
    seqIOwrite (siOut,
		siIn->idLen ? sqioId(siIn) : 0,
		siIn->descLen ? sqioDesc(siIn) : 0,
		siIn->seqLen, sqioSeq(siIn),
		siIn->isQual ? sqioQual(siIn) : 0) ;
  seqIOclose (siIn) ;
  seqIOclose (siOut) ;

  if (isVerbose)
    { fprintf (stderr, "written %lld sequences to file type %s, total length %lld, max length %lld\n",
	       siOut->nSeq, seqIOtypeName[siOut->type], siOut->totSeqLen, siOut->maxSeqLen) ;
  
      timeTotal (stderr) ;
    }
}

/****************/
