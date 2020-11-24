/*  File: composition.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2018
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Aug  7 13:46 2020 (rd109)
 * Created: Sun Nov 11 17:21:40 2018 (rd109)
 *-------------------------------------------------------------------
 */

#include "seqio.h"
#include "array.h"
#include <stdbool.h>
#include <ctype.h>
#include <math.h>

static int lengthBins = 20 ;

void usage (void)
{
  fprintf (stderr, "Usage: composition [opts] <filename>\n") ;
  fprintf (stderr, "  will read fasta, fastq, bam/sam/cram, 1code, custom-binary.  Use filename '-' for stdin (not 1code binary)\n") ;
  fprintf (stderr, "  options:\n") ;
  fprintf (stderr, "    -b : show base counts\n") ;
  fprintf (stderr, "    -q : show quality counts\n") ;
  fprintf (stderr, "    -t : show time and memory used\n") ;
  fprintf (stderr, "    -l : show length distribution in up to %d quadratic bins\n", lengthBins) ;
}

int main (int argc, char *argv[])
{
  --argc ; ++argv ;
  U64   *totBase = 0, *totQual = 0 ;
  bool  isTime = false ;
  Array lengthCount = 0, lengthSum = 0 ;

  if (!argc) usage () ;

  while (argc && **argv == '-' && strcmp (*argv, "-"))
    if (!strcmp (*argv, "-b")) { totBase = new0 (256, U64) ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-q")) { totQual = new0 (256, U64) ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-t")) { isTime = true ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-l"))
      { lengthCount = arrayCreate (10000, int) ;
	lengthSum = arrayCreate (10000, U64) ;
	--argc ; ++argv ;
      }
    else usage () ;
  
  if (isTime) timeUpdate (stdout) ;
  
  SeqIO *si = seqIOopenRead (*argv, 0, true) ;
  if (!si) die ("failed to open sequence file %s\n", *argv) ;

  U64 lenMin = 0, lenMax = 0, totLen = 0, n = 0 ;
  int i ;
  while (seqIOread (si))
    { char *s = sqioSeq(si), *e = s + si->seqLen ;
      if (totBase) while (s < e) ++totBase[*s++] ;
      totLen += si->seqLen ;
      if (si->seqLen > lenMax) lenMax = si->seqLen ;
      if (!lenMin || si->seqLen < lenMin) lenMin = si->seqLen ;
      if (lengthCount)
	{ i = 10.*sqrt((double)si->seqLen) ;
	  ++array(lengthCount, i, int) ; array(lengthSum, i, U64) += si->seqLen ;
	}
      if (totQual && si->isQual)
	{ char *q = sqioQual(si), *e = q + si->seqLen ;
	  while (q < e) ++totQual[*q++] ;
	}
    }
  printf ("%s file, %llu sequences >= 0, %llu total, %.2f average, %llu min, %llu max\n",
	  seqIOtypeName[si->type], si->nSeq, totLen, totLen / (double) si->nSeq, lenMin, lenMax) ;
  if (totBase)
    { U64 totUnprint = 0 ;
      printf ("bases\n") ;
      for (i = 0 ; i < 256 ; ++i)
	if (totBase[i])
	  { if (isprint(i)) printf ("  %c %llu %4.1f %%\n", i, totBase[i], totBase[i]*100.0/totLen) ;
	    else totUnprint += totBase[i] ;
	  }
      if (totUnprint) printf (" unprintable %llu %4.1f %%\n", totUnprint, totUnprint*100.0/totLen) ;
      free (totBase) ;
    }

  if (totQual && si->isQual)
    { printf ("qualities\n") ;
      U64 sum = 0 ;
      for (i = 0 ; i < 256 ; ++i)
	{ sum += totQual[i] ;
	  if (totQual[i]) printf (" %3d %llu %4.1f %% %5.1f %%\n",
				  i, totQual[i], totQual[i]*100.0/totLen, sum*100.0/totLen) ;
	}
      free (totQual) ;
    }

  if (lengthSum)
    { if (lenMin < lenMax)
	{ U64 tot50 = 0 ;
	  for (i = 0 ; i < arrayMax(lengthSum) && tot50 < 0.5*totLen ; ++i) 
	    tot50 += arr(lengthSum, i, U64) ;
	  printf ("approximate N50 %d\n", (i*(i+1))/100) ;
	  printf ("length distribution (quadratic bins)\n") ;
	  int s = 0 ;
	  int d = arrayMax(lengthCount) / 20 ;
	  for (i = 0 ; i < arrayMax(lengthCount) ; ++i)
	    { s += arr(lengthCount, i, int) ;
	      if (s && !((arrayMax(lengthCount)-1-i) % d))
		{ printf ("  %d\t%d\n", (i*i)/100, s) ;
		  s = 0 ;
		}
	    }
	}
      arrayDestroy (lengthCount) ; arrayDestroy (lengthSum) ;
    }
  
  if (isTime) timeTotal (stdout) ;
  seqIOclose (si) ;
}

/****************/
