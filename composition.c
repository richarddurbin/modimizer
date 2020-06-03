/*  File: composition.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2018
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Jul 22 21:07 2019 (rd109)
 * Created: Sun Nov 11 17:21:40 2018 (rd109)
 *-------------------------------------------------------------------
 */

#include "seqio.h"
#include <ctype.h>

int main (int argc, char *argv[])
{
  --argc ; ++argv ;

  timeUpdate (stdout) ;

  if (!argc) *--argv = "-" ;	/* minor abuse, so that default is to read stdin */
  
  SeqIO *si = seqIOopenRead (*argv, 0, TRUE) ;
  if (!si) die ("failed to open sequence file %s\n", *argv) ;

  U64 totLen = 0, n = 0 ;
  U64 *totBase = new0 (256, U64) ;
  U64 *totQual = new0 (256, U64) ;
  while (seqIOread (si))
    { char *s = sqioSeq(si), *e = s + si->seqLen ;
      while (s < e) ++totBase[*s++] ;
      totLen += si->seqLen ;
      if (si->isQual)
	{ char *q = sqioQual(si), *e = q + si->seqLen ;
	  while (q < e) ++totQual[*q++] ;
	}
    }
  printf ("%s file, %llu sequences >= 0, %llu total, %.2f average\n",
	  seqIOtypeName[si->type], si->nSeq, totLen, totLen / (double) si->nSeq) ;
  int i ;
  U64 totUnprint = 0 ;
  printf ("bases\n") ;
  for (i = 0 ; i < 256 ; ++i)
    if (totBase[i])
      { if (isprint(i)) printf ("  %c %llu %4.1f %%\n", i, totBase[i], totBase[i]*100.0/totLen) ;
	else totUnprint += totBase[i] ;
      }
  if (totUnprint) printf (" unprintable %llu %4.1f %%\n", totUnprint, totUnprint*100.0/totLen) ;

  if (si->isQual)
    { printf ("qualities\n") ;
      U64 sum = 0 ;
      for (i = 0 ; i < 256 ; ++i)
	{ sum += totQual[i] ;
	  if (totQual[i]) printf (" %3d %llu %4.1f %% %5.1f %%\n",
				  i, totQual[i], totQual[i]*100.0/totLen, sum*100.0/totLen) ;
	}
    }
  
  seqIOclose (si) ;
  timeTotal (stdout) ;
}

/****************/
