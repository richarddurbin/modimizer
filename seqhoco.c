/*  File: seqhoco.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2020
 *-------------------------------------------------------------------
 * Description: comp
 * Exported functions:
 * HISTORY:
 * Last edited: Dec  7 23:58 2020 (rd109)
 * Created: Thu Jun  4 23:40:57 2020 (rd109)
 *-------------------------------------------------------------------
 */

#include "seqio.h"
#include <ctype.h>

int main (int argc, char *argv[])
{
  --argc ; ++argv ;
  if (!argc) *--argv = "-" ; // standard trick to read stdin if no argument
  
  SeqIO *siIn = seqIOopenRead (*argv, dna2textConv, false) ;
  if (!siIn) die ("failed to read sequence file %s", *argv) ;
  //  dna2textConv['N'] = 0 ; dna2textConv['n'] = 0 ; // why this?
  SeqIO *siOut = seqIOopenWrite ("-z", FASTA, dna2textConv, 0) ;
  if (!siOut) die ("failed to open stdio to write compressed fasta output") ;
  while (seqIOread (siIn) && siIn->seqLen)
    { char *t, *s ;
      t = s = sqioSeq(siIn) ;
      int n = siIn->seqLen ;
      while (n-- > 0) if (toupper(*++s) != toupper(*t)) *++t = *s ;
      seqIOwrite (siOut, sqioId(siIn), 0, t-sqioSeq(siIn)+1, sqioSeq(siIn), 0) ;
    }
  seqIOclose (siIn) ;
  seqIOclose (siOut) ;
}

// GCGGTTTGAGTGAnNCGAGACGAGAcgcgcccctcccaCGCGGGGAAGGG
//     GCGTGAGTGAnCGAGACGAGAcgcgctcaCGCGAG
