/*  File: modtype.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2021
 *-------------------------------------------------------------------
 * Description: idea is to type structural variants by kmers spanning the breakpoints
 * Exported functions:
 * HISTORY:
 * Last edited: Aug  4 14:55 2021 (rd109)
 * Created: Mon Jun 21 11:31:08 2021 (rd109)
 *-------------------------------------------------------------------
 */

#include "modset.h"
#include "ONElib.h"
#define BAMIO		// need to be able to read CRAM files
#include "seqio.h"

FILE *outFile ;
int   numThreads = 1 ;
bool  isVerbose = false ;

void usage (void)
{ fprintf (stderr, "Usage: modtype OPTIONS <reference> <sitefile> <samplefile>\n") ;
  fprintf (stderr, "  -v | --verbose : toggle verbose mode\n") ;
  fprintf (stderr, "  -t | --threads <number of threads for parallel ops> [%d]\n", numThreads) ;
  fprintf (stderr, "  -o | --output <output filename> : '-' for stdout\n") ;
  exit (1) ;
}

/* Strategy
   At each site in <sites> record the preceding and succeeding 31mer = triggers, using ModSet
   For each sample, go through all the reads 
       for every trigger hit build table of successor counts 
       collapse to major variants - perhaps by 3mer or 4mer distribution
       somehow summarise these to give genotypes
 */

/************************************************************/

char *schemaText =
  "1 3 def 1 0  schema for modtype\n"
  ".\n"
  "P 3 var                    variant file\n"
  "S 3 ins                    insertion file\n"
  "G c 2 3 INT 6 STRING          chromosome\n"
  "O I 2 3 INT 3 INT             insertion between left_pos and right_pos\n"
  "D A 1 4 CHAR                  0 for ref ancestral, 1 for alt ancestral\n"
  "D G 1 6 STRING                genotype: 0, 1 or 2 as a char per sample\n"
  "D K 2 4 CHAR 3 DNA            L|R, reference kmer up to left/right position\n"
  "D k 2 4 CHAR 3 DNA            L|R, insertion kmer following left/right position\n"
  "D L 1 8 INT_LIST              per sample numbers of left breakpoint insertion spans\n"
  "D R 1 8 INT_LIST              per sample numbers of right breakpoint insertion spans\n"
  "D F 1 8 INT_LIST              per sample numbers of reference spans\n"
  ".\n"
  "P 3 smp                    sample file\n"
  "O N 1 6 STRING                sample name\n"
  "D F 1 6 STRING                filename\n"
  "D C 1 4 REAL                  coverage\n"           
  ".\n"
  "P 3 nul                    empty file - comments only\n" ;

OneSchema *schema ;		/* make this from schemaText in main() */

/************************************************************/

#define TOPBIT  0x80000000	/* set for FORWARD orientation, unset for reverse orientation */
#define TOPMASK 0x7fffffff

typedef struct {
  DICT  *names ;       		/* chromosome/contig names */
  Array  seq ; 			/* of char* */
  Array  len ;			/* of int = sequence lengths */
} Reference ;

Reference *ref ;

typedef struct {
  int chrom ;
  U32 leftPos, rightPos ;
  U64 pre, post ;
} Site ;

Array sites ;			/* of Site */

typedef struct {
  char*  fileName ;
  double coverage ;
} Sample ;

typedef struct {
  DICT  *names ;
  Array sample ;		/* of Sample */
} SampleSet ;

SampleSet *samples ;

/************************************************************/

Reference *referenceRead (char* fileName)
{
  SeqIO *si = seqIOopenRead (fileName, dna2textConv, false) ;
  if (!si) die ("failed to open reference sequence file %s", fileName) ;
  Reference *ref = new0 (1, Reference) ;
  ref->names = dictCreate (64) ;
  ref->seq = arrayCreate (64, char*) ;
  ref->len = arrayCreate (64, int) ;
  int totLen = 0 ;
  while (seqIOread (si))
    { int i ;
      if (!dictAdd (ref->names, sqioId(si), &i))
	die ("duplicate sequence name %s in reference", sqioId(si)) ;
      array(ref->len, i, int) = si->seqLen ;
      array(ref->seq, i, char*) = new0 (si->seqLen+1, char) ;
      memcpy (arr(ref->seq, i, char*), sqioSeq(si), si->seqLen) ;
      totLen += si->seqLen ;
    }
  fprintf (stderr, "  reference read %d sequences total length %d from %s\n",
	   arrayMax(ref->len), totLen, fileName) ;
  seqIOclose (si) ;
  return ref ;
}

/************************************************************/

Array sitesRead (char* fileName, Reference *ref)
{
  OneFile *vf = oneFileOpenRead (fileName, schema, "ins", 1) ;
  if (!vf) die ("failed to open sites file %s", fileName) ;
  Array a = arrayCreate (256, Site) ;
  int chrom, cmax ;
  while (oneReadLine (vf))
    switch (vf->lineType)
      {
      case 'c':
	if (!dictFind (ref->names, oneString(vf), &chrom))
	  die ("bad contig/chrom name %s at line %d in %s", oneString(vf), vf->line, fileName) ;
	cmax = arr(ref->len,chrom,int) ;
	break ;
      case 'I':
	{ Site *s = arrayp(a, arrayMax(a), Site) ;
	  s->chrom = chrom ; s->leftPos = oneInt(vf,0) ; s->rightPos = oneInt(vf,1) ;
	  if (s->leftPos >= s->rightPos)
	    die ("positions out of order at line %d in site file %s", vf->line, fileName) ;
	  if (s->leftPos < 0)
	    die ("left position %d at line %d in %s is < 0", s->leftPos, vf->line, fileName) ;
	  if (s->rightPos > cmax)
	    die ("right position %d at line %d in %s is > %d",
		 s->rightPos, vf->line, fileName, cmax) ;
	}
	break ;
      }

  oneFileClose (vf) ;
  return a ;
}

/************************************************************/

SampleSet *samplesRead (char* fileName)
{
  SampleSet *ss = new (1, SampleSet) ;
  ss->names = dictCreate (256) ;
  ss->sample = arrayCreate (256, Sample) ;
  Sample *s ;

  OneFile *vf = oneFileOpenRead (fileName, schema, "smp", 1) ;
  if (!vf) die ("failed to open samples file %s", fileName) ;
  while (oneReadLine (vf))
    switch (vf->lineType)
      {
      case 'N':
	{ int k ;
	  if (!dictAdd (ss->names, oneString(vf), &k))
	    die ("duplicate sample name %s", oneString(vf)) ;
	  s = arrayp (ss->sample, k, Sample) ;
	}
	break ;
      case 'F':
	s->fileName = strdup (oneString(vf)) ;
	break ;
      }
  oneFileClose (vf) ;

  fprintf (stderr, "read %d samples from %s\n", dictMax(ss->names), fileName) ;
  
  return ss ;
}

/************************************************************/

int main (int argc, char *argv[])
{
  --argc ; ++argv ;		/* eat program name */

  outFile = stdout ;
  schema = oneSchemaCreateFromText (schemaText) ;
  
  timeUpdate (stdout) ;		/* initialise timer */
#ifdef OMP
  numThreads = omp_get_max_threads () ;
  omp_set_num_threads (numThreads) ;
#endif

  if (!argc) usage () ;

  while (argc > 3) {
    if (**argv != '-')
      die ("option/command %s does not start with '-': run without arguments for usage", *argv) ;
    //fprintf (stderr, "COMMAND %s", *argv) ;
    //for (i = 1 ; i < argc && *argv[i] != '-' ; ++i) fprintf (stderr, " %s", argv[i]) ;
    //fputc ('\n', stderr) ;
    
#define ARGMATCH(x,y,n)	((!strcmp (*argv, x) || (!strcmp (*argv,y))) && argc >= n && (argc -= n, argv += n))
    if (ARGMATCH("-v","--verbose",1)) isVerbose = !isVerbose ;
    else if (ARGMATCH("-t","--threads",2))
      {
#ifdef OMP
	numThreads = atoi(argv[-1]) ;
	if (numThreads > omp_get_max_threads ()) numThreads = omp_get_max_threads () ;
	omp_set_num_threads (numThreads) ;
#else
	fprintf (stderr, "  can't set thread number - not compiled with OMP\n") ;
#endif
      }
    else if (ARGMATCH("-o","--output",2))
      { if (!strcmp (argv[-1], "-"))
	  outFile = stdout ;
	else if (!(outFile = fopen (argv[-1], "w")))
	  { fprintf (stderr, "can't open output file %s - resetting to stdout\n", argv[-1]) ;
	    outFile = stdout ;
	  }
      }
    else die ("unkown command %s - run without arguments for usage", *argv) ;

    timeUpdate (outFile) ;
  }

  if (argc != 3) die ("missing three file names after options - run without args for usage") ;

  ref = referenceRead (*argv++) ;
  sites = sitesRead (*argv++, ref) ;
  samples = samplesRead (*argv++) ;
  
  fprintf (outFile, "total resources used: ") ; timeTotal (outFile) ;
  if (outFile != stdout) { printf ("total resources used: ") ; timeTotal (stdout) ; }
}

/************* end of file ************/
