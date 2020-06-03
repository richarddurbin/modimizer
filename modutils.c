/*  File: modutils.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2018
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Feb 17 10:13 2019 (rd109)
 * Created: Wed Nov 14 22:31:47 2018 (rd109)
 *-------------------------------------------------------------------
 */

#include "modset.h"
#include "seqio.h"

FILE *outFile ;
BOOL isVerbose = FALSE ;

static int addSequence (Modset *ms, char *s, int len) /* return number of hashes */
{
  int nHash = 0 ;
  SeqhashRCiterator *mi = modRCiterator (ms->hasher, s, len) ;
  U64 hash ; int pos ;
  while (modRCnext (mi, &hash, &pos, 0))
    { U32 index = modsetIndexFind (ms, hash, TRUE) ;
      U16 *di = &ms->depth[index] ; ++*di ; if (!*di) *di = U16MAX ;
      ++nHash ;
    }
  return nHash ;
}

static BOOL addSequenceFile (Modset *ms, char *filename, BOOL is10x)
{
  char *seq ;			/* ignore the name for now */
  int len ;
  U64 nSeq = 0, totLen = 0, totHash = 0 ;

  dna2indexConv['N'] = dna2indexConv['n'] = 0 ; /* to get 2-bit encoding */
  SeqIO *si = seqIOopenRead (filename, dna2indexConv, FALSE) ; /* false for no qualities */
  if (!si) return FALSE ;
  while (seqIOread (si))
    { ++nSeq ; totLen += si->seqLen ;
      if (is10x && (nSeq & 0x1)) totHash += addSequence (ms, sqioSeq(si)+23, si->seqLen-23) ;
      else totHash += addSequence (ms, sqioSeq(si), si->seqLen) ;
    }
  seqIOclose (si) ;
  fprintf (outFile, "added %llu sequences total length %llu total hashes %llu, new max %u\n",
	   nSeq, totLen, totHash, ms->max) ;
  return TRUE ;
}

void depthHistogram (Modset *ms, FILE *f)
{
  Array h = arrayCreate (256, U32) ;
  U32 i ;
  for (i = 1 ; i <= ms->max ; ++i)
    if (ms->depth[i] < arrayMax(h)) ++arr(h,ms->depth[i],U32) ; /* more efficient to check */
    else ++array(h,ms->depth[i],U32) ;
  for (i = 0 ; i < arrayMax(h) ; ++i)
    if (arr(h,i,U32)) fprintf (f, "DP\t%u\t%u\n", i, arr(h,i,U32)) ;
  arrayDestroy (h) ;
}

void reportDepths (Modset *ms, Array ma, FILE *f)
{
  U32 i, j, index ;
  for (i = 1 ; i <= ms->max ; ++i)
    { fprintf (f, "MH\t%llx\t%d\t%u", ms->value[i], msCopy(ms,i), ms->depth[i]) ;
      for (j = 0 ; j < arrayMax(ma) ; ++j)
	if ((index = modsetIndexFind (arr(ma,j,Modset*), ms->value[i], FALSE)))
	  fprintf (f, "\t%u", arr(ma,j,Modset*)->depth[index]) ;
	else
	  fprintf (f, "\t0") ;
      fputc ('\n', f) ;
    }
}

void usage (void)
{ fprintf (stderr, "Usage: modutils <commands>\n") ;
  fprintf (stderr, "Commands are executed in order - set parameters before using them!\n") ;
  fprintf (stderr, "  -v | --verbose : toggle verbose mode\n") ;
  fprintf (stderr, "  -o | --output <output filename> : '-' for stdout\n") ;
  fprintf (stderr, "  -c | --modcreate table_bits{28} kmer{19} mod{31} seed{17}: can truncate parameters\n") ;
  fprintf (stderr, "  -r | --read <mod file>\n") ;
  fprintf (stderr, "  -w | --write <mod file>\n") ;
  fprintf (stderr, "  -a | --add <read file> : add kmers from read file\n") ;
  fprintf (stderr, "  -x | --add10x <10x read file> : add kmers from 10x read file\n") ;
  fprintf (stderr, "  -m | --merge <mod file> : add kmers from read file; writes depths\n") ;
  fprintf (stderr, "  -p | --prune <min> <max> : remove mod entries < min or >= max\n") ;
  fprintf (stderr, "  -s | --setcopy <copy1min> <copy2min> <copyMmin> : reset mod copy\n") ;
  fprintf (stderr, "  -sM | --setcopyM <copyMmin> : set copyM if depth > copyMmin\n") ;
  fprintf (stderr, "  -H | --hist <outfile> : print depth histogram\n") ;
  fprintf (stderr, "  -d | --depth <outfile> <mod file>* : print depth per mod [also in other files]\n") ;
  fprintf (stderr, "command -c or -r must come before other commands from -w onwards\n") ;
  fprintf (stderr, "read files can be fasta or fastq, gzipped or not\n") ;
  fprintf (stderr, "example usage\n") ;
  fprintf (stderr, "  modutils -c 30 19 31 17 -a XR1.fa.gz -a XR2.fa.gz -w X.mod\n") ;
  fprintf (stderr, "  modutils -c 30 19 31 17 -a YR1.fa.gz -a YR2.fa.gz -w Y.mod\n") ;
  fprintf (stderr, "  modutils -r X.mod -m Y.mod -w XY1.mod -H XY.his\n") ;
  fprintf (stderr, "then look at histogram XY.his and decide on thresholds, then\n") ;
  fprintf (stderr, "  modutils -r XY1.mod -p 5 200 -s 10 50 100 -w XY2.mod\n") ;
  fprintf (stderr, "  modutils -r XY2.mod -d XY.depths X.mod Y.mod\n") ;
  fprintf (stderr, "XY.depths will have columns: hash, depth_in_XY2, depth_inX, depth_in_Y\n") ;
}

int main (int argc, char *argv[])
{
  --argc ; ++argv ;		/* eat program name */
  if (!argc) usage () ;

  outFile = stdout ;
  timeUpdate (stdout) ;		/* initialise timer */

  Modset *ms = 0 ;
  int i ;			/* generically useful variables */
  FILE *f ;

  while (argc) {
    if (**argv != '-')
      die ("option/command %s does not start with '-': run without arguments for usage", *argv) ;
    fprintf (stderr, "COMMAND %s", *argv) ;
    for (i = 1 ; i < argc && *argv[i] != '-' ; ++i) fprintf (stderr, " %s", argv[i]) ;
    fputc ('\n', stderr) ;
    
#define ARGMATCH(x,y,n)	((!strcmp (*argv, x) || (!strcmp (*argv,y))) && argc >= n && (argc -= n, argv += n))
    if (ARGMATCH("-v","--verbose",1)) isVerbose = !isVerbose ;
    else if (ARGMATCH("-o","--output",2))
      { if (!strcmp (argv[-1], "-"))
	  outFile = stdout ;
	else if (!(outFile = fopen (argv[-1], "w")))
	  { fprintf (stderr, "can't open output file %s - resetting to stdout\n", argv[-1]) ;
	    outFile = stdout ;
	  }
      }
    else if (!ms && ARGMATCH("-c","--create",1))
      { int B = 28, k = 19, w = 31, s = 17 ;
	if (argc && **argv != '-')
	  { if (!(B = atoi(*argv)) || B < 20 || B > 34) die ("bad modbuild B %s", *argv) ;
	    if (--argc && **++argv != '-')
	      { if (!(k = atoi(*argv)) || k < 1) die ("bad modbuild k %s", *argv) ;
		if (--argc && **++argv != '-')
		  { if (!(w = atoi(*argv)) || k < 1) die ("bad modbuild w %s", *argv) ;
		    if (--argc && **++argv != '-')
		      { if (!(s = atoi(*argv))) die ("bad modbuild w %s", *argv) ;
			--argc ; ++argv ;
		      }
		  }
	      }
	  }
	Seqhash *sh = seqhashCreate (k, w, s) ;
	seqhashReport (sh, outFile) ;
	ms = modsetCreate (sh, B, 0) ;
      }
    else if (!ms && ARGMATCH("-r","--read",2))
      { if (!(f = fopen (argv[-1], "r"))) die ("failed to open mod file %s", argv[-1]) ;
	ms = modsetRead (f) ;
	fclose (f) ;
	modsetSummary (ms, outFile) ;
      }
    else if (ms && ARGMATCH("-w","--write",2))
      { if (!(f = fopen (argv[-1], "w"))) die ("failed to open mod file %s", argv[-1]) ;
	modsetWrite (ms, f) ;
	fclose (f) ;
      }
    else if (ms && ARGMATCH("-p","--prune",3))
      { modsetDepthPrune (ms, atoi(argv[-2]), atoi(argv[-1])) ;
	modsetSummary (ms, outFile) ;
      }
    else if (ms && ARGMATCH("-s","--setcopy",4))
      { int copy1min = atoi(argv[-3]), copy2min = atoi(argv[-2]), copyMmin = atoi(argv[-1]) ;
	U32 u ;
	for (u = 1 ; u <= ms->max ; ++u)
	  if (ms->depth[u] < copy1min) msSetCopy0(ms,u) ;
	  else if (ms->depth[u] < copy2min) msSetCopy1(ms,u) ;
	  else if (ms->depth[u] < copyMmin) msSetCopy2(ms,u) ;
	  else msSetCopyM(ms,u) ;
	modsetSummary (ms, outFile) ;
      }
    else if (ms && ARGMATCH("-sM","--setcopyM",2))
      { int copyMmin = atoi(argv[-1]) ;
	U32 u ; for (u = 1 ; u <= ms->max ; ++u) if (ms->depth[u] >= copyMmin) msSetCopyM(ms,u) ;
	modsetSummary (ms, outFile) ;
      }
    else if (ms && ARGMATCH("-a","--add",2))
      { if (!addSequenceFile (ms, argv[-1], FALSE))
	  die ("failed to open sequence file %s", argv[-1]) ;
	modsetSummary (ms, outFile) ;
      }
    else if (ms && ARGMATCH("-x","--add10x",2))
      { if (!addSequenceFile (ms, argv[-1], TRUE))
	  die ("failed to open sequence file %s", argv[-1]) ;
	modsetSummary (ms, outFile) ;
      }
    else if (ms && ARGMATCH ("-m","--merge",2))
      { if (!(f = fopen (argv[-1], "r"))) die ("failed to open mod file %s", argv[-1]) ;
	Modset *ms2 = modsetRead (f) ;
	modsetSummary (ms2, outFile) ;
	fclose (f) ;
	if (!modsetMerge (ms, ms2))
	  fprintf (stderr, "modset %s incompatible with current - unable to merge\n", argv[-1]) ;
	modsetDestroy (ms2) ;
	modsetSummary (ms, outFile) ;
      }
    else if (ms && ARGMATCH ("-H","--hist",2))
      { if (!(f = fopen (argv[-1], "w"))) die ("failed to open histogram file %s", argv[-1]) ;
       depthHistogram (ms, f) ;
       fclose (f) ;
      }
    else if (ms && ARGMATCH ("-d","--depths",2))
      { FILE *fd ;
	if (!(fd = fopen (argv[-1], "w"))) die ("failed to open depths file %s", argv[-1]) ;
	Array ma = arrayCreate (32, Modset*) ;
	while (argc && **argv != '-')
	  { if (!(f = fopen (*argv, "r"))) die ("failed to open mod file %s", *argv) ;
	    array(ma,arrayMax(ma),Modset*) = modsetRead (f) ;
	    modsetSummary (arr(ma,arrayMax(ma)-1,Modset*), outFile) ;
	    fclose (f) ;
	    --argc ; ++argv ;
	  }
	reportDepths (ms, ma, fd) ;
	arrayDestroy (ma) ;
	fclose (fd) ;
      }
    else die ("unknown command %s - run without arguments for usage", *argv) ;

    timeUpdate (outFile) ;
  }

  fprintf (outFile, "total resources used: ") ; timeTotal (outFile) ;
  if (outFile != stdout) { printf ("total resources used: ") ; timeTotal (stdout) ; }
}

/************* end of file ************/
