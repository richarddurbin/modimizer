# makefile for modimizer, developed on Richard's Mac.

#CFLAGS= -O3
CFLAGS= -g				# for debugging
#CFLAGS= -03 -DOMP -fopenmp		# for OMP parallelisation - doesn't compile on Mac

all: modmap modasm modutils composition # mod10x fq2b 

clean:
	$(RM) *.o *~ modmap modasm modutils composition # mod10x fq2b 

### object files

UTILS_OBJS=hash.o dict.o array.o utils.o
UTILS_HEADERS=utils.h array.h dict.h hash.h
$(UTILS_OBJS): utils.h $(UTILS_HEADERS)

seqhash.o: seqhash.h

modset.c: modset.h

seqio.o: seqio.h

### programs

modmap: modmap.c seqio.o seqhash.o modset.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ -lm -lz

modasm: modasm.c seqio.o seqhash.o modset.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ -lm -lz

modutils: modutils.c seqio.o seqhash.o modset.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ -lm -lz

composition: composition.c seqio.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ -lz

#fq2b: fq2b.c $(UTILS_OBJS)
#	$(CC) $(CFLAGS) $^ -o $@ -lz

#mod10x: mod10x.c seqio.o seqhash.o $(UTILS_OBJS)
#	$(CC) $(CFLAGS) $^ -o $@ -lm

### end of file
