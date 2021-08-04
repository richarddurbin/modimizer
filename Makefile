# makefile for modimizer, developed on Richard's Mac.

#CFLAGS= -O3
CFLAGS= -g -target arm64-apple-macos11				# for debugging
#CFLAGS= -03 -DOMP -fopenmp		# for OMP parallelisation - doesn't compile on Mac

ALL=modmap modasm modutils composition seqconvert seqhoco modrep modtype # mod10x fq2b 

DESTDIR=~/bin

all: $(ALL) # mod10x fq2b 

install:
	cp $(ALL) $(DESTDIR)

clean:
	$(RM) *.o *~ $(ALL) # mod10x fq2b 
	\rm -r *.dSYM

### object files

UTILS_OBJS=hash.o dict.o array.o utils.o
UTILS_HEADERS=utils.h array.h dict.h hash.h
$(UTILS_OBJS): utils.h $(UTILS_HEADERS)

ONE_DIR = ../vgp-tools/Core
HTS_DIR = $(PWD)/../htslib
SEQIO_OPTS = -DONEIO -DBAMIO -I$(HTS_DIR)/htslib/
SEQIO_LIBS = -L$(ONE_DIR) -lONE -L$(HTS_DIR) -Wl,-rpath $(HTS_DIR) -lhts -lm -lbz2 -llzma -lcurl -lz 
# the "-Wl,-rpath $(HTS_DIR)" incantation is needed for local dynamic linking if htslib is not installed centrally

seqhash.o: seqhash.h

modset.c: modset.h

seqio.o: seqio.c seqio.h 
	$(CC) $(CFLAGS) $(SEQIO_OPTS) -c $^

ONElib.o: ONElib.c ONElib.h 
	$(CC) $(CFLAGS) -c $^

### programs

modmap: modmap.c seqio.o seqhash.o modset.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(SEQIO_LIBS)

modrep: modrep.c seqio.o seqhash.o modset.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(SEQIO_LIBS)

modasm: modasm.c seqio.o seqhash.o modset.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(SEQIO_LIBS)

modutils: modutils.c seqio.o seqhash.o modset.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(SEQIO_LIBS)

composition: composition.c seqio.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(SEQIO_LIBS)

seqconvert: seqconvert.c seqio.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(SEQIO_LIBS)

seqhoco: seqhoco.c seqio.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(SEQIO_LIBS)

modtype: modtype.c seqio.o $(UTILS_OBJS) ONElib.o
	$(CC) $(CFLAGS) $^ -o $@ $(SEQIO_LIBS)

#fq2b: fq2b.c $(UTILS_OBJS)
#	$(CC) $(CFLAGS) $^ -o $@ 

#mod10x: mod10x.c seqio.o seqhash.o $(UTILS_OBJS)
#	$(CC) $(CFLAGS) $^ -o $@ -lm

### end of file
