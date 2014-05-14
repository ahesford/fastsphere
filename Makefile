CC= /usr/local/bin/clang
RM= rm -f
LD= $(CC)

OPTFLAGS= -fopenmp -O3 -march=core2
DFLAGS=

ARCHLIBS= -framework Accelerate
ARCHFLAGS= -D_MACOSX -flax-vector-conversions

CFLAGS= $(OPTFLAGS) $(ARCHFLAGS) $(DFLAGS) -I/usr/local/include
LFLAGS= $(OPTFLAGS) $(ARCHFLAGS) $(DFLAGS) -L/usr/local/lib

LIBS= -lgsl -lfftw3_threads -lfftw3

OBJS= config.o fsht.o init.o scatmat.o farfield.o ptsrc.o \
      spbessel.o spreflect.o translator.o util.o

FASTSPHERE= fastsphere
SPHEREPIX= spherepix

fastsphere: $(OBJS) fastsphere.o spherepix.o
	$(LD) $(LFLAGS) -o $(FASTSPHERE) fastsphere.o $(OBJS) $(LIBS) $(ARCHLIBS)
	$(LD) $(LFLAGS) -o $(SPHEREPIX) spherepix.o $(OBJS) $(LIBS) $(ARCHLIBS)
	@echo "Finished building $(FASTSPHERE) and $(SPHEREPIX)"

habis: ARCHLIBS= -L$(ATLAS_DIR)/lib -L$(GSL_DIR)/lib -L$(FFTW_DIR)/lib \
	-llapack -lptf77blas -lptcblas -latlas -lgfortran
habis: ARCHFLAGS= -D_ATLAS -I$(ATLAS_DIR)/include \
	-I$(FFTW_DIR)/include -I$(GSL_DIR)/include
habis: OPTFLAGS= -fopenmp -O2 -march=native -mtune=native
habis: CC= gcc
habis: LD= gfortran
habis: fastsphere

clean:
	$(RM) $(FASTSPHERE) $(SPHEREPIX) $(OBJS) fastsphere.o spherepix.o *.core core

.SUFFIXES: .o .c

.c.o:
	$(CC) $(CFLAGS) -o $@ -c $<
