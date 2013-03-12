CC= llvm-gcc
RM= rm -f
LD= $(CC)

OPTFLAGS= -fopenmp -O3 -Xarch_x86_64 -march=core2 -Xarch_i386 -march=prescott
DFLAGS=

ARCHLIBS= -framework Accelerate
ARCHFLAGS= -D_MACOSX -arch x86_64 -arch i386

CFLAGS= $(OPTFLAGS) $(ARCHFLAGS) $(DFLAGS) -I/usr/local/include
LFLAGS= $(OPTFLAGS) $(ARCHFLAGS) $(DFLAGS) -L/usr/local/lib

LIBS= -lgsl -lfftw3_threads -lfftw3

OBJS= config.o fsht.o init.o scatmat.o farfield.o spbessel.o \
      shrotate.o spreflect.o translator.o util.o

FASTSPHERE= fastsphere
SPHEREPIX= spherepix

fastsphere: $(OBJS) fastsphere.o spherepix.o
	$(LD) $(LFLAGS) -o $(FASTSPHERE) fastsphere.o $(OBJS) $(LIBS) $(ARCHLIBS)
	$(LD) $(LFLAGS) -o $(SPHEREPIX) spherepix.o $(OBJS) $(LIBS) $(ARCHLIBS)
	@echo "Finished building $(FASTSPHERE) and $(SPHEREPIX)"

ultra: ARCHLIBS= -L$(ATLAS_DIR)/lib -L$(GSL_DIR)/lib -L$(FFTW_DIR)/lib \
	-llapack -lptf77blas -lptcblas -latlas -lgfortran
ultra: ARCHFLAGS= -D_ATLAS -I$(ATLAS_DIR)/include \
	-I$(FFTW_DIR)/include -I$(GSL_DIR)/include
ultra: OPTFLAGS= -fopenmp -O2 -march=native -mtune=native
ultra: CC= gcc
ultra: LD= gfortran
ultra: fastsphere

darwin32: ARCHFLAGS= -D_MACOSX -arch i386
darwin32: fastsphere

clean:
	$(RM) $(FASTSPHERE) $(SPHEREPIX) $(OBJS) fastsphere.o spherepix.o *.core core

.SUFFIXES: .o .c

.c.o:
	$(CC) $(CFLAGS) -o $@ -c $<
