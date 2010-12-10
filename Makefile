CC= llvm-gcc
RM= rm -f
LD= $(CC)

OPTFLAGS= -fopenmp -O3 -Xarch_x86_64 -march=core2 -Xarch_i386 -march=prescott
DFLAGS=

ARCHLIBS= -framework Accelerate
ARCHFLAGS= -D_MACOSX -arch x86_64 -arch i386

CFLAGS= $(OPTFLAGS) $(ARCHFLAGS) $(DFLAGS) -I/opt/local/include -I/usr/local/include
LFLAGS= $(OPTFLAGS) $(ARCHFLAGS) $(DFLAGS) -L/opt/local/lib -L/usr/local/lib

LIBS= -lgsl -lfftw3_threads -lfftw3

OBJS= config.o fsht.o init.o scatmat.o farfield.o spbessel.o \
      shrotate.o spreflect.o translator.o util.o

FASTSPHERE= fastsphere
SPHEREPIX= spherepix

fastsphere: $(OBJS) fastsphere.o spherepix.o
	$(LD) $(LFLAGS) -o $(FASTSPHERE) fastsphere.o $(OBJS) $(LIBS) $(ARCHLIBS)
	$(LD) $(LFLAGS) -o $(SPHEREPIX) spherepix.o $(OBJS) $(LIBS) $(ARCHLIBS)
	@echo "Finished building $(FASTSPHERE) and $(SPHEREPIX)"

bsd: ARCHLIBS= -lalapack_r -lptf77blas -lptcblas -latlas_r -lgfortran
bsd: ARCHFLAGS= -D_FREEBSD -L/usr/local/lib/gcc45
bsd: OPTFLAGS= -fopenmp -O2 -march=opteron -mtune=opteron
bsd: CC= gcc
bsd: LD= gfortran
bsd: fastsphere

darwin32: ARCHFLAGS= -D_MACOSX -arch i386
darwin32: fastsphere

clean:
	$(RM) $(FASTSPHERE) $(SPHEREPIX) $(OBJS) fastsphere.o spherepix.o *.core core

.SUFFIXES: .o .c

.c.o:
	$(CC) $(CFLAGS) -o $@ -c $<
