CC= gcc
RM= rm -f
LD= $(CC)

OPTFLAGS= -fopenmp -O2 -march=opteron -mtune=opteron

CFLAGS= $(OPTFLAGS) $(ARCHFLAGS) -I/opt/local/include -I/usr/local/include
LFLAGS= $(OPTFLAGS) $(ARCHFLAGS) -L/opt/local/lib -L/usr/local/lib

LIBS= -lgsl -lfftw3_threads -lfftw3
ARCHLIBS= -lalapack_r -lptf77blas -lptcblas -latlas_r -lgfortran
ARCHFLAGS= -D_FREEBSD -L/usr/local/lib/gcc45

OBJS= config.o fsht.o init.o scatmat.o farfield.o spbessel.o \
      shrotate.o spreflect.o translator.o util.o

FASTSPHERE= fastsphere
SPHEREPIX= spherepix

fastsphere: $(OBJS) fastsphere.o spherepix.o
	$(LD) $(LFLAGS) -o $(FASTSPHERE) fastsphere.o $(OBJS) $(LIBS) $(ARCHLIBS)
	$(LD) $(LFLAGS) -o $(SPHEREPIX) spherepix.o $(OBJS) $(LIBS) $(ARCHLIBS)

darwin: OPTFLAGS= -fopenmp -O3 -Xarch_x86_64 -march=core2 -Xarch_i386 -march=prescott
darwin: ARCHLIBS= -framework Accelerate
darwin: ARCHFLAGS= -D_MACOSX -arch x86_64 -arch i386
darwin: CC= llvm-gcc
darwin: fastsphere
	@echo "Built universal binary on Darwin."

darwin32: OPTFLAGS= -fopenmp -O3 -march=prescott
darwin32: ARCHLIBS= -framework Accelerate
darwin32: ARCHFLAGS= -D_MACOSX -arch i386
darwin32: CC= llvm-gcc
darwin32: fastsphere
	@echo "Built 32-bit binary on Darwin."

clean:
	$(RM) $(FASTSPHERE) $(SPHEREPIX) $(OBJS) fastsphere.o spherepix.o *.core core

.SUFFIXES: .o .c

.c.o:
	$(CC) $(CFLAGS) -o $@ -c $<
