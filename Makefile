CC= gcc
RM= rm -f
LD= gfortran

OPTFLAGS= -fopenmp -O2 -march=opteron -mtune=opteron

CFLAGS= $(OPTFLAGS) -I/opt/local/include -I/usr/local/include
LFLAGS= $(OPTFLAGS) -L/opt/local/lib -L/usr/local/lib -L../spherepack31 -L../gmres

LIBS= -lgmres -lspherepack -lgsl -lfftw3_threads -lfftw3
ARCHLIBS= -alapack_r -lptf77blas -lptcblas -latlas_r

OBJS= config.o fastsphere.o fsht.o init.o scatmat.o farfield.o spbessel.o \
      shrotate.o shtranslate.o spreflect.o translator.o util.o

FASTSPHERE= fastsphere

fastsphere: $(OBJS)
	$(LD) $(LFLAGS) -o $(FASTSPHERE) $(OBJS) $(LIBS) $(ARCHLIBS)

darwin: OPTFLAGS= -fopenmp -O3 -march=nocona -mtune=nocona -arch x86_64 -arch i386
darwin: ARCHLIBS= -framework Accelerate
darwin: $(OBJS)
	@echo "Building universal binary on Darwin."
	$(LD) $(LFLAGS) -o $(FASTSPHERE) $(OBJS) $(LIBS) $(ARCHLIBS)

clean:
	$(RM) $(FASTSPHERE) $(OBJS) *.core core

.SUFFIXES: .o .c

.c.o:
	$(CC) $(CFLAGS) -o $@ -c $<
