CC= gcc
RM= rm -f
LD= $(CC)

OPTFLAGS= -fopenmp -O2 -march=opteron -mtune=opteron

CFLAGS= $(OPTFLAGS) -I/opt/local/include -I/usr/local/include
LFLAGS= $(OPTFLAGS) -L/opt/local/lib -L/usr/local/lib -L../spherepack31 -L../gmres

LIBS= -lgmres -lspherepack -lgsl -lfftw3_threads -lfftw3 -lgfortran -lm
ARCHLIBS= -alapack_r -lptf77blas -lptcblas -latlas_r

OBJS= config.o fsht.o init.o scatmat.o farfield.o spbessel.o \
      shrotate.o spreflect.o translator.o util.o

FASTSPHERE= fastsphere
SPHEREPIX= spherepix

fastsphere: $(OBJS) fastsphere.o spherepix.o
	$(LD) $(LFLAGS) -o $(FASTSPHERE) fastsphere.o $(OBJS) $(LIBS) $(ARCHLIBS)
	$(LD) $(LFLAGS) -o $(SPHEREPIX) spherepix.o $(OBJS) $(LIBS) $(ARCHLIBS)

darwin: OPTFLAGS= -fopenmp -O3 -march=nocona -mtune=nocona -arch x86_64 -arch i386
darwin: ARCHLIBS= -framework Accelerate
darwin: $(OBJS) fastsphere.o spherepix.o
	@echo "Building universal binary on Darwin."
	$(LD) $(LFLAGS) -o $(FASTSPHERE) fastsphere.o $(OBJS) $(LIBS) $(ARCHLIBS)
	$(LD) $(LFLAGS) -o $(SPHEREPIX) spherepix.o $(OBJS) $(LIBS) $(ARCHLIBS)

darwin32: OPTFLAGS= -fopenmp -O3 -march=nocona -mtune=nocona -arch i386
darwin32: ARCHLIBS= -framework Accelerate
darwin32: $(OBJS) fastsphere.o spherepix.o
	@echo "Building universal binary on Darwin."
	$(LD) $(LFLAGS) -o $(FASTSPHERE) fastsphere.o $(OBJS) $(LIBS) $(ARCHLIBS)
	$(LD) $(LFLAGS) -o $(SPHEREPIX) spherepix.o $(OBJS) $(LIBS) $(ARCHLIBS)

clean:
	$(RM) $(FASTSPHERE) $(SPHEREPIX) $(OBJS) fastsphere.o spherepix.o *.core core

.SUFFIXES: .o .c

.c.o:
	$(CC) $(CFLAGS) -o $@ -c $<
