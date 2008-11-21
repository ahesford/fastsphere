FF= gfortran-mp-4.2
CC= gcc-mp-4.2
AR= ar
RM= rm -f

OPTFLAGS= -fopenmp -O -march=prescott -mtune=prescott
FFLAGS= $(OPTFLAGS)
CFLAGS= $(OPTFLAGS) -I/opt/local/include

LFLAGS= $(OPTFLAGS) -L/opt/local/lib -L../spherepack31 -L../gmres
LIBS= -lgmres -lspherepack -lgsl -lfftw3_threads -lfftw3 -framework Accelerate

TROBJS= trtest.o shtranslate.o shrotate.o spbessel.o util.o
OBJS= config.o fastsphere.o fsht.o init.o scatmat.o farfield.o spbessel.o \
      shrotate.o shtranslate.o spreflect.o translator.o util.o

TRTEST= trtest
FASTSPHERE= fastsphere

fastsphere: $(OBJS)
	$(FF) $(LFLAGS) -o $(FASTSPHERE) $(OBJS) $(LIBS)

trtest: $(TROBJS)
	$(FF) $(LFLAGS) -o $(TRTEST) $(TROBJS) $(LIBS)

all: fastsphere trtest

clean:
	$(RM) $(FASTSPHERE) $(OBJS) *.core core \
		$(SHTEST) $(SHOBJS) $(TRTEST) $(TROBJS) 

.SUFFIXES: .o .f .c

.f.o:
	$(FF) $(FFLAGS) -o $@ -c $<

.c.o:
	$(CC) $(CFLAGS) -o $@ -c $<
