FF= /usr/local/bin/gfortran
CC= /usr/local/bin/gcc
AR= ar
RM= rm -f

OPTFLAGS= -fopenmp -O -march=nocona -mtune=nocona -m64
FFLAGS= $(OPTFLAGS)
CFLAGS= $(OPTFLAGS) -I/opt/local/include

LFLAGS= $(OPTFLAGS) -L/opt/local/lib -L../gmres
LIBS= -lgmres -lgsl -lfftw3_threads -lfftw3 -framework Accelerate

OBJS= config.o fastsphere.o fsht.o init.o scatmat.o farfield.o spbessel.o \
      shrotate.o shtranslate.o spreflect.o translator.o util.o

FASTSPHERE= fastsphere

fastsphere: $(OBJS)
	$(FF) $(LFLAGS) -o $(FASTSPHERE) $(OBJS) $(LIBS)

clean:
	$(RM) $(FASTSPHERE) $(OBJS) *.core core $(SHTEST) $(SHOBJS)

.SUFFIXES: .o .f .c

.f.o:
	$(FF) $(FFLAGS) -o $@ -c $<

.c.o:
	$(CC) $(CFLAGS) -o $@ -c $<
