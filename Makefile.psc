FF= ifort
CC= icc
AR= ar
RM= rm -f

MKLPATH=/opt/intel/mkl/10.0.2.018/lib/64
OPTFLAGS= -O3 -openmp
FFLAGS= $(OPTFLAGS)
CFLAGS= $(OPTFLAGS) -I/usr/local/packages/fftw3/include -I../gsl/include

LFLAGS= $(OPTFLAGS) -nofor-main -L/usr/local/packages/fftw3/lib \
	-L../gsl/lib -L../gmres -L$(MKLPATH)
LIBS= -lgmres -lgsl -lgslcblas -lfftw3_threads -lfftw3 \
      -lmkl_lapack64 -lmkl_ipf -lmkl

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
