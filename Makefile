FF= gfortran-mp-4.2
CC= gcc-mp-4.2
AR= ar
RM= rm -f

OPTFLAGS= -fopenmp -O -march=prescott -mtune=prescott
FFLAGS= $(OPTFLAGS)
CFLAGS= $(OPTFLAGS) -I/opt/local/include -I/usr/local/include

LFLAGS= $(OPTFLAGS) -L/opt/local/lib -L/usr/local/lib -L../spherepack31 \
	-L../libamos -L../gmres
LIBS= -lgmres -lamos -lspherepack -lgsl -lfftw3 -framework Accelerate

SHOBJS= fsht.o shtest.o util.o
TROBJS= fsht.o translator.o trtest.o spbessel.o util.o
OBJS= config.o fastsphere.o fsht.o init.o scatmat.o farfield.o \
      spbessel.o spreflect.o translator.o util.o

SHTEST= shtest
TRTEST= trtest
FASTSPHERE= fastsphere

default: fastsphere shtest trtest

fastsphere: $(OBJS)
	$(FF) $(LFLAGS) -o $(FASTSPHERE) $(OBJS) $(LIBS)

shtest: $(SHOBJS)
	$(FF) $(LFLAGS) -o $(SHTEST) $(SHOBJS) $(LIBS)

trtest: $(TROBJS)
	$(FF) $(LFLAGS) -o $(TRTEST) $(TROBJS) $(LIBS)

clean:
	$(RM) $(FASTSPHERE) $(OBJS) *.core core \
		$(SHTEST) $(SHOBJS) $(TRTEST) $(TROBJS) 

.SUFFIXES: .o .f .c

.f.o:
	$(FF) $(FFLAGS) -o $@ -c $<

.c.o:
	$(CC) $(CFLAGS) -o $@ -c $<
