FF?= gfortran
CC?= gcc
AR?= ar
RM?= rm -f

FFLAGS= -O2 -march=nocona -mtune=nocona
CFLAGS= -O2 -march=nocona -mtune=nocona -I/opt/local/include
#FFLAGS= -O2 -mcpu=G4 -mtune=G4
#CFLAGS= -O2 -mcpu=G4 -mtune=G4 -I/opt/local/include

LFLAGS= -L/opt/local/lib -L../spherepack31 -L../libamos
LIBS= -lspherepack -lgsl -lfftw3 -lamos

SHOBJS= fsht.o shtest.o util.o
TROBJS= fsht.o translator.o trtest.o spbessel.o util.o
OBJS= spreflect.o config.o init.o fastsphere.o

SHTEST= shtest
TRTEST= trtest

default: $(OBJS) shtest trtest

shtest: $(SHOBJS)
	$(FF) $(LFLAGS) -o $(SHTEST) $(SHOBJS) $(LIBS)

trtest: $(TROBJS)
	$(FF) $(LFLAGS) -o $(TRTEST) $(TROBJS) $(LIBS)

clean:
	$(RM) $(SHTEST) $(TRTEST) $(SHOBJS) $(TROBJS) $(OBJS) *.core core

.SUFFIXES: .o .f .c

.f.o:
	$(FF) $(FFLAGS) -o $@ -c $<

.c.o:
	$(CC) $(CFLAGS) -o $@ -c $<
