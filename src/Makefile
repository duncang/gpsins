
INCL =-I/usr/realtime/include -I/usr/src/linux/include -I/opt/local/include -I/usr/local/include -I../../novatelstore/src
CFLAGS =-Wall -Wstrict-prototypes -pipe -g
LIBS = -L/usr/realtime/lib -L/opt/local/lib -L/usr/local/lib -pthread  -lgsl -lgslcblas

#USE_RTAI=1
#LIBS += -llxrt

all: gpsins fixmqueue gps_single

fixmqueue: fixmqueue.c
	gcc -o fixmqueue $(INCL) $(CLFAGS) $(LIBS) fixmqueue.c

gpsins: gpsins.c gpsins.h gps.o matrix.o
	gcc -o gpsins  $(INCL) $(CFLAGS) $(LIBS) gpsins.c gps.o matrix.o

gpsstore: gpsstore.c gpsins.h novatel.h
	gcc -o gpsstore $(INCL) $(CFLAGS) -ansi -pthread -lgsl \
	-lgslcblas gpsstore.c 

gps_single: gps_single.c  gps.o matrix.o 
	gcc -o gps_single $(INCL) $(CFLAGS) -pthread  \
	gps_single.c ../../novatelstore/src/novatel.o gps.o matrix.o \
	$(LIBS)


gps.o: gps.c gps.h
	gcc -o gps.o -c gps.c  $(CFLAGS) $(INCL)

matrix.o: matrix.c matrix.h
	gcc -o matrix.o -c matrix.c $(CFLAGS)

quat.o: quat.c quat.h
	gcc -o quat.o -c quat.c $(CFLAGS)

clean:
	-rm *~ *.o gpsins gpsstore gps_single
	
doc: gpsins.h gpsins.c 
	doxygen

      
