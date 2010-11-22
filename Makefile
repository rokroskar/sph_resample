#
# Makefile for SKID (Spline Kernel Interpolative Denmax).
#
# on an SGI R4400 add -mips2 to CFLAGS
# on an SGI R8000 (TFP) add -mips4 -O3 to CFLAGS
#
CC=gcc
#CFLAGS	=	-O2
#Following added for gcc
#CFLAGS = -O3 -funroll-loops
CFLAGS = -g 
LIBS	=   -lm

default: sph_rsmpl 

clean:
	rm -f *.o

sph_rsmpl: main.o kd.o smooth1.o 
	$(CC) $(CFLAGS) -o sph_rsmpl main.o kd.o smooth1.o $(LIBS)

main.o: main.c kd.h smooth1.h 

kd.o: kd.c kd.h tipsydefs.h 

smooth1.o: smooth1.c kd.h smooth1.h


