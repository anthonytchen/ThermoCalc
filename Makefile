# C-compiler
CC=g++


# libraries
#LIBSPATH=/home/tc0008/sundials/lib
#LIBS=-L/usr/local/lib -lgsl -lgslcblas -lm -pthread -lsundials_cvode -lsundials_nvecserial
#LIBS=-L/usr/local/lib -lm -pthread -lsundials_cvode -lsundials_nvecserial
LIBS=-lgsl
#INC=-I/usr/local/include/

#CFLAGS=-O2 -ffast-math -Wno-deprecated -c
#CFLAGS=-g -O0 -c -Wno-deprecated
CFLAGS=-g -O0 -c -Wno-write-strings

# linker
LD=g++

#LFLAGS=-O2 -ffast-math -Wno-deprecated
#LFLAGS=-g -O0 -Wno-deprecated
LFLAGS=-g -O0 -Wno-write-strings


all:	except.o arg.o cal_pH

clean:
	rm -f *.o *.lnk

### executables ##
cal_pH: cal_pH.cpp arg.o except.o
	$(LD) $(LFLAGS)  cal_pH.cpp arg.o except.o $(LIBS) $(INC) -o cal_pH


### objective files ##
except.o: except.cpp except.h
	$(CC) $(CFLAGS) except.cpp $(INC) -o except.o
arg.o: arg.h arg.c
	$(CC) $(CFLAGS) arg.c $(INC) -o arg.o
