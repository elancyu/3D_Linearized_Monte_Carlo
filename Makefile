LIBS = -lm
CFLAGS = -Wall -O3 -g -I/usr/local/include
LFLAGS = -L/HOME/jyuaq/RTMC/RTMC
LD = cc $(LFLAGS) -fopenmp
CC = gcc -c -fopenmp

OBJS = main.o\
       Builder.o\
       Calculation.o\
	   FreeMem.o\
	   GlobalVariables.o\
	   Rand.o\
	   
	   
main: $(OBJS)
	$(LD) -o fSLMC_serial $(OBJS) $(LIBS)
	 
.c.o:
	$(CC) $(CFLAGS) $<

>PHONY: clean
clean:
	rm ${OBJS} *.dat *.log