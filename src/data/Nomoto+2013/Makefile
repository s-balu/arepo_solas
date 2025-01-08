# Makefile for CELib.

.SUFFIX: .s

CC = gcc 
#CFLAGS = -O3 -h list=a #-fPIC
CFLAGS = -std=c99 -O3 -g -Wall -ffast-math  #-fPIC
#CFLAGS = -std=c99 -O3 -Wall -ffast-math  #-fPIC
LDFLAGS =  -lm # -lefence
OBJS = main.c

all:
	$(CC) $(CFLAGS) ./main.c -o ngen.out $(LDFLAGS) 

clean:
	rm *.o *.s *.a *.swp ngen.out -f

.s.o:
	$(CC) $(CFLAGS) -c -g $^ 

.c.s:
	$(CC) $(CFLAGS) -S $^ 

tag:
	ctags *.[c,h] 

main.o: main.c Makefile
