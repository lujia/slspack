all: lib exam

include ../Makefile.inc

lib:
	@(cd ../src; make)

exam.o: exam.c ../src/libslspack.a

clean:
	@rm -f *.o *.log exam
