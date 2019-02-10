
.SUFFIXES: .cxx

default: lib
all: lib

include ../Makefile.inc
include Makefile.dep

SRC_CXX   = $(wildcard *.cxx)
DIR_CXX   = $(notdir $(SRC_CXX))
OBJS_CXX  = $(patsubst %.cxx,%.o,$(DIR_CXX))

SRC_C   = $(wildcard *.c)
DIR_C   = $(notdir $(SRC_C))
OBJS_C  = $(patsubst %.c,%.o,$(DIR_C))

OBJS     := $(OBJS_CXX) $(OBJS_C)

lib: libslspack.a

libslspack.a: ${OBJS}
	${AR} -ru $@ $(OBJS)
	${RANLIB} $@

.PHONY: clean

clean:
	@rm -fv libslspack.a
	@rm -fv *.o 
