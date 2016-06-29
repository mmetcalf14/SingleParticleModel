ifeq "$(ARCH)" ""
ARCH = x86_64
endif

ifeq "$(OS)" ""
OS = $(shell uname -s)
endif

ifeq "$(NODE)" ""
NODE = $(shell uname -n)
endif

EIGENINC = /Volumes/Files/GitHub/eigen

LAPACK = -lblas -llapack -lm
LAPACK_OMP = $(LAPACK)
CC = clang++ -O3 -m64 -std=c++11 -stdlib=libc++ -I./ -I$(EIGENINC)

LIB := CorrelationMatrix.o SingleParticle_Hamiltonian_MLM.o

.PHONY: all clean

all: sp.app

sp.app: sp.o $(LIB)
	$(CC) -o $@ $< $(LIB)

sp.o: main.cpp
	$(CC) -c -o $@ $<

CorrelationMatrix.o: CorrelationMatrix.cpp
	$(CC) -c -o $@ $<

SingleParticle_Hamiltonian_MLM.o: SingleParticle_Hamiltonian_MLM.cpp
	$(CC) -c -o $@ $<

clean:
	@rm -rf *.o
