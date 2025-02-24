SOURCE1=src/deedsBCV0.cu

ifeq ($(SLOW),1)
	OPT =-O
else
	OPT =-O3 -g -Xcompiler -fopenmp -Xcompiler -mavx2 -Xcompiler -msse4.2
endif

.PHONY: target

all: deeds

deeds: $(SOURCE1) Makefile
	nvcc $(SOURCE1) -I src -lz -o deedsBCV -std=c++11 $(OPT)



clean:
	rm -f deedsBCV

