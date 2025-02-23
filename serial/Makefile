SOURCE1=src/deedsBCV0.cpp

ifeq ($(SLOW),1)
	OPT =-O
else
	OPT =-O3 -g -mavx2 -msse4.2
endif

.PHONY: target

all: deeds

deeds: $(SOURCE1) Makefile
	g++ $(SOURCE1) -I src -lz -o deedsBCV -std=c++11 $(OPT)

clean:
	rm -f deedsBCV

