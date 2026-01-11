# GNU-like toolchain assumed

ifeq (icpc,${CXX})
    HOWSTRICT ?= -std=gnu++98 -ansi -Wall
else
    HOWSTRICT ?= -std=c++98 -ansi -pedantic -Wall -Wextra -Wno-long-long
endif
HOWFAST   ?= -g -O2 -DNDEBUG
PRECISION ?= -DREAL=double
CXXFLAGS  ?= $(HOWSTRICT) $(HOWFAST) $(PRECISION)
VERSION   := $(shell git describe --always --dirty 2>/dev/null || echo 'unknown')

all:     zohar example test ar6 arsel faber1986 collomb2009 lorenz

CC = $(CXX) # Force compilation and linking with C++ compiler

zohar.o:   zohar.cpp   ar.hpp
zohar:     zohar.o

example.o: example.cpp ar.hpp
example:   example.o

test.o:    test.cpp    ar.hpp
test:      test.o

ar6.o:   ar6.cpp   ar.hpp
ar6:     ar6.o

arsel.o:   arsel.cpp   ar.hpp
	$(CXX) $(CXXFLAGS) -DARSEL_VERSION='"$(VERSION)"' -DARSEL_CXXFLAGS='"$(CXXFLAGS)"' -c -o arsel.o arsel.cpp

arsel:     arsel.o

faber1986.o:  faber1986.cpp
faber1986:    faber1986.o

collomb2009.o:  collomb2009.cpp
collomb2009:    collomb2009.o

lorenz.o:  lorenz.cpp
lorenz:    lorenz.o

clean:
	rm -f example zohar test ar6 arsel collomb2009 faber1986 lorenz *.o

# Some test cases from http://paulbourke.net/miscellaneous/ar/
check: zohar example test
	@echo
	./zohar
	@echo
	./example
	@echo
	./test test0.coeff test0.dat
	@echo
	./test test1.coeff test1.dat
	@echo
	./test test2.coeff test2.dat
	@echo
	./test test3.coeff test3.dat
	@echo
	./test --subtract-mean rhoe.coeff rhoe.dat
	@echo

# Run quite a bit of random data through the test routines
stress: SHELL=/bin/bash                      # Bash required here
stress: TIME=command time --verbose          # Command timer to use
stress: ORDER=0 0 0 0 0 0 0 0 0 0 0 0 0 0 0  # Fit an AR(15) model
stress: COUNT=5000000                        # Number of samples to generate
stress: RAND=/dev/urandom                    # Random source to use
stress: test
	@printf "Fitting model to %g samples from %s...\n\n" $(COUNT) $(RAND)
	$(TIME) ./test --subtract-mean <(echo $(ORDER)) <(od -tu1 -vAn -N$(COUNT) $(RAND))

###################################################################
# Expose functionality as a Python module called 'ar' when possible
###################################################################
PYTHON ?= $(shell which python)
ifneq "$(PYTHON)" ""

all:   ar.so
ar.so: ar-python.cpp ar.hpp setup.py
	ARSEL_VERSION="$(VERSION)" $(PYTHON) setup.py build_ext --inplace

python-check: ar.so
	$(PYTHON) test_ar.py

check: python-check

clean: python-clean
python-clean:
	rm -rf ar.so ar.cpython*.so build

endif

###################################################################
# Build LaTeX-based write ups as PDFs whenever latexmk is available
###################################################################
LATEXMK ?= $(shell which latexmk)
ifneq "$(LATEXMK)" ""

all:      writeups
WRITEUPS := $(patsubst %.tex,%.pdf,$(wildcard *.tex))
writeups: $(WRITEUPS)

%.pdf : %.tex
	$(LATEXMK)    -dvi- -ps- -pdf $<
	$(LATEXMK) -c -dvi- -ps- -pdf $<

clean: writeups-clean
writeups-clean:
	rm -f $(WRITEUPS)

endif
