# GNU-compatible toolchain assumed

HOWSTRICT ?= -std=c++98 -ansi -pedantic -Wall -Wextra -Werror -Wfatal-errors -Wno-long-long
HOWFAST   ?= -g -O3 -DNDEBUG -funsafe-math-optimizations # Not finite-math-only!
# HOWFAST ?= -g -O0 -fno-unsafe-math-optimizations -D_GLIBCXX_DEBUG
CXXFLAGS  ?= $(HOWSTRICT) $(HOWFAST)

all:     zohar example test arsel

CC = $(CXX) # Force compilation and linking with C++ compiler

zohar.o:   zohar.cpp   ar.hpp
zohar:     zohar.o

example.o: example.cpp ar.hpp
example:   example.o

test.o:    test.cpp    ar.hpp
test:      test.o

arsel.o:   arsel.cpp   ar.hpp
arsel:     arsel.o

clean:
	rm -f example zohar test arsel *.o

# Some test cases from http://paulbourke.net/miscellaneous/ar/
check: all
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
	@$(TIME) ./test --subtract-mean <(echo $(ORDER)) <(od -tu1 -vAn -N$(COUNT) $(RAND))

# Expose functionality through GNU Octave when mkoctfile available
MKOCTFILE ?= $(shell which mkoctfile)
ifneq "$(MKOCTFILE)" ""

all:      octfiles
OCTFILES := $(patsubst %-octfile.cpp,%.oct,$(wildcard *-octfile.cpp))
octfiles: $(OCTFILES)

%.oct : %-octfile.cpp ar.hpp
	env "CXXFLAGS=$(HOWFAST)" $(MKOCTFILE) -v $< -o $@

clean: octfiles-clean
octfiles-clean:
	rm -f $(OCTFILES)

endif

# Expose functionality as a Python module called 'ar' when possible
PYTHON ?= $(shell which python)
ifneq "$(PYTHON)" ""

all:   ar.so
ar.so: ar-python.cpp ar.hpp setup.py
	$(PYTHON) setup.py build_ext --inplace --build-temp python-build

clean: python-clean
python-clean:
	rm -rf ar.so python-build

endif
