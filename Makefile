HOWSTRICT ?= -std=c++98 -ansi -pedantic -Wall -Wextra -Werror -Wfatal-errors
CXXFLAGS  ?= $(HOWSTRICT) -g -O3 -DNDEBUG
#  CXXFLAGS ?= $(HOWSTRICT) -g -O0 -fno-unsafe-math-optimizations -D_GLIBCXX_DEBUG

all:     zohar example test arsel

zohar:   zohar.cpp   ar.hpp
example: example.cpp ar.hpp
test:    test.cpp    ar.hpp
arsel:   arsel.cpp   ar.hpp

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
stress: COUNT=1000000                        # Number of samples to generate
stress: RAND=/dev/urandom                    # Random source to use
stress: test
	@printf "Fitting model to %g samples from %s...\n\n" $(COUNT) $(RAND)
	@$(TIME) ./test --subtract-mean <(echo $(ORDER)) <(od -tu1 -vAn -N$(COUNT) $(RAND))

# Expose functionality through GNU Octave when mkoctfile available
MKOCTFILE ?= $(shell which mkoctfile)
ifdef MKOCTFILE

all: octfiles

octfiles: arsel.oct arcov.oct

octfiles-clean:
	rm -f *.oct

clean: octfiles-clean

%.oct : %-octfile.cpp ar.hpp
	$(MKOCTFILE) -g $< -o $@

endif
