STRICTER ?= -Werror -Wall -Wextra -std=c++98 -ansi -pedantic
CXXFLAGS ?= $(STRICTER) -g -O3 -DNDEBUG
#  CXXFLAGS ?= $(STRICTER) -g -O0 -fno-unsafe-math-optimizations -D_GLIBCXX_DEBUG

all:     zohar example test

zohar:   zohar.cpp   burg.hpp
example: example.cpp burg.hpp
test:    test.cpp    burg.hpp
clean:
	rm -f example zohar test

# Some test cases from http://paulbourke.net/miscellaneous/ar/
check: all
	@echo
	./zohar
	@echo
	./example
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
