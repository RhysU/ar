CXXFLAGS ?= -Wall -O3 # -O0 -fno-unsafe-math-optimizations
all:     example test
example: example.cpp burg.hpp
test:    test.cpp    burg.hpp
clean:
	rm -f example test

# Test cases from http://paulbourke.net/miscellaneous/ar/
check: test
	@echo
	./test test1.coeff     test1.dat
	@echo
	./test test2.coeff     test2.dat
	@echo
	./test test3.coeff     test3.dat
	@echo
