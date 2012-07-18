CXXFLAGS ?= -g -Wall -O3 # -O0 -fno-unsafe-math-optimizations -D_GLIBCXX_DEBUG

all:     example zohar test

example: example.cpp burg.hpp
test:    test.cpp    burg.hpp
zohar:   zohar.cpp   burg.hpp
clean:
	rm -f example zohar test

# Test cases from http://paulbourke.net/miscellaneous/ar/
check: all
	@echo
	./example
	@echo
	./zohar
	@echo
	./test test1.coeff     test1.dat
	@echo
	./test test2.coeff     test2.dat
	@echo
	./test test3.coeff     test3.dat
	@echo
