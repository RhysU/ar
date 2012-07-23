CXXFLAGS ?= -g -Wall -O3 # -O0 -fno-unsafe-math-optimizations -D_GLIBCXX_DEBUG

all:     zohar example test

zohar:   zohar.cpp   burg.hpp
example: example.cpp burg.hpp
test:    test.cpp    burg.hpp
clean:
	rm -f example zohar test

# Test cases from http://paulbourke.net/miscellaneous/ar/
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
