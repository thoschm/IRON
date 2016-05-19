all:
	g++ -std=c++11 -O3 -march=native -o examples/demo src/demo.cpp -Iinclude -I/usr/include/eigen3/ -Wall -pedantic -Wno-unused-private-field

