# tensor.hpp
	I created this to have a nice interface for multidimensional matrices:
	when I have to access the element of a 3D matrix I want just to have to write mat(i, j, k).

# genMatrix.cpp
	I use this for creating random matrices.
	To compile: `clang++ -std=c++14 genMatrix.cpp -o genMatrix`.
	To create a 2048x2048 random matrices you have to run: `genMatrix 2048x2048 > outfile`

# test2D.cpp
	I implemented here 2D*2D matrix block multiplication.
	I just learnt here the block multiplication.

# main.cpp
	This is my solution for assignment number 1.
