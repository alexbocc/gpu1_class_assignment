/*
 * test2D.cpp
 *
 *  Created on: Jun 14, 2016
 *      Author: Attila Bagoly <battila93@gmail.com>
 */


#include "tensor.hpp"
#include <iostream>
#include <fstream>
#include <iterator>
#include <cmath>
#include <chrono>



using namespace std;


template<typename T> matrix<T, 2> mult2D(const matrix<T,2>& m1, const matrix<T, 2>& m2)
{
	auto M1N1 = m1.getDimension(0),
			 M1N2 = m1.getDimension(1),
			 M2N1 = m2.getDimension(0),
			 M2N2 = m2.getDimension(1);

	if (M1N2!=M2N1){
		cerr << "Matrices not compatibile" << endl;
		return matrix<T,2>();
	}

	matrix<T, 2> M(M1N1, M2N2);

	size_t i, j, k;
	T sum;
	for(i=0;i<M1N1;i++){
		for(j=0;j<M2N2;j++){
			sum = 0;
			for(k=0;k<M1N2;k++){
				sum += m1(i, k)*m2(k, j);
			}
			M(i, j) = sum;
		}
	}
	return M;
}

template<typename T> matrix<T, 2> mult2DBlock(const matrix<T,2>& m1, const matrix<T, 2>& m2, size_t b)
{
	auto M1N1 = m1.getDimension(0),
			 M1N2 = m1.getDimension(1),
			 M2N1 = m2.getDimension(0),
			 M2N2 = m2.getDimension(1);

	if (M1N2!=M2N1){
		cerr << "Matrices not compatibile" << endl;
		return matrix<T,2>();
	}

	matrix<T, 2> M(M1N1, M2N2);

	size_t i, j, k, bi, bj, bk, ii, jj, kk, i0, j0, k0;
	T sum;

	size_t bS1 = M1N1/b,
				  bS = M1N2/b,
			   bS2 = M2N2/b;

	for(bi=0;bi<bS1;bi++){
		i0 = bi*b;
		for(bj=0;bj<bS2;bj++){
			j0 = bj*b;
			for(bk=0;bk<bS;bk++){
				k0 = bk*b;
				for(i=0;i<b;i++){
					ii = i0+i;
					for(j=0;j<b;j++){
						jj = j0+j;
						sum = 0;
						for(k=0;k<b;k++){
							sum += m1(ii, k+k0)*m2(k+k0, jj);
						}
						M(ii, jj) += sum;
					}
				}
			}
		}
	}
	return M;
}

int main(int argc, char** argv)
{
	if (argc!=2){
		cerr << "You have to pass a block size as parameter"<<endl;
		return -1;
	}
	matrix<double, 2> m1 = readMatrix<double,2>("mat2d_1");
	matrix<double, 2> m2 = readMatrix<double, 2>("mat2d_2");

	auto t0 = chrono::high_resolution_clock::now();

	matrix<double, 2> m = mult2D<double>(m1, m2);

	auto t1 = chrono::high_resolution_clock::now();

	matrix<double, 2> mB = mult2DBlock<double>(m1, m2, atoi(argv[1]));

	auto t2 = chrono::high_resolution_clock::now();


	auto N = m.getDimension(0),
			 M = m.getDimension(1);
	float sum = 0;
	for(auto i=0, j=0;i<N;i++){
		for(j=0; j<M;j++){
			sum += fabs(m(i, j)-mB(i, j));
		}
	}
	cout <<  "Sum of abs of difference: " << sum << endl;

	cout << "Non block version time: " << chrono::duration_cast<chrono::milliseconds>(t1-t0).count() << " ms" << endl;
	cout << "Block version time (b="<< argv[1] <<"): " << chrono::duration_cast<chrono::milliseconds>(t2-t1).count() << " ms" << endl;


	return 0;
}
