/*
 * genMatrix.cpp
 *
 *  Created on: Jun 14, 2016
 *      Author: Attila Bagoly <battila93@gmail.com>
 */

#include "tensor.hpp"
#include <iostream>
#include <random>

using namespace std;

int main(int argc, char**argv)
{
  if (argc<2){
    cerr << "You must pass dimension (eg. 2x2, 3x3, 3x3x3)" << endl;
    return 1;
  }
  string tmp(argv[1]);
  stringstream stream(tmp);
  vector<int> dim;
  while(getline(stream, tmp, 'x')){
    dim.push_back(stoi(tmp));
  }
  if (dim.size()!=2 && dim.size()!=3){
    cerr << "Just 2D and 3D tensors are supported!" << endl;
    return 1;
  }
  std::random_device engine;
  std::uniform_real_distribution<float> uniform(-100, 100);
  if (dim.size()==2){
    matrix<float, 2> mat(dim[0], dim[1]);
    for(std::size_t i=0, j;i<dim[0];i++){
      for(j=0;j<dim[1];j++){
        mat(i, j) = uniform(engine);
      }
    }
    cout << mat;
  }
  if (dim.size()==3){
    matrix<float, 3> mat(dim[0], dim[1], dim[2]);
    for(std::size_t i=0, j, k;i<dim[0];i++){
      for(j=0;j<dim[1];j++){
        for(k=0;k<dim[2];k++){
          mat(i, j, k) = uniform(engine);
        }
      }
    }
    cout << mat;
  }

  return 0;
}
