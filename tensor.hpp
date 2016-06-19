/*
 * tensor.hpp
 *
 *  Created on: Jun 14, 2016
 *      Author: Attila Bagoly <battila93@gmail.com>
 */

#ifndef __TENSOR_HPP__
#define __TENSOR_HPP__

#include <vector>
#include <cstddef>
#include <type_traits>
#include <algorithm>
#include <ostream>
#include <istream>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>



template<typename T> std::enable_if_t<std::is_same<int, T>::value, T> to_number(const std::string& s){
  return std::stoi(s);
}
template<typename T> std::enable_if_t<std::is_same<float, T>::value, T> to_number(const std::string& s){
  return std::stof(s);
}
template<typename T> std::enable_if_t<std::is_same<double, T>::value, T> to_number(const std::string& s){
  return std::stod(s);
}


template<typename T, std::size_t DIMENSION> class matrix {
private:
  std::vector<std::size_t> dim;
  std::vector<T> mat;

public:
  matrix(){}
  template<typename... Args> matrix(Args... dims) : dim({static_cast<std::size_t>(dims)...}){
    static_assert(sizeof...(dims)==DIMENSION, "Number of arguments must be equal to dimension");
    mat.resize( std::accumulate(dim.begin(), dim.end(), 1, std::multiplies<std::size_t>()) );
  }
  matrix(const matrix& src) = default;
  matrix(matrix&& src): dim(std::move(src.dim)), mat(std::move(src.mat)){};

  matrix& operator=(const matrix& m) = default;
  matrix<T, DIMENSION>& operator=(matrix&& src){
    mat = std::move(src.mat);
    dim = std::move(src.dim);
  }

  template<typename... Args> void resize(Args ... dims){
    static_assert(sizeof...(dims)==DIMENSION, "Number of arguments must be equal to dimension");
    size_t dim_idx = 0;
    dim.resize(DIMENSION);
    for(auto& d: {static_cast<std::size_t>(dims)...}){
      dim[dim_idx++] = d;
    }
    mat.resize( std::accumulate(dim.begin(), dim.end(), 1, std::multiplies<std::size_t>()) );
  }

  std::size_t getDimension(int i) const{
    return dim[i];
  }

  void makeZeroMatrix(){
    std::fill(mat.begin(), mat.end(), 0);
  }


  // 2D matrix
  T& operator()(std::size_t i1, std::size_t i2){
    static_assert(DIMENSION==2, "Not 2D matrix!");
    return mat[i1+dim[0]*i2];
  }
  const T& operator()(std::size_t i1, std::size_t i2) const{
    static_assert(DIMENSION==2, "Not 2D matrix!");
    return mat[i1+dim[0]*i2];
  }


  // 3D matrix
  T& operator()(std::size_t i1, std::size_t i2, std::size_t i3){
    static_assert(DIMENSION==3, "Not 3D matrix!");
    return mat[i1+dim[0]*i2+dim[1]*dim[0]*i3];
  }
  const T& operator()(std::size_t i1, std::size_t i2, std::size_t i3) const{
    static_assert(DIMENSION==3, "Not 3D matrix!");
    return mat[i1+dim[0]*i2+dim[1]*dim[0]*i3];
  }

  /**
    * output operaotrs
  **/
  template<typename TT=T> friend std::ostream& operator<<(std::ostream& o, typename std::enable_if_t<DIMENSION==2, const matrix<TT, 2>>& t){
    o << t.dim[0] << " " << t.dim[1] << std::endl;
    for(size_t i=0, j;i<t.dim[0];i++){
      for(j=0;j<t.dim[1];j++){
        o << t(i, j) <<",";
      }
      o<<std::endl;
    }
    return o;
  }

  template<typename TT=T>friend std::ostream& operator<<(std::ostream& o, typename std::enable_if_t<DIMENSION==3, const matrix<TT, 3>>& t){
    o << t.dim[0] << " " << t.dim[1] << " " << t.dim[2] << std::endl;
    for(size_t i=0, j, k;i<t.dim[0];i++){
      for(j=0;j<t.dim[1];j++){
        for(k=0;k<t.dim[2];k++){
            o << t(i, j, k) <<",";
        }
        o << ";";
      }
      o<<std::endl;
    }
    return o;
  }

  /**
   * input operators
  **/

  template<typename TT=T> friend std::istream& operator>>(std::istream& in, typename std::enable_if_t<DIMENSION==2, matrix<TT,2>>& mat){
    std::string tmp;
    std::size_t row = 0,
                col = 0;
    while(std::getline(in, tmp)){
      if (tmp.size()<=0) continue;
      if (row>=mat.dim[0]){
        std::cerr << "more row in file then what was set up in first line!" << std::endl;
        return in;
      }
      col = 0;
      std::stringstream sstream(tmp);
      while(std::getline(sstream, tmp, ',')){
        if (col>=mat.dim[1]){
          std::cerr << "more element in row then what was set up in first line!" << std::endl;
          return in;
        }
        mat(row, col++) = to_number<T>(tmp);
      }
      row++;
    }
    return in;
  }

  template<typename TT=T> friend std::istream& operator>>(std::istream& in, typename std::enable_if_t<DIMENSION==3, matrix<TT,3>>& mat){
    std::string tmp;
    std::size_t d1, d2, d3;
    d1 = 0;
    while(std::getline(in, tmp)){
      if (tmp.size()<=0) continue;
      if (d1>=mat.dim[0]){
        std::cerr << "more row in file then what was set up in first line!" << std::endl;
        return in;
      }
      d2 = 0;
      std::stringstream sstream(tmp);
      while(std::getline(sstream, tmp, ';')){
        if (d2>=mat.dim[1]){
          std::cerr << "more element in row then what was set up in first line!" << std::endl;
          return in;
        }
        d3 = 0;
        std::stringstream sstream2(tmp);
        while(std::getline(sstream2, tmp, ',')){
          if (d3>=mat.dim[2]){
            std::cerr << "more element in z direction then what was set up in first line!" << std::endl;
            return in;
          }
          mat(d1,d2, d3++) = to_number<T>(tmp);
        }
        d2++;
      }
      d1++;
    }
    return in;
  }
};


/**
 * Read matrix from file.
 */

template<typename T, std::size_t D> std::enable_if_t<D==2, matrix<T,D>> readMatrix(const std::string& file_name){
	matrix<T, D> mat;
	std::ifstream in(file_name.c_str(), std::ios::in);
	if (!in.is_open()){
		std::cerr << "Wasn't able to read " << file_name << " file!" << std::endl;
		return mat;
	}
  std::size_t d1, d2;
  in >> d1 >> d2;
  mat.resize(d1, d2);
	in >> mat;
	return mat;
}

template<typename T, std::size_t D> std::enable_if_t<D==3, matrix<T,D>> readMatrix(const std::string& file_name){
	matrix<T, D> mat;
	std::ifstream in(file_name.c_str(), std::ios::in);
	if (!in.is_open()){
		std::cerr << "Wasn't able to read " << file_name << " file!" << std::endl;
		return mat;
	}
  std::size_t d1, d2, d3;
  in >> d1 >> d2 >> d3;
  mat.resize(d1, d2, d3);
	in >> mat;
	return mat;
}


#endif
