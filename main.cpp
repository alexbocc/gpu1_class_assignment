/*
 * main.cpp
 *
 *  Created on: Jun 18, 2016
 *      Author: Attila Bagoly <battila93@gmail.com>
 */

#include "tensor.hpp"
#include <chrono>

using namespace std;

template<typename T> matrix<T, 3> multiply(const matrix<T, 3>& A, const matrix<T, 2>& m1, const matrix<T, 2>& m2, const matrix<T, 2>& m3)
{
  auto AN1 = A.getDimension(0),
       AN2 = A.getDimension(1),
       AN3 = A.getDimension(2);

  auto m1N1 = m1.getDimension(0),
       m1N2 = m1.getDimension(1),
       m2N1 = m2.getDimension(0),
       m2N2 = m2.getDimension(1),
       m3N1 = m3.getDimension(0),
       m3N2 = m3.getDimension(1);

  matrix<T, 3> M(m1N1, m2N1, m3N1);

  if (AN1!=m1N2 || AN2!=m2N2 || AN3!=m3N2){
    cerr << "Dimension mismatch! " << endl;
    return M;
  }

  size_t i, j, k, x, y, z;
  T sum;

  for(i=0;i<m1N1;i++){
    for(j=0;j<m2N1;j++){
      for(k=0;k<m3N1;k++){
        sum = 0;
        for(x=0;x<AN1;x++){
          for(y=0;y<AN2;y++){
            for(z=0;z<AN3;z++){
              sum += A(x, y, z)*m1(i, x)*m2(j, y)*m3(k, z);
            }
          }
        }
        M(i, j, k) = sum;
      }
    }
  }
  return M;
}

template<typename T> matrix<T, 3> multiplyBlock0(const matrix<T, 3>& A, const matrix<T, 2>& m1, const matrix<T, 2>& m2, const matrix<T, 2>& m3, unsigned int b)
{
  auto AN1 = A.getDimension(0),
       AN2 = A.getDimension(1),
       AN3 = A.getDimension(2);

  auto m1N1 = m1.getDimension(0),
       m1N2 = m1.getDimension(1),
       m2N1 = m2.getDimension(0),
       m2N2 = m2.getDimension(1),
       m3N1 = m3.getDimension(0),
       m3N2 = m3.getDimension(1);

  matrix<T, 3> M(m1N1, m2N1, m3N1);

  if (AN1!=m1N2 || AN2!=m2N2 || AN3!=m3N2){
    cerr << "Dimension mismatch! " << endl;
    return M;
  }

  T sum;
  size_t bi, bj, bk, bm, bn, bl, i, j, k, x, y, z, i0, j0, k0, m0, n0, l0, ii, jj, kk;
  size_t bS1 = m1N1/b,
         bS2 = m2N1/b,
         bS3 = m3N1/b,
         bS4 = AN1/b,
         bS5 = AN2/b,
         bS6 = AN3/b;

  for(bi=0;bi<bS1;++bi){
    i0 = bi * b;
    for(bj=0;bj<bS2;++bj){
      j0 = bj * b;
      for(bk=0;bk<bS3;++bk){
        k0 = bk * b;
        for(bm=0;bm<bS4;++bm){
          m0 = bm * b;
          for(bn=0;bn<bS5;++bn){
            n0 = bn * b;
            for(bl=0;bl<bS6;++bl){
              l0 = bl * b;
              for(i=0;i<b;++i){
                ii = i0 + i;
                for(j=0;j<b;++j){
                  jj = j0 + j;
                  for(k=0;k<b;++k){
                    kk = k0 + k;
                    sum = 0;
                    for(x=0;x<b;++x){
                      for(y=0;y<b;++y){
                        for(z=0;z<b;++z){
                          sum += A(x+m0, y+n0, z+l0)*m1(ii, x+m0)*m2(jj, y+n0)*m3(kk, z+l0);
                        }
                      }
                    }
                    M(ii, jj, kk) += sum;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return M;
}



template<typename T> matrix<T, 3> multiplyBlock(const matrix<T, 3>& A, const matrix<T, 2>& m1, const matrix<T, 2>& m2, const matrix<T, 2>& m3, unsigned int b)
{
  auto AN1 = A.getDimension(0),
       AN2 = A.getDimension(1),
       AN3 = A.getDimension(2);

  auto m1N1 = m1.getDimension(0),
       m1N2 = m1.getDimension(1),
       m2N1 = m2.getDimension(0),
       m2N2 = m2.getDimension(1),
       m3N1 = m3.getDimension(0),
       m3N2 = m3.getDimension(1);

  matrix<T, 3> M;

  if (AN1!=m1N2 || AN2!=m2N2 || AN3!=m3N2){
    cerr << "Dimension mismatch! " << endl;
    return M;
  }

  T sum;

  size_t bi, bj, bk, bm, i, j, k, i0, j0, k0, m0, ii, jj, kk, x, y, z, xx, yy, zz, s;

  M.resize(m1N1,AN2, AN3);

  size_t bS1 = m1N1/b,
         bS2 = AN2/b,
         bS3 = AN3/b,
         bS  = AN1/b;
  for(bi=0;bi<bS1;++bi){
    i0 = bi * b;
    for(bj=0;bj<bS2;++bj){
      j0 = bj * b;
      for(bk=0;bk<bS3;++bk){
        k0 = bk * b;
        for(bm=0;bm<bS;++bm){
          m0  = bm * b;
          for(i=0;i<b;++i){
            ii = i0 + i;
            for(y=0;y<b;++y){
              yy = j0 + y;
              for(z=0;z<b;++z){
                zz = k0 + z;

                sum = 0;
                for(s=0;s<b;++s){
                  sum += A(s+m0, yy, zz)*m1(ii, s+m0);
                }
                M(ii, yy, zz) += sum;

              }
            }
          }
        }
      }
    }
  }

  matrix<T,3> A1(m1N1, m2N1, AN3);
  bS1 = m1N1/b;
  bS2 = m2N1/b;
  bS3 = AN3/b;
  bS  = AN2/b;
  for(bi=0;bi<bS1;++bi){
    i0 = bi * b;
    for(bj=0;bj<bS2;++bj){
      j0 = bj * b;
      for(bk=0;bk<bS3;++bk){
        k0 = bk * b;
        for(bm=0;bm<bS;++bm){
          m0  = bm * b;
          for(i=0;i<b;++i){
            ii = i0 + i;
            for(j=0;j<b;++j){
              jj = j0 + j;
              for(z=0;z<b;++z){
                zz = k0 + z;

                sum = 0;
                for(s=0;s<b;++s){
                  sum += M(ii, s+m0, zz)*m2(jj, s+m0);
                }
                A1(ii,jj, zz) += sum;
              }
            }
          }
        }
      }
    }
  }

  M.resize(m1N1, m2N1, m3N1);
  M.makeZeroMatrix();
  bS1 = m1N1/b;
  bS2 = m2N1/b;
  bS3 = m3N1/b;
  bS  = AN3/b;
  for(bi=0;bi<bS1;++bi){
    i0 = bi * b;
    for(bj=0;bj<bS2;++bj){
      j0 = bj * b;
      for(bk=0;bk<bS3;++bk){
        k0 = bk * b;
        for(bm=0;bm<bS;++bm){
          m0  = bm * b;
          for(i=0;i<b;++i){
            ii = i0 + i;
            for(j=0;j<b;++j){
              jj = j0 + j;
              for(k=0;k<b;++k){
                kk = k0 + k;

                sum = 0;
                for(s=0;s<b;++s){
                  sum += A1(ii, jj, s+m0)*m3(kk, s+m0);
                }
                M(ii,jj, kk) += sum;
              }
            }
          }
        }
      }
    }
  }

  return M;
}



int main(int argc, char** argv)
{
  if (argc<3){
		cerr << "You have to pass a block size as parameter and block method number (0=naiv, 1=advanced)"<<endl;
		return -1;
	}

  matrix<double, 3> A  = readMatrix<double, 3>("mat3d");
  matrix<double, 2> m1 = readMatrix<double, 2>("mat2d_1");
  matrix<double, 2> m2 = readMatrix<double, 2>("mat2d_2");
  matrix<double, 2> m3 = readMatrix<double, 2>("mat2d_3");



  auto t0 = chrono::high_resolution_clock::now();

  matrix<double, 3> res2 = atoi(argv[2])==0 ? multiplyBlock0<double>(A, m1,m2, m3, atoi(argv[1]))
                                            : multiplyBlock<double>(A, m1, m2, m3, atoi(argv[1]));

  auto t1 = chrono::high_resolution_clock::now();

  if (argc==4 && string(argv[3])=="test"){
    matrix<double, 3> res1 = multiply<double>(A, m1, m2, m3);

    auto t2 = chrono::high_resolution_clock::now();

    matrix<double, 3> diff(res1.getDimension(0), res1.getDimension(1), res1.getDimension(2));

    auto N = res1.getDimension(0),
         M = res1.getDimension(1),
         O = res1.getDimension(2);
    double sum = 0;
    for(auto i=0, j=0, k=0;i<N;i++){
      for(j=0; j<M;j++){
        for(k=0;k<O;k++){
          diff(i, j, k) = fabs(res1(i, j, k) - res2(i, j, k));
          sum += diff(i, j, k);
        }
      }
    }
    ofstream out2("difference", ios::out);
    out2 << diff;
    out2.close();

    cout <<  "Sum of abs of difference: " << sum << endl;
    cout << "Non block version time: " << chrono::duration_cast<chrono::milliseconds>(t2-t1).count() << " ms" << endl;
    cout << "Block version time (b="<< argv[1] <<"): " << chrono::duration_cast<chrono::milliseconds>(t1-t0).count() << " ms" << endl;
  } else {
    cout << "Block size = "<< argv[1] <<"; time = " << chrono::duration_cast<chrono::milliseconds>(t1-t0).count() << " ms" << endl;
  }

  ofstream out1("result", ios::out);

  out1 << res2;

  out1.close();

  return 0;
}
