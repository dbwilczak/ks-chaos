#ifndef DW_PROJECTIONS_H__
#define DW_PROJECTIONS_H__
#include "capd/capdlib.h"

/**
 * This file contains auxiliary routines that are responsible for embeding and projecting vectors and matrices from/to Poincare sections.
 * Note, the integration od the infinite-dimensional ODE is performed in the full space, while the sections are affine subspaces.
 */ 
 
/// Applies symmetry S(a1,a2,a3,...) = (-a1,a2,-a3,...) to a vector
template <class V>
V symVector(V v){
  for(int i=0;i<v.dimension();i+=2)
    v[i] = -v[i];
  return v;
}

/// Applies the same symmetry to a matrix. That is it computes S*m
template<class M>
M symMatrix(M m){
  for(int i=0;i<m.numberOfRows();i+=2)
    m.row(i) *= -1.;
  return m;
}

/// Removes first coordinate from a finite-dimensional vector 
template<class V>
V projectVector(const V& v){
  return V(v.dimension()-1,v.begin()+1);
}

/// Extends a finite-dimensional vector by adding zero as the first coordinate
template<class V>
V embedVector(V u){
  V result (u.dimension()+1);
  for(int i=0;i<u.dimension();++i)
    result[i+1] = u[i];
  return result;
}

/// Project finite-dimensional matrix by removing first row and first column
template<class M>
M projectMatrix(const M& m){
  int D = m.numberOfRows();
  M r(D-1,D-1);
  for(int i=1;i<D;++i)
    for(int j=1;j<D;++j)
      r(i,j) = m(i+1,j+1);
  return r;
}

/// Enlarge matrix by addind leading row and column that contain zeroes, except r(1,1)=1
template<class M>
M embedMatrix(const M& m){
  int D = m.numberOfRows();
  M r(D+1,D+1);
  r(1,1) = 1.;
  for(int i=1;i<=D;++i)
    for(int j=1;j<=D;++j)
      r(i+1,j+1) = m(i,j);
  return r;
}

/// Intresects each pair of coefficients A_{ij} and A_{ji} of an interval matrix 
/// and sets the results back to A_{ij} and A_{ji}
template<class M>
void symmetrize(M& A){
 for(int i=1;i<=A.numberOfRows();++i){
    for(int j=i+1;j<=A.numberOfRows();++j){
      intersection(A(i,j),A(j,i),A(i,j));
      A(j,i) = A(i,j);
    }
  }
}

/// Assuming M is bigger than the dimensions of m1
/// the routine sets martix m2 to a square matrix with m2_{ij} = m1_{ij}
/// and adding 1 on the diagonal
template<class M1, class M2>
void expandMatrix(const M1& m1, M2& m2, int M){
  m2 = M2::Identity(M);
  for(int i=1;i<=m1.numberOfRows();++i)
    for(int j=1;j<=m1.numberOfColumns();++j)
      m2(i,j) = m1(i,j);
}

/// Computes martix |A_{ij}|
template<class M>
M absMatrix(M A){
  for(int i=1;i<=A.numberOfRows();++i)
    for(int j=1;j<=A.numberOfColumns();++j)
      A(i,j) = capd::abs(A(i,j));
  return A;
}

/// Given a square matrix 'm', resize it to dimension 'dim' by adding 1 on the diagonal.
template<class M>
M resizeMatrix(const M& m, int dim){
  int i,j;
  int D = m.numberOfRows();
  M r(dim,dim);
  int d = capd::min(dim,D);
  for(i=1;i<=d;++i)
    for(j=1;j<=d;++j)
      r(i,j) = m(i,j);
  for(i=d+1;i<=dim;++i)
    r(i,i) = 1.;
  return r;
}

/// Given a vector u, resize it to dimension dim by adding zeroes to the tail
template<class V>
V resizeVector(const V& u, int dim){
  int d = capd::min(dim,(int)u.dimension());
  V r(dim);
  for(int i=0;i<d;++i)
    r[i] = u[i];
  return r;
}

/// Given a matrix Q, normalize its columns.
template<class M>
void normalizeColumns(M& Q){
  for(int i=0;i<Q.numberOfColumns();++i)
    Q.column(i).normalize();
}

/// Perform QR-decomopistion of Q and store the result again in Q.
template<class M>
void orthonormalizeColumns(M& Q){
  M A = Q;
  M R = Q;
  capd::matrixAlgorithms::QR_decompose(A,Q,R);
}

/// Perform QR-decomopistion of Q and store the result again in Q.
template<class M>
void orthonormalizeColumns(M& Q, int d){
  M A = Q;
  orthonormalizeColumns(Q);
  for(int i=0;i<d;++i){
    Q.column(i) = A.column(i);
    Q.column(i).normalize();
  }
}

#endif
