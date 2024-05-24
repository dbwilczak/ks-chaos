#ifndef DW___HSET_H__
#define DW___HSET_H__

#include <fstream>
#include <vector>
#include "capd/capdlib.h"
#include "include/projections.h"

/**
 * HSet is a data structure which represents a set on the section as
 *   X = c + B*r
 * where
 * c - is the origin of the coordinate system
 * B - is a linear coordinate system
 * r - is a box centered at the origin
 * We assume that r[0] = 0 and this is the 'normal variable'.
 * The matrix inv_B is a rigorous inverse of B.
 * 
 * The set is embeded into Poincare section by an orthogonal mapping E, i.e.
 *   Y = c + E*B*r
 * The matrix E maps first coordinate of B*r into ortogonal direction to the section.
 * The remaining columns of E form an orthogonal basis on the section.
 * 
 * The matrix inv_E is a rigorous inverse of E.
 *  
 */ 
struct HSet{
  capd::IVector c;
  capd::IMatrix B, inv_B;
  capd::IVector r;
  capd::IMatrix E, inv_E;
  double S;
  
  capd::IAffineSection section;
  int iterations;

  /// This constructor initializes an h-set of dimension m.
  /// All coordinates should be initialized later.
  HSet(int m=1) 
    : c(m), B(m,m), inv_B(m,m), r(m), 
      E(capd::IMatrix::Identity(m)), inv_E(capd::IMatrix::Identity(m)),
      section(capd::IVector(m),capd::IVector(m))
  {}

  /// This method sets a Poincare section in which the h-set is embedded. 
  /// The section is computed as an affine plane attached art the point _c,
  /// which is orthogonal (wrt finite dimensional l_2 inner product) 
  /// to the vector _n.
  void setCenterAndSection(capd::LDVector _c, capd::LDVector _n){
    this->c = capd::vectalg::convertObject<capd::IVector>(_c);
    capd::IVector n = capd::vectalg::convertObject<capd::IVector>(_n);
    this->E.setToIdentity();
    this->E.column(0) = n;
    orthonormalizeColumns(this->E);
    this->inv_E = Transpose(this->E);
    this->section = capd::IAffineSection(this->c,n);
  }

  /// This constructor reads coordinates of an h-set from a file.
  HSet(const char* filename) : HSet(1){
    capd::DMatrix _E,_B;
    capd::DVector _c;
    std::ifstream in(filename);
    auto reader = [&in](auto& t) { while(in.get()!='='){} in >> t;};
    reader(S); reader(_c); reader(_E); reader(_B); reader(r); reader(iterations);
    in.close();  

    this->c = capd::vectalg::convertObject<capd::IVector>(_c);
    this->E = capd::vectalg::convertObject<capd::IMatrix>(_E);
    this->B = capd::vectalg::convertObject<capd::IMatrix>(_B);

    this->inv_E = capd::matrixAlgorithms::krawczykInverse(this->E);
    this->inv_B = capd::matrixAlgorithms::krawczykInverse(this->B);

    this->section = capd::IAffineSection(this->c,this->E.column(0));  
  } 

  /// This constructor reads coordinates of an h-set from a file.
  /// Then it embeds it into (higher) dimension M.
  HSet(const char*filename, int M, interval q) : HSet(filename) {
    int m = c.dimension();
    c = resizeVector(c,M);
    r = resizeVector(r,M);
    B = resizeMatrix(B,M);
    inv_B = resizeMatrix(inv_B,M);
    E = resizeMatrix(E,M);
    inv_E = resizeMatrix(inv_E,M);
    
    /// now we have to updated tail in r and the constant S
    for(int i=m+1;i<=M;++i)
      r[i-1] = S*interval(-1,1)/power(q,i);
    S = (S/power(q,M-m)).rightBound();
  } 

  /// shifts the h-set and the section to a new 'origin' _c
  void setCenter(capd::LDVector _c){
    this->c = capd::vectalg::convertObject<capd::IVector>(_c);
    this->section.setOrigin(this->c);
  }

  /// sets new shape matrix of the h-set
  void setCoordinateSystem(capd::IMatrix _B){
    this->B = _B;
    this->inv_B = capd::matrixAlgorithms::krawczykInverse(B);
  }

  /// sets new shape matrix of the h-set
  void setCoordinateSystem(capd::LDMatrix _B){
    setCoordinateSystem(capd::vectalg::convertObject<capd::IMatrix>(_B));
  }
};

#endif

