/**
 * This class defines data types for rigorous integrator odf the KS equation.
 * Note, that the CAPD library provides a general class that computes Poincare map of an abstract system.
 * Parameters of this class are
 * - Solver - that integrates the system
 * - Section - there are three types, but we use AffineSection
 *
 * In capdDynSys/pdes directory of the CAPD library we have implemented classes
 * - PdeSolver - it integrates an abstract Pde provided the vector field is defined
 * - PdeAffineSection - class that represents affine Poincare section
 *
 * This directory also contains code specific for the KS equation
 * - OneDimKSSineVectorField - defines vector field for the KS equation (parameter to PdeSolver) and implements all the bound presented in the article.
 */ 
#ifndef _IKS_POINCARE_MAP_
#define _IKS_POINCARE_MAP_

#include "include/typedefs.h"
#include "projections.h"

// ###############################################
/**
 * This class is a wrapper for the (rigorous) Poincare map from the CAPD library.
 * It provides algorithms for computation of Poincare map of KS PDE in given coordinates.
 */ 
struct IKSPoincareMap {
  /// constructor 
  /// @param [in] _m - dimension of the projection
  /// @param [in] nu - viscosity
  /// @param [in] order - order of the numerical method
  /// @param [in] dissEncDim - index of first dissipative coordinate. Used in the algorithm that finds an enclosure.
  /// @param [in] _q - geometric decay of the tail
  /// The default Poincare section is set as x[0] = 0
  IKSPoincareMap(int dissEncDim, int _m, interval nu, int order, double tol, interval _q)
    : m(_m),                                               
      vectorField(nu,_m,dissEncDim),
      solver(vectorField,order),
      section(GeometricBound(_m,0.,_q),GeometricBound(_m,0.,_q)),
      pm(solver,section),
      q(_q)
  {
    solver.setAbsoluteTolerance(tol);
    solver.setRelativeTolerance(tol);
    GeometricBound n(m,0.,q);
    n[0] = 1.;
    section.setNormalVector(n);
  }
  
  /// This routine computes n-th iteration of the Poincare map.
  /// Result is represented in the coordinate system Q centered at u, i.e.
  /// Q*(P^n(x+C*r0)-u) is returned
  /// @param [in] x - center of the set
  /// @param [in] C - shape matrix of the set
  /// @param [in] r - the set in coordinates C
  /// @param [in] Q - linear coordinate system in target section
  /// @param [in] u - origin of the coordinate system in target section
  /// @param [in] n - positive integer
  /// @param [out] returnTime - bound on return time to the section
  /// @return bound for Q*(P^n(x+C*r+tail)-u)
  GeometricBound computeImage(GeometricBound x, IMatrix C, IVector r, IMatrix Q, IVector u, int n, interval& returnTime ){
    C0DoubletonSetGeometricTail initialCondition(x,C,r);

    returnTime=0.;
    /// Note: the number of explicit coordinates in h-set does not need to be the same as
    ///       M = the number of explicit coordinates used by rigorous integrator PdeSolver
    Q = resizeMatrix(Q,m);
    u = resizeVector(u,m);
    /// call algorithm from CAPD
    return pm(initialCondition,toSeries(u,0.),Q,returnTime,n);
  }

  /// This routine computes n-th iteration of the Poincare map.
  /// Result is represented in the coordinate system Q centered at u, i.e.
  /// Q*(P^n(x+C*r0)-u) is returned
  /// @param [in] x - center of the set
  /// @param [in] C - shape matrix of the set
  /// @param [in] r - the set in coordinates C 
  /// @param [in] S - constant used to define tail
  /// @param [in] Q - linear coordinate system in target section
  /// @param [in] u - origin of the coordinate system in target section
  /// @param [in] n - positive integer
  /// @param [out] returnTime - bound on return time to the section
  /// @return bound for Q*(P^n(x+C*r+tail)-u)
  GeometricBound computeImage(IVector x, IMatrix C, IVector r, double S, IMatrix Q, IVector u, int n, interval& returnTime ){
    return computeImage(toSeries(x,S),C,r,Q,u,n,returnTime);
  }

  /// This routine computes derivative of n-th iteration of Poincare map.
  /// D P^n(x+C*r+tail)
  /// @param [in] x - center of the set
  /// @param [in] C - shape matrix of the set
  /// @param [in] r - the set in coordinates C
  /// @param [in] n - positive integer
  /// @param [out] returnTime - bound on return time to the section
  /// @param [out] Py - bound on operator norms of rows in the block D_{xy}P^n 
  /// @param [out] Pyx - bound on operator norm of the block D_{yx}P^n
  /// @return bound for D_{xx}P^n
  IMatrix computeDerivative(GeometricBound x, IMatrix C, IVector r, int n, IMatrix Q, interval& returnTime, IVector& Py, interval& Pyx){
    IMatrix DP(r.dimension(),r.dimension());
	  
    GeometricBound X = x;
    X.getExplicitCoefficients() += C*r;

    IVector tmp = r;
    split(tmp,r);
    x.getExplicitCoefficients() += C*tmp;
    typedef C1DoubletonSetGeometricTail::FiniteDimensionalBaseSet C1Set;
    C1Set set = C1Set(C1Set::C0BaseSet(x.getExplicitCoefficients(),C,r),C1Set::C1BaseSet(Q));
    C1DoubletonSetGeometricTail initialCondition(X,set);
    initialCondition.printLog = false;
    auto u = pm(initialCondition,DP,returnTime,n);
    auto P = pm.getVectorField()(0.,u);
    auto& secDerEnc = pm.getSectionDerivativesEnclosure();
    DP = section.computeDP(u,P,DP,secDerEnc.encDy,secDerEnc.encDyx,Py,Pyx);
    return DP;
  } 

  /// This routine computes derivative of n-th iteration of Poincare map.
  /// D P^n(x+C*r+tail)
  /// @param [in] x - center of the set
  /// @param [in] C - shape matrix of the set
  /// @param [in] r - the set in coordinates C
  /// @param [in] n - positive integer
  /// @param [in] Q - new coordinate system
  /// @param [in] invQ - some matrix, prossibly invQ = Q^{-1} but not required
  /// @param [out] returnTime - bound on return time to the section
  /// @param [out] Py - bound on operator norms of rows in the block D_{xy}(invQ*P^n*Q)
  /// @param [out] Pyx - bound on operator norm of the block D_{yx}(invQ*P^n*Q)
  /// @return bound for invQ*D_{xx}P^n*Q
  IMatrix computeDerivative(GeometricBound x, IMatrix C, IVector r, int n, IMatrix Q, IMatrix invQ, interval& returnTime, IVector& Jxy, interval& Jyx, interval& Jyy){
    IVector Py = IVector(m+1);
    interval Pyx;
    IMatrix DP = this->computeDerivative(x,C,r,n,Q,returnTime,Py,Pyx);
    IVector tmp = absMatrix(invQ)*IVector(Py.dimension()-1,Py.begin()); 
    IMatrix Jxx = invQ*DP;
    Jxy = IVector(m-1,tmp.begin()+1);
    Jyy = Py[DP.numberOfRows()];
    Jyx = Pyx;
    //Jyx = IEuclNorm()(Q)*Pyx;
    return Jxx;
  } 

  /// set new Poincare section
  void setSection(const IAffineSection& s){   
    this->section.setOrigin(toSeries(s.getOrigin(),0.));
    this->section.setNormalVector(toSeries(s.getNormalVector(),0.));
  }

private:
  /// This method creates an instance of GeometricBound from an interval vector x and constant S.
  /// The leading coeffiients of GeometricBound are taken from x.
  /// The tail is set to S*q^{-i}, for i>dimension(x).
  GeometricBound toSeries(IVector x, double S){
    GeometricBound result(m,S,q);
    unsigned i;
    for(i=0;i<x.dimension();++i) result[i] = x[i];
    for(;i<m;++i) result[i] = interval(-S,S)/power(interval(q),i);
    return result;
  }

  int m;
  OneDimKSSineVectorField vectorField;
  PdeSolver solver;
  PdeAffineSection section;
  PdePoincareMap pm;
  interval q;
};


/// Some parameters remain constant in computation. We set them in Factory class to simplify construction of PdeSolvers
struct IKSPoincareMapFactory{
  int m, dissEncDim, order;
  double tolerance;
  interval nu, q;
  IAffineSection section;
  
  IKSPoincareMapFactory(int _dissEncDim, int _m, interval _nu, int _order, double _tolerance, interval _q) 
    : m(_m ), dissEncDim(_dissEncDim), order(_order), tolerance(_tolerance), nu(_nu), q(_q),
      section(IAffineSection(IVector(m), IVector(m)))
  {}
  
  IKSPoincareMap* createPoincareMap(const IAffineSection& s) const{
    std::unique_lock<std::mutex> lock(poincareMapConstructorGuard);
    IKSPoincareMap* r = new IKSPoincareMap(dissEncDim,m,nu,order,tolerance,q);
    r->setSection(s);
    return r;
  }

  IKSPoincareMap* createPoincareMap() const{
    std::unique_lock<std::mutex> lock(poincareMapConstructorGuard);
    IKSPoincareMap* r = new IKSPoincareMap(dissEncDim,m,nu,order,tolerance,q);
    r->setSection(section);
    return r;
  }
private:
  static std::mutex poincareMapConstructorGuard;
};
std::mutex IKSPoincareMapFactory::poincareMapConstructorGuard;

#endif
