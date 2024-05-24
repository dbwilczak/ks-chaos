#ifndef DW__LDKS_POINCARE_MAP__
#define DW__LDKS_POINCARE_MAP__

// ###############################################
/**
 * This class is a wrapper for the (nonrigorous) Poincare map from the CAPD library.
 * It implements D-dimensional Galerkin projection of the KS PDE.
 * It is used for simulation and for finding approximate objects used rigorous computation.
 */ 

#include "capd/capdlib.h"

using namespace capd;
using autodiff::Node;

template<class V>
V getNormalToSection(int D){
  const static typename V::ScalarType c[100] = {typename V::ScalarType(1.)};
  return V(D,c);
}

// ###############################################

void ksVectorField(
    Node _time, 
    Node in[],   // array of input variables
    int D,       // dimension of domain
    Node out[],  // array of result variables
    int ,        // dimension of counterdomain - here ignored since the same as D
    Node p[],    // array of parameters
    int          // number of parameters
  )
{
  for(int k=1;k<=D;++k)
  {
    int n;
    Node temp(0.);
    int s = capd::min(D-k,(k-1)/2);

    for(n=s+1;n<=(k-1)/2;++n)
      temp -= in[n-1]*in[k-n-1];

    for(n=1;n<=s;++n)
      temp += in[n-1]*(in[n+k-1] - in[k-n-1]);

    for(n=s+1;n<=D-k;++n)
      temp += in[n-1]*in[n+k-1];


    if(k%2==0)
      temp -= Node(0.5)*sqr(in[k/2-1]);
    Node k2(k*k);
    out[k-1] = Node(2*k)*temp + k2*(Node(1)-p[0]*k2)*in[k-1];
  }
}

// ###############################################

void ksPoincareSection(
    Node _time, 
    Node in[],   // array of input variables
    int D,       // dimension of domain
    Node out[],  // array of result variables
    int ,        // dimension of counterdomain - here ignored since equals to 1
    Node p[],    // array of parameters
    int          // number of parameters
  )
{
  out[0] = p[D];
  for(int i=0;i<D;++i) out[0] += p[i]*in[i];
}


struct LDKSPoincareMap {
  /// constructor 
  /// @param [in] D - dimension of the projection
  /// @param [in] nu - viscosity
  /// @param [in] order - order of the numerical (Taylor) method
  /// The default Poincare section is set as x[0] = 0
  LDKSPoincareMap(int D, long double nu, int order=10)
    : vectorField(ksVectorField,D,D,1),
      solver(vectorField,order),
      section(LDVector(D),LDVector(D)),
      pm(solver,section/*,poincare::PlusMinus*/)
  {
    vectorField.setParameter(0,nu);
    LDVector n(D);
    n[0] = 1.;
    section.setNormalVector(n);
  }

  /// This routine sets a new affine Poincare section 
  /// which passes trough the point x
  /// and which is a hyperplane orthogonal to the vector field at x.
  /// @param [in]  x - a vector at which the section is set.
  /// @return - normal vector to the section (normalized vector field at x)
  LDVector setOrthogonalSectionAt(LDVector x){
    LDVector n = vectorField(x);
    n.normalize();
    section = LDAffineSection(x,n);
    return n;
  }

  /// This routine computes n-th iteration of the Poincare map, i.e.
  /// P^n(x) is returned
  /// @param [in] x - a point (not need to be on the section)
  /// @param [in] n - positive integer
  /// @return n-th itersection of the trajectory of x with the section
  LDVector operator()(LDVector x, int n=1){
    for(int i=0;i<n;++i)
      x = pm(x);
    return x;
  }

  /// This routine computes n-th iteration of the Poincare map
  /// and derivative of the flow at x. It returns
  /// P^n(x) and DPhi^n(x)
  /// @param [in] x - a point (not need to be on the section)
  /// @param [in] n - positive integer
  /// @param [out] D - a matrix tha will be set to DPhi^n(x)
  /// @return n-th itersection of the trajectory of x with the section
  LDVector flowDerivative(const LDVector& x, LDMatrix& D, int n=1){
    LDVector r = x;
    LDMatrix M(x.dimension(),x.dimension());
    D.setToIdentity();
    for(int i=0;i<n;++i){
      r = pm(r,M); 
      D = M*D;
    }
    return r;
  }

  /// This routine computes n-th iteration of the Poincare map
  /// and its derivative at x. It returns
  /// P^n(x) and DP^n(x)
  /// @param [in] x - a point (not need to be on the section)
  /// @param [in] n - positive integer
  /// @param [out] D - a matrix tha will be set to DP^n(x)
  /// @return n-th itersection of the trajectory of x with the section
  LDVector operator()(const LDVector& x, LDMatrix& D, int n=1){
    LDVector r = flowDerivative(x,D,n);
    D = pm.computeDP(r,D);
    return r;
  }

  LDMap vectorField;
  LDTaylor solver;
  LDAffineSection section;
  LDPoincareMap pm;
};

#endif
