#include <iostream>
#include <fstream>
#include <sstream>

// Algorithms specific for this program
#include "include/tictac.h"            /// measures time of computation
#include "include/IKSPoincareMap.h"    /// rigorous computation of Poincare map for the KS equation
#include "include/LDKSPoincareMap.h"   /// nonrigorous computation of Poincare map of Galerkin projection
const double q = 1.5;

// ###################################################################
/// nonrigorous routine. Finds an approximate symmetric stable PO by simple iteration
LDVector findStableSymmetricPeriodicOrbit(double nu, LDVector u){
  LDKSPoincareMap pm(u.dimension(),nu,8);
  double r = 0;
  do{
    LDVector v = pm(u);
    for(int i=0;i<u.dimension();i+=2)
      v[i] = -v[i];
    r = (u-v).euclNorm();
    u = v;
  } while(r>1e-12); 
  u[0] = 0.;//  u[2] = -u[0];
  return u;
}

// ###################################################################
/// nonrigorous routine. Find approximate eigenvectors of S\circ P
/// where S is the symmetry at approximate fixed point.
LDMatrix findDiagonalCoordsAtSymmetricPO(double nu, LDVector u){
  int D = u.dimension();
  LDKSPoincareMap pm(D,nu,6);
  LDMatrix A(D,D);
  LDVector v = pm(u,A);
  /// apply the reflectional symmetry S
  A = symMatrix(A);
  LDMatrix B = projectMatrix(A);
  
  int d = B.numberOfRows();
  LDVector rV(d), iV(d);         // real and imaginary parts of eigenvalues
  LDMatrix rVec(d,d), iVec(d,d); // real and imaginary parts of eigenvectors
  capd::alglib::computeEigenvaluesAndEigenvectors(B,rV,iV,rVec,iVec);
  for(int i=0;i<d;++i) {
    rVec.column(i).normalize();
    if(iV[i]!=0.){
      rVec.column(i) = iVec.column(i+1);
      ++i;
    }
  }
  /// embed again in the full dimension
  /// adding normalized vector field at u as the first column  
  A = embedMatrix(rVec);
  A.column(0) = pm.vectorField(u);
  return A;
}

/// Auxililary rotuine.
/// It finds a stable periodic point by simple iteration.
/// Then it computes an approximate Jordan basis at this point. 
/// This will be input for rigorous computaiton.
std::pair<IVector,IMatrix> initPerPoint(interval nu, int m){
  // find approximate fixed point, here is a candidate
  long double fp_0127[100]  = {0.3880963864168981863369,1.320636180607189233806,-0.38809638641689818612,-0.3441693503991267324797,0.1064018499502311431793,0.03224480081440828437787,-0.01530749744245004089276,-0.001967432611581025603906,0.001665882888158957331689,2.792702856471190139358e-05,-0.0001474119864200247229072,1.04190055590547967658e-05,1.110652829460937019944e-05,-1.765010587333960545915e-06,-7.135295085766261252834e-07,1.93381798328456584331e-07,3.746523970866851961057e-08,-1.731069590744421978845e-08,-1.319849112872784146994e-09,1.351029049304772606711e-09};
  LDVector u(m,fp_0127);
  // refine it by a simple iteration
  std::cout << "Finding candidate for a stable periodic orbit by simple iteration ...\n" << std::flush;
  u = findStableSymmetricPeriodicOrbit(nu.rightBound(),u);

  // find approximate derivative of S\circ P, where S is the reflectional symmetry    
  LDMatrix A = findDiagonalCoordsAtSymmetricPO(nu.rightBound(),u);
  return {vectalg::convertObject<IVector>(u),vectalg::convertObject<IMatrix>(A)};
}


/// auxiliary routine that prints double in a Latex friendly format
std::string writeDouble(double c){
  if(c!=0){
    int n = floorl(log10(fabs(c)));
    std::ostringstream out;
    out.precision(17);
    out << c*pow(10,-n);
    if(n!=0) out << "\\cdot 10^{" << n << "}";
    return out.str();
  } else{
    return "0"; 
  }
}

/// auxiliary routine. Writes an approximate periodic point as a Latex array - to be included in the article
void writeVector(std::ostream& out, LDVector u){
  out.precision(17);
  out << "\\begin{array}{l|l}\n";
  for(int i=0;i<u.dimension();++i){
    out << "u^0_{" <<(i+1) << "}=" << writeDouble(u[i]);
    out << (i%2 ? "\\\\\n" : " & ");
  }
  out << "\n\\end{array}\n";
}

/// auxiliary routine. Writes computed bounds in a Latex format to be included in the article
void writeData(std::ostream& out, IVector u, IVector r, double Cu, double Cr){
  out.precision(17);
  out << "\\begin{array}{|c|c|c|}\n";
  out << "i & u  & r \\\\\n";
  for(int i=0;i<u.dimension();++i){
    // i+2 because indexing in the article starts from 1 and we skip the section one coordinate a_1=0 
    out << (i+2) << " & " << writeDouble(abs(u[i]).rightBound()) << " & " << writeDouble(abs(r[i]).rightBound()) << "\\\\\n";
  }
  out << "C & " << writeDouble(Cu) << " & " << writeDouble(Cr) << "\\\\\n";
  out << "\\end{array}\n"; 
}

// ###################################################################
/// auxiliarty routine. Writes computed C^0 bounds to a given stream and in Latex format
void logStablePOC0Data(std::ostream& out, const IVector& pi_r, const IVector& pi_u, const GeometricBound& r, const GeometricBound& x, bool validated){
  writeData(out,pi_u,pi_r,x.getConstant().rightBound(),r.getConstant().rightBound());  
  LOGGER(out,validated);
}

// ###################################################################
/// auxiliarty routine. Writes computed C^1 bounds to a given stream
void logStablePOC1Data(std::ostream& out, interval PxxNorm, interval PxyNorm, interval PyxNorm, interval PyyNorm, interval norm, bool validated){
  LOGGER(out,PxxNorm);
  LOGGER(out,PxyNorm);
  LOGGER(out,PyxNorm);
  LOGGER(out,PyyNorm);
  LOGGER(out,norm);
  LOGGER(out,validated);
}

// ###################################################################
/// @return invA*(P(x0+A*r)-y0)
GeometricBound computeImage(IKSPoincareMap& pm, GeometricBound x0, IMatrix A, IVector r, GeometricBound y0, IMatrix invA){
  interval returnTime = 0.;
  auto t = tic(cout,"Computing image");
  GeometricBound result = pm.computeImage(x0,A,r,resizeMatrix(invA,x0.dimension()),y0,1/*period*/,returnTime);  
  tac(cout,t);
  result[0] = 0;    
  return result;
}
  
// ###################################################################
/// @param [in] A is an approximate matrix that diagonalizes S*DP(x0) (on the section)
///             and the first column is the flowdir direction
/// @param [in] x0 =(x1,...,xM,0,...) is an approximate fixed point for Galerkin projection
bool validateStablePO(IKSPoincareMap& pm, GeometricBound x0, IMatrix A, double inflateFactor){
  writeVector(cout,vectalg::convertObject<LDVector>(x0.getExplicitCoefficients()));

  bool validated;
  GeometricBound y0 = symVector(x0); // y0 = S(x0)
  IMatrix invSA = matrixAlgorithms::krawczykInverse(symMatrix(A));  // invSA = (S*A)^{-1}  
  
  // preliminary computation pm(x0), here r = 1e-9. We want to find a candidate set K for validation.
  // See article for details.
  std::cout << "\nPreliminary computation of a candidate set K for further validation ...\n" << std::flush;
  GeometricBound r0 = computeImage(pm,x0,A,IVector(A.numberOfRows()),y0,invSA);

  // now we set candidate for validation of a fixed point. Blow up r0.
  GeometricBound x = x0;
  x.setConstant(r0.getConstant()*1.25);
  IVector u = r0.getExplicitCoefficients()*(inflateFactor*interval(-1,1));

  // here we compute r = invSA*(P(x+A*u)-y0), where y0 = S(x0) and A diagonalizes S*DP(x0)
  // that is the result is represented as y0 + (SA)*r
  std::cout << "\nCompute PS(K) and check the inclusion PS(K)\\subset K ...\n" << std::flush;
  GeometricBound r = computeImage(pm,x,A,IVector(A.numberOfRows(),u.begin()),y0,invSA);
  
  // We have to check that SP(x0 + A*u)\subset x0 + A*u
  // From the above we have
  // SP(x0+A*u) \subset S(y0 + (SA)*r) = x0 + A*r
  // It suficces to check the inclusion r\subset u 
  // First we take the projection to the section (remove the 0-th coordinate)
  IVector pi_r = IVector(r.getExplicitCoefficients().dimension()-1,r.getExplicitCoefficients().begin()+1);
  IVector pi_u = IVector(u.dimension()-1,u.begin()+1);
  validated = subsetInterior(pi_r,pi_u) and r.getConstant() < x.getConstant();

  logStablePOC0Data(std::cout,pi_r,pi_u,r,x,validated);

  if(!validated) return false;
  
  // If validated, start C^1 computation.
  // We have to check that the norm of SDP(x) is less than 1.
  // Since S is an involution, we have
  // A^-1 SDP(x) A = invSA*DP(x)*A
  // Here we compute invSA*DP(x)*A for x\in X:=x0 + A*r.
  // From previous computation we know that the fixed point belongs to X
  x = x0;
  x.setConstant(r.getConstant());
  
  interval returnTime = 0.;
  IVector Jxy;
  interval Jyx, Jyy;
  
  auto t = tic(cout,"Computing derivative");
  IMatrix Jxx = pm.computeDerivative(x,A,r.getExplicitCoefficients(),1,A,invSA,returnTime,Jxy,Jyx,Jyy); 
  tac(cout,t);
  
  interval PxxNorm = IEuclNorm()(projectMatrix(Jxx));
  interval PxyNorm = IEuclNorm()(Jxy).rightBound();
  interval PyxNorm = Jyx;
  interval PyyNorm = Jyy.rightBound();
  interval norm =  capd::max(PxxNorm + PxyNorm,PyxNorm + PyyNorm).rightBound(); 
  validated = (norm<1.);

  logStablePOC1Data(std::cout,PxxNorm,PxyNorm,PyxNorm,PyyNorm,norm,validated);
  
  return validated; 
}

// ###################################################################

int main(int argc, char** argv) {
  try {
    cout.precision(17);

    /// parameter of the system
    const interval nu = interval(127)/1000.;

    /// Control data
    const int m = 15;                // number of explicit coefficients
    const int order = 2;             // order of Taylor method
    const double tol = 1e-7;         // tolerance per time step of the PdeSolver
    const int dissipativeEncDim = 9; // index of first coordinate on which dissipative enclosure is used
    double inflateFactor = 3;
    // other parameters:  m=18, order=3, tol=1e-11, inflateFactor = 10
    
    // find an approximate periodic point by simple iteration
    // A is an approximate Jordan basis of a m-dimensional block of SDP(x0)
    auto [_x0,A] = initPerPoint(nu,m);
    
    // create an instance of Poincare map
    GeometricBound x0(m,0.,q,_x0.begin());
    IKSPoincareMap pm(dissipativeEncDim,m,nu,order,tol,q);
    
    // call validation routine  
    validateStablePO(pm,x0,A,inflateFactor);

  } catch (exception& e) {
    cout << e.what() << endl;
  }
}
