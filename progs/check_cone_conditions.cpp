// STL libraries
#include <iostream>
#include <fstream>
#include <sstream>

// Algorithms specific for this program
#include "include/tictac.h"                /// for time of computation
#include "include/ComputeDerivativeTask.h" /// asynchronous computation of derivative of Poincare map
#include "include/HSet.h"                  /// data structure that represents an h-set

// Constants used in the computation.
const int g_order = 3;                                    /// Order of the PdeSolver.
const interval g_nu = interval(1212)/interval(10000);     /// parameter of the KS PDE

const double g_q = 1.5;             /// geometric decay of all tails
const int g_dissEncDimension = 5;   /// index of first variable on which dissipative enclosure acts
const int g_m = 18;                 /// dimension of HSets - see data structure HSet

struct HSetWithCones : public HSet{
  IMatrix P;
  capd::IMatrix Q;

  HSetWithCones(const char*filename, int M, interval q, const char* coordSystemFilename, double q1, double q2) 
    : HSet(filename,M,q), Q(-IMatrix::Identity(M-1))
  {
    Q(1,1) = q1;
    Q(2,2) = q2;
    LDMatrix A;
    ifstream in(coordSystemFilename);
    in >> A;
    P = resizeMatrix(vectalg::convertObject<IMatrix>(A),g_m);
  }
};

threading::ThreadPool threadPool(thread::hardware_concurrency());

// ###################################################################

struct CheckConeCondition {

  CheckConeCondition(const HSetWithCones& s1, const HSetWithCones& s2, double tol)
    : s1(s1), s2(s2), 
      pmFactory(g_dissEncDimension,g_m,g_nu,g_order,tol,g_q)
  {
    pmFactory.section = s2.section;
    out.precision(17);
  }
  
  const HSetWithCones s1, s2;

  IMatrix Jxx;    
  IVector Jxy;
  interval Jyy, Jyx;
  
  IKSPoincareMapFactory pmFactory;
  std::vector<ComputeDerivativeTask*> tasks;
  std::ostringstream out; // for logging the output

  /// This function creates asynchronous tasks that compute derivative of Poincare map
  /// @param [in] g1, g2 - split initial set s1 into g1*g2 pieces on main coordinares
  void initTasks(int g1, int g2){
    tasks.clear();
    GeometricBound x = GeometricBound(s1.S,g_q,s1.c);
    
    IMatrix A = s1.E*s1.B;
    A.column(0).clear(); // this is the normal direction to section
    IMatrix Q = s1.E*s1.P;
    IMatrix invQ = capd::matrixAlgorithms::krawczykInverse(s2.E*s2.P); 

    for(int i=0;i<g1;++i){
      IVector r = s1.r;
      r[1] = s1.r[1].left() + interval(i,i+1)*(s1.r[1].right()-s1.r[1].left())/g1;
      for(int j=0;j<g2;++j){
        r[2] = s1.r[2].left() + interval(j,j+1)*(s1.r[2].right()-s1.r[2].left())/g2;
        tasks.push_back(new ComputeDerivativeTask(pmFactory,x,A,r,Q,invQ,s1.iterations));
      }
    }
  }

  /// After all tasks complete computation, we merge results from these computation.
  /// Final result will be stored in members
  void mergeData()
  {  
    this->Jxx = tasks[0]->Jxx;
    this->Jxy = tasks[0]->Jxy;
    this->Jyx = tasks[0]->Jyx;
    this->Jyy = tasks[0]->Jyy;
    for(unsigned i=1;i<tasks.size();++i){
      intervalHull(tasks[i]->Jxx,this->Jxx,this->Jxx);
      intervalHull(tasks[i]->Jxy,this->Jxy,this->Jxy);
      this->Jyx = intervalHull(tasks[i]->Jyx,this->Jyx);
      this->Jyy = intervalHull(tasks[i]->Jyy,this->Jyy);
    }
  }
  
  /// Start asynchronous computation of derivative, wait for tasks and merge data.
  /// Bounds on derivative are stored in the file
  void computeDerivative(int g1, int g2, const char* filename){
    initTasks(g1,g2);
    for(auto* t : tasks) threadPool.process(t);
    for(auto* t : tasks) t->join();
    mergeData();
    for(auto* t : tasks) delete t;
    if(filename){
      ofstream out(filename);
      out.precision(17);
      out << Jxx << endl;
      out << Jxy << endl;
      out << Jyx << endl;
      out << Jyy << endl;
      out.close();
    }
  }
  
  interval minEigenvalueByGershgorin(IMatrix K){
    interval minEigenvalue;
    for(int i=1;i<=K.numberOfRows();++i)
    {
      interval lambda = 0.;
      for(int j=1;j<=K.numberOfColumns();++j)
        if(i!=j)
          lambda += capd::abs(K(i,j));
      lambda = K(i,i) + interval(-1,1)*lambda;
      out << "lambda["<<i<<"]=" << lambda << endl;
      minEigenvalue = (i==1) ? lambda : capd::min(lambda,minEigenvalue);
    }
    return minEigenvalue;
  }

  /// Assume Jxx, Jxy, Jyx, Jyy are computed.
  /// This funciton checks if the cone condition is then satisfied with matrixes s1.Q and s2.Q.
  bool checkConeCondition() {
    IMatrix J = projectMatrix(this->Jxx);
    IMatrix absQ1 = absMatrix(s1.Q);

    interval D = 1. - (absQ1*Jxy)*Jxy - sqr(Jyy);
    interval C = Jyx*Jyy;
    for(int k=1;k<=Jxy.dimension();++k)
      C += absQ1(k,k)*IEuclNorm()(IVector(J.row(k-1)))*Jxy[k-1];

    IMatrix QF =  Transpose(J)*s1.Q*J-s2.Q;
    symmetrize(QF);
    matrixAlgorithms::symMatrixDiagonalize(QF,J);
    interval A = minEigenvalueByGershgorin(J) - sqr(Jyx);
    bool result = A>0 and D>0 and A*D-sqr(C)>0;

    LOGGER(out,sqr(Jyx));
    LOGGER(out,A);
    LOGGER(out,D);
    LOGGER(out,C);
    LOGGER(out,A*D-sqr(C));
    if(!isSingular(C)) LOGGER(out,A*D/sqr(C));
    LOGGER(out,result);
    return result;
  }
};

/// run validation for parameter value 0.1212 and numerically observed hyperbolic period orbit
void validateConeCondition(HSetWithCones s1, HSetWithCones s2, int g1, int g2, double tol, const char* filename, const char* description){
  CheckConeCondition ccc(s1,s2,tol);
  ccc.computeDerivative(g1,g2,filename);     
  ccc.checkConeCondition();
  cout << description << endl << ccc.out.str() << endl;
}

// ###################################################################

int main(int argc, char** argv) {
  try {
    cout.precision(17);

    HSetWithCones N1("hsets/N1.txt",g_m,g_q,"data/coordinateSystemN1.dat",1.0,-0.5);
    HSetWithCones M("hsets/M.txt",g_m,g_q,"data/coordinateSystemM.dat",2.0,-1.0);

    HSetWithCones N2("hsets/N2.txt",g_m,g_q,"data/coordinateSystemN2.dat",1.0,-1.0);
    HSetWithCones K1("hsets/K1.txt",g_m,g_q,"data/coordinateSystemK1.dat",1.0,-1.0);
    HSetWithCones K2("hsets/K2.txt",g_m,g_q,"data/coordinateSystemK2.dat",1.0,-1.0);
    HSetWithCones K3("hsets/K3.txt",g_m,g_q,"data/coordinateSystemK3.dat",4.0,-2.0);

    /// Start clock
    auto tictac = tic(std::cout,"Check cone condition\n");

    // create an automatically run tasks
      std::thread ccTasks[] ={
        std::thread(validateConeCondition,N1,M ,2,1,1e-9,"output/N1.dat","N1=>M"),
        std::thread(validateConeCondition,M ,N1,1,1,1e-9,"output/M.dat","M=>N1"),
        std::thread(validateConeCondition,N2,K1,1,1,1e-9,"output/N2.dat","N2=>K1"),
        std::thread(validateConeCondition,K1,K2,1,1,1e-9,"output/K1.dat","K1=>K2"),
        std::thread(validateConeCondition,K2,K3,6,1,1e-9,"output/K2.dat","K2=>K3"),    
        std::thread(validateConeCondition,K3,N2,1,1,1e-9,"output/K3.dat","K3=>N2")
        };
      // wait for tasks  
      for(auto& t : ccTasks) t.join();
    
    /// Stop clock and print time of computation
    tac(std::cout,tictac);

    threadPool.interrupt();
  } catch (exception& e) {
    cout << e.what();
  }
}
