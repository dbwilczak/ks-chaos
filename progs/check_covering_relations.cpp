// Standard libraries and the CAPD library
#include <iostream>
#include <fstream>
#include <sstream>
#include "capd/auxil/ConfigFileReader.h"

// Algorithms specific for this program
#include "include/tictac.h"           /// measures time of computation
#include "include/ComputeImageTask.h" /// asynchronous computation of Poincare map
#include "include/HSet.h"             /// data structure that represents an h-set

// Synchronize writting to screen by many tasks
std::mutex outputGuard;

// Constants used in the computation.
const int g_order = 4;                                    /// Order of the PdeSolver.
const interval g_nu = interval(1212)/interval(10000);     /// parameter of the KS PDE

const double g_q = 1.5;             /// geometric decay o tails
const int g_dissEncDimension = 9;   /// index of first variable on which dissipative enclosure acts
const int g_m = 15;                 /// dimension of HSets - see structure HSet

// #####################################################################

// All h-sets that appear in the article
// Coordinates are read from files hsets/*.txt
HSet *N1, *N2, *M, *K1, *K2, *K3; /// sets centered at approximate periodic points a1, a2
std::vector<HSet*> h12;           /// heteroclinic chain a1->a2
std::vector<HSet*> h21;           /// heteroclinic chain a2->a1

/// Read coordinates of all h-sets from files (precomputed data)
void initHSets(){
  N1 = new HSet("hsets/N1.txt");
  N2 = new HSet("hsets/N2.txt");
  M  = new HSet("hsets/M.txt");
  K1 = new HSet("hsets/K1.txt");
  K2 = new HSet("hsets/K2.txt");
  K3 = new HSet("hsets/K3.txt");
  for(int i=0;i<=10;++i){
    ostringstream filename1, filename2;
    filename1 << "hsets/N" << i << "_12.txt";
    filename2 << "hsets/N" << i << "_21.txt";
    h12.push_back(new HSet(filename1.str().c_str()));
    h21.push_back(new HSet(filename2.str().c_str()));
  }
}

// #####################################################################
/**
 * Let s1 be an HSet. This class is a task that asynchronously computes images of 
 * - left face of the HSet s1
 * - right face of the HSet s1
 * - entire set s1
 * 
 * After all tasks are completed, this object checks if the inequalities required for covering relation s1=>s2 are satisfied.
 */ 
class CheckCoveringRelation{
public:
  static auxil::ConfigFileReader reader;
  
  IKSPoincareMapFactory pmFactory = IKSPoincareMapFactory(g_dissEncDimension,0,g_nu,g_order,0.0,g_q);

  CheckCoveringRelation(threading::ThreadPool& _pool, HSet* _s1, HSet* _s2, std::string _d)
    : pool(_pool), s1(*_s1), s2(*_s2), description(_d)
  {
    
    /// reader reads parameters of PdeSolver from the config file.
    /// These are: m = number of explicitly stored coorfinates in integration of PDE and tol=acceptable tolerance per one step of integration.
    /// See file proofData/config.dat
    /// Remaining parameters of the solvers are fixed for all covering relations (see global variables)
    pmFactory.m = reader.get<int>(_d+".M");
    pmFactory.tolerance = reader.get<double>(_d+".tol");

    /// C is the shape matrix of the set expressed in the coordinates on the section.
    /// and it is equal to product E*B, where E is the shape matrix of s1 and B are local orthogonal coordiantes on the section 
    IMatrix C = s1.E * s1.B;
    /// Q is the cooridinate system on the target section, in which the image wil be represented.
    /// If s1=s2 then Q=C^{-1} and we compute C^{-1}*P*C.
    /// Otherwise, we compute C2^{-1}*P*C1
    IMatrix Q = s2.inv_B * s2.inv_E;

    /// create task for computation of P(set)
    this->set = new ComputeImageTask(s1.c,C,s1.r,Q,s2.c,s1.S,s1.iterations,s2.section,pmFactory);
    
    /// create task for computation of P(rightFace(set))
    IVector r = s1.r;
    r[1] = r[1].rightBound();
    this->right = new ComputeImageTask(s1.c,C,r,Q,s2.c,s1.S,s1.iterations,s2.section,pmFactory);

    /// create task for computation of P(leftFace(set))
    r = s1.r;
    r[1] = r[1].leftBound();
    this->left = new ComputeImageTask(s1.c,C,r,Q,s2.c,s1.S,s1.iterations,s2.section,pmFactory);

    /// Start computation: execute method 'run' in a separate thread.
    this->t = new std::thread(&CheckCoveringRelation::run,std::ref(*this));
  }
  
  ~CheckCoveringRelation(){
    delete left;
    delete right;
    delete set;
  }

  void join(){
    if(t==nullptr)
      throw std::runtime_error("CheckCoveringRelation error: main thread has not been started!");
    t->join();
    delete t;
  }
  
  bool result() { return validated; }

private:
  /// This is the main function called asynchronously from the constructor
  void run(){
    // start tasks that compute Poincare map on left, right face and the set
    pool.process(left);
    pool.process(right);
    pool.process(set);
    // wait until they complete
    left->join();
    right->join();
    set->join();

    this->check();    // check inclusions for covering relation
    this->report();   // print report to stdout
  }

  /// auxiliary function, which checks if computed image 't' satisfies predicate 'condition'
  template<class Condition>
  static bool checkCondition(ComputeImageTask* t, Condition condition){
    return (t->status == ComputeImageTask::Status::Completed) and condition(t->result);
  }

  /// GeometricSeries x has M-explicit coefficients, while infinite tail is bounded by S*q^{-i}.
  /// This function refines the constant S to the new dimension N<=M by bounding explicit coefficients N+1,...,M to the form new_S*q^{-i}
  static double getConstant(const GeometricBound& x, int N){
    double S = x.getConstant().rightBound();
    for(int i=N;i<x.dimension();++i){
      interval t = abs(x[i])*power(interval(g_q),i);
      S = capd::max(S,t.rightBound());
    }
    return S;
  }

  /// This function checks inclusions for the covering relation s1=>s2
  bool check(){
    IVector r = s2.r; 
    double S = s2.S;
    /// lambda function (predicate), which checks if the projection of 'im' onto stable coordinates is inside the interior of r  
    auto accross = [r,S](const GeometricBound& im){
      IVector pi_im = IVector(r.dimension()-2,im.getExplicitCoefficients().begin()+2);
      IVector pi_r = IVector(r.dimension()-2,r.begin()+2);
      /// Check inequality on finite number of coordinates and compare constants that give size of geometric tail.
      /// The constant in the image should be smaller than that of the set to be covered.
      return subsetInterior(pi_im,pi_r) and getConstant(im,r.dimension()) < S;
    };
    
    /// lambda function (predicate), which cheks if the projection onto unique unstable coordinate is ouside of r[1].
    auto mapaway = [r](const GeometricBound& im){
      return abs(im[1])>abs(r[1]);
    };
    
    /// here we check the conditions for covering relation and return the result
    this->validatedLeft = checkCondition(left,mapaway);
    this->validatedRight = checkCondition(right,mapaway);
    this->validatedAcross = checkCondition(set,accross);
    /// Moreover two faces of s1 must be mapped into oposite sides of s2.r[1].
    /// It suffices to check that they have oposite signs.
    this->validatedToOpositeSides = left->result[1]*right->result[1] < 0.;
    this->validated = 
          this->validatedLeft and
          this->validatedRight and
          this->validatedAcross and
          this->validatedToOpositeSides;
    return this->validated;
  }
  
  /// This method prints to the standard output log from the computation and checking the inequalities
  void report(){
    std::ostringstream out;
    out.precision(17);
    out << "############ " << description << " ################" << std::endl;
    if(this->validated)
      out << description << " checked" << std::endl;
    else
      out << description << " failed" << std::endl;
    /// projections onto Poincare section, i.e. we remove first irrelevant coordinate.
    int D = s2.c.dimension()-1;
    IVector imR(D,right->result.getExplicitCoefficients().begin()+1);
    IVector imL(D,left->result.getExplicitCoefficients().begin()+1);
    IVector imS(D,set->result.getExplicitCoefficients().begin()+1);
    out << "\n1. First (unstable) coordinate of right and left images must be mapped into oposite sides of the first coordinate of the box to be covered."
        << "\n2. All but first coordinates of the 'set image' must be mapped into corresponding coordinates of the box to be covered."
        << "\n3. The constant that gives uniform geometric bound on the tail in 'set image' must be smaller than that of target HSet."
        << "\n\nbox to be covered:\n" << IVector(D,s2.r.begin()+1)
        << "\n\nrightImage:\n" << imR
        << "\n\nleftImage:\n" << imL
        << "\n\nsetImage:\n"  << imS
        << "\n\nBoundOnTail(setImage)  : "  << getConstant(set->result,s2.r.dimension()) << " * " << g_q << "^{-i}"
        << "\nBoundOnTail(TargetHSet): "  << s2.S  << " * " << g_q << "^{-i}"
        << "\n\nBound on return time: " << set->returnTime
        << std::endl;

    std::unique_lock<std::mutex> lock(outputGuard);
    std::cout << out.str() << std::flush;
  }
  
  threading::ThreadPool& pool;
  /// check covering between s1 and s2
  HSet& s1, s2;
  /// text description of covering relation
  std::string description;
  bool validated = false;
  bool validatedLeft = false;
  bool validatedRight = false;
  bool validatedAcross = false;
  bool validatedToOpositeSides = false;
  std::thread* t = nullptr;
  
  ComputeImageTask *left, *right, *set;
};
auxil::ConfigFileReader CheckCoveringRelation::reader = auxil::ConfigFileReader("data/config.dat");

// ###############################################

int main(int argc,char* argv[]){
  cout.precision(16);

  try{
    int numberOfThreads = thread::hardware_concurrency();
    if(argc>1)  numberOfThreads = atoi(argv[1]);
    threading::ThreadPool pool(numberOfThreads);
    initHSets();
    
    ostringstream out;
    out << "\nThreads = " << numberOfThreads << endl; 
    out << "Dimension = " << g_m << endl;
    out << "Order of PdeSolver = " << g_order << endl;
    

    /// Start clock
    auto tictac = tic(std::cout,out.str().c_str());

      /// An array of tasks that validate covering relations
      vector<CheckCoveringRelation*> covrelList;

      /// Create and run tasks for validation of covering relations.
      
      /// Heteroclinic chain a2->a1
      HSet* h2 = new HSet(*N2);
      h2->iterations = 4;
      covrelList.push_back(new CheckCoveringRelation(pool, h2,      h21[0], "N2->N0_21"));
      covrelList.push_back(new CheckCoveringRelation(pool, h21[0],  h21[1], "N0_21->N1_21"));
      covrelList.push_back(new CheckCoveringRelation(pool, h21[1],  h21[2], "N1_21->N2_21"));
      covrelList.push_back(new CheckCoveringRelation(pool, h21[2],  h21[3], "N2_21->N3_21"));
      covrelList.push_back(new CheckCoveringRelation(pool, h21[3],  h21[4], "N3_21->N4_21"));
      covrelList.push_back(new CheckCoveringRelation(pool, h21[4],  h21[5], "N4_21->N5_21"));
      covrelList.push_back(new CheckCoveringRelation(pool, h21[5],  h21[6], "N5_21->N6_21"));
      covrelList.push_back(new CheckCoveringRelation(pool, h21[6],  h21[7], "N6_21->N7_21"));
      covrelList.push_back(new CheckCoveringRelation(pool, h21[7],  h21[8], "N7_21->N8_21"));
      covrelList.push_back(new CheckCoveringRelation(pool, h21[8],  h21[9], "N8_21->N9_21"));       
      covrelList.push_back(new CheckCoveringRelation(pool, h21[9],  h21[10],"N9_21->N10_21"));
      covrelList.push_back(new CheckCoveringRelation(pool, h21[10], N1,     "N10_21->N1"));

      /// Heteroclinic chain a1->a2
      HSet* h1 = new HSet(*N1);
      h1->iterations = 3;
      covrelList.push_back(new CheckCoveringRelation(pool, h1,      h12[0], "N1->N0_12")); 
      covrelList.push_back(new CheckCoveringRelation(pool, h12[0],  h12[1], "N0_12->N1_12"));
      covrelList.push_back(new CheckCoveringRelation(pool, h12[1],  h12[2], "N1_12->N2_12"));
      covrelList.push_back(new CheckCoveringRelation(pool, h12[2],  h12[3], "N2_12->N3_12"));
      covrelList.push_back(new CheckCoveringRelation(pool, h12[3],  h12[4], "N3_12->N4_12"));
      covrelList.push_back(new CheckCoveringRelation(pool, h12[4],  h12[5], "N4_12->N5_12"));
      covrelList.push_back(new CheckCoveringRelation(pool, h12[5],  h12[6], "N5_12->N6_12"));
      covrelList.push_back(new CheckCoveringRelation(pool, h12[6],  h12[7], "N6_12->N7_12"));
      covrelList.push_back(new CheckCoveringRelation(pool, h12[7],  h12[8], "N7_12->N8_12"));
      covrelList.push_back(new CheckCoveringRelation(pool, h12[8],  h12[9], "N8_12->N9_12"));
      covrelList.push_back(new CheckCoveringRelation(pool, h12[9],  h12[10],"N9_12->N10_12"));
      covrelList.push_back(new CheckCoveringRelation(pool, h12[10], N2,     "N10_12->N2"));

      /// Periodic loops
      covrelList.push_back(new CheckCoveringRelation(pool, N1, M,  "N1->M"));
      covrelList.push_back(new CheckCoveringRelation(pool, M,  N1, "M->N1"));
      
      covrelList.push_back(new CheckCoveringRelation(pool, N2, K1, "N2->K1"));
      covrelList.push_back(new CheckCoveringRelation(pool, K1, K2, "K1->K2"));
      covrelList.push_back(new CheckCoveringRelation(pool, K2, K3, "K2->K3"));
      covrelList.push_back(new CheckCoveringRelation(pool, K3, N2, "K3->N2"));

      /// Wait for tasks 
      for(CheckCoveringRelation* t : covrelList) t->join();

      /// Print the final result
      bool validated = true;
      for(CheckCoveringRelation* t : covrelList)
        validated = validated and t->result();
      cout << "\nFinal result: " << boolalpha << validated << endl;

    /// Stop clock and print time of computation
    tac(std::cout,tictac);
    
    /// Release memory
    for(CheckCoveringRelation* t : covrelList) delete t;
    for(auto* p : h12) delete p;
    for(auto* p : h21) delete p;
    delete h1, h2, N1, N2, M, K1, K2, K3;
    
  }catch(exception& e){
    cout << "Exception: " << e.what() << endl;
  }
}
