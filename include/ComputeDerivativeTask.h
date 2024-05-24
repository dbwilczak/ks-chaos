#ifndef DW__COMPUTE_DERIVATIVE_TASK_H__
#define DW__COMPUTE_DERIVATIVE_TASK_H__

#include "capd/threading/ThreadPool.h"
#include "include/IKSPoincareMap.h"

/// This class is responsible for asynchronous computation of derivative of Poincare map.
/// It runs in a spearate thread from a thread pool, that is the method 'run' is called async.
/// After the task completes computation, the task changes its status from New to Completed
/// and the results is stored in members Dp, Py, Pyx.
/// If the computation fails, the status is set to Failure;

struct ComputeDerivativeTask : public threading::Task{
  /// Control data
  const IKSPoincareMapFactory& pmFactory;
  int iterations; 
  
  /// Input data, the set is represented as (x+A*r,tail), where tail is of the form |a_k|<=S*q^{-k}
  GeometricBound x;
  IMatrix A; 
  IVector r;
  IMatrix Q, invQ;
  
  /// Output in coordinate system A, that is
  IMatrix Jxx;    /// Jxx = invQ*DxxP*Q
  IVector Jxy;
  interval Jyy, Jyx;

  interval returnTime = 0.;
  enum Status { New, Completed, Failure };
  Status status = Status::New;
  
  /// Constructor saves initial condition for integration
  ComputeDerivativeTask(const IKSPoincareMapFactory& _pmFactory, GeometricBound _x, IMatrix _A, IVector _r, IMatrix Q, IMatrix _invQ, int _iterations) 
    : pmFactory(_pmFactory), iterations(_iterations),
      x(_x), A(_A), r(_r), Q(Q), invQ(_invQ),
      Jxx(A.numberOfRows(),A.numberOfColumns()), Jxy(A.numberOfColumns())
  {}

  /// computes invA*D_{xx}P(x0+A*r)*A
  /// and norms of blocks DxyP, DyxP, DyyP
  void run(unsigned threadId){
    std::unique_ptr<IKSPoincareMap> pm;

    try{
      /// create an instance of IKSPoincare map
      pm.reset(pmFactory.createPoincareMap());

      /// here we call the routine which computes derivative of Poincare map. 
      Jxx = pm->computeDerivative(x,A,r,iterations,Q,invQ,returnTime,Jxy,Jyx,Jyy);
      
      /// we mark the status as completed.
      /// If the routine was able to integrate the set of initial conditions until its intersection with Poincare section (no exception caught)
      /// then we set status of computation as Completed. We will check if the computed bound satisfies required inequalities if all tasks finish computation.
      this->status = Status::Completed;
    } catch(std::exception& e){
        cout << "Exception caught:\n" << e.what() << endl;
        status = Status::Failure;
    }
  }    
}; // struct ComputeDerivativeTask

#endif
