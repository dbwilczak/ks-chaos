#ifndef DW__COMPUTE_IMAGE_TASK_H_
#define DW__COMPUTE_IMAGE_TASK_H_

/** 
  * ComputeImageTask is responsible for asynchronous computation of Poincare map on a given domain.
  * The image of Poincare map is represented in an affine coordinate system.
  * 
  * Computation is executed asychrnonously in a thread from the thread pool.
  * After receiving the result the status is changed to Completed 
  * and the result os stored in the member 'result'.
  * If the rigorous integration fails, the status is changed to Failure.
  * 
  * This class implements an abstract class Task from the CAPD library. 
  * Objects of this type can be executed as tasks in a ThreadPool.
  * The method 'run' is executed asychrnonously in a thread from the thread pool.
  * 
*/
#include "capd/threading/ThreadPool.h"
#include "include/IKSPoincareMap.h"

class ComputeImageTask: public threading::Task{
public:
  /// Task input data
  /// x+C*r - set on which we compute Poincare map
  const IVector x; const IMatrix C; const IVector r;
  /// y - > Q*(P^iterations(y)-u) - coordinate system in which result is represented
  const IMatrix Q; const IVector u;
  
  const int iterations;         /// number of iterations of Poincare map
  const IAffineSection section; /// Poincare section at u
  const double S;               /// constant describing geometric tail S*q^{-k}
    
  /// Control data
  const IKSPoincareMapFactory& pmFactory;
  
  /// Task results
  interval returnTime;           /// bound for return time
  GeometricBound result; /// bound for the result Q*(P^iterations(x+C*r)-u)
  /// Status 'Completed' means, that the rutine which computes Poincare map returned the result without throwing an exception - see method 'run' below.
  enum Status { New, Completed, Failure };
  Status status = Status::New;

  ComputeImageTask(const ComputeImageTask& ) = delete;

  /// Constructor
  ComputeImageTask(
    IVector _x, IMatrix _C, IVector _r, 
    IMatrix _Q, IVector _u, 
    double _S,
    int _iterations, 
    IAffineSection _section, 
    const IKSPoincareMapFactory& _pmFactory
  ) :
    x(_x), C(_C), r(_r), 
    Q(_Q), u(_u), 
    S(_S),
    iterations(_iterations), 
    section(_section),
    result(_x.dimension(),0.,_pmFactory.q),
    pmFactory(_pmFactory)
  {}

  /// This method is executed asynchronously by a thread from the ThreadPool
  virtual void run(unsigned){
    std::unique_ptr<IKSPoincareMap> pm;
    try{
      pm.reset(pmFactory.createPoincareMap(section));

      int D = x.dimension();
      IVector x0(D), dr(D);
      split(r,x0,dr);
      /// here we call the routine which computes Poincare map. 
      this->result = pm->computeImage(x+C*x0,C,dr,S,Q,u,iterations,returnTime);
      
      /// we mark the status as completed.
      /// If the routine was able to integrate the set of initial conditions until its intersection with Poincare section (no exception caught)
      /// then we set status of computation as Completed. We will check if the computed bound satisfies required inequalities if all tasks finish computation.
      /// See class below CheckCoveringRelation
      this->status = Status::Completed;
    }catch(std::exception& e){
      std::cout << "Exception: " << e.what() << std::endl;
      status = Status::Failure;
    }
  }
};


#endif
