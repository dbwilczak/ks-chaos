# KS-CHAOS: suplement to the article 
### Daniel Wilczak and Piotr Zgliczy&#324;ski, Validated integration of variational equations for a class of dissipative PDEs

## Table of content

- [Introduction](#introduction)
- [Numerical constants used in the proof](#constants)
- [Description of h-sets](#hsets)
- [List of h-sets that appear in the computation](#hsets)
- [The C++-17 programs](#program)
  - [Required resources](#requirements)
  - [Compilation of the program](#gcc)
  - [Existence of a stable symmetric periodic orbit](#stablepo)
  - [Existence of symbolic dynamics](#symdyn)
  - [Existence of connecting orbits](#conecondition)
  - [Source code of the program](#code)
  
## Introduction<a name="introduction"></a>

The Kuramoto-Sivashinsky PDE (referred in the sequel as the KS equation), after imposing odd and periodic boundary conditions, can be rewritten as an infinite-dimensional ODE

> $`a_k' = k^2(1-\nu k^2)a_k - k\sum_{n=1}^{k-1}a_na_{k-n}+ 2k\sum_{n=1}^\infty a_na_{k+n}`$,

where $`a_k, k\geq 1`$ are the Fourier coefficients of the solutions to the KS PDE. Symbolic dynamics for the above system is verified by means of the method of covering relations (see article for details).

We constructed two sequences of covering relations along numerically observed heteroclinic chains joining two selected  periodic solutions in both directions.
Then, using rigorous integrator for the KS equation, we were able to check assumptions of abstract theorems, which guarantee the existence of symbolic dynamics, connecting orbits and stable periodic orbits.

The purpose of this document is to provide
- an information about numerical constants and settings used in the computation
- definition of h-sets and Poincare sections used in the computation
- an instruction of how one can compile and run the program
  
## Numerical constants used in the proof <a name="constants"></a>
The following constants were used in our computation. Some of them are mentioned in the article.
- $`q = 3/2`$ - this is geometric decay ratio used to define infinite tails (see data structures <b>GBound</b> and <b>HSet</b> in the article)</li>
- $`m = 15`$ -number of leading coordinates of an h-set represented explicitly in the computer memory (see data structure <b>HSet</b> in the article)</li>

## Description of h-sets <a name="hsets"></a>
In order to reduce overestimation in rigorous computation of Poincare maps, each h-set is embedded into a different Poincare section, which is almost orthogonal to the vector field near the centre of this h-set (an approximate heteroclinic point).

Each h-set that appears in the computation is described by the following data:
- $`S\geq 0`$ - constant used to define geometrically decaying tail $`|a_k|\leq Sq^{-k}`$ for $`k\geq m`$
- $`a = (a_1,\ldots,a_m)`$ - this is the centre of an h-set (an approximate periodic or heteroclinic point), represented as&nbsp; $`m-`$dimensional vector with $`a_1=0`$<br/>
    *Note that approximate periodic and heteroclinic points have been found on the section* 
    > $`\Pi = \{ a_1=0 \ \text{and}\ a_1'>0\}`$
- $`E`$  - $`m\times m`$ - almost orthogonal matrix (up to floating point accuracy) used to define Poincare section. First column of the matrix is close to flow direction at $`a`$, while the remaining columns span $`(m-1)`$-dimensional subspace of the section. 
Formally, the section is given by 
> $`\Pi = \{ (x_k)_{k=1}^\infty | \pi_1 E*(\pi_{\leq m}x-a) = 0 \}`$

- $`B`$ - $`m\times m`$ - shape matrix of an h-set. First row and first column of this matrix is set to $`(1,0,...,0)`$ just to simplify further notation. 
The remaining $`(m-1)\times(m-1)`$ block defines an affine change of coordinates.

- $`r = (r_1,\ldots,r_m`$ - an interval vector centred at zero, with $`r_1=[0,0]`$. It defines size of an h-set in the main $`(m-1)`$ coordinates on the section

The above data structure $`[S,a,E,B,r]`$ defines an h-set $`H`$ by
> $`\pi_{\leq m}H = a + E*B*r`$<br/>
> $`H_k = [-S,S]q^{-k}`$ for $`k>m`$
  

## List of h-sets that appear in the computation <a name="listhsets"></a>

Below we give list of links to&nbsp; 28 text files. Each file contains a record $`[S,a,E,B,r]`$ of data, as described above.
File names correspond to notation used in the article.

### H-sets at periodic points

[**N1.txt**](hsets/N1.txt) - the h-set centred at an approximate period two point $`a^1`$ of Poincare map on the section $`a_1=0`$

[**N2.txt**](hsets/N2.txt) - the h-set centred at an approximate period four point $`a^2`$ of Poincare map on the section $`a_1=0`$

[**M.txt**](hsets/M.txt) - the h-set centred at an approximate period two point $`P(a^1)`$ of Poincare map on the section $`a_1=0`$

[**K1.txt**](hsets/K1.txt) - the h-set centred at an approximate period four point $`P(a^2)`$ of Poincare map on the section $`a_1=0`$

[**K2.txt**](hsets/K2.txt) - the h-set centred at an approximate period four point $`P^2(a^2)`$ of Poincare map on the section $`a_1=0`$

[**K3.txt**](hsets/K3.txt) - the h-set centred at an approximate period four point $`P^3(a^2)`$ of Poincare map on the section $`a_1=0`$

### H-sets along $`a^1\to a^2`$ heteroclinic chain
  - [**N0_12.txt**](hsets/N0_12.txt) - the h-set $`N^0_{1\to2}`$
  - [**N1_12.txt**](hsets/N1_12.txt) - the h-set $`N^1_{1\to2}`$
  - [**N2_12.txt**](hsets/N2_12.txt) - the h-set $`N^2_{1\to2}`$
  - [**N3_12.txt**](hsets/N3_12.txt) - the h-set $`N^3_{1\to2}`$
  - [**N4_12.txt**](hsets/N4_12.txt) - the h-set $`N^4_{1\to2}`$
  - [**N5_12.txt**](hsets/N5_12.txt) - the h-set $`N^5_{1\to2}`$
  - [**N6_12.txt**](hsets/N6_12.txt) - the h-set $`N^6_{1\to2}`$
  - [**N7_12.txt**](hsets/N7_12.txt) - the h-set $`N^7_{1\to2}`$
  - [**N8_12.txt**](hsets/N8_12.txt) - the h-set $`N^8_{1\to2}`$
  - [**N9_12.txt**](hsets/N9_12.txt) - the h-set $`N^9_{1\to2}`$
  - [**N10_12.txt**](hsets/N10_12.txt) - the h-set $`N^{10}_{1\to2}`$

### H-sets along $`a^2\to a^1`$ heteroclinic chain
  - [**N0_21.txt**](hsets/N0_21.txt) - the h-set $`N^0_{2\to1}`$
  - [**N1_21.txt**](hsets/N1_21.txt) - the h-set $`N^1_{2\to1}`$
  - [**N2_21.txt**](hsets/N2_21.txt) - the h-set $`N^2_{2\to1}`$
  - [**N3_21.txt**](hsets/N3_21.txt) - the h-set $`N^3_{2\to1}`$
  - [**N4_21.txt**](hsets/N4_21.txt) - the h-set $`N^4_{2\to1}`$
  - [**N5_21.txt**](hsets/N5_21.txt) - the h-set $`N^5_{2\to1}`$
  - [**N6_21.txt**](hsets/N6_21.txt) - the h-set $`N^6_{2\to1}`$
  - [**N7_21.txt**](hsets/N7_21.txt) - the h-set $`N^7_{2\to1}`$
  - [**N8_21.txt**](hsets/N8_21.txt) - the h-set $`N^8_{2\to1}`$
  - [**N9_21.txt**](hsets/N9_21.txt) - the h-set $`N^9_{2\to1}`$
  - [**N10_21.txt**](hsets/N10_21.txt) - the h-set $`N^{10}_{2\to1}`$

## The C++-17 programs <a name="program"></a>

### Required resources <a name="requirements"></a>
- The programs are parallelized by means of the C++17 threading.
- Each of the above programs should not exceed 100 minutes of computation on a modern personal computer having 16 cores.
- None of the programs require large RAM memory -&nbsp; **1GB** is more than enough.

### Compilation of the program <a name="gcc"></a>
The program can be compiled and run on any computer running linux or under ubuntu console for MS Windows (WSL2).

The program is based on the [**CAPD library**](http://capd.ii.uj.edu.pl) which provides in particular
- interval arithmetics
- data structures and algorithms of linear algebra
- data structures and algorithms for automatic differentiation
- data structures and algorithms for rigorous integration of ODEs
- data structures and algorithms for rigorous computation of Poincare maps

In an empty directory: clone source code of the CAPD library and compile it
- clone source code of the CAPD library and compile it  
    ```
      git clone https://github.com/CAPDGroup/CAPD
      cd CAPD && mkdir build && cd build 
      cmake .. -DCAPD_ENABLE_MULTIPRECISION=false -DCAPD_BUILD_EXAMPLES=false
      make -j
      cd ../../
    ```
- after succesfull compilation of the CAPD library clone source code of the main programs and compile them
    ```
      git clone https://github.com/dbwilczak/ks-chaos
      cd ks-chaos && mkdir dep && mkdir obj
      make -j
    ```

- After succesfull compilation three executables are created:
  - **check_stable_po**
  - **check_cone_conditions** and 
  - **check_covering_relations**

### Existence of a stable symmetric periodic orbit<a name="stablepo"></a>
The executable ```check_stable_po``` is built from the following source
> [**progs/check_stable_po.cpp**](progs/check_stable_po.cpp)


- The program runs within few minutes in a single thread.
- It does not create output files.
- During the computation it prints to the standard output computed approximate periodic orbit, bounds from the validation of the existence of this periodic orbit and, finally, bounds on the norm of the derivative of Poincare map.

Here are two sample outputs of the program (as reported in the article):
- [**output/outStablePeriodicOrbitDim15**](output/outStablePeriodicOrbitDim15) - validation in effective dimension 15 and larger tolerance $`10^{-7}`$ of the PDE solver<br/>
- [**output/outStablePeriodicOrbitDim18**](output/outStablePeriodicOrbitDim18) - validation in effective dimension 18 and smaller tolerance $`10^{-11}`$ of the PDE solver<br/>

By default the program runs computation in dimension 15 (much faster) - see <b>main</b> function in [**progs/check_stable_po.cpp**](progs/check_stable_po.cpp) for details.



### <a name="symdyn"></a>Existence of symbolic dynamics.
The executable ```check_covering_relations``` is built from the following source
>  [**progs/check_covering_relations.cpp**](progs/check_covering_relations.cpp)

The program creates and runs **3x30=90 concurrent tasks** for verification of the existence of the following covering relations.

-  $`N^1\Longrightarrow M \Longrightarrow N^1`$
-  $`N^2\Longrightarrow K_1\Longrightarrow K_2\Longrightarrow K_3\Longrightarrow N^2`$
-  $`N^1\Longrightarrow N^0_{1\to2}\Longrightarrow => N^1_{1\to2}\Longrightarrow\cdots\Longrightarrow N^{10}_{1\to2}\Longrightarrow N^1`$
-  $`N^2\Longrightarrow N^0_{2\to1}\Longrightarrow => N^1_{2\to1}\Longrightarrow\cdots\Longrightarrow N^{10}_{2\to1}\Longrightarrow N^2`$

**Comments**
- The program automatically recognizes the number of available cores and runs all 90 tasks in a thread pool.
- The program should finish computation within 100 minutes on a modern CPU with 16 cores.
- The program does not create output files.
- Computed bounds and some logs from the verification of the covering relations are printed to the standard output during the computation. The order of output data is not determined (depends on scheduling of tasks).

**Here is a sample output of the program**
>  [**output/outCoveringRelations**](output/outCoveringRelations)



### <a name="conecondition"></a>Existence of connecting orbits
The executable ```check_cone_conditions``` is built from the following source
>  [**progs/check_cone_conditions.cpp**](progs/check_cone_conditions.cpp)

The program creates and runs **12 concurrent tasks** tasks for verification of cone conditions for covering relations
-  $`N^1\Longrightarrow M \Longrightarrow N^1`$
-  $`N^2\Longrightarrow K_1\Longrightarrow K_2\Longrightarrow K_3\Longrightarrow N^2`$</b></tt></blockquote>

**Comments**
- The program automatically recognizes the number of available cores and runs all 12 tasks in a thread pool.
- The program should finish computation within an hour on a modern CPU with at least 12 cores.
- Computed bounds and some logs from the verification of cone conditions  are printed to the standard output during the computation. The order of output data is not determined (depends on scheduling of tasks).
- The program creates output files, which contain computed bounds on the derivative of Poincare map expressed in coordinates of the corresponding h-sets - as described in the article. They are located in the 
  [**output**](output) directory:

  - [**output/N1.dat**](output/N1.dat)
  - [**output/N2.dat**](output/N2.dat)
  - [**output/M.dat**](output/M.dat)
  - [**output/K1.dat**](output/K1.dat)
  - [**output/K2.dat**](output/K2.dat)
  - [**output/K3.dat**](output/K3.dat)

**Here is a sample output of the program:**
> [**output/outConeConditions**](output/outConeConditions)



## <a name="code"></a>Source code of the program.
Here we give an overview on the code specific for the computer-assisted proof of symbolic dynamics and connecting orbits in the KS equations. 
The program consists ot three executables. They use several algorithms and auxiliary routines implemented in the following files:
- [**include/typedefs.h**](include/typedefs.h)<br/>
    main header file which includes the CAPD library and provides type definitions for several template types
- [**include/tictac.h**](include/tictac.h)<br/>
    an auxiliary file which contains two global functions used to measure and print time of computation
- [**include/projections.h**](include/projections.h)<br/>
    an auxiliary file which implements short routines - projections and embeddings from/to Poincare sections
- [**include/HSet.h**](include/HSet.h)<br/>
    data structure that represents an h-set
- [**include/LDKSPoincareMap.h**](include/KSPoincareMap.h)<br/>
    a wrapper for the classes **LDOdeSolver* and **LDPoincareMap** from the [**CAPD library**](http://capd.ii.uj.edu.pl). This class is used for nonrigrous integration (in long double precision), for instance to find a candidate of a stable periodic orbit by simple iteration.
- [**include/IKSPoincareMap.h**](include/IPoincareMap.h)<br/> 
    a wrapper for the classes **PdeSolver** and **IPoincareMap** from the [**CAPD library**](http://capd.ii.uj.edu.pl). This class is used in rigorous integration. It provides in particular algorithms for computation of Poincare maps and derivatives in given coordinate systems in domain and codomain.
- [**include/ComputeImageTask.h**](include/ComputeImageTask.h)<br/>
    a task that can be executed in a thread pool and which computes a bound on the Poincare map between sections on a given domain. It uses class **IPoincareMap** to integrate the equation.
- [**include/ComputeDerivativeTask.h**](include/ComputeDerivativeTask.h)<br/>
    a task that can be executed in a thread pool and which computes a bound on the derivative of Poincare map between sections on a given domain. It uses class **IPoincareMap** to integrate the variational equation.

Based on the above general classes and routines, the three executables are written that realize computer assisted proofs of the existence of a stable, symmetric periodic orbit, symbolic dynamics and connecting orbits between periodic orbits
- [**progs/check_stable_po.cpp**](progs/check_stable_po.cpp)
- [**progs/check_covering_relations.cpp**](progs/check_covering_relations.cpp)
- [**progs/check_cone_conditions.cpp**](progs/check_cone_conditions.cpp)
