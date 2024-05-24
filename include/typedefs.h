/**
 * This header files provides names of data structures and classes used in the program.
 * All types are defined in the CAPD as template classes. 
 * Here we just provide tempalte parameters and skip namespaces
 */ 
#ifndef DW__MAIN_H__
#define DW__MAIN_H__

#define LOGGER(out,x) out << "====>" << (#x) << "=" <<  (x) << "\n" << std::flush;
#define LOGGERN(out,x) out << "====>" << (#x) << "=" <<  (x) <<s "\n\n" << std::flush;

#include "capd/capdlib.h"
#include "capd/poincare/PoincareMap.hpp"
#include "capd/poincare/TimeMap.hpp"
#include "capd/poincare/SectionDerivativesEnclosure.h"
#include "capd/pdes/GeometricBound.h"
#include "capd/pdes/DissipativeVectorField.h"
#include "capd/pdes/PdeCurve.h"
#include "capd/pdes/OneDimKSSineVectorField.h"
#include "capd/pdes/PdeSolver.h"
#include "capd/pdes/C1DoubletonSetGeometricTail.h"
#include "capd/pdes/C0DoubletonSetGeometricTail.h"
#include "capd/pdes/C0HODoubletonSetGeometricTail.h"
#include "capd/pdes/PdeCoordinateSection.h"
#include "capd/pdes/PdeAffineSection.h"

using namespace capd;
using namespace std;

typedef capd::pdes::GeometricBound<capd::interval> GeometricBound;
using capd::pdes::OneDimKSSineVectorField;

template class capd::poincare::SectionDerivativesEnclosure<IMatrix>;
typedef capd::pdes::PdeSolver< GeometricBound > PdeSolver;
typedef capd::poincare::TimeMap<PdeSolver> PdeTimeMap;

typedef capd::pdes::PdeAbstractSection<GeometricBound,IMatrix> PdeAbstractSection;
typedef capd::pdes::PdeCoordinateSection<GeometricBound,IMatrix> PdeCoordinateSection;
typedef capd::pdes::PdeAffineSection<GeometricBound,IMatrix> PdeAffineSection;
typedef capd::poincare::PoincareMap<PdeSolver,PdeAbstractSection> PdePoincareMap;

typedef capd::dynset::FactorReorganization<capd::dynset::PartialQRWithPivoting<5> >  PartialQR;
typedef capd::dynset::FactorReorganization<capd::dynset::PartialQRWithPivoting<8> >  C1PartialQR;
typedef capd::pdes::C0DoubletonSetGeometricTail<capd::dynset::C0DoubletonSet<capd::IMatrix,PartialQR>> C0DoubletonSetGeometricTail;
typedef capd::pdes::C1DoubletonSetGeometricTail<C1PartialQR> C1DoubletonSetGeometricTail;

#endif

