// @(#)root/mathmore:$Id$
// Author: L. Moneta  03/2011

 /**********************************************************************
  *                                                                    *
  * Copyright (c) 2004 ROOT Foundation,  CERN/PH-SFT                   *
  *                                                                    *
  * This library is free software; you can redistribute it and/or      *
  * modify it under the terms of the GNU General Public License        *
  * as published by the Free Software Foundation; either version 2     *
  * of the License, or (at your option) any later version.             *
  *                                                                    *
  * This library is distributed in the hope that it will be useful,    *
  * but WITHOUT ANY WARRANTY; without even the implied warranty of     *
  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU   *
  * General Public License for more details.                           *
  *                                                                    *
  * You should have received a copy of the GNU General Public License  *
  * along with this library (see file COPYING); if not, write          *
  * to the Free Software Foundation, Inc., 59 Temple Place, Suite      *
  * 330, Boston, MA 02111-1307 USA, or contact the author.             *
  *                                                                    *
  **********************************************************************/

// Header file for class GSLMultiRootFinder
//

#ifndef ROOT_Math_GSLMultiRootFinder
#define ROOT_Math_GSLMultiRootFinder



#include "Math/IFunction.h"

#include "Math/WrappedFunction.h"

#include <vector>
#include <utility>
#include <iostream>

namespace ROOT {
namespace Math {


   class GSLMultiRootBaseSolver;

   /** @defgroup MultiRoot  Multidimensional ROOT finding 
       Classes for finding the roots of a multi-dimensional system. 
       @ingroup NumAlgo
   */

  /**
     Class for  Multidimensional root finding algorithms bassed on GSL. This class is used to solve a
     non-linear system of equations:

     f1(x1,....xn) = 0
     f2(x1,....xn) = 0
     ..................
     fn(x1,....xn) = 0

     See the GSL <A HREF="http://www.gnu.org/software/gsl/manual/html_node/Multidimensional-Root_002dFinding.html"> online manual</A> for
     information on the GSL MultiRoot finding algorithms

     The available GSL algorithms require the derivatives of the supplied functions or not (they are
     computed internally by GSL). In the first case the user needs to provide a list of multidimensional functions implementing the
     gradient interface (ROOT::Math::IMultiGradFunction) while in the second case it is enough to supply a list of
     functions impelmenting the ROOT::Math::IMultiGenFunction interface.
     The available algorithms requiring derivatives (see also the GSL
     <A HREF="http://www.gnu.org/software/gsl/manual/html_node/Algorithms-using-Derivatives.html">documentation</A> )
     are the followings:
     <ul>
         <li><tt>ROOT::Math::GSLMultiRootFinder::kHybridSJ</tt>  with name <i>"HybridSJ"</i>: modified Powell's hybrid
     method as implemented in HYBRJ in MINPACK
         <li><tt>ROOT::Math::GSLMultiRootFinder::kHybridJ</tt>  with name <i>"HybridJ"</i>: unscaled version of the
     previous algorithm</li>
         <li><tt>ROOT::Math::GSLMultiRootFinder::kNewton</tt>  with name <i>"Newton"</i>: Newton method </li>
         <li><tt>ROOT::Math::GSLMultiRootFinder::kGNewton</tt>  with name <i>"GNewton"</i>: modified Newton method </li>
     </ul>
     The algorithms without derivatives (see also the GSL
     <A HREF="http://www.gnu.org/software/gsl/manual/html_node/Algorithms-without-Derivatives.html">documentation</A> )
     are the followings:
     <ul>
         <li><tt>ROOT::Math::GSLMultiRootFinder::kHybridS</tt>  with name <i>"HybridS"</i>: same as HybridSJ but using
     finate difference approximation for the derivatives</li>
         <li><tt>ROOT::Math::GSLMultiRootFinder::kHybrid</tt>  with name <i>"Hybrid"</i>: unscaled version of the
     previous algorithm</li>
         <li><tt>ROOT::Math::GSLMultiRootFinder::kDNewton</tt>  with name <i>"DNewton"</i>: discrete Newton algorithm </li>
         <li><tt>ROOT::Math::GSLMultiRootFinder::kBroyden</tt>  with name <i>"Broyden"</i>: Broyden algorithm </li>
     </ul>

     @ingroup MultiRoot
  */


 class GSLMultiRootFinder {

 public:

   /**
      enumeration specifying the types of GSL multi root finders
      requiring the derivatives

   */
    enum EDerivType {
       kHybridSJ,
       kHybridJ,
       kNewton,
       kGNewton
    };
    /**
       enumeration specifying the types of GSL multi root finders
       which do not require the derivatives

    */
    enum EType {
       kHybridS,
       kHybrid,
       kDNewton,
       kBroyden
    };



    /// create a multi-root finder based on an algorithm not requiring function derivative
    GSLMultiRootFinder(EType type);

    /// create a multi-root finder based on an algorithm requiring function derivative
    GSLMultiRootFinder(EDerivType type);

    /*
      create a multi-root finder using a string.
      The names are those defined in the GSL manuals
      after having remived the GSL prefix (gsl_multiroot_fsolver).
      Default algorithm  is "hybrids" (without derivative).
    */
    GSLMultiRootFinder(const char * name = 0);

    /// destructor
    virtual ~GSLMultiRootFinder();

 private:
    // usually copying is non trivial, so we make this unaccessible
    GSLMultiRootFinder(const GSLMultiRootFinder &);
    GSLMultiRootFinder & operator = (const GSLMultiRootFinder &);

 public:

    /// set the type for an algorithm without derivatives
    void SetType(EType type) {
       fType = type; fUseDerivAlgo = false;
    }

    /// set the type of algorithm using derivatives
    void SetType(EDerivType type) {
       fType = type; fUseDerivAlgo = true;
    }

    /// set the type using a string
    void SetType(const char * name);

    /*
       add the list of functions f1(x1,..xn),...fn(x1,...xn). The list must contain pointers of
       ROOT::Math::IMultiGenFunctions. The method requires the
       the begin and end of the list iterator.
       The list can be any stl container or a simple array of  ROOT::Math::IMultiGenFunctions* or
       whatever implementing an iterator.
       If using a derivative type algorithm the function pointers must implement the
       ROOT::Math::IMultiGradFunction interface
    */
    template<class FuncIterator>
    bool SetFunctionList( FuncIterator begin, FuncIterator end) {
       bool ret = true;
       for (FuncIterator itr = begin; itr != end; ++itr) {
          const ROOT::Math::IMultiGenFunction * f = *itr;
          // Using bitwise operator &= require the operand to be a bool
          // to have the intended effect here.
          ret &= (AddFunction( *f) != 0);
       }
       return ret;
    }

    /*
      add (set) a single function fi(x1,...xn) which is part of the system of
       specifying the begin and end of the iterator.
       If using a derivative type algorithm the function must implement the
       ROOT::Math::IMultiGradFunction interface
       Return the current number of function in the list and 0 if failed to add the function
     */
    int AddFunction( const ROOT::Math::IMultiGenFunction & func);

    /// same method as before but using any function implementing
    /// the operator(), so can be wrapped in a IMultiGenFunction interface
    template <class Function>
    int AddFunction( Function & f, int ndim) {
       // no need to care about lifetime of wfunc. It will be cloned inside AddFunction
       WrappedMultiFunction<Function &> wfunc(f, ndim);
       return AddFunction(wfunc);
    }

    /**
       return the number of sunctions set in the class.
       The number must be equal to the dimension of the functions
     */
    unsigned  int Dim() const { return fFunctions.size(); }

    /// clear list of functions
    void Clear();

    /// return the root X values solving the system
    const double * X() const;

    /// return the function values f(X) solving the system
    /// i.e. they must be close to zero at the solution
    const double * FVal() const;

    /// return the last step size
    const double * Dx() const;


    /**
       Find the root starting from the point X;
       Use the number of iteration and tolerance if given otherwise use
       default parameter values which can be defined by
       the static method SetDefault...
    */
    bool Solve(const double * x,  int maxIter = 0, double absTol = 0, double relTol = 0);

    /// Return number of iterations
    int Iterations() const {
       return fIter;
    }

    /// Return the status of last root finding
    int Status() const { return fStatus; }

    /// Return the algorithm name used for solving
    /// Note the name is available only after having called solved
    /// Otherwise an empyty string is returned
    const char * Name() const;

    /*
       set print level
       level = 0  quiet (no messages print)
             = 1  print only the result
             = 3  max debug. Print result at each iteration
    */
    void SetPrintLevel(int level) { fPrintLevel = level; }

    /// return the print level
    int PrintLevel() const { return fPrintLevel; }


    //-- static methods to set configurations

    /// set tolerance (absolute and relative)
    /// relative tolerance is only use to verify the convergence
    /// do it is a minor parameter
    static void SetDefaultTolerance(double abstol, double reltol = 0 );

    /// set maximum number of iterations
    static void SetDefaultMaxIterations(int maxiter);

    /// print iteration state
    void PrintState(std::ostream & os = std::cout);


 protected:

    // return type given a name
    std::pair<bool,int> GetType(const char * name);
    // clear list of functions
    void ClearFunctions();


 private:

    int fIter;           // current numer of iterations
    int fStatus;         // current status
    int fPrintLevel;     // print level

    // int fMaxIter;        // max number of iterations
    // double fAbsTolerance;  // absolute tolerance
    // double fRelTolerance;  // relative tolerance
    int fType;            // type of algorithm
    bool fUseDerivAlgo; // algorithm using derivative

    GSLMultiRootBaseSolver * fSolver;
    std::vector<ROOT::Math::IMultiGenFunction *> fFunctions;   //! transient Vector of the functions


 };

   // use typedef for most sensible name
   typedef GSLMultiRootFinder MultiRootFinder;

} // namespace Math
} // namespace ROOT


#endif /* ROOT_Math_GSLMultiRootFinder */
