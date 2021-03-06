#ifndef ROBUST_HPP
#define ROBUST_HPP 1

#include "iterativeMethod.hpp"
#include "bisection.hpp"
#include "newton.hpp"

namespace RootFinding {
  /*!
    \class Robust
    \brief Find the root of a non-linear equation using a robust method 
    \author Daniele A. Di Pietro
    \date 8-26-2005
  */
  class Robust : public IterativeMethod {
  public:
    //! Coarse method type
    typedef Bisection coarseT;
    //! Fine method type
    typedef Newton fineT;
    /** @name Constructors
     */
    //@{
    //! Constructor
    Robust(real, int, checkT = INCREMENT, real = 5e4);
    //@}
    /** @name Accessors
     */
    //@{
    //! Return the number of iterations of the coarse method
    inline int nitcoarse() { return M_coarse.nit(); }
    //! Return the number to iterations of the fine method
    inline int nitfine() { return M_fine.nit(); }
    //@}
    /** @name Methods
     */
    //@{
    //! Apply the method
    real apply(real, real, fctptr, fctptr);
    //@}
    //! Display the results
    friend std::ostream& operator<<(std::ostream&, Robust&);
  private:
	//! The ratio between coarse and fine method tolerance
	real M_cfratio;
    //! The coarse method
    coarseT M_coarse;
    //! The fine method
    fineT M_fine;
  };
}

#endif
