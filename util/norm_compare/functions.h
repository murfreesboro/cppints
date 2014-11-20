/**
 * \file   functions.h
 * \brief  This file contains the variety mathematical functions 
 * \author Fenglai Liu and Jing Kong
 * \date   June, 2011
 * \note
 * Oct. 2011: add the fm function for the integral codes
 */
#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include "libgen.h"
#include <boost/math/special_functions/binomial.hpp>  // special functions in math
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/gamma.hpp>
using namespace boost::math;

namespace functions {

	/**
	 * factorial of n!
	 */
	inline Double getFactorial(const Double& n) {
		return factorial<Double>(n);
	};

	/**
	 * BionomialCoe returns C^{l}_{n}
	 */
	inline Double bionomialCoe(const Double& n, const Double& l) {
		return binomial_coefficient<Double>(n, l);
	};

	/**
	 * DoubleFac returns N!!
	 */
	inline Double doubleFac(const Double& n) {
		if (fabs(n)<THRESHOLD_MATH || fabs(n+ONE)<THRESHOLD_MATH) {
			return ONE; //case for n = -1 and n = 0
		}else{
			return double_factorial<Double>(n);
		}
	}; 


	// max and min function are provided in the STD library


	/**
	 * GammaFunInt calculate the Gamma integral below:
	 * \int_-infinity^+infinity  x^n*exp(-alpha*x^2) dx
	 */
   inline Double gammaFunInt(const Int& n, const Double& alpha) {
		if (n%2 == 0) {
			if (n == 0) {
				return sqrt(PI/alpha);
			} else {
				Double x  = sqrt(PI/alpha);
				Double ta = TWO*alpha;
				for(Int i=n/2; i>0; i--) {
					x = x*(TWO*i-1)/ta;
				}
				return x;
			}
		} else {
			return ZERO;
		}
	};

	/**
	 * This functions simply calculates the Fm(a,m) function, which is shown below:
	 * \f$ Fm(x,a,m) = \int^{x}_{0} e^{-at^{2}}t^{2m} dt \f$
	 * such function can be transformed finally into the incomplete Gamma function:
	 * \f$ Fm(x,a,m) = \frac{1}{2a^{m+1/2}} \int^{ax^{2}}_{0} e^{-x} x^{m-1/2} dx \f$
	 */
	inline Double fm(const Double& x,  const Double& a, const Int& m) {
		if (fabs(a)<THRESHOLD_MATH) return (pow(x,2*m+1)/(2*m+1));
		Double f12 = ONE/TWO;
		Double p   = ONE/(TWO*pow(a,m+f12));
		Double upperLimit = a*pow(x,TWO);
		return tgamma_lower(m+f12,upperLimit)*p;
	};

	/**
	 * Based on the fm function, we can get its polynomial form which is used in 
	 * the NAI:
	 * \f$ \int^{x}_{0} e^{-av^{2}} (x^{2}-v^{2})^{L^{'}}v^{2(L-2L^{'})} dv \f$
	 */
	inline Double fmNAI(const Double& x, const Double& a, const Int& Lp, const Int& L) {
		Double r  = ZERO;
		Double x2 = x*x;
		for(Int i=0;i<=Lp;i++){
			Double p = bionomialCoe(Lp,i)*pow(MINUS_ONE,i)*pow(x2,Lp-i);
			Int    e = i+L-2*Lp;
			r += p*fm(x,a,e);
		}
		return r;
	};

	/**
	 * integration for the expression that:
	 * \f$ \int^{1}_{0} (1-t^{2})^{L} dt \f$
	 */
	inline Double oneMinusT2L(const Int& L) {
		Double r = ZERO;
		for(Int i=0;i<=L;i++){
			Double p = bionomialCoe(L,i)*pow(MINUS_ONE,i);
			r += p/(2*i+1); 
		}
		return r;
	};

}

#endif
