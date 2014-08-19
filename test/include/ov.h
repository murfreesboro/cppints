#ifndef OV_H
#define OV_H
#include "libgen.h"

namespace ov {

	/**
	 * This function is used to calculate the fk in the form like below:
	 *
	 * \f$(X + a)^{l1}(X + b)^{l2} = \sum_{k=0}^{l1+l2}X^{k}fk(a,b,l1,l2)\f$
	 *
	 * This expression is used in the Gaussian Primitive Product Theorem
	 * A, B are atom centers, P is the new center created in the theorem
	 * l are the original angular momentums for the two primitives of 1 and 2
	 * See the document of Gaussian Primitive Product Theorem for more information
	 * \param k   order for angular momentum of x 
	 * \param l1  exponent for PAx
	 * \param l2  exponent for PBx
	 * \param x1  PAx
	 * \param x2  PBx
	 * \return    fk(PAx,PBx,l1,l2)
	 */
	Double Fk(const Int& k, const Int& l1, const Int& l2, 
			const Double& x1, const Double& x2); 

	/**
	 * This function is used to calculate the overlap integral 
	 * over two Gaussian primitive functions i and j
	 *
	 * \param d         : the product for the coefficients of two primitives. Already
	 *                    mulplied with pre-factor if it has
	 * \param li, mi, ni: angular momentum for i
	 * \param lj, mj, nj: angular momentum for j
	 * \param alpha     : the sum of exponents between two primitives
	 * \param isSameCenter: wether the two primitives are share the same center?
	 */
	Double overlapIntegral(
			const Double& d,   const Double& alpha, 
			const Double& PAx, const Double& PAy, const Double& PAz,  
			const Double& PBx, const Double& PBy, const Double& PBz,  
			const Int& li, const Int& mi, const Int& ni,
			const Int& lj, const Int& mj, const Int& nj, 
			bool isSameCenter);

	/**
	 * This function is used to calculate the overlap integral 
	 * over two contracted shells: (i|j)
	 * here shell data should be pure shell, no composite shell
	 * allowed
	 * \param Li   : angular momentum for shell i
	 * \param Lj   : angular momentum for shell j
	 * \param inp  : contraction degree for shell i
	 * \param icoe : coefficients for shell i
	 * \param iexp : exponetial factors for shell i
	 * \param A    : center for shell i
	 * \param jnp  : contraction degree for shell j
	 * \param jcoe : coefficients for shell j
	 * \param jexp : exponetial factors for shell j
	 * \param B    : center for shell j
	 * \return abcd: result integrals
	 */
	void directov(const Int& Li, const Int& Lj, const Int& inp, const Double* icoe, 
			const Double* iexp, const Double* A, const Int& jnp, const Double* jcoe, 
			const Double* jexp, const Double* B, Double* abcd);

	/**
	 * main function to perform two overlap integral comparison
	 * the two shell's data are given without angular momentum
	 * \param maxl : the maximum angular momentum for testing
	 * \param inp  : contraction degree for shell i
	 * \param icoe : coefficients for shell i
	 * \param iexp : exponetial factors for shell i
	 * \param A    : center for shell i
	 * \param jnp  : contraction degree for shell j
	 * \param jcoe : coefficients for shell j
	 * \param jexp : exponetial factors for shell j
	 * \param B    : center for shell j
	 */
	void ov_test(const Int& maxL,const Int& inp,const vector<Double>& icoe, 
			const vector<Double>& iexp, const Double* A, const Int& jnp, 
			const vector<Double>& jcoe, const vector<Double>& jexp, 
			const Double* B);
}

#endif
