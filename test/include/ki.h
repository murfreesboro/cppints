#ifndef KI_H
#define KI_H
#include "libgen.h"

namespace ki {

	/**
	 * This function is used to calculate the kinetic integral 
	 * over two Gaussian primitive functions i and j
	 *
	 * \param d         : the product for the coefficients of two primitives. Already
	 *                    mulplied with pre-factor if it has
	 * \param li, mi, ni: angular momentum for i
	 * \param lj, mj, nj: angular momentum for j
	 * \param ai         : the exponent for i
	 * \param aj         : the exponent for j
	 * \param isSameCenter: wether the two primitives are share the same center?
	 */
	Double kineticIntegral(
			const Double& d,   const Double& ai,  const Double& aj, 
			const Double& PAx, const Double& PAy, const Double& PAz,  
			const Double& PBx, const Double& PBy, const Double& PBz,  
			const Int& li, const Int& mi, const Int& ni,
			const Int& lj, const Int& mj, const Int& nj, 
			bool isSameCenter);

	/**
	 * This function is used to calculate the kinetic integral 
	 * over two contracted shells: (i|nabla|j)
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
	void directki(const Int& Li, const Int& Lj, const Int& inp, const Double* icoe, 
			const Double* iexp, const Double* A, const Int& jnp, const Double* jcoe, 
			const Double* jexp, const Double* B, Double* abcd);

	/**
	 * main function to perform kinetic integral comparison
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
	void ki_test(const Int& maxL,const Int& inp,const vector<Double>& icoe, 
			const vector<Double>& iexp, const Double* A, const Int& jnp, 
			const vector<Double>& jcoe, const vector<Double>& jexp, 
			const Double* B);
}

#endif
