#ifndef MOM_H
#define MOM_H
#include "libgen.h"

//
// MOM integrals could be derived directly from 
// three body integrals. simply, the gamma is 
// set to zero, we can get the MOM integrals
//

namespace mom {

	/**
	 * This function is used to calculate the moment integral 
	 * over two contracted shells: (i|\nu|j)
	 * here shell data should be pure shell, no composite shell
	 * allowed
	 * \param Li   : angular momentum for shell i
	 * \param Lj   : angular momentum for shell j
	 * \param Lk   : angular momentum for the momentum operator
	 * \param inp  : contraction degree for shell i
	 * \param icoe : coefficients for shell i
	 * \param iexp : exponetial factors for shell i
	 * \param A    : center for shell i
	 * \param jnp  : contraction degree for shell j
	 * \param jcoe : coefficients for shell j
	 * \param jexp : exponetial factors for shell j
	 * \param B    : center for shell j
	 * \param C    : the center of operator
	 * \return abcd: result integrals
	 */
	void directmom(const Int& Li, const Int& Lj, const Int& Lk, 
			const Int& inp, const Double* icoe, const Double* iexp, const Double* A, 
			const Int& jnp, const Double* jcoe, const Double* jexp, const Double* B, 
			const Double* C, Double* abcd);

	/**
	 * main function to perform moment integral comparison
	 * \param maxl : the maximum angular momentum for testing
	 * \param inp  : contraction degree for shell i
	 * \param icoe : coefficients for shell i
	 * \param iexp : exponetial factors for shell i
	 * \param A    : center for shell i
	 * \param jnp  : contraction degree for shell j
	 * \param jcoe : coefficients for shell j
	 * \param jexp : exponetial factors for shell j
	 * \param B    : center for shell j
	 * \param C    : center moment operator
	 */
	void mom_test(const Int& maxL, const Int& auxL,
			const Int& inp, const vector<Double>& icoe, const vector<Double>& iexp, 
			const Double* A, const Int& jnp, const vector<Double>& jcoe, 
			const vector<Double>& jexp, const Double* B, const Double* C);
}

#endif
