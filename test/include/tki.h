#ifndef TKI_H
#define TKI_H
#include "libgen.h"

namespace tki {

	/**
	 * perform three body kinetic integral calculations over the gaussian
	 * primitive functions
	 */
	Double threeki(const Double& alpha, const Double& beta, const Double& gamma,
			const Double* A, const Double* B, const Double* C,
			const Int& li, const Int& mi, const Int& ni,
			const Int& lj, const Int& mj, const Int& nj, 
			const Int& lk, const Int& mk, const Int& nk);

	/**
	 * This function is used to calculate the kinetic integral 
	 * over three contracted shells: (ij|k)
	 * here shell data should be pure shell, no composite shell
	 * allowed
	 * \param Li   : angular momentum for shell i
	 * \param Lj   : angular momentum for shell j
	 * \param Lk   : angular momentum for shell k
	 * \param inp  : contraction degree for shell i
	 * \param icoe : coeffcient array for shell i
	 * \param iexp : exponetial factors for shell i
	 * \param A    : center for shell i
	 * \param jnp  : contraction degree for shell j
	 * \param jcoe : coeffcient array for shell j
	 * \param jexp : exponetial factors for shell j
	 * \param B    : center for shell k
	 * \param knp  : contraction degree for shell k
	 * \param kcoe : coeffcient array for shell k
	 * \param kexp : exponetial factors for shell k
	 * \param C    : center for shell k
	 * \return abcd: result integrals
	 */
	void directtki(const Int& Li, const Int& Lj, const Int& Lk, 
			const Int& inp, const Double* icoe, const Double* iexp, const Double* A, 
			const Int& jnp, const Double* jcoe, const Double* jexp, const Double* B, 
			const Int& knp, const Double* kcoe, const Double* kexp, const Double* C, 
			Double* abcd);

	/**
	 * main function to perform three kinetic integral comparison
	 * the shell's data are given without angular momentum
	 * \param maxl : the maximum angular momentum for testing
	 * \param inp  : contraction degree for shell i
	 * \param icoe : coefficient array for shell  i
	 * \param iexp : exponetial factors for shell i
	 * \param A    : center for shell i
	 * \param jnp  : contraction degree for shell j
	 * \param jcoe : coefficient array for shell  j
	 * \param jexp : exponetial factors for shell j
	 * \param B    : center for shell j
	 * \param knp  : contraction degree for shell k
	 * \param kcoe : coefficient array for shell  k
	 * \param kexp : exponetial factors for shell k
	 * \param C    : center for shell k
	 */
	void tki_test(const Int& maxL,
			const Int& inp,const vector<Double>& icoe,const vector<Double>& iexp,const Double* A,
			const Int& jnp,const vector<Double>& jcoe,const vector<Double>& jexp,const Double* B,
			const Int& knp,const vector<Double>& kcoe,const vector<Double>& kexp,const Double* C);
}

#endif
