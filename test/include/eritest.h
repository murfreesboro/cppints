#ifndef ERITEST_H
#define ERITEST_H
#include "libgen.h"

//
// Here we note that the real calculation here
// is depending on the eri.cpp and eri.h provided
// by libint package
// the direct calculation does not have any 
// coefficients incolved, therefore all 
// coefficients will be set to one in HGP part
//

namespace eritest {

	//
	// define the work type for ERI
	// basically, 2 body and 3 body eri
	// and normal 4 body eri
	//
	const Int NORMAL_ERI     = 1;
	const Int TWO_BODY_ERI   = 2;
	const Int THREE_BODY_ERI = 3;

	/**
	 * This function is used to calculate the eri  
	 * over four contracted shells: (ij|kl)
	 * here shell data should be pure shell, no composite shell
	 * allowed
	 * \param Li   : angular momentum for shell i
	 * \param Lj   : angular momentum for shell j
	 * \param Lk   : angular momentum for shell k
	 * \param Ll   : angular momentum for shell l
	 * \param inp  : contraction degree for shell i
	 * \param icoe : coefficient array  for shell i
	 * \param iexp : exponetial factors for shell i
	 * \param A    : center for shell i
	 * \param jnp  : contraction degree for shell j
	 * \param jcoe : coefficient array  for shell j
	 * \param jexp : exponetial factors for shell j
	 * \param B    : center for shell j
	 * \param knp  : contraction degree for shell k
	 * \param kcoe : coefficient array  for shell k
	 * \param kexp : exponetial factors for shell k
	 * \param C    : center for shell k
	 * \param lnp  : contraction degree for shell l
	 * \param lcoe : coefficient array  for shell l
	 * \param lexp : exponetial factors for shell l
	 * \param D    : center for shell l
	 * \return abcd: result integrals
	 */
	void directeri(const Int& Li, const Int& Lj, const Int& Lk, const Int& Ll,
			const Int& inp, const Double* icoe, const Double* iexp, const Double* A, 
			const Int& jnp, const Double* jcoe, const Double* jexp, const Double* B, 
			const Int& knp, const Double* kcoe, const Double* kexp, const Double* C, 
			const Int& lnp, const Double* lcoe, const Double* lexp, const Double* D, 
			Double* abcd);

	/**
	 * main function to perform ERI comparison
	 * this is without coefficient information
	 * \param maxl     : the maximum angular momentum for testing
	 * \param maxl     : the maximum angular momentum used for aux 
	 *                   shell dimension for testing
	 * \param workType : normal ERI, 2/3 body ERI?
	 * \param inp      : contraction degree for shell i
	 * \param icoe     : coefficient array  for shell i
	 * \param iexp     : exponetial factors for shell i
	 * \param A        : center for shell i
	 * \param jnp      : contraction degree for shell j
	 * \param jcoe     : coefficient array  for shell j
	 * \param jexp     : exponetial factors for shell j
	 * \param B        : center for shell j
	 * \param knp      : contraction degree for shell k
	 * \param kcoe     : coefficient array  for shell k
	 * \param kexp     : exponetial factors for shell k
	 * \param C        : center for shell k
	 * \param lnp      : contraction degree for shell l
	 * \param lcoe     : coefficient array  for shell l
	 * \param lexp     : exponetial factors for shell l
	 * \param D        : center for shell l
	 */
	void eri_test(const Int& maxL, const Int& auxMaxL, const Int& workType,
			const Int& inp,const vector<Double>& icoe,const vector<Double>& iexp, const Double* A, 
			const Int& jnp,const vector<Double>& jcoe,const vector<Double>& jexp, const Double* B,
			const Int& knp,const vector<Double>& kcoe,const vector<Double>& kexp, const Double* C, 
			const Int& lnp,const vector<Double>& lcoe,const vector<Double>& lexp, const Double* D);
}

#endif
