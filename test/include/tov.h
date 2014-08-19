#ifndef TOV_H
#define TOV_H
#include "libgen.h"

//
// for the tov direct calculation, we note that the use of 
// threeOverlapIxyz and threeOverlap can not give the correct
// results. Do not know why, and we tried a lot in debugging 
// the codes. Therefore, we finally called the threeov by
// transform the three ov into two ov. Now everyhing is happy now.
//
// Therefore, do not call the function of threeOverlap.
// present here only for archive purpose.
//

namespace tov {

	/**
	 * this function is used to calculate each integral piece of three center overlap 
	 * integral(on x, y or z)
	 * see the paper below, equation 17
	 * Efficient recursive computation of molecular integrals 
	 * over Cartesian Gaussian functions 
	 * Obara, S. and Saika, A.
	 * J. CHEM. Phys. 84, 3963, 1986
	 */
	Double threeOverlapIxyz(const Double& alpha, const Double& beta, const Double& gamma,
			const Int& na, const Int& nb, const Int& nc, const Double& G, const Double& A,
			const Double& B, const Double& C);

	/**
	 * this function is used to calculate three center integrals: (ij|k)
	 * implementation is following the same paper above
	 */
	Double threeOverlap(const Double& alpha, const Double& beta, const Double& gamma,
			const Double* A, const Double* B, const Double* C,
			const Int& li, const Int& mi, const Int& ni,
			const Int& lj, const Int& mj, const Int& nj, 
			const Int& lk, const Int& mk, const Int& nk);

	/**
	 * this is another way to perform three body overlap integral calculation
	 * basically, it transform the three body integral into the two body
	 * integral
	 */
	Double threeov(const Double& alpha, const Double& beta, const Double& gamma,
			const Double* A, const Double* B, const Double* C,
			const Int& li, const Int& mi, const Int& ni,
			const Int& lj, const Int& mj, const Int& nj, 
			const Int& lk, const Int& mk, const Int& nk);

	/**
	 * This function is used to calculate the overlap integral 
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
	void directtov(const Int& Li, const Int& Lj, const Int& Lk, 
			const Int& inp, const Double* icoe, const Double* iexp, const Double* A, 
			const Int& jnp, const Double* jcoe, const Double* jexp, const Double* B, 
			const Int& knp, const Double* kcoe, const Double* kexp, const Double* C, 
			Double* abcd);

	/**
	 * main function to perform three overlap integral comparison
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
	void tov_test(const Int& maxL,
			const Int& inp,const vector<Double>& icoe,const vector<Double>& iexp,const Double* A,
			const Int& jnp,const vector<Double>& jcoe,const vector<Double>& jexp,const Double* B,
			const Int& knp,const vector<Double>& kcoe,const vector<Double>& kexp,const Double* C);
}

#endif
