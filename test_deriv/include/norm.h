#ifndef NORM_H
#define NORM_H
#include "libgen.h"

namespace norm {

	/**
	 * This function is used to generate the normalization factor
	 * for angular momentum of (L,0,0)
	 * \param iShell the L
	 * \param nPrim  length of c and e
	 * \param c      coefficients of shell data
	 * \param e      exponetial factors of shell data
	 * \return       normalization factor for basis set with (L,0,0)
	 */
	Double normCartBasisSets(const Int& iShell, const Int& nPrim, 
			const Double* c, const Double* e);

	/**
	 * this function is used to get the ratio between the normalization
	 * factor of N_{L00}(lx=L,ly=0,lz=0) and N_{lmn}(lx=l,ly=m,lz=n)
	 * basically, here what we compute is the ratio of N_{lmn}/N_{L00}
	 */
	void scaleNormFactors(const Int& Lmin, const Int& Lmax, Double* N); 

	/**
	 * normalize the coefficients array for the given angular momentum 
	 * of lmin and lmax
	 * the input c array should contain all of original
	 * \param e exponential factor array
	 * \return c the normalized coefficients array
	 */
	void normCoe(const Int& lmin, const Int& lmax, const vector<Double>& e, 
			const vector<Double>& oric, vector<Double>& c);

	/**
	 * generate the coe2 array for the HGP part of codes
	 * \param l1min, l1max: angular momentum for shell 1
	 * \param l2min, l2max: angular momentum for shell 2
	 * \param inp,jnp     : number of primitive functions for the shell 1/shell 2
	 * \param c1          : coefficient array for shell 1
	 * \param c2          : coefficient array for shell 2
	 * \return coe2       : result coefficient pair array
	 */
	void makeC2(const Int& l1min, const Int& l1max, 
		const Int& l2min, const Int& l2max, const Int& inp, const Int& jnp,
		const vector<Double>& c1, const vector<Double>& c2, vector<Double>& coe2);

	/**
	 * for the given result integrals (1 body array), scale it with with normalized
	 * ratio factors
	 */
	void scale1BodyInts(const Int& LBra1Min, const Int& LBra1Max, Double* result); 

	/**
	 * for the given result integrals (2 body array), scale it with with normalized
	 * ratio factors
	 */
	void scale2BodyInts(const Int& LBra1Min, const Int& LBra1Max, 
		const Int& LBra2Min, const Int& LBra2Max, Double* result); 

	/**
	 * for the given result integrals (3 body array), scale it with with normalized
	 * ratio factors
	 */
	void scale3BodyInts(const Int& LBra1Min, const Int& LBra1Max, 
			const Int& LBra2Min, const Int& LBra2Max, 
			const Int& LKet1Min, const Int& LKet1Max, Double* result);

	/**
	 * for the given result integrals (4 body array), scale it with with normalized
	 * ratio factors
	 */
	void scale4BodyInts(const Int& LBra1Min, const Int& LBra1Max, 
			const Int& LBra2Min, const Int& LBra2Max, 
			const Int& LKet1Min, const Int& LKet1Max, 
			const Int& LKet2Min, const Int& LKet2Max, Double* result);

}

#endif
