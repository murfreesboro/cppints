#ifndef EXPR12_H
#define EXPR12_H
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

namespace expr12test {

	/**
	 * main function to perform the integral comparison
	 * we treat the last shell as dummy, and as we set omega
	 * value to be huge; then the integral becomes three body overlap
	 *
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
	void expr12_test(const Int& maxL, 
		const Int& inp, const vector<Double>& icoe, const vector<Double>& iexp, const Double* A, 
		const Int& jnp, const vector<Double>& jcoe, const vector<Double>& jexp, const Double* B, 
		const Int& knp, const vector<Double>& kcoe, const vector<Double>& kexp, const Double* C);

}


#endif
