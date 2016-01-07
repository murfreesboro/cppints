#ifndef NAI_H
#define NAI_H
#include "libgen.h"

namespace nai {

	/**
	 * This function is used to calculate the integral part of NAI 
	 * over two Gaussian primitive functions i and j
	 * \param P         : the shell pair center
	 * \param C         : the center for the nuclear
	 * \param i,j,k     : the exponent order
	 * \param alpha     : sum of exponents
	 */
	Double nuclearIntegralKernel(const Double& alpha, const Double CP[3], 
			const Int& i, const Int& j, const Int& k);

	/**
	 * This function is used to calculate the nuclear attraction integral 
	 * over two Gaussian primitive functions i and j
	 *
	 * \param d         : the product for the coefficients of two primitives. Already
	 *                    mulplied with pre-factor if it has
	 * \param C         : coordinates for the nuclear
	 * \param P         : coordinates for the shell center
	 * \param li, mi, ni: angular momentum for i
	 * \param lj, mj, nj: angular momentum for j
	 * \param alpha     : the sum of exponents between two primitives
	 * \param isSameCenter: wether the two primitives are share the same center?
	 */
	Double nuclearIntegral(
			const Double& d,   const Double& alpha, const Double* C, const Int* Atomic,
			const Double& PAx, const Double& PAy, const Double& PAz,  
			const Double& PBx, const Double& PBy, const Double& PBz,  
			const Double& Px, const Double& Py, const Double& Pz,  
			const Int& li, const Int& mi, const Int& ni,
			const Int& lj, const Int& mj, const Int& nj, 
			const Int& nAtoms, bool isSameCenter);

	/**
	 * This function is used to calculate the NAI derivatives
	 * over two contracted shells: (i|j)
	 * here shell data should be pure shell, no composite shell
	 * allowed
	 * \param Li    : angular momentum for shell i
	 * \param Lj    : angular momentum for shell j
	 * \param natoms: number of atoms for the passing N and AN
	 * \param inp   : contraction degree for shell i
	 * \param icoe  : coefficients for shell i
	 * \param iexp  : exponetial factors for shell i
	 * \param A     : center for shell i
	 * \param jnp   : contraction degree for shell j
	 * \param jcoe  : coefficients for shell j
	 * \param jexp  : exponetial factors for shell j
	 * \param B     : center for shell j
	 * \param N     : nuclear centers
	 * \param AN    : atomic numbers array
	 * \param pos   : derivatives position
	 * \param dir   : derivatives direction
	 * \return abcd : result integrals
	 */
	void directnai_d1(const Int& Li, const Int& Lj, const Int& natoms, 
			const Int& inp, const Double* icoe, const Double* iexp, const Double* A, 
			const Int& jnp, const Double* jcoe, const Double* jexp, const Double* B, 
			const Double* N, const Int* AN, const Int& pos, const Int& dir, Double* abcd);

	/**
	 * This function is used to calculate the NAI second derivatives
	 * over two contracted shells: (i|j) 
	 * here shell data should be pure shell, no composite shell
	 * allowed
	 * \param Li    : angular momentum for shell i
	 * \param Lj    : angular momentum for shell j
	 * \param natoms: number of atoms for the passing N and AN
	 * \param inp   : contraction degree for shell i
	 * \param icoe  : coefficients for shell i
	 * \param iexp  : exponetial factors for shell i
	 * \param A     : center for shell i
	 * \param jnp   : contraction degree for shell j
	 * \param jcoe  : coefficients for shell j
	 * \param jexp  : exponetial factors for shell j
	 * \param B     : center for shell j
	 * \param N     : nuclear centers
	 * \param AN    : atomic numbers array
	 * \param pos   : derivatives position
	 * \param dir   : derivatives direction
	 * \return abcd : result integrals
	 */
	void directnai_d2(const Int& Li, const Int& Lj, const Int& natoms, 
			const Int& inp, const Double* icoe, const Double* iexp, const Double* A, 
			const Int& jnp, const Double* jcoe, const Double* jexp, const Double* B, 
			const Double* N, const Int* AN, const Int& pos, const Int& dir, Double* abcd);

	/**
	 * This function is used to calculate the NAI second derivatives
	 * over two contracted shells: (i|j)
	 * here shell data should be pure shell, no composite shell
	 * allowed
	 * \param Li    : angular momentum for shell i
	 * \param Lj    : angular momentum for shell j
	 * \param natoms: number of atoms for the passing N and AN
	 * \param inp   : contraction degree for shell i
	 * \param icoe  : coefficients for shell i
	 * \param iexp  : exponetial factors for shell i
	 * \param A     : center for shell i
	 * \param jnp   : contraction degree for shell j
	 * \param jcoe  : coefficients for shell j
	 * \param jexp  : exponetial factors for shell j
	 * \param B     : center for shell j
	 * \param N     : nuclear centers
	 * \param AN    : atomic numbers array
	 * \param pos   : derivatives position
	 * \param dir   : derivatives direction
	 * \return abcd : result integrals
	 */
	void directnai_d2(const Int& Li, const Int& Lj, const Int& natoms, 
			const Int& inp, const Double* icoe, const Double* iexp, const Double* A, 
			const Int& jnp, const Double* jcoe, const Double* jexp, const Double* B, 
			const Double* N, const Int* AN, const Int& pos1, const Int& pos2, 
			const Int& dir1, const Int& dir2, Double* abcd);

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
	 * \param N    : nuclear centers coordinates
	 * \param atomicNumbers: atomic numbers array
	 */
	void nai_test(const string& derivFile, const Int& maxL, const Int& order, 
		const Int& inp,const vector<Double>& icoe, const vector<Double>& iexp, const Double* A, 
		const Int& jnp,const vector<Double>& jcoe, const vector<Double>& jexp, const Double* B, 
		const vector<Double>& N, const vector<Int>& atomicNumbers);

}

#endif
