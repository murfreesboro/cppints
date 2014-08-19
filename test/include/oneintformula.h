/**
 * \file  oneintformula.h
 * \brief This file contains the math formulas for single electron integral
 *        finally used only for debugging purpose
 * \author Fenglai Liu and Jing Kong
 * \date   June, 2011
 */
#ifndef ONEINTFORMULA_H
#define ONEINTFORMULA_H

#include "functions.h"
using namespace functions;

namespace oneint {

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
	inline Double Fk(const Int& k, const Int& l1, const Int& l2, 
			const Double& x1, const Double& x2) {
		if (l1 > 0 && l2 > 0) {
			Int low   = max(0,k-l2);
			Int upper = min(k,l1);
			Double sum = ZERO;
			for (Int i=low;i<=upper;i++) {
				sum += pow(x1,l1-i)*pow(x2,l2-k+i)*bionomialCoe(l1,i)*bionomialCoe(l2,k-i);
			}
			return sum;
		} else if (l1 >0 && l2 == 0) {
			return pow(x1,l1-k)*bionomialCoe(l1,k);
		} else if (l1 ==0 && l2 > 0) {
			return pow(x2,l2-k)*bionomialCoe(l2,k);
		} else if (l1 ==0 && l2 == 0) {
			return ONE;
		} else {
			crash(true, "In the Fk function the value for l1 and l2 are invalid");
			return ZERO;
		}
	};


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
	inline Double overlapIntegral(
			const Double& d,   const Double& alpha, 
			const Double& PAx, const Double& PAy, const Double& PAz,  
			const Double& PBx, const Double& PBy, const Double& PBz,  
			const Int& li, const Int& mi, const Int& ni,
			const Int& lj, const Int& mj, const Int& nj, 
			bool isSameCenter) {

		if (isSameCenter) {
			return d*gammaFunInt(li+lj,alpha)*gammaFunInt(mi+mj,alpha)*gammaFunInt(ni+nj,alpha);
		}else{
			Double sum = ZERO;
			for (Int i=0; i<=(Int)(li+lj)/2; i++) {
				for (Int j=0; j<=(Int)(mi+mj)/2; j++) {
					for (Int k=0; k<=(Int)(ni+nj)/2; k++) {
						Double fk = Fk(2*i,li,lj,PAx,PBx)*Fk(2*j,mi,mj,PAy,PBy)*Fk(2*k,ni,nj,PAz,PBz);
						Double g  = gammaFunInt(2*i,alpha)*gammaFunInt(2*j,alpha)*gammaFunInt(2*k,alpha);
						sum += fk*g; 
					}
				}
			}
			return sum*d;
		}
	};

	/**
	 * this function is used to calculate each integral piece of three center overlap 
	 * integral(on x, y or z)
	 * see the OS 1986 paper equation 17
	 */
	inline Double threeOverlapIxyz(const Double& alpha, const Double& beta, const Double& gamma,
			const Int& na, const Int& nb, const Int& nc, const Double& G, const Double& A,
			const Double& B, const Double& C) {
		Double zeta   = alpha + beta + gamma;
		Double pidz   = pow(PI/zeta,0.5E0);
		Double oned2z = 0.5E0/zeta;
		//cout << "oned2z " << oned2z << endl;
		Double result = 0;
		for(Int ka=0; ka<=na; ka++) {
			for(Int kb=0; kb<=nb; kb++) {
				for(Int kc=0; kc<=nc; kc++) {
					Int k=ka+kb+kc;
					if (k%2 != 0) continue;
					Double GA = pow((G-A),na-ka);
					Double GB = pow((G-B),nb-kb);
					Double GC = pow((G-C),nc-kc);
					Double in = doubleFac(k-1)*pow(oned2z,k);
					cout << "k index " << ka << kb << kc << endl;
					//cout << "bionomialCoe  " << bionomialCoe(na,ka) << " "
					//	<< bionomialCoe(nb,kb) << " " << bionomialCoe(nc,kc) << endl;
					//cout << "GA " << GA << endl;
					//cout << "GB " << GB << endl;
					//cout << "GC " << GC << endl;
					//cout << "k " << k << endl;
					//cout << "integral " << in << endl;
					Double r = bionomialCoe(na,ka)*bionomialCoe(nb,kb)*bionomialCoe(nc,kc)*GA*GB*GC*in;
					//cout << "r " << r << endl;
					result += r;
				}
			}
		}
		return result*pidz;
	};

	/**
	 * this function is used to calculate three center integrals: (ij|k)
	 */
	inline Double threeOverlap(const Double& alpha, const Double& beta, const Double& gamma,
			const Double* A, const Double* B, const Double* C,
			const Int& li, const Int& mi, const Int& ni,
			const Int& lj, const Int& mj, const Int& nj, 
			const Int& lk, const Int& mk, const Int& nk) {
		Double zab = alpha + beta;
		Double zeta= zab + gamma;
		Double Px  = (alpha*A[0] + beta*B[0])/zab;
		Double Py  = (alpha*A[1] + beta*B[1])/zab;
		Double Pz  = (alpha*A[2] + beta*B[2])/zab;
		Double Gx  = (zab*Px     + gamma*C[0])/zeta;
		Double Gy  = (zab*Py     + gamma*C[1])/zeta;
		Double Gz  = (zab*Pz     + gamma*C[2])/zeta;
		Double AB2 = (A[0]-B[0])*(A[0]-B[0]) + (A[1]-B[1])*(A[1]-B[1]) + (A[2]-B[2])*(A[2]-B[2]);
		Double PC2 = (Px  -C[0])*(Px  -C[0]) + (Py  -C[1])*(Py  -C[1]) + (Pz  -C[2])*(Pz  -C[2]);
		Double kappa = exp(-(alpha*beta/zab)*AB2)*exp(-(zab*gamma/zeta)*PC2);
		cout << "l " << li << lj << lk << endl;
		cout << "m " << mi << mj << mk << endl;
		cout << "n " << ni << nj << nk << endl;
		//cout << "exp " << alpha << beta << gamma << endl;
		cout << "Ix " << endl;
		Double Ix = threeOverlapIxyz(alpha,beta,gamma,li,lj,lk,Gx,A[0],B[0],C[0]);
		cout << "Iy " << endl;
		Double Iy = threeOverlapIxyz(alpha,beta,gamma,mi,mj,mk,Gy,A[1],B[1],C[1]);
		cout << "Iz " << endl;
		Double Iz = threeOverlapIxyz(alpha,beta,gamma,ni,nj,nk,Gz,A[2],B[2],C[2]);
		//cout << "I " << kappa << " " << Ix << " " << Iy << " " << Iz << endl;
		return kappa*Ix*Iy*Iz;
	};

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
	inline Double kineticIntegral(
			const Double& d,   const Double& ai,  const Double& aj, 
			const Double& PAx, const Double& PAy, const Double& PAz,  
			const Double& PBx, const Double& PBy, const Double& PBz,  
			const Int& li, const Int& mi, const Int& ni,
			const Int& lj, const Int& mj, const Int& nj, 
			bool isSameCenter) {

		Double r = ZERO;

		// Ix
		if (li-1>=0){
			r += MINUS_TWO*aj*li*overlapIntegral(d,ai+aj,PAx,PAy,PAz,PBx,PBy,PBz,
					li-1,mi,ni,lj+1,mj,nj,isSameCenter);
		}
		if (lj-1>=0){
			r += MINUS_TWO*ai*lj*overlapIntegral(d,ai+aj,PAx,PAy,PAz,PBx,PBy,PBz,
					li+1,mi,ni,lj-1,mj,nj,isSameCenter);
		}
		if (li-1>=0 && lj-1>=0){
			r += li*lj*overlapIntegral(d,ai+aj,PAx,PAy,PAz,PBx,PBy,PBz,
					li-1,mi,ni,lj-1,mj,nj,isSameCenter);
		}
		r += FOUR*ai*aj*overlapIntegral(d,ai+aj,PAx,PAy,PAz,PBx,PBy,PBz,
				li+1,mi,ni,lj+1,mj,nj,isSameCenter);


		// Iy
		if (mi-1>=0){
			r += MINUS_TWO*aj*mi*overlapIntegral(d,ai+aj,PAx,PAy,PAz,PBx,PBy,PBz,
					li,mi-1,ni,lj,mj+1,nj,isSameCenter);
		}
		if (mj-1>=0){
			r += MINUS_TWO*ai*mj*overlapIntegral(d,ai+aj,PAx,PAy,PAz,PBx,PBy,PBz,
					li,mi+1,ni,lj,mj-1,nj,isSameCenter);
		}
		if (mi-1>=0 && mj-1>=0){
			r += mi*mj*overlapIntegral(d,ai+aj,PAx,PAy,PAz,PBx,PBy,PBz,
					li,mi-1,ni,lj,mj-1,nj,isSameCenter);
		}
		r += FOUR*ai*aj*overlapIntegral(d,ai+aj,PAx,PAy,PAz,PBx,PBy,PBz,
				li,mi+1,ni,lj,mj+1,nj,isSameCenter);


		// Iz
		if (ni-1>=0){
			r += MINUS_TWO*aj*ni*overlapIntegral(d,ai+aj,PAx,PAy,PAz,PBx,PBy,PBz,
					li,mi,ni-1,lj,mj,nj+1,isSameCenter);
		}
		if (nj-1>=0){
			r += MINUS_TWO*ai*nj*overlapIntegral(d,ai+aj,PAx,PAy,PAz,PBx,PBy,PBz,
					li,mi,ni+1,lj,mj,nj-1,isSameCenter);
		}
		if (ni-1>=0 && nj-1>=0){
			r += ni*nj*overlapIntegral(d,ai+aj,PAx,PAy,PAz,PBx,PBy,PBz,
					li,mi,ni-1,lj,mj,nj-1,isSameCenter);
		}
		r += FOUR*ai*aj*overlapIntegral(d,ai+aj,PAx,PAy,PAz,PBx,PBy,PBz,
				li,mi,ni+1,lj,mj,nj+1,isSameCenter);

		r *= ONE/TWO;
		return r;

	};


	/**
	 * This function is used to calculate the integral part of NAI 
	 * over two Gaussian primitive functions i and j
	 * \param P         : the shell pair center
	 * \param C         : the center for the nuclear
	 * \param i,j,k     : the exponent order
	 * \param alpha     : sum of exponents
	 */
	inline Double nuclearIntegralKernel(const Double& alpha, const Double CP[3], 
			const Int& i, const Int& j, const Int& k) {

		Double r   = ZERO;
		Double CPx = CP[0];
		Double CPy = CP[1];
		Double CPz = CP[2];
		Double PC  = sqrt(CPx*CPx+CPy*CPy+CPz*CPz);
		Int    L   = i+j+k;

		bool isSameCenter = (fabs(PC-ZERO) < THRESHOLD_MATH) ? true : false ;
		//cout << " is same center?? " << isSameCenter << endl;
		if(isSameCenter) {
			if (i%2 == 0 && j%2 == 0 && k%2 == 0) {
				Int Lp    = (i+j+k)/2;
				Double c1 = TWO*PI/(pow(alpha,Lp+1)*pow(TWO,Lp));
				Double c2 = doubleFac(i-1)*doubleFac(j-1)*doubleFac(k-1);
				Double in = oneMinusT2L(Lp); 
				r = c1*c2*in;
			}
		}else{
			for (Int ip=0; ip<=(Int)(i/2); ip++) {
				for (Int jp=0; jp<=(Int)(j/2); jp++) {
					for (Int kp=0; kp<=(Int)(k/2); kp++) {

						// prefactors in the integral
						Int Lp    = ip+jp+kp;
						Double c1 = TWO*PI/(pow(alpha,Lp+1)*pow(TWO,Lp));
						Double c2 = bionomialCoe(i,2*ip)*bionomialCoe(j,2*jp)*bionomialCoe(k,2*kp);
						Double c3 = pow(CPx,i-2*ip)*pow(CPy,j-2*jp)*pow(CPz,k-2*kp);
						Double c4 = doubleFac(2*ip-1)*doubleFac(2*jp-1)*doubleFac(2*kp-1);
						Double c5 = ONE/pow(PC,2*(L-Lp)+1);
						Double K  = c1*c2*c3*c4*c5;

						// the integral over s
						Double in = fmNAI(PC,alpha,Lp,L);
						r += K*in;
					}
				}
			}
		}
		//printf("nuclear integral core: %f\n", r);
		return r;
	};



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
	inline Double nuclearIntegral(
			const Double& d,   const Double& alpha, const Double* C, const Int* Atomic,
			const Double& PAx, const Double& PAy, const Double& PAz,  
			const Double& PBx, const Double& PBy, const Double& PBz,  
			const Double& Px, const Double& Py, const Double& Pz,  
			const Int& li, const Int& mi, const Int& ni,
			const Int& lj, const Int& mj, const Int& nj, 
			const Int& nAtoms, bool isSameCenter) {

		Double result = ZERO;

		// loop over all of nuclears
		for(Int n=0; n<nAtoms; n++) {
			Double CP[3];
			CP[0] = C[3*n  ] - Px;
			CP[1] = C[3*n+1] - Py;
			CP[2] = C[3*n+2] - Pz;
			Double Z = MINUS_ONE*Atomic[n];

			// real calculation
			if (isSameCenter) {
				result += Z*d*nuclearIntegralKernel(alpha,CP,li+lj,mi+mj,ni+nj); 
			}else{
				Double sum = ZERO;
				for (Int i=0; i<=li+lj; i++) {
					for (Int j=0; j<=mi+mj; j++) {
						for (Int k=0; k<=ni+nj; k++) {
							Double fk = Fk(i,li,lj,PAx,PBx)*Fk(j,mi,mj,PAy,PBy)*Fk(k,ni,nj,PAz,PBz);
							Double g  = nuclearIntegralKernel(alpha,CP,i,j,k); 
							sum += fk*g; 
						}
					}
				}
				result += Z*sum*d;
			}
		}
		return result;
	};


}

#endif
