#include "libgen.h"
#include "angmomlist.h"
#include "shellprop.h"
#include "functions.h"
#include "ov.h"
#include "norm.h"
#include "nai.h"
using namespace shellprop;
using namespace functions;
using namespace norm;
using namespace ov;
using namespace nai;

extern void hgp_os_nai(const LInt& LCode, const UInt& inp2, const UInt& nAtoms, 
		const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, 
		const Double* A, const Double* B, const Double* N, const UInt* Z, Double* abcd);

Double nai::nuclearIntegralKernel(const Double& alpha, const Double CP[3], 
		const Int& i, const Int& j, const Int& k) 
{

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
}

Double nai::nuclearIntegral(
		const Double& d,   const Double& alpha, const Double* C, const Int* Atomic,
		const Double& PAx, const Double& PAy, const Double& PAz,  
		const Double& PBx, const Double& PBy, const Double& PBz,  
		const Double& Px, const Double& Py, const Double& Pz,  
		const Int& li, const Int& mi, const Int& ni,
		const Int& lj, const Int& mj, const Int& nj, 
		const Int& nAtoms, bool isSameCenter) 
{

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
}

void nai::directnai(const Int& Li, const Int& Lj, const Int& natoms, 
		const Int& inp, const Double* icoe, const Double* iexp, const Double* A, 
		const Int& jnp, const Double* jcoe, const Double* jexp, const Double* B, 
		const Double* N, const Int* AN, Double* abcd)
{
	// number of basis sets
	Int nBra1Bas = getCartBas(Li,Li);
	Int nBra2Bas = getCartBas(Lj,Lj);

	// loop over each basis set
	Double AB2 = (A[0]-B[0])*(A[0]-B[0])+(A[1]-B[1])*(A[1]-B[1])+(A[2]-B[2])*(A[2]-B[2]);
	bool isSameCenter = false;
	if (AB2 <THRESHOLD_MATH) isSameCenter = true;
	Int count = 0;
	for(Int jBas=0; jBas<nBra2Bas; jBas++) {
		for(Int iBas=0; iBas<nBra1Bas; iBas++) {

			// get the angular momentum
			Int l1,m1,n1,l2,m2,n2;
			getlmn(Li, Li, iBas, l1, m1, n1);
			getlmn(Lj, Lj, jBas, l2, m2, n2);

			// form the primitives data
			Double nai = ZERO;
			for(Int jp=0; jp<jnp; jp++) {
				for(Int ip=0; ip<inp; ip++) {

					// prefactors etc.
					Double ia   = iexp[ip];
					Double ja   = jexp[jp];
					Double ic   = icoe[ip];
					Double jc   = jcoe[jp];
					Double alpla= ia+ja; 
					Double ab   = -ia*ja/alpla;
					Double pref = exp(ab*AB2)*ic*jc;

					// form P point according to the 
					// Gaussian pritimive product theorem
					Double adab= ia/alpla; 
					Double bdab= ja/alpla; 
					Double Px  = A[0]*adab + B[0]*bdab;
					Double Py  = A[1]*adab + B[1]*bdab;
					Double Pz  = A[2]*adab + B[2]*bdab;
					Double PAx = Px - A[0];
					Double PAy = Py - A[1];
					Double PAz = Pz - A[2];
					Double PBx = Px - B[0];
					Double PBy = Py - B[1];
					Double PBz = Pz - B[2];

					// now it's ready to call for direct calculation
					Double r = nuclearIntegral(pref, alpla, N, AN, PAx, PAy, PAz, PBx, PBy, PBz,  
							Px, Py, Pz, l1, m1, n1, l2, m2, n2, natoms, isSameCenter);
					nai += r;
				}
			}

			// push it to the abcd
			abcd[count] = nai;
			count++;
		}
	}
}

void nai::nai_test(const Int& maxL, const Int& inp, const vector<Double>& icoe, 
		const vector<Double>& iexp, const Double* A, const Int& jnp, 
		const vector<Double>& jcoe, const vector<Double>& jexp, 
		const Double* B, const vector<Double>& N, const vector<Int>& atomicNumbers)
{
	//
	// form the data for hgp calculation
	// we note that coefficient array is complicated,
	// will be depending on the angular momentum cases
	// so do it later
	//
	Double AB2 = (A[0]-B[0])*(A[0]-B[0])+(A[1]-B[1])*(A[1]-B[1])+(A[2]-B[2])*(A[2]-B[2]);
	Int np2 = inp*jnp;
	vector<Double> iexp2(np2,ZERO);
	vector<Double> fbra(np2,ZERO);
	vector<Double> P(3*np2,ZERO);
	Int count = 0;
	for(Int jp=0; jp<jnp; jp++) {
		for(Int ip=0; ip<inp; ip++) {

			// prefactors etc.
			Double ia   = iexp[ip];
			Double ja   = jexp[jp];
			Double alpla= ia+ja; 
			Double ab   = -ia*ja/alpla;
			Double pref = exp(ab*AB2)*pow(PI/alpla,1.5E0);
			iexp2[count]= ONE/alpla;
			fbra[count] = pref;

			// form P point according to the 
			// Gaussian pritimive product theorem
			Double adab = ia/alpla; 
			Double bdab = ja/alpla; 
			Double Px   = A[0]*adab + B[0]*bdab;
			Double Py   = A[1]*adab + B[1]*bdab;
			Double Pz   = A[2]*adab + B[2]*bdab;
			P[3*count+0]= Px;
			P[3*count+1]= Py;
			P[3*count+2]= Pz;
			count++;
		}
	}

	//
	// since the input atomic numbers is in int type
	// we make a copy and do it in UInt
	//
	vector<UInt> Z(atomicNumbers.size());
	for(Int i=0; i<atomicNumbers.size(); i++) {
		Z[i] = atomicNumbers[i];
	}

	// now it's real work to do ov test
	// we note that the bra1's L >= bra2's L
	//
	cout << "**************************************************************" << endl;
	cout << " two body NAI test:" << endl;
	cout << "**************************************************************" << endl;
	Int nAtoms = atomicNumbers.size();
	for(Int jpos=0; jpos<MAX_SHELL_TYPES; jpos++) {
		for(Int ipos=jpos; ipos<MAX_SHELL_TYPES; ipos++) {

			// now decode the L
			Int iLmin, iLmax;
			decodeL(SHELL_ANG_MOM_CODE[ipos],iLmin,iLmax);
			Int jLmin, jLmax;
			decodeL(SHELL_ANG_MOM_CODE[jpos],jLmin,jLmax);
			if (iLmax>maxL || jLmax>maxL) continue;

			// now form the L code
			LInt LCode = codeSQ(SHELL_ANG_MOM_CODE[ipos],SHELL_ANG_MOM_CODE[jpos]);

			// normalize the coefficient array
			Int nLBra1 = iLmax-iLmin+1;
			Int nLBra2 = jLmax-jLmin+1;
			vector<Double> bra1Coe(nLBra1*inp);
			vector<Double> bra2Coe(nLBra2*jnp);
			normCoe(iLmin, iLmax, iexp, icoe, bra1Coe);
			normCoe(jLmin, jLmax, jexp, jcoe, bra2Coe);

			// generate the shell pair for coe
			vector<Double> braCoePair(nLBra1*nLBra2*np2);
			makeC2(iLmin,iLmax,jLmin,jLmax,inp,jnp,bra1Coe,bra2Coe,braCoePair);

			// result array
			Int nBra1Bas = getCartBas(iLmin,iLmax);
			Int nBra2Bas = getCartBas(jLmin,jLmax);
			vector<Double> result(nBra1Bas*nBra2Bas);
			UInt np2_    = static_cast<UInt>(np2);
			UInt nAtoms_ = static_cast<UInt>(nAtoms);
			hgp_os_nai(LCode,np2_,nAtoms_,&braCoePair.front(),&iexp2.front(), 
					&fbra.front(),&P.front(),A,B,&N.front(),
					&Z.front(),&result.front());

			// now let's directly calculate the nai
			Int d1 = getCartBas(iLmin,iLmax);
			for(Int Lj=jLmin; Lj<=jLmax; Lj++) {
				for(Int Li=iLmin; Li<=iLmax; Li++) {

					// initilize the results
					Int nBra1Bas = getCartBas(Li,Li);
					Int nBra2Bas = getCartBas(Lj,Lj);
					vector<Double> abcd(nBra1Bas*nBra2Bas);

					// get the coefficient array
					const Double* iCArray = &bra1Coe[(Li-iLmin)*inp];
					const Double* jCArray = &bra2Coe[(Lj-jLmin)*jnp];

					// now do direct nai
					directnai(Li,Lj,nAtoms,inp,iCArray,&iexp.front(),A,jnp,jCArray, 
							&jexp.front(),B,&N.front(),&atomicNumbers.front(),&abcd.front());

					// now it's comparison
					for(Int jBas=0; jBas<nBra2Bas; jBas++) {
						for(Int iBas=0; iBas<nBra1Bas; iBas++) {
							Int iBasIndex = getBasOffset(iLmin,Li,iBas);
							Int jBasIndex = getBasOffset(jLmin,Lj,jBas);
							Int index = iBasIndex+jBasIndex*d1;
							Double v1 = hgp[index];
							Double v2 = abcd[iBas+jBas*nBra1Bas];
							if (fabs(v1-v2)>THRESH) {
								cout << "Bra1's L: " << iLmin << " " << iLmax << endl;
								cout << "Bra2's L: " << jLmin << " " << jLmax << endl;
								cout << "result did not match: " << endl;
								printf("difference  : %-16.10f\n", fabs(v1-v2));
								printf("hgp    value: %-16.10f\n", v1);
								printf("direct value: %-16.10f\n", v2);
								Int l1,m1,n1,l2,m2,n2;
								getlmn(Li, Li, iBas, l1, m1, n1);
								getlmn(Lj, Lj, jBas, l2, m2, n2);
								cout << l1 << m1 << n1 << l2 <<m2 << n2 << endl;
							}
						}
					}
				}
			}
		}
	}
}


