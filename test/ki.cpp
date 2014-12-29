#include "libgen.h"
#include "angmomlist.h"
#include "shellprop.h"
#include "functions.h"
#include "norm.h"
#include "ov.h"
#include "ki.h"
using namespace shellprop;
using namespace functions;
using namespace norm;
using namespace ov;
using namespace ki;

extern void hgp_os_kinetic(const LInt& LCode, const UInt& inp2, const Double* icoe, 
		const Double* iexp, const Double* iexpdiff, const Double* ifac, 
		const Double* P, const Double* A, const Double* B, Double* abcd);

Double ki::kineticIntegral(
		const Double& d,   const Double& ai,  const Double& aj, 
		const Double& PAx, const Double& PAy, const Double& PAz,  
		const Double& PBx, const Double& PBy, const Double& PBz,  
		const Int& li, const Int& mi, const Int& ni,
		const Int& lj, const Int& mj, const Int& nj, 
		bool isSameCenter) 
{

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

}

void ki::directki(const Int& Li, const Int& Lj, const Int& inp, const Double* icoe, 
		const Double* iexp, const Double* A, const Int& jnp, const Double* jcoe, 
		const Double* jexp, const Double* B, Double* abcd)
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
			Double ki = ZERO;
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
					Double r = kineticIntegral(pref,ia,ja, PAx, PAy, PAz, PBx, PBy, PBz,  
							l1, m1, n1, l2, m2, n2, isSameCenter);
					ki += r;
				}
			}

			// push it to the abcd
			abcd[count] = ki;
			count++;
		}
	}
}

void ki::ki_test(const Int& maxL, const Int& inp, const vector<Double>& icoe, 
		const vector<Double>& iexp, const Double* A, const Int& jnp, 
		const vector<Double>& jcoe, const vector<Double>& jexp, 
		const Double* B)
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
	vector<Double> iexpdiff(np2,ZERO);
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
			iexpdiff[count]= ia-ja;

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

	// now it's real work to do ov test
	// we note that the bra1's L >= bra2's L
	//
	cout << "**************************************************************" << endl;
	cout << " two body kinetic integrals test:" << endl;
	cout << "**************************************************************" << endl;
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
			hgp_os_kinetic(LCode,np2_,&braCoePair.front(),&iexp2.front(),&iexpdiff.front(), 
					&fbra.front(),&P.front(),A,B,&result.front());

			// now let's directly calculate the ov
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

					// now do direct ov
					directki(Li,Lj,inp,iCArray,&iexp.front(),A,jnp,jCArray, 
							&jexp.front(),B,&abcd.front());

					// now it's comparison
					for(Int jBas=0; jBas<nBra2Bas; jBas++) {
						for(Int iBas=0; iBas<nBra1Bas; iBas++) {
							Int iBasIndex = getBasOffset(iLmin,Li,iBas);
							Int jBasIndex = getBasOffset(jLmin,Lj,jBas);
							Int index = iBasIndex+jBasIndex*d1;
							Double v1 = result[index];
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


