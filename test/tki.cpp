#include "libgen.h"
#include "angmomlist.h"
#include "shellprop.h"
#include "functions.h"
#include "norm.h"
#include "tov.h"
#include "tki.h"
using namespace shellprop;
using namespace functions;
using namespace norm;
using namespace tov;
using namespace tki;

extern void hgp_os_threebodyki(const LInt& LCode, const UInt& inp2, const UInt& jnp2, 
		const Double* icoe, const Double* iexp, const Double* iexpdiff, 
		const Double* ifac,  const Double* P, const Double* A, const Double* B, 
		const Double* jcoe, const Double* jexp, const Double* C, Double* abcd);

Double tki::threeki(const Double& alpha, const Double& beta, const Double& gamma,
		const Double* A, const Double* B, const Double* C,
		const Int& li, const Int& mi, const Int& ni,
		const Int& lj, const Int& mj, const Int& nj, 
		const Int& lk, const Int& mk, const Int& nk) 
{

	Double r = ZERO;

	// Ix
	
	// A and C
	if (li-1>=0 && lk-1>=0){
		r += li*lk*threeov(alpha,beta,gamma,A,B,C,li-1,mi,ni,lj,mj,nj,lk-1,mk,nk); 
	}
	if (lk-1>=0){
		r += MINUS_TWO*alpha*lk*threeov(alpha,beta,gamma,A,B,C,li+1,mi,ni,lj,mj,nj,lk-1,mk,nk); 
	}
	if (li-1>=0){
		r += MINUS_TWO*gamma*li*threeov(alpha,beta,gamma,A,B,C,li-1,mi,ni,lj,mj,nj,lk+1,mk,nk); 
	}
	r += FOUR*alpha*gamma*threeov(alpha,beta,gamma,A,B,C,li+1,mi,ni,lj,mj,nj,lk+1,mk,nk);

	// B and C
	if (lj-1>=0 && lk-1>=0){
		r += lj*lk*threeov(alpha,beta,gamma,A,B,C,li,mi,ni,lj-1,mj,nj,lk-1,mk,nk); 
	}
	if (lk-1>=0){
		r += MINUS_TWO*beta*lk*threeov(alpha,beta,gamma,A,B,C,li,mi,ni,lj+1,mj,nj,lk-1,mk,nk); 
	}
	if (lj-1>=0){
		r += MINUS_TWO*gamma*lj*threeov(alpha,beta,gamma,A,B,C,li,mi,ni,lj-1,mj,nj,lk+1,mk,nk); 
	}
	r += FOUR*beta*gamma*threeov(alpha,beta,gamma,A,B,C,li,mi,ni,lj+1,mj,nj,lk+1,mk,nk);


	// Iy
	
	// A and C
	if (mi-1>=0 && mk-1>=0){
		r += mi*mk*threeov(alpha,beta,gamma,A,B,C,li,mi-1,ni,lj,mj,nj,lk,mk-1,nk); 
	}
	if (mk-1>=0){
		r += MINUS_TWO*alpha*mk*threeov(alpha,beta,gamma,A,B,C,li,mi+1,ni,lj,mj,nj,lk,mk-1,nk); 
	}
	if (mi-1>=0){
		r += MINUS_TWO*gamma*mi*threeov(alpha,beta,gamma,A,B,C,li,mi-1,ni,lj,mj,nj,lk,mk+1,nk); 
	}
	r += FOUR*alpha*gamma*threeov(alpha,beta,gamma,A,B,C,li,mi+1,ni,lj,mj,nj,lk,mk+1,nk);

	// B and C
	if (mj-1>=0 && mk-1>=0){
		r += mj*mk*threeov(alpha,beta,gamma,A,B,C,li,mi,ni,lj,mj-1,nj,lk,mk-1,nk); 
	}
	if (mk-1>=0){
		r += MINUS_TWO*beta*mk*threeov(alpha,beta,gamma,A,B,C,li,mi,ni,lj,mj+1,nj,lk,mk-1,nk); 
	}
	if (mj-1>=0){
		r += MINUS_TWO*gamma*mj*threeov(alpha,beta,gamma,A,B,C,li,mi,ni,lj,mj-1,nj,lk,mk+1,nk); 
	}
	r += FOUR*beta*gamma*threeov(alpha,beta,gamma,A,B,C,li,mi,ni,lj,mj+1,nj,lk,mk+1,nk);


	// Iz
	
	// A and C
	if (ni-1>=0 && nk-1>=0){
		r += ni*nk*threeov(alpha,beta,gamma,A,B,C,li,mi,ni-1,lj,mj,nj,lk,mk,nk-1); 
	}
	if (nk-1>=0){
		r += MINUS_TWO*alpha*nk*threeov(alpha,beta,gamma,A,B,C,li,mi,ni+1,lj,mj,nj,lk,mk,nk-1); 
	}
	if (ni-1>=0){
		r += MINUS_TWO*gamma*ni*threeov(alpha,beta,gamma,A,B,C,li,mi,ni-1,lj,mj,nj,lk,mk,nk+1); 
	}
	r += FOUR*alpha*gamma*threeov(alpha,beta,gamma,A,B,C,li,mi,ni+1,lj,mj,nj,lk,mk,nk+1);

	// B and C
	if (nj-1>=0 && nk-1>=0){
		r += nj*nk*threeov(alpha,beta,gamma,A,B,C,li,mi,ni,lj,mj,nj-1,lk,mk,nk-1); 
	}
	if (nk-1>=0){
		r += MINUS_TWO*beta*nk*threeov(alpha,beta,gamma,A,B,C,li,mi,ni,lj,mj,nj+1,lk,mk,nk-1); 
	}
	if (nj-1>=0){
		r += MINUS_TWO*gamma*nj*threeov(alpha,beta,gamma,A,B,C,li,mi,ni,lj,mj,nj-1,lk,mk,nk+1); 
	}
	r += FOUR*beta*gamma*threeov(alpha,beta,gamma,A,B,C,li,mi,ni,lj,mj,nj+1,lk,mk,nk+1);


	r *= ONE/TWO;
	return r;
}

void tki::directtki(const Int& Li, const Int& Lj, const Int& Lk, 
		const Int& inp, const Double* icoe, const Double* iexp, const Double* A, 
		const Int& jnp, const Double* jcoe, const Double* jexp, const Double* B, 
		const Int& knp, const Double* kcoe, const Double* kexp, const Double* C, Double* abcd)
{
	// number of basis sets
	Int nBra1Bas = getCartBas(Li,Li);
	Int nBra2Bas = getCartBas(Lj,Lj);
	Int nKet1Bas = getCartBas(Lk,Lk);

	// loop over each basis set
	Int count = 0;
	for(Int kBas=0; kBas<nKet1Bas; kBas++) {
		for(Int jBas=0; jBas<nBra2Bas; jBas++) {
			for(Int iBas=0; iBas<nBra1Bas; iBas++) {

				// get the angular momentum
				Int l1,m1,n1,l2,m2,n2,l3,m3,n3;
				getlmn(Li, Li, iBas, l1, m1, n1);
				getlmn(Lj, Lj, jBas, l2, m2, n2);
				getlmn(Lk, Lk, kBas, l3, m3, n3);

				// form the primitives data
				Double ki = ZERO;
				for(Int kp=0; kp<knp; kp++) {
					for(Int jp=0; jp<jnp; jp++) {
						for(Int ip=0; ip<inp; ip++) {
							Double ia = iexp[ip];
							Double ja = jexp[jp];
							Double ka = kexp[kp];
							Double ic = icoe[ip];
							Double jc = jcoe[jp];
							Double kc = kcoe[kp];
							Double r  = threeki(ia,ja,ka,A,B,C,  
									l1,m1,n1,l2,m2,n2,l3,m3,n3);
							ki += ic*jc*kc*r;
						}
					}
				}

				// push it to the abcd
				abcd[count] = ki;
				count++;
			}
		}
	}
}

void tki::tki_test(const Int& maxL, 
		const Int& inp, const vector<Double>& icoe, const vector<Double>& iexp, const Double* A, 
		const Int& jnp, const vector<Double>& jcoe, const vector<Double>& jexp, const Double* B,
		const Int& knp, const vector<Double>& kcoe, const vector<Double>& kexp, const Double* C)
{
	//
	// form the data for hgp calculation
	// bra side
	//
	Double AB2 = (A[0]-B[0])*(A[0]-B[0])+(A[1]-B[1])*(A[1]-B[1])+(A[2]-B[2])*(A[2]-B[2]);
	Int inp2 = inp*jnp;
	vector<Double> iexp2(inp2,ZERO);
	vector<Double> iexpdiff(inp2,ZERO);
	vector<Double> fbra(inp2,ZERO);
	vector<Double> P(3*inp2,ZERO);
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

	//
	// form the data for hgp calculation
	// ket side
	//
	count = 0;
	vector<Double> jexp2(knp,ZERO);
	for(Int kp=0; kp<knp; kp++) {
		Double ia = kexp[kp];
		jexp2[count]= ONE/ia;
		count++;
	}

	// now it's real work to do ki test
	cout << "**************************************************************" << endl;
	cout << " three body kinetic integrals test:" << endl;
	cout << "**************************************************************" << endl;
	for(Int kpos=0; kpos<MAX_SHELL_TYPES; kpos++) {
		for(Int jpos=0; jpos<MAX_SHELL_TYPES; jpos++) {
			for(Int ipos=jpos; ipos<MAX_SHELL_TYPES; ipos++) {

				// now decode the L
				Int iLmin, iLmax;
				decodeL(SHELL_ANG_MOM_CODE[ipos],iLmin,iLmax);
				Int jLmin, jLmax;
				decodeL(SHELL_ANG_MOM_CODE[jpos],jLmin,jLmax);
				Int kLmin, kLmax;
				decodeL(SHELL_ANG_MOM_CODE[kpos],kLmin,kLmax);
				if (iLmax>maxL || jLmax>maxL || kLmax>maxL) continue;

				// debug test
				//if (iLmin != iLmax) continue;
				//if (jLmin != jLmax) continue;
				//if (kLmin != kLmax) continue;
				//Int totalL = iLmax+jLmax+kLmax;
				//if (totalL>=3) continue;
				//if (iLmin != 1 && iLmax != 1) continue;
				//if (jLmin != 1 && jLmax != 1) continue;
				//if (kLmin != 0 && kLmax != 0) continue;

				// now form the L code
				LInt LCode = codeSQ(SHELL_ANG_MOM_CODE[ipos],
						SHELL_ANG_MOM_CODE[jpos],SHELL_ANG_MOM_CODE[kpos]);

				// normalize the coefficient array
				Int nLBra1 = iLmax-iLmin+1;
				Int nLBra2 = jLmax-jLmin+1;
				Int nLKet1 = kLmax-kLmin+1;
				vector<Double> bra1Coe(nLBra1*inp);
				vector<Double> bra2Coe(nLBra2*jnp);
				vector<Double> ket1Coe(nLKet1*knp);
				normCoe(iLmin, iLmax, iexp, icoe, bra1Coe);
				normCoe(jLmin, jLmax, jexp, jcoe, bra2Coe);
				normCoe(kLmin, kLmax, kexp, kcoe, ket1Coe);

				// generate the shell pair for coe
				vector<Double> braCoePair(nLBra1*nLBra2*inp2);
				makeC2(iLmin,iLmax,jLmin,jLmax,inp,jnp,bra1Coe,bra2Coe,braCoePair);

				// result array
				Int nBra1Bas = getCartBas(iLmin,iLmax);
				Int nBra2Bas = getCartBas(jLmin,jLmax);
				Int nKet1Bas = getCartBas(kLmin,kLmax);
				vector<Double> result(nBra1Bas*nBra2Bas*nKet1Bas);
				UInt inp2_   = static_cast<UInt>(inp2);
				UInt knp_    = static_cast<UInt>(knp);
				hgp_os_threebodyki(LCode,inp2_,knp_,&braCoePair.front(),&iexp2.front(), 
						&iexpdiff.front(),&fbra.front(),&P.front(),A,B,&ket1Coe.front(),
						&jexp2.front(),C,&result.front());

				// now let's directly calculate the ov
				Int nCarBas1 = getCartBas(iLmin,iLmax);
				Int nCarBas2 = getCartBas(jLmin,jLmax);
				Int d1       = nCarBas1;
				Int d2       = nCarBas1*nCarBas2;
				for(Int Lk=kLmin; Lk<=kLmax; Lk++) {
					for(Int Lj=jLmin; Lj<=jLmax; Lj++) {
						for(Int Li=iLmin; Li<=iLmax; Li++) {

							// initilize the results
							Int nBra1Bas = getCartBas(Li,Li);
							Int nBra2Bas = getCartBas(Lj,Lj);
							Int nKet1Bas = getCartBas(Lk,Lk);
							vector<Double> abcd(nBra1Bas*nBra2Bas*nKet1Bas);

							// get the coefficient array
							const Double* iCArray = &bra1Coe[(Li-iLmin)*inp];
							const Double* jCArray = &bra2Coe[(Lj-jLmin)*jnp];
							const Double* kCArray = &ket1Coe[(Lk-kLmin)*knp];

							// now do direct ov
							directtki(Li,Lj,Lk,inp,iCArray,&iexp.front(),A,jnp,jCArray,
									&jexp.front(),B,knp,kCArray,&kexp.front(),C,&abcd.front());

							// now it's comparison
							for(Int kBas=0; kBas<nKet1Bas; kBas++) {
								for(Int jBas=0; jBas<nBra2Bas; jBas++) {
									for(Int iBas=0; iBas<nBra1Bas; iBas++) {
										Int iBasIndex = getBasOffset(iLmin,Li,iBas);
										Int jBasIndex = getBasOffset(jLmin,Lj,jBas);
										Int kBasIndex = getBasOffset(kLmin,Lk,kBas);
										Int index = iBasIndex+jBasIndex*d1+kBasIndex*d2;
										Double v1 = result[index];
										Double v2 = abcd[iBas+jBas*nBra1Bas+kBas*nBra1Bas*nBra2Bas];
										if (fabs(v1-v2)>THRESH) {
											cout << "Bra1's L: " << iLmin << " " << iLmax << endl;
											cout << "Bra2's L: " << jLmin << " " << jLmax << endl;
											cout << "Ket1's L: " << kLmin << " " << kLmax << endl;
											cout << "result did not match: " << endl;
											printf("difference  : %-16.10f\n", fabs(v1-v2));
											printf("hgp    value: %-16.10f\n", v1);
											printf("direct value: %-16.10f\n", v2);
											Int l1,m1,n1,l2,m2,n2,l3,m3,n3;
											getlmn(Li, Li, iBas, l1, m1, n1);
											getlmn(Lj, Lj, jBas, l2, m2, n2);
											getlmn(Lk, Lk, kBas, l3, m3, n3);
											cout << l1 << m1 << n1 << endl;
											cout << l2 << m2 << n2 << endl;
											cout << l3 << m3 << n3 << endl;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}


