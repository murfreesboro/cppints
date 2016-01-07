#include "libgen.h"
#include "angmomlist.h"
#include "shellprop.h"
#include "functions.h"
#include "norm.h"
#include "tov.h"
#include "localmemscr.h"
#include "mom.h"
using namespace shellprop;
using namespace functions;
using namespace localmemscr;
using namespace norm;
using namespace tov;
using namespace mom;

extern void hgp_os_mom(const LInt& LCode, const UInt& inp2, const Double* icoe, 
		const Double* iexp, const Double* ifac, const Double* P, 
		const Double* A, const Double* B, const Double* C, Double* abcd);

void mom::directmom(const Int& Li, const Int& Lj, const Int& Lk, 
		const Int& inp, const Double* icoe, const Double* iexp, const Double* A, 
		const Int& jnp, const Double* jcoe, const Double* jexp, const Double* B, 
		const Double* C, Double* abcd)
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
				// here we will call three ov directly
				// for mom integrals, the gamma is 0
				Double ov = ZERO;
				for(Int jp=0; jp<jnp; jp++) {
					for(Int ip=0; ip<inp; ip++) {
						Double ia = iexp[ip];
						Double ja = jexp[jp];
						Double ka = ZERO;
						Double ic = icoe[ip];
						Double jc = jcoe[jp];
						Double r  = threeov(ia,ja,ka,A,B,C,  
								l1,m1,n1,l2,m2,n2,l3,m3,n3);
						ov += ic*jc*r;
					}
				}

				// push it to the abcd
				abcd[count] = ov;
				count++;
			}
		}
	}
}

void mom::mom_test(const Int& maxL, const Int& auxL, 
		const Int& inp, const vector<Double>& icoe, const vector<Double>& iexp, const Double* A, 
		const Int& jnp, const vector<Double>& jcoe, const vector<Double>& jexp, const Double* B, 
		const Double* C)
{
	//
	// form the data for hgp calculation
	// bra side
	//
	Double AB2 = (A[0]-B[0])*(A[0]-B[0])+(A[1]-B[1])*(A[1]-B[1])+(A[2]-B[2])*(A[2]-B[2]);
	Int inp2 = inp*jnp;
	vector<Double> iexp2(inp2,ZERO);
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
	// set up the local mem scr
	// the length is set according to maxL = 5
	// auxMaxL = 5
	//
	//LocalMemScr scr(92610);

	// now it's real work to do ov test
	cout << "**************************************************************" << endl;
	cout << " moment overlap integrals test:" << endl;
	cout << "**************************************************************" << endl;

	// here the Loper is closely related to the code in cppints
	for(Int LOper=1; LOper<=auxL; LOper++) {
		for(Int jpos=0; jpos<MAX_SHELL_TYPES; jpos++) {
			for(Int ipos=jpos; ipos<MAX_SHELL_TYPES; ipos++) {

				// now decode the L
				Int iLmin, iLmax;
				decodeL(SHELL_ANG_MOM_CODE[ipos],iLmin,iLmax);
				Int jLmin, jLmax;
				decodeL(SHELL_ANG_MOM_CODE[jpos],jLmin,jLmax);
				if (iLmax>maxL || jLmax>maxL) continue;

				// now form the L code
				LInt LCode = codeSQ(SHELL_ANG_MOM_CODE[ipos], 
						SHELL_ANG_MOM_CODE[jpos], LOper);

				// normalize the coefficient array
				Int nLBra1 = iLmax-iLmin+1;
				Int nLBra2 = jLmax-jLmin+1;
				vector<Double> bra1Coe(nLBra1*inp);
				vector<Double> bra2Coe(nLBra2*jnp);
				normCoe(iLmin, iLmax, iexp, icoe, bra1Coe);
				normCoe(jLmin, jLmax, jexp, jcoe, bra2Coe);

				// generate the shell pair for coe
				vector<Double> braCoePair(nLBra1*nLBra2*inp2);
				makeC2(iLmin,iLmax,jLmin,jLmax,inp,jnp,bra1Coe,bra2Coe,braCoePair);

				// result array
				Int nBra1Bas = getCartBas(iLmin,iLmax);
				Int nBra2Bas = getCartBas(jLmin,jLmax);
				Int nKet1Bas = getCartBas(LOper,LOper);
				UInt inp2_   = static_cast<UInt>(inp2);
				vector<Double> result(nBra1Bas*nBra2Bas*nKet1Bas);
				hgp_os_mom(LCode,inp2_,&braCoePair.front(),&iexp2.front(),
						&fbra.front(),&P.front(),A,B,C,&result.front());

				// reset the scr
				//scr.reset();

				// now let's directly calculate the ov
				Int nCarBas1 = getCartBas(iLmin,iLmax);
				Int nCarBas2 = getCartBas(jLmin,jLmax);
				Int d1       = nCarBas1;
				Int d2       = nCarBas1*nCarBas2;
				for(Int Lj=jLmin; Lj<=jLmax; Lj++) {
					for(Int Li=iLmin; Li<=iLmax; Li++) {

						// initilize the results
						Int nBra1Bas = getCartBas(Li,Li);
						Int nBra2Bas = getCartBas(Lj,Lj);
						vector<Double> abcd(nBra1Bas*nBra2Bas*nKet1Bas);

						// get the coefficient array
						const Double* iCArray = &bra1Coe[(Li-iLmin)*inp];
						const Double* jCArray = &bra2Coe[(Lj-jLmin)*jnp];

						// now do direct ov
						directmom(Li,Lj,LOper,inp,iCArray,&iexp.front(),A,jnp,
								jCArray,&jexp.front(),B,C,&abcd.front());

						// now it's comparison
						for(Int kBas=0; kBas<nKet1Bas; kBas++) {
							for(Int jBas=0; jBas<nBra2Bas; jBas++) {
								for(Int iBas=0; iBas<nBra1Bas; iBas++) {
									Int iBasIndex = getBasOffset(iLmin,Li,iBas);
									Int jBasIndex = getBasOffset(jLmin,Lj,jBas);
									Int kBasIndex = getBasOffset(LOper,LOper,kBas);
									Int index = iBasIndex+jBasIndex*d1+kBasIndex*d2;
									Double v1 = result[index];
									Double v2 = abcd[iBas+jBas*nBra1Bas+kBas*nBra1Bas*nBra2Bas];
									if (fabs(v1-v2)>THRESH) {
										cout << "Bra1's L: " << iLmin << " " << iLmax << endl;
										cout << "Bra2's L: " << jLmin << " " << jLmax << endl;
										cout << "L Oper  : " << LOper << endl;
										cout << "result did not match: " << endl;
										printf("difference  : %-16.10f\n", fabs(v1-v2));
										printf("hgp    value: %-16.10f\n", v1);
										printf("direct value: %-16.10f\n", v2);
										Int l1,m1,n1,l2,m2,n2,l3,m3,n3;
										getlmn(Li, Li, iBas, l1, m1, n1);
										getlmn(Lj, Lj, jBas, l2, m2, n2);
										getlmn(LOper, LOper, kBas, l3, m3, n3);
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


