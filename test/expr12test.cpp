#include "libgen.h"
#include "angmomlist.h"
#include "shellprop.h"
#include "functions.h"
#include "norm.h"
#include "localmemscr.h"
#include "tov.h"
#include "expr12test.h"
using namespace shellprop;
using namespace functions;
using namespace norm;
using namespace tov;
using namespace localmemscr;
using namespace expr12test;

extern void hgp_os_expr12(const LInt& LCode, const UInt& inp2, const UInt& jnp2, const Double& thresh, 
		const Double& omega, const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, 
		const Double* A, const Double* B, const Double* jcoe, const Double* jexp, const Double* jfac, 
		const Double* Q, const Double* C, const Double* D, Double* abcd, LocalMemScr& scr);

void expr12test::expr12_test(const Int& maxL, 
		const Int& inp, const vector<Double>& icoe, const vector<Double>& iexp, const Double* A, 
		const Int& jnp, const vector<Double>& jcoe, const vector<Double>& jexp, const Double* B, 
		const Int& knp, const vector<Double>& kcoe, const vector<Double>& kexp, const Double* C) 
{
	//
	// form the data for hgp calculation
	// firstly it's the 
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

	// let's define the dummy center and it's data
	// the dummy shell will be in the last dimension of the 
	// shell quartet
	Int lnp = 1;
	Double lcoe[1];
	lcoe[0] = ONE;
	Double lexp[1];
	lexp[0] = ZERO;
	Double D[3];
	D[0] = 0.0E0;
	D[1] = 0.0E0;
	D[2] = 0.0E0;

	// now it's ket side data
	Double CD2 = (C[0]-D[0])*(C[0]-D[0])+(C[1]-D[1])*(C[1]-D[1])+(C[2]-D[2])*(C[2]-D[2]);
	Int jnp2 = knp*lnp;
	vector<Double> jexp2(jnp2,ZERO);
	vector<Double> fket(jnp2,ZERO);
	vector<Double> Q(3*jnp2,ZERO);
	count = 0;
	for(Int lp=0; lp<lnp; lp++) {
		for(Int kp=0; kp<knp; kp++) {

			// prefactors etc.
			Double ia   = kexp[kp];
			Double ja   = lexp[lp];
			Double alpla= ia+ja; 
			Double ab   = -ia*ja/alpla;
			Double pref = exp(ab*CD2)*pow(PI/alpla,1.5E0);
			jexp2[count]= ONE/alpla;
			fket[count] = pref;

			// form P point according to the 
			// Gaussian pritimive product theorem
			Double adab = ia/alpla; 
			Double bdab = ja/alpla; 
			Double Qx   = C[0]*adab + D[0]*bdab;
			Double Qy   = C[1]*adab + D[1]*bdab;
			Double Qz   = C[2]*adab + D[2]*bdab;
			Q[3*count+0]= Qx;
			Q[3*count+1]= Qy;
			Q[3*count+2]= Qz;
			count++;
		}
	}

	//
	// set up the local mem scr 
	// the information is depending on maxL = 5
	// auxMaxL = 5
	// for <HH|HH> shell quartet, the maximum memory
	// length is 194481
	//
	LocalMemScr scr(1944810);

	//
	// real work begins
	// 
	cout << "**************************************************************" << endl;
	cout << "* four body exp(-omega*r12^2) integral test:" << endl;
	cout << "* the last dimension is always the dummy shell " << endl;
	cout << "**************************************************************" << endl;
	for(Int n2=0; n2<MAX_SHELL_PAIR_NUMBER; n2++) {
		for(Int n1=n2; n1<MAX_SHELL_PAIR_NUMBER; n1++) {

			// the last shell must be S
			Int LBra  = SHELL_PAIR_ORDER_ARRAY[n1];
			Int LKet  = SHELL_PAIR_ORDER_ARRAY[n2];
			Int iLmin = -1;
			Int iLmax = -1;
			Int jLmin = -1;
			Int jLmax = -1;
			decodeSQ(LBra,iLmin,iLmax,jLmin,jLmax);
			if (iLmax>maxL) continue;
			if (jLmax>maxL) continue;
			Int kLmin = -1;
			Int kLmax = -1;
			Int lLmin = -1;
			Int lLmax = -1;
			decodeSQ(LKet,kLmin,kLmax,lLmin,lLmax);
			if (kLmax>maxL) continue;
			if (lLmax != S) continue;
			Int L1 = codeL(iLmin,iLmax);
			Int L2 = codeL(jLmin,jLmax);
			Int L3 = codeL(kLmin,kLmax);
			Int L4 = codeL(lLmin,lLmax);

			// now form the L code
			LInt LCode = codeSQ(L1,L2,L3,L4);

			// normalize the coefficient array
			Int nLBra1 = iLmax-iLmin+1;
			Int nLBra2 = jLmax-jLmin+1;
			Int nLKet1 = kLmax-kLmin+1;
			Int nLKet2 = lLmax-lLmin+1;
			vector<Double> bra1Coe(nLBra1*inp);
			vector<Double> bra2Coe(nLBra2*jnp);
			vector<Double> ket1Coe(nLKet1*knp);
			vector<Double> ket2Coe(lnp,ONE);
			normCoe(iLmin, iLmax, iexp, icoe, bra1Coe);
			normCoe(jLmin, jLmax, jexp, jcoe, bra2Coe);
			normCoe(kLmin, kLmax, kexp, kcoe, ket1Coe);

			// generate the shell pair for coe
			vector<Double> braCoePair(nLBra1*nLBra2*inp2);
			makeC2(iLmin,iLmax,jLmin,jLmax,inp,jnp,bra1Coe,bra2Coe,braCoePair);
			vector<Double> ketCoePair(nLKet1*nLKet2*jnp2);
			makeC2(kLmin,kLmax,lLmin,lLmax,knp,lnp,ket1Coe,ket2Coe,ketCoePair);

			// result array
			// set omega to be a very large number
			Int nBra1Bas = getCartBas(iLmin,iLmax);
			Int nBra2Bas = getCartBas(jLmin,jLmax);
			Int nKet1Bas = getCartBas(kLmin,kLmax);
			Int nKet2Bas = getCartBas(lLmin,lLmax);
			vector<Double> result(nBra1Bas*nBra2Bas*nKet1Bas*nKet2Bas);
			Double omega = 1000000.0E0;
			hgp_os_expr12(LCode,inp2,jnp2,INTEGRAL_THRESH,omega,
					&braCoePair.front(),&iexp2.front(),&fbra.front(),&P.front(),A,B, 
					&ketCoePair.front(),&jexp2.front(),&fket.front(),&Q.front(),C,D, 
					&result.front(),scr);

			// clear the scr
			scr.reset();

			// now let's directly calculate the three body ov
			// for large omega value, it should become the 
			// three body overlap integral
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
						directtov(Li,Lj,Lk,inp,iCArray,&iexp.front(),A,jnp,jCArray,
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
									if (fabs(v1-v2)>1.0E-7) {
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


