#include "libgen.h"
#include "angmomlist.h"
#include "shellprop.h"
#include "functions.h"
#include "norm.h"
#include "eri.h"
#include "derivinfor.h"
#include "localmemscr.h"
#include "eritest.h"
using namespace shellprop;
using namespace functions;
using namespace norm;
using namespace localmemscr;
using namespace derivinfor;
using namespace eritest;

#ifdef ORDER1
void hgp_os_eri_d1(const LInt& LCode, const UInt& inp2, const UInt& jnp2, 
		const Double& thresh, const Double& pMax, const Double& omega, const Double* icoe, 
		const Double* iexp, const Double* iexpdiff, const Double* ifac, 
		const Double* P, const Double* A, const Double* B, const Double* jcoe, 
		const Double* jexp, const Double* jexpdiff, const Double* jfac, 
		const Double* Q, const Double* C, const Double* D, Double* abcd, LocalMemScr& scr);
#endif

#ifdef ORDER2
void hgp_os_eri_d2(const LInt& LCode, const UInt& inp2, const UInt& jnp2, 
		const Double& thresh, const Double& pMax, const Double& omega, const Double* icoe, 
		const Double* iexp, const Double* iexpdiff, const Double* ifac, 
		const Double* P, const Double* A, const Double* B, const Double* jcoe, 
		const Double* jexp, const Double* jexpdiff, const Double* jfac, 
		const Double* Q, const Double* C, const Double* D, Double* abcd, LocalMemScr& scr);
#endif

void eritest::directeri_d1(const Int& Li, const Int& Lj, const Int& Lk, const Int& Ll,
		const Int& inp, const Double* icoe, const Double* iexp, const Double* A, 
		const Int& jnp, const Double* jcoe, const Double* jexp, const Double* B, 
		const Int& knp, const Double* kcoe, const Double* kexp, const Double* C, 
		const Int& lnp, const Double* lcoe, const Double* lexp, const Double* D, 
		const Int& pos, const Int& dir, Double* abcd)
{
	// number of basis sets
	Int nBra1Bas = getCartBas(Li,Li);
	Int nBra2Bas = getCartBas(Lj,Lj);
	Int nKet1Bas = getCartBas(Lk,Lk);
	Int nKet2Bas = getCartBas(Ll,Ll);

	// re-define the center
	double AA[3];
	AA[0] = A[0];
	AA[1] = A[1];
	AA[2] = A[2];
	double BB[3];
	BB[0] = B[0];
	BB[1] = B[1];
	BB[2] = B[2];
	double CC[3];
	CC[0] = C[0];
	CC[1] = C[1];
	CC[2] = C[2];
	double DD[3];
	DD[0] = D[0];
	DD[1] = D[1];
	DD[2] = D[2];

	// loop over each basis set
	Int count = 0;
	for(Int lBas=0; lBas<nKet2Bas; lBas++) {
		for(Int kBas=0; kBas<nKet1Bas; kBas++) {
			for(Int jBas=0; jBas<nBra2Bas; jBas++) {
				for(Int iBas=0; iBas<nBra1Bas; iBas++) {

					// get the angular momentum
					Int ll1,mm1,nn1,ll2,mm2,nn2;
					Int ll3,mm3,nn3,ll4,mm4,nn4;
					getlmn(Li, Li, iBas, ll1, mm1, nn1);
					getlmn(Lj, Lj, jBas, ll2, mm2, nn2);
					getlmn(Lk, Lk, kBas, ll3, mm3, nn3);
					getlmn(Ll, Ll, lBas, ll4, mm4, nn4);

					// form the primitives data
					Double result = ZERO;
					for(Int lp=0; lp<lnp; lp++) {
						for(Int kp=0; kp<knp; kp++) {
							for(Int jp=0; jp<jnp; jp++) {
								for(Int ip=0; ip<inp; ip++) {

									// we do not have coefficient in eri calculation
									// therefore only fetch the exponentials
									// also for either float/double type, because the 
									// function of eri only support double; therefore
									// we use all double type here
									double ia = iexp[ip];
									double ja = jexp[jp];
									double ka = kexp[kp];
									double la = lexp[lp];
									double ic = icoe[ip];
									double jc = jcoe[jp];
									double kc = kcoe[kp];
									double lc = lcoe[lp];

									// here make a copy of original l,m,n values
									Int l1 = ll1;
									Int l2 = ll2;
									Int l3 = ll3;
									Int l4 = ll4;
									Int m1 = mm1;
									Int m2 = mm2;
									Int m3 = mm3;
									Int m4 = mm4;
									Int n1 = nn1;
									Int n2 = nn2;
									Int n3 = nn3;
									Int n4 = nn4;

									// first term, do derivatives on x^ly^mz^n*exp(-alpha*r^2)
									// this will be, for example if on x;
									// -l*x^l-1y^mz^n*exp(-alpha*r^2)
									int newL = 1;
									if (pos == BRA1) {
										if (dir == DERIV_X) {
											newL = l1;
											l1 = l1-1;
										}else if (dir == DERIV_Y) {
											newL = m1;
											m1 = m1-1;
										}else {
											newL = n1;
											n1 = n1-1;
										}
									}else if (pos == BRA2) {
										if (dir == DERIV_X) {
											newL = l2;
											l2 = l2-1;
										}else if (dir == DERIV_Y) {
											newL = m2;
											m2 = m2-1;
										}else {
											newL = n2;
											n2 = n2-1;
										}
									}else if (pos == KET1) {
										if (dir == DERIV_X) {
											newL = l3;
											l3 = l3-1;
										}else if (dir == DERIV_Y) {
											newL = m3;
											m3 = m3-1;
										}else {
											newL = n3;
											n3 = n3-1;
										}
									}else{
										if (dir == DERIV_X) {
											newL = l4;
											l4 = l4-1;
										}else if (dir == DERIV_Y) {
											newL = m4;
											m4 = m4-1;
										}else {
											newL = n4;
											n4 = n4-1;
										}
									}
									double r1 = ZERO;
									if (newL>=1) {
										double r = eri(l1,m1,n1,ia,AA,
												l2,m2,n2,ja,BB,
												l3,m3,n3,ka,CC,
												l4,m4,n4,la,DD, 0);
										r1 = -1.0E0*newL*ic*jc*kc*lc*r;
									}

									// second term, do derivatives on x^ly^mz^n*exp(-alpha*r^2)
									// this will be, for example if on x;
									// 2*alpha*x^l+1y^mz^n*exp(-alpha*r^2)
									// copy the original l,m,n values first
									l1 = ll1;
									l2 = ll2;
									l3 = ll3;
									l4 = ll4;
									m1 = mm1;
									m2 = mm2;
									m3 = mm3;
									m4 = mm4;
									n1 = nn1;
									n2 = nn2;
									n3 = nn3;
									n4 = nn4;

									// we note that this term always exist
									double alpha = ZERO;
									if (pos == BRA1) {
										alpha = ia;
										if (dir == DERIV_X) {
											l1 = l1+1;
										}else if (dir == DERIV_Y) {
											m1 = m1+1;
										}else {
											n1 = n1+1;
										}
									}else if (pos == BRA2) {
										alpha = ja;
										if (dir == DERIV_X) {
											l2 = l2+1;
										}else if (dir == DERIV_Y) {
											m2 = m2+1;
										}else {
											n2 = n2+1;
										}
									}else if (pos == KET1) {
										alpha = ka;
										if (dir == DERIV_X) {
											l3 = l3+1;
										}else if (dir == DERIV_Y) {
											m3 = m3+1;
										}else {
											n3 = n3+1;
										}
									}else{
										alpha = la;
										if (dir == DERIV_X) {
											l4 = l4+1;
										}else if (dir == DERIV_Y) {
											m4 = m4+1;
										}else {
											n4 = n4+1;
										}
									}
									double r2 = ZERO;
									double tmp_r2 = eri(l1,m1,n1,ia,AA,
											l2,m2,n2,ja,BB,
											l3,m3,n3,ka,CC,
											l4,m4,n4,la,DD, 0);
									r2 = 2.0E0*alpha*ic*jc*kc*lc*tmp_r2;

									// now it's final result
									result += r1 + r2;
								}
							}
						}
					}

					// push it to the abcd
					abcd[count] = result;
					count++;
				}
			}
		}
	}
}

//  this is for derivatives on the same position and same direction
void eritest::directeri_d2(const Int& Li, const Int& Lj, const Int& Lk, const Int& Ll,
		const Int& inp, const Double* icoe, const Double* iexp, const Double* A, 
		const Int& jnp, const Double* jcoe, const Double* jexp, const Double* B, 
		const Int& knp, const Double* kcoe, const Double* kexp, const Double* C, 
		const Int& lnp, const Double* lcoe, const Double* lexp, const Double* D, 
		const Int& pos, const Int& dir, Double* abcd)
{
	// number of basis sets
	Int nBra1Bas = getCartBas(Li,Li);
	Int nBra2Bas = getCartBas(Lj,Lj);
	Int nKet1Bas = getCartBas(Lk,Lk);
	Int nKet2Bas = getCartBas(Ll,Ll);

	// re-define the center
	double AA[3];
	AA[0] = A[0];
	AA[1] = A[1];
	AA[2] = A[2];
	double BB[3];
	BB[0] = B[0];
	BB[1] = B[1];
	BB[2] = B[2];
	double CC[3];
	CC[0] = C[0];
	CC[1] = C[1];
	CC[2] = C[2];
	double DD[3];
	DD[0] = D[0];
	DD[1] = D[1];
	DD[2] = D[2];

	// loop over each basis set
	Int count = 0;
	for(Int lBas=0; lBas<nKet2Bas; lBas++) {
		for(Int kBas=0; kBas<nKet1Bas; kBas++) {
			for(Int jBas=0; jBas<nBra2Bas; jBas++) {
				for(Int iBas=0; iBas<nBra1Bas; iBas++) {

					// get the angular momentum
					Int ll1,mm1,nn1,ll2,mm2,nn2;
					Int ll3,mm3,nn3,ll4,mm4,nn4;
					getlmn(Li, Li, iBas, ll1, mm1, nn1);
					getlmn(Lj, Lj, jBas, ll2, mm2, nn2);
					getlmn(Lk, Lk, kBas, ll3, mm3, nn3);
					getlmn(Ll, Ll, lBas, ll4, mm4, nn4);

					// form the primitives data
					Double result = ZERO;
					for(Int lp=0; lp<lnp; lp++) {
						for(Int kp=0; kp<knp; kp++) {
							for(Int jp=0; jp<jnp; jp++) {
								for(Int ip=0; ip<inp; ip++) {

									// we do not have coefficient in eri calculation
									// therefore only fetch the exponentials
									// also for either float/double type, because the 
									// function of eri only support double; therefore
									// we use all double type here
									double ia = iexp[ip];
									double ja = jexp[jp];
									double ka = kexp[kp];
									double la = lexp[lp];
									double ic = icoe[ip];
									double jc = jcoe[jp];
									double kc = kcoe[kp];
									double lc = lcoe[lp];

									// here make a copy of original l,m,n values
									Int l1 = ll1;
									Int l2 = ll2;
									Int l3 = ll3;
									Int l4 = ll4;
									Int m1 = mm1;
									Int m2 = mm2;
									Int m3 = mm3;
									Int m4 = mm4;
									Int n1 = nn1;
									Int n2 = nn2;
									Int n3 = nn3;
									Int n4 = nn4;

									// first term, do derivatives on x^ly^mz^n*exp(-alpha*r^2)
									// for example, if derivatives is on x; first on exp then
									// on the angular part, it gives
									// (-2*alpha*(l+1))*x^ly^mz^n*exp(-alpha*r^2)
									// this term always exist
									double alpha = ZERO;
									Int newL = 1;
									if (pos == BRA1) {
										alpha = ia;
										if (dir == DERIV_X) {
											newL = l1;
										}else if (dir == DERIV_Y) {
											newL = m1;
										}else {
											newL = n1;
										}
									}else if (pos == BRA2) {
										alpha = ja;
										if (dir == DERIV_X) {
											newL = l2;
										}else if (dir == DERIV_Y) {
											newL = m2;
										}else {
											newL = n2;
										}
									}else if (pos == KET1) {
										alpha = ka;
										if (dir == DERIV_X) {
											newL = l3;
										}else if (dir == DERIV_Y) {
											newL = m3;
										}else {
											newL = n3;
										}
									}else{
										alpha = la;
										if (dir == DERIV_X) {
											newL = l4;
										}else if (dir == DERIV_Y) {
											newL = m4;
										}else {
											newL = n4;
										}
									}
									double tmp_r1 = eri(l1,m1,n1,ia,AA,
											l2,m2,n2,ja,BB,
											l3,m3,n3,ka,CC,
											l4,m4,n4,la,DD, 0);
									double r1 = -2.0E0*alpha*(newL+1)*ic*jc*kc*lc*tmp_r1;

									// second term, do derivatives on x^ly^mz^n*exp(-alpha*r^2)
									// for example, if derivatives is on x; firstly on the 
									// angular part then on the exp term, it gives
									// (-2*alpha*l)*x^ly^mz^n*exp(-alpha*r^2)
									// this term may not exist if l = 0
									double r2 = -2.0E0*alpha*newL*ic*jc*kc*lc*tmp_r1;

									// third term, do derivatives on x^ly^mz^n*exp(-alpha*r^2)
									// for example, if derivatives is on x; 
									// this is all on the angular part. it gives:
									// l*(l-1)x^l-2y^mz^n*exp(-alpha*r^2)
									newL = 1;
									if (pos == BRA1) {
										if (dir == DERIV_X) {
											newL = l1;
											l1 = l1-2;
										}else if (dir == DERIV_Y) {
											newL = m1;
											m1 = m1-2;
										}else {
											newL = n1;
											n1 = n1-2;
										}
									}else if (pos == BRA2) {
										if (dir == DERIV_X) {
											newL = l2;
											l2 = l2-2;
										}else if (dir == DERIV_Y) {
											newL = m2;
											m2 = m2-2;
										}else {
											newL = n2;
											n2 = n2-2;
										}
									}else if (pos == KET1) {
										if (dir == DERIV_X) {
											newL = l3;
											l3 = l3-2;
										}else if (dir == DERIV_Y) {
											newL = m3;
											m3 = m3-2;
										}else {
											newL = n3;
											n3 = n3-2;
										}
									}else{
										if (dir == DERIV_X) {
											newL = l4;
											l4 = l4-2;
										}else if (dir == DERIV_Y) {
											newL = m4;
											m4 = m4-2;
										}else {
											newL = n4;
											n4 = n4-2;
										}
									}
									double r3 = ZERO;
									if (newL>1) {
										double tmp_r3 = eri(l1,m1,n1,ia,AA,
												l2,m2,n2,ja,BB,
												l3,m3,n3,ka,CC,
												l4,m4,n4,la,DD, 0);
										r3 = newL*(newL-1)*ic*jc*kc*lc*tmp_r3;
									}

									// fourth term, do derivatives on x^ly^mz^n*exp(-alpha*r^2)
									// and all on exp term. this will be, for example if on x;
									// 4*alpha^2*x^l+2y^mz^n*exp(-alpha*r^2)
									// copy the original l,m,n values first
									l1 = ll1;
									l2 = ll2;
									l3 = ll3;
									l4 = ll4;
									m1 = mm1;
									m2 = mm2;
									m3 = mm3;
									m4 = mm4;
									n1 = nn1;
									n2 = nn2;
									n3 = nn3;
									n4 = nn4;

									// we note that this term always exist
									if (pos == BRA1) {
										if (dir == DERIV_X) {
											l1 = l1+2;
										}else if (dir == DERIV_Y) {
											m1 = m1+2;
										}else {
											n1 = n1+2;
										}
									}else if (pos == BRA2) {
										if (dir == DERIV_X) {
											l2 = l2+2;
										}else if (dir == DERIV_Y) {
											m2 = m2+2;
										}else {
											n2 = n2+2;
										}
									}else if (pos == KET1) {
										if (dir == DERIV_X) {
											l3 = l3+2;
										}else if (dir == DERIV_Y) {
											m3 = m3+2;
										}else {
											n3 = n3+2;
										}
									}else{
										if (dir == DERIV_X) {
											l4 = l4+2;
										}else if (dir == DERIV_Y) {
											m4 = m4+2;
										}else {
											n4 = n4+2;
										}
									}
									double tmp_r4 = eri(l1,m1,n1,ia,AA,
											l2,m2,n2,ja,BB,
											l3,m3,n3,ka,CC,
											l4,m4,n4,la,DD, 0);
									double r4 = 4.0E0*alpha*alpha*ic*jc*kc*lc*tmp_r4;

									// now it's final result
									result += r1 + r2 + r3 + r4;
								}
							}
						}
					}

					// push it to the abcd
					abcd[count] = result;
					count++;
				}
			}
		}
	}
}

void eritest::directeri_d2(const Int& Li, const Int& Lj, const Int& Lk, const Int& Ll,
		const Int& inp, const Double* icoe, const Double* iexp, const Double* A, 
		const Int& jnp, const Double* jcoe, const Double* jexp, const Double* B, 
		const Int& knp, const Double* kcoe, const Double* kexp, const Double* C, 
		const Int& lnp, const Double* lcoe, const Double* lexp, const Double* D, 
		const Int& pos1, const Int& pos2, const Int& dir1, const int& dir2, Double* abcd)
{
	// number of basis sets
	Int nBra1Bas = getCartBas(Li,Li);
	Int nBra2Bas = getCartBas(Lj,Lj);
	Int nKet1Bas = getCartBas(Lk,Lk);
	Int nKet2Bas = getCartBas(Ll,Ll);

	// re-define the center
	double AA[3];
	AA[0] = A[0];
	AA[1] = A[1];
	AA[2] = A[2];
	double BB[3];
	BB[0] = B[0];
	BB[1] = B[1];
	BB[2] = B[2];
	double CC[3];
	CC[0] = C[0];
	CC[1] = C[1];
	CC[2] = C[2];
	double DD[3];
	DD[0] = D[0];
	DD[1] = D[1];
	DD[2] = D[2];

	// loop over each basis set
	Int count = 0;
	for(Int lBas=0; lBas<nKet2Bas; lBas++) {
		for(Int kBas=0; kBas<nKet1Bas; kBas++) {
			for(Int jBas=0; jBas<nBra2Bas; jBas++) {
				for(Int iBas=0; iBas<nBra1Bas; iBas++) {

					// get the angular momentum
					Int ll1,mm1,nn1,ll2,mm2,nn2;
					Int ll3,mm3,nn3,ll4,mm4,nn4;
					getlmn(Li, Li, iBas, ll1, mm1, nn1);
					getlmn(Lj, Lj, jBas, ll2, mm2, nn2);
					getlmn(Lk, Lk, kBas, ll3, mm3, nn3);
					getlmn(Ll, Ll, lBas, ll4, mm4, nn4);

					// form the primitives data
					Double result = ZERO;
					for(Int lp=0; lp<lnp; lp++) {
						for(Int kp=0; kp<knp; kp++) {
							for(Int jp=0; jp<jnp; jp++) {
								for(Int ip=0; ip<inp; ip++) {

									// we do not have coefficient in eri calculation
									// therefore only fetch the exponentials
									// also for either float/double type, because the 
									// function of eri only support double; therefore
									// we use all double type here
									double ia = iexp[ip];
									double ja = jexp[jp];
									double ka = kexp[kp];
									double la = lexp[lp];
									double ic = icoe[ip];
									double jc = jcoe[jp];
									double kc = kcoe[kp];
									double lc = lcoe[lp];

									// here make a copy of original l,m,n values
									Int l1 = ll1;
									Int l2 = ll2;
									Int l3 = ll3;
									Int l4 = ll4;
									Int m1 = mm1;
									Int m2 = mm2;
									Int m3 = mm3;
									Int m4 = mm4;
									Int n1 = nn1;
									Int n2 = nn2;
									Int n3 = nn3;
									Int n4 = nn4;

									// first term 4*alpha*beta*I(l+1,l'+1)
									double alpha = ZERO;
									double beta  = ZERO;
									for(Int i=0; i<2; i++) {

										// set the pos and dir
										Int pos = pos1;
										Int dir = dir1;
										if (i == 1) {
											pos = pos2;
											dir = dir2;
										}

										// set the exponential value
										double val = ZERO;

										// now according to the pos and dir
										// change the value
										if (pos == BRA1) {
											val = ia;
											if (dir == DERIV_X) {
												l1 = l1+1;
											}else if (dir == DERIV_Y) {
												m1 = m1+1;
											}else {
												n1 = n1+1;
											}
										}else if (pos == BRA2) {
											val = ja;
											if (dir == DERIV_X) {
												l2 = l2+1;
											}else if (dir == DERIV_Y) {
												m2 = m2+1;
											}else {
												n2 = n2+1;
											}
										}else if (pos == KET1) {
											val = ka;
											if (dir == DERIV_X) {
												l3 = l3+1;
											}else if (dir == DERIV_Y) {
												m3 = m3+1;
											}else {
												n3 = n3+1;
											}
										}else{
											val = la;
											if (dir == DERIV_X) {
												l4 = l4+1;
											}else if (dir == DERIV_Y) {
												m4 = m4+1;
											}else {
												n4 = n4+1;
											}
										}

										// finally assign value
										if (i == 0) {
											alpha = val;
										}else{
											beta  = val;
										}
									}
									double tmp_r1 = eri(l1,m1,n1,ia,AA,
											l2,m2,n2,ja,BB,
											l3,m3,n3,ka,CC,
											l4,m4,n4,la,DD, 0);
									double r1 = 4.0E0*alpha*beta*ic*jc*kc*lc*tmp_r1;

									// second term 
									// l*l'*I(l-1,l'-1)
									l1 = ll1;
									l2 = ll2;
									l3 = ll3;
									l4 = ll4;
									m1 = mm1;
									m2 = mm2;
									m3 = mm3;
									m4 = mm4;
									n1 = nn1;
									n2 = nn2;
									n3 = nn3;
									n4 = nn4;
									Int newL1 = 0;
									Int newL2 = 0;
									for(Int i=0; i<2; i++) {

										// set the pos and dir
										Int pos = pos1;
										Int dir = dir1;
										if (i == 1) {
											pos = pos2;
											dir = dir2;
										}

										// set the value
										Int newL = 0;

										// now according to the pos and dir
										// change the value
										if (pos == BRA1) {
											if (dir == DERIV_X) {
												newL = l1;
												l1 = l1-1;
											}else if (dir == DERIV_Y) {
												newL = m1;
												m1 = m1-1;
											}else {
												newL = n1;
												n1 = n1-1;
											}
										}else if (pos == BRA2) {
											if (dir == DERIV_X) {
												newL = l2;
												l2 = l2-1;
											}else if (dir == DERIV_Y) {
												newL = m2;
												m2 = m2-1;
											}else {
												newL = n2;
												n2 = n2-1;
											}
										}else if (pos == KET1) {
											if (dir == DERIV_X) {
												newL = l3;
												l3 = l3-1;
											}else if (dir == DERIV_Y) {
												newL = m3;
												m3 = m3-1;
											}else {
												newL = n3;
												n3 = n3-1;
											}
										}else{
											if (dir == DERIV_X) {
												newL = l4;
												l4 = l4-1;
											}else if (dir == DERIV_Y) {
												newL = m4;
												m4 = m4-1;
											}else {
												newL = n4;
												n4 = n4-1;
											}
										}

										//assign newL
										if (i == 0) {
											newL1 = newL;
										}else{
											newL2 = newL;
										}
									}
									double r2 = ZERO;
									if (newL1>=1 && newL2>=1) {
										double tmp_r2 = eri(l1,m1,n1,ia,AA,
												l2,m2,n2,ja,BB,
												l3,m3,n3,ka,CC,
												l4,m4,n4,la,DD, 0);
										r2 = newL1*newL2*ic*jc*kc*lc*tmp_r2;
									}

									// third term 
									// -2*alpha*l'*I(l+1,l'-1)
									// l' is just the newL2
									l1 = ll1;
									l2 = ll2;
									l3 = ll3;
									l4 = ll4;
									m1 = mm1;
									m2 = mm2;
									m3 = mm3;
									m4 = mm4;
									n1 = nn1;
									n2 = nn2;
									n3 = nn3;
									n4 = nn4;
									for(Int i=0; i<2; i++) {

										// set the pos and dir
										Int pos = pos1;
										Int dir = dir1;
										if (i == 1) {
											pos = pos2;
											dir = dir2;
										}

										// also, incremental value is different
										// first one is 1, second one is -1
										Int inc = 1;
										if (i == 1) {
											inc = -1;
										}

										// now according to the pos and dir
										// change the value
										if (pos == BRA1) {
											if (dir == DERIV_X) {
												l1 = l1+inc;
											}else if (dir == DERIV_Y) {
												m1 = m1+inc;
											}else {
												n1 = n1+inc;
											}
										}else if (pos == BRA2) {
											if (dir == DERIV_X) {
												l2 = l2+inc;
											}else if (dir == DERIV_Y) {
												m2 = m2+inc;
											}else {
												n2 = n2+inc;
											}
										}else if (pos == KET1) {
											if (dir == DERIV_X) {
												l3 = l3+inc;
											}else if (dir == DERIV_Y) {
												m3 = m3+inc;
											}else {
												n3 = n3+inc;
											}
										}else{
											if (dir == DERIV_X) {
												l4 = l4+inc;
											}else if (dir == DERIV_Y) {
												m4 = m4+inc;
											}else {
												n4 = n4+inc;
											}
										}
									}
									double r3 = ZERO;
									if (newL2>=1) {
										double tmp_r3 = eri(l1,m1,n1,ia,AA,
												l2,m2,n2,ja,BB,
												l3,m3,n3,ka,CC,
												l4,m4,n4,la,DD, 0);
										r3 = -2.0E0*alpha*newL2*ic*jc*kc*lc*tmp_r3;
									}

									// fourth term 
									// -2*beta*l*I(l-1,l'+1)
									// l is just the newL1
									l1 = ll1;
									l2 = ll2;
									l3 = ll3;
									l4 = ll4;
									m1 = mm1;
									m2 = mm2;
									m3 = mm3;
									m4 = mm4;
									n1 = nn1;
									n2 = nn2;
									n3 = nn3;
									n4 = nn4;
									for(Int i=0; i<2; i++) {

										// set the pos and dir
										Int pos = pos1;
										Int dir = dir1;
										if (i == 1) {
											pos = pos2;
											dir = dir2;
										}

										// also, incremental value is different
										// first one is -1, second one is +1
										Int inc = -1;
										if (i == 1) {
											inc = 1;
										}

										// now according to the pos and dir
										// change the value
										if (pos == BRA1) {
											if (dir == DERIV_X) {
												l1 = l1+inc;
											}else if (dir == DERIV_Y) {
												m1 = m1+inc;
											}else {
												n1 = n1+inc;
											}
										}else if (pos == BRA2) {
											if (dir == DERIV_X) {
												l2 = l2+inc;
											}else if (dir == DERIV_Y) {
												m2 = m2+inc;
											}else {
												n2 = n2+inc;
											}
										}else if (pos == KET1) {
											if (dir == DERIV_X) {
												l3 = l3+inc;
											}else if (dir == DERIV_Y) {
												m3 = m3+inc;
											}else {
												n3 = n3+inc;
											}
										}else{
											if (dir == DERIV_X) {
												l4 = l4+inc;
											}else if (dir == DERIV_Y) {
												m4 = m4+inc;
											}else {
												n4 = n4+inc;
											}
										}
									}
									double r4 = ZERO;
									if (newL1>=1) {
										double tmp_r4 = eri(l1,m1,n1,ia,AA,
												l2,m2,n2,ja,BB,
												l3,m3,n3,ka,CC,
												l4,m4,n4,la,DD, 0);
										r4 = -2.0E0*beta*newL1*ic*jc*kc*lc*tmp_r4;
									}

									// now it's final result
									result += r1 + r2 + r3 + r4;
								}
							}
						}
					}

					// push it to the abcd
					abcd[count] = result;
					count++;
				}
			}
		}
	}
}

void eritest::eri_test(const string& derivFile, 
		const Int& maxL,const Int& auxMaxL,         const Int& order,           const Int& workType,
		const Int& inp, const vector<Double>& icoe, const vector<Double>& iexp, const Double* A, 
		const Int& jnp, const vector<Double>& jcoe, const vector<Double>& jexp, const Double* B, 
		const Int& knp, const vector<Double>& kcoe, const vector<Double>& kexp, const Double* C, 
		const Int& lnp, const vector<Double>& lcoe, const vector<Double>& lexp, const Double* D)
{
	//
	// form the data for hgp calculation
	// firstly it's the 
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
			Double diff = ia-ja; 
			Double ab   = -ia*ja/alpla;
			Double pref = exp(ab*AB2)*pow(PI/alpla,1.5E0);
			iexp2[count]= ONE/alpla;
			iexpdiff[count]= diff;
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

	// now it's ket side data
	Double CD2 = (C[0]-D[0])*(C[0]-D[0])+(C[1]-D[1])*(C[1]-D[1])+(C[2]-D[2])*(C[2]-D[2]);
	Int jnp2 = knp*lnp;
	vector<Double> jexp2(jnp2,ZERO);
	vector<Double> jexpdiff(jnp2,ZERO);
	vector<Double> fket(jnp2,ZERO);
	vector<Double> Q(3*jnp2,ZERO);
	count = 0;
	for(Int lp=0; lp<lnp; lp++) {
		for(Int kp=0; kp<knp; kp++) {

			// prefactors etc.
			Double ia   = kexp[kp];
			Double ja   = lexp[lp];
			Double alpla= ia+ja; 
			Double diff = ia-ja; 
			Double ab   = -ia*ja/alpla;
			Double pref = exp(ab*CD2)*pow(PI/alpla,1.5E0);
			jexp2[count]= ONE/alpla;
			jexpdiff[count]= diff;
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
	//
	Int len = 4000000;
	if (order == 2) len =100000000;
	LocalMemScr scr(len);

	//
	// real work begins
	// we note, that the work has to be in different branch
	// according to the ERI work type
	if (workType == TWO_BODY_ERI) {
		cout << "**************************************************************" << endl;
		cout << "two body electron repulsion integrals test:" << endl;
		cout << "**************************************************************" << endl;
		for(Int jpos=0; jpos<MAX_SHELL_TYPES; jpos++) {
			for(Int ipos=jpos; ipos<MAX_SHELL_TYPES; ipos++) {

				// now decode the L
				Int iLmin, iLmax;
				Int L1 = SHELL_ANG_MOM_CODE[ipos];
				decodeL(L1,iLmin,iLmax);
				Int kLmin, kLmax;
				Int L2 = SHELL_ANG_MOM_CODE[jpos];
				decodeL(L2,kLmin,kLmax);

				// now additional angular momentum
				Int jLmin = S; 
				Int jLmax = S;
				Int lLmin = S; 
				Int lLmax = S;

				// test with maxL
				if (iLmax>auxMaxL || kLmax>auxMaxL) continue;

				// now form the L code
				LInt LCode = codeSQ(L1,S,L2,S);

				// get the deriv information first
				DerivInforArray deriv(order,LCode);
				deriv.readInformation(derivFile);

				// normalize the coefficient array
				Int nLBra1 = iLmax-iLmin+1;
				Int nLBra2 = 1;
				Int nLKet1 = kLmax-kLmin+1;
				Int nLKet2 = 1;
				vector<Double> bra1Coe(nLBra1*inp);
				vector<Double> bra2Coe(jcoe);
				vector<Double> ket1Coe(nLKet1*knp);
				vector<Double> ket2Coe(lcoe);
				normCoe(iLmin, iLmax, iexp, icoe, bra1Coe);
				//normCoe(jLmin, jLmax, jexp, jcoe, bra2Coe);
				normCoe(kLmin, kLmax, kexp, kcoe, ket1Coe);
				//normCoe(lLmin, lLmax, lexp, lcoe, ket2Coe);

				// generate the shell pair for coe
				vector<Double> braCoePair(nLBra1*nLBra2*inp2);
				makeC2(iLmin,iLmax,jLmin,jLmax,inp,jnp,bra1Coe,bra2Coe,braCoePair);
				vector<Double> ketCoePair(nLKet1*nLKet2*jnp2);
				makeC2(kLmin,kLmax,lLmin,lLmax,knp,lnp,ket1Coe,ket2Coe,ketCoePair);

				// result array
				Int nBra1Bas = getCartBas(iLmin,iLmax);
				Int nKet1Bas = getCartBas(kLmin,kLmax);
				Int nDeriv   = deriv.getLenDerivInforArray();
				Int nTotalDeriv = deriv.getTotalNumDeriv();
				vector<Double> result(nBra1Bas*nKet1Bas*nTotalDeriv);
				UInt inp2_ = static_cast<UInt>(inp2);
				UInt jnp2_ = static_cast<UInt>(jnp2);
				Double pmax = 1.0E0;
				Double omega = 0.0E0;
				if (order == 1) {
#ifdef ORDER1
					hgp_os_eri_d1(LCode,inp2_,jnp2_,INTEGRAL_THRESH,pmax,omega,
							&braCoePair.front(),&iexp2.front(),&iexpdiff.front(),&fbra.front(),&P.front(),A,B, 
							&ketCoePair.front(),&jexp2.front(),&jexpdiff.front(),&fket.front(),&Q.front(),C,D, 
							&result.front(),scr);
#endif
				}else{
#ifdef ORDER2
					hgp_os_eri_d2(LCode,inp2_,jnp2_,INTEGRAL_THRESH,pmax,omega,
							&braCoePair.front(),&iexp2.front(),&iexpdiff.front(),&fbra.front(),&P.front(),A,B, 
							&ketCoePair.front(),&jexp2.front(),&jexpdiff.front(),&fket.front(),&Q.front(),C,D, 
							&result.front(),scr);
#endif
				}

				// clear the scr
				scr.reset();

				// now let's directly calculate the eri derivatives
				// loop over the number of derivatives
				Int derivIndex = 0;
				for(Int iDeriv=0; iDeriv<nDeriv; iDeriv++) {

					// get the deriv information
					const DerivInfor& derivInfor = deriv.getDerivInfor(iDeriv);
					Int pos1, pos2;
					derivInfor.getDerivPos(pos1,pos2);
					for(Int iDir=0; iDir<derivInfor.getDerivDirLen(); iDir++) {
						Int dir1, dir2;
						derivInfor.getDerivDirection(iDir,dir1,dir2);

						// get the pointer of result
						const Double* hgp = &result[nBra1Bas*nKet1Bas*derivIndex];

						Int nCarBas1 = getCartBas(iLmin,iLmax);
						Int nCarBas2 = getCartBas(jLmin,jLmax);
						Int nCarBas3 = getCartBas(kLmin,kLmax);
						Int d1       = nCarBas1;
						Int d2       = nCarBas1*nCarBas2;
						Int d3       = nCarBas1*nCarBas2*nCarBas3;
						for(Int Lk=kLmin; Lk<=kLmax; Lk++) {
							for(Int Li=iLmin; Li<=iLmax; Li++) {

								// additional two dimensions are S
								Int Lj = S;
								Int Ll = S;

								// initilize the results
								Int nBra1Bas = getCartBas(Li,Li);
								Int nBra2Bas = getCartBas(Lj,Lj);
								Int nKet1Bas = getCartBas(Lk,Lk);
								Int nKet2Bas = getCartBas(Ll,Ll);
								vector<Double> abcd(nBra1Bas*nBra2Bas*nKet1Bas*nKet2Bas);

								// coefficient array
								const Double* iC = &bra1Coe[(Li-iLmin)*inp];
								const Double* jC = &bra2Coe[(Lj-jLmin)*jnp];
								const Double* kC = &ket1Coe[(Lk-kLmin)*knp];
								const Double* lC = &ket2Coe[(Ll-lLmin)*lnp];

								// now do direct eri
								if (order == 1) {
									directeri_d1(Li,Lj,Lk,Ll,
											inp,iC,&iexp.front(),A,
											jnp,jC,&jexp.front(),B,
											knp,kC,&kexp.front(),C,
											lnp,lC,&lexp.front(),D,
											pos1,dir1,&abcd.front());
								}else if (order == 2) {
									if (pos1 == pos2 && dir1 == dir2) {
										directeri_d2(Li,Lj,Lk,Ll,
												inp,iC,&iexp.front(),A,
												jnp,jC,&jexp.front(),B,
												knp,kC,&kexp.front(),C,
												lnp,lC,&lexp.front(),D,
												pos1,dir1,&abcd.front());
									}else{
										directeri_d2(Li,Lj,Lk,Ll,
												inp,iC,&iexp.front(),A,
												jnp,jC,&jexp.front(),B,
												knp,kC,&kexp.front(),C,
												lnp,lC,&lexp.front(),D,
												pos1,pos2,dir1,dir2,&abcd.front());
									}
								}

								// now it's comparison
								for(Int l=0; l<nKet2Bas; l++) {
									for(Int k=0; k<nKet1Bas; k++) {
										for(Int j=0; j<nBra2Bas; j++) {
											for(Int i=0; i<nBra1Bas; i++) {

												// we need to get the index for the HGP 
												Int iBasIndex = getBasOffset(iLmin,Li,i);
												Int jBasIndex = getBasOffset(jLmin,Lj,j);
												Int kBasIndex = getBasOffset(kLmin,Lk,k);
												Int lBasIndex = getBasOffset(lLmin,Ll,l);
												Int index = iBasIndex+jBasIndex*d1+kBasIndex*d2+lBasIndex*d3;
												Double v1 = hgp[index];
												Double v2 = abcd[i+j*nBra1Bas+k*nBra1Bas*nBra2Bas+
													l*nBra1Bas*nBra2Bas*nKet1Bas];
												if (fabs(v1-v2)>THRESH) {
													cout << "derivative position  1: " << stringPos(pos1) << endl;
													cout << "derivative direction 1: " << stringDir(dir1) << endl;
													if (pos2>0) {
														cout << "derivative position  2: " << stringPos(pos2) << endl;
														cout << "derivative direction 2: " << stringDir(dir2) << endl;
													}
													cout << "Bra1's L: " << iLmin << " " << iLmax << endl;
													cout << "Bra2's L: " << S     << " " << S     << endl;
													cout << "Ket1's L: " << jLmin << " " << jLmax << endl;
													cout << "Ket2's L: " << S     << " " << S     << endl;
													cout << "result did not match: " << endl;
													printf("difference  : %-16.10f\n", fabs(v1-v2));
													printf("hgp    value: %-16.10f\n", v1);
													printf("direct value: %-16.10f\n", v2);
													Int l0,m0,n0,l1,m1,n1,l2,m2,n2,l3,m3,n3;
													getlmn(Li, Li, i, l0, m0, n0);
													getlmn(Lj, Lj, j, l1, m1, n1);
													getlmn(Lk, Lk, k, l2, m2, n2);
													getlmn(Ll, Ll, l, l3, m3, n3);
													cout << l0 << m0 << n0 << endl;
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

						// increment the deriv index
						derivIndex++;
					}
				}
			}
		}
	}else if (workType == THREE_BODY_ERI) {
		cout << "**************************************************************" << endl;
		cout << "three body electron repulsion integrals test:" << endl;
		cout << "**************************************************************" << endl;
		for(Int kpos=0; kpos<MAX_SHELL_TYPES; kpos++) {
			for(Int jpos=0; jpos<MAX_SHELL_TYPES; jpos++) {
				for(Int ipos=jpos; ipos<MAX_SHELL_TYPES; ipos++) {

					// now decode the L
					Int iLmin, iLmax;
					Int L1 = SHELL_ANG_MOM_CODE[ipos];
					decodeL(L1,iLmin,iLmax);
					Int jLmin, jLmax;
					Int L2 = SHELL_ANG_MOM_CODE[jpos];
					decodeL(L2,jLmin,jLmax);
					Int kLmin, kLmax;
					Int L3 = SHELL_ANG_MOM_CODE[kpos];
					decodeL(L3,kLmin,kLmax);
					Int lLmin = S; 
					Int lLmax = S;
					Int L4 = S;

					// test with max L
					if (iLmax>maxL || jLmax>maxL) continue;
					if (kLmax>auxMaxL) continue;

					// now form the L code
					LInt LCode = codeSQ(L1,L2,L3,L4);

					// get the deriv information first
					DerivInforArray deriv(order,LCode);
					deriv.readInformation(derivFile);

					// normalize the coefficient array
					Int nLBra1 = iLmax-iLmin+1;
					Int nLBra2 = jLmax-jLmin+1;
					Int nLKet1 = kLmax-kLmin+1;
					Int nLKet2 = lLmax-lLmin+1;
					vector<Double> bra1Coe(nLBra1*inp);
					vector<Double> bra2Coe(nLBra2*jnp);
					vector<Double> ket1Coe(nLKet1*knp);
					vector<Double> ket2Coe(lcoe);
					normCoe(iLmin, iLmax, iexp, icoe, bra1Coe);
					normCoe(jLmin, jLmax, jexp, jcoe, bra2Coe);
					normCoe(kLmin, kLmax, kexp, kcoe, ket1Coe);
					//normCoe(lLmin, lLmax, lexp, lcoe, ket2Coe);

					// generate the shell pair for coe
					vector<Double> braCoePair(nLBra1*nLBra2*inp2);
					makeC2(iLmin,iLmax,jLmin,jLmax,inp,jnp,bra1Coe,bra2Coe,braCoePair);
					vector<Double> ketCoePair(nLKet1*nLKet2*jnp2);
					makeC2(kLmin,kLmax,lLmin,lLmax,knp,lnp,ket1Coe,ket2Coe,ketCoePair);

					// result array
					Int nBra1Bas = getCartBas(iLmin,iLmax);
					Int nBra2Bas = getCartBas(jLmin,jLmax);
					Int nKet1Bas = getCartBas(kLmin,kLmax);
					Int nKet2Bas = getCartBas(lLmin,lLmax);
					Int nDeriv   = deriv.getLenDerivInforArray();
					Int nTotalDeriv = deriv.getTotalNumDeriv();
					vector<Double> result(nBra1Bas*nBra2Bas*nKet1Bas*nKet2Bas*nTotalDeriv);
					Double pmax = 1.0E0;
					Double omega = 0.0E0;
					if (order == 1) {
#ifdef ORDER1
						hgp_os_eri_d1(LCode,inp2,jnp2,INTEGRAL_THRESH,pmax,omega,
								&braCoePair.front(),&iexp2.front(),&iexpdiff.front(),&fbra.front(),&P.front(),A,B, 
								&ketCoePair.front(),&jexp2.front(),&jexpdiff.front(),&fket.front(),&Q.front(),C,D, 
								&result.front(),scr);
#endif
					}else{
#ifdef ORDER2
						hgp_os_eri_d2(LCode,inp2,jnp2,INTEGRAL_THRESH,pmax,omega,
								&braCoePair.front(),&iexp2.front(),&iexpdiff.front(),&fbra.front(),&P.front(),A,B, 
								&ketCoePair.front(),&jexp2.front(),&jexpdiff.front(),&fket.front(),&Q.front(),C,D, 
								&result.front(),scr);
#endif
					}

					// clear the scr
					scr.reset();

					// now let's directly calculate the eri derivatives
					// loop over the number of derivatives
					Int derivIndex = 0;
					for(Int iDeriv=0; iDeriv<nDeriv; iDeriv++) {

						// get the deriv information
						const DerivInfor& derivInfor = deriv.getDerivInfor(iDeriv);
						Int pos1, pos2;
						derivInfor.getDerivPos(pos1,pos2);
						for(Int iDir=0; iDir<derivInfor.getDerivDirLen(); iDir++) {
							Int dir1, dir2;
							derivInfor.getDerivDirection(iDir,dir1,dir2);

							// get the pointer of result
							const Double* hgp = &result[nBra1Bas*nBra2Bas*nKet1Bas*nKet2Bas*derivIndex];

							Int nCarBas1 = getCartBas(iLmin,iLmax);
							Int nCarBas2 = getCartBas(jLmin,jLmax);
							Int nCarBas3 = getCartBas(kLmin,kLmax);
							Int d1       = nCarBas1;
							Int d2       = nCarBas1*nCarBas2;
							Int d3       = nCarBas1*nCarBas2*nCarBas3;
							for(Int Ll=lLmin; Ll<=lLmax; Ll++) {
								for(Int Lk=kLmin; Lk<=kLmax; Lk++) {
									for(Int Lj=jLmin; Lj<=jLmax; Lj++) {
										for(Int Li=iLmin; Li<=iLmax; Li++) {

											// initilize the results
											Int nBra1Bas = getCartBas(Li,Li);
											Int nBra2Bas = getCartBas(Lj,Lj);
											Int nKet1Bas = getCartBas(Lk,Lk);
											Int nKet2Bas = getCartBas(Ll,Ll);
											vector<Double> abcd(nBra1Bas*nBra2Bas*nKet1Bas*nKet2Bas);

											// coefficient array
											const Double* iC = &bra1Coe[(Li-iLmin)*inp];
											const Double* jC = &bra2Coe[(Lj-jLmin)*jnp];
											const Double* kC = &ket1Coe[(Lk-kLmin)*knp];
											const Double* lC = &ket2Coe[(Ll-lLmin)*lnp];

											// now do direct eri
											if (order == 1) {
												directeri_d1(Li,Lj,Lk,Ll,
														inp,iC,&iexp.front(),A,
														jnp,jC,&jexp.front(),B,
														knp,kC,&kexp.front(),C,
														lnp,lC,&lexp.front(),D,
														pos1,dir1,&abcd.front());
											}else if (order == 2) {
												if (pos1 == pos2 && dir1 == dir2) {
													directeri_d2(Li,Lj,Lk,Ll,
															inp,iC,&iexp.front(),A,
															jnp,jC,&jexp.front(),B,
															knp,kC,&kexp.front(),C,
															lnp,lC,&lexp.front(),D,
															pos1,dir1,&abcd.front());
												}else{
													directeri_d2(Li,Lj,Lk,Ll,
															inp,iC,&iexp.front(),A,
															jnp,jC,&jexp.front(),B,
															knp,kC,&kexp.front(),C,
															lnp,lC,&lexp.front(),D,
															pos1,pos2,dir1,dir2,&abcd.front());
												}
											}

											// now it's comparison
											for(Int l=0; l<nKet2Bas; l++) {
												for(Int k=0; k<nKet1Bas; k++) {
													for(Int j=0; j<nBra2Bas; j++) {
														for(Int i=0; i<nBra1Bas; i++) {

															// we need to get the index for the HGP 
															Int iBasIndex = getBasOffset(iLmin,Li,i);
															Int jBasIndex = getBasOffset(jLmin,Lj,j);
															Int kBasIndex = getBasOffset(kLmin,Lk,k);
															Int lBasIndex = getBasOffset(lLmin,Ll,l);
															Int index = iBasIndex+jBasIndex*d1+kBasIndex*d2+lBasIndex*d3;
															Double v1 = hgp[index];
															Double v2 = abcd[i+j*nBra1Bas+k*nBra1Bas*nBra2Bas+
																l*nBra1Bas*nBra2Bas*nKet1Bas];
															if (fabs(v1-v2)>THRESH) {
																cout << "derivative position  1: " << stringPos(pos1) << endl;
																cout << "derivative direction 1: " << stringDir(dir1) << endl;
																if (pos2>0) {
																	cout << "derivative position  2: " << stringPos(pos2) << endl;
																	cout << "derivative direction 2: " << stringDir(dir2) << endl;
																}
																cout << "Bra1's L: " << iLmin << " " << iLmax << endl;
																cout << "Bra2's L: " << jLmin << " " << jLmax << endl;
																cout << "Ket1's L: " << kLmin << " " << kLmax << endl;
																cout << "Ket2's L: " << lLmin << " " << lLmax << endl;
																cout << "result did not match: " << endl;
																printf("difference  : %-16.10f\n", fabs(v1-v2));
																printf("hgp    value: %-16.10f\n", v1);
																printf("direct value: %-16.10f\n", v2);
																Int l0,m0,n0,l1,m1,n1,l2,m2,n2,l3,m3,n3;
																getlmn(Li, Li, i, l0, m0, n0);
																getlmn(Lj, Lj, j, l1, m1, n1);
																getlmn(Lk, Lk, k, l2, m2, n2);
																getlmn(Ll, Ll, l, l3, m3, n3);
																cout << l0 << m0 << n0 << endl;
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

							// increase the deriv index
							derivIndex++;
						}
					}
				}
			}
		}
	}else{
		cout << "**************************************************************" << endl;
		cout << "four body electron repulsion integrals test:" << endl;
		cout << "**************************************************************" << endl;
		for(Int n2=0; n2<MAX_SHELL_PAIR_NUMBER; n2++) {
			for(Int n1=n2; n1<MAX_SHELL_PAIR_NUMBER; n1++) {

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
				if (lLmax>maxL) continue;

				// next we have to form the LCode for single shell
				Int L1 = codeL(iLmin,iLmax);
				Int L2 = codeL(jLmin,jLmax);
				Int L3 = codeL(kLmin,kLmax);
				Int L4 = codeL(lLmin,lLmax);

				// now form the L code
				LInt LCode = codeSQ(L1,L2,L3,L4);

				// get the deriv information first
				DerivInforArray deriv(order,LCode);
				deriv.readInformation(derivFile);

				// normalize the coefficient array
				Int nLBra1 = iLmax-iLmin+1;
				Int nLBra2 = jLmax-jLmin+1;
				Int nLKet1 = kLmax-kLmin+1;
				Int nLKet2 = lLmax-lLmin+1;
				vector<Double> bra1Coe(nLBra1*inp);
				vector<Double> bra2Coe(nLBra2*jnp);
				vector<Double> ket1Coe(nLKet1*knp);
				vector<Double> ket2Coe(nLKet2*lnp);
				normCoe(iLmin, iLmax, iexp, icoe, bra1Coe);
				normCoe(jLmin, jLmax, jexp, jcoe, bra2Coe);
				normCoe(kLmin, kLmax, kexp, kcoe, ket1Coe);
				normCoe(lLmin, lLmax, lexp, lcoe, ket2Coe);

				// generate the shell pair for coe
				vector<Double> braCoePair(nLBra1*nLBra2*inp2);
				makeC2(iLmin,iLmax,jLmin,jLmax,inp,jnp,bra1Coe,bra2Coe,braCoePair);
				vector<Double> ketCoePair(nLKet1*nLKet2*jnp2);
				makeC2(kLmin,kLmax,lLmin,lLmax,knp,lnp,ket1Coe,ket2Coe,ketCoePair);

				// result array
				Int nBra1Bas = getCartBas(iLmin,iLmax);
				Int nBra2Bas = getCartBas(jLmin,jLmax);
				Int nKet1Bas = getCartBas(kLmin,kLmax);
				Int nKet2Bas = getCartBas(lLmin,lLmax);
				Int nDeriv   = deriv.getLenDerivInforArray();
				Int nTotalDeriv = deriv.getTotalNumDeriv();
				vector<Double> result(nBra1Bas*nBra2Bas*nKet1Bas*nKet2Bas*nTotalDeriv);

				// now call hgp here
				Double pmax = 1.0E0;
				Double omega = 0.0E0;
				if (order == 1) {
#ifdef ORDER1
					hgp_os_eri_d1(LCode,inp2,jnp2,INTEGRAL_THRESH,pmax,omega,
							&braCoePair.front(),&iexp2.front(),&iexpdiff.front(),&fbra.front(),&P.front(),A,B, 
							&ketCoePair.front(),&jexp2.front(),&jexpdiff.front(),&fket.front(),&Q.front(),C,D, 
							&result.front(),scr);
#endif
				}else{
#ifdef ORDER2
					hgp_os_eri_d2(LCode,inp2,jnp2,INTEGRAL_THRESH,pmax,omega,
							&braCoePair.front(),&iexp2.front(),&iexpdiff.front(),&fbra.front(),&P.front(),A,B, 
							&ketCoePair.front(),&jexp2.front(),&jexpdiff.front(),&fket.front(),&Q.front(),C,D, 
							&result.front(),scr);
#endif
				}

				// clear the scr
				scr.reset();

				// now let's directly calculate the eri derivatives
				// loop over the number of derivatives
				Int derivIndex = 0;
				for(Int iDeriv=0; iDeriv<nDeriv; iDeriv++) {

					// get the deriv information
					const DerivInfor& derivInfor = deriv.getDerivInfor(iDeriv);
					Int pos1, pos2;
					derivInfor.getDerivPos(pos1,pos2);
					for(Int iDir=0; iDir<derivInfor.getDerivDirLen(); iDir++) {
						Int dir1, dir2;
						derivInfor.getDerivDirection(iDir,dir1,dir2);

						// get the pointer of result
						const Double* hgp = &result[nBra1Bas*nBra2Bas*nKet1Bas*nKet2Bas*derivIndex];

						// now let's generate the eri deriv
						Int nCarBas1 = getCartBas(iLmin,iLmax);
						Int nCarBas2 = getCartBas(jLmin,jLmax);
						Int nCarBas3 = getCartBas(kLmin,kLmax);
						Int d1       = nCarBas1;
						Int d2       = nCarBas1*nCarBas2;
						Int d3       = nCarBas1*nCarBas2*nCarBas3;
						for(Int Ll=lLmin; Ll<=lLmax; Ll++) {
							for(Int Lk=kLmin; Lk<=kLmax; Lk++) {
								for(Int Lj=jLmin; Lj<=jLmax; Lj++) {
									for(Int Li=iLmin; Li<=iLmax; Li++) {

										// initilize the results
										Int nBra1Bas = getCartBas(Li,Li);
										Int nBra2Bas = getCartBas(Lj,Lj);
										Int nKet1Bas = getCartBas(Lk,Lk);
										Int nKet2Bas = getCartBas(Ll,Ll);
										vector<Double> abcd(nBra1Bas*nBra2Bas*nKet1Bas*nKet2Bas);

										// coefficient array
										const Double* iC = &bra1Coe[(Li-iLmin)*inp];
										const Double* jC = &bra2Coe[(Lj-jLmin)*jnp];
										const Double* kC = &ket1Coe[(Lk-kLmin)*knp];
										const Double* lC = &ket2Coe[(Ll-lLmin)*lnp];

										// now do direct eri
										if (order == 1) {
											directeri_d1(Li,Lj,Lk,Ll,
													inp,iC,&iexp.front(),A,
													jnp,jC,&jexp.front(),B,
													knp,kC,&kexp.front(),C,
													lnp,lC,&lexp.front(),D,
													pos1,dir1,&abcd.front());
										}else if (order == 2) {
											if (pos1 == pos2 && dir1 == dir2) {
												directeri_d2(Li,Lj,Lk,Ll,
														inp,iC,&iexp.front(),A,
														jnp,jC,&jexp.front(),B,
														knp,kC,&kexp.front(),C,
														lnp,lC,&lexp.front(),D,
														pos1,dir1,&abcd.front());
											}else{
												directeri_d2(Li,Lj,Lk,Ll,
														inp,iC,&iexp.front(),A,
														jnp,jC,&jexp.front(),B,
														knp,kC,&kexp.front(),C,
														lnp,lC,&lexp.front(),D,
														pos1,pos2,dir1,dir2,&abcd.front());
											}
										}

										// now it's comparison
										for(Int l=0; l<nKet2Bas; l++) {
											for(Int k=0; k<nKet1Bas; k++) {
												for(Int j=0; j<nBra2Bas; j++) {
													for(Int i=0; i<nBra1Bas; i++) {

														// we need to get the index for the HGP 
														Int iBasIndex = getBasOffset(iLmin,Li,i);
														Int jBasIndex = getBasOffset(jLmin,Lj,j);
														Int kBasIndex = getBasOffset(kLmin,Lk,k);
														Int lBasIndex = getBasOffset(lLmin,Ll,l);
														Int index = iBasIndex+jBasIndex*d1+kBasIndex*d2+lBasIndex*d3;
														Double v1 = hgp[index];
														Double v2 = abcd[i+j*nBra1Bas+k*nBra1Bas*nBra2Bas+
															l*nBra1Bas*nBra2Bas*nKet1Bas];
														if (fabs(v1-v2)>THRESH) {
															cout << "derivative position  1: " << stringPos(pos1) << endl;
															cout << "derivative direction 1: " << stringDir(dir1) << endl;
															if (pos2>0) {
																cout << "derivative position  2: " << stringPos(pos2) << endl;
																cout << "derivative direction 2: " << stringDir(dir2) << endl;
															}
															cout << "Bra1's L: " << iLmin << " " << iLmax << endl;
															cout << "Bra2's L: " << jLmin << " " << jLmax << endl;
															cout << "Ket1's L: " << kLmin << " " << kLmax << endl;
															cout << "Ket2's L: " << lLmin << " " << lLmax << endl;
															cout << "result did not match: " << endl;
															printf("difference  : %-16.10f\n", fabs(v1-v2));
															printf("hgp    value: %-16.10f\n", v1);
															printf("direct value: %-16.10f\n", v2);
															Int l0,m0,n0,l1,m1,n1,l2,m2,n2,l3,m3,n3;
															getlmn(Li, Li, i, l0, m0, n0);
															getlmn(Lj, Lj, j, l1, m1, n1);
															getlmn(Lk, Lk, k, l2, m2, n2);
															getlmn(Ll, Ll, l, l3, m3, n3);
															cout << l0 << m0 << n0 << endl;
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

						// increase the deriv index
						derivIndex++;
					}
				}
			}
		}
	}
}


