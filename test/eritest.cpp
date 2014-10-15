#include "libgen.h"
#include "angmomlist.h"
#include "shellprop.h"
#include "functions.h"
#include "norm.h"
#include "eri.h"
#include "localmemscr.h"
#include "eritest.h"
using namespace shellprop;
using namespace functions;
using namespace norm;
using namespace localmemscr;
using namespace eritest;

extern bool hgp_os_eri(const LInt& LCode, const UInt& inp2, const UInt& jnp2, 
		const Double& thresh, const Double* icoe, const Double* iexp, const Double* ifac, 
		const Double* P, const Double* A, const Double* B, const Double* jcoe, 
		const Double* jexp,const Double* jfac, const Double* Q, const Double* C, 
		const Double* D, Double* abcd, LocalMemScr& scr);

void eritest::directeri(const Int& Li, const Int& Lj, const Int& Lk, const Int& Ll,
		const Int& inp, const Double* icoe, const Double* iexp, const Double* A, 
		const Int& jnp, const Double* jcoe, const Double* jexp, const Double* B, 
		const Int& knp, const Double* kcoe, const Double* kexp, const Double* C, 
		const Int& lnp, const Double* lcoe, const Double* lexp, const Double* D, 
		Double* abcd)
{
	// number of basis sets
	Int nBra1Bas = getCartBas(Li,Li);
	Int nBra2Bas = getCartBas(Lj,Lj);
	Int nKet1Bas = getCartBas(Lk,Lk);
	Int nKet2Bas = getCartBas(Ll,Ll);

	// loop over each basis set
	Int count = 0;
	for(Int lBas=0; lBas<nKet2Bas; lBas++) {
		for(Int kBas=0; kBas<nKet1Bas; kBas++) {
			for(Int jBas=0; jBas<nBra2Bas; jBas++) {
				for(Int iBas=0; iBas<nBra1Bas; iBas++) {

					// get the angular momentum
					Int l1,m1,n1,l2,m2,n2;
					Int l3,m3,n3,l4,m4,n4;
					getlmn(Li, Li, iBas, l1, m1, n1);
					getlmn(Lj, Lj, jBas, l2, m2, n2);
					getlmn(Lk, Lk, kBas, l3, m3, n3);
					getlmn(Ll, Ll, lBas, l4, m4, n4);

					// form the primitives data
					Double result = ZERO;
					for(Int lp=0; lp<lnp; lp++) {
						for(Int kp=0; kp<knp; kp++) {
							for(Int jp=0; jp<jnp; jp++) {
								for(Int ip=0; ip<inp; ip++) {

#ifdef WITH_SINGLE_PRECISION
									// we do not have coefficient in eri calculation
									// therefore only fetch the exponentials
									double ia = iexp[ip];
									double ja = jexp[jp];
									double ka = kexp[kp];
									double la = lexp[lp];
									double ic = icoe[ip];
									double jc = jcoe[jp];
									double kc = kcoe[kp];
									double lc = lcoe[lp];

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

									// now do calculation here
									double r = eri(l1,m1,n1,ia,AA,
								         			l2,m2,n2,ja,BB,
											         l3,m3,n3,ka,CC,
											         l4,m4,n4,la,DD, 0);
									result += ic*jc*kc*lc*r;
#else
									// we do not have coefficient in eri calculation
									// therefore only fetch the exponentials
									Double ia = iexp[ip];
									Double ja = jexp[jp];
									Double ka = kexp[kp];
									Double la = lexp[lp];
									Double ic = icoe[ip];
									Double jc = jcoe[jp];
									Double kc = kcoe[kp];
									Double lc = lcoe[lp];

									// now do calculation here
									Double r = eri(l1,m1,n1,ia,A,
								         			l2,m2,n2,ja,B,
											         l3,m3,n3,ka,C,
											         l4,m4,n4,la,D, 0);
									result += ic*jc*kc*lc*r;
#endif
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

void eritest::eri_test(const Int& maxL, const Int& auxMaxL, const Int& workType,
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
	// we note, that the work has to be in different branch
	// according to the ERI work type
	// 
	cout << "**************************************************************" << endl;
	cout << "four body electron repulsion integrals test:" << endl;
	cout << "**************************************************************" << endl;
	if (workType == TWO_BODY_ERI) {
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
				vector<Double> result(nBra1Bas*nKet1Bas);
				UInt inp2_ = static_cast<UInt>(inp2);
				UInt jnp2_ = static_cast<UInt>(jnp2);
				hgp_os_eri(LCode,inp2_,jnp2_,INTEGRAL_THRESH,
						&braCoePair.front(),&iexp2.front(),&fbra.front(),&P.front(),A,B, 
						&ketCoePair.front(),&jexp2.front(),&fket.front(),&Q.front(),C,D, 
						&result.front(),scr);

				// clear the scr
				scr.reset();

				// now let's directly calculate the eri
				Int offset = 0;
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

						// now do direct ov
						directeri(Li,Lj,Lk,Ll,
								inp,iC,&iexp.front(),A,
								jnp,jC,&jexp.front(),B,
								knp,kC,&kexp.front(),C,
								lnp,lC,&lexp.front(),D,&abcd.front());

						// now it's comparison
						const Double* hgp = &result[offset];
						for(Int l=0; l<nKet2Bas; l++) {
							for(Int k=0; k<nKet1Bas; k++) {
								for(Int j=0; j<nBra2Bas; j++) {
									for(Int i=0; i<nBra1Bas; i++) {
										Double v1 = hgp[i+j*nBra1Bas+k*nBra1Bas*nBra2Bas+
											l*nBra1Bas*nBra2Bas*nKet1Bas];
										Double v2 = abcd[i+j*nBra1Bas+k*nBra1Bas*nBra2Bas+
											l*nBra1Bas*nBra2Bas*nKet1Bas];
										if (fabs(v1-v2)>THRESH) {
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

						// increase the offset
						offset += nBra1Bas*nBra2Bas*nKet1Bas*nKet2Bas;
					}
				}
			}
		}
	}else if (workType == THREE_BODY_ERI) {
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
					vector<Double> result(nBra1Bas*nBra2Bas*nKet1Bas*nKet2Bas);
					hgp_os_eri(LCode,inp2,jnp2,INTEGRAL_THRESH,
							&braCoePair.front(),&iexp2.front(),&fbra.front(),&P.front(),A,B, 
							&ketCoePair.front(),&jexp2.front(),&fket.front(),&Q.front(),C,D, 
							&result.front(),scr);

					// clear the scr
					scr.reset();

					// now let's directly calculate the eri
					Int offset = 0;
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

									// now do direct ov
									directeri(Li,Lj,Lk,Ll,
											inp,iC,&iexp.front(),A,
											jnp,jC,&jexp.front(),B,
											knp,kC,&kexp.front(),C,
											lnp,lC,&lexp.front(),D,&abcd.front());

									// now it's comparison
									const Double* hgp = &result[offset];
									for(Int l=0; l<nKet2Bas; l++) {
										for(Int k=0; k<nKet1Bas; k++) {
											for(Int j=0; j<nBra2Bas; j++) {
												for(Int i=0; i<nBra1Bas; i++) {
													Double v1 = hgp[i+j*nBra1Bas+k*nBra1Bas*nBra2Bas+
														l*nBra1Bas*nBra2Bas*nKet1Bas];
													Double v2 = abcd[i+j*nBra1Bas+k*nBra1Bas*nBra2Bas+
														l*nBra1Bas*nBra2Bas*nKet1Bas];
													if (fabs(v1-v2)>THRESH) {
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

									// increase the offset
									offset += nBra1Bas*nBra2Bas*nKet1Bas*nKet2Bas;
								}
							}
						}
					}
				}
			}
		}
	}else{
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
				vector<Double> result(nBra1Bas*nBra2Bas*nKet1Bas*nKet2Bas);

				// now call hgp here
				hgp_os_eri(LCode,inp2,jnp2,INTEGRAL_THRESH,
						&braCoePair.front(),&iexp2.front(),&fbra.front(),&P.front(),A,B, 
						&ketCoePair.front(),&jexp2.front(),&fket.front(),&Q.front(),C,D, 
						&result.front(),scr);

				// clear the scr
				scr.reset();

				// now let's directly calculate the eri
				Int offset = 0;
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

								// now do direct ov
								directeri(Li,Lj,Lk,Ll,
										inp,iC,&iexp.front(),A,
										jnp,jC,&jexp.front(),B,
										knp,kC,&kexp.front(),C,
										lnp,lC,&lexp.front(),D,&abcd.front());

								// now it's comparison
								const Double* hgp = &result[offset];
								for(Int l=0; l<nKet2Bas; l++) {
									for(Int k=0; k<nKet1Bas; k++) {
										for(Int j=0; j<nBra2Bas; j++) {
											for(Int i=0; i<nBra1Bas; i++) {
												Double v1 = hgp[i+j*nBra1Bas+k*nBra1Bas*nBra2Bas+
													l*nBra1Bas*nBra2Bas*nKet1Bas];
												Double v2 = abcd[i+j*nBra1Bas+k*nBra1Bas*nBra2Bas+
													l*nBra1Bas*nBra2Bas*nKet1Bas];
												if (fabs(v1-v2)>THRESH) {
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

								// increase the offset
								offset += nBra1Bas*nBra2Bas*nKet1Bas*nKet2Bas;
							}
						}
					}
				}
			}
		}
	}
}


