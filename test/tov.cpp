#include "libgen.h"
#include "angmomlist.h"
#include "shellprop.h"
#include "functions.h"
#include "norm.h"
#include "ov.h"
#include "localmemscr.h"
#include "tov.h"
using namespace shellprop;
using namespace functions;
using namespace norm;
using namespace ov;
using namespace localmemscr;
using namespace tov;

void hgp_os_threebodyoverlap(const LInt& LCode, const UInt& inp2, const UInt& jnp2, 
		const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, 
		const Double* A, const Double* B, const Double* jcoe, const Double* jexp, 
		const Double* C, Double* abcd, LocalMemScr& scr);

Double tov::threeOverlapIxyz(const Double& alpha, const Double& beta, const Double& gamma,
		const Int& na, const Int& nb, const Int& nc, const Double& G, const Double& A,
		const Double& B, const Double& C) 
{
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
				cout << "bionomialCoe  " << bionomialCoe(na,ka) << " "
					<< bionomialCoe(nb,kb) << " " << bionomialCoe(nc,kc) << endl;
				cout << "GA " << GA << endl;
				cout << "GB " << GB << endl;
				cout << "GC " << GC << endl;
				cout << "k " << k << endl;
				cout << "integral " << in << endl;
				Double r = bionomialCoe(na,ka)*bionomialCoe(nb,kb)*bionomialCoe(nc,kc)*GA*GB*GC*in;
				cout << "r " << r << endl;
				result += r;
			}
		}
	}
	return result*pidz;
}

Double tov::threeOverlap(const Double& alpha, const Double& beta, const Double& gamma,
		const Double* A, const Double* B, const Double* C,
		const Int& li, const Int& mi, const Int& ni,
		const Int& lj, const Int& mj, const Int& nj, 
		const Int& lk, const Int& mk, const Int& nk) 
{
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
	Double Ix = threeOverlapIxyz(alpha,beta,gamma,li,lj,lk,Gx,A[0],B[0],C[0]);
	cout << "Ix " << endl;
	Double Iy = threeOverlapIxyz(alpha,beta,gamma,mi,mj,mk,Gy,A[1],B[1],C[1]);
	cout << "Iy " << endl;
	Double Iz = threeOverlapIxyz(alpha,beta,gamma,ni,nj,nk,Gz,A[2],B[2],C[2]);
	cout << "Iz " << endl;
	//cout << "I " << kappa << " " << Ix << " " << Iy << " " << Iz << endl;
	return kappa*Ix*Iy*Iz;
}

Double tov::threeov(const Double& alpha, const Double& beta, const Double& gamma,
		const Double* A, const Double* B, const Double* C,
		const Int& li, const Int& mi, const Int& ni,
		const Int& lj, const Int& mj, const Int& nj, 
		const Int& lk, const Int& mk, const Int& nk) 
{
	// firstly calculate the prefactors etc.
	// this is for AB pair
	// P is the new center combined by A and B primitives
	Double zab = alpha + beta;
	Double Px  = (alpha*A[0] + beta*B[0])/zab;
	Double Py  = (alpha*A[1] + beta*B[1])/zab;
	Double Pz  = (alpha*A[2] + beta*B[2])/zab;
	Double PAx = Px - A[0];
	Double PAy = Py - A[1];
	Double PAz = Pz - A[2];
	Double PBx = Px - B[0];
	Double PBy = Py - B[1];
	Double PBz = Pz - B[2];
	Double AB2 = (A[0]-B[0])*(A[0]-B[0]) + (A[1]-B[1])*(A[1]-B[1]) + (A[2]-B[2])*(A[2]-B[2]);
	Double pre1= exp(-(alpha*beta/zab)*AB2);
	bool ABSameCenter = false;
	if (AB2<THRESHOLD_MATH) ABSameCenter = true;

	// now it's for PC pair
	// G is the new center combined by P and C primitives
	Double zeta= zab + gamma;
	Double Gx  = (zab*Px + gamma*C[0])/zeta;
	Double Gy  = (zab*Py + gamma*C[1])/zeta;
	Double Gz  = (zab*Pz + gamma*C[2])/zeta;
	Double GPx = Gx - Px;
	Double GPy = Gy - Py;
	Double GPz = Gz - Pz;
	Double GCx = Gx - C[0];
	Double GCy = Gy - C[1];
	Double GCz = Gz - C[2];
	Double PC2 = (Px-C[0])*(Px-C[0]) + (Py-C[1])*(Py-C[1]) + (Pz-C[2])*(Pz-C[2]);
	Double pre2= exp(-(zab*gamma/zeta)*PC2);
	bool PCSameCenter = false;
	if (PC2<THRESHOLD_MATH) PCSameCenter = true;

	if (ABSameCenter) {
		Double ov = overlapIntegral(pre2,zeta,GPx,GPy,GPz,GCx,GCy,GCz,  
				li+lj,mi+mj,ni+nj,lk,mk,nk,PCSameCenter);
		return ov*pre1;
	}else{
		Double sum = ZERO;
		for (Int i=0; i<=(Int)(li+lj); i++) {
			for (Int j=0; j<=(Int)(mi+mj); j++) {
				for (Int k=0; k<=(Int)(ni+nj); k++) {
					Double fk = Fk(i,li,lj,PAx,PBx)*Fk(j,mi,mj,PAy,PBy)*Fk(k,ni,nj,PAz,PBz);
					Double ov = overlapIntegral(pre2,zeta,GPx,GPy,GPz,GCx,GCy,GCz,  
							i,j,k,lk,mk,nk,PCSameCenter);
					sum += fk*ov; 
				}
			}
		}
		return sum*pre1;
	}
}

void tov::directtov(const Int& Li, const Int& Lj, const Int& Lk, 
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
				Double ov = ZERO;
				for(Int kp=0; kp<knp; kp++) {
					for(Int jp=0; jp<jnp; jp++) {
						for(Int ip=0; ip<inp; ip++) {
							Double ia = iexp[ip];
							Double ja = jexp[jp];
							Double ka = kexp[kp];
							Double ic = icoe[ip];
							Double jc = jcoe[jp];
							Double kc = kcoe[kp];
							Double r  = threeov(ia,ja,ka,A,B,C,  
									l1,m1,n1,l2,m2,n2,l3,m3,n3);
							ov += ic*jc*kc*r;
						}
					}
				}

				// push it to the abcd
				abcd[count] = ov;
				count++;
			}
		}
	}
}

void tov::tov_test(const Int& maxL, 
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

	//
	// set up the local mem scr
	// the length is set according to maxL = 5
	// auxMaxL = 5
	//
	LocalMemScr scr(92610);

	// now it's real work to do ov test
	cout << "**************************************************************" << endl;
	cout << " three body overlap integrals test:" << endl;
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
				hgp_os_threebodyoverlap(LCode,inp2_,knp_,&braCoePair.front(),&iexp2.front(), 
						&fbra.front(),&P.front(),A,B,&ket1Coe.front(),&jexp2.front(),
						C,&result.front(), scr);

				// reset the scr
				scr.reset();

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
										Double v1 = hgp[index];
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


