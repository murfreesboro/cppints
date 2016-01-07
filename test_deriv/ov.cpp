#include "libgen.h"
#include "angmomlist.h"
#include "shellprop.h"
#include "functions.h"
#include "derivinfor.h"
#include "norm.h"
#include "ov.h"
using namespace shellprop;
using namespace functions;
using namespace derivinfor;
using namespace norm;
using namespace ov;

#ifdef ORDER1
void hgp_os_twobodyoverlap_d1(const LInt& LCode, 
		const UInt& inp2, const Double* icoe, const Double* iexp, const Double* iexpdiff, 
		const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd);
#endif

#ifdef ORDER2
void hgp_os_twobodyoverlap_d2(const LInt& LCode, 
		const UInt& inp2, const Double* icoe, const Double* iexp, const Double* iexpdiff, 
		const Double* ifac, const Double* P, const Double* A, const Double* B, Double* abcd);
#endif

Double ov::Fk(const Int& k, const Int& l1, const Int& l2, 
			const Double& x1, const Double& x2) 
{
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
}

Double ov::overlapIntegral(
			const Double& d,   const Double& alpha, 
			const Double& PAx, const Double& PAy, const Double& PAz,  
			const Double& PBx, const Double& PBy, const Double& PBz,  
			const Int& li, const Int& mi, const Int& ni,
			const Int& lj, const Int& mj, const Int& nj, 
			bool isSameCenter) 
{

	if (isSameCenter) {
		return d*gammaFunInt(li+lj,alpha)*gammaFunInt(mi+mj,alpha)*gammaFunInt(ni+nj,alpha);
	}else{
		//
		// Ix
		//
		Double Ix  = ZERO;
		for (Int i=0; i<=(Int)(li+lj)/2; i++) {
			Double fk = Fk(2*i,li,lj,PAx,PBx);
			Double g  = gammaFunInt(2*i,alpha);
			Ix += fk*g; 
		}

		//
		// Iy
		//
		Double Iy  = ZERO;
		for (Int i=0; i<=(Int)(mi+mj)/2; i++) {
			Double fk= Fk(2*i,mi,mj,PAy,PBy);
			Double g = gammaFunInt(2*i,alpha);
			Iy += fk*g; 
		}

		//
		// Iz
		//
		Double Iz  = ZERO;
		for (Int i=0; i<=(Int)(ni+nj)/2; i++) {
			Double fk= Fk(2*i,ni,nj,PAz,PBz);
			Double g = gammaFunInt(2*i,alpha);
			Iz += fk*g; 
		}

		// now it's result
		return d*Ix*Iy*Iz;

	}
}

void ov::directov_d1(const Int& Li, const Int& Lj, 
		const Int& inp, const Double* icoe, const Double* iexp, const Double* A, 
		const Int& jnp, const Double* jcoe, const Double* jexp, const Double* B, 
		const Int& pos, const Int& dir, Double* abcd)
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
			Int ll1,mm1,nn1,ll2,mm2,nn2;
			getlmn(Li, Li, iBas, ll1, mm1, nn1);
			getlmn(Lj, Lj, jBas, ll2, mm2, nn2);

			// form the primitives data
			Double ov = ZERO;
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

					// here make a copy of original l,m,n values
					Int l1 = ll1;
					Int l2 = ll2;
					Int m1 = mm1;
					Int m2 = mm2;
					Int n1 = nn1;
					Int n2 = nn2;

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
					}
					Double r1 = ZERO;
					if (newL>=1) {
						Double r = overlapIntegral(pref,alpla, PAx, PAy, PAz, PBx, PBy, PBz,  
								l1, m1, n1, l2, m2, n2, isSameCenter);
						r1 = -1.0E0*newL*r;
					}

					// second term, do derivatives on x^ly^mz^n*exp(-alpha*r^2)
					// this will be, for example if on x;
					// 2*alpha*x^l+1y^mz^n*exp(-alpha*r^2)
					// copy the original l,m,n values first
					l1 = ll1;
					l2 = ll2;
					m1 = mm1;
					m2 = mm2;
					n1 = nn1;
					n2 = nn2;

					// we note that this term always exist
					Double alpha = ZERO;
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
					}
					Double tmp_r2 = overlapIntegral(pref,alpla, PAx, PAy, PAz, PBx, PBy, PBz,  
							l1, m1, n1, l2, m2, n2, isSameCenter);
					Double r2 = 2.0E0*alpha*tmp_r2;

					// now adding all things together
					ov += r1 + r2;
				}
			}

			// push it to the abcd
			abcd[count] = ov;
			count++;
		}
	}
}

void ov::directov_d2(const Int& Li, const Int& Lj, 
		const Int& inp, const Double* icoe, const Double* iexp, const Double* A, 
		const Int& jnp, const Double* jcoe, const Double* jexp, const Double* B, 
		const Int& pos, const Int& dir, Double* abcd)
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
			Int ll1,mm1,nn1,ll2,mm2,nn2;
			getlmn(Li, Li, iBas, ll1, mm1, nn1);
			getlmn(Lj, Lj, jBas, ll2, mm2, nn2);

			// form the primitives data
			Double ov = ZERO;
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

					// here make a copy of original l,m,n values
					Int l1 = ll1;
					Int l2 = ll2;
					Int m1 = mm1;
					Int m2 = mm2;
					Int n1 = nn1;
					Int n2 = nn2;

					// first term, do derivatives on x^ly^mz^n*exp(-alpha*r^2)
					// for example, if derivatives is on x; first on exp then
					// on the angular part, it gives
					// (-2*alpha*(l+1))*x^ly^mz^n*exp(-alpha*r^2)
					// this term always exist
					Double alpha = ZERO;
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
					}
					Double tmp_r1 = overlapIntegral(pref,alpla, PAx, PAy, PAz, PBx, PBy, PBz,  
							l1, m1, n1, l2, m2, n2, isSameCenter);
					Double r1 = -2.0E0*alpha*(newL+1)*tmp_r1;

					// second term, do derivatives on x^ly^mz^n*exp(-alpha*r^2)
					// for example, if derivatives is on x; firstly on the 
					// angular part then on the exp term, it gives
					// (-2*alpha*l)*x^ly^mz^n*exp(-alpha*r^2)
					// this term may not exist if l = 0
					Double r2 = -2.0E0*alpha*newL*tmp_r1;

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
					}
					Double r3 = 0.0E0;
					if (newL>1) {
						Double tmp_r3 = overlapIntegral(pref,alpla, PAx, PAy, PAz, PBx, PBy, PBz,  
								l1, m1, n1, l2, m2, n2, isSameCenter);
						r3 = newL*(newL-1)*tmp_r3;
					}

					// fourth term, do derivatives on x^ly^mz^n*exp(-alpha*r^2)
					// and all on exp term. this will be, for example if on x;
					// 4*alpha^2*x^l+2y^mz^n*exp(-alpha*r^2)
					// copy the original l,m,n values first
					l1 = ll1;
					l2 = ll2;
					m1 = mm1;
					m2 = mm2;
					n1 = nn1;
					n2 = nn2;

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
					}
					double tmp_r4 = overlapIntegral(pref,alpla, PAx, PAy, PAz, PBx, PBy, PBz,  
							l1, m1, n1, l2, m2, n2, isSameCenter);
					double r4 = 4.0E0*alpha*alpha*tmp_r4;

					// now adding all things together
					ov += r1 + r2 + r3 + r4;
				}
			}

			// push it to the abcd
			abcd[count] = ov;
			count++;
		}
	}
}

void ov::directov_d2(const Int& Li, const Int& Lj, 
		const Int& inp, const Double* icoe, const Double* iexp, const Double* A, 
		const Int& jnp, const Double* jcoe, const Double* jexp, const Double* B, 
		const Int& pos1, const Int& pos2, const Int& dir1, const Int& dir2, Double* abcd)
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
			Int ll1,mm1,nn1,ll2,mm2,nn2;
			getlmn(Li, Li, iBas, ll1, mm1, nn1);
			getlmn(Lj, Lj, jBas, ll2, mm2, nn2);

			// form the primitives data
			Double ov = ZERO;
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

					// here make a copy of original l,m,n values
					Int l1 = ll1;
					Int l2 = ll2;
					Int m1 = mm1;
					Int m2 = mm2;
					Int n1 = nn1;
					Int n2 = nn2;

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
						}

						// finally assign value
						if (i == 0) {
							alpha = val;
						}else{
							beta  = val;
						}
					}
					Double tmp_r1 = overlapIntegral(pref,alpla, PAx, PAy, PAz, PBx, PBy, PBz,  
							l1, m1, n1, l2, m2, n2, isSameCenter);
					Double r1 = 4.0E0*alpha*beta*tmp_r1;

					// here make a copy of original l,m,n values
					l1 = ll1;
					l2 = ll2;
					m1 = mm1;
					m2 = mm2;
					n1 = nn1;
					n2 = nn2;

					// second term: l*l'*I(l-1,l'-1)
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
						Double tmp_r2 = overlapIntegral(pref,alpla, PAx, PAy, PAz, PBx, PBy, PBz,  
								l1, m1, n1, l2, m2, n2, isSameCenter);
						r2 = newL1*newL2*tmp_r2;
					}

					// third term 
					// -2*alpha*l'*I(l+1,l'-1)
					// l' is just the newL2
					l1 = ll1;
					l2 = ll2;
					m1 = mm1;
					m2 = mm2;
					n1 = nn1;
					n2 = nn2;
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
						}
					}
					double r3 = ZERO;
					if (newL2>=1) {
						double tmp_r3 = overlapIntegral(pref,alpla, PAx, PAy, PAz, PBx, PBy, PBz,  
								l1, m1, n1, l2, m2, n2, isSameCenter);
						r3 = -2.0E0*alpha*newL2*tmp_r3;
					}

					// fourth term 
					// -2*beta*l*I(l-1,l'+1)
					// l is just the newL1
					l1 = ll1;
					l2 = ll2;
					m1 = mm1;
					m2 = mm2;
					n1 = nn1;
					n2 = nn2;
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
						}
					}
					double r4 = ZERO;
					if (newL1>=1) {
						Double tmp_r4 = overlapIntegral(pref,alpla, PAx, PAy, PAz, PBx, PBy, PBz,  
								l1, m1, n1, l2, m2, n2, isSameCenter);
						r4 = -2.0E0*beta*newL1*tmp_r4;
					}

					// now adding all things together
					ov += r1 + r2 + r3 + r4;
				}
			}

			// push it to the abcd
			abcd[count] = ov;
			count++;
		}
	}
}

void ov::ov_test(const string& derivFile, const Int& maxL, const Int& order,
		const Int& inp, const vector<Double>& icoe, const vector<Double>& iexp, const Double* A, 
		const Int& jnp, const vector<Double>& jcoe, const vector<Double>& jexp, const Double* B)
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

	// now it's real work to do ov test
	// we note that the bra1's L >= bra2's L
	//
	cout << "**************************************************************" << endl;
	cout << " two body overlap integrals test:" << endl;
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

			// get the deriv information first
			DerivInforArray deriv(order,LCode);
			deriv.readInformation(derivFile);

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
			Int nDeriv   = deriv.getLenDerivInforArray();
			Int nTotalDeriv = deriv.getTotalNumDeriv();
			vector<Double> result(nBra1Bas*nBra2Bas*nTotalDeriv);
			UInt np2_    = static_cast<UInt>(np2);
			if (order == 1) {
#ifdef ORDER1
			hgp_os_twobodyoverlap_d1(LCode,np2_,&braCoePair.front(),&iexp2.front(), 
					&iexpdiff.front(),&fbra.front(),&P.front(),A,B,&result.front());
#endif
			}else{
#ifdef ORDER2
			hgp_os_twobodyoverlap_d2(LCode,np2_,&braCoePair.front(),&iexp2.front(), 
					&iexpdiff.front(),&fbra.front(),&P.front(),A,B,&result.front());
#endif
			}

			// now let's directly calculate the ov derivatives
			// loop over the number of derivatives
			Int derivIndex = 0;
			for(Int iDeriv=0; iDeriv<nDeriv; iDeriv++) {

				// get the deriv information
				const DerivInfor& derivInfor = deriv.getDerivInfor(iDeriv);
				Int pos1, pos2;
				derivInfor.getDerivPos(pos1,pos2);
				for(Int iDir=0; iDir<derivInfor.getDerivDirLen(); iDir++) {
					Int dir1,dir2;
					derivInfor.getDerivDirection(iDir,dir1,dir2);

					// get the pointer of result
					const Double* hgp = &result[nBra1Bas*nBra2Bas*derivIndex];

					// now direct calculation
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
							if (order == 1) {
								directov_d1(Li,Lj,inp,iCArray,&iexp.front(),A,jnp,jCArray, 
										&jexp.front(),B,pos1,dir1,&abcd.front());
							}else if (order == 2) {
								if (pos1 == pos2 && dir1 == dir2) {
									directov_d2(Li,Lj,inp,iCArray,&iexp.front(),A,jnp,jCArray, 
											&jexp.front(),B,pos1,dir1,&abcd.front());
								}else{
									directov_d2(Li,Lj,inp,iCArray,&iexp.front(),A,jnp,jCArray, 
											&jexp.front(),B,pos1,pos2,dir1,dir2,&abcd.front());
								}
							}

							// now it's comparison
							for(Int jBas=0; jBas<nBra2Bas; jBas++) {
								for(Int iBas=0; iBas<nBra1Bas; iBas++) {
									Int iBasIndex = getBasOffset(iLmin,Li,iBas);
									Int jBasIndex = getBasOffset(jLmin,Lj,jBas);
									Int index = iBasIndex+jBasIndex*d1;
									Double v1 = hgp[index];
									Double v2 = abcd[iBas+jBas*nBra1Bas];
									if (fabs(v1-v2)>THRESH) {
										cout << "derivative position  1 " << stringPos(pos1) << endl;
										cout << "derivative direction 1 " << stringDir(dir1) << endl;
										if (pos2>0) {
											cout << "derivative position  2: " << stringPos(pos2) << endl;
											cout << "derivative direction 2: " << stringDir(dir2) << endl;
										}
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

					// increment the deriv index
					derivIndex++;
				}
			}
		}
	}
}


