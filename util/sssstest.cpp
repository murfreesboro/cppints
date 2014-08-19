#include<cstddef>
#include "constants.h"
// common C head files used in the program
#include<cstdlib>
#include<cstdio>

// external math functions 
#include<cmath>
#include <boost/math/special_functions/gamma.hpp>
#include<string> 
#include<vector>
#include<iostream>
#include<fstream>
#include<iomanip>
#include<string> 
#include<list>
#include<set>
#include<map>
#include<iterator>
#include<algorithm>
#define THRESHOLD_MATH  0.000000000000001E0
using namespace std;
using namespace boost::math;
typedef size_t        UInt;
typedef double        Double;

Double fm(const Double& a, const UInt& m); 

void newSSSS(const UInt& maxM, const UInt& inp2, const Double* iexp, const Double* icoe, 
		const Double* ifac, 
		const UInt& jnp2, const Double* jexp, const Double* jcoe, const Double* jfac, 
		const Double* A, const Double* B, const Double* P,
		const Double* C, const Double* D, const Double* Q, Double* ssss);

void oldSSSS(const UInt& maxM, const UInt& iNP2, const Double* iCoe, const Double* iExp, const Double* A, const Double* B, const UInt& jNP2, const Double* jCoe, const Double* jExp, const Double* C, const Double* D);

Double fm(const Double& a, const UInt& m) {
	Double x = 1.0E0;
	if (fabs(a)<1.0E-14) return (pow(x,2*m+1)/(2*m+1));
	Double f12 = 1.0E0/2.0E0;
	Double p   = 1.0E0/(2.0E0*pow(a,m+f12));
	Double upperLimit = a*pow(x,2.0E0);
	return tgamma_lower(m+f12,upperLimit)*p;
}

void newSSSS(const UInt& maxM, const UInt& inp2, const Double* iexp, const Double* icoe, 
		const Double* ifac, 
		const UInt& jnp2, const Double* jexp, const Double* jcoe, const Double* jfac, 
		const Double* A, const Double* B, const Double* P,
		const Double* C, const Double* D, const Double* Q) 
{

	vector<Double> result(maxM+1,0.0E0);
	for(UInt ip2=0; ip2<inp2; ip2++) {
		Double onedz = iexp[ip2];
		Double ic2   = icoe[ip2];
		Double fbra  = ifac[ip2];
		Double oned2z= 0.5E0*onedz;
		UInt offsetP = 3*ip2;
		Double PX    = P[offsetP  ];
		Double PY    = P[offsetP+1];
		Double PZ    = P[offsetP+2];
		Double PAX   = PX - A[0];
		Double PAY   = PY - A[1];
		Double PAZ   = PZ - A[2];
		Double PBX   = PX - B[0];
		Double PBY   = PY - B[1];
		Double PBZ   = PZ - B[2];
		for(UInt jp2=0; jp2<jnp2; jp2++) {
			Double onede = jexp[jp2];
			Double jc2   = jcoe[jp2];
			Double fket  = jfac[jp2];
			Double oned2e= 0.5E0*onedz;
			UInt offsetQ= 3*jp2;
			Double QX    = Q[offsetQ  ];
			Double QY    = Q[offsetQ+1];
			Double QZ    = Q[offsetQ+2];
			Double QCX   = QX - C[0];
			Double QCY   = QY - C[1];
			Double QCZ   = QZ - C[2];
			Double QDX   = QX - D[0];
			Double QDY   = QY - D[1];
			Double QDZ   = QZ - D[2];
			Double rho   = 1.0E0/(onedz+onede);
			Double sqrho = sqrt(rho);
			Double WX    = rho*(PX*onede + QX*onedz);
			Double WY    = rho*(PY*onede + QY*onedz);
			Double WZ    = rho*(PZ*onede + QZ*onedz);
			Double oned2k= 0.5E0*rho*onede*onedz;
			Double WPX   = WX - PX;
			Double WPY   = WY - PY;
			Double WPZ   = WZ - PZ;
			Double rhod2zsq = rho*oned2z*onedz;
			Double WQX   = WX - QX;
			Double WQY   = WY - QY;
			Double WQZ   = WZ - QZ;
			Double rhod2esq = rho*oned2e*onede;
			Double prefactor = ic2*jc2*fbra*fket;
			Double PQ2   = (PX-QX)*(PX-QX)+(PY-QY)*(PY-QY)+(PZ-QZ)*(PZ-QZ);
			Double u     = rho*PQ2;
			Double squ   = sqrt(u);

			// here we need to set the ssss integrals
			Double I_SSSS_M0  = 0.0E0;
			Double I_SSSS_M1  = 0.0E0;
			Double I_SSSS_M2  = 0.0E0;
			Double I_SSSS_M3  = 0.0E0;
			Double I_SSSS_M4  = 0.0E0;
			Double I_SSSS_M5  = 0.0E0;
			Double I_SSSS_M6  = 0.0E0;
			Double I_SSSS_M7  = 0.0E0;
			Double I_SSSS_M8  = 0.0E0;
			Double I_SSSS_M9  = 0.0E0;
			Double I_SSSS_M10 = 0.0E0;
			Double I_SSSS_M11 = 0.0E0;

			//
			// we use this way to generate the (SS|SS)^{m}
			// maxM=0, we just use erf function
			// if maxM<=10 and maxM>0, then mainly we use up recursive relation
			// else we will calculate fm directly and use down recursive relation
			// see the doc for more information
			//
			if (maxM==0) {
				I_SSSS_M0 = (prefactor*sqrho/squ)*erf(squ);
				//cout << "new ssss prefactor " << (prefactor*sqrho)*TWOOVERSQRTPI << endl;
				result[0] += I_SSSS_M0;
			}else if (maxM<=10) {

				//
				// now here for maxM<=10 and maxM>0
				// we also have two situations
				// 1  if u <=1; use power series to get f_{max}(u), then go down recursive way;
				// 2  if u > 1; calculate erf(u), then go up recursive way;
				//

				if (u<=1.0E0) {

					// calculate (SS|SS)^{Mmax}
					Double u2  = 2.0E0*u;
					I_SSSS_M8  = 1.0E0+u2*ONEOVER41; 
					I_SSSS_M8  = 1.0E0+u2*ONEOVER39*I_SSSS_M8; 
					I_SSSS_M8  = 1.0E0+u2*ONEOVER37*I_SSSS_M8; 
					I_SSSS_M8  = 1.0E0+u2*ONEOVER35*I_SSSS_M8; 
					I_SSSS_M8  = 1.0E0+u2*ONEOVER33*I_SSSS_M8; 
					I_SSSS_M8  = 1.0E0+u2*ONEOVER31*I_SSSS_M8; 
					I_SSSS_M8  = 1.0E0+u2*ONEOVER29*I_SSSS_M8; 
					I_SSSS_M8  = 1.0E0+u2*ONEOVER27*I_SSSS_M8; 
					I_SSSS_M8  = 1.0E0+u2*ONEOVER25*I_SSSS_M8; 
					I_SSSS_M8  = 1.0E0+u2*ONEOVER23*I_SSSS_M8; 
					I_SSSS_M8  = 1.0E0+u2*ONEOVER21*I_SSSS_M8; 
					I_SSSS_M8  = 1.0E0+u2*ONEOVER19*I_SSSS_M8; 
					I_SSSS_M8  = ONEOVER17*I_SSSS_M8;
					Double eu  = exp(-u);
					Double f   = TWOOVERSQRTPI*prefactor*sqrho*eu;
					I_SSSS_M8  = f*I_SSSS_M8;

					// now use down recursive relation to get 
					// rest of (SS|SS)^{m}
					I_SSSS_M7  = ONEOVER15*(u2*I_SSSS_M8+f);
					I_SSSS_M6  = ONEOVER13*(u2*I_SSSS_M7+f);
					I_SSSS_M5  = ONEOVER11*(u2*I_SSSS_M6+f);
					I_SSSS_M4  = ONEOVER9*(u2*I_SSSS_M5+f);
					I_SSSS_M3  = ONEOVER7*(u2*I_SSSS_M4+f);
					I_SSSS_M2  = ONEOVER5*(u2*I_SSSS_M3+f);
					I_SSSS_M1  = ONEOVER3*(u2*I_SSSS_M2+f);
					I_SSSS_M0  = ONEOVER1*(u2*I_SSSS_M1+f);

					result[0] += I_SSSS_M0;
					result[1] += I_SSSS_M1;
					result[2] += I_SSSS_M2;
					result[3] += I_SSSS_M3;
					result[4] += I_SSSS_M4;
					result[5] += I_SSSS_M5;
					result[6] += I_SSSS_M6;
					result[7] += I_SSSS_M7;
					result[8] += I_SSSS_M8;
				} else {

					//cout << "up recursive relation " << endl;
					// calculate (SS|SS)^{0}
					I_SSSS_M0     = (prefactor*sqrho/squ)*erf(squ);

					// now use up recursive relation
					Double oneO2u = 0.5E0/u;
					Double eu     = exp(-u);
					Double f      = TWOOVERSQRTPI*prefactor*sqrho*eu;
					I_SSSS_M1     = oneO2u*(1.0E0*I_SSSS_M0-f);
					I_SSSS_M2     = oneO2u*(3.0E0*I_SSSS_M1-f);
					I_SSSS_M3     = oneO2u*(5.0E0*I_SSSS_M2-f);
					I_SSSS_M4     = oneO2u*(7.0E0*I_SSSS_M3-f);
					I_SSSS_M5     = oneO2u*(9.0E0*I_SSSS_M4-f);
					I_SSSS_M6     = oneO2u*(11.0E0*I_SSSS_M5-f);
					I_SSSS_M7     = oneO2u*(13.0E0*I_SSSS_M6-f);
					I_SSSS_M8     = oneO2u*(15.0E0*I_SSSS_M7-f);

					result[0] += I_SSSS_M0;
					result[1] += I_SSSS_M1;
					result[2] += I_SSSS_M2;
					result[3] += I_SSSS_M3;
					result[4] += I_SSSS_M4;
					result[5] += I_SSSS_M5;
					result[6] += I_SSSS_M6;
					result[7] += I_SSSS_M7;
					result[8] += I_SSSS_M8;

					//
					// The algorithm below is actually not numerically
					// good. and I think erf function may not more 
					// expensive than calculating O*(2m-1)!!/|PQ|
					// therefore in our codes we do not do it
					//
					/*

					//cout << "far end expansion" << endl;
					// calculate (SS|SS)^{Mmax} through O*(2m-1)!!/|PQ|
					Double f0     = prefactor*sqrho;
					I_SSSS_M8     = f0/squ;
					Double oneO2u = 0.5E0/u;
					I_SSSS_M8    *= oneO2u; 
					I_SSSS_M8    *= 3.0E0*oneO2u; 
					I_SSSS_M8    *= 5.0E0*oneO2u; 
					I_SSSS_M8    *= 7.0E0*oneO2u; 
					I_SSSS_M8    *= 9.0E0*oneO2u; 
					I_SSSS_M8    *= 11.0E0*oneO2u; 
					I_SSSS_M8    *= 13.0E0*oneO2u; 
					I_SSSS_M8    *= 15.0E0*oneO2u; 

					//Double x = TWOOVERSQRTPI*prefactor*sqrho*fm(u,8);
					//cout << "difference between his I_SSSS_M8 and standard" << 
					//	fabs(I_SSSS_M8-x) << endl;

					// now use down recursive relation to get 
					// rest of (SS|SS)^{m}
					Double u2     = 2.0E0*u;
					Double f      = TWOOVERSQRTPI*f0*exp(-u);
					I_SSSS_M7     = ONEOVER15*(u2*I_SSSS_M8+f);
					I_SSSS_M6     = ONEOVER13*(u2*I_SSSS_M7+f);
					I_SSSS_M5     = ONEOVER11*(u2*I_SSSS_M6+f);
					I_SSSS_M4     = ONEOVER9*(u2*I_SSSS_M5+f);
					I_SSSS_M3     = ONEOVER7*(u2*I_SSSS_M4+f);
					I_SSSS_M2     = ONEOVER5*(u2*I_SSSS_M3+f);
					I_SSSS_M1     = ONEOVER3*(u2*I_SSSS_M2+f);
					I_SSSS_M0     = u2*I_SSSS_M1+f;

					cout << "u is " << u << endl;
					Double x = TWOOVERSQRTPI*prefactor*sqrho*fm(u,0);
					cout << "difference between his I_SSSS_M0 and standard in down recursive " << 
					fabs(I_SSSS_M0-x) << endl;
					x = TWOOVERSQRTPI*prefactor*sqrho*fm(u,8);
					cout << "difference between his I_SSSS_M8 and standard in down recursive " << 
					fabs(I_SSSS_M8-x) << endl;
					x = TWOOVERSQRTPI*prefactor*sqrho*fm(u,7);
					cout << "difference between his I_SSSS_M7 and standard in down recursive " << 
					fabs(I_SSSS_M7-x) << endl;
					x = TWOOVERSQRTPI*prefactor*sqrho*fm(u,6);
					cout << "difference between his I_SSSS_M6 and standard in down recursive " << 
					fabs(I_SSSS_M6-x) << endl;
					x = TWOOVERSQRTPI*prefactor*sqrho*fm(u,5);
					cout << "difference between his I_SSSS_M5 and standard in down recursive " << 
					fabs(I_SSSS_M5-x) << endl;
					x = TWOOVERSQRTPI*prefactor*sqrho*fm(u,4);
					cout << "difference between his I_SSSS_M4 and standard in down recursive " << 
					fabs(I_SSSS_M4-x) << endl;
					x = TWOOVERSQRTPI*prefactor*sqrho*fm(u,3);
					cout << "difference between his I_SSSS_M3 and standard in down recursive " << 
					fabs(I_SSSS_M3-x) << endl;
					x = TWOOVERSQRTPI*prefactor*sqrho*fm(u,2);
					cout << "difference between his I_SSSS_M2 and standard in down recursive " << 
					fabs(I_SSSS_M2-x) << endl;
					x = TWOOVERSQRTPI*prefactor*sqrho*fm(u,1);
					cout << "difference between his I_SSSS_M1 and standard in down recursive " << 
					fabs(I_SSSS_M1-x) << endl;

					result[0] += I_SSSS_M0;
					result[1] += I_SSSS_M1;
					result[2] += I_SSSS_M2;
					result[3] += I_SSSS_M3;
					result[4] += I_SSSS_M4;
					result[5] += I_SSSS_M5;
					result[6] += I_SSSS_M6;
					result[7] += I_SSSS_M7;
					result[8] += I_SSSS_M8;
					*/
				}

			}else{

				// calculate (SS|SS)^{Mmax} with incomplete gamma function
				Double f0     = TWOOVERSQRTPI*prefactor*sqrho;
				I_SSSS_M11    = f0*fm(u,11);

				// now use down recursive relation to get 
				// rest of (SS|SS)^{m}
				Double u2     = 2.0E0*u;
				Double f      = f0*exp(-u);
				I_SSSS_M10    = ONEOVER21*(u2*I_SSSS_M11+f);
				I_SSSS_M9     = ONEOVER19*(u2*I_SSSS_M10+f);
				I_SSSS_M8     = ONEOVER17*(u2*I_SSSS_M9+f);
				I_SSSS_M7     = ONEOVER15*(u2*I_SSSS_M8+f);
				I_SSSS_M6     = ONEOVER13*(u2*I_SSSS_M7+f);
				I_SSSS_M5     = ONEOVER11*(u2*I_SSSS_M6+f);
				I_SSSS_M4     = ONEOVER9*(u2*I_SSSS_M5+f);
				I_SSSS_M3     = ONEOVER7*(u2*I_SSSS_M4+f);
				I_SSSS_M2     = ONEOVER5*(u2*I_SSSS_M3+f);
				I_SSSS_M1     = ONEOVER3*(u2*I_SSSS_M2+f);
				I_SSSS_M0     = ONEOVER1*(u2*I_SSSS_M1+f);

				result[0] += I_SSSS_M0;
				result[1] += I_SSSS_M1;
				result[2] += I_SSSS_M2;
				result[3] += I_SSSS_M3;
				result[4] += I_SSSS_M4;
				result[5] += I_SSSS_M5;
				result[6] += I_SSSS_M6;
				result[7] += I_SSSS_M7;
				result[8] += I_SSSS_M8;
				result[9] += I_SSSS_M9;
				result[10] += I_SSSS_M10;
				result[11] += I_SSSS_M11;
			}
		}
	}
	for(UInt n=0; n<=maxM; n++) {
		printf("M is: %d  integral: %-18.14f\n", n, result[n]);
	}
}

void oldSSSS(const UInt& maxM, const UInt& iNP2, const Double* iCoe, const Double* iExp, const Double* A, const Double* B, const UInt& jNP2, const Double* jCoe, const Double* jExp, const Double* C, const Double* D) 
{
	vector<Double> result(maxM+1,0.0E0);
	for(UInt ip=0; ip<iNP2; ip++) {

		Double alpha = iExp[2*ip];
		Double beta  = iExp[2*ip+1];
		Double zeta  = alpha+beta;
		Double onedz = 1.0E0/zeta;
		Double ic    = iCoe[ip];
		Double adz   = alpha*onedz;
		Double bdz   = beta*onedz;
		Double PX    = A[0]*adz + B[0]*bdz;
		Double PY    = A[1]*adz + B[1]*bdz;
		Double PZ    = A[2]*adz + B[2]*bdz;
		Double PAX   = PX - A[0];
		Double PAY   = PY - A[1];
		Double PAZ   = PZ - A[2];
		Double pidz = PI*onedz;
		Double I_S_S_bra = ic*pow(pidz,1.5E0);
		Double oned2z = 0.5E0*onedz;
		Double oned2zsq  = oned2z*onedz;
		for(UInt jp=0; jp<jNP2; jp++) {

			Double gamma = jExp[2*jp];
			Double delta = jExp[2*jp+1];
			Double eta   = gamma+delta;
			Double onede = 1.0E0/eta;
			Double jc    = jCoe[jp];
			Double gde   = gamma*onede;
			Double dde   = delta*onede;
			Double QX    = C[0]*gde + D[0]*dde;
			Double QY    = C[1]*gde + D[1]*dde;
			Double QZ    = C[2]*gde + D[2]*dde;
			Double QCX   = QX - C[0];
			Double QCY   = QY - C[1];
			Double QCZ   = QZ - C[2];
			Double PQX   = PX - QX;
			Double PQY   = PY - QY;
			Double PQZ   = PZ - QZ;
			Double PQ2   = PQX*PQX+PQY*PQY+PQZ*PQZ;
			Double kappa = zeta+eta;
			Double onedk = 1.0E0/kappa;
			Double rho   = zeta*eta*onedk;
			Double zdk   = zeta*onedk;
			Double edk   = eta*onedk;
			Double WX    = PX*zdk + QX*edk;
			Double WY    = PY*zdk + QY*edk;
			Double WZ    = PZ*zdk + QZ*edk;
			Double oned2k= 0.5E0*onedk;
			Double WPX   = WX - PX;
			Double WPY   = WY - PY;
			Double WPZ   = WZ - PZ;
			Double rhod2zsq = rho*oned2zsq;
			Double WQX   = WX - QX;
			Double WQY   = WY - QY;
			Double WQZ   = WZ - QZ;
			Double oned2e   = 0.5E0*onede;
			Double oned2esq = oned2e*onede;
			Double rhod2esq = rho*oned2esq;
			Double pide  = PI*onede;
			Double rhodpi= rho/PI;
			Double factor= jc*pow(pide,1.5E0)*pow(rhodpi,0.5E0);
			Double u     = rho*PQ2;
			Double prefactor = 2.0E0*I_S_S_bra*factor;
			//cout << "old ssss prefactor " << prefactor << endl;

			for(UInt n=0; n<=maxM; n++) {
				result[n] += prefactor*fm(u,n);
			}
		}
	}
	for(UInt n=0; n<=maxM; n++) {
		printf("M is: %d  integral: %-18.14f\n", n, result[n]);
	}
}


int main() {

	//
	// lambda is used to adjust the distance betwen P and Q
	//
	Double lambda = 2.0;
	UInt ntest = 3;


	//
	// form the basic shell data - arbitrary
	//
	Double exp1[3];
	Double c1[3];
	exp1[0] = 18.7311370;
	exp1[1] = 2.8253937;
	exp1[2] = 0.6401217;
	c1[0]   = 0.03349460;
	c1[1]   = 0.23472695;
	c1[2]   = 0.81375733;
	Double A[3];
	A[0]    = 0.0;
	A[1]    = 0.0;
	A[2]    = 0.0;

	Double exp2[3];
	Double c2[3];
	exp2[0] = 7.8682724;
	exp2[1] = 1.8812885;
	exp2[2] = 0.5442493;
	c2[0]   = 0.0689991;
	c2[1]   = 0.3164240;
	c2[2]   = 0.7443083;
	Double B[3];
	B[0]    = 1.0;
	B[1]    = 0.0;
	B[2]    = 0.0;

	Double exp3[3];
	Double c3[3];
	exp3[0] = 5.8682724;
	exp3[1] = 2.8812885;
	exp3[2] = 0.7442493;
	c3[0]   = 0.0149991;
	c3[1]   = 0.6164240;
	c3[2]   = 0.5443083;
	Double C[3];
	C[0]    = 0.0;
	C[1]    = 1.0*lambda;
	C[2]    = 0.0;

	Double exp4[3];
	Double c4[3];
	exp4[0] = 8.78682724;
	exp4[1] = 4.11112885;
	exp4[2] = 0.2442493;
	c4[0]   = 0.0449991;
	c4[1]   = 0.264240;
	c4[2]   = 0.2443083;
	Double D[3];
	D[0]    = 1.0;
	D[1]    = 1.0*lambda;
	D[2]    = 0.0;


	//
	// now do the SSSS in old way
	//
	bool doOldSSSS = true;
	if (doOldSSSS) {

		// AB pairs
		UInt iNP2 = 9;
		vector<Double> iExp;
		vector<Double> iCoe;
		iExp.reserve(2*iNP2);
		iCoe.reserve(iNP2);
		Double AB2 = pow(A[0]-B[0],2)+pow(A[1]-B[1],2)+pow(A[2]-B[2],2);
		for(UInt j=0;j<3;j++) {
			for(UInt i=0;i<3;i++) {
				Double ia    = exp1[i];
				Double ja    = exp2[j];
				Double ic    = c1[i];
				Double jc    = c2[j];
				Double alpha = ia+ja; 
				Double ab    = -ia*ja/alpha;
				Double pref  = exp(ab*AB2);
				Double d     = ic*jc*pref;
				iExp.push_back(ia);
				iExp.push_back(ja);
				iCoe.push_back(d);
			}
		}

		// CD pairs
		UInt jNP2 = 9;
		vector<Double> jExp;
		vector<Double> jCoe;
		jExp.reserve(2*jNP2);
		jCoe.reserve(jNP2);
		Double CD2 = pow(C[0]-D[0],2)+pow(C[1]-D[1],2)+pow(C[2]-D[2],2);
		for(UInt j=0;j<3;j++) {
			for(UInt i=0;i<3;i++) {
				Double ia    = exp3[i];
				Double ja    = exp4[j];
				Double ic    = c3[i];
				Double jc    = c4[j];
				Double alpha = ia+ja; 
				Double ab    = -ia*ja/alpha;
				Double pref  = exp(ab*CD2);
				Double d     = ic*jc*pref;
				jExp.push_back(ia);
				jExp.push_back(ja);
				jCoe.push_back(d);
			}
		}

		// now do the work
		for(UInt i=0; i<ntest; i++) {

			UInt maxM = 8;
			if (i == 1) maxM = 11;
			if (i == 2) maxM = 0;
			oldSSSS(maxM,iNP2,&iCoe.front(),&iExp.front(),A,B, 
					jNP2,&jCoe.front(),&jExp.front(),C,D); 
		}
	}

	//
	// now do the SSSS in new way
	//
	bool donewSSSS = true;
	if (donewSSSS) {

		// AB pairs
		UInt inp2 = 9;
		vector<Double> iexp;
		vector<Double> icoe;
		vector<Double> ifac;
		vector<Double> P;
		iexp.reserve(inp2);
		icoe.reserve(inp2);
		ifac.reserve(inp2);
		P.reserve(3*inp2);
		Double AB2 = pow(A[0]-B[0],2)+pow(A[1]-B[1],2)+pow(A[2]-B[2],2);
		for(UInt j=0;j<3;j++) {
			for(UInt i=0;i<3;i++) {
				Double ia    = exp1[i];
				Double ja    = exp2[j];
				Double ic    = c1[i];
				Double jc    = c2[j];
				Double alpha = 1/(ia+ja); 
				iexp.push_back(alpha);
				icoe.push_back(ic*jc);
				Double ab    = -ia*ja/(ia+ja);
				Double pref  = exp(ab*AB2);
				Double d     = pref*pow(PI*alpha,1.5E0);
				ifac.push_back(d);
				Double adz   = ia*alpha;
				Double bdz   = ja*alpha;
				Double PX    = A[0]*adz + B[0]*bdz;
				Double PY    = A[1]*adz + B[1]*bdz;
				Double PZ    = A[2]*adz + B[2]*bdz;
				P.push_back(PX);
				P.push_back(PY);
				P.push_back(PZ);
			}
		}

		// CD pairs
		UInt jnp2 = 9;
		vector<Double> jexp;
		vector<Double> jcoe;
		vector<Double> jfac;
		vector<Double> Q;
		jexp.reserve(jnp2);
		jcoe.reserve(jnp2);
		jfac.reserve(jnp2);
		Q.reserve(3*jnp2);
		Double CD2 = pow(C[0]-D[0],2)+pow(C[1]-D[1],2)+pow(C[2]-D[2],2);
		for(UInt j=0;j<3;j++) {
			for(UInt i=0;i<3;i++) {
				Double ia    = exp3[i];
				Double ja    = exp4[j];
				Double ic    = c3[i];
				Double jc    = c4[j];
				Double alpha = 1/(ia+ja); 
				jexp.push_back(alpha);
				jcoe.push_back(ic*jc);
				Double ab    = -ia*ja/(ia+ja);
				Double pref  = exp(ab*CD2);
				Double d     = pref*pow(PI*alpha,1.5E0);
				jfac.push_back(d);
				Double adz   = ia*alpha;
				Double bdz   = ja*alpha;
				Double QX    = C[0]*adz + D[0]*bdz;
				Double QY    = C[1]*adz + D[1]*bdz;
				Double QZ    = C[2]*adz + D[2]*bdz;
				Q.push_back(QX);
				Q.push_back(QY);
				Q.push_back(QZ);
			}
		}

		// now do the work
		for(UInt i=0; i<ntest; i++) {

			UInt maxM = 8;
			if (i == 1) maxM = 11;
			if (i == 2) maxM = 0;
			newSSSS(maxM,inp2,&iexp.front(),&icoe.front(),&ifac.front(), 
					jnp2,&jexp.front(),&jcoe.front(),&jfac.front(),A,B,&P.front(),
					C,D,&Q.front()); 
		}
	}
	return 0;
}

