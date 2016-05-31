//
// this is for testing the timing of the given individual integrals
//

#include "libgen.h"
#include <boost/algorithm/string.hpp>   // string handling
#include <boost/lexical_cast.hpp>
#include <boost/timer.hpp>
using namespace boost;

extern void hgp_os_eri_sp_sp_sp_sp(const UInt& inp2, const UInt& jnp2, const Double& pMax, 
		const Double& omega, const Double* icoe, const Double* iexp, const Double* ifac, 
		const Double* P, const Double* A, const Double* B, const Double* jcoe, const Double* jexp, 
		const Double* jfac, const Double* Q, const Double* C, const Double* D, Double* abcd);

extern void hgp_os_eri_sp_sp_sp_sp_d1(const UInt& inp2, const UInt& jnp2, const Double& pMax, 
		const Double& omega, const Double* icoe, const Double* iexp, const Double* iexpdiff, 
		const Double* ifac, const Double* P, const Double* A, const Double* B, const Double* jcoe, 
		const Double* jexp, const Double* jexpdiff, const Double* jfac, const Double* Q, const Double* C, 
		const Double* D, Double* abcd);

Int main(int argc, char* argv[])
{
	/////////////////////////////////////////////////////////////////////////////
	// setting the shell data
	// all of shell data are in default to be composite shell
	// which contains two sub-shells
	/////////////////////////////////////////////////////////////////////////////
	
	// shell 1
	Int inp = 2;
	vector<Double> iexp(inp);
	vector<Double> icoe(2*inp);
	iexp[0] = 0.1812885;
	iexp[1] = 0.0442493;
	icoe[0] = 0.1193324;
	icoe[1] = 0.1608542;
	icoe[2] = 0.8434564;
	icoe[3] = 0.0689991;
	Double A[3];
	A[0]    = 1.0;
	A[1]    = 0.0;
	A[2]    = 0.0;

	// shell 2
	Int jnp = 2;
	vector<Double> jexp(jnp);
	vector<Double> jcoe(2*jnp);
	jexp[0] = 0.1439130;
	jexp[1] = 0.0914760;
	jcoe[0] = 6.139703E-02;
	jcoe[1] = 3.061130E-01;
	jcoe[2] = 1.154890;
	jcoe[3] = 0.1891265;
	Double B[3];
	B[0]    = 0.0;
	B[1]    = 1.0;
	B[2]    = 0.0;

	// shell 3
	Int knp = 2;
	vector<Double> kexp(knp);
	vector<Double> kcoe(2*knp);
	kexp[0] = 0.03628970;
	kexp[1] = 0.14786010;
	kcoe[0] = 0.09996723;
	kcoe[1] = 0.39951283;
	kcoe[2] = 0.70011547;
	kcoe[3] = 0.15591627;
	Double C[3];
	C[0]    = 0.0;
	C[1]    = 0.0;
	C[2]    = 1.0;

	// shell 4
	Int lnp = 2;
	vector<Double> lexp(lnp);
	vector<Double> lcoe(2*lnp);
	lexp[0] = 0.01982050;
	lexp[1] = 0.16906180;
	lcoe[0] = 0.06723;
	lcoe[1] = 0.51283;
	lcoe[2] = 0.30011547;
	lcoe[3] = 0.5591627;
	Double D[3];
	D[0]    = 0.0;
	D[1]    = 1.0;
	D[2]    = 1.0;

	/////////////////////////////////////////////////////////////////////////////
	// setting the shell data
	// all of shell data are in default to be composite shell
	// which contains two sub-shells
	/////////////////////////////////////////////////////////////////////////////
	//
	// form the data for hgp calculation
	// firstly it's the bra side
	//
	Double AB2 = (A[0]-B[0])*(A[0]-B[0])+(A[1]-B[1])*(A[1]-B[1])+(A[2]-B[2])*(A[2]-B[2]);
	Int inp2 = inp*jnp;
	vector<Double> iexp2(inp2,ZERO);
	vector<Double> fbra(inp2,ZERO);
	vector<Double> P(3*inp2,ZERO);
	vector<Double> iexpdiff(inp2,ZERO);
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
			fbra[count] = pref;
			iexpdiff[count]= diff;

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
			fket[count] = pref;
			jexpdiff[count]= diff;

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

	/////////////////////////////////////////////////////////////////////////////
	// now it's the formal timing compare
	/////////////////////////////////////////////////////////////////////////////
	
	// set the angular momentum of the results
	// here the user need to do that
	Int lmin1 = 0;
	Int lmax1 = 1;
	Int lmin2 = 0;
	Int lmax2 = 1;
	Int lmin3 = 0;
	Int lmax3 = 1;
	Int lmin4 = 0;
	Int lmax4 = 1;
	Int nBas1 = ((lmax1+1)*(lmax1+2)*(lmax1+3)-lmin1*(lmin1+1)*(lmin1+2))/6;
	Int nBas2 = ((lmax2+1)*(lmax2+2)*(lmax2+3)-lmin2*(lmin2+1)*(lmin2+2))/6;
	Int nBas3 = ((lmax3+1)*(lmax3+2)*(lmax3+3)-lmin3*(lmin3+1)*(lmin3+2))/6;
	Int nBas4 = ((lmax4+1)*(lmax4+2)*(lmax4+3)-lmin4*(lmin4+1)*(lmin4+2))/6;

	// make the coefficients pair
	// bra side
	Int nL1 = lmax1-lmin1+1;
	Int nL2 = lmax2-lmin2+1;
	count = 0;
	vector<Double> braCoePair(inp*jnp*nL1*nL2,ZERO);
	for(Int j=0; j<nL2; j++) {
		for(Int i=0; i<nL1; i++) {
			for(Int jp=0; jp<jnp; jp++) {
				for(Int ip=0; ip<inp; ip++) {
					Double ic = icoe[ip+i*inp];
					Double jc = jcoe[jp+j*jnp];
					braCoePair[count] = ic*jc;
					count++;
				}
			}
		}
	}

	// ket side
	nL1 = lmax3-lmin3+1;
	nL2 = lmax4-lmin4+1;
	count = 0;
	vector<Double> ketCoePair(knp*lnp*nL1*nL2,ZERO);
	for(Int j=0; j<nL2; j++) {
		for(Int i=0; i<nL1; i++) {
			for(Int jp=0; jp<lnp; jp++) {
				for(Int ip=0; ip<knp; ip++) {
					Double ic = kcoe[ip+i*inp];
					Double jc = lcoe[jp+j*jnp];
					ketCoePair[count] = ic*jc;
					count++;
				}
			}
		}
	}

	// now set up the result vectors
	// N is the total number of running times
	Int N = 1000000;
	Double pMax   = 1.0E0;
	Double omega  = 0.0E0;
	vector<Double> result1(nBas1*nBas2*nBas3*nBas4);
	timer t1;
	for(Int i=0; i<N; i++) {
		hgp_os_eri_sp_sp_sp_sp(inp2,jnp2,pMax,omega,
				&braCoePair.front(),&iexp2.front(),&fbra.front(),&P.front(),A,B, 
				&ketCoePair.front(),&jexp2.front(),&fket.front(),&Q.front(),C,D, 
				&result1.front());
	}
	printf("hgp sp_sp_sp_sp energy time consuming: %-14.7f\n", t1.elapsed());

	// this is another integral file
	vector<Double> result2(nBas1*nBas2*nBas3*nBas4*9); // dimension is responsible by the user
	timer t2;
	for(Int i=0; i<N; i++) {
		hgp_os_eri_sp_sp_sp_sp_d1(inp2,jnp2,pMax,omega,
				&braCoePair.front(),&iexp2.front(),&iexpdiff.front(),&fbra.front(),&P.front(),A,B, 
				&ketCoePair.front(),&jexp2.front(),&jexpdiff.front(),&fket.front(),&Q.front(),C,D, 
				&result2.front());
	}
	printf("hgp sp_sp_sp_sp gradient time consuming: %-14.7f\n", t2.elapsed());

	return 0;

}
