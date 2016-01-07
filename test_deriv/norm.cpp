#include "libgen.h"
#include "angmomlist.h"
#include "shellprop.h"
#include "functions.h"
#include "norm.h"
using namespace shellprop;
using namespace functions;
using namespace norm;

Double norm::normCartBasisSets(const Int& iShell, const Int& nPrim, const Double* c, 
		const Double* e)
{
	Int l,m,n;
	l = iShell;
	m = n = 0;
	Double sum = ZERO;
	for(Int i=0; i<nPrim; i++) {
		for(Int j=i; j<nPrim; j++) {
			Double exp2 = e[i]+e[j];
			Double coe  = c[i]*c[j];
			Double t    = ONE;
			if (exp2<THRESHOLD_MATH) {
				cout << "For the angular momentum " << iShell << endl;
				cout << "The sum of exp2 is too small: " << exp2 << endl;
				crash(true, "It's impossible to calculate Cartesian normalization factor");
			}
			t = coe*gammaFunInt(2*l,exp2)*gammaFunInt(2*m,exp2)*gammaFunInt(2*n,exp2);
			if (j==i) {
				sum += t;
			} else {
				sum += TWO*t;
			}
		}
	}
	Double N = ONE/sqrt(sum);
	return N;
}

void norm::scaleNormFactors(const Int& Lmin,const Int& Lmax,Double* N) 
{
	Int i = 0;
	for(Int iShell=Lmin;iShell<=Lmax;iShell++) {
		Int nbas= getCartBas(iShell,iShell);
		for(Int iBas=0;iBas<nbas;iBas++) {
			Int l = 0;
			Int m = 0;
			Int n = 0;
			getlmn(iShell,iShell,iBas,l,m,n);
			//Double n1 = sqrt(doubleFac(TWO*l-ONE)*doubleFac(TWO*m-ONE)*doubleFac(TWO*n-ONE));
			//Double n2 = sqrt(doubleFac(TWO*iShell-ONE));

			//
			// because we have to avoid overflow of double for large
			// L in double factor calculation, now we change into
			// new way
			//
			Double x = ONE;
			Int L = iShell;
			for(Int count=0; count<L; count++) {
				Int l0 = l-count;
				if (l0<1) l0 = 1;
				Int m0 = m-count;
				if (m0<1) m0 = 1;
				Int n0 = n-count;
				if (n0<1) n0 = 1;
				Int L0 = L-count;
				if (L0<1) L0 = 1;
				x = x*(TWO*l0-ONE)*(TWO*m0-ONE)*(TWO*n0-ONE)/(TWO*L0-ONE);
			}
			x = ONE/x;
			N[i] = sqrt(x);
			i++;
		}
	}
}

void norm::normCoe(const Int& lmin, const Int& lmax, 
		const vector<Double>& e, const vector<Double>& oric, vector<Double>& c)
{
	// for the coefficient array, we have to do more treatment. The
	// treatment below is that each Gaussian primitives also need
	// normalized! the final normalization step will be done later
	Int nPrim = e.size();
	for(Int iShell=lmin;iShell<=lmax;iShell++) {
		Int np = (iShell-lmin)*nPrim; 
		for(Int i=0; i<nPrim; i++){
			Double coe  = oric[np+i];
			Double s    = pow(e[i],(iShell+THREE/TWO)/TWO);
			c[np+i]     = coe*s;
		}
	}

	// normalization for basis set
	// here we assume that all of basis sets are in type of x^Ly^0z^0
	// later in the transformation calculation we will scale it back
	// see shelltransscale.h
	for(Int iShell=lmin;iShell<=lmax;iShell++) {
		Int np = (iShell-lmin)*nPrim; 
		Double norm = normCartBasisSets(iShell, nPrim, &c[np], &e[0]);
		for(Int i=0; i<nPrim; i++){
			c[np+i] = c[np+i]*norm;
		}
	}
}

void norm::makeC2(const Int& l1min, const Int& l1max, 
		const Int& l2min, const Int& l2max, const Int& inp, const Int& jnp,
		const vector<Double>& c1, const vector<Double>& c2, vector<Double>& coe2)
{
	// firstly do additional check
	// we should only have SP composite shells
	if (l1min == 0 && l1max > 1) {
		crash(true, "Only SP shell is allowed here in makec2, input shell 1 is illegal");
	}
	if (l2min == 0 && l2max > 1) {
		crash(true, "Only SP shell is allowed here in makec2, input shell 2 is illegal");
	}

	// now let's form the coefficient pairs
	// the input c1 and c2 should be default 
	// in SP shell form
	//
	// we could have four situations below:
	// 1  shell 1 is composite and shell 2 is composite
	// 2  shell 1 is composite and shell 2 is single   
	// 3  shell 1 is single    and shell 2 is composite
	// 4  shell 1 is single    and shell 2 is single   
	//
	Int nL1 = l1max-l1min+1;
	Int nL2 = l2max-l2min+1;
	Int count = 0;
	coe2.assign(inp*jnp*nL1*nL2,ZERO);
	for(Int j=0; j<nL2; j++) {
		for(Int i=0; i<nL1; i++) {
			for(Int jp=0; jp<jnp; jp++) {
				for(Int ip=0; ip<inp; ip++) {
					Double ic = c1[ip+i*inp];
					Double jc = c2[jp+j*jnp];
					coe2[count] = ic*jc;
					count++;
				}
			}
		}
	}
}

void norm::scale1BodyInts(const Int& LBra1Min, const Int& LBra1Max, Double* result) 
{

	// do we return it?
	// for L<=P, we do nothing here
	if (LBra1Max<=1) return;

	// first step, get the normalization ratios
	Int nBra1Bas= getCartBas(LBra1Min,LBra1Max);
	vector<Double> bra1Ratio(nBra1Bas,ONE);
	if (LBra1Max>1) {
		scaleNormFactors(LBra1Min,LBra1Max,&bra1Ratio.front());
	}	

	// now let's scale the result
	for(Int iBas=0; iBas<nBra1Bas; iBas++) {
		Double N = bra1Ratio[iBas];
		result[iBas] *= N;
	}
}

void norm::scale2BodyInts(const Int& LBra1Min, const Int& LBra1Max, 
		const Int& LBra2Min, const Int& LBra2Max, Double* result) 
{

	// do we return it?
	// for L<=P, we do nothing here
	if (LBra1Max<=1 && LBra2Max<=1) return;

	// first step, get the normalization ratios
	Int nBra1Bas= getCartBas(LBra1Min,LBra1Max);
	vector<Double> bra1Ratio(nBra1Bas,ONE);
	if (LBra1Max>1) {
		scaleNormFactors(LBra1Min,LBra1Max,&bra1Ratio.front());
	}	
	Int nBra2Bas= getCartBas(LBra2Min,LBra2Max);
	vector<Double> bra2Ratio(nBra2Bas,ONE);
	if (LBra2Max>1) {
		scaleNormFactors(LBra2Min,LBra2Max,&bra2Ratio.front());
	}	

	// now let's scale the result
	for(Int jBas=0; jBas<nBra2Bas; jBas++) {
		for(Int iBas=0; iBas<nBra1Bas; iBas++) {
			Double N = bra1Ratio[iBas]*bra2Ratio[jBas];
			result[iBas+jBas*nBra1Bas] *= N;
		}
	}
}

void norm::scale3BodyInts(const Int& LBra1Min, const Int& LBra1Max, 
		const Int& LBra2Min, const Int& LBra2Max, 
		const Int& LKet1Min, const Int& LKet1Max, Double* result) 
{

	// do we return it?
	// for L<=P, we do nothing here
	if (LBra1Max<=1 && LBra2Max<=1 && LKet1Max<=1) return;

	// first step, get the normalization ratios
	Int nBra1Bas= getCartBas(LBra1Min,LBra1Max);
	vector<Double> bra1Ratio(nBra1Bas,ONE);
	if (LBra1Max>1) {
		scaleNormFactors(LBra1Min,LBra1Max,&bra1Ratio.front());
	}	
	Int nBra2Bas= getCartBas(LBra2Min,LBra2Max);
	vector<Double> bra2Ratio(nBra2Bas,ONE);
	if (LBra2Max>1) {
		scaleNormFactors(LBra2Min,LBra2Max,&bra2Ratio.front());
	}	
	Int nKet1Bas= getCartBas(LKet1Min,LKet1Max);
	vector<Double> ket1Ratio(nKet1Bas,ONE);
	if (LKet1Max>1) {
		scaleNormFactors(LKet1Min,LKet1Max,&ket1Ratio.front());
	}	

	// now let's scale the result
	for(Int kBas=0; kBas<nKet1Bas; kBas++) {
		for(Int jBas=0; jBas<nBra2Bas; jBas++) {
			for(Int iBas=0; iBas<nBra1Bas; iBas++) {
				Double N = bra1Ratio[iBas]*bra2Ratio[jBas]*ket1Ratio[kBas];
				result[iBas+jBas*nBra1Bas+kBas*nBra1Bas*nBra2Bas] *= N;
			}
		}
	}
}

void norm::scale4BodyInts(const Int& LBra1Min, const Int& LBra1Max, 
		const Int& LBra2Min, const Int& LBra2Max, 
		const Int& LKet1Min, const Int& LKet1Max, 
		const Int& LKet2Min, const Int& LKet2Max, Double* result) 
{

	// do we return it?
	// for L<=P, we do nothing here
	if (LBra1Max<=1 && LBra2Max<=1 && LKet1Max<=1 && LKet2Max<=1) return;

	// first step, get the normalization ratios
	Int nBra1Bas= getCartBas(LBra1Min,LBra1Max);
	vector<Double> bra1Ratio(nBra1Bas,ONE);
	if (LBra1Max>1) {
		scaleNormFactors(LBra1Min,LBra1Max,&bra1Ratio.front());
	}	
	Int nBra2Bas= getCartBas(LBra2Min,LBra2Max);
	vector<Double> bra2Ratio(nBra2Bas,ONE);
	if (LBra2Max>1) {
		scaleNormFactors(LBra2Min,LBra2Max,&bra2Ratio.front());
	}	
	Int nKet1Bas= getCartBas(LKet1Min,LKet1Max);
	vector<Double> ket1Ratio(nKet1Bas,ONE);
	if (LKet1Max>1) {
		scaleNormFactors(LKet1Min,LKet1Max,&ket1Ratio.front());
	}	
	Int nKet2Bas= getCartBas(LKet2Min,LKet2Max);
	vector<Double> ket2Ratio(nKet2Bas,ONE);
	if (LKet2Max>1) {
		scaleNormFactors(LKet2Min,LKet2Max,&ket2Ratio.front());
	}	

	// now let's scale the result
	for(Int lBas=0; lBas<nKet2Bas; lBas++) {
		for(Int kBas=0; kBas<nKet1Bas; kBas++) {
			for(Int jBas=0; jBas<nBra2Bas; jBas++) {
				for(Int iBas=0; iBas<nBra1Bas; iBas++) {
					Double N = bra1Ratio[iBas]*bra2Ratio[jBas]*ket1Ratio[kBas]*ket2Ratio[lBas];
					result[iBas+jBas*nBra1Bas+kBas*nBra1Bas*nBra2Bas+
						lBas*nBra1Bas*nBra2Bas*nKet1Bas] *= N;
				}
			}
		}
	}
}
