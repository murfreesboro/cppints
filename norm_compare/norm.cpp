#include "libgen.h"
#include "functions.h"
using namespace functions;

Double normCartBasisSets(const Int& iShell, const Int& nPrim, const Double* c, 
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

Int main()
{
	// set shell data
	Int np = 1;
	vector<Double> jexp(np);
	vector<Double> jcoe(np);
	jexp[0] = 0.03;
	//jexp[1] = 0.5439130;
	//jexp[2] = 0.0914760;
	jcoe[0] = 1.0;
	//jcoe[1] = 1.0;
	//jcoe[2] = 1.0;

	// see the normalization factor
	Int n = 6;
	Double normS = normCartBasisSets(0,np,&jcoe.front(),&jexp.front());
	for(Int i=1; i<=n; i++) {
		Double n = normCartBasisSets(i,np,&jcoe.front(),&jexp.front());
		cout << "order " << i << " ratio " << n/normS << endl;
	}

	return 0;

}
