//
// we only list the latest check up for the codes
//
// 1 8th 2016:
//
// full test is made after we revise code in 
// github commit 10734ecd72e2fcaa45b2a8cc1933d237b7e0ac97
//
// the following modules has been well tested:
// two and three body KI, OV
// NAI
// all ERI
// MOM (auxL up to 5)
// are well tested, for L = 4, auxL = 5
// the test are performed for both single precision and 
// double precision
// threshold value for single is given as 1.0*10^-5
// threshold value for double is given as 1.0*10^-10
// all of tests for "double" are below the given thresh 
// value.
// see the sing.log and double.log for the test results.
//

#include "libgen.h"
#include "ov.h"
#include "ki.h"
#include "nai.h"
#include "tov.h"
#include "tki.h"
#include "eritest.h"
#include "mom.h"
#include "expr12test.h"
#include <boost/algorithm/string.hpp>   // string handling
#include <boost/lexical_cast.hpp>
using namespace ov;
using namespace tov;
using namespace tki;
using namespace ki;
using namespace nai;
using namespace eritest;
using namespace mom;
using namespace expr12test;

Int main(int argc, char* argv[])
{
	/////////////////////////////////////////////////////////////////////////////
	// setting the testing jobs
	// we note that the maxL must be in accordance with the hgp codes
	// so as the auxMaxL, order etc.
	/////////////////////////////////////////////////////////////////////////////
	bool testOV  = false;
	bool testKI  = false;
	bool testNAI = false;
	bool testERI = false;
	bool testTOV = false;
	bool testTKI = false;
	bool testMOM = false;
	bool test2BodyERI = false;
	bool test3BodyERI = false;
	bool testEXPR12   = false;

	// settings we need to further processed
	// here is default value
	Int order   = 0;
	Int maxL    = 4;
	Int auxMaxL = 5;
	Int momL    = 5;

	// parse the input parameter
	for(Int i=1; i<argc; i++) {
		string com = argv[i];
		boost::algorithm::to_lower(com);
		if (com == "ov"   ) testOV  = true;
		if (com == "tov"  ) testTOV = true;
		if (com == "tki"  ) testTKI = true;
		if (com == "ki"   ) testKI  = true;
		if (com == "nai"  ) testNAI = true;
		if (com == "eri"  ) testERI = true;
		if (com == "mom"  ) testMOM = true;
		if (com == "2eri" ) test2BodyERI = true;
		if (com == "3eri" ) test3BodyERI = true;
		if (com == "expr12") testEXPR12  = true;
	}

	// now print out the input information
	cout << "********job testing list********" << endl;
	cout << "integral testing for energy calculation" << endl;
	cout << "maximum angular momentum for integrals is " << maxL << endl;
	cout << "maximum aux angular momentum for integrals is " << auxMaxL << endl;
	cout << "maximum angular momentum for MOM integrals is " << momL << endl;
	if (testOV) cout <<  "two body overlap" << endl;
	if (testTOV) cout << "three body overlap" << endl;
	if (testKI) cout  << "two body kinetic" << endl;
	if (testTKI) cout << "three body kinetic" << endl;
	if (testNAI) cout << "two body nuclear attraction" << endl;
	if (testERI) cout << "normal ERI" << endl;
	if (test2BodyERI) cout << "two body ERI" << endl;
	if (test3BodyERI) cout << "three body ERI" << endl;
	if (testMOM) cout << "momentum integrals" << endl;
	if (testEXPR12) cout << "expr12 integrals" << endl;
	cout << "********************************" << endl;

	/////////////////////////////////////////////////////////////////////////////
	// now setting the shell data
	// all of shell data are in default to be composite shell
	// which contains two sub-shells
	/////////////////////////////////////////////////////////////////////////////
	
	// shell 1
	Int inp = 2;
	if (order > 0) inp = 1;
	vector<Double> iexp(inp);
	vector<Double> icoe(2*inp);
	if (order == 0) {
		iexp[0] = 1.8812885;
		iexp[1] = 0.5442493;
		icoe[0] = 0.1193324;
		icoe[1] = 0.1608542;
		icoe[2] = 0.8434564;
		icoe[3] = 0.0689991;
	}else{
		iexp[0] = 0.1442493;
		icoe[0] = 1.0;
		icoe[1] = 1.0;
	}
	Double A[3];
	A[0]    = 1.0;
	A[1]    = 0.0;
	A[2]    = 0.0;

	// shell 2
	Int jnp = 2;
	if (order > 0) jnp = 1;
	vector<Double> jexp(jnp);
	vector<Double> jcoe(2*jnp);
	if (order == 0) {
		jexp[0] = 0.5439130;
		jexp[1] = 0.0914760;
		jcoe[0] = 6.139703E-02;
		jcoe[1] = 3.061130E-01;
		jcoe[2] = 1.154890;
		jcoe[3] = 0.1891265;
	}else{
		jexp[0] = 0.0914760;
		jcoe[0] = 1.0;
		jcoe[1] = 1.0;
	}
	Double B[3];
	B[0]    = 0.0;
	B[1]    = 1.0;
	B[2]    = 0.0;

	// shell 3
	Int knp = 2;
	if (order > 0) knp = 1;
	vector<Double> kexp(knp);
	vector<Double> kcoe(2*knp);
	if (order == 0) {
		kexp[0] = 0.63628970;
		kexp[1] = 0.14786010;
		kcoe[0] = -0.09996723;
		kcoe[1] = 0.39951283;
		kcoe[2] = 0.70011547;
		kcoe[3] = 0.15591627;
	}else{
		kexp[0] = 0.04808870;
		kcoe[0] = 1.0;
		kcoe[1] = 1.0;
	}
	Double C[3];
	C[0]    = 0.0;
	C[1]    = 0.0;
	C[2]    = 1.0;

	// shell 4
	Int lnp = 2;
	if (order > 0) lnp = 1;
	vector<Double> lexp(lnp);
	vector<Double> lcoe(2*lnp);
	if (order == 0) {
		lexp[0] = 0.51982050;
		lexp[1] = 0.16906180;
		lcoe[0] = 0.06723;
		lcoe[1] = 0.51283;
		lcoe[2] = 0.30011547;
		lcoe[3] = 0.5591627;
	}else{
		lexp[0] = 0.23695610;
		lcoe[0] = 1.0;
		lcoe[1] = 1.0;
	}
	Double D[3];
	D[0]    = 0.0;
	D[1]    = 1.0;
	D[2]    = 1.0;

	// create the dummy shell 
	// the most important thing for the 
	// shell is that its exp is all zero
	Int dnp = 1;
	vector<Double> dexp(dnp,ZERO);
	vector<Double> dcoe(dnp,ONE);
	Double dummy[3];
	dummy[0] = 0.0;
	dummy[1] = 0.0;
	dummy[2] = 0.0;

	// nuclear centers as well as atomic numbers
	Int nAtoms = 3;
	vector<Double> coord(3*nAtoms);
	vector<Int> atomicNumbers(nAtoms);
	coord[0] = 0.0;
	coord[1] = 0.0;
	coord[2] = 1.0;
	coord[3] = 0.0;
	coord[4] = 1.0;
	coord[5] = 0.0;
	coord[6] = 1.0;
	coord[7] = 0.0;
	coord[8] = 0.0;
	atomicNumbers[0] = 1;
	atomicNumbers[1] = 2;
	atomicNumbers[2] = 3;

	/////////////////////////////////////////////////////////////////////////////
	//                  now this is the real working section
	/////////////////////////////////////////////////////////////////////////////
	if (testOV) {
		ov_test(maxL,inp,icoe,iexp,A,jnp,jcoe,jexp,B);
	}
	if (testTOV) {
		tov_test(maxL,inp,icoe,iexp,A,jnp,jcoe,jexp,B,knp,kcoe,kexp,C);
	}
	if (testTKI) {
		tki_test(maxL,inp,icoe,iexp,A,jnp,jcoe,jexp,B,knp,kcoe,kexp,C);
	}
	if (testMOM) {

		// the center for moment operator is set to 0
		Double cen[3];
		cen[0] = ZERO;
		cen[1] = ZERO;
		cen[2] = ZERO;
		mom_test(maxL,momL,inp,icoe,iexp,A,jnp,jcoe,jexp,B,cen);
	}
	if (testKI) {
		ki_test(maxL,inp,icoe,iexp,A,jnp,jcoe,jexp,B);
	}
	if (testNAI) {
		nai_test(maxL,inp,icoe,iexp,A,jnp,jcoe,jexp,B,coord,atomicNumbers);
	}
	if (testERI) {
		eri_test(maxL,auxMaxL,NORMAL_ERI,
				inp,icoe,iexp,A,
				jnp,jcoe,jexp,B,
				knp,kcoe,kexp,C,
				lnp,lcoe,lexp,D);
	}
	if (test2BodyERI) {
		eri_test(maxL,auxMaxL,TWO_BODY_ERI,
				inp,icoe,iexp,A,
				dnp,dcoe,dexp,dummy,
				knp,kcoe,kexp,C,
				dnp,dcoe,dexp,dummy);
	}
	if (test3BodyERI) {
		eri_test(maxL,auxMaxL,THREE_BODY_ERI,
				inp,icoe,iexp,A,
				jnp,jcoe,jexp,B,
				knp,kcoe,kexp,C,
				dnp,dcoe,dexp,dummy);
	}
	if (testEXPR12) {
		expr12_test(maxL,
				inp,icoe,iexp,A,
				jnp,jcoe,jexp,B,
				knp,kcoe,kexp,C);
	}
	cout << "all of jobs finished" << endl;

	return 0;

}
