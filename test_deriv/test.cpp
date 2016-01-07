//
// we only list the latest check up for the codes
//
// in the version e7a49ba58de5e920953877092628ac6758ba5d82
// we checked the gradient for NAI, two body overlap, KI
// and normal ERI for maxL = 4
//
// all of calculation are passed. the result is in gradient_double.log
//
// in version ac6270130fac825d3df86e2139af960310e9f793
// we checked the second derivatives for NAI, two body overlap, KI
// and normal ERI for maxL = 3
//
// all of calculation are passed. the result is in hessian_double.log
//

#include "libgen.h"
#include "ov.h"
#include "ki.h"
#include "nai.h"
#include "eritest.h"
#include <boost/algorithm/string.hpp>   // string handling
#include <boost/lexical_cast.hpp>
using namespace ov;
using namespace ki;
using namespace nai;
using namespace eritest;

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
	bool test2BodyERI = false;
	bool test3BodyERI = false;

	// settings you may need to change here
	// here is default value
	Int order   = 1;
#ifdef ORDER2
	order       = 2;
#endif
	Int maxL    = 4;
	if (order == 2) maxL = 3;
	Int auxMaxL = 5;
	Int momL    = 5;

	// parse the input parameter
	for(Int i=1; i<argc; i++) {
		string com = argv[i];
		boost::algorithm::to_lower(com);
		if (com == "ov"   ) testOV  = true;
		if (com == "ki"   ) testKI  = true;
		if (com == "nai"  ) testNAI = true;
		if (com == "eri"  ) testERI = true;
		if (com == "2eri" ) test2BodyERI = true;
		if (com == "3eri" ) test3BodyERI = true;
	}

	// now print out the input information
	cout << "********job testing list********" << endl;
	if (order == 0) cout << "integral testing for energy calculation" << endl;
	if (order == 1) cout << "integral testing for gradient calculation" << endl;
	if (order == 2) cout << "integral testing for Hessian calculation" << endl;
	cout << "maximum angular momentum for integrals is " << maxL << endl;
	cout << "maximum aux angular momentum for integrals is " << auxMaxL << endl;
	cout << "maximum angular momentum for MOM integrals is " << momL << endl;
	if (testOV) cout <<  "two body overlap" << endl;
	if (testKI) cout  << "two body kinetic" << endl;
	if (testNAI) cout << "two body nuclear attraction" << endl;
	if (testERI) cout << "normal ERI" << endl;
	if (test2BodyERI) cout << "two body ERI" << endl;
	if (test3BodyERI) cout << "three body ERI" << endl;
	cout << "********************************" << endl;

	/////////////////////////////////////////////////////////////////////////////
	// now setting the shell data
	// all of shell data are in default to be composite shell
	// which contains two sub-shells
	/////////////////////////////////////////////////////////////////////////////
	
	// shell 1
	Int inp = 1;
	vector<Double> iexp(inp);
	vector<Double> icoe(2*inp);
	iexp[0] = 0.1442493;
	icoe[0] = 1.0;
	icoe[1] = 1.0;
	Double A[3];
	A[0]    = 1.0;
	A[1]    = 0.0;
	A[2]    = 0.0;

	// shell 2
	Int jnp = 1;
	vector<Double> jexp(jnp);
	vector<Double> jcoe(2*jnp);
	jexp[0] = 0.0914760;
	jcoe[0] = 1.0;
	jcoe[1] = 1.0;
	Double B[3];
	B[0]    = 0.0;
	B[1]    = 1.0;
	B[2]    = 0.0;

	// shell 3
	Int knp = 1;
	vector<Double> kexp(knp);
	vector<Double> kcoe(2*knp);
	kexp[0] = 0.04808870;
	kcoe[0] = 1.0;
	kcoe[1] = 1.0;
	Double C[3];
	C[0]    = 0.0;
	C[1]    = 0.0;
	C[2]    = 1.0;

	// shell 4
	Int lnp = 1;
	vector<Double> lexp(lnp);
	vector<Double> lcoe(2*lnp);
	lexp[0] = 0.23695610;
	lcoe[0] = 1.0;
	lcoe[1] = 1.0;
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
		string derivFile = "deriv_infor_twobodyoverlap_1.txt"; 
		if (order == 2) {
			derivFile = "deriv_infor_twobodyoverlap_2.txt"; 
		}
		ov_test(derivFile,maxL,order,inp,icoe,iexp,A,jnp,jcoe,jexp,B);
	}
	if (testKI) {
		string derivFile = "deriv_infor_kinetic_1.txt"; 
		if (order == 2) {
			derivFile = "deriv_infor_kinetic_2.txt"; 
		}
		ki_test(derivFile,maxL,order,inp,icoe,iexp,A,jnp,jcoe,jexp,B);
	}
	if (testNAI) {
		string derivFile = "deriv_infor_nai_1.txt"; 
		if (order == 2) {
			derivFile = "deriv_infor_nai_2.txt"; 
		}
		nai_test(derivFile,maxL,order,inp,icoe,iexp,A,jnp,jcoe,jexp,B,coord,atomicNumbers);
	}
	if (testERI) {
		string derivFile = "deriv_infor_eri_1.txt"; 
		if (order == 2) {
			derivFile = "deriv_infor_eri_2.txt"; 
		}
		eri_test(derivFile,maxL,auxMaxL,order,NORMAL_ERI,
				inp,icoe,iexp,A,
				jnp,jcoe,jexp,B,
				knp,kcoe,kexp,C,
				lnp,lcoe,lexp,D);
	}
	if (test2BodyERI) {
		string derivFile = "deriv_infor_eri_1.txt"; 
		if (order == 2) {
			derivFile = "deriv_infor_eri_2.txt"; 
		}
		eri_test(derivFile,maxL,auxMaxL,order,TWO_BODY_ERI,
				inp,icoe,iexp,A,
				dnp,dcoe,dexp,dummy,
				knp,kcoe,kexp,C,
				dnp,dcoe,dexp,dummy);
	}
	if (test3BodyERI) {
		string derivFile = "deriv_infor_eri_1.txt"; 
		if (order == 2) {
			derivFile = "deriv_infor_eri_2.txt"; 
		}
		eri_test(derivFile,maxL,auxMaxL,order,THREE_BODY_ERI,
				inp,icoe,iexp,A,
				jnp,jcoe,jexp,B,
				knp,kcoe,kexp,C,
				dnp,dcoe,dexp,dummy);
	}

	cout << "all of jobs are finished " << endl;
	return 0;

}
