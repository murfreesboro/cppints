//
// CPPINTS: A C++ Program to Generate Analytical Integrals Based on Gaussian
// Form Primitive Functions
//
// Copyright (C) 2015 Fenglai Liu
// This softare uses the MIT license as below:
//
//	Permission is hereby granted, free of charge, to any person obtaining 
//	a copy of this software and associated documentation files (the "Software"), 
//	to deal in the Software without restriction, including without limitation 
//	the rights to use, copy, modify, merge, publish, distribute, sublicense, 
//	and/or sell copies of the Software, and to permit persons to whom the Software 
//	is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all 
// copies or substantial portions of the Software.
//						    
//	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
//	INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
//	PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
//	FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
//	ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
//
#include "boost/lexical_cast.hpp"
#include <boost/algorithm/string.hpp>   // string handling
#include <boost/filesystem.hpp>
#include "inttype.h"
#include "rr.h"
#include "printing.h"
#include "sqintsprint.h"
#include "sqints.h"
using boost::lexical_cast;
using namespace boost::filesystem;
using namespace boost;
using namespace printing;
using namespace inttype;
using namespace rr;
using namespace sqintsprint;
using namespace sqints;

void SQInts::headPrinting(ofstream& file) const 
{
	// first part, the comment of software license
	string line = "//";
	printLine(0,line,file);
	line = "// ";
	printLine(0,line,file);
	line = "// This code is generated from CPPINTS, a C++ program to generate the ";
	printLine(0,line,file);
	line = "// analytical integrals based on Gaussian form primitive functions. ";
	printLine(0,line,file);
   line = "// Copyright (C) 2015 Fenglai Liu ";
	printLine(0,line,file);
   line = "// This softare uses the MIT license as below: ";
	printLine(0,line,file);
   line = "// ";
	printLine(0,line,file);
   line = "// Permission is hereby granted, free of charge, to any person obtaining ";
	printLine(0,line,file);
   line = "// a copy of this software and associated documentation files (the \"Software\"), ";
	printLine(0,line,file);
   line = "// to deal in the Software without restriction, including without limitation ";
	printLine(0,line,file);
   line = "// the rights to use, copy, modify, merge, publish, distribute, sublicense, ";
	printLine(0,line,file);
   line = "// and/or sell copies of the Software, and to permit persons to whom the Software ";
	printLine(0,line,file);
   line = "// is furnished to do so, subject to the following conditions: ";
	printLine(0,line,file);
   line = "// ";
	printLine(0,line,file);
   line = "// The above copyright notice and this permission notice shall be included in all ";
	printLine(0,line,file);
   line = "// copies or substantial portions of the Software. ";
	printLine(0,line,file);
   line = "// ";					    
	printLine(0,line,file);
   line = "// THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, ";
	printLine(0,line,file);
   line = "// INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR ";
	printLine(0,line,file);
   line = "// PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE ";
	printLine(0,line,file);
   line = "// FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR ";
	printLine(0,line,file);
	line = "// OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER  ";
	printLine(0,line,file);
	line = "// DEALINGS IN THE SOFTWARE. ";
	printLine(0,line,file);
	line = "// ";
	printLine(0,line,file);
	line = "// ";
	printLine(0,line,file);
	file << endl;

	// do we use vector or scr etc.?
	bool usescr = infor.withSCRVec();
	bool useDoubleVec = infor.withDoubleVec();
	bool withTBBVec = infor.useTBBVec();
	bool withBoostGamma = infor.useBoostGamma();

	// now second part, it's the include file
	line = "#include \"constants.h\""; 
	printLine(0,line,file);
	line = "#include <cstddef>"; 
	printLine(0,line,file);
	line = "#include <math.h>"; 
	printLine(0,line,file);
	if (withBoostGamma) {
		line = "#include <boost/math/special_functions/gamma.hpp>";
		printLine(0,line,file);
	}
	if (usescr) {
		line = "#include \"localmemscr.h\""; 
		printLine(0,line,file);
		line = "using namespace localmemscr;";
		printLine(0,line,file);
	}else {
		if (useDoubleVec){
			line = "#include <vector>"; 
			printLine(0,line,file);
			if (withTBBVec) {
				line = "#include \"tbb/scalable_allocator.h\"";
				printLine(0,line,file);
			}
		}
	}
	file << endl;

	// now do the typedef work
	// however, the information here is defined in the localmemscr
	// therefore, we do not repeat it if we use LocalMemScr
	if (! usescr) {
		line = "typedef int             Int;";
		printLine(0,line,file);
		line = "typedef size_t          UInt;";
		printLine(0,line,file);
		line = "#ifdef WITH_SINGLE_PRECISION";
		printLine(0,line,file);
		line = "typedef float           Double;";
		printLine(0,line,file);
		line = "#define THRESHOLD_MATH  0.0000001";
		printLine(0,line,file);
		line = "#else";
		printLine(0,line,file);
		line = "typedef double          Double;";
		printLine(0,line,file);
		line = "#define THRESHOLD_MATH  0.00000000000001";
		printLine(0,line,file);
		line = "#endif";
		printLine(0,line,file);
		file << endl;
	}

	// now if we use vector, we also do typedef here
	if (useDoubleVec) {
		line = "// typedef the vector type so that they all have the same type of DoubleVec";
		if (withTBBVec) {
			line = "typedef std::vector<Double,tbb::scalable_allocator<Double> >   DoubleVec; ";
		}else{
			line = "typedef std::vector<Double>   DoubleVec; ";
		}
		file << endl;
	}

	// print out variable comments
	line = "//";
	printLine(0,line,file);
	line = "//  here below is a list of variables used in the program";
	printLine(0,line,file);
	line = "//";
	printLine(0,line,file);
	line = "//  alpha is the bra1's exponent";
	printLine(0,line,file);
	line = "//  beta  is the bra2's exponent";
	printLine(0,line,file);
	line = "//  gamma is the ket1's exponent";
	printLine(0,line,file);
	line = "//  delta is the ket2's exponent";
	printLine(0,line,file);
	line = "//  A is the nuclear center for bra1";
	printLine(0,line,file);
	line = "//  B is the nuclear center for bra2";
	printLine(0,line,file);
	line = "//  C is the nuclear center for ket1";
	printLine(0,line,file);
	line = "//  D is the nuclear center for ket2";
	printLine(0,line,file);
	line = "//  P is the new center after bra1 combined with bra2";
	printLine(0,line,file);
	line = "//  Q is the new center after ket1 combined with ket2";
	printLine(0,line,file);
	line = "//  W is the new center after P combined with Q";
	printLine(0,line,file);
	line = "//  thresh value is threshold to perform significance check on primitive integrals";
	printLine(0,line,file);
	line = "//  pMax is maximum value of corresponding density matrix block, used for ERI";
	printLine(0,line,file);
	line = "//";
	printLine(0,line,file);
	line = "//  variables:";
	printLine(0,line,file);
	line = "//";
	printLine(0,line,file);
	line = "//  zeta      = alpha + beta";
	printLine(0,line,file);
	line = "//  eta       = gamma + delta";
	printLine(0,line,file);
	line = "//  oned2z    = 1/(2*zeta)";
	printLine(0,line,file);
	line = "//  oned2e    = 1/(2*eta)";
	printLine(0,line,file);
	line = "//  onedz     = 1/zeta";
	printLine(0,line,file);
	line = "//  onede     = 1/eta";
	printLine(0,line,file);
	line = "//  kappa     = zeta + eta";
	printLine(0,line,file);
	line = "//  onedk     = 1/kappa";
	printLine(0,line,file);
	line = "//  oned2zeta = 1/(2*(alpha+beta+gamma))";
	printLine(0,line,file);
	line = "//  xi        = alpha*beta*onedz";
	printLine(0,line,file);
	line = "//  twoxi     = 2*alpha*beta*onedz";
	printLine(0,line,file);
	line = "//  rho       = zeta*eta*onedk";
	printLine(0,line,file);
	line = "//  rhod2zsq  = rho/(2*zeta*zeta)";
	printLine(0,line,file);
	line = "//  rhod2esq  = rho/(2*eta*eta)";
	printLine(0,line,file);
	line = "//  adz       = alpha*onedz";
	printLine(0,line,file);
	line = "//  bdz       = beta*onedz";
	printLine(0,line,file);
	line = "//  gde       = gamma*onede";
	printLine(0,line,file);
	line = "//  gde       = delta*onede";
	printLine(0,line,file);
	line = "//";
	printLine(0,line,file);
	line = "//  input parameters based on primitive functions pair:";
	printLine(0,line,file);
	line = "//";
	printLine(0,line,file);
	line = "//  bra side shell pair is index as i";
	printLine(0,line,file);
	line = "//  inp2  is the number of primitive pairs";
	printLine(0,line,file);
	line = "//  iexp  is the array of 1/(alpha+beta)";
	printLine(0,line,file);
	line = "//  icoe  is the array of ic_bra1*ic_bra2";
	printLine(0,line,file);
	line = "//  ifac  is the array of pre-factor on bra side ";
	printLine(0,line,file);
	line = "//  for (SS|SS)^{m} etc. type of integrals";
	printLine(0,line,file);
	line = "//  ket side shell pair is index as j";
	printLine(0,line,file);
	line = "//  jnp2  is the number of primitive pairs";
	printLine(0,line,file);
	line = "//  jexp  is the array of 1/(gamma+delta)";
	printLine(0,line,file);
	line = "//  jcoe  is the array of jc_ket1*jc_ket2";
	printLine(0,line,file);
	line = "//  jfac  is the array of pre-factor on ket side ";
	printLine(0,line,file);
	line = "//  for (SS|SS)^{m} etc. type of integrals";
	printLine(0,line,file);
	line = "//";
	printLine(0,line,file);
	file << endl;
}

string SQInts::getArgList() const
{
	string arg;
	int intOperator = infor.getOper();
	switch(intOperator) {
		case TWOBODYOVERLAP:
			arg = "const UInt& inp2, const Double* icoe, "
				"const Double* iexp, const Double* ifac, const Double* P, " 
				"const Double* A, const Double* B, Double* abcd";
			break;
		case MOM:
			arg = "const UInt& inp2, const Double* icoe, "
				"const Double* iexp, const Double* ifac, const Double* P, " 
				"const Double* A, const Double* B, const Double* C, Double* abcd";
			break;
		case THREEBODYOVERLAP:
			arg = "const UInt& inp2, const UInt& jnp2, const Double* icoe, "
				"const Double* iexp, const Double* ifac, const Double* P, "
				"const Double* A, const Double* B, const Double* jcoe, "
				"const Double* jexp, const Double* C, Double* abcd";
			break;
		case THREEBODYKI:
			arg = "const UInt& inp2, const UInt& jnp2, const Double* icoe, "
				"const Double* iexp, const Double* iexpdiff, "
				"const Double* ifac, const Double* P, "
				"const Double* A, const Double* B, const Double* jcoe, "
				"const Double* jexp, const Double* C, Double* abcd";
			break;
		case NAI:
			arg = "const UInt& inp2, const UInt& nAtoms, const Double* icoe, " 
				"const Double* iexp, const Double* ifac, const Double* P, "
				"const Double* A, const Double* B, const Double* N, const UInt* Z, " 
				"Double* abcd";
			break;
		case ESP:
			arg = "const UInt& inp2, const UInt& nGrids, const Double* icoe, " 
				"const Double* iexp, const Double* ifac, const Double* P, "
				"const Double* A, const Double* B, const Double* R, " 
				"Double* abcd";
			break;
		case ERI:
			arg = "const UInt& inp2, const UInt& jnp2, const Double& thresh, const Double& pMax, "
				"const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, "
				"const Double* A, const Double* B, const Double* jcoe, "
				"const Double* jexp, const Double* jfac, const Double* Q, "
				"const Double* C, const Double* D, Double* abcd";
			break;
		case KINETIC:
			arg = "const UInt& inp2, const Double* icoe, "
				"const Double* iexp, const Double* iexpdiff, const Double* ifac, " 
				"const Double* P, const Double* A, const Double* B, Double* abcd";
			break;
		default:
			crash(true, "Invalid operator passed in getArgList");
			break;
	}

	// finally, consider that whether we have the scr class add in?
	if (infor.withSCRVec()) {
		arg = arg + ", LocalMemScr& scr";
	}

	return arg;
}

void SQInts::doCoreHRR(const int& iSide, const vector<ShellQuartet>& inputList, 
		vector<ShellQuartet>& bottomSQList) const 
{
	// determine the side
	int side = infor.get2edSide();
	string fileName = infor.getFileName();
	if (iSide == 0) {
		fileName += ".hrr2";
	}else if (iSide == 1) {
		side = infor.get1stSide();
		fileName += ".hrr1";
	}

	// is it worthy to do HRR?
	bool weDoHRR = false;
	if (side>0) {
		for(int iSQ=0; iSQ<(int)inputList.size(); iSQ++) {
			const ShellQuartet& sq = inputList[iSQ];
			if (sq.canDoHRR(side)) {
				weDoHRR = true;
				break;
			}
		}
	}

	// testing code
#ifdef SQINTS_DEBUG
	cout << "Shall we do the HRR for the given input shell quartet? " << weDoHRR << endl;
	for(int iSQ=0; iSQ<(int)inputList.size(); iSQ++) {
		cout << inputList[iSQ].getName() << endl;
	}
#endif

	// now do the hrr work
	// if we can do the hrr work, we print out 
	// the hrr codes and push the result bottom
	// sq in the result list
	// else we just push the input sq list
	// into the result list for further work
	if (weDoHRR) {
		RR hrr(HRR,inputList,infor,side);
		hrr.print(infor,fileName);
		hrr.getBottomSQList(bottomSQList); 
	}else{
		for(int iSQ=0; iSQ<(int)inputList.size(); iSQ++) {
			const ShellQuartet& sq = inputList[iSQ];
			vector<ShellQuartet>::iterator it = find(bottomSQList.begin(),bottomSQList.end(),sq);
			if (it == bottomSQList.end()) {
				bottomSQList.push_back(sq);
			}
		}
	}
}

void SQInts::doHRR(vector<ShellQuartet>& hrrResultSQList) const
{
	// shall we do HRR work?
	if (infor.hasHrr()) {

		// set up the print 
		const vector<ShellQuartet>& inputSQList = infor.getInputSQList();
		SQIntsPrint sqintsPrint(HRR,infor,inputSQList);

		// print out the necessary variables for HRR
		// i == 0 is second side, i == 1 is first side
		for(int i=0; i<2; i++) {
			string fileName = infor.getFileName();
			if (i == 0) {
				int secondSide = infor.get2edSide();
				if (secondSide<0) continue;
				fileName += ".hrr2";
				sqintsPrint.printHRRSideVar(secondSide,fileName);
			}else if (i == 1) {
				int firstSide = infor.get1stSide();
				if (firstSide<0) continue;
				fileName += ".hrr1";
				sqintsPrint.printHRRSideVar(firstSide,fileName);
			}
		}

		// we note, that each sq in the inputSQList has a 
		// indepedent HRR process. Therefore, for each
		// single sq in the composite sq we will do its 
		// HRR one by one
		for(int iSQ=0; iSQ<(int)inputSQList.size(); iSQ++) {

			// initial sq
			const ShellQuartet& sq = inputSQList[iSQ];

			// form the input vector for hrr
			vector<ShellQuartet> inputList;
			inputList.reserve(1);
			inputList.push_back(sq);

			// hrr result for this sq
			vector<ShellQuartet> hrrResult;

			// now do the work
			int iSide = 0;
			doCoreHRR(iSide,inputList,hrrResult); 

			// testing code
#ifdef SQINTS_DEBUG
			cout << "hrr result sq list from second side:" << endl;
			for(int i=0; i<(int)hrrResult.size(); i++) {
				cout << hrrResult[i].getName() << endl;
			}
#endif
			// 
			// now let's do the first side
			//
			inputList = hrrResult;
			hrrResult.clear();

			// now do the work
			iSide = 1;
			doCoreHRR(iSide,inputList,hrrResult); 

			// now we need to push the hrrResult into the final 
			// result list
			for(int i=0; i<(int)hrrResult.size(); i++) {
				const ShellQuartet& sq = hrrResult[i];
				vector<ShellQuartet>:: iterator it=find(hrrResultSQList.begin(),
						hrrResultSQList.end(), sq);
				if (it == hrrResultSQList.end()) hrrResultSQList.push_back(sq);
			}
		}

		// testing code
#ifdef SQINTS_DEBUG
		cout << "final hrr result sq list:" << endl;
		for(int iSQ=0; iSQ<(int)hrrResultSQList.size(); iSQ++) {
			cout << hrrResultSQList[iSQ].getName() << endl;
		}
#endif

	}else{

		// if we do not do the HRR process, then the input shell quartet
		// is carried out into the VRR process
		const vector<ShellQuartet>& inputSQList = infor.getInputSQList();
		hrrResultSQList = inputSQList;
	}

	// finally, we have to see that whether we are going to generate
	// independent hrr cpp file
	// the core cpp code is already contained in the hrr1 and hrr2 file
	// here what we need to complete the whole task is the function head
	// information
	if (infor.splitCPPFile() && infor.hasHrr()) {

		// set up the print 
		// here we need to set up the calling of HRR
		// therefore we will use VRR's result
		// which is actually the above hrrResultSQList
		SQIntsPrint sqintsPrint(HRR,infor,hrrResultSQList);

		// get the file name and open the stream
		string fileName = infor.getFileName();
		fileName += ".hrr_head";
		ofstream cppHead;
		cppHead.open(fileName.c_str());

		// print out function type information to this file
		string functionName = infor.getFuncName(); 
		functionName = functionName + "_hrr";
		string argList = sqintsPrint.getHRRArgList();
		string line = "void " + functionName + "( " + argList + " )";
		printLine(0,line,cppHead);
		cppHead.close();

		// now we also need the string that how can we call the function
		// in the code
		fileName  = infor.getFileName();
		fileName += ".hrr_func";
		ofstream cppCode;
		cppCode.open(fileName.c_str());

		// print out function line to this file
		string argListInCall = sqintsPrint.transformArgList(argList);
		line = functionName + "( " + argListInCall + " );";
		printLine(2,line,cppCode);
		cppCode << endl;
		cppCode.close();
	}
}

void SQInts::doVRR(const vector<ShellQuartet>& hrrResultSQList) const 
{
	// set up the sqints print
	int vrr_method = infor.getVRRMethod();
	int oper = infor.getOper();
	SQIntsPrint sqintsPrint(vrr_method,infor,hrrResultSQList);

	// VRR head
	string fileName = infor.getFileName();
	string varFile  = fileName + ".vrr";
	if (infor.splitCPPFile()) varFile = fileName + ".vrr_var";
	sqintsPrint.printVRRHead(varFile);

	// do we need to add in bottom integral generations 
	// into vrr file
	string vrrFile  = fileName + ".vrr";
	if (infor.splitCPPFile() && complicatedBottomIntegrals(oper)) {
		sqintsPrint.printBottomIntegrals(vrrFile);
	}

	// if this is the pure S integral, then no actually RR process
	// needed. Therefore, we will return here
	if (infor.areAllBottomSQ()) return;

	// we need to convert the hrr result sq list into the vrr input
	// that is to say, drop all of division inforamtion
	vector<ShellQuartet> vrrInputList;
	for(int i=0; i<(int)hrrResultSQList.size(); i++) {
		ShellQuartet sq(hrrResultSQList[i]);
		sq.destroyDivision();
		vector<ShellQuartet>:: iterator it=find(vrrInputList.begin(),
				vrrInputList.end(), sq);
		if (it == vrrInputList.end()) vrrInputList.push_back(sq);
	}

	// now let's do VRR process
	RR vrr(vrr_method,vrrInputList,infor,NULL_POS);
	vrr.print(infor,vrrFile);

	// perform contraction for VRR if necessary
	sqintsPrint.vrrContraction(vrrFile);

	// finally, we have to see that whether we are going to generate
	// independent vrr cpp file
	// the core cpp code is already contained in the vrr file
	// here what we need to complete the whole task is the function head
	// information
	if (infor.splitCPPFile()) {

		// get the file name and open the stream
		string fileName = infor.getFileName();
		fileName += ".vrr_head";
		ofstream cppHead;
		cppHead.open(fileName.c_str());

		// print out function type information to this file
		string functionName = infor.getFuncName(); 
		functionName = functionName + "_vrr";
		string argList = sqintsPrint.getVRRArgList();
		string returnType = "void ";
		string line = returnType + functionName + "( " + argList + " )";
		printLine(0,line,cppHead);
		cppHead.close();

		// we have to get the nspace for the function code calling
		int intOperator = infor.getOper();
		int nSpace = getNSpaceByOper(intOperator);

		// now we also need the string that how can we call the function
		// in the code
		fileName  = infor.getFileName();
		fileName += ".vrr_func";
		ofstream cppCode;
		cppCode.open(fileName.c_str());

		// print out function line to this file
		// also we need to complete the VRR part
		string argListInCall = sqintsPrint.transformArgList(argList);
		line = functionName + "( " + argListInCall + " );";
		printLine(nSpace,line,cppCode);
		cppCode << endl;
		sqintsPrint.printVRREnd(cppCode);
		cppCode << endl;
		cppCode.close();
	}
}

void SQInts::doRR() const
{
	vector<ShellQuartet> hrrResultSQList;
	doHRR(hrrResultSQList); 
	doVRR(hrrResultSQList); 
}

bool SQInts::isFileExist() const 
{
	bool withTmpWorkDir = false;
	string file = infor.getFileName(withTmpWorkDir);
	string fileName = file + ".cpp";

	// now look that whether the file exists?
	path cppFile(fileName.c_str());
	if (exists(cppFile)) {
		return true;
	}
	return false;
}

void SQInts::assembleWorkingCPPFile(bool isHRR) const
{
	//
	// when we do file split for large code generation
	// project, we will assemble the working cpp file pieces
	// together here
	//
	bool withTmpDir = false;
	string cppFile = infor.getFileName(withTmpDir);
	if (isHRR) {
		cppFile = cppFile + "_hrr.cpp";
	}else{
		cppFile = cppFile + "_vrr.cpp";
	}
	ofstream RR;
	RR.open(cppFile.c_str());

	// generate the head of file
	headPrinting(RR);

	// function type information
	string ext = ".vrr_head";
	if (isHRR) ext = ".hrr_head";
	string rrHeadFile = infor.getFileName();
	rrHeadFile = rrHeadFile + ext;
	path headFile(rrHeadFile.c_str());
	if (! exists(headFile)) {
		cout << "Missing head file " << rrHeadFile << endl;
		crash(true, "Why the head file does not exist???");
	}
	ifstream rrHead;
	rrHead.open(rrHeadFile.c_str(),ios::in);
	string line;
	getline(rrHead,line);
	RR << line << endl;
	rrHead.close();
	line = "{";
	printLine(0,line,RR);
	RR << endl;

	// now we are going to generate the main body of 
	// the codes, namely, we just copy the contents
	// of files
	vector<string> fileExtList;
	int nFiles = 1;
	if (isHRR) {
		nFiles = 2;
		fileExtList.reserve(nFiles);
		fileExtList.push_back(".hrr1");
		fileExtList.push_back(".hrr2");
	}else{
		fileExtList.reserve(nFiles);
		fileExtList.push_back(".vrr");
	}

	// now let's handle the working piece of file
	for(int i=0; i<nFiles; i++) {

		// get the file name
		// we note, that the HRR file may not exist
		// then we just continue
		string ext = fileExtList[i];
		string rrFile = infor.getFileName();
		rrFile = rrFile + ext;
		if (! exists(rrFile)) continue;
		ifstream rrWork;
		rrWork.open(rrFile.c_str(),ios::in);
		while(getline(rrWork,line)) {
			RR << line << endl;
		}
		rrWork.close();
	}

	// now finalize the cpp file
	RR << "}" << endl;
	RR.close();
}

void SQInts::assembleTopCPPFile() const
{
	// open the top cpp file
	bool withTmpDir = false;
	string cppFile = infor.getFileName(withTmpDir);
	cppFile = cppFile + ".cpp";
	ofstream CPP;
	CPP.open(cppFile.c_str());

	// generate the head of file
	headPrinting(CPP);

	// now work is different depending on that
	// whether we do split file or not
	// we just use HRR to see whether we have 
	// the file split requirement
	bool doFileSplit = infor.splitCPPFile();

	// if it's in file split mode, we need to 
	// include the function type information
	if (doFileSplit) {

		// include two head files
		vector<string> extList;
		extList.reserve(2);
		extList.push_back(".vrr_head");
		extList.push_back(".hrr_head");

		// now make the function head
		for(int i=0; i<2; i++) {

			// get the file name
			string ext = extList[i];
			string func = infor.getFileName();
			func = func + ext;

			// we note, that the file may not exist
			if (! exists(func)) continue;
			ifstream head;
			string line;
			head.open(func.c_str(),ios::in);
			getline(head,line);

			// now print it 
			line = "extern " + line + ";";
			CPP << line << endl;
			CPP << endl;
			head.close();
		}
	}

	// function type information
	// now all functions are in void type
	string functionName = infor.getFuncName();
	string argList = getArgList();
	string returnType = "void ";
	//if(sigCheck(infor.getOper())) returnType = "bool ";
	string line = returnType + functionName + "( " + argList + " )";
	printLine(0,line,CPP);
	line = "{";
	printLine(0,line,CPP);
	CPP << endl;

	// before VRR part, let's see that whether we need additional
	// loop etc. over the whole VRR/HRR body
	// for example, ESP etc. need a loop of grids
	// we do it here
	int oper = infor.getOper();
	if (oper == ESP) {
		line = "// loop over grid points ";
		printLine(2,line,CPP);
		line = "for(UInt iGrid=0; iGrid<nGrids; iGrid++) {";
		printLine(2,line,CPP);
		CPP << endl;
	}

	// generate the VRR and HRR part
	vector<string> fileExtList;
	fileExtList.reserve(3);
	if (doFileSplit) {
		fileExtList.push_back(".vrr_var");
		fileExtList.push_back(".vrr_func");
		fileExtList.push_back(".hrr_func");
	}else{
		fileExtList.push_back(".vrr");
		fileExtList.push_back(".hrr1");
		fileExtList.push_back(".hrr2");
	}

	// now let's handle the working piece of file
	for(int i=0; i<3; i++) {

		// get the file name
		string ext = fileExtList[i];
		string rrFile = infor.getFileName();
		rrFile = rrFile + ext;

		// we note, that the file may not exist
		if (! exists(rrFile)) continue;
		ifstream rrWork;
		rrWork.open(rrFile.c_str(),ios::in);
		while(getline(rrWork,line)) {
			CPP << line << endl;
		}
		rrWork.close();
	}

	// we do not do it anymore
	// at the HRR end we need to return true 
	// for the case of using fmt function
	//if(sigCheck(infor.getOper())) {
	//	CPP << endl;
	//	CPP << "  // for shell quartets using fmt function, we do significance check" << endl;
	//	CPP << "  return true; " << endl;
	//}

	// for ESP etc. we need an additional } to close the loop body
	if (oper == ESP) {
		CPP << "  }" << endl;
	}

	// now finalize the cpp file
	CPP << "}" << endl;
	CPP.close();
}

void SQInts::codeGeneration() const
{
	// generate the tmp folder
	// the work dir name should be same with
	// the one used in function of getProjectFileDir
	int intOperator = infor.getOper();
	int jobOrder    = infor.getJobOrder();
	string workDir = infor.getProjectTmpFileDir(intOperator,jobOrder);
	path tmpWorkDir(workDir.c_str());
	create_directory(tmpWorkDir);

	// firstly, let's finish the RR process
	// so to generate the tmp files
	// namely, vrr, hrr1 and/or hrr2
	// all of these are in the work dir
	doRR();

	// now let's assemble the vrr cpp file
	// as well as hrr cpp file if split file
	// demanded
	// first is VRR, second is HRR
	if (infor.splitCPPFile()){
		assembleWorkingCPPFile(false);
		if (infor.hasHrr()) assembleWorkingCPPFile(true);
	}

	// now assemble the top cpp file
	assembleTopCPPFile();

	// finally, clean the tmp file
	remove_all(tmpWorkDir);
}

