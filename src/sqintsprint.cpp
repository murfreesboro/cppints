//
// CPPINTS: A C++ Program to Generate Analytical Integrals Based on Gaussian
// Form Primitive Functions
//
// Copyright (C) 2012-2014 Fenglai Liu
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
#include "printing.h"
#include "sqintsinfor.h"
#include "derivinfor.h"
#include "integral.h"
#include "inttype.h"
#include "shell.h"
#include "basis.h"
#include "shellsymbol.h"
#include "boost/lexical_cast.hpp"
#include <boost/algorithm/string.hpp>   // string handling
#include "sqintsprint.h"
using namespace printing;
using namespace sqintsinfor;
using namespace inttype;
using namespace integral;
using namespace shell;
using namespace basis;
using namespace derivinfor;
using boost::lexical_cast;
using namespace boost;
using namespace sqintsprint;

//
// !!!! second print section: bottom integrals printing
// We note, that the bottom integral generation could be easy
// or complicated. We will do it here
//
void SQIntsPrint::printMOMBottomIntegrals(const string& name) const
{
	// we need to determine the maximum momentum order
	int maxMOMOrder = 0;
	const vector<ShellQuartet>& inputSQList = infor.getInputSQList();
	for(int iSQ=0; iSQ<(int)inputSQList.size(); iSQ++) {
		const Shell& s = inputSQList[iSQ].getShell(KET1);
		int L = s.getL();
		if (L>maxMOMOrder) maxMOMOrder = L;
	}
	maxMOMOrder += infor.getJobOrder();

	// what is the nspace for the code?
	// if everything is put into one cpp file,
	// it's 4
	// else it's in a vrr cpp file, then it's 2
	int nSpace = 4;
	if (infor.splitCPPFile()) nSpace = 2;

	// create the file 
	ofstream file;
	file.open(name.c_str(),std::ofstream::app);

	// based on the two body overlap, we further 
	// print more about for the mom integrals
	file << endl;
	string line = "// ";
	printLine(nSpace,line,file);
	line = "// now create the bottom integrals for momentum";
	printLine(nSpace,line,file);
	line = "// ";
	printLine(nSpace,line,file);
	line = "Double I_MOM_S_S_S   = I_TWOBODYOVERLAP_S_S;";
	printLine(nSpace,line,file);

	//
	// all of MOM bottom integrals are in the form of 
	// S_S_x^(l,m,n)
	// x is an arbitrary basis set
	// the recursive relation is given as:
	// S_S_x = PC_i*S_S_(x-li) + c*oned2z*S_S_(x-2*li)
	// we will manually code it up
	//
	for(int iOrder=1; iOrder<=maxMOMOrder; iOrder++) {
		Shell s(iOrder);
		int nBasis = s.getBasisSetNumber();
		for(int iBas=0; iBas<nBasis; iBas++) {

			// now get the lmn
			int l=-1;
			int m=-1;
			int n=-1;
			s.getBasisSetFromIndex(iBas,l,m,n);

			// now take a look that who is the smallest one?
			// l, m or n?
			// we note, that the smallest but should be at 
			// least 1
			int max1 = std::max(l,m);
			int pos  = std::max(max1,n);
			if (n<pos && n>0) pos = n;
			if (m<pos && m>0) pos = m;
			if (l<pos && l>0) pos = l;

			// now judge the direction
			string direction = "Z";
			if (pos == m) direction = "Y";
			if (pos == l) direction = "X";
			string PC = "PC" + direction;

			// name of this integral
			string momIntName = "I_MOM_S_S";
			Basis bas(l,m,n);
			string lhs = momIntName + "_" + bas.getName();

			// now it's the terms in RHS
			string rhs1;
			string rhs2;
			bool hasSecondTerm = false;
			string coef;
			if (direction == "X") {
				Basis bas1(l-1,m,n);
				rhs1 = momIntName + "_" + bas1.getName();
				coef = lexical_cast<string>(l-1);
				if (l-2>=0) {
					hasSecondTerm = true;
					Basis bas2(l-2,m,n);
					rhs2 = momIntName + "_" + bas2.getName();
				}
			}else if (direction == "Y") {
				Basis bas1(l,m-1,n);
				rhs1 = momIntName + "_" + bas1.getName();
				coef = lexical_cast<string>(m-1);
				if (m-2>=0) {
					hasSecondTerm = true;
					Basis bas2(l,m-2,n);
					rhs2 = momIntName + "_" + bas2.getName();
				}
			}else{
				Basis bas1(l,m,n-1);
				rhs1 = momIntName + "_" + bas1.getName();
				coef = lexical_cast<string>(n-1);
				if (n-2>=0) {
					hasSecondTerm = true;
					Basis bas2(l,m,n-2);
					rhs2 = momIntName + "_" + bas2.getName();
				}
			}

			// now let's compose the final integral expression
			string line = "Double ";
			line = line + lhs + " = ";
			line = line + PC + "*" + rhs1;
			if (hasSecondTerm) {
				if (coef != "1") {
					line = line + " + " + coef;
					line = line + "*oned2z*" + rhs2;
				}else{
					line = line + " + oned2z*" + rhs2;
				}
			}
			line = line + ";";
			printLine(nSpace,line,file);
		}
	}

	// finally, close the file
	file.close();
}

void SQIntsPrint::printBottomIntegrals(const string& name) const
{
	int oper = infor.getOper();
	if (oper == MOM) {
		printMOMBottomIntegrals(name);
	}else{
		crash(true, "In printBottomIntegrals the operator is not supported");
	}
}

//
// !!!! third print section: the head of the function
// arguments list etc.
// they are very useful when we use file split
//
string SQIntsPrint::getVRRArgList() const
{

	///////////////////////////////////////////////////////
	//             get the status of the VRR             //
	///////////////////////////////////////////////////////
	
	bool hasRROnBRA1 = hasVRROnVar(BRA1);
	bool hasRROnBRA2 = hasVRROnVar(BRA2);
	bool hasRROnKET1 = hasVRROnVar(KET1);
	bool hasRROnKET2 = hasVRROnVar(KET2);
	bool hasRROnBRA  = hasRROnBRA1 || hasRROnBRA2;
	bool hasRROnKET  = hasRROnKET1 || hasRROnKET2;
	int  maxLSum     = getMaxLSum();
	int  intOperator = infor.getOper();

	///////////////////////////////////////////////////////
	//                now get the arg list               //
	//                firstly it's input                 //
	//                      variables                    //
	///////////////////////////////////////////////////////
	string arg;

	// coefficients is the first part 
	// we note, it's only that this is composite shell quartet,
	// we need to bring the coefficients into the subroutine
	// so that to form contraction in the end stage of VRR
	// else the contraction of coefficients is done at the 
	// beginning
	if (infor.isComSQ()) {
		for(int i=0; i<2; i++) {

			// get the coe length
			int nCoeArray = 0;
			if (i == 0) {
				nCoeArray = infor.getCoeArrayLength(BRA);
			}else{
				nCoeArray = infor.getCoeArrayLength(KET);
			}

			// now print 
			if (nCoeArray >= 1) {
				if (i == 0) {
					arg = arg + "const Double& ic2, ";
					for(int p=1; p<nCoeArray; p++) {
						string lhs  = "ic2_" + lexical_cast<string>(p);
						string line = "const Double& " + lhs + ", ";
						arg  = arg + line;
					}
				}else{
					arg = arg + "const Double& jc2, ";
					for(int p=1; p<nCoeArray; p++) {
						string lhs  = "jc2_" + lexical_cast<string>(p);
						string line = "const Double& " + lhs + ", ";
						arg  = arg + line;
					}
				}
			}
		}
	}

	// now let's deal with the variables used in VRR
	// two body overlap we note, that mom integrals 
	// is also put into the VRR body part
	if (intOperator == TWOBODYOVERLAP || intOperator == MOM) {

		// VRR variable
		arg = arg + "const Double& oned2z, ";
		if (hasRROnBRA1) {
			arg = arg + "const Double& PAX, const Double& PAY, const Double& PAZ, ";
		}
		if (hasRROnBRA2) {
			arg = arg + "const Double& PBX, const Double& PBY, const Double& PBZ, ";
		}

		// additional var for MOM
		// we note that this is always needed for mom integrals
		if (intOperator == MOM) {
			arg = arg + "const Double& PCX, const Double& PCY, const Double& PCZ, ";
		}

		// bottom integral
		arg = arg + "const Double& I_TWOBODYOVERLAP_S_S, ";
	}

	// kinetic integrals
	if (intOperator == KINETIC) {

		// VRR variable
		arg = arg + "const Double& oned2z, const Double& adz, const Double& bdz, "
				"const Double& twoxi, ";
		if (hasRROnBRA1) {
			arg = arg + "const Double& PAX, const Double& PAY, const Double& PAZ, ";
		}
		if (hasRROnBRA2) {
			arg = arg + "const Double& PBX, const Double& PBY, const Double& PBZ, ";
		}

		// bottom integral
		arg = arg + "const Double& I_TWOBODYOVERLAP_S_S, const Double& I_KINETIC_S_S, ";
	}

	// three body overlap
	if (intOperator == THREEBODYOVERLAP) {

		// VRR variable
		arg = arg + "const Double& oned2zeta, ";
		if (hasRROnBRA1) {
			arg = arg + "const Double& GAX, const Double& GAY, const Double& GAZ, ";
		}
		if (hasRROnBRA2) {
			arg = arg + "const Double& GBX, const Double& GBY, const Double& GBZ, ";
		}
		if (hasRROnKET1) {
			arg = arg + "const Double& GCX, const Double& GCY, const Double& GCZ, ";
		}

		// bottom integral
		arg = arg + "const Double& I_THREEBODYOVERLAP_S_S_S, ";
	}

	// three body kinetic energy
	if (intOperator == THREEBODYKI) {

		// VRR variable
		// we need to add in all of three exp factors
		// alpha beta nad gamma
		// also, since the passed in sqlist is empty
		// therefore we will add in
		arg = arg + "const Double& oned2zeta, ";
		arg = arg + "const Double& alpha, const Double& beta, const Double& gamma, ";
		arg = arg + "const Double& GAX, const Double& GAY, const Double& GAZ, ";
		arg = arg + "const Double& GBX, const Double& GBY, const Double& GBZ, ";
		arg = arg + "const Double& GCX, const Double& GCY, const Double& GCZ, ";

		// bottom integral
		arg = arg + "const Double& I_THREEBODYOVERLAP_S_S_S, ";
	}

	// nai
	if (intOperator == NAI) {

		// VRR variable
		arg = arg + "const Double& oned2z, ";
		arg = arg + "const Double& PNX, const Double& PNY, const Double& PNZ, ";
		if (hasRROnBRA1) {
			arg = arg + "const Double& PAX, const Double& PAY, const Double& PAZ, ";
		}
		if (hasRROnBRA2) {
			arg = arg + "const Double& PBX, const Double& PBY, const Double& PBZ, ";
		}

		// bottom integral
		for(int m=0; m<=maxLSum; m++) {
			string name = getBottomIntName(m,intOperator);
			arg = arg + "const Double& " + name + ", ";
		}
	}

	// esp
	if (intOperator == ESP) {

		// VRR variable
		arg = arg + "const Double& oned2z, ";
		arg = arg + "const Double& PRX, const Double& PRY, const Double& PRZ, ";
		if (hasRROnBRA1) {
			arg = arg + "const Double& PAX, const Double& PAY, const Double& PAZ, ";
		}
		if (hasRROnBRA2) {
			arg = arg + "const Double& PBX, const Double& PBY, const Double& PBZ, ";
		}

		// bottom integral
		for(int m=0; m<=maxLSum; m++) {
			string name = getBottomIntName(m,intOperator);
			arg = arg + "const Double& " + name + ", ";
		}
	}

	// eri
	if (intOperator == ERI) {

		// VRR variable
		arg = arg + "const Double& oned2k, ";

		if (hasRROnBRA) {
			arg = arg + "const Double& oned2z, const Double& rhod2zsq, ";
			arg = arg + "const Double& WPX, const Double& WPY, const Double& WPZ, ";
		}
		if (hasRROnKET) {
			arg = arg + "const Double& oned2e, const Double& rhod2esq, ";
			arg = arg + "const Double& WQX, const Double& WQY, const Double& WQZ, ";
		}
		if (hasRROnBRA1) {
			arg = arg + "const Double& PAX, const Double& PAY, const Double& PAZ, ";
		}
		if (hasRROnBRA2) {
			arg = arg + "const Double& PBX, const Double& PBY, const Double& PBZ, ";
		}
		if (hasRROnKET1) {
			arg = arg + "const Double& QCX, const Double& QCY, const Double& QCZ, ";
		}
		if (hasRROnKET2) {
			arg = arg + "const Double& QDX, const Double& QDY, const Double& QDZ, ";
		}

		// bottom integral
		for(int m=0; m<=maxLSum; m++) {
			string name = getBottomIntName(m,intOperator);
			arg = arg + "const Double& " + name + ", ";
		}
	}

	///////////////////////////////////////////////////////
	//                now get the arg list               //
	//                secondly it's output               //
	//                      variables                    //
	///////////////////////////////////////////////////////
	bool hasResultSQ = false;
	for(int iSQ=0; iSQ<(int)rrSQList.size(); iSQ++) {

		// here we have to test that whether the given sq is real result
		const ShellQuartet& sq = rrSQList[iSQ];
		if (infor.isResult(sq)) {
			hasResultSQ = true;
			continue;
		}

		// now include it in the argument list
		string line = "Double* " + sq.getName();
		line = line + ", ";
		arg = arg + line;
	}

	// if the input shell quartet contains the final result,
	// we need to include the abcd array
	if (hasResultSQ) {
		arg = arg + "Double* abcd";
	}

	// add in the localmemscr use
	int vrrMethod = infor.getVRRMethod();
	bool vrrUseArray = infor.withArrayIndex(vrrMethod);
	if (vrrUseArray && infor.useSCRVec()) {
		if (hasResultSQ) {
			arg = arg + ", LocalMemScr& scr";
		}else{
			arg = arg + "LocalMemScr& scr";
		}
	}

	// finally we have to check the ","
	// if there's no result abcd, we will have 
	// arglist with additional ,
	for(int pos=arg.size()-1; pos>=0; pos--) {
		char c = arg[pos];
		
		// for digit or word, then this is the variable name
		// we stop here
		if (isalnum(c)) break;

		// we change the additional ' to space
		if (c == ',') arg[pos] = ' ';
	}
	return arg;
}

string SQIntsPrint::getHRRArgList() const
{
	// firstly, push in the A B centers etc.
	string arg;
	int firstSide  = infor.get1stSide();
	int secondSide = infor.get2edSide();
	bool hasBra = false;
	bool hasKet = false;
	if (firstSide == BRA  || secondSide == BRA)  hasBra = true;
	if (firstSide == BRA1 || secondSide == BRA1) hasBra = true;
	if (firstSide == BRA2 || secondSide == BRA2) hasBra = true;
	if (firstSide == KET  || secondSide == KET)  hasKet = true;
	if (firstSide == KET1 || secondSide == KET1) hasKet = true;
	if (firstSide == KET2 || secondSide == KET2) hasKet = true;
	if (hasBra) arg = arg + "const Double* A, const Double* B, ";
	if (hasKet) arg = arg + "const Double* C, const Double* D, ";

	// secondly, push in the input variables
	// we note, that the list also contains the 
	// shell quartet that could be finished in VRR step
	// we will check it here
	for(int iSQ=0; iSQ<(int)rrSQList.size(); iSQ++) {

		// here we have to test that whether the given sq is real result
		const ShellQuartet& sq = rrSQList[iSQ];
		if (infor.isResult(sq)) continue;

		// now push it in
		arg  = arg + "const Double* " + sq.getName() + ", ";
	}

	// finally, the output vectors
	// which has aunique name of "abcd"
	arg  = arg + "Double* abcd";

	// add in the localmemscr use
	bool hrrUseArray = infor.withArrayIndex(HRR);
	if (hrrUseArray && infor.useSCRVec()) {
		arg = arg + ", LocalMemScr& scr";
	}

	return arg;
}

string SQIntsPrint::transformArgList(string arg) const
{
	// the input arg is the arg get for VRR/HRR function
	// we need to drop the const Double& etc. so that 
	// to use it in the place to call it
	LineParse lp(arg);
	string result;
	for(int i=0; i<lp.getNPieces(); i++) {
		
		// we begin to determine that whether this is variable name
		string val = lp.findValue(i);
		if (val.find("const")!=std::string::npos) continue;
		if (val.find("vector")!=std::string::npos) continue;
		if (val.find("Double")!=std::string::npos) continue;
		if (val.find("LocalMemScr")!=std::string::npos) continue;

		// if the name is turned out to be a shell quartet name,
		// then this is must be VRR's output result or HRR input
		// result array
		// since we all use pointer, we need to pass in pointer inside
		if (val.find("SQ")!=std::string::npos) {

			// see that whether the name has ' additional to the name or not
			int len = val.length();
			if (val[len-1] == ',') {
				val[len-1] = '[';
				val.append("0],");
				val = "&" + val;
			}else{
				val = "&" + val + "[0]";
			}
		};

		// now this should be good enough
		// the value should be the variable name now
		// also the comma should be also within it
		result = result + val;
	}
	return result;
}

//
// !!!! fourth print section prepare the VRR process:
// this part is to print out the VRR head for future VRR process
// 1  vector declaration; 
// 2  print out comments for variables;
// 3  initial the VRR variables;
// 4  generate the SSSS^{m} integrals if necessary
//
bool SQIntsPrint::hasVRROnVar(const int& var) const
{
	for(int iSQ=0; iSQ<(int)rrSQList.size(); iSQ++) {
		if (rrSQList[iSQ].isNonSShell(var)) return true;
	}
	return false;
}

int SQIntsPrint::nTotalInts() const
{
	int nInts = 0;
	for(int iSQ=0; iSQ<(int)rrSQList.size(); iSQ++) {
		nInts += rrSQList[iSQ].getNInts();
	}
	return nInts;
}

int SQIntsPrint::getMaxLSum() const
{
	int maxL = 0;
	for(int iSQ=0; iSQ<(int)rrSQList.size(); iSQ++) {
		int LSum = rrSQList[iSQ].getLSum();
		if (LSum > maxL) maxL = LSum;
	}
	return maxL+infor.getJobOrder();
}

string SQIntsPrint::getBottomIntName(const int& m, const int& oper) const
{
	//
	// Here we note that we should count in all of 
	// integral bodies 
	// therefore we use the rr order
	//
	// on the other hand, such name is used only in
	// VRR mode. so we do not have any division infor
	// the division infor will be destroyed for any
	// sq passed into vrr
	//
	string name = "I";
	string oper_name = getOperStringName(oper);
	name = name + "_" + oper_name;
	int nBody = getRROrder(oper);
	if (nBody == 1) {
		name = name + "_S";
	}else if (nBody == 2) {
		name = name + "_S_S";
	}else if (nBody == 3) {
		name = name + "_S_S_S";
	}else {
		name = name + "_S_S_S_S";
	}
	if (m>0) {
		name = name + "_M" + lexical_cast<string>(m);
	}
	return name;
}

void SQIntsPrint::fmtIntegralsGeneration(const int& maxLSum, 
		const int& oper, const int& nSpace, ofstream& file) const
{
	// get the information from infor
	int m_limit = infor.M_limit;
	int fmt_error = infor.fmt_error;

	//
	// to generate the (SS|SS)^{m} in terms of f_{m}(t) function
	// maxM=0, we just use erf function
	// if maxM<=M_limit and maxM>0, then mainly we use up recursive relation
	// else we will calculate fm directly and use down recursive relation
	// see the doc for more information
	//
	// We note, that the code below applies for both ERI and NAI etc. situations
	// where the (SSSS)^{m} could be expressed as: O*rho^{0.5}*fm(t)
	//
	if (maxLSum==0) {

		// 
		// for (SS|SS)^{0}
		// we need to consider the case that u == 0
		//
		string name = getBottomIntName(0,oper);
		string line = "Double " + name + " = 0.0E0;"; 
		printLine(nSpace,line,file);
		line = "if (fabs(u)<THRESHOLD_MATH) {";
		printLine(nSpace,line,file);
		line = name + " = prefactor*sqrho*TWOOVERSQRTPI;";
		printLine(nSpace+2,line,file);
		line = "}else{";
		printLine(nSpace,line,file);
		line = "#ifdef WITH_SINGLE_PRECISION";
		printLine(0,line,file);
		line = name + " = (prefactor*sqrho/squ)*erf(squ);";
		printLine(nSpace+2,line,file);
		line = "#else ";
		printLine(0,line,file);
		line = name + " = (prefactor*sqrho/squ)*erfVal;";
		printLine(nSpace+2,line,file);
		line = "#endif";
		printLine(0,line,file);
		line = "}";
		printLine(nSpace,line,file);
		file << endl;

	}else {

		string mlimit   = lexical_cast<string>(m_limit);
		string fmterror = lexical_cast<string>(fmt_error);

		// comments
		file << endl;
		string line ="//";
		printLine(nSpace,line,file);
		line = "//";
		printLine(nSpace,line,file);
		line = "// now here for maxM>0 to compute the infamous incomplete Gamma function f_{m}(u)";
		printLine(nSpace,line,file);
		line = "// the implementation is divided in two situations:";
		printLine(nSpace,line,file);
		line = "// 1  if u <=1.8; use power series to get f_{Mmax}(u), then use down recursive"; 
		printLine(nSpace,line,file);
		line = "//    relation to get the rest of incomplete Gamma functions;";
		printLine(nSpace,line,file);
		line = "// 2  for u >1.8 and M <= " + mlimit + " we calculate erf(u), then use up recursive";
		printLine(nSpace,line,file);
		line = "//    relation to calculate the rest of results";
		printLine(nSpace,line,file);
		line = "// 3  for u> 1.8 and M >  " + mlimit + " we calculate f_{Mmax}(u) then use down ";
		printLine(nSpace,line,file);
		line = "//    recursive relation to get rest of incomplete Gamma functions ";
		printLine(nSpace,line,file);
		line = "// The above procedure is tested for u between 0 to 40 with step length 1.0E-6"; 
		printLine(nSpace,line,file);
		line = "// (or 1.0E-5 for float double data), for up recursive relation it shows the error"; 
		printLine(nSpace,line,file);
		line = "// within 1.0E-12 (for M_limit = 12 or error within 1.0E-6 for float type of data";
		printLine(nSpace,line,file);
		line = "// For the polynomial expansion and down recursive procedure the error is within ";
		printLine(nSpace,line,file);
		line = "// 1.0E-14. All of the testing details please refer to the fmt_test folder";
		printLine(nSpace,line,file);
		line = "// ";
		printLine(nSpace,line,file);
		line = "// There's one thing need to note for up recursive process. We found that the up";
		printLine(nSpace,line,file);
		line = "// recursive procedure is only stable for maxM<=" + mlimit + " and u>1.8 with double";
		printLine(nSpace,line,file);
		line = "// precision data, single precision data will lose accuracy quickly so the result";
		printLine(nSpace,line,file);
		line = "// for single precision calculation is not doable. Therefore if the \"WITH_SINGLE_PRECISION\"";
		printLine(nSpace,line,file);
		line = "// is defined, then for erf function calculation as well as up recursive";
		printLine(nSpace,line,file);
		line = "// process we will use the double type of data";
		printLine(nSpace,line,file);
		line = "// ";
		printLine(nSpace,line,file);
		line = "//";
		printLine(nSpace,line,file);
		file << endl;

		// list the name of the integral
		for(int m=0; m<=maxLSum; m++) {
			string name = getBottomIntName(m,oper);
			line = "Double " + name + "  = 0.0E0;"; 
			printLine(nSpace,line,file);
		}
		file << endl;

		if (maxLSum<=m_limit) {

			// we note, that for the float type of variable
			// we need the double type of variable to 
			// store the accuracy
			line = "// declare double type of results to store the accuracy";
			printLine(nSpace,line,file);
			line = "#ifdef WITH_SINGLE_PRECISION";
			printLine(0,line,file);
			for(int m=0; m<=maxLSum; m++) {
				string name = getBottomIntName(m,oper);
				name = name + "_d";
				line = "double " + name + "  = 0.0E0;"; 
				printLine(nSpace,line,file);
			}
			line = "#endif";
			printLine(0,line,file);
			file << endl;
		}

		// now calculating the (SS|SS)^{Mmax} in terms of power series
		// the power series is 18 terms
		// however, nTerms is 17 since the index starts from 0, then to 17
		// totally 18 terms
		int nTerms = 17;
		line="if (u<=1.8E0) {";
		printLine(nSpace,line,file);
		file << endl;
		line = "// calculate (SS|SS)^{Mmax}";
		printLine(nSpace+2,line,file);
		line = "// use 18 terms power series to expand the (SS|SS)^{Mmax}";
		printLine(nSpace+2,line,file);
		line = "Double u2 = 2.0E0*u;";
		printLine(nSpace+2,line,file);
		string coe = "ONEOVER" + lexical_cast<string>(2*(maxLSum+nTerms)+1);
		string name= getBottomIntName(maxLSum,oper);
		line = name + " = " + "1.0E0+u2*" + coe + ";";
		printLine(nSpace+2,line,file);
		for(int m=maxLSum+nTerms-1; m>maxLSum; m--) {
			coe  = "ONEOVER" + lexical_cast<string>(2*m+1);
			line = name + " = " + "1.0E0+u2*" + coe + "*" + name + ";";
			printLine(nSpace+2,line,file);
		}
		coe  = "ONEOVER" + lexical_cast<string>(2*maxLSum+1);
		line = name + " = " + coe + "*" + name + ";";
		printLine(nSpace+2,line,file);
		line = "Double eu = exp(-u);";
		printLine(nSpace+2,line,file);
		line = "Double f  = TWOOVERSQRTPI*prefactor*sqrho*eu;";
		printLine(nSpace+2,line,file);
		line = name + "  = " + "f*" + name + ";";
		printLine(nSpace+2,line,file);
		file << endl;

		// do down recursive relation
		line = "// now use down recursive relation to get"; 
		printLine(nSpace+2,line,file);
		line = "// rest of (SS|SS)^{m}";
		printLine(nSpace+2,line,file);
		for(int m=maxLSum-1; m>=0; m--) {
			coe  = "ONEOVER" + lexical_cast<string>(2*m+1);
			name= getBottomIntName(m,oper);
			string pInt = getBottomIntName(m+1,oper); 
			line = name + "  = " + coe + "*(u2*" + pInt + "+f);";
			printLine(nSpace+2,line,file);
		}
		file << endl;

		// now calculating the (SS|SS)^{0} for u > 1.8
		line="}else{";
		printLine(nSpace,line,file);

		if (maxLSum<=m_limit) {

			// now begin the float var part
			line = "#ifdef WITH_SINGLE_PRECISION";
			printLine(0,line,file);
			file << endl;

			// SSSS variables
			line = "// recompute the variable in terms of double accuracy"; 
			printLine(nSpace+2,line,file);
			line = "double u_d     = u;";
			printLine(nSpace+2,line,file);
			line = "double rho_d   = rho;";
			printLine(nSpace+2,line,file);
			line = "double fac_d   = prefactor;";
			printLine(nSpace+2,line,file);
			line = "double sqrho_d = sqrt(rho_d);";
			printLine(nSpace+2,line,file);
			line = "double squ_d   = sqrt(u_d);";
			printLine(nSpace+2,line,file);
			file << endl;

			line = "// use erf function to get (SS|SS)^{0}"; 
			printLine(nSpace+2,line,file);
			name = getBottomIntName(0,oper);
			name = name + "_d";
			line = "if (fabs(u_d)<THRESHOLD_MATH) {";
			printLine(nSpace+2,line,file);
			line = name + " = fac_d*sqrho_d*TWOOVERSQRTPI;";
			printLine(nSpace+4,line,file);
			line = "}else{";
			printLine(nSpace+2,line,file);
			line = name + " = (fac_d*sqrho_d/squ_d)*erf(squ_d);";
			printLine(nSpace+4,line,file);
			line = "}";
			printLine(nSpace+2,line,file);
			file << endl;

			// do up recursive relation
			line = "// now use up recursive relation to get"; 
			printLine(nSpace+2,line,file);
			line = "// rest of (SS|SS)^{m}";
			printLine(nSpace+2,line,file);
			line = "double oneO2u  = 0.5E0/u_d;";
			printLine(nSpace+2,line,file);
			line = "double eu      = exp(-u_d);";
			printLine(nSpace+2,line,file);
			line = "double f       = TWOOVERSQRTPI*fac_d*sqrho_d*eu;";
			printLine(nSpace+2,line,file);
			for(int m=1; m<=maxLSum; m++) {
				coe  = lexical_cast<string>(2*m-1) + ".0E0";
				name = getBottomIntName(m,oper);
				name = name + "_d";
				string pInt = getBottomIntName(m-1,oper); 
				pInt = pInt + "_d";
				line = name + " = oneO2u*(" + coe + "*" + pInt + "-f);";
				printLine(nSpace+2,line,file);
			}
			file << endl;

			// now let's write the result back
			line = "// write the double result back to the float var";
			printLine(nSpace+2,line,file);
			for(int m=0; m<=maxLSum; m++) {
				string name = getBottomIntName(m,oper);
				string name_d = name + "_d";
				line = name + " = static_cast<Double>(" + name_d + ");"; 
				printLine(nSpace+2,line,file);
			}
			file << endl;

			// now switch to the double part of codes
			line = "#else";
			printLine(0,line,file);
			file << endl;

			line = "// use erf function to get (SS|SS)^{0}"; 
			printLine(nSpace+2,line,file);
			name = getBottomIntName(0,oper);
			line = "if (fabs(u)<THRESHOLD_MATH) {";
			printLine(nSpace+2,line,file);
			line = name + " = prefactor*sqrho*TWOOVERSQRTPI;";
			printLine(nSpace+4,line,file);
			line = "}else{";
			printLine(nSpace+2,line,file);
			line = name + " = (prefactor*sqrho/squ)*erfVal;";
			printLine(nSpace+4,line,file);
			line = "}";
			printLine(nSpace+2,line,file);
			file << endl;

			// do up recursive relation
			line = "// now use up recursive relation to get"; 
			printLine(nSpace+2,line,file);
			line = "// rest of (SS|SS)^{m}";
			printLine(nSpace+2,line,file);
			line = "Double oneO2u = 0.5E0/u;";
			printLine(nSpace+2,line,file);
			line = "Double eu     = exp(-u);";
			printLine(nSpace+2,line,file);
			line = "Double f      = TWOOVERSQRTPI*prefactor*sqrho*eu;";
			printLine(nSpace+2,line,file);
			for(int m=1; m<=maxLSum; m++) {
				coe  = lexical_cast<string>(2*m-1) + ".0E0";
				name = getBottomIntName(m,oper);
				string pInt = getBottomIntName(m-1,oper); 
				line = name + " = oneO2u*(" + coe + "*" + pInt + "-f);";
				printLine(nSpace+2,line,file);
			}
			file << endl;

			// this should be the end of maxM<=10 and u>1.1
			line = "#endif";
			printLine(0,line,file);
			file << endl;

		}else{

			// comments
			file << endl;
			string line ="//";
			printLine(nSpace+2,line,file);
			line = "// now here for maxM>M_limit";
			printLine(nSpace+2,line,file);
			line = "// use external function to calculate f_{Mmax}(t)"; 
			printLine(nSpace+2,line,file);
			line = "// then use down recursive relation to get others";
			printLine(nSpace+2,line,file);
			line = "//";
			printLine(nSpace+2,line,file);

			// now calculating the (SS|SS)^{Mmax} with incomplete gamma function
			file << endl;
			line = "// calculate (SS|SS)^{Mmax} with incomplete gamma function";
			printLine(nSpace+2,line,file);
			line = "// currently we use boost library for calculation";
			printLine(nSpace+2,line,file);
			string name= getBottomIntName(maxLSum,oper);
			string maxOrder = lexical_cast<string>(maxLSum);
			line = "if (fabs(u)<THRESHOLD_MATH) {";
			printLine(nSpace+2,line,file);
			line = name + " = 1.0E0/(2.0E0*" + maxOrder +"+1.0E0);";
			printLine(nSpace+4,line,file);
			line = "}else{";
			printLine(nSpace+2,line,file);
			line = name + " = (0.5E0/squ)*boost::math::tgamma_lower(" + maxOrder + "+0.5E0,u);";
			printLine(nSpace+4,line,file);
			line = "Double oneOveru = 1.0E0/u;";
			printLine(nSpace+4,line,file);
			line = "for(UInt i=0; i<" + maxOrder +"; i++) {";
			printLine(nSpace+4,line,file);
			line = name + " = " + name + "*oneOveru;";
			printLine(nSpace+6,line,file);
			line = "}";
			printLine(nSpace+4,line,file);
			line = "}";
			printLine(nSpace+2,line,file);
			line = "Double f = TWOOVERSQRTPI*prefactor*sqrho;";
			printLine(nSpace+2,line,file);
			line = name + " *= f;"; 
			printLine(nSpace+2,line,file);

			// do down recursive relation
			file << endl;
			line = "// now use down recursive relation to get"; 
			printLine(nSpace+2,line,file);
			line = "// rest of (SS|SS)^{m}";
			printLine(nSpace+2,line,file);
			line = "Double u2 = 2.0E0*u;";
			printLine(nSpace+2,line,file);
			line = "Double eu = exp(-u);";
			printLine(nSpace+2,line,file);
			line = "f = f*eu;";
			printLine(nSpace+2,line,file);
			for(int m=maxLSum-1; m>=0; m--) {
				string coe  = "ONEOVER" + lexical_cast<string>(2*m+1);
				name= getBottomIntName(m,oper);
				string pInt = getBottomIntName(m+1,oper); 
				line = name + "  = " + coe + "*(u2*" + pInt + "+f);";
				printLine(nSpace+2,line,file);
			}
		}

		// now finish the code section for M > 0
		line="}";
		printLine(nSpace,line,file);
	}
	file << endl;
}

void SQIntsPrint::fmtIntegralsTest(const int& maxLSum, 
		const int& oper, const int& nSpace, ofstream& file) const
{
	// for L == 0 we do not repeat it
	if (maxLSum == 0) return;

	// for other case
	file << endl;
	string line = "// ";
	printLine(nSpace,line,file);
	line = "// we use (SS|Oper|SS)^{0} to testify the magnitude order of integrals";
	printLine(nSpace,line,file);
	line = "// for the given integral based on Gaussian primitive functions.";
	printLine(nSpace,line,file);
	line = "// Because (SS|Oper|SS)^{0} and (SS|Oper|SS)^{m}(m>0) only differs with";
	printLine(nSpace,line,file);
	line = "// f_{m}(t) function, and in most of the cases f_{m}(t) > f_{m+1}(t); ";
	printLine(nSpace,line,file);
	line = "// hence we only test the case of (SS|Oper|SS)^{0} ";
	printLine(nSpace,line,file);
	line = "// ";
	printLine(nSpace,line,file);
	line = "// question about coefficients:";
	printLine(nSpace,line,file);
	line = "// the coefficient already contains normalization factor for the given L.";
	printLine(nSpace,line,file);
	line = "// In general, for a given primitive function; if it's exponent <1; as L";
	printLine(nSpace,line,file);
	line = "// goes larger, the normalization factor is smaller. If exponent > 1, ";
	printLine(nSpace,line,file);
	line = "// as L goes larger, normalization factor becomes larger. This trend";
	printLine(nSpace,line,file);
	line = "// in general complys with the understanding for normalization.";
	printLine(nSpace,line,file);
	line = "// On the other hand, in the recursive generation process, the following ";
	printLine(nSpace,line,file);
	line = "// integrals becomes larger and larger no doubt. Therefore from numerical";
	printLine(nSpace,line,file);
	line = "// point of view, if the coefficients > 1 then we can not judge the ";
	printLine(nSpace,line,file);
	line = "// integrals based on (SS|Oper|SS)^{0}, since normalization can not scale ";
	printLine(nSpace,line,file);
	line = "// it back. Therefore, if coefficients > 1 we do not do significance test";
	printLine(nSpace,line,file);
	line = "// ";
	printLine(nSpace,line,file);
	line = "#ifndef WITH_SINGLE_PRECISION";
	printLine(0,line,file);
	line = "if (fabs(ic2*jc2)<1.0E0) {"; 
	if (oper == NAI || oper == ESP) {
		line = "if (fabs(ic2)<1.0E0) {"; 
	}
	printLine(nSpace,line,file);
	string name = getBottomIntName(0,oper);
	name = name + "_IntegralTest";
	line = "Double " + name + " = 0.0E0;"; 
	printLine(nSpace+2,line,file);
	line = "if (fabs(u)<THRESHOLD_MATH) {";
	printLine(nSpace+2,line,file);
	line = name + " = pref*sqrho*TWOOVERSQRTPI;";
	printLine(nSpace+4,line,file);
	line = "}else{";
	printLine(nSpace+2,line,file);
	line = name + " = (pref*sqrho/squ)*erfVal;";
	printLine(nSpace+4,line,file);
	line = "}";
	printLine(nSpace+2,line,file);
	if (oper == ERI) {
		line = "Double prim2Thresh = thresh/(inp2*jnp2);";
	}else if (oper == NAI || oper == ESP) {
		line = "Double prim2Thresh = thresh/inp2;";
	}
	printLine(nSpace+2,line,file);

	// here we need to know the limit of fmt error
	int fmt_error = infor.fmt_error;
	string fmterror = "1.0E-12";
	if (fmt_error != 12) {
		if (fmt_error == 13) {
			fmterror = "1.0E-13";
		}else{
			crash(true,"the fmt error is not supported here in the fmtIntegralsTest function");
		}
	}
	line = "Double thresh_integralTest = prim2Thresh > " + fmterror + " ? prim2Thresh : " + fmterror + ";";
	printLine(nSpace+2,line,file);

	// here it's the code to judge the significance 
	line = "if(fabs(" + name + ")<thresh_integralTest) continue;";
	printLine(nSpace+2,line,file);

	// ok, fmt test end here
	line = "}";
	printLine(nSpace,line,file);
	line = "#endif";
	printLine(0,line,file);

	// this boolean variable is to see whether we will
	// do the following HRR part
	line = "isSignificant = true;";
	printLine(nSpace,line,file);
	file << endl;
}

void SQIntsPrint::printVRRHead(const string& name) const
{
	// create the file 
	ofstream file;
	file.open(name.c_str(),std::ofstream::app);

	// now is the vrr result 
	// here the sqlist refers to the VRR's input
	// must be the input rrSQList
	if (rrSQList.size() > 0) {
		vrrResultStatement(file);
	}

	// set up significance test if the file employs fmt function
	// however, if it's all bottom integrals we do not need to do it
	// just return true 
	if(sigCheck(infor.getOper()) && ! infor.areAllBottomSQ()) {
		string line = "// initialize the significance check for VRR part ";
		printLine(2,line,file);
		line = "bool isSignificant = false;";
		printLine(2,line,file);
		file << endl;
	}

	// now go to the given function
	int oper = infor.getOper();
	switch(oper) {
		case TWOBODYOVERLAP:
			printTwoBodyOverlapHead(file);
			break;
		case THREEBODYOVERLAP:
			printThreeBodyOverlapHead(file);
			break;
		case MOM:
			printMOMHead(file);
			break;
		case NAI:
			printNAIHead(file);
			break;
		case ESP:
			printESPHead(file);
			break;
		case ERI:
			printERIHead(file);
			break;
		case KINETIC:
			printKineticHead(file);
			break;
		case THREEBODYKI:
			printThreeBodyKIHead(file);
			break;
		default:
			crash(true,"Illegal operator type in print head function");
			break;
	}

	// now close the file
	file.close();

	//
	// do we need to calculate the bottom integrals in a complicated way?
	// this may be demanded for the operator itself
	// however, only if we do not split file
	//
	if (! infor.splitCPPFile() && complicatedBottomIntegrals(oper)) {
		printBottomIntegrals(name);
	}

	//
	// if the input shell quartets are all bottom type 
	// of shell quartets, then we just print them here
	//
	if (infor.areAllBottomSQ()) {

		//
		// re-open the stream
		//
		ofstream file;
		file.open(name.c_str(),std::ofstream::app);

		// retreive the operator infor
		// so that we set the number of space
		int nSpace = getNSpaceByOper(oper);

		//
		// this function is suitable for the situation that 
		// the sqlist (which should be the VRR part input),
		// are all S type of integrals.
		// here the situations are divided into several cases
		// case 1: the bottom integral are just pure S
		// for this case, the sq only has one integral
		// case 2: there will be a group of bottom integrals
		// for example, the MOM integrals
		// case 3:  the result abcd may have additional offset
		// other than caused by the situation that number of 
		// sq > 1. 
		//
		// in case that we may have multiple sq in the 
		// input shell quartet list, we have offset
		// to get the correct ressult index for abcd
		//

		// detect that wether we have compilcated offset for result?
		// we note that if it's all bottom sq, then the rrsqlist 
		// is same with the sqlist got from infor class
		// therefore it's safe to call ntotalInts function
		string additionalOffset;
		bool hasAdditionalOffset = resultIntegralHasAdditionalOffset(oper);
		if (hasAdditionalOffset) {
			int nInts = nTotalInts();
			additionalOffset = determineAdditionalOffset(oper,nInts);
		}

		// now do the contraction here
		int offset = 0;
		const vector<ShellQuartet>& sqlist = infor.getInputSQList();
		for(int iSQ=0; iSQ<(int)sqlist.size(); iSQ++) {
			const ShellQuartet& sq = sqlist[iSQ];
			int nInts = sq.getNInts();
			for(int i=0; i<nInts; i++) {
				Integral I(sq,i);
				string name = I.getName();
				string line;
				if (hasAdditionalOffset) {
					line = "abcd[" + additionalOffset + "+" + 
						lexical_cast<string>(i+offset) += "] += " + name + ";";
				}else{
					line = "abcd[" + lexical_cast<string>(i+offset) += "] += " + name + ";";
				}
				printLine(nSpace,line,file);
			}
			offset += nInts;
		}

		// in this case, let's close the VRR too
		printVRREnd(file);

		// now thing is done
		file.close();
	}
}

void SQIntsPrint::vrrResultStatement(ofstream& myfile) const 
{
	// set the nSpace
	// consider the additional offset for ESP etc.
	int nSpace = 2;
	if (resultIntegralHasAdditionalOffset(infor.getOper())) {
		nSpace += 2;
	}

	//
	// this is the variable declare head
	//
	string line = "//";
	printLine(nSpace,line,myfile);
	line = "// declare the variables as result of VRR process";
	printLine(nSpace,line,myfile);
	line = "//";
	printLine(nSpace,line,myfile);

	//
	// now print out the input rr sq
	//
	for(int iSQ=0; iSQ<(int)rrSQList.size(); iSQ++) {

		// the first thing is that we need to judge whether this 
		// sq is the result
		const ShellQuartet& sq = rrSQList[iSQ];
		if (infor.isResult(sq)) continue;

		// now print the VRR's input
		bool isHRRwithArrayIndex = infor.withArrayIndex(HRR);
		if (isHRRwithArrayIndex) {
			string name = sq.getName();
			int nInts   = sq.getNInts();
			string arrayType = infor.getArrayType();
			string declare   = infor.getArrayDeclare(lexical_cast<string>(nInts));
			string line      = arrayType + name + declare;
			printLine(nSpace,line,myfile);
		}else{
			int nInts = sq.getNInts();
			for(int i=0; i<nInts; i++) {
				Integral I(sq,i);
				string name = I.getName();
				string line = "Double " + name + " = 0.0E0;";
				printLine(nSpace,line,myfile);
			}
		}
	}
	myfile << endl;

}

void SQIntsPrint::printTwoBodyOverlapHead(ofstream& file) const 
{
	///////////////////////////////////////////////////////
	//     information for the given shell quartets      //
	///////////////////////////////////////////////////////
	bool hasRROnBRA1 = hasVRROnVar(BRA1);
	bool hasRROnBRA2 = hasVRROnVar(BRA2);
	bool hasRROnBRA  = hasRROnBRA1 || hasRROnBRA2;

	// for composite shell, we will handle the coefficients
	// in the final step
	bool comSQ = infor.isComSQ();
	int nBraCoeArray = infor.getCoeArrayLength(BRA);

	///////////////////////////////////////////////////////
	//                     bra side                      //
	///////////////////////////////////////////////////////
	string line = "for(UInt ip2=0; ip2<inp2; ip2++) {";
	printLine(2,line,file);

	// coefficients and exponents
	// exponents is only useful for RR work
	line = "Double ic2   = icoe[ip2];";
	printLine(4,line,file);
	for(int i=1; i<nBraCoeArray; i++) {
		string array= "icoe";
		string lhs  = "ic2_" + lexical_cast<string>(i);
		string pos  = "ip2+"  + lexical_cast<string>(i) + "*" + "inp2";
		string rhs  = array  + "[" + pos + "];"; 
		line = "Double " + lhs + " = " + rhs;
		printLine(4,line,file);
	}

	// prefactors for the SSSS integrals
	line = "Double fbra  = ifac[ip2];";
	printLine(4,line,file);

	// if we do RR on BRA side, we need 
	// them
	if (hasRROnBRA) {
		line = "Double onedz = iexp[ip2];";
		printLine(4,line,file);
		line = "Double oned2z= 0.5E0*onedz;";
		printLine(4,line,file);
		line = "UInt offsetP  = 3*ip2;";
		printLine(4,line,file);
		line = "Double PX    = P[offsetP  ];";
		printLine(4,line,file);
		line = "Double PY    = P[offsetP+1];";
		printLine(4,line,file);
		line = "Double PZ    = P[offsetP+2];";
		printLine(4,line,file);
	}

	// if we do RR on BRA1
	if (hasRROnBRA1) {
		line = "Double PAX   = PX - A[0];";
		printLine(4,line,file);
		line = "Double PAY   = PY - A[1];";
		printLine(4,line,file);
		line = "Double PAZ   = PZ - A[2];";
		printLine(4,line,file);
	}

	// if we do RR on BRA2
	if (hasRROnBRA2) {
		line = "Double PBX   = PX - B[0];";
		printLine(4,line,file);
		line = "Double PBY   = PY - B[1];";
		printLine(4,line,file);
		line = "Double PBZ   = PZ - B[2];";
		printLine(4,line,file);
	}

	// now let's go to generate the S integrals
	if (! comSQ) {
		line = "Double I_TWOBODYOVERLAP_S_S = ic2*fbra;";
		printLine(4,line,file);
	}else{
		line = "Double I_TWOBODYOVERLAP_S_S = fbra;";
		printLine(4,line,file);
	}

	// whether we pass it
	line = "if(fabs(I_TWOBODYOVERLAP_S_S)<THRESHOLD_MATH) continue;";
	printLine(4,line,file);
	file << endl;
}

void SQIntsPrint::printMOMHead(ofstream& file) const 
{
	///////////////////////////////////////////////////////
	//     information for the given shell quartets      //
	//     in default, since we need PCX etc.            //
	//     therefore we need the P point infor           //
	///////////////////////////////////////////////////////
	bool hasRROnBRA1 = hasVRROnVar(BRA1);
	bool hasRROnBRA2 = hasVRROnVar(BRA2);
	bool hasRROnBRA  = true;

	// for composite shell, we will handle the coefficients
	// in the final step
	bool comSQ = infor.isComSQ();
	int nBraCoeArray = infor.getCoeArrayLength(BRA);

	///////////////////////////////////////////////////////
	//                     bra side                      //
	///////////////////////////////////////////////////////
	string line = "for(UInt ip2=0; ip2<inp2; ip2++) {";
	printLine(2,line,file);

	// coefficients and exponents
	// exponents is only useful for RR work
	line = "Double ic2   = icoe[ip2];";
	printLine(4,line,file);
	for(int i=1; i<nBraCoeArray; i++) {
		string array= "icoe";
		string lhs  = "ic2_" + lexical_cast<string>(i);
		string pos  = "ip2+"  + lexical_cast<string>(i) + "*" + "inp2";
		string rhs  = array  + "[" + pos + "];"; 
		line = "Double " + lhs + " = " + rhs;
		printLine(4,line,file);
	}

	// prefactors for the SSSS integrals
	line = "Double fbra  = ifac[ip2];";
	printLine(4,line,file);

	// in default, we need this to print out the bottom integrals
	if (hasRROnBRA) {
		line = "Double onedz = iexp[ip2];";
		printLine(4,line,file);
		line = "Double oned2z= 0.5E0*onedz;";
		printLine(4,line,file);
		line = "UInt offsetP  = 3*ip2;";
		printLine(4,line,file);
		line = "Double PX    = P[offsetP  ];";
		printLine(4,line,file);
		line = "Double PY    = P[offsetP+1];";
		printLine(4,line,file);
		line = "Double PZ    = P[offsetP+2];";
		printLine(4,line,file);
		line = "Double PCX   = PX - C[0];";
		printLine(4,line,file);
		line = "Double PCY   = PY - C[1];";
		printLine(4,line,file);
		line = "Double PCZ   = PZ - C[2];";
		printLine(4,line,file);
	}

	// if we do RR on BRA1
	if (hasRROnBRA1) {
		line = "Double PAX   = PX - A[0];";
		printLine(4,line,file);
		line = "Double PAY   = PY - A[1];";
		printLine(4,line,file);
		line = "Double PAZ   = PZ - A[2];";
		printLine(4,line,file);
	}

	// if we do RR on BRA2
	if (hasRROnBRA2) {
		line = "Double PBX   = PX - B[0];";
		printLine(4,line,file);
		line = "Double PBY   = PY - B[1];";
		printLine(4,line,file);
		line = "Double PBZ   = PZ - B[2];";
		printLine(4,line,file);
	}

	// now let's go to generate the S integrals
	if (! comSQ) {
		line = "Double I_TWOBODYOVERLAP_S_S = ic2*fbra;";
		printLine(4,line,file);
	}else{
		line = "Double I_TWOBODYOVERLAP_S_S = fbra;";
		printLine(4,line,file);
	}

	// whether we pass it
	line = "if(fabs(I_TWOBODYOVERLAP_S_S)<THRESHOLD_MATH) continue;";
	printLine(4,line,file);
	file << endl;
}

void SQIntsPrint::printKineticHead(ofstream& file) const 
{
	///////////////////////////////////////////////////////
	//     information for the given shell quartets      //
	///////////////////////////////////////////////////////
	bool hasRROnBRA1 = hasVRROnVar(BRA1);
	bool hasRROnBRA2 = hasVRROnVar(BRA2);
	bool hasRROnBRA  = hasRROnBRA1 || hasRROnBRA2;

	// for composite shell, we will handle the coefficients
	// in the final step
	bool comSQ = infor.isComSQ();
	int nBraCoeArray = infor.getCoeArrayLength(BRA);

	///////////////////////////////////////////////////////
	//                     bra side                      //
	///////////////////////////////////////////////////////
	string line = "Double AB2 = (A[0]-B[0])*(A[0]-B[0])+(A[1]-B[1])*(A[1]-B[1])"
		"+(A[2]-B[2])*(A[2]-B[2]);";
	printLine(2,line,file);
	line = "for(UInt ip2=0; ip2<inp2; ip2++) {";
	printLine(2,line,file);

	// coefficients and exponents
	line = "Double ic2   = icoe[ip2];";
	printLine(4,line,file);
	for(int i=1; i<nBraCoeArray; i++) {
		string array= "icoe";
		string lhs  = "ic2_" + lexical_cast<string>(i);
		string pos  = "ip2+"  + lexical_cast<string>(i) + "*" + "inp2";
		string rhs  = array  + "[" + pos + "];"; 
		line = "Double " + lhs + " = " + rhs;
		printLine(4,line,file);
	}

	// exponents
	line = "Double onedz = iexp[ip2];";
	printLine(4,line,file);
	line = "Double zeta  = 1.0E0/onedz;";
	printLine(4,line,file);
	line = "Double zdiff = iexpdiff[ip2];";
	printLine(4,line,file);
	line = "Double alpha = 0.5E0*(zeta+zdiff);";
	printLine(4,line,file);
	line = "Double beta  = 0.5E0*(zeta-zdiff);";
	printLine(4,line,file);
	line = "Double xi    = alpha*beta*onedz;";
	printLine(4,line,file);
	line = "Double twoxi = 2.0E0*xi;";
	printLine(4,line,file);

	// prefactors for the SSSS integrals
	line = "Double fbra  = ifac[ip2];";
	printLine(4,line,file);

	// if we do RR on BRA side, we need 
	// them
	if (hasRROnBRA) {
		line = "Double oned2z= 0.5E0*onedz;";
		printLine(4,line,file);
		line = "Double adz   = alpha*onedz;";
		printLine(4,line,file);
		line = "Double bdz   = beta*onedz;";
		printLine(4,line,file);
		line = "UInt offsetP  = 3*ip2;";
		printLine(4,line,file);
		line = "Double PX    = P[offsetP  ];";
		printLine(4,line,file);
		line = "Double PY    = P[offsetP+1];";
		printLine(4,line,file);
		line = "Double PZ    = P[offsetP+2];";
		printLine(4,line,file);
	}

	// if we do RR on BRA1
	if (hasRROnBRA1) {
		line = "Double PAX   = PX - A[0];";
		printLine(4,line,file);
		line = "Double PAY   = PY - A[1];";
		printLine(4,line,file);
		line = "Double PAZ   = PZ - A[2];";
		printLine(4,line,file);
	}

	// if we do RR on BRA2
	if (hasRROnBRA2) {
		line = "Double PBX   = PX - B[0];";
		printLine(4,line,file);
		line = "Double PBY   = PY - B[1];";
		printLine(4,line,file);
		line = "Double PBZ   = PZ - B[2];";
		printLine(4,line,file);
	}

	// now let's go to generate the S integrals
	if (! comSQ) {
		line = "Double I_KINETIC_S_S = ic2*fbra*xi*(3.0E0-twoxi*AB2);";
		printLine(4,line,file);
		line = "Double I_TWOBODYOVERLAP_S_S = ic2*fbra;";
		printLine(4,line,file);
	}else{
		line = "Double I_KINETIC_S_S = fbra*xi*(3.0E0-twoxi*AB2);";
		printLine(4,line,file);
		line = "Double I_TWOBODYOVERLAP_S_S = fbra;";
		printLine(4,line,file);
	}

	// whether we pass it
	line = "if(fabs(I_KINETIC_S_S)<THRESHOLD_MATH) continue;";
	printLine(4,line,file);
	file << endl;
}

void SQIntsPrint::printNAIHead(ofstream& file) const 
{
	///////////////////////////////////////////////////////
	//     information for the given shell quartets      //
	///////////////////////////////////////////////////////
	bool hasRROnBRA1 = hasVRROnVar(BRA1);
	bool hasRROnBRA2 = hasVRROnVar(BRA2);
	bool hasRROnBRA  = hasRROnBRA1 || hasRROnBRA2;
	int  maxLSum     = getMaxLSum();

	// for composite shell, we will handle the coefficients
	// in the final step
	bool comSQ = infor.isComSQ();
	int nBraCoeArray = infor.getCoeArrayLength(BRA);

	///////////////////////////////////////////////////////
	//                     bra side                      //
	///////////////////////////////////////////////////////
	string line = "for(UInt ip2=0; ip2<inp2; ip2++) {";
	printLine(2,line,file);

	// coefficients 
	line = "Double ic2   = icoe[ip2];";
	printLine(4,line,file);
	for(int i=1; i<nBraCoeArray; i++) {
		string array= "icoe";
		string lhs  = "ic2_" + lexical_cast<string>(i);
		string pos  = "ip2+"  + lexical_cast<string>(i) + "*" + "inp2";
		string rhs  = array  + "[" + pos + "];"; 
		line = "Double " + lhs + " = " + rhs;
		printLine(4,line,file);
	}

	// exponents
	line = "Double onedz = iexp[ip2];";
	printLine(4,line,file);
	line = "Double rho   = 1.0E0/onedz;";
	printLine(4,line,file);
	line = "Double sqrho = sqrt(rho);";
	printLine(4,line,file);

	// prefactors for the SSSS integrals
	line = "Double fbra  = ifac[ip2];";
	printLine(4,line,file);

	// if we do RR on BRA side, we need 
	// them
	if (hasRROnBRA) {
		line = "Double oned2z= 0.5E0*onedz;";
		printLine(4,line,file);
	}

	// we need the new center
	line = "UInt offsetP  = 3*ip2;";
	printLine(4,line,file);
	line = "Double PX    = P[offsetP  ];";
	printLine(4,line,file);
	line = "Double PY    = P[offsetP+1];";
	printLine(4,line,file);
	line = "Double PZ    = P[offsetP+2];";
	printLine(4,line,file);

	// if we do RR on BRA1
	if (hasRROnBRA1) {
		line = "Double PAX   = PX - A[0];";
		printLine(4,line,file);
		line = "Double PAY   = PY - A[1];";
		printLine(4,line,file);
		line = "Double PAZ   = PZ - A[2];";
		printLine(4,line,file);
	}

	// if we do RR on BRA2
	if (hasRROnBRA2) {
		line = "Double PBX   = PX - B[0];";
		printLine(4,line,file);
		line = "Double PBY   = PY - B[1];";
		printLine(4,line,file);
		line = "Double PBZ   = PZ - B[2];";
		printLine(4,line,file);
	}

	// now loop over the nuclear center
   line = "for(UInt iAtom=0; iAtom<nAtoms; iAtom++) {";
	printLine(4,line,file);
	line = "Double PNX   = PX - N[iAtom*3  ];";
	printLine(6,line,file);
	line = "Double PNY   = PY - N[iAtom*3+1];";
	printLine(6,line,file);
	line = "Double PNZ   = PZ - N[iAtom*3+2];";
	printLine(6,line,file);
	line = "Double PN2   = PNX*PNX+PNY*PNY+PNZ*PNZ;";
	printLine(6,line,file);
	line = "Double charge= Z[iAtom];";
	printLine(6,line,file);
	line = "Double u     = rho*PN2;";
	printLine(6,line,file);
	line = "Double squ   = sqrt(u);";
	printLine(6,line,file);
	line = "#ifndef WITH_SINGLE_PRECISION";
	printLine(0,line,file);
	line = "Double erfVal= erf(squ);";
	printLine(6,line,file);
	line = "#endif";
	printLine(0,line,file);
	if (! comSQ) {
		line = "Double prefactor = -ic2*charge*fbra;";
		printLine(6,line,file);
	}else{
		line = "Double prefactor = -charge*fbra;";
		printLine(6,line,file);
	}

	// now calculate the bottom integrals
	fmtIntegralsGeneration(maxLSum,NAI,6,file);
}

void SQIntsPrint::printESPHead(ofstream& file) const 
{
	///////////////////////////////////////////////////////
	//     information for the given shell quartets      //
	///////////////////////////////////////////////////////
	bool hasRROnBRA1 = hasVRROnVar(BRA1);
	bool hasRROnBRA2 = hasVRROnVar(BRA2);
	bool hasRROnBRA  = hasRROnBRA1 || hasRROnBRA2;
	int  maxLSum     = getMaxLSum();

	// for composite shell, we will handle the coefficients
	// in the final step
	bool comSQ = infor.isComSQ();
	int nBraCoeArray = infor.getCoeArrayLength(BRA);

	///////////////////////////////////////////////////////
	//                     bra side                      //
	///////////////////////////////////////////////////////
	string line = "for(UInt ip2=0; ip2<inp2; ip2++) {";
	printLine(4,line,file);

	// coefficients 
	line = "Double ic2   = icoe[ip2];";
	printLine(6,line,file);
	for(int i=1; i<nBraCoeArray; i++) {
		string array= "icoe";
		string lhs  = "ic2_" + lexical_cast<string>(i);
		string pos  = "ip2+"  + lexical_cast<string>(i) + "*" + "inp2";
		string rhs  = array  + "[" + pos + "];"; 
		line = "Double " + lhs + " = " + rhs;
		printLine(6,line,file);
	}

	// exponents
	line = "Double onedz = iexp[ip2];";
	printLine(6,line,file);
	line = "Double rho   = 1.0E0/onedz;";
	printLine(6,line,file);
	line = "Double sqrho = sqrt(rho);";
	printLine(6,line,file);

	// prefactors for the SSSS integrals
	line = "Double fbra  = ifac[ip2];";
	printLine(6,line,file);

	// if we do RR on BRA side, we need 
	// them
	if (hasRROnBRA) {
		line = "Double oned2z= 0.5E0*onedz;";
		printLine(6,line,file);
	}

	// we need the new center
	line = "UInt offsetP  = 3*ip2;";
	printLine(6,line,file);
	line = "Double PX    = P[offsetP  ];";
	printLine(6,line,file);
	line = "Double PY    = P[offsetP+1];";
	printLine(6,line,file);
	line = "Double PZ    = P[offsetP+2];";
	printLine(6,line,file);

	// if we do RR on BRA1
	if (hasRROnBRA1) {
		line = "Double PAX   = PX - A[0];";
		printLine(6,line,file);
		line = "Double PAY   = PY - A[1];";
		printLine(6,line,file);
		line = "Double PAZ   = PZ - A[2];";
		printLine(6,line,file);
	}

	// if we do RR on BRA2
	if (hasRROnBRA2) {
		line = "Double PBX   = PX - B[0];";
		printLine(6,line,file);
		line = "Double PBY   = PY - B[1];";
		printLine(6,line,file);
		line = "Double PBZ   = PZ - B[2];";
		printLine(6,line,file);
	}

	// now loop over the nuclear center
	line = "Double PRX   = PX - R[iGrid*3  ];";
	printLine(6,line,file);
	line = "Double PRY   = PY - R[iGrid*3+1];";
	printLine(6,line,file);
	line = "Double PRZ   = PZ - R[iGrid*3+2];";
	printLine(6,line,file);
	line = "Double PR2   = PRX*PRX+PRY*PRY+PRZ*PRZ;";
	printLine(6,line,file);
	line = "Double u     = rho*PR2;";
	printLine(6,line,file);
	line = "Double squ   = sqrt(u);";
	printLine(6,line,file);
	line = "#ifndef WITH_SINGLE_PRECISION";
	printLine(0,line,file);
	line = "Double erfVal= erf(squ);";
	printLine(6,line,file);
	line = "#endif";
	printLine(0,line,file);
	if (! comSQ) {
		line = "Double prefactor = ic2*fbra;";
		printLine(6,line,file);
	}else{
		line = "Double prefactor = fbra;";
		printLine(6,line,file);
	}

	// now calculate the bottom integrals
	fmtIntegralsGeneration(maxLSum,ESP,6,file);
}

void SQIntsPrint::printThreeBodyOverlapHead(ofstream& file) const 
{
	///////////////////////////////////////////////////////
	//     information for the given shell quartets      //
	///////////////////////////////////////////////////////
	bool hasRROnBRA1 = hasVRROnVar(BRA1);
	bool hasRROnBRA2 = hasVRROnVar(BRA2);
	bool hasRROnKET1 = hasVRROnVar(KET1);
	bool hasRR       = hasRROnBRA1  || hasRROnBRA2 || hasRROnKET1;

	// for composite shell, we will handle the coefficients
	// in the final step
	bool comSQ = infor.isComSQ();
	int nBraCoeArray = infor.getCoeArrayLength(BRA);
	int nKetCoeArray = infor.getCoeArrayLength(KET); 

	///////////////////////////////////////////////////////
	//                     bra side                      //
	///////////////////////////////////////////////////////
	string line = "for(UInt ip2=0; ip2<inp2; ip2++) {";
	printLine(2,line,file);

	// coefficients 
	line = "Double ic2   = icoe[ip2];";
	printLine(4,line,file);
	for(int i=1; i<nBraCoeArray; i++) {
		string array= "icoe";
		string lhs  = "ic2_" + lexical_cast<string>(i);
		string pos  = "ip2+"  + lexical_cast<string>(i) + "*" + "inp2";
		string rhs  = array  + "[" + pos + "];"; 
		line = "Double " + lhs + " = " + rhs;
		printLine(4,line,file);
	}

	// exponents
	line = "Double onedz = iexp[ip2];";
	printLine(4,line,file);

	// prefactors for the SSSS integrals
	line = "Double fbra  = ifac[ip2];";
	printLine(4,line,file);

	// we need the new center even for the SSS integral
	line = "UInt offsetP  = 3*ip2;";
	printLine(4,line,file);
	line = "Double PX    = P[offsetP  ];";
	printLine(4,line,file);
	line = "Double PY    = P[offsetP+1];";
	printLine(4,line,file);
	line = "Double PZ    = P[offsetP+2];";
	printLine(4,line,file);

	// now loop over the primitives on ket1
	line = "for(UInt jp2=0; jp2<jnp2; jp2++) {";
	printLine(4,line,file);
	line = "Double onede = jexp[jp2];";
	printLine(6,line,file);
	line = "Double jc2   = jcoe[jp2];";
	printLine(6,line,file);
	for(int i=1; i<nKetCoeArray; i++) {
		string array= "jcoe";
		string lhs  = "jc2_" + lexical_cast<string>(i);
		string pos  = "jp2+"  + lexical_cast<string>(i) + "*" + "jnp2";
		string rhs  = array  + "[" + pos + "];"; 
		line = "Double " + lhs + " = " + rhs;
		printLine(6,line,file);
	}

	// now combine the bra and ket part together 
	// to generate rho etc. - this is needed for S integral
	line = "Double rho   = 1.0E0/(onedz+onede);";
	printLine(6,line,file);
	line = "Double PC2   = (PX-C[0])*(PX-C[0])+(PY-C[1])*(PY-C[1])+(PZ-C[2])*(PZ-C[2]);";
	printLine(6,line,file);

	// now combine the bra and ket part into the new center 
	// G. This is needed when bra/ket needs RR work
	if (hasRR) {
		line = "Double GX    = rho*(PX*onede + C[0]*onedz);";
		printLine(6,line,file);
		line = "Double GY    = rho*(PY*onede + C[1]*onedz);";
		printLine(6,line,file);
		line = "Double GZ    = rho*(PZ*onede + C[2]*onedz);";
		printLine(6,line,file);
	}

	// now create the RR coefficients
	if (hasRR) {
		line = "Double oned2zeta = 0.5E0*rho*onede*onedz;";
		printLine(6,line,file);
	}

	// if we do RR on BRA1
	if (hasRROnBRA1) {
		line = "Double GAX   = GX - A[0];";
		printLine(6,line,file);
		line = "Double GAY   = GY - A[1];";
		printLine(6,line,file);
		line = "Double GAZ   = GZ - A[2];";
		printLine(6,line,file);
	}

	// if we do RR on BRA2
	if (hasRROnBRA2) {
		line = "Double GBX   = GX - B[0];";
		printLine(6,line,file);
		line = "Double GBY   = GY - B[1];";
		printLine(6,line,file);
		line = "Double GBZ   = GZ - B[2];";
		printLine(6,line,file);
	}

	// if we do RR on KET1
	if (hasRROnKET1) {
		line = "Double GCX   = GX - C[0];";
		printLine(6,line,file);
		line = "Double GCY   = GY - C[1];";
		printLine(6,line,file);
		line = "Double GCZ   = GZ - C[2];";
		printLine(6,line,file);
	}

	// now assemble overlap for bra part
	line = "Double p     = pow(onede*rho,1.5E0);";
	printLine(6,line,file);
	if (comSQ) {
		line = "Double I_THREEBODYOVERLAP_S_S_S = p*fbra*exp(-rho*PC2);";
		printLine(6,line,file);
	}else{
		line = "Double I_THREEBODYOVERLAP_S_S_S = ic2*jc2*p*fbra*exp(-rho*PC2);";
		printLine(6,line,file);
	}
	line = "if(fabs(I_THREEBODYOVERLAP_S_S_S)<THRESHOLD_MATH) continue;";
	printLine(6,line,file);
	file << endl;
}

void SQIntsPrint::printThreeBodyKIHead(ofstream& file) const 
{
	///////////////////////////////////////////////////////
	//     information for the given shell quartets      //
	//     for three body kinetic energy, RR will be     //
	//     done for all of its three positions           //
	///////////////////////////////////////////////////////
	bool hasRROnBRA1 = true;
	bool hasRROnBRA2 = true;
	bool hasRROnKET1 = true;
	bool hasRR       = true;

	// for composite shell, we will handle the coefficients
	// in the final step
	bool comSQ = infor.isComSQ();
	int nBraCoeArray = infor.getCoeArrayLength(BRA);
	int nKetCoeArray = infor.getCoeArrayLength(KET); 

	///////////////////////////////////////////////////////
	//                     bra side                      //
	///////////////////////////////////////////////////////
	string line = "for(UInt ip2=0; ip2<inp2; ip2++) {";
	printLine(2,line,file);

	// coefficients 
	line = "Double ic2   = icoe[ip2];";
	printLine(4,line,file);
	for(int i=1; i<nBraCoeArray; i++) {
		string array= "icoe";
		string lhs  = "ic2_" + lexical_cast<string>(i);
		string pos  = "ip2+"  + lexical_cast<string>(i) + "*" + "inp2";
		string rhs  = array  + "[" + pos + "];"; 
		line = "Double " + lhs + " = " + rhs;
		printLine(4,line,file);
	}

	// exponents
	// we need both alpha and beta exponents
	line = "Double onedz = iexp[ip2];";
	printLine(4,line,file);
	line = "Double zeta  = 1.0E0/onedz;";
	printLine(4,line,file);
	line = "Double zdiff = iexpdiff[ip2];";
	printLine(4,line,file);
	line = "Double alpha = 0.5E0*(zeta+zdiff);";
	printLine(4,line,file);
	line = "Double beta  = 0.5E0*(zeta-zdiff);";
	printLine(4,line,file);

	// prefactors for the SSSS integrals
	line = "Double fbra  = ifac[ip2];";
	printLine(4,line,file);

	// we need the new center even for the SSS integral
	line = "UInt offsetP  = 3*ip2;";
	printLine(4,line,file);
	line = "Double PX    = P[offsetP  ];";
	printLine(4,line,file);
	line = "Double PY    = P[offsetP+1];";
	printLine(4,line,file);
	line = "Double PZ    = P[offsetP+2];";
	printLine(4,line,file);

	// now loop over the primitives on ket1
	line = "for(UInt jp2=0; jp2<jnp2; jp2++) {";
	printLine(4,line,file);
	line = "Double onede = jexp[jp2];";
	printLine(6,line,file);
	line = "Double gamma = 1.0E0/onede;";
	printLine(6,line,file);
	line = "Double jc2   = jcoe[jp2];";
	printLine(6,line,file);
	for(int i=1; i<nKetCoeArray; i++) {
		string array= "jcoe";
		string lhs  = "jc2_" + lexical_cast<string>(i);
		string pos  = "jp2+"  + lexical_cast<string>(i) + "*" + "jnp2";
		string rhs  = array  + "[" + pos + "];"; 
		line = "Double " + lhs + " = " + rhs;
		printLine(6,line,file);
	}

	// now combine the bra and ket part together 
	// to generate rho etc. - this is needed for S integral
	line = "Double rho   = 1.0E0/(onedz+onede);";
	printLine(6,line,file);
	line = "Double PC2   = (PX-C[0])*(PX-C[0])+(PY-C[1])*(PY-C[1])+(PZ-C[2])*(PZ-C[2]);";
	printLine(6,line,file);

	// now combine the bra and ket part into the new center 
	// G. This is needed when bra/ket needs RR work
	if (hasRR) {
		line = "Double GX    = rho*(PX*onede + C[0]*onedz);";
		printLine(6,line,file);
		line = "Double GY    = rho*(PY*onede + C[1]*onedz);";
		printLine(6,line,file);
		line = "Double GZ    = rho*(PZ*onede + C[2]*onedz);";
		printLine(6,line,file);
	}

	// now create the RR coefficients
	if (hasRR) {
		line = "Double oned2zeta = 0.5E0*rho*onede*onedz;";
		printLine(6,line,file);
	}

	// if we do RR on BRA1
	if (hasRROnBRA1) {
		line = "Double GAX   = GX - A[0];";
		printLine(6,line,file);
		line = "Double GAY   = GY - A[1];";
		printLine(6,line,file);
		line = "Double GAZ   = GZ - A[2];";
		printLine(6,line,file);
	}

	// if we do RR on BRA2
	if (hasRROnBRA2) {
		line = "Double GBX   = GX - B[0];";
		printLine(6,line,file);
		line = "Double GBY   = GY - B[1];";
		printLine(6,line,file);
		line = "Double GBZ   = GZ - B[2];";
		printLine(6,line,file);
	}

	// if we do RR on KET1
	if (hasRROnKET1) {
		line = "Double GCX   = GX - C[0];";
		printLine(6,line,file);
		line = "Double GCY   = GY - C[1];";
		printLine(6,line,file);
		line = "Double GCZ   = GZ - C[2];";
		printLine(6,line,file);
	}

	// now assemble overlap for bra part
	line = "Double p     = pow(onede*rho,1.5E0);";
	printLine(6,line,file);
	if (comSQ) {
		line = "Double I_THREEBODYOVERLAP_S_S_S = p*fbra*exp(-rho*PC2);";
		printLine(6,line,file);
	}else{
		line = "Double I_THREEBODYOVERLAP_S_S_S = ic2*jc2*p*fbra*exp(-rho*PC2);";
		printLine(6,line,file);
	}
	line = "if(fabs(I_THREEBODYOVERLAP_S_S_S)<THRESHOLD_MATH) continue;";
	printLine(6,line,file);
	file << endl;

	// now compute the SSS integrals
	if (infor.areAllBottomSQ()) {
		line = "Double f   = (alpha+beta)*gamma/(alpha+beta+gamma);";
		printLine(6,line,file);
		line = "Double f1x = 2.0E0*(alpha*GAX+beta*GBX)*gamma*GCX+f;";
		printLine(6,line,file);
		line = "Double f1y = 2.0E0*(alpha*GAY+beta*GBY)*gamma*GCY+f;";
		printLine(6,line,file);
		line = "Double f1z = 2.0E0*(alpha*GAZ+beta*GBZ)*gamma*GCZ+f;";
		printLine(6,line,file);
		line = "Double I_THREEBODYKI_S_S_S = (f1x+f1y+f1z)*I_THREEBODYOVERLAP_S_S_S;";
		printLine(6,line,file);
		file << endl;
	}
}

void SQIntsPrint::printERIHead(ofstream& file) const 
{
	///////////////////////////////////////////////////////
	//     information for the given shell quartets      //
	///////////////////////////////////////////////////////
	bool hasRROnBRA1 = hasVRROnVar(BRA1);
	bool hasRROnBRA2 = hasVRROnVar(BRA2);
	bool hasRROnKET1 = hasVRROnVar(KET1);
	bool hasRROnKET2 = hasVRROnVar(KET2);
	bool hasRROnBRA  = hasRROnBRA1 || hasRROnBRA2;
	bool hasRROnKET  = hasRROnKET1 || hasRROnKET2;
	bool hasRR       = hasRROnBRA  || hasRROnKET;
	int  maxLSum     = getMaxLSum();

	// for composite shell, we will handle the coefficients
	// in the final step
	bool comSQ = infor.isComSQ();
	int nBraCoeArray = infor.getCoeArrayLength(BRA);
	int nKetCoeArray = infor.getCoeArrayLength(KET); 

	///////////////////////////////////////////////////////
	//                     bra side                      //
	///////////////////////////////////////////////////////
	string line = "for(UInt ip2=0; ip2<inp2; ip2++) {";
	printLine(2,line,file);

	// coefficients and exponents
	line = "Double onedz = iexp[ip2];";
	printLine(4,line,file);
	line = "Double ic2   = icoe[ip2];";
	printLine(4,line,file);
	for(int i=1; i<nBraCoeArray; i++) {
		string array= "icoe";
		string lhs  = "ic2_" + lexical_cast<string>(i);
		string pos  = "ip2+"  + lexical_cast<string>(i) + "*" + "inp2";
		string rhs  = array  + "[" + pos + "];"; 
		line = "Double " + lhs + " = " + rhs;
		printLine(4,line,file);
	}

	// prefactors for the SSSS integrals
	line = "Double fbra  = ifac[ip2];";
	printLine(4,line,file);

	// if we do RR on BRA side, we need 
	// them
	if (hasRROnBRA) {
		line = "Double oned2z= 0.5E0*onedz;";
		printLine(4,line,file);
	}

	// however, even for S integral calculation, we need
	// P for calculating |PQ|
	line = "UInt offsetP  = 3*ip2;";
	printLine(4,line,file);
	line = "Double PX    = P[offsetP  ];";
	printLine(4,line,file);
	line = "Double PY    = P[offsetP+1];";
	printLine(4,line,file);
	line = "Double PZ    = P[offsetP+2];";
	printLine(4,line,file);

	// if we do RR on BRA1
	if (hasRROnBRA1) {
		line = "Double PAX   = PX - A[0];";
		printLine(4,line,file);
		line = "Double PAY   = PY - A[1];";
		printLine(4,line,file);
		line = "Double PAZ   = PZ - A[2];";
		printLine(4,line,file);
	}

	// if we do RR on BRA2
	if (hasRROnBRA2) {
		line = "Double PBX   = PX - B[0];";
		printLine(4,line,file);
		line = "Double PBY   = PY - B[1];";
		printLine(4,line,file);
		line = "Double PBZ   = PZ - B[2];";
		printLine(4,line,file);
	}

	///////////////////////////////////////////////////////
	//                     ket side                      //
	///////////////////////////////////////////////////////
	line = "for(UInt jp2=0; jp2<jnp2; jp2++) {";
	printLine(4,line,file);

	// coefficients and exponents
	line = "Double onede = jexp[jp2];";
	printLine(6,line,file);
	line = "Double jc2   = jcoe[jp2];";
	printLine(6,line,file);
	for(int i=1; i<nKetCoeArray; i++) {
		string array= "jcoe";
		string lhs  = "jc2_" + lexical_cast<string>(i);
		string pos  = "jp2+"  + lexical_cast<string>(i) + "*" + "jnp2";
		string rhs  = array  + "[" + pos + "];"; 
		line = "Double " + lhs + " = " + rhs;
		printLine(6,line,file);
	}

	// prefactors for the SSSS integrals
	line = "Double fket  = jfac[jp2];";
	printLine(6,line,file);

	// we need Q coordinates
	line = "UInt offsetQ  = 3*jp2;";
	printLine(6,line,file);
	line = "Double QX    = Q[offsetQ  ];";
	printLine(6,line,file);
	line = "Double QY    = Q[offsetQ+1];";
	printLine(6,line,file);
	line = "Double QZ    = Q[offsetQ+2];";
	printLine(6,line,file);

	// based on bra and ket part, generate the prefactor
	// as well as other things to form (SS|SS)^{m} 
	// integrals
	if (comSQ) {
		line = "Double pref      = fbra*fket;";
		printLine(6,line,file);
		line = "Double prefactor = pref;";
		printLine(6,line,file);
	}else{
		line = "Double pref      = fbra*fket;";
		printLine(6,line,file);
		line = "Double prefactor = ic2*jc2*pref;";
		printLine(6,line,file);
	}
	line = "Double rho   = 1.0E0/(onedz+onede);";
	printLine(6,line,file);
	line = "Double sqrho = sqrt(rho);";
	printLine(6,line,file);
	line = "Double PQ2   = (PX-QX)*(PX-QX)+(PY-QY)*(PY-QY)+(PZ-QZ)*(PZ-QZ);";
	printLine(6,line,file);
	line = "Double u     = rho*PQ2;";
	printLine(6,line,file);
	line = "Double squ   = sqrt(u);";
	printLine(6,line,file);
	line = "#ifndef WITH_SINGLE_PRECISION";
	printLine(0,line,file);
	line = "Double erfVal= erf(squ);";
	printLine(6,line,file);
	line = "#endif";
	printLine(0,line,file);

	// now test the significance of integrals
	fmtIntegralsTest(maxLSum,ERI,6,file);
	file << endl;

	// if we do RR on KET1
	if (hasRROnKET1) {
		line = "Double QCX   = QX - C[0];";
		printLine(6,line,file);
		line = "Double QCY   = QY - C[1];";
		printLine(6,line,file);
		line = "Double QCZ   = QZ - C[2];";
		printLine(6,line,file);
	}

	// if we do RR on KET2
	if (hasRROnKET2) {
		line = "Double QDX   = QX - D[0];";
		printLine(6,line,file);
		line = "Double QDY   = QY - D[1];";
		printLine(6,line,file);
		line = "Double QDZ   = QZ - D[2];";
		printLine(6,line,file);
	}


	// now combine the bra and ket part into the new center 
	// W. This is needed when bra/ket needs RR work
	if (hasRR) {
		line = "Double WX    = rho*(PX*onede + QX*onedz);";
		printLine(6,line,file);
		line = "Double WY    = rho*(PY*onede + QY*onedz);";
		printLine(6,line,file);
		line = "Double WZ    = rho*(PZ*onede + QZ*onedz);";
		printLine(6,line,file);
		line = "Double oned2k= 0.5E0*rho*onede*onedz;";
		printLine(6,line,file);
	}

	// if BRA part is both S integral, we do not need it
	// used only in RR
	if (hasRROnBRA) {
		line = "Double WPX   = WX - PX;";
		printLine(6,line,file);
		line = "Double WPY   = WY - PY;";
		printLine(6,line,file);
		line = "Double WPZ   = WZ - PZ;";
		printLine(6,line,file);
		line = "Double rhod2zsq = rho*oned2z*onedz;";
		printLine(6,line,file);
	}

	// if KET part is both S integral, we do not need it
	// used only in RR on ket side
	if (hasRROnKET) {
		line = "Double WQX   = WX - QX;";
		printLine(6,line,file);
		line = "Double WQY   = WY - QY;";
		printLine(6,line,file);
		line = "Double WQZ   = WZ - QZ;";
		printLine(6,line,file);
		line = "Double oned2e= 0.5E0*onede;";
		printLine(6,line,file);
		line = "Double rhod2esq= rho*oned2e*onede;";
		printLine(6,line,file);
	}
	file << endl;

	// now let's go to generate the S integrals
	fmtIntegralsGeneration(maxLSum,ERI,6,file);
}

//
// !!!! fifth print section: VRR contraction
// 1  we need to contract the VRR results into final results;
// 2  end the VRR process
//

void SQIntsPrint::vrrContraction(const string& fname) const
{
	// let's go to see whether this is for contraction on composite shells
	bool comSQ = infor.isComSQ();

	// determine that how many space should be given for each line printing
	int oper = infor.getOper();
	int nSpace = getNSpaceByOper(oper);
	if (infor.splitCPPFile()) nSpace = 2;

	// create the RR file
	// this is already in a mod of a+
	ofstream myfile;
	myfile.open(fname.c_str(),std::ofstream::app);

	// firstly, determine the number of sections
	// for contraction
	// currently many section situation just corresponding
	// to the case involving derivatives
	int nSect = 1;
	if (needIntDeriv(oper)) {
		nSect = getDerivNumberFromOper(oper,infor.getJobOrder());
	}
	
	// loop over sections first
	for(int isec=0; isec<nSect; isec++) {

		// check whether this is integral derivatives?
		int var = NULL_POS;
		if (needIntDeriv(oper)) {
			int derOrder = getDerivOrderFromOper(oper,infor.getJobOrder());
			var = getDerivVar(isec,derOrder);
		}

		// detect that wether we have compilcated offset for result?
		string additionalOffset;
		bool hasAdditionalOffset = resultIntegralHasAdditionalOffset(oper);
		if (hasAdditionalOffset) {
			int nInts = nTotalInts();
			additionalOffset = determineAdditionalOffset(oper,nInts);
		}

		//
		// now loop over the RR sq list
		// which is also the VRR's output results
		// there's one thing to note, that since the VRR
		// result must be in full integral list (that is 
		// to say, integral index is same with array index).
		// therefore, it's safe to convert the index into 
		// integral here. That is also why we do it freely
		// here
		// for more information, please refer to rrints.h
		//
		for(int iSQ=0; iSQ<(int)rrSQList.size(); iSQ++) {
			const ShellQuartet& sq = rrSQList[iSQ];

			//
			// now print out the comment for this section's 
			// contraction
			//
			myfile << endl;
			string line = "/************************************************************";
			printLine(nSpace,line,myfile);
			line = " * shell quartet name: " + sq.getName();
			printLine(nSpace,line,myfile);
			line = " * doing contraction work for VRR part ";
			printLine(nSpace,line,myfile);
			line = " ************************************************************/";
			printLine(nSpace,line,myfile);

			// according to the shell quartet, get the coe 
			// offset
			int ic2Offset = -1; 
			int jc2Offset = -1;
			if (comSQ) {
				infor.getCoeOffset(sq,ic2Offset,jc2Offset);
#ifdef SQINTS_DEBUG
				if (ic2Offset == -1 && jc2Offset == -1) {
					crash(true, "incorrect getCoeOffset result");
				}
#endif
			}

			// now let's form the coefficient part
			// we note that as long as this is composite 
			// shell quartet, that all of c2 offsets are >= 0
			// however, it's possible that ket part does not 
			// exist. Then the jc2Offset is set to -1
			string ic = "ic2";
			if (ic2Offset>0) ic += "_" + lexical_cast<string>(ic2Offset);
			string jc = "jc2";
			if (jc2Offset>0) jc += "_" + lexical_cast<string>(jc2Offset);
			string coe = ic;
			if (jc2Offset>=0) {
				coe += "*" + jc;
			}

			// is it the result shell quartet?
			// in that case, we will use name of "abcd" rather than
			// the shell quartet name
			bool isResult = infor.isResult(sq);
			int  sqOffset = 0;
			if (comSQ && isResult) sqOffset = infor.getOffset(sq);

			// now step into contraction
			string name = sq.getName();
			int nInts   = sq.getNInts();
			for(int i=0; i<nInts; i++) {

				//
				// LHS name
				// lhs will be the VRR result possibily with modifiers
				// such as division and exponent infor
				//
				Integral I(sq,i);
				string lhsName = sq.getName();
				if (! infor.withArrayIndex(HRR)) {
					lhsName = I.getName();
				}
				if (isResult) {
					lhsName = "abcd";
				}

				// index for lhs
				string lhsIndex;
				if (isResult || infor.withArrayIndex(HRR)) {

					// determine the offset
					int offset = i;
					if (isResult) offset = sqOffset + i;

					// counting the additional offset
					if (hasAdditionalOffset) {
						lhsIndex  = "[" + additionalOffset + "+" + lexical_cast<string>(offset) + "]";
					}else{
						lhsIndex  = "[" + lexical_cast<string>(offset) + "]";
					}
				}

				// finish lhs
				string lhs = lhsName + lhsIndex;

				//
				// RHS name
				// for rhs, we need to destroy the modifier information
				// also the name is the one used in RR process
				// we need to determine the name in RR process
				//
				// On the other hand, we note, that sometimes this sq 
				// may be the bottom integral (for example, high L in MOM case)
				// it uses the integral form rather than array form
				//
				I.destroyDivision();
				string rhsName = I.formVarName(rrType,MODULE_RESULT,var,comSQ);
				int vrr_method = infor.getVRRMethod();
				bool vrrWithArrayIndex = infor.withArrayIndex(vrr_method);
				if (vrrWithArrayIndex && ! sq.isSTypeSQ()) {
					ShellQuartet newSQ(sq);
					newSQ.destroyDivision();
					rhsName = newSQ.formArrayName(rrType,MODULE_RESULT,var,comSQ);
				}

				// add in coefficients to the RHS
				if (comSQ) rhsName = coe + "*" + rhsName; 

				// add in array index
				string rhs = rhsName;
				if (vrrWithArrayIndex && ! sq.isSTypeSQ()) {
					string rhs = rhsName + "[" + lexical_cast<string>(i) + "]";
				}

				// form the code
				string line = lhs + " += " + rhs + ";";
				printLine(nSpace,line,myfile);
			}
		}
	}

	// form the end of the vrr file
	// for split cpp file case, we will
	// do it when forming the function call code
	if (! infor.splitCPPFile()) {
		printVRREnd(myfile);
	}

	// now close the whole file
	myfile.close();

}

void SQIntsPrint::printVRREnd(ofstream& myfile) const
{
	// determine that how many space should be given for each line printing
	int oper = infor.getOper();
	int nSpace = getNSpaceByOper(oper);

	// for ESP etc. that nSpace gives the indentation for the whole cpp function
	// here we only close the VRR section
	// so consider to -2
	int nSpaceStop = 2;
	if (resultIntegralHasAdditionalOffset(oper)) nSpaceStop += 2;

	// finally, we need to add braket closure to the vrr body
	string line = "}";
	for(int iSpace= nSpace-2; iSpace>=nSpaceStop; iSpace = iSpace - 2) {
		printLine(iSpace,line,myfile);
	}

	// if it's using significance check and all bottom integral
	// we do nothing and just return 
	if(sigCheck(infor.getOper()) && infor.areAllBottomSQ()) {
		return;
	}

	// if we have significance test, then we may
	// test the significance first
	// if VRR step all of integrals are omitted,
	// then we do not need to do HRR accordingly
	if(sigCheck(infor.getOper())) {
		myfile << endl;
		line = "/************************************************************";
		printLine(2,line,myfile);
		line = " * let's see the significance test result";
		printLine(2,line,myfile);
		line = " * if the VRR step is insignificant, we do not need to do HRR";
		printLine(2,line,myfile);
		line = " ************************************************************/";
		printLine(2,line,myfile);

		// the handling of significance at the bottom of VRR 
		// this is to see whether we do HRR part
		// for the case that VRR/HRR are in a loop 
		// we can not return, just use continue
		line = "if (! isSignificant) return false;";
		if (resultIntegralHasAdditionalOffset(infor.getOper())) {
			line = "if (! isSignificant) continue;";
		}
		printLine(2,line,myfile);
		myfile << endl;
	}	
}

//
// !!!! sixth print section: HRR variable printing
//

void SQIntsPrint::printHRRSideVar(const int& side, const string& fileName) const
{
	// open file
	ofstream myfile;
	myfile.open(fileName.c_str(),std::ofstream::app);

	// set nSpace
	// additionally, for HRR if it's inside additional loop like ESP etc.
	// we need to consider add more nSpace
	int nSpace = 2;
	if (resultIntegralHasAdditionalOffset(infor.getOper())) {
		nSpace += 2;
	}
	
	// some comments
	myfile << endl;
	string line = "/************************************************************";
	printLine(nSpace,line,myfile);
	line = " * initilize the HRR steps : build the AB/CD variables";
	printLine(nSpace,line,myfile);
	line = " ************************************************************/";
	printLine(nSpace,line,myfile);

	// is it AB or CD side?
	if (side == BRA1 || side == BRA2 || side == BRA) {
		line = "Double ABX = A[0] - B[0];";
		printLine(nSpace,line,myfile);
		line = "Double ABY = A[1] - B[1];";
		printLine(nSpace,line,myfile);
		line = "Double ABZ = A[2] - B[2];";
		printLine(nSpace,line,myfile);
	}else{
		line = "Double CDX = C[0] - D[0];";
		printLine(nSpace,line,myfile);
		line = "Double CDY = C[1] - D[1];";
		printLine(nSpace,line,myfile);
		line = "Double CDZ = C[2] - D[2];";
		printLine(nSpace,line,myfile);
	}
	myfile << endl;
	myfile.close();
}

