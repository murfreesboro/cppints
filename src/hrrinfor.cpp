//
// CPPINTS: A C++ Program to Generate Analytical Integrals Based on Gaussian
// Form Primitive Functions
//
// Copyright (C) 2015 The State University of New York at Buffalo
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
#include<algorithm>
#include "printing.h"
#include "integral.h"
#include "inttype.h"
#include "shell.h"
#include "basis.h"
#include "shellsymbol.h"
#include "sqintsinfor.h"
#include "rr.h"
#include "boost/lexical_cast.hpp"
#include <boost/algorithm/string.hpp>   // string handling
#include "vrrinfor.h"
using namespace printing;
using namespace inttype;
using namespace integral;
using namespace shell;
using namespace basis;
using namespace sqintsinfor;
using namespace rr;
using boost::lexical_cast;
using namespace boost;
using namespace vrrinfor;

void HRRInfor::declareArray() const 
{
	// if it's not in file split, just return
	// we do not need the array declare
	if (! hrrFileSplit) return;

	// set the nSpace
	int nSpace = 2;
	if (resultIntegralHasAdditionalOffset(oper)) {
		nSpace += 2;
	}

	//
	// If we do not have last section, which means
	// HRR is the beginning and end section for all
	// of the cpp file, then we do not need declare
	//
	if (nextSection == NULL_POS) {
		return;
	}

	// now open the file
	string varFileName = infor.getWorkFuncName(false,section);
	ofstream varfile;
	varfile.open (varFileName.c_str(),std::ofstream::app);
	varfile << endl;
	string line;
	line = "/************************************************************";
	printLine(nSpace,line,varfile);
	line = " * declare the HRR result shell quartets";
	printLine(nSpace,line,varfile);
	line = " ************************************************************/";
	printLine(nSpace,line,varfile);

	// now print out all of shell quartets
	string arrayType = infor.getArrayType();
	for(int iSQ=0; iSQ<(int)outputSQList.size(); iSQ++) {
		const ShellQuartet& sq = outputSQList[iSQ];
		int status = outputSQStatus[iSQ];
		if (! inArrayStatus(status)) continue;
		int nInts = outputSQIntNumList[iSQ];
		string arrayName = sq.formArrayName(rrType);
		string nLHSInts  = lexical_cast<string>(nInts);
		string declare   = infor.getArrayDeclare(nLHSInts);
		line = arrayType + arrayName + declare;
		printLine(nSpace,line,varfile);
	}

	// now close the file
	varfile << endl;
	varfile.close();

				string line = "Double " + name + " = 0.0E0;";
				printLine(nSpace,line,myfile);
			}
		}
	}

	//
	// we need to declare the VRR results in case that VRR contraction is split
	//
	if (vrrContSplit) {

		// print head
		myfile << endl;
		string line = "//";
		printLine(nSpace,line,myfile);
		line = "// declare VRR results in array form, for the case that VRR and ";
		printLine(nSpace,line,myfile);
		line = "// contraction is split out, so contraction will be done after VRR ";
		printLine(nSpace,line,myfile);
		line = "//";
		printLine(nSpace,line,myfile);

		//
		// now print out the output rr sq list
		// they will be passed to the next module
		//
		for(int iSQ=0; iSQ<(int)vrrSQList.size(); iSQ++) {
			const ShellQuartet& sq  = vrrSQList[iSQ];
			if (sq.isSTypeSQ()) {
				Integral I(sq,0);
				string name = I.formVarName(VRR);
				string line = "Double " + name + " = 0.0E0;";
				printLine(nSpace,line,myfile);
			}else{
				const set<int>& intList = solvedIntList[iSQ];
				string name      = sq.formArrayName(VRR);
				string arrayType = getArrayType();
				int nInts        = intList.size();
				string declare   = getArrayDeclare(lexical_cast<string>(nInts));
				string line      = arrayType + name + declare;
				printLine(nSpace,line,myfile);
			}
		}
	}

	//
	// finally, for VRR split into multiple files; we need to declare
	// the array form of input/output shell quartets for the VRR
	// sub functions
	//
	// we note, that this only applies to the case that has more than
	// one sub file record. For the one sub file record case, the 
	// function does not have input (bottom integrals direclty computed),
	// and the output is just the VRR results vrrSQList (in vrrContSplit),
	// or the outputSQList (if contraction is done inside the function). 
	//
	if (subFilesList.size() > 1) {
		for(int iSubFile=0; iSubFile<(int)subFilesList.size()-1; iSubFile++) {
			const vector<ShellQuartet>& output = subFilesList[iSubFile].getLHSSQList();
			const vector<int>& status = subFilesList[iSubFile].getLHSSQStatus();
			for(int iSQ=0; iSQ<(int)output.size(); iSQ++) {
				const ShellQuartet& sq = output[iSQ];
				if(status[iSQ] != FUNC_INOUT_SQ) continue;

				// let's see whether it has been already declared
				bool byPass = false;
				vector<ShellQuartet>::const_iterator it = find(vrrSQList.begin(),vrrSQList.end(),sq);
				if (it != vrrSQList.end()) byPass = true;

				// now print
				if (! byPass) {
					if (sq.isSTypeSQ()) {
						Integral I(sq,0);
						string name = I.formVarName(VRR);
						string line = "Double " + name + " = 0.0E0;";
						printLine(nSpace,line,myfile);
					}else{
						int nInts        = subFilesList[iSubFile].getLHSSQIntNum(sq);
						string name      = sq.formArrayName(VRR);
						string arrayType = getArrayType();
						string declare   = getArrayDeclare(lexical_cast<string>(nInts));
						string line      = arrayType + name + declare;
						printLine(nSpace,line,myfile);
					}
				}
			}
		}
	}

	// now close the result statment part
	myfile << endl;
}

void VRRInfor::printVRRHead(const SQIntsInfor& infor) const
{
	// derive the name of VRR head
	string name = infor.getWorkFuncName(false,VRR_HEAD);

	// create the file 
	ofstream file;
	file.open(name.c_str(),std::ofstream::out);

	// firstly we need to see that wether we need a loop on the top
	if (oper == ESP) {
		string line = "// loop over grid points ";
		printLine(2,line,file);
		line = "for(UInt iGrid=0; iGrid<nGrids; iGrid++) {";
		printLine(2,line,file);
		file << endl;
	}

	// let's check that whether the operator is with error function
	// form, which is, operator is erf(r12)/r12
	if (withErf(oper)) {
		string line = "// check that whether we use erf(r12)/r12 form operator ";
		printLine(2,line,file);
		line = "bool withErfR12 = false;";
		printLine(2,line,file);
		line = "if (fabs(omega)>THRESHOLD_MATH) withErfR12 = true;";
		printLine(2,line,file);
		file << endl;
	}

	//
	// we print the vrr declaration variables here
	//
	printResultStatement(file);

	// set the nSpace
	int nSpace = 2;
	if (resultIntegralHasAdditionalOffset(oper) && ! vrrInFileSplit) {
		nSpace += 2;
	}

	// set up significance test if the file employs fmt function
	// however, if it's all bottom integrals we do not need to do it
	// just return true 
	if(sigCheck(oper) && ! isLastSection()) {
		string line = "// initialize the significance check for VRR part ";
		printLine(nSpace,line,file);
		line = "// this will determine that whether we skip the following part ";
		printLine(nSpace,line,file);
		line = "bool isSignificant = false;";
		printLine(nSpace,line,file);
		file << endl;
	}

	// now go to the given function
	switch(oper) {
		case TWOBODYOVERLAP:
			printTwoBodyOverlapHead(file,infor);
			break;
		case THREEBODYOVERLAP:
			printThreeBodyOverlapHead(file,infor);
			break;
		case MOM:
			printMOMHead(file,infor);
			break;
		case NAI:
			printNAIHead(file,infor);
			break;
		case ESP:
			printESPHead(file,infor);
			break;
		case ERI:
			printERIHead(file,infor);
			break;
		case EXPR12:
			printEXPR12Head(file,infor);
			break;
		case KINETIC:
			printKineticHead(file,infor);
			break;
		case THREEBODYKI:
			printThreeBodyOverlapHead(file,infor);
			break;
		default:
			crash(true,"Illegal operator type in print head function");
			break;
	}

	// now close the file
	file.close();

	//
	// do we need to calculate the bottom integrals in a complicated way?
	// right now this is only demanded by the MOM integrals
	//
	if (! vrrInFileSplit && oper == MOM) {
		const vector<ShellQuartet>& inputSQList = infor.getInputSQList();
		printMOMBottomIntegrals(name,inputSQList);
	}
}

void VRRInfor::printTwoBodyOverlapHead(ofstream& file,const SQIntsInfor& infor) const 
{
	///////////////////////////////////////////////////////
	//     information for the given shell quartets      //
	///////////////////////////////////////////////////////
	bool hasRROnBRA1 = hasVRROnVar(BRA1);
	bool hasRROnBRA2 = hasVRROnVar(BRA2);
	bool hasRROnBRA  = hasRROnBRA1 || hasRROnBRA2;
	bool withExpFac  = infor.withExpFac();

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

	// for gradient calculation, we need the alpha and beta value
	if (hasRROnBRA || withExpFac) {
		line = "Double onedz = iexp[ip2];";
		printLine(4,line,file);
		if (withExpFac) {
			line = "Double zeta  = 1.0E0/onedz;";
			printLine(4,line,file);
			line = "Double zdiff = iexpdiff[ip2];";
			printLine(4,line,file);
			line = "Double alpha = 0.5E0*(zeta+zdiff);";
			printLine(4,line,file);
			line = "Double beta  = 0.5E0*(zeta-zdiff);";
			printLine(4,line,file);
		}
	}

	// if we do RR on BRA side, we need 
	// them
	if (hasRROnBRA) {
		line = "Double oned2z= 0.5E0*onedz;";
		printLine(4,line,file);
		line = "UInt offsetP = 3*ip2;";
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
		line = "Double I_TWOBODYOVERLAP_S_S_vrr = ic2*fbra;";
		printLine(4,line,file);
	}else{
		line = "Double I_TWOBODYOVERLAP_S_S_vrr = fbra;";
		printLine(4,line,file);
	}

	// whether we pass it
	line = "if(fabs(I_TWOBODYOVERLAP_S_S_vrr)<THRESHOLD_MATH) continue;";
	printLine(4,line,file);
	file << endl;
}

void VRRInfor::printMOMHead(ofstream& file,const SQIntsInfor& infor) const 
{
	///////////////////////////////////////////////////////
	//     information for the given shell quartets      //
	//     in default, since we need PCX etc.            //
	//     therefore we need the P point infor           //
	///////////////////////////////////////////////////////
	bool hasRROnBRA1 = hasVRROnVar(BRA1);
	bool hasRROnBRA2 = hasVRROnVar(BRA2);
	bool hasRROnBRA  = true;
	bool withExpFac  = infor.withExpFac();

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
	if (hasRROnBRA || withExpFac) {
		line = "Double onedz = iexp[ip2];";
		printLine(4,line,file);
		if (withExpFac) {
			line = "Double zeta  = 1.0E0/onedz;";
			printLine(4,line,file);
			line = "Double zdiff = iexpdiff[ip2];";
			printLine(4,line,file);
			line = "Double alpha = 0.5E0*(zeta+zdiff);";
			printLine(4,line,file);
			line = "Double beta  = 0.5E0*(zeta-zdiff);";
			printLine(4,line,file);
		}
		line = "Double oned2z= 0.5E0*onedz;";
		printLine(4,line,file);
		line = "UInt offsetP = 3*ip2;";
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
		line = "Double I_TWOBODYOVERLAP_S_S_vrr = ic2*fbra;";
		printLine(4,line,file);
	}else{
		line = "Double I_TWOBODYOVERLAP_S_S_vrr = fbra;";
		printLine(4,line,file);
	}

	// whether we pass it
	line = "if(fabs(I_TWOBODYOVERLAP_S_S_vrr)<THRESHOLD_MATH) continue;";
	printLine(4,line,file);
	file << endl;
}

void VRRInfor::printKineticHead(ofstream& file, const SQIntsInfor& infor) const 
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
		line = "UInt offsetP = 3*ip2;";
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
		line = "Double I_KINETIC_S_S_vrr = ic2*fbra*xi*(3.0E0-twoxi*AB2);";
		printLine(4,line,file);
		line = "Double I_TWOBODYOVERLAP_S_S_vrr = ic2*fbra;";
		printLine(4,line,file);
	}else{
		line = "Double I_KINETIC_S_S_vrr = fbra*xi*(3.0E0-twoxi*AB2);";
		printLine(4,line,file);
		line = "Double I_TWOBODYOVERLAP_S_S_vrr = fbra;";
		printLine(4,line,file);
	}

	// whether we pass it
	line = "if(fabs(I_KINETIC_S_S_vrr)<THRESHOLD_MATH) continue;";
	printLine(4,line,file);
	file << endl;
}

void VRRInfor::printNAIHead(ofstream& file,const SQIntsInfor& infor) const 
{
	///////////////////////////////////////////////////////
	//     information for the given shell quartets      //
	///////////////////////////////////////////////////////
	bool hasRROnBRA1 = hasVRROnVar(BRA1);
	bool hasRROnBRA2 = hasVRROnVar(BRA2);
	bool hasRROnBRA  = hasRROnBRA1 || hasRROnBRA2;
	int  maxLSum     = getMaxLSum();
	bool withExpFac  = infor.withExpFac();

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
	if (withExpFac) {
		line = "Double zeta  = rho;";
		printLine(4,line,file);
		line = "Double zdiff = iexpdiff[ip2];";
		printLine(4,line,file);
		line = "Double alpha = 0.5E0*(zeta+zdiff);";
		printLine(4,line,file);
		line = "Double beta  = 0.5E0*(zeta-zdiff);";
		printLine(4,line,file);
	}
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
	line = "UInt offsetP = 3*ip2;";
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

void VRRInfor::printESPHead(ofstream& file, const SQIntsInfor& infor) const 
{
	///////////////////////////////////////////////////////
	//     information for the given shell quartets      //
	///////////////////////////////////////////////////////
	bool hasRROnBRA1 = hasVRROnVar(BRA1);
	bool hasRROnBRA2 = hasVRROnVar(BRA2);
	bool hasRROnBRA  = hasRROnBRA1 || hasRROnBRA2;
	int  maxLSum     = getMaxLSum();
	bool withExpFac  = infor.withExpFac();

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
	if (withExpFac) {
		line = "Double zeta  = rho;";
		printLine(6,line,file);
		line = "Double zdiff = iexpdiff[ip2];";
		printLine(6,line,file);
		line = "Double alpha = 0.5E0*(zeta+zdiff);";
		printLine(6,line,file);
		line = "Double beta  = 0.5E0*(zeta-zdiff);";
		printLine(6,line,file);
	}
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
	line = "UInt offsetP = 3*ip2;";
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

void VRRInfor::printThreeBodyOverlapHead(ofstream& file,const SQIntsInfor& infor) const 
{
	///////////////////////////////////////////////////////
	//     information for the given shell quartets      //
	///////////////////////////////////////////////////////
	bool hasRROnBRA1 = hasVRROnVar(BRA1);
	bool hasRROnBRA2 = hasVRROnVar(BRA2);
	bool hasRROnKET1 = hasVRROnVar(KET1);
	bool hasRR       = hasRROnBRA1  || hasRROnBRA2 || hasRROnKET1;
	bool withExpFac  = infor.withExpFac();

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

	// with exp fac?
	if (withExpFac) {
		line = "Double zeta  = 1.0E0/onedz;";
		printLine(4,line,file);
		line = "Double zdiff = iexpdiff[ip2];";
		printLine(4,line,file);
		line = "Double alpha = 0.5E0*(zeta+zdiff);";
		printLine(4,line,file);
		line = "Double beta  = 0.5E0*(zeta-zdiff);";
		printLine(4,line,file);
	}

	// prefactors for the SSSS integrals
	line = "Double fbra  = ifac[ip2];";
	printLine(4,line,file);

	// we need the new center even for the SSS integral
	line = "UInt offsetP = 3*ip2;";
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
	if (withExpFac) {
		line = "Double gamma = 1.0E0/onede;";
		printLine(6,line,file);
	}
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
		line = "Double I_THREEBODYOVERLAP_S_S_S_vrr = p*fbra*exp(-rho*PC2);";
		printLine(6,line,file);
	}else{
		line = "Double I_THREEBODYOVERLAP_S_S_S_vrr = ic2*jc2*p*fbra*exp(-rho*PC2);";
		printLine(6,line,file);
	}
	line = "if(fabs(I_THREEBODYOVERLAP_S_S_S_vrr)<THRESHOLD_MATH) continue;";
	printLine(6,line,file);
	file << endl;
}

void VRRInfor::printERIHead(ofstream& file, const SQIntsInfor& infor) const 
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
	bool withExpFac  = infor.withExpFac();

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

	// with exp fac?
	if (withExpFac) {
		line = "Double zeta  = 1.0E0/onedz;";
		printLine(4,line,file);
		line = "Double zdiff = iexpdiff[ip2];";
		printLine(4,line,file);
		line = "Double alpha = 0.5E0*(zeta+zdiff);";
		printLine(4,line,file);
		line = "Double beta  = 0.5E0*(zeta-zdiff);";
		printLine(4,line,file);
	}

	// bra side coefficients
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
	line = "UInt offsetP = 3*ip2;";
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

	// with exp fac?
	if (withExpFac) {
		line = "Double eta   = 1.0E0/onede;";
		printLine(6,line,file);
		line = "Double ediff = jexpdiff[jp2];";
		printLine(6,line,file);
		line = "Double gamma = 0.5E0*(eta+ediff);";
		printLine(6,line,file);
		line = "Double delta = 0.5E0*(eta-ediff);";
		printLine(6,line,file);
	}

	// ket side coefficients
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

	// let's do testing here first
	// now test the significance of integrals
	fmtIntegralsTest(maxLSum,ERI,6,file);
	file << endl;

	// continue to generate variables 
	line = "UInt offsetQ  = 3*jp2;";
	printLine(6,line,file);
	line = "Double QX    = Q[offsetQ  ];";
	printLine(6,line,file);
	line = "Double QY    = Q[offsetQ+1];";
	printLine(6,line,file);
	line = "Double QZ    = Q[offsetQ+2];";
	printLine(6,line,file);
	line = "Double rho   = 1.0E0/(onedz+onede);";
	printLine(6,line,file);
	line = "Double sqrho = sqrt(rho);";
	printLine(6,line,file);
	line = "Double PQ2   = (PX-QX)*(PX-QX)+(PY-QY)*(PY-QY)+(PZ-QZ)*(PZ-QZ);";
	printLine(6,line,file);

	// here set up the u
	line = "Double u     = rho*PQ2;";
	printLine(6,line,file);
	line = "if (withErfR12) u = PQ2/(1.0E0+1.0E0/(omega*omega)+1.0E0/rho);";
	printLine(6,line,file);
	line = "Double squ   = sqrt(u);";
	printLine(6,line,file);

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

	// we may also need to correct the bottom integral
	// if in error function form
	setupErfPrefactors(maxLSum,ERI,6,file);
}

void VRRInfor::printEXPR12Head(ofstream& file, const SQIntsInfor& infor) const 
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
	bool withExpFac  = infor.withExpFac();

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

	// with exp fac?
	if (withExpFac) {
		line = "Double zeta  = 1.0E0/onedz;";
		printLine(4,line,file);
		line = "Double zdiff = iexpdiff[ip2];";
		printLine(4,line,file);
		line = "Double alpha = 0.5E0*(zeta+zdiff);";
		printLine(4,line,file);
		line = "Double beta  = 0.5E0*(zeta-zdiff);";
		printLine(4,line,file);
	}

	// bra side coefficients
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
	line = "UInt offsetP = 3*ip2;";
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

	// with exp fac?
	if (withExpFac) {
		line = "Double eta   = 1.0E0/onede;";
		printLine(6,line,file);
		line = "Double ediff = jexpdiff[jp2];";
		printLine(6,line,file);
		line = "Double gamma = 0.5E0*(eta+ediff);";
		printLine(6,line,file);
		line = "Double delta = 0.5E0*(eta-ediff);";
		printLine(6,line,file);
	}

	// ket side coefficients
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

	// continue to generate variables 
	line = "UInt offsetQ  = 3*jp2;";
	printLine(6,line,file);
	line = "Double QX    = Q[offsetQ  ];";
	printLine(6,line,file);
	line = "Double QY    = Q[offsetQ+1];";
	printLine(6,line,file);
	line = "Double QZ    = Q[offsetQ+2];";
	printLine(6,line,file);
	line = "Double rho   = 1.0E0/(onedz+onede);";
	printLine(6,line,file);

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
		line = "Double odorho= omega/(rho+omega);";
		printLine(6,line,file);
		line = "Double od2k  = oned2k*odorho;";
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
		line = "Double orhod2z2 = rhod2zsq*odorho;";
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
		line = "Double orhod2e2= rhod2esq*odorho;";
		printLine(6,line,file);
	}

	// now compute the bottom integral
	file << endl;
	line = "// now begin compute the bottom integral";
	line = "Double rhodorho= rho/(rho+omega);";
	printLine(6,line,file);
	file << endl;

	// here we need to do something for the debugging
	line = "// if operator is normalized, the bottom integral will";
	printLine(6,line,file);
	line = "// multiply (omega/PI)^{3/2}, this is only used when omega is very large";
	printLine(6,line,file);
	line = "// so that the exp(-omega*r12^2) becomes delta function, such multiplier";
	printLine(6,line,file);
	line = "// ensures that the bottom integral is not zero";
	printLine(6,line,file);
	line = "#ifdef  DEBUG_EXPR12";
	printLine(0,line,file);
	line = "prefactor = prefactor*pow(1.0E0/PI,1.5E0)*pow(rho*(omega/(rho+omega)),1.5E0);";
	printLine(6,line,file);
	line = "#else";
	printLine(0,line,file);
	line = "prefactor = prefactor*pow(rhodorho,1.5E0);";
	printLine(6,line,file);
	line = "#endif";
	printLine(0,line,file);

	// continue the bottom integral calculation
	line = "Double PQ2     = (PX-QX)*(PX-QX)+(PY-QY)*(PY-QY)+(PZ-QZ)*(PZ-QZ);";
	printLine(6,line,file);
	line = "Double expFac  = exp(-omega*rhodorho*PQ2);";
	printLine(6,line,file);
	line = "Double I_EXPR12_S_S_S_S_vrr = prefactor*expFac;";
	printLine(6,line,file);
	line = "if(fabs(I_EXPR12_S_S_S_S_vrr)<THRESHOLD_MATH) continue;";
	printLine(6,line,file);
	file << endl;
}

//////////////////////////////////////////////////////////////////////////
//           @@@@ printing VRR contraction part of code                 //
//////////////////////////////////////////////////////////////////////////
void VRRInfor::contraction(const SQIntsInfor& infor, 
		const vector<ShellQuartet>& sqlist, const int& fileIndex) const
{
	// let's go to see whether this is for contraction on composite shells
	bool comSQ = infor.isComSQ();

	// determine that how many space should be given for each line printing
	int nSpace = getNSpaceByOper(oper);

	// detect that wether we have compilcated offset for result?
	// we note, that the nInts should be total number of integrals
	// in terms of final results, so we use infor to give the right number
	string additionalOffset;
	bool hasAdditionalOffset = resultIntegralHasAdditionalOffset(oper);
	if (hasAdditionalOffset) {
		int nInts = infor.nInts();
		additionalOffset = determineAdditionalOffset(oper,nInts);
	}

	// also consider the nSpace
	// for the non file split mode, we need 
	// to increment the nSpace
	if (hasAdditionalOffset) nSpace += 2;
	if (vrrInFileSplit || vrrContSplit) nSpace = 2;

	// get the file name
	// for contraction and VRR together, the module name is VRR
	// else it's VRR_CONT
	int module = VRR;
	if(vrrContSplit) module = VRR_CONT;
	string filename = infor.getWorkFuncName(false,module,fileIndex);

	// create the RR file
	// this is already in a mod of a+
	ofstream myfile;
	myfile.open(filename.c_str(),std::ofstream::app);

	// for VRR and contraction split case,
	// here we will transform the input sq into the variable
	// form for contraction
	if (vrrContSplit) {
		for(int iSQ=0; iSQ<(int)sqlist.size(); iSQ++) {

			// get the corresponding sq in vrr result
			const ShellQuartet& sq = sqlist[iSQ];
			int pos = -1;
			for(int iSQ2=0; iSQ2<(int)vrrSQList.size(); iSQ2++) {
				const ShellQuartet& sq2 = vrrSQList[iSQ2];
				if (sq2 == sq) {
					pos = iSQ2;
					break;
				}
			}

			// check pos
			if (pos < 0) {
				crash(true,"something wrong in VRR contraction working function");
			}

			// now let's print
			const set<int>& intList = solvedIntList[pos];

			//
			// now print out the comment for this section's 
			// contraction
			//
			myfile << endl;
			string line = "/************************************************************";
			printLine(nSpace,line,myfile);
			line = " * transform the array form of integral into variables: " + sq.getName();
			printLine(nSpace,line,myfile);
			line = " ************************************************************/";
			printLine(nSpace,line,myfile);

			// now work begins
			pos  = -1; 
			string arrayName = sq.formArrayName(VRR);
			for(set<int>::const_iterator it=intList.begin(); it != intList.end(); ++it) {

				// form the lhs
				pos++;
				string rhs = arrayName + "[" + boost::lexical_cast<string>(pos) + "]";

				//
				// LHS name
				// lhs will be the VRR result possibily with modifiers
				// such as division and exponent infor
				//
				int intIndex = *it;
				Integral I(sq,intIndex);
				string lhs = I.formVarName(VRR);

				// form the code
				string line = lhs + " = " + rhs + ";";
				printLine(nSpace,line,myfile);
			}
		}
	}

	// set up some vector, to hold the output list
	vector<ShellQuartet> outputList;
	outputList.reserve(100);

	//
	// now loop over the input sq list
	// which is also the VRR's output results
	//
	for(int iSQ=0; iSQ<(int)outputSQList.size(); iSQ++) {

		// see whether the sq is in the sqlist?
		// the sqlist is the output of the given sub file
		// we note, that it's without the multipler information
		bool hasIt = false;
		const ShellQuartet& sq  = outputSQList[iSQ];
		ShellQuartet newSQ(sq);
		newSQ.destroyMultipliers();
		for(int iSQ2=0; iSQ2<(int)sqlist.size(); iSQ2++) {
			const ShellQuartet& sq2 = sqlist[iSQ2];
			if (sq2 == newSQ) {
				hasIt = true;
				break;
			}
		}
		if(! hasIt) continue;

		// update the output list
		if (vrrContSplit) {
			outputList.push_back(sq);
		}

		// get the sq and it's corresponding integral list
		const set<int>& intList = outputIntList[iSQ];
		int status = outputSQStatus[iSQ];
		int nInts = intList.size();
		int nTotalInts = sq.getNInts();
		int diff = nTotalInts - nInts;

		// whether the LHS is array or not?
		// whether the LHS is final result?
		bool lhsUseArray = inArrayStatus(status);
		bool isResult    = isGlobalResult(status);

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
		line = " * totally " + lexical_cast<string>(diff) + " integrals are omitted ";
		printLine(nSpace,line,myfile);
		line = " ************************************************************/";
		printLine(nSpace,line,myfile);

		// according to the shell quartet, get the coe 
		// offset
		int ic2Offset = -1; 
		int jc2Offset = -1;
		if (comSQ) {
			infor.getCoeOffset(sq,ic2Offset,jc2Offset);
			if (ic2Offset == -1 && jc2Offset == -1) {
				crash(true, "incorrect getCoeOffset result");
			}
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

		// let's form the code for coefs
		// adding in the composite shell case if possible
		string coefsName = sq.getName() + "_coefs";
		string coefsCodeLHS = "Double " + coefsName + " = ";
		string coefsCodeRHS;
		if (comSQ) {
			coefsCodeRHS = coe;
		}

		// do we have exponential factors add in?
		if (sq.withExpFac()) {
			string expFactors = sq.getExpFacMultiplers();
			if (comSQ) {
				coefsCodeRHS = coefsCodeRHS + "*" + expFactors;
			}else{
				coefsCodeRHS = expFactors;
			}
		}

		// now print out the whole coefs
		bool withModifier = comSQ || sq.withExpFac();
		if (withModifier) {
			line = coefsCodeLHS + coefsCodeRHS + ";";
			printLine(nSpace,line,myfile);
		}

		// now step into contraction
		// pos is the position for the integral in the array form
		int pos  = -1; 
		for(set<int>::const_iterator it=intList.begin(); it != intList.end(); ++it) {

			// get the intIndex
			// also form the position for the integral
			// in the array
			int intIndex = *it;
			pos++;

			//
			// LHS name
			// lhs will be the VRR result possibily with modifiers
			// such as division and exponent infor
			//
			Integral I(sq,intIndex);
			string lhsName = sq.getName();
			if (! lhsUseArray) {
				lhsName = I.getName();
			}
			if (isResult) {
				lhsName = "abcd";
			}

			// index for lhs
			string lhs;
			if (isResult || lhsUseArray) {

				// determine the offset
				int offset = pos;
				if (isResult) {
					offset = infor.getOffset(sq,pos);
				}

				// counting the additional offset
				// if we apply additional offset to the VRR result,
				// it must be the final results
				string lhsIndex;
				if (hasAdditionalOffset && isResult) {
					lhsIndex  = "[" + additionalOffset + "+" + lexical_cast<string>(offset) + "]";
				}else{
					lhsIndex  = "[" + lexical_cast<string>(offset) + "]";
				}
				lhs = lhsName + lhsIndex;
			}else{
				lhs = lhsName;
			}


			//
			// RHS name
			// for rhs, we need to destroy the modifier information
			// for both composite shell quartet(division)
			// as well as the exponential factors
			//
			// Because the VRR always uses the variable form, therefore
			// in the normal digestion, the rhs is always in variable rather
			// than array form
			//
			I.destroyMultipliers();
			string rhs = I.formVarName(VRR);
			if (withModifier) rhs = coefsName + "*" + rhs; 

			// form the code
			string line = lhs + " += " + rhs + ";";
			printLine(nSpace,line,myfile);
		}
	}

	// now close the whole file
	myfile.close();

	// finally for the contraction and VRR split case,
	// we need to generate a contraction function statement
	if (vrrContSplit) {

		// let's get the argument
		string arg;
		arg.reserve(1000);

		// this is the composite shell coefficients
		if (comSQ) {
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

		// add in the exponential factors
		if (infor.withExpFac()) {
			int nBody = getOperOrder(oper);
			if (nBody == 1) {
				arg = arg + "const Double& alpha, ";
			}else if (nBody == 2) {
				arg = arg + "const Double& alpha, const Double& beta, ";
			}else if (nBody == 3) {
				arg = arg + "const Double& alpha, const Double& beta, const Double& gamma, ";
			}else {
				arg = arg + "const Double& alpha, const Double& beta, const Double& gamma, const Double& delta, ";
			}
		}

		// now it's input shell quartet
		for(int iSQ=0; iSQ<(int)sqlist.size(); iSQ++) {
			const ShellQuartet& sq = sqlist[iSQ];
			string line = "const Double* " + sq.formArrayName(VRR);
			line = line + ", ";
			arg = arg + line;
		}

		// finally it's output
		bool hasResultSQ = false;
		for(int iSQ=0; iSQ<(int)outputList.size(); iSQ++) {
			const ShellQuartet& sq = outputList[iSQ];
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

		// now form the function statement line
		string name = infor.getWorkFuncName(true,VRR_CONT,fileIndex);
		string func = "void " + name + "(" + arg + ");";

		// open file stream, write the result
		string f    = infor.getWorkFuncName(false,VRR_CONT_STATEMENT);
		ofstream file;
		file.open(f.c_str(),std::ofstream::app);

		// add some comments
		file << endl;
		string line = "// VRR contraction function: " + 
			boost::lexical_cast<string>(fileIndex);
		printLine(0,line,file);
		printLine(0,func,file);
		file << endl;

		// now close file
		file.close();
	}
}

void VRRInfor::vrrContraction(const SQIntsInfor& infor) const
{
	// now deal with the VRR contraction split case
	if(vrrContSplit) {

		// initialization
		int fileIndex = 1;
		int nParam = 0;
		int nLHS   = 0;
		vector<ShellQuartet> sqlist;
		sqlist.reserve(200);

		// now let's loop over the VRR result
		sqlist.clear();
		for(int iSQ=0; iSQ<(int)vrrSQList.size(); iSQ++) {

			// now add in sqlist
			const ShellQuartet& sq = vrrSQList[iSQ];
			sqlist.push_back(sq);
			nParam += 1;
			const set<int>& intList = solvedIntList[iSQ];
			nLHS += intList.size();

			// let's see how many output sq corresponding to this one
			for(int iSQ2=0; iSQ2<(int)outputSQList.size(); iSQ2++) {
				const ShellQuartet& sq2  = outputSQList[iSQ2];
				ShellQuartet newSQ(sq2);
				newSQ.destroyMultipliers();
				if (sq == newSQ) {
					nParam += 1;
					nLHS += intList.size();
				}
			}

			// let's see whether we begin to print it
			if (nLHS>nLHSForVRRSplit || nParam>maxParaForFunction || iSQ == int(vrrSQList.size()-1)) {

				// print the contraction file
				contraction(infor,sqlist,fileIndex);

				// now clear result
				sqlist.clear();
				nLHS = 0;
				nParam = 0;
				fileIndex += 1;
			}
		}

		// now finish work, let's return
		return;
	}

	// now this is the case that VRR and contraction are doing together,
	// however VRR is split into different parts 
	if(vrrInFileSplit) {

		// set up shell quartet vector
		vector<ShellQuartet> sqlist;
		sqlist.reserve(200);

		// do contraction in terms of the VRR part
		for(int iSubFile=0; iSubFile<(int)subFilesList.size(); iSubFile++) {

			// got the section output shell quartets for this sub file
			sqlist.clear();
			const vector<ShellQuartet>& lhs = subFilesList[iSubFile].getLHSSQList();
			for(int iSQ=0; iSQ<(int)lhs.size(); iSQ++) {

				// whether this sq has been already been considered
				const ShellQuartet& sq = lhs[iSQ];
				vector<ShellQuartet>::const_iterator it2 = find(sqlist.begin(),sqlist.end(),sq);
				if (it2 != sqlist.end()) continue;

				// now check this output sq is the local module result?
				vector<ShellQuartet>::const_iterator it = find(vrrSQList.begin(),vrrSQList.end(),sq);
				if (it != vrrSQList.end()) {
						sqlist.push_back(sq);
				}
			}

			// if this sub file does not have any VRR module output, just by pass
			if (sqlist.size() == 0) continue;

			// file index
			int subFileIndex = iSubFile + 1;

			// print contraction to the given sub file
			contraction(infor,sqlist,subFileIndex);
		}

		// now return
		return;
	}

	// now it's not VRR split, nor VRR contraction split,
	// we just do normal contraction
	int fileIndex = -1;
	contraction(infor,vrrSQList,fileIndex);
}

//////////////////////////////////////////////////////////////////////////
//                     @@@@ form sub files for VRR                      //
//////////////////////////////////////////////////////////////////////////
int VRRInfor::contrationCount(const ShellQuartet& sq) const
{
	// does this sq appears in the VRR result?
	vector<ShellQuartet>::const_iterator it = find(vrrSQList.begin(),vrrSQList.end(),sq);
	if (it == vrrSQList.end()) {
		return 0;
	}
	
	// now let's count
	int nCont = 0;
	for(int iSQ=0; iSQ<(int)outputSQList.size(); iSQ++) {
		const ShellQuartet& sq2  = outputSQList[iSQ];
		ShellQuartet newSQ(sq2);
		newSQ.destroyMultipliers();
		if (sq == newSQ) {
			const set<int>& intList = outputIntList[iSQ];
			nCont += intList.size();
		}
	}

	// now return 
	return nCont;
}

void VRRInfor::formSubFiles(bool onlyOneSubFile, const SQIntsInfor& infor, const RR& vrr) 
{
	// set up a working copy of sub file record
	SubFileRecord record(VRR);
	record.init();
	int nLHS = 0;

	// for operator using fmt function, we are able to count
	// the number of function parameters as a simulation
	int nCurrentMSQ  = 0;
	int nPreviousMSQ = 0;
	int mVal = -1;

	// for operator using fmt, we can also estimate the 
	// local contractions. This is the contraction work
	// to summarize the calculation into the function output
	// array
	int nLocalCont   = 0;

	//
	// here we will form the sub file according to the printing order of rrsq list
	// basically, it's reverse order of rrsqlist
	//
	// here we can not consider the function parameter number for each sub file
	// this is because the sub file input/output can be only formed when we have
	// all of sub files formed
	//
	// only an estimation is made for oper using fmt function
	//
	const list<RRSQ>& rrsqList = vrr.getRRSQList();
	for(list<RRSQ>::const_reverse_iterator it=rrsqList.rbegin(); it!=rrsqList.rend(); ++it) {

		// add this RRSQ
		record.updateRRSQ(*it);

		// count the LHS 
		const list<int>& LHS = it->getLHSIndexArray();
		const ShellQuartet& sq = it->getLHSSQ();
		nLHS += LHS.size();

		// if it's required that with only one sub file
		// then we omit all of following steps
		if (onlyOneSubFile) continue;

		// we count in the local contraction for the 
		// shell quartets with m > 0
		// m == 0 possibly could be the VRR results
		// we only count in local contraction for current M
		if (useFmt(oper) && sq.getM() > 0) {
			nLocalCont += LHS.size();
		}

		// let's see the m value situation
		// we can do this is because all of sq shares the 
		// same M value are grouped together
		if (useFmt(oper)) {
			// is this the first M value defined?
			// else we may do some statistics
			int m = sq.getM();
			if (mVal<0) {
				mVal = m;
				nCurrentMSQ = 1;
			}else{
				// do we have m value changed?
				if (mVal != m) {
					// exchange the previous sq counting with current one
					// also reset mVal
					// restart the current m value sq counting
					mVal = m;
					nPreviousMSQ = nCurrentMSQ;
					nCurrentMSQ  = 1;

					// also reset the local contraction work
					// when M value changes
					// we try to count all of LHS for current M value
					// this is estimation of local contraction
					nLocalCont   = 0;
				}else{
					nCurrentMSQ  = nCurrentMSQ + 1;
				}
			}
		}

		// do we have contraction work for this rrsq?
		// for intetgrals with fmt we have a simple trick
		// because only mvalue = 0 are those VRR outputs
		bool countingContraction = true;
		if (useFmt(oper)) {
			int m = sq.getM();
			if (m>0) countingContraction = false;
		}

		// count in the possible contraction work
		if (countingContraction) {
			if (vrrContSplit) {
				vector<ShellQuartet>::const_iterator it = find(vrrSQList.begin(),vrrSQList.end(),sq);
				if (it != vrrSQList.end()) {
					nLHS += LHS.size();
				}
			}else{
				const ShellQuartet& sq = it->getLHSSQ();
				int nCont = contrationCount(sq);
				nLHS += nCont;
			}
		}

		// for oper using fmt, let's see it's estimation of 
		// function parameters
		// for ERI part, the subroutine will approximately
		// use previous M value shell quartets as well as 
		// current M value shell quartets, if the subroutine
		// is not the first one
		//
		// the first one will take all of bottom integrals,
		// so it's not using the previous M value sq counting 
		int nFuncParam = 0;
		if (useFmt(oper)) {
			if (subFilesList.size() == 0) {
				nFuncParam = nCurrentMSQ;
			}else{
				nFuncParam = nPreviousMSQ + nCurrentMSQ;
			}
		}

		// do we reach the sub file limit?
		// if so we add it to the record list,
		// and everything restarts
		if (nLHS+nLocalCont>nLHSForVRRSplit || nFuncParam>maxParaForFunction) {
			subFilesList.push_back(record);
			record.clear();
			nLHS = 0;
		}
	}

	// we need to add in the last record
	// if it has some content inside
	// this is for the case that we reach
	// the end of rrsqlist
	if (nLHS>0) {
		subFilesList.push_back(record);
	}

	// before forming the input/output shell quartets for the sub file list,
	// we need to know which one is the module output
	// it's only when we do VRR contraction split we need to identify them
	// because in this case we will do contraction in different files
	if (vrrContSplit) {
		for(int iSub=0; iSub<(int)subFilesList.size(); iSub++) {
			SubFileRecord& subFile = subFilesList[iSub];
			subFile.updateModuleOutput(infor,vrrSQList);
		}
	}

	// we need to form the input/output parameter 
	// for each sub file
	for(int iSub=0; iSub<(int)subFilesList.size(); iSub++) {
		SubFileRecord& record1 = subFilesList[iSub];
		for(int jSub=iSub+1; jSub<(int)subFilesList.size(); jSub++) {
			SubFileRecord& record2 = subFilesList[jSub];

			// now update the output for record 1
			record1.updateOutput(record2);

			// also update the input for record 2
			record2.updateInput(record1);
		}
	}

	// now let's print out the function statement
	for(int iSub=0; iSub<(int)subFilesList.size(); iSub++) {

		// get the function parameter
		const SubFileRecord& record = subFilesList[iSub];
		string arg = getVRRArgList(iSub,infor,record);

		// now form the function statement line
		int fileIndex = iSub + 1;
		string name = infor.getWorkFuncName(true,VRR,fileIndex);
		string func = "void " + name + "(" + arg + ");";

		// open file stream, write the result
		string f    = infor.getWorkFuncName(false,VRR_FUNC_STATEMENT);
		ofstream file;
		file.open(f.c_str(),std::ofstream::app);

		// add some comments
		file << endl;
		string line = "// VRR working function: " + 
			boost::lexical_cast<string>(fileIndex);
		printLine(0,line,file);
		printLine(0,func,file);
		file << endl;

		// now close file
		file.close();
	}
}

string VRRInfor::getVRRArgList(const int& subFileIndex, 
		const SQIntsInfor& infor, const SubFileRecord& record) const
{

	// get the status of the VRR    
	bool hasRROnBRA1 = hasVRROnVar(BRA1);
	bool hasRROnBRA2 = hasVRROnVar(BRA2);
	bool hasRROnKET1 = hasVRROnVar(KET1);
	bool hasRROnKET2 = hasVRROnVar(KET2);
	bool hasRROnBRA  = hasRROnBRA1 || hasRROnBRA2;
	bool hasRROnKET  = hasRROnKET1 || hasRROnKET2;
	int  maxLSum     = getMaxLSum();
	int  intOperator = oper;
	bool withExpFac  = infor.withExpFac();

	// reserve space for the output
	string arg;
	arg.reserve(2000);

	// coefficients is the first part 
	// we note, it's only that this is composite shell quartet,
	// we need to bring the coefficients into the subroutine
	// so that to form contraction in the end stage of VRR
	// else the contraction of coefficients is done in the 
	// pure contraction function
	//
	// also when VRR contraction is split into different
	// work, the coefficients is only used for contraction
	// therefore we will handled it over there
	//
	if (! vrrContSplit) {
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

		// for gradient calcualtion, we also need alpha and beta
		if (withExpFac) {
			arg = arg + "const Double& alpha, const Double& beta, ";
		}

		// bottom integral
		arg = arg + "const Double& I_TWOBODYOVERLAP_S_S_vrr, ";
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

		// for gradient calcualtion, we also need alpha and beta
		if (withExpFac) {
			arg = arg + "const Double& alpha, const Double& beta, ";
		}

		// bottom integral
		arg = arg + "const Double& I_TWOBODYOVERLAP_S_S_vrr, const Double& I_KINETIC_S_S_vrr, ";
	}

	// three body overlap
	if (intOperator == THREEBODYOVERLAP || intOperator == THREEBODYKI) {

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

		// for gradient calcualtion, we also need alpha and beta etc.
		if (withExpFac) {
			arg = arg + "const Double& alpha, const Double& beta, const Double& gamma, ";
		}

		// bottom integral
		arg = arg + "const Double& I_THREEBODYOVERLAP_S_S_S_vrr, ";
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

		// for gradient calcualtion, we also need alpha and beta
		if (withExpFac) {
			arg = arg + "const Double& alpha, const Double& beta, ";
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

		// for gradient calcualtion, we also need alpha and beta
		if (withExpFac) {
			arg = arg + "const Double& alpha, const Double& beta, ";
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

		// for gradient calcualtion, we also need alpha and beta
		if (withExpFac) {
			arg = arg + "const Double& alpha, const Double& beta, const Double& gamma, const Double& delta, ";
		}

		// bottom integral
		int m1 = -1;
		int m2 = -1;
		record.getMValueLimit(m1,m2);
		for(int m=m1; m<=m2; m++) {
			string name = getBottomIntName(m,intOperator);
			arg = arg + "const Double& " + name + ", ";
		}
	}

	// expr12
	if (intOperator == EXPR12) {

		// VRR variable
		arg = arg + "const Double& od2k, ";

		if (hasRROnBRA) {
			arg = arg + "const Double& oned2z, const Double& orhod2z2, ";
			arg = arg + "const Double& WPX, const Double& WPY, const Double& WPZ, ";
		}
		if (hasRROnKET) {
			arg = arg + "const Double& oned2e, const Double& orhod2e2, ";
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

		// for gradient calcualtion, we also need alpha and beta
		if (withExpFac) {
			arg = arg + "const Double& alpha, const Double& beta, const Double& gamma, const Double& delta, ";
		}

		// bottom integral
		arg = arg + "const Double& I_EXPR12_S_S_S_S_vrr, ";
	}

	// now let's consider the input shell quartets
	// it's only avaiable for the 2ed or higher sub file record
	// since the first one will use the bottom integrals
	if (subFileIndex > 0) {
		const vector<ShellQuartet>& rhs = record.getRHSSQList();
		const vector<int>&    rhsStatus = record.getRHSSQStatus();
		for(int iSQ=0; iSQ<(int)rhs.size(); iSQ++) {
			if (rhsStatus[iSQ] != FUNC_INOUT_SQ) continue;
			const ShellQuartet& sq = rhs[iSQ];
			string line = "const Double* " + sq.formArrayName(VRR) + ", ";
			arg = arg + line;
		}
	}

	// now it's the output shell quartets
	// firstly we need to consider the output as function output
	const vector<ShellQuartet>& lhs = record.getLHSSQList();
	const vector<int>&    lhsStatus = record.getLHSSQStatus();
	for(int iSQ=0; iSQ<(int)lhs.size(); iSQ++) {
		if (lhsStatus[iSQ] != FUNC_INOUT_SQ) continue;
		const ShellQuartet& sq = lhs[iSQ];
		string line = "Double* " + sq.formArrayName(VRR) + ", ";
		arg = arg + line;
	}

	// for the case that VRR and contraction are doing together
	// we need to pass the VRR module output to the function, too
	if (! vrrContSplit) {

		// let's see whether we have VRR module result in terms of
		// function output, because we already count in the sq in
		// FUNC_INOUT_SQ above, we will not do them here
		vector<ShellQuartet> sqlist;
		sqlist.reserve(200);
		hasABCD = false;
		for(int iSQ=0; iSQ<(int)lhs.size(); iSQ++) {
			if (lhsStatus[iSQ] == GLOBAL_RESULT_SQ) hasABCD = true;
			if (lhsStatus[iSQ] != FUNC_INOUT_SQ) continue;
			const ShellQuartet& sq = lhs[iSQ];
			vector<ShellQuartet>::const_iterator it = find(vrrSQList.begin(),vrrSQList.end(),sq);
			if (it == vrrSQList.end()) continue;
			sqlist.push_back(sq);
		}

		// now let's get the corresponding output list
		if (sqlist.size()>0) {

			// get the output sq list
			vector<ShellQuartet> outputList;
			outputList.reserve(200);
			for(int iSQ=0; iSQ<(int)outputSQList.size(); iSQ++) {
				const ShellQuartet& sq2  = outputSQList[iSQ];
				ShellQuartet newSQ(sq2);
				newSQ.destroyMultipliers();
				vector<ShellQuartet>::const_iterator it = find(sqlist.begin(),sqlist.end(),newSQ);
				if (it == sqlist.end()) continue;
				outputList.push_back(sq2);
			}

			// now let's print it
			for(int iSQ=0; iSQ<(int)outputList.size(); iSQ++) {
				const ShellQuartet& sq = outputList[iSQ];
				string line = "Double* " + sq.getName() + ", ";
				arg = arg + line;
			}
		}

		// now let's see whether we have global result
		if (hasABCD) {
			arg = arg + "Double* abcd";
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

	// basically these are the input stuff
	return arg;
}

void VRRInfor::subFilesForming(const SQIntsInfor& infor, const RR& vrr)
{
	// let's count how many LHS for VRR and contraction part in total
	int nLHSCon = 0;
	int nVRRLHS = vrr.countLHSIntNumbers();
	for(int iSQ=0; iSQ<(int)outputSQList.size(); iSQ++) {
		const set<int>& intList = outputIntList[iSQ];
		nLHSCon += intList.size();
	}

	// do we have VRR split? 
	int nLHS = nVRRLHS + nLHSCon;
	if (nLHS>=nLHSForVRRSplit) {
		vrrInFileSplit = true;
	}else{
		vrrInFileSplit = false;
	}

	// do we have VRR contraction split?
	// for integrals with we need to also count int the number 
	// of bottom integrals
	//
	// we only do VRR contraction split when VRR is in spliting
	//
	vrrContSplit = false;
	int nOutputSQ = outputSQList.size();
	if (useFmt(oper)) {
		int maxAng = getMaxLSum() + 1;
		nOutputSQ += maxAng;
	}
	if (vrrInFileSplit && nOutputSQ>maxParaForFunction) {
		vrrContSplit = true;
	}

	// do we have only one sub file?
	bool onlyOneSubFile = false;
	if (vrrInFileSplit) {
		if (nLHS<nLHSForVRRSplit*(1+VRRSplitCoefs)) onlyOneSubFile = true;
	}

	// now let's form the sub files
	if (vrrInFileSplit) {
		formSubFiles(onlyOneSubFile,infor,vrr);
	}
}

//////////////////////////////////////////////////////////////////////////
//                     @@@@ contructors etc.                            //
//////////////////////////////////////////////////////////////////////////
VRRInfor::VRRInfor(const SQIntsInfor& infor, const RR& vrr):Infor(infor),vrrInFileSplit(false),
	vrrContSplit(false),nextSection(infor.nextSection(VRR)),oper(infor.getOper()),
	vrrSQList(vrr.getRRResultSQList()),solvedIntList(vrr.getRRUnsolvedIntList()),
	outputSQList(vrr.getRRResultSQList()),outputIntList(vrr.getRRUnsolvedIntList()),
	outputSQStatus(outputSQList.size(),VARIABLE_SQ)  
{
	vrr.updateSQIntListForVRR(vrrSQList,solvedIntList); 
	for(int iSQ=0; iSQ<(int)outputSQList.size(); iSQ++) {
		const ShellQuartet& sq = outputSQList[iSQ];
		if (sq.isSTypeSQ()) {
			outputSQStatus[iSQ] = BOTTOM_SQ;
		}else if (infor.isResult(sq)) {
			outputSQStatus[iSQ] = GLOBAL_RESULT_SQ;
		}
	}
}

void VRRInfor::updateOutputSQInArray(const vector<ShellQuartet>& moduleInput)
{
	for(int iSQ=0; iSQ<(int)outputSQList.size(); iSQ++) {
		const ShellQuartet& sq = outputSQList[iSQ];
		if (outputSQStatus[iSQ] != VARIABLE_SQ) continue;
		vector<ShellQuartet>::const_iterator it = find(moduleInput.begin(),moduleInput.end(),sq);
		if (it != moduleInput.end()) {
			outputSQStatus[iSQ] = ARRAY_SQ;
		}
	}
}