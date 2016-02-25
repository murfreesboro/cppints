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
#include "boost/lexical_cast.hpp"
#include <boost/algorithm/string.hpp>   // string handling
#include <boost/filesystem.hpp>
#include "inttype.h"
#include "derivinfor.h"
#include "rr.h"
#include "nonrr.h"
#include "printing.h"
#include "vrrinfor.h"
#include "hrrinfor.h"
#include "nonrrinfor.h"
#include "sqints.h"
using boost::lexical_cast;
using namespace boost::filesystem;
using namespace boost;
using namespace derivinfor;
using namespace printing;
using namespace inttype;
using namespace rr;
using namespace nonrr;
using namespace vrrinfor;
using namespace hrrinfor;
using namespace nonrrinfor;
using namespace sqints;

void SQInts::intCodeGeneration()
{
	// this is the output shell quartet list and unsolved 
	// integral list passing from one code section to another
	vector<ShellQuartet> outputSQList; 
	vector<set<int> > unsolvedList; 

	// initilize the output sq list and unsolved integral list 
	// for deriv job, we initilize it from the derivSQList
	// else we will do it from inputSQList in SQIntsInfor 
	if (infor.getJobOrder() > 0) {
		outputSQList = infor.getDerivSQList();
	}else{
		outputSQList = infor.getInputSQList();
	}
	for(int iSQ=0; iSQ<(int)outputSQList.size(); iSQ++) {
		set<int> list;
		outputSQList[iSQ].getIntegralList(list);
		unsolvedList.push_back(list);
	}

	///////////////////////////////////////////////////////////////////////
	// some note before we proceed the real work:                        //
	// the code generation is in reverse order of the real cpp file.     //
	// the result is generated first, then to it's previous section;     //
	// finally the generation will stop at the VRR section               //
	//                                                                   //
	//      %%%             CODE GENERATION STEP                         //
	///////////////////////////////////////////////////////////////////////

	////////////////////////////
	//   derivative section   //
	////////////////////////////
	NONRR derivJob(outputSQList,unsolvedList);
	size_t nLHSDeriv = 0;
	if (infor.getJobOrder() > 0) {

		// generate the integrals
		derivJob.buildRRSQList();

		// add the deriv section to infor
		int moduleName = DERIV;
		infor.appendCodeSection(moduleName);

		// let's carry the shell quartets to the next section
		// here the unsolved list must be in integral index
		// form, so do it before array index transformation
		outputSQList.clear();
		unsolvedList.clear();
		outputSQList = derivJob.getBottomSQList();
		unsolvedList = derivJob.getBottomIntList();

		// update nLHS
		nLHSDeriv = (size_t)derivJob.countLHSIntNumbers();
	}
	NONRRInfor derivJobInfor(infor,derivJob);
	infor.updateWithArray(derivJobInfor.fileSplit());

	////////////////////////////
	// possible NON-RR section//
	////////////////////////////
	NONRR nonRRJob(outputSQList,unsolvedList);
	size_t nLHSNonRR = 0;
	if (isNONRROper(infor.getOper())) {

		// build the rrsq
		nonRRJob.buildRRSQList();

		// let's carry the shell quartets to the next section
		// here the unsolved list must be in integral index
		// form, so do it before array index transformation
		outputSQList.clear();
		unsolvedList.clear();
		outputSQList = nonRRJob.getBottomSQList();
		unsolvedList = nonRRJob.getBottomIntList();

		// add the non-RR section to infor
		int moduleName = NON_RR;
		infor.appendCodeSection(moduleName);

		// update nLHS
		nLHSNonRR  = (size_t)nonRRJob.countLHSIntNumbers();
	}
	NONRRInfor nonRRJobInfor(infor,nonRRJob);
	infor.updateWithArray(nonRRJobInfor.fileSplit());

	////////////////////////////
	// second HRR section     //
	////////////////////////////
	RR hrr2(HRR2,HRR,outputSQList,unsolvedList);
	size_t nLHSHRR2 = 0;
	int firstSide  = NULL_POS;
	int secondSide = NULL_POS;
	if (infor.hasHRR()) {

		// determine the first and second side
		hrr2.sideDeterminationInHRR(firstSide,secondSide);

		// now do the second side HRR
		// for generated code, second side is performed after
		// first side; so in generation of code second side 
		// is prior to the first side
		if (secondSide != NULL_POS) {

			// generat rrsq
			hrr2.generateRRSQList(secondSide);

			// let's carry the shell quartets to the next section
			// here the unsolved list must be in integral index
			// form, so do it before array index transformation
			outputSQList.clear();
			unsolvedList.clear();
			outputSQList = hrr2.getHRRBottomSQList();
			unsolvedList = hrr2.getHRRBottomIntList();

			// add the HRR2 section to infor
			int moduleName = HRR2;
			infor.appendCodeSection(moduleName);

			// update nLHS
			nLHSHRR2   = (size_t)hrr2.countLHSIntNumbers();
		}
	}
	HRRInfor HRR2JobInfor(infor,hrr2);
	infor.updateWithArray(HRR2JobInfor.fileSplit());

	////////////////////////////
	//  first HRR section     //
	////////////////////////////
	RR hrr1(HRR1,HRR,outputSQList,unsolvedList);
	size_t nLHSHRR1 = 0;
	if (firstSide != NULL_POS) {

		// now do the first side
		hrr1.generateRRSQList(firstSide);

		// now rewrite the shell quartet list
		outputSQList.clear();
		unsolvedList.clear();
		outputSQList = hrr1.getHRRBottomSQList();
		unsolvedList = hrr1.getHRRBottomIntList();

		// add the HRR1 section to infor
		int moduleName = HRR1;
		infor.appendCodeSection(moduleName);

		// update nLHS
		nLHSHRR1   = (size_t)hrr1.countLHSIntNumbers();
	}
	HRRInfor HRR1JobInfor(infor,hrr1);
	infor.updateWithArray(HRR1JobInfor.fileSplit());

	////////////////////////////
	//      VRR section       //
	////////////////////////////
	RR vrr(VRR,infor.getVRRMethod(),outputSQList,unsolvedList);
	vrr.generateRRSQList(NULL_POS);

	// update nLHS
	size_t nLHSVRR = (size_t)vrr.countLHSIntNumbers();

	// add the VRR section to infor
	int moduleName = VRR;
	infor.appendCodeSection(moduleName);

	// now build the vrr infor
	VRRInfor vrrInfor(infor,vrr);
	infor.updateWithArray(vrrInfor.fileSplit());

	///////////////////////////////////////////////////////////////////////
	// now let's analyze the situation of file split. All modules has    //
	// been set up, and now it's time to see whether we need to make     //
	// some changes                                                      //
	//      %%%%            FILE SPLIT DETERMINATION STEP                //
	///////////////////////////////////////////////////////////////////////

	// do we need to re-evalate the file split status
	// for all of codes?
	// we do this when the main driver cpp file becomes too large
	size_t nLHSMainCPP = 0;
	const vector<int>& secInfor = infor.getSectionInfor(); 
	for(int iSec=0; iSec<(int)secInfor.size(); iSec++) {
		int sec = secInfor[iSec];
		if (sec == DERIV) {
			if (! derivJobInfor.fileSplit()) {
				nLHSMainCPP = nLHSMainCPP + nLHSDeriv;
			}
		}else if (sec == NON_RR) {
			if (! nonRRJobInfor.fileSplit()) {
				nLHSMainCPP = nLHSMainCPP + nLHSNonRR;
			}
		}else if (sec == HRR2) {
			if (! HRR2JobInfor.fileSplit()) {
				nLHSMainCPP = nLHSMainCPP + nLHSHRR2;
			}
		}else if (sec == HRR1) {
			if (! HRR1JobInfor.fileSplit()) {
				nLHSMainCPP = nLHSMainCPP + nLHSHRR1;
			}
		}else if (sec == VRR) {
			if (! vrrInfor.fileSplit()) {
				nLHSMainCPP = nLHSMainCPP + nLHSVRR;
			}
		}
	}

	// do we need to split the single cpp file?
	bool doSplit = false;
	double nTotalLHSCoef = infor.totalLHSScaleFac;
	int nTotalLHSLimit   = infor.nTotalLHSLimit;
	if (nLHSMainCPP>(size_t)(nTotalLHSLimit*(1+nTotalLHSCoef))) {
		doSplit = true;
		infor.updateWithArray(true);
	}

	// now let's update the infor for sections
	// from the last section to VRR
	if (doSplit) {
		size_t nLHS = nLHSMainCPP;
		for(int iSec=0; iSec<(int)secInfor.size(); iSec++) {

			// now let's update
			int sec = secInfor[iSec];
			if (sec == DERIV) {
				if (! derivJobInfor.fileSplit()) {
					derivJobInfor.updateFileSplit();
					nLHS = nLHS - nLHSDeriv;
				}
			}else if (sec == NON_RR) {
				if (! nonRRJobInfor.fileSplit()) {
					nonRRJobInfor.updateFileSplit();
					nLHS = nLHS - nLHSNonRR;
				}
			}else if (sec == HRR2) {
				if (! HRR2JobInfor.fileSplit()) {
					HRR2JobInfor.updateFileSplit();
					nLHS = nLHS - nLHSHRR2;
				}
			}else if (sec == HRR1) {
				if (! HRR1JobInfor.fileSplit()) {
					HRR1JobInfor.updateFileSplit();
					nLHS = nLHS - nLHSHRR1;
				}
			}else if (sec == VRR) {
				if (! vrrInfor.fileSplit()) {
					vrrInfor.updateFileSplit();
					nLHS = nLHS - nLHSVRR;
				}
			}

			// do we get to a break?
			if (nLHS<(size_t)(nTotalLHSLimit*(1+nTotalLHSCoef))) {
				break;
			}
		}
	}

	///////////////////////////////////////////////////////////////////////
	//  based on the file split status, we can figure out the input      //
	//  and output shell quartet status, in array or var?                //
	//    %%%%    SET UP MODULE INPUT/OUPUT STATUS                       //
	///////////////////////////////////////////////////////////////////////

	// firstly for VRR, update it's output when it's not in file split
	if (! vrrInfor.fileSplit()) {
		for(int iSec=0; iSec<(int)secInfor.size(); iSec++) {
			int sec = secInfor[iSec];
			if (sec == VRR) break;
			if (sec == DERIV && derivJobInfor.fileSplit()) {
				const vector<ShellQuartet>& output = derivJobInfor.getOutputSQList();
				vrrInfor.updateOutputSQInArray(output);
			}else if (sec == NON_RR && nonRRJobInfor.fileSplit()) {
				const vector<ShellQuartet>& output = nonRRJobInfor.getOutputSQList();
				vrrInfor.updateOutputSQInArray(output);
			}else if (sec == HRR2 && HRR2JobInfor.fileSplit()) {
				const vector<ShellQuartet>& output = HRR2JobInfor.getOutputSQList();
				vrrInfor.updateOutputSQInArray(output);
			}else if (sec == HRR1 && HRR1JobInfor.fileSplit()) {
				const vector<ShellQuartet>& output = HRR1JobInfor.getOutputSQList();
				vrrInfor.updateOutputSQInArray(output);
			}
		}
	}else{
		vrrInfor.subFilesForming(infor,vrr);
	}

	// form output sq which are in array status
	vector<ShellQuartet> sqlist;
	sqlist.reserve(500);
	const vector<ShellQuartet>& vrrOutputSQ = vrrInfor.getOutputSQList();
	const vector<int>& vrrOutputSQStatus    = vrrInfor.getOutputSQStatus();
	for(int iSQ=0; iSQ<(int)vrrOutputSQ.size(); iSQ++) {
		if (inArrayStatus(vrrOutputSQStatus[iSQ])) {
			sqlist.push_back(vrrOutputSQ[iSQ]);
		}
	}

	// now VRR output sq has fixed, reversely update all of other modules
	// which has VRR module results. We note, that we only do it to these modules
	// when they are not in file split status 
	if (sqlist.size() > 0) {
		for(int iSec=0; iSec<(int)secInfor.size(); iSec++) {
			int sec = secInfor[iSec];
			if (sec == VRR) break;
			if (sec == DERIV && ! derivJobInfor.fileSplit()) {
				derivJobInfor.updateInputSQInArray(sqlist);
			}else if (sec == NON_RR && ! nonRRJobInfor.fileSplit()) {
				nonRRJobInfor.updateInputSQInArray(sqlist);
			}else if (sec == HRR2 && ! HRR2JobInfor.fileSplit()) {
				HRR2JobInfor.updateInputSQInArray(sqlist);
			}else if (sec == HRR1 && ! HRR1JobInfor.fileSplit()) {
				HRR1JobInfor.updateInputSQInArray(sqlist);
			}
		}
	}

	// now let's do HRR1, we may not have HRR
	if (infor.hasSection(HRR1)) {

		// update the module output from all of following modules
		if (! HRR1JobInfor.fileSplit()) {
			for(int iSec=0; iSec<(int)secInfor.size(); iSec++) {
				int sec = secInfor[iSec];
				if (sec == HRR1) break;
				if (sec == DERIV && derivJobInfor.fileSplit()) {
					const vector<ShellQuartet>& output = derivJobInfor.getOutputSQList();
					HRR1JobInfor.updateOutputSQInArray(output);
				}else if (sec == NON_RR && nonRRJobInfor.fileSplit()) {
					const vector<ShellQuartet>& output = nonRRJobInfor.getOutputSQList();
					HRR1JobInfor.updateOutputSQInArray(output);
				}else if (sec == HRR2 && HRR2JobInfor.fileSplit()) {
					const vector<ShellQuartet>& output = HRR2JobInfor.getOutputSQList();
					HRR1JobInfor.updateOutputSQInArray(output);
				}
			}
		}else{
			HRR1JobInfor.formSubFiles(infor,hrr1);
		}

		// form output sq which are in array status
		sqlist.clear();
		const vector<ShellQuartet>& outputSQ = HRR1JobInfor.getOutputSQList();
		const vector<int>& outputSQStatus    = HRR1JobInfor.getOutputSQStatus();
		for(int iSQ=0; iSQ<(int)outputSQ.size(); iSQ++) {
			if (inArrayStatus(outputSQStatus[iSQ])) {
				sqlist.push_back(outputSQ[iSQ]);
			}
		}

		// now HRR1 output sq has fixed, reversely update all of other modules
		// which has HRR1 module results. We note, that we only do it to these modules
		// when they are not in file split status 
		for(int iSec=0; iSec<(int)secInfor.size(); iSec++) {
			int sec = secInfor[iSec];
			if (sec == HRR1) break;
			if (sec == DERIV && ! derivJobInfor.fileSplit()) {
				derivJobInfor.updateInputSQInArray(sqlist);
			}else if (sec == NON_RR && ! nonRRJobInfor.fileSplit()) {
				nonRRJobInfor.updateInputSQInArray(sqlist);
			}else if (sec == HRR2 && ! HRR2JobInfor.fileSplit()) {
				HRR2JobInfor.updateInputSQInArray(sqlist);
			}
		}
	}

	// now let's do HRR2, we may not have HRR2
	if (infor.hasSection(HRR2)) {

		// update the module output from all of following modules
		if (! HRR2JobInfor.fileSplit()) {
			for(int iSec=0; iSec<(int)secInfor.size(); iSec++) {
				int sec = secInfor[iSec];
				if (sec == HRR2) break;
				if (sec == DERIV && derivJobInfor.fileSplit()) {
					const vector<ShellQuartet>& output = derivJobInfor.getOutputSQList();
					HRR2JobInfor.updateOutputSQInArray(output);
				}else if (sec == NON_RR && nonRRJobInfor.fileSplit()) {
					const vector<ShellQuartet>& output = nonRRJobInfor.getOutputSQList();
					HRR2JobInfor.updateOutputSQInArray(output);
				}
			}
		}else{
			HRR2JobInfor.formSubFiles(infor,hrr2);
		}

		// form output sq which are in array status
		sqlist.clear();
		const vector<ShellQuartet>& outputSQ = HRR2JobInfor.getOutputSQList();
		const vector<int>& outputSQStatus    = HRR2JobInfor.getOutputSQStatus();
		for(int iSQ=0; iSQ<(int)outputSQ.size(); iSQ++) {
			if (inArrayStatus(outputSQStatus[iSQ])) {
				sqlist.push_back(outputSQ[iSQ]);
			}
		}

		// now HRR2 output sq has fixed, reversely update all of other modules
		// which has HRR2 module results. We note, that we only do it to these modules
		// when they are not in file split status 
		if (sqlist.size() > 0) {
			for(int iSec=0; iSec<(int)secInfor.size(); iSec++) {
				int sec = secInfor[iSec];
				if (sec == HRR2) break;
				if (sec == DERIV && ! derivJobInfor.fileSplit()) {
					derivJobInfor.updateInputSQInArray(sqlist);
				}else if (sec == NON_RR && ! nonRRJobInfor.fileSplit()) {
					nonRRJobInfor.updateInputSQInArray(sqlist);
				}
			}
		}
	}

	// now let's do NON-RR
	if (infor.hasSection(NON_RR)) {

		// update the module output for DERIV module
		if (! nonRRJobInfor.fileSplit()) {
			if (infor.hasSection(DERIV) && derivJobInfor.fileSplit()) {
				const vector<ShellQuartet>& output = derivJobInfor.getOutputSQList();
				nonRRJobInfor.updateOutputSQInArray(output);
			}
		}else{
			nonRRJobInfor.formSubFiles(infor,nonRRJob);
		}

		// form output sq which are in array status
		sqlist.clear();
		const vector<ShellQuartet>& outputSQ = nonRRJobInfor.getOutputSQList();
		const vector<int>& outputSQStatus    = nonRRJobInfor.getOutputSQStatus();
		for(int iSQ=0; iSQ<(int)outputSQ.size(); iSQ++) {
			if (inArrayStatus(outputSQStatus[iSQ])) {
				sqlist.push_back(outputSQ[iSQ]);
			}
		}

		// we only need to consider the DERIV module
		if (sqlist.size() > 0) {
			if (infor.hasSection(DERIV) && ! derivJobInfor.fileSplit()) {
				derivJobInfor.updateInputSQInArray(sqlist);
			}
		}
	}

	// for deriv job we do not need to update the input and output
	// all of input shell quartet updating are finished
	// and output sq must be global results
	if (infor.hasSection(DERIV) && derivJobInfor.fileSplit()) {
		derivJobInfor.formSubFiles(infor,derivJob);
	}

	///////////////////////////////////////////////////////////////////////
	//               %%%%     PRINT OUT THE CODES                        //
	///////////////////////////////////////////////////////////////////////
	vrr.vrrPrint(infor,vrrInfor);
	if (infor.hasSection(HRR1)) {
		hrr1.hrrPrint(infor,HRR1JobInfor);
	}
	if (infor.hasSection(HRR2)) {
		hrr2.hrrPrint(infor,HRR2JobInfor);
	}
	if (infor.hasSection(NON_RR)) {
		nonRRJob.print(infor,nonRRJobInfor);
	}
	if (infor.hasSection(DERIV)) {
		derivJob.print(infor,derivJobInfor);
	}
}

void SQInts::assembleCPPFiles() const
{
	// open the top cpp file
	string cppFile = infor.getWorkFuncName(false,NULL_POS,-1,true);
	ofstream CPP;
	CPP.open(cppFile.c_str());

	// generate the head of file
	infor.headPrinting(CPP);

	// now append all of function prototype to the head
	// particularly, because the vrr contraction may also
	// in a seperate function, therefore we need to check it
	vector<int> modules(infor.getSectionInfor());
	for(int iFunc=modules.size()-1; iFunc>=0; iFunc--) {
		int module = modules[iFunc];
		appendFuncPrototype(module,CPP);
	}
	if (hasFileDefined(VRR_CONT)) {
		appendFuncPrototype(VRR_CONT,CPP);
	}

	// now it's the function name
	string func = infor.getFuncName();
	func = "void " + func;
	string arg  = infor.getArgList();
	string line = func + "(" + arg + ")";
	CPP << line << endl;
	CPP << "{" << endl;

	//
	// now we begin to generate code starting from VRR 
	// section, until the final section is met
	//

	//////////////////////////////
	//       VRR section        //
	//////////////////////////////

	// generate the VRR head
	// append it to the main cpp file
	int moduleName = VRR_HEAD;
	appendFile(moduleName,CPP);

	// now it's the main body of VRR
	moduleName = VRR;
	if (hasFileDefined(moduleName)) {
		formWorkFile(moduleName);
		formFunctionCall(moduleName,CPP);
	}else{
		appendFile(moduleName,CPP);
	}

	// contraction part
	// for contraction if it's not in file split
	// mode, it will be with VRR code
	moduleName = VRR_CONT;
	if (hasFileDefined(moduleName)) {
		formWorkFile(moduleName);
		formFunctionCall(moduleName,CPP);
	}

	// now let's close the VRR part
	// determine the number space for printing etc.
	int oper = infor.getOper();
	int nSpace = getNSpaceByOper(oper);
	int nSpaceStop = 2;
	if (resultIntegralHasAdditionalOffset(oper)) nSpaceStop += 2;

	// finally, we need to add braket closure to the vrr body
	line = "}";
	for(int iSpace= nSpace-2; iSpace>=nSpaceStop; iSpace = iSpace - 2) {
		printLine(iSpace,line,CPP);
	}

	// if we have significance test, then we may
	// test the significance first
	// if VRR step all of integrals are omitted,
	// then we do not need to do HRR accordingly
	if(sigCheck(oper) && ! infor.isLastSection(VRR)) {
		CPP << endl;
		line = "/************************************************************";
		printLine(nSpaceStop,line,CPP);
		line = " * let's see the significance test result. if VRR result is";
		printLine(nSpaceStop,line,CPP);
		line = " * insignificant, there's no need to do following codes";
		printLine(nSpaceStop,line,CPP);
		line = " ************************************************************/";
		printLine(nSpaceStop,line,CPP);

		// the handling of significance at the bottom of VRR 
		// this is to see whether we do HRR part
		// for the case that VRR/HRR are in a loop 
		// we can not return, just use continue
		line = "if (! isSignificant) return;";
		if (resultIntegralHasAdditionalOffset(oper)) {
			line = "if (! isSignificant) continue;";
		}
		printLine(nSpaceStop,line,CPP);
	}	

	//////////////////////////////
	//       HRR1 section etc.  //
	//////////////////////////////
	for(int iFunc=modules.size()-1; iFunc>=0; iFunc--) {
		int module = modules[iFunc];
		if (module == VRR || module == VRR_CONT) continue;
		if (hasFileDefined(module)) {

			// form the work file
			formWorkFile(module);

			// append the main module file if we have
			// this contains the variables declare or 
			// the results declare etc.
			appendFile(module,CPP,false);

			// form the function call to the result cpp file
			formFunctionCall(module,CPP);
		}else{
			appendFile(module,CPP);
		}
	}

	// now finalize the cpp file
	if (oper == ESP) {
		CPP << "  }" << endl;
		CPP << "}" << endl;
	}else{
		CPP << "}" << endl;
	}
	CPP.close();
}

bool SQInts::isFileExist() const 
{
	string fileName = infor.getWorkFuncName(false,NULL_POS,-1,true);
	path cppFile(fileName.c_str());
	if (exists(cppFile)) {
		return true;
	}
	return false;
}

bool SQInts::hasFileDefined(int moduleName) const
{
	// we will check t he function statement
	// if it has, then we have the files for this section
	int module = -1;
	if (moduleName == VRR) {
		module = VRR_FUNC_STATEMENT;
	}else if (moduleName == VRR_CONT) {
		module = VRR_CONT_STATEMENT;
	}else if (moduleName == HRR1) {
		module = HRR1_FUNC_STATEMENT;
	}else if (moduleName == HRR2) {
		module = HRR2_FUNC_STATEMENT;
	}else if (moduleName == NON_RR) {
		module = NON_RR_FUNC_STATEMENT;
	}else if (moduleName == DERIV) {
		module = DERIV_FUNC_STATEMENT;
	}else {
		crash(true,"the input module name is invalid in SQInts::hasFileDefined");
	}
	bool onlyWithFuncName = false;
	bool inFinalDir = false;
	string pro = infor.getWorkFuncName(onlyWithFuncName,module,-1,inFinalDir);
	path p(pro.c_str());
	if (exists(p)) {
		return true;
	}
	return false;
}

void SQInts::appendFuncPrototype(int moduleName, ofstream& CPP) const
{
	// deriv the function prototype
	// get the name
	int module = -1;
	if (moduleName == VRR) {
		module = VRR_FUNC_STATEMENT;
	}else if (moduleName == VRR_CONT) {
		module = VRR_CONT_STATEMENT;
	}else if (moduleName == HRR1) {
		module = HRR1_FUNC_STATEMENT;
	}else if (moduleName == HRR2) {
		module = HRR2_FUNC_STATEMENT;
	}else if (moduleName == NON_RR) {
		module = NON_RR_FUNC_STATEMENT;
	}else if (moduleName == DERIV) {
		module = DERIV_FUNC_STATEMENT;
	}else {
		crash(true,"the input module name is invalid in SQInts::appendFuncPrototype");
	}

	// now append the prototype to the main cpp file
	// we note  that if the module is not in file split mode, we do not
	// have the prototype file
	// so we check it's existence
	bool onlyWithFuncName = false;
	bool inFinalDir = false;
	string pro = infor.getWorkFuncName(onlyWithFuncName,module,-1,inFinalDir);
	path p(pro.c_str());
	if (! exists(p)) {
		return;
	}

	// now do the work
	ifstream PRO;
	PRO.open(pro.c_str(),ios::in);
	string line;
	while(getline(PRO,line)) {
		CPP << line << endl;
	}
	PRO.close();
}

void SQInts::formFunctionCall(int moduleName, ofstream& CPP) const
{
	// deriv the function prototype
	// get the name
	int module = -1;
	if (moduleName == VRR) {
		module = VRR_FUNC_STATEMENT;
	}else if (moduleName == VRR_CONT) {
		module = VRR_CONT_STATEMENT;
	}else if (moduleName == HRR1) {
		module = HRR1_FUNC_STATEMENT;
	}else if (moduleName == HRR2) {
		module = HRR2_FUNC_STATEMENT;
	}else if (moduleName == NON_RR) {
		module = NON_RR_FUNC_STATEMENT;
	}else if (moduleName == DERIV) {
		module = DERIV_FUNC_STATEMENT;
	}else {
		crash(true,"the input module name is invalid in SQInts::appendFuncPrototype");
	}

	// now append the prototype to the main cpp file
	// we note  that if the module is not in file split mode, we do not
	// have the prototype file
	// so we check it's existence
	bool onlyWithFuncName = false;
	bool inFinalDir = false;
	string pro = infor.getWorkFuncName(onlyWithFuncName,module,-1,inFinalDir);
	path p(pro.c_str());
	if (! exists(p)) {
		cout << "the module name " << moduleName << endl;
		crash(true, "the prototype file does not exist in formFunctionCall of sqints");
	}

	// now get the function statement
	vector<string> prototype;
	prototype.reserve(200);
	ifstream PRO;
	PRO.open(pro.c_str(),ios::in);
	string line;
	while(getline(PRO,line)) {

		// get the function statement
		if (line.find("void") == std::string::npos) continue;
		string statement = line;

		// now get rid of the type information etc.
		LineParse lp(statement);
		string result;
		result.reserve(1000);
		for(int i=0; i<lp.getNPieces(); i++) {

			// we begin to determine that whether this is variable name
			string val = lp.findValue(i);

			// check whether this is the function name
			if (val.find("(")!=std::string::npos) {
				int pos = val.find("(");
				string v = val.substr(0,pos+1);
				result = result + v;
				continue;
			}

			// omit these values
			if (val.find("const")!=std::string::npos) continue;
			if (val.find("vector")!=std::string::npos) continue;
			if (val.find("Double")!=std::string::npos) continue;
			if (val.find("LocalMemScr")!=std::string::npos) continue;
			if (val.find("void")!=std::string::npos) continue;

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
			}

			// append it to the result
			result = result + val;
		}
		prototype.push_back(result);
	}
	PRO.close();

	// determine the nspace
	int nSpace = 2;
	int oper = infor.getOper();
	if (moduleName == VRR || moduleName == VRR_CONT) {
		nSpace = getNSpaceByOper(oper);
	}
	if (resultIntegralHasAdditionalOffset(oper)) {
		nSpace += 2;
	}

	// finally, append these function calls to the main file
	for(int i=0; i<(int)prototype.size(); i++) {
		string func = prototype[i];
		printLine(nSpace,func,CPP);
		CPP << endl;
	}
	CPP << endl;
}

void SQInts::formWorkFile(int moduleName) const 
{
	// let's see how many sub files we have for this given module
	// except VRR, VRR only have one module
	int iFile = 1;
	int nFiles = 0;
	bool onlyWithFuncName = false;
	bool inFinalDir = false;
	while(true) {
		string fileName = infor.getWorkFuncName(onlyWithFuncName,moduleName,iFile,inFinalDir);
		path cppFile(fileName.c_str());
		if (exists(cppFile)) {
			nFiles++;
		}else{
			break;
		}
		iFile++;
	}

	// we should not have 0 sub files
	// because it already in file split mode
	if (nFiles == 0) {
		crash(true,"why there are no subfiles in SQInts::formWorkFile");
	}

	// deriv the function prototype
	// get the name
	int module = -1;
	if (moduleName == VRR) {
		module = VRR_FUNC_STATEMENT;
	}else if (moduleName == VRR_CONT) {
		module = VRR_CONT_STATEMENT;
	}else if (moduleName == HRR1) {
		module = HRR1_FUNC_STATEMENT;
	}else if (moduleName == HRR2) {
		module = HRR2_FUNC_STATEMENT;
	}else if (moduleName == NON_RR) {
		module = NON_RR_FUNC_STATEMENT;
	}else if (moduleName == DERIV) {
		module = DERIV_FUNC_STATEMENT;
	}else {
		crash(true,"the input module name is invalid in SQInts::formWorkFile");
	}

	// now let's deriv the prototype
	vector<string> prototype;
	prototype.reserve(nFiles);
	string pro = infor.getWorkFuncName(onlyWithFuncName,module,-1,inFinalDir);
	path p(pro.c_str());
	if (! exists(p)) {
		cout << "Missing function prototype file " << pro << endl;
		crash(true, "Why the function prototype file does not exist??? we need it in formWorkFile");
	}
	ifstream PRO;
	PRO.open(pro.c_str(),ios::in);
	string line;
	while(getline(PRO,line)) {
		if (line.find("void") == std::string::npos) continue;
		if (line.find(";") != std::string::npos) {
			int pos = line.find(";");
			line[pos] = ' ';
		}
		prototype.push_back(line);
	}
	PRO.close();

	// also check the prototype number
	if (prototype.size() != (size_t)nFiles) {
		cout << "module name " << moduleName << endl;
		cout << "prototype size " << prototype.size() << endl;
		cout << "nFiles " << nFiles << endl;
		crash(true,"the prototype number does not match nFiles in SQInts::formWorkFile");
	}

	// now let's form every subfile
	inFinalDir = true;
	for(iFile = 1; iFile<=nFiles; iFile++) {

		// get the work file name
		int fileIndex = iFile;
		string fileName = infor.getWorkFuncName(onlyWithFuncName,moduleName,fileIndex,inFinalDir);

		// now let's open it
		ofstream CPP;
		CPP.open(fileName.c_str(),ios::out);

		// firstly create the include part
		infor.headPrinting(CPP);

		// now print the function prototyoe
		string pro = prototype[iFile-1];
		CPP << pro << endl;

		// add in the "{" so that to start the main body
		CPP << "{" << endl;

		// now let's open the work file
		string work = infor.getWorkFuncName(onlyWithFuncName,moduleName,fileIndex,false);

		// check existance
		path p(work.c_str());
		if (! exists(p)) {
			cout << "Missing work file " << work << endl;
			crash(true, "Why the main body code file does not exist??? we need it in formWorkFile");
		}

		// now let's open it
		ifstream workFile;
		string line;
		workFile.open(work.c_str(),ios::in);
		while(getline(workFile,line)) {
			CPP << line << endl;
		}
		workFile.close();

		// close the cpp file
		CPP << "}" << endl;
		CPP.close();
	}
}

void SQInts::appendFile(int moduleName, ofstream& CPP, bool checkExist) const 
{
	// according to the module and status,
	// get the work file name
	string file = infor.getWorkFuncName(false,moduleName);

	// check existance
	// if we do not have this file, we do not do a crash
	// just return
	path p(file.c_str());
	if (! exists(p)) {
		if (checkExist) {
			cout << "Missing file " << file << endl;
			crash(true, "Why the file does not exist??? we need it in appendFile function of sqints");
		}else{
			return;
		}
	}

	// now let's open it
	ifstream workFile;
	workFile.open(file.c_str(),ios::in);
	string line;
	while(getline(workFile,line)) {
		CPP << line << endl;
	}
	workFile.close();
}

void SQInts::codeGeneration() 
{
	// generate the tmp folder
	// the work dir name should be same with
	// the one used in function of getProjectFileDir
	int intOperator = infor.getOper();
	int jobOrder    = infor.getJobOrder();
	string workDir = infor.getProjectTmpFileDir(intOperator,jobOrder);
	path tmpWorkDir(workDir.c_str());
	create_directory(tmpWorkDir);

	// generate the codes
	intCodeGeneration();

	// now assemble cpp files
	assembleCPPFiles();

	// this is debugging codes
	//const vector<int> shellCodes = infor.getShellCodeArray();
	//if (intOperator == THREEBODYKI && shellCodes[0] == 3 && shellCodes[1] == 3 && shellCodes[2] == 3) {
	//	crash(true, "job complete");
	//}

	// finally, clean the tmp file
	remove_all(tmpWorkDir);
}

