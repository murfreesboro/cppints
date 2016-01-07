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

	//
	// some note before we proceed the real work:
	// the code generation is in reverse order of the real cpp file.
	// the result is generated first, then to it's previous section;
	// finally the generation will stop at the VRR section
	//


	////////////////////////////
	//   derivative section   //
	////////////////////////////
	if (infor.getJobOrder() > 0) {

		// generate the integrals
		NONRR derivJob(outputSQList,unsolvedList);
		derivJob.buildRRSQList();

		// add the deriv section to infor
		int moduleName = DERIV;
		infor.appendCodeSection(moduleName);

		// let's carry the shell quartets to the next section
		// here the unsolved list must be in integral index
		// form, so do it before array index transformation
		outputSQList.clear();
		unsolvedList.clear();
		outputSQList = derivJob.getResultSQList();
		unsolvedList = derivJob.getResultIntSQList();

		// let's see whether we do it in file split mode?
		int nDerivLHS = (int)derivJob.countLHSIntNumbers();
		bool inFileSplit = false;
		if (nDerivLHS>infor.nLHSForDerivSplit) inFileSplit = true;
		infor.setFileSplitMode(moduleName,inFileSplit);
		if (infor.withArrayIndex(moduleName)) { 
			derivJob.arrayIndexTransformation(infor);
		}

		// generate the code
		derivJob.derivPrint(infor);

		// now let's update the VRR/HRR shell quartet list must be in array form
		if (infor.fileSplit(moduleName)) {
			infor.updateVRRSQListInArray(outputSQList); 
			infor.updateHRRSQListInArray(outputSQList); 
		}
	}

	////////////////////////////
	// possible NON-RR section//
	////////////////////////////
	if (isNONRROper(infor.getOper())) {

		// build the rrsq
		NONRR nonRRJob(outputSQList,unsolvedList);
		nonRRJob.buildRRSQList();

		// let's carry the shell quartets to the next section
		// here the unsolved list must be in integral index
		// form, so do it before array index transformation
		outputSQList.clear();
		unsolvedList.clear();
		outputSQList = nonRRJob.getResultSQList();
		unsolvedList = nonRRJob.getResultIntSQList();

		// add the non-RR section to infor
		int moduleName = NON_RR;
		infor.appendCodeSection(moduleName);

		// let's see whether we do it in file split mode?
		int nNonRRLHS = (int)nonRRJob.countLHSIntNumbers();
		bool inFileSplit = false;
		if (nNonRRLHS>infor.nLHSForNonRRSplit) inFileSplit = true;
		infor.setFileSplitMode(moduleName,inFileSplit);
		nonRRJob.arrayIndexTransformation(infor);

		// generate the code
		nonRRJob.nonRRPrint(infor);

		// now let's update the VRR/HRR shell quartet list must be in array form
		if (infor.fileSplit(moduleName)) {
			infor.updateVRRSQListInArray(outputSQList); 
			infor.updateHRRSQListInArray(outputSQList); 
		}
	}

	////////////////////////////
	//      HRR section       //
	////////////////////////////
	if (infor.hasHRR()) {

		// determine the first and second side
		RR hrr2(HRR2,HRR,outputSQList,unsolvedList);
		int firstSide  = NULL_POS;
		int secondSide = NULL_POS;
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

			// let's see whether we do it in file split mode?
			int nHRRLHS = (int)hrr2.countLHSIntNumbers();
			bool inFileSplit = false;
			if (nHRRLHS>infor.nLHSForHRR2Split) inFileSplit = true;
			infor.setFileSplitMode(moduleName,inFileSplit);
			hrr2.arrayIndexTransformation(infor);

			// print the code
			hrr2.hrrPrint(infor);

			// now let's update the VRR shell quartet list must be in array form
			// we will pick up those who can not do HRR work for the both sides
			if (infor.fileSplit(moduleName)) {
				infor.updateVRRSQListInArray(outputSQList); 
			}
		}

		// do you have work on the first side?
		if (firstSide != NULL_POS) {

			// now do the first side
			RR hrr1(HRR1,HRR,outputSQList,unsolvedList);
			hrr1.generateRRSQList(firstSide);

			// now rewrite the shell quartet list
			outputSQList.clear();
			unsolvedList.clear();
			outputSQList = hrr1.getHRRBottomSQList();
			unsolvedList = hrr1.getHRRBottomIntList();

			// add the HRR1 section to infor
			int moduleName = HRR1;
			infor.appendCodeSection(moduleName);

			// let's see whether we do it in file split mode?
			int nHRRLHS = (int)hrr1.countLHSIntNumbers();
			bool inFileSplit = false;
			if (nHRRLHS>infor.nLHSForHRR1Split) inFileSplit = true;
			infor.setFileSplitMode(moduleName,inFileSplit);
			hrr1.arrayIndexTransformation(infor);

			// print it
			hrr1.hrrPrint(infor);

			// now let's update the VRR shell quartet list must be in array form
			// we will pick up those who can not do HRR work for the both sides
			if (infor.fileSplit(moduleName)) {
				infor.updateVRRSQListInArray(outputSQList); 
			}
		}
	}

	////////////////////////////
	//      VRR section       //
	////////////////////////////
	RR vrr(VRR,infor.getVRRMethod(),outputSQList,unsolvedList);
	vrr.generateRRSQList(NULL_POS);

	// add the VRR section to infor
	int moduleName = VRR;
	infor.appendCodeSection(moduleName);

	// now build the vrr infor
	VRRInfor vrrInfor(infor,vrr);

	// update infor
	infor.updateVRRInfor(vrrInfor);

	// generate the VRR head
	vrrInfor.printVRRHead(infor);

	// generate the VRR code section
	vrr.vrrPrint(infor,vrrInfor);

	// now finally let's do contraction
	vrrInfor.vrrContraction(infor);
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
	if (infor.fileSplit(VRR_CONT)) {
		appendFuncPrototype(VRR_CONT,CPP);
	}
	vector<int> modules(infor.getSectionInfor());
	for(int iFunc=modules.size()-1; iFunc>=0; iFunc--) {
		int module = modules[iFunc];
		appendFuncPrototype(module,CPP);
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
	if (infor.fileSplit(moduleName)) {
		formWorkFile(moduleName);
		formFunctionCall(moduleName,CPP);
	}else{
		appendFile(moduleName,CPP);
	}

	// contraction part
	// for contraction if it's not in file split
	// mode, it will be with VRR code
	moduleName = VRR_CONT;
	if (infor.fileSplit(moduleName)) {
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
		if (infor.fileSplit(module)) {

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
		for(int i=0; i<lp.getNPieces(); i++) {

			// we begin to determine that whether this is variable name
			string val = lp.findValue(i);
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
	if (moduleName == VRR) {
		nFiles = 1;
	}else{
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
		// for VRR, we do not have sub files
		// so reset the iFile to be -1
		int fileIndex = iFile;
		if (moduleName == VRR) fileIndex = -1;
		string fileName = infor.getWorkFuncName(onlyWithFuncName,moduleName,fileIndex,inFinalDir);

		// now let's open it
		ofstream CPP;
		CPP.open(fileName.c_str(),ios::out);

		// firstly create the include part
		infor.headPrinting(CPP);

		// now print the function prototyoe
		string pro;
		if (moduleName == VRR) {
			pro = prototype[0];
		}else{
			pro = prototype[iFile-1];
		}
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
	// here we have an exception for the VRR variable statement
	// in fact, we may not have this file generated
	// so if we do not have this file, we do not do a crash
	// just return
	path p(file.c_str());
	if (checkExist) {
		if (! exists(p)) {
			cout << "Missing file " << file << endl;
			crash(true, "Why the file does not exist??? we need it in appendFile function of sqints");
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

