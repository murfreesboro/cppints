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
#include "sqintsinfor.h"
#include "rr.h"
#include "boost/lexical_cast.hpp"
#include <boost/algorithm/string.hpp>   // string handling
#include "hrrinfor.h"
using namespace printing;
using namespace inttype;
using namespace integral;
using namespace sqintsinfor;
using namespace rr;
using boost::lexical_cast;
using namespace boost;
using namespace hrrinfor;

void HRRInfor::declareArray(const SQIntsInfor& infor) const 
{
	//
	// we note that even if the given HRR is not in file split,
	// we may also need to declare array results. This is because
	// the array form result may be passed into following module
	// like DERIV
	//

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
	line = " * declare the HRR result shell quartets in array form";
	printLine(nSpace,line,varfile);
	line = " ************************************************************/";
	printLine(nSpace,line,varfile);

	// now print out all of shell quartets
	string arrayType = getArrayType();
	for(int iSQ=0; iSQ<(int)outputSQList.size(); iSQ++) {

		// test that whether the sq is not HRR output?
		// if not, this may be the previous HRR output,
		// or even HRR output
		// it's just passed from other modules like deriv
		const ShellQuartet& sq = outputSQList[iSQ];
		if (! sq.canDoHRR(side)) continue;

		// whether it's in array? we only declare array form
		int status = outputSQStatus[iSQ];
		if (! inArrayStatus(status)) continue;

		// now this is the declare part
		int nInts = outputSQIntNumList[iSQ];
		string arrayName = sq.formArrayName(section);
		string nLHSInts  = boost::lexical_cast<string>(nInts);
		string declare   = getArrayDeclare(nLHSInts);
		line = arrayType + arrayName + declare;
		printLine(nSpace,line,varfile);
	}

	//
	// finally, for HRR split into multiple files; we need to declare
	// the array form of input/output shell quartets for the HRR
	// sub functions
	//
	// we note, that this only applies to the case that has more than
	// one sub file record. 
	//
	if (subFilesList.size() > 1) {
		for(int iSubFile=0; iSubFile<(int)subFilesList.size()-1; iSubFile++) {
			const vector<ShellQuartet>& output = subFilesList[iSubFile].getLHSSQList();
			const vector<int>& status = subFilesList[iSubFile].getLHSSQStatus();
			for(int iSQ=0; iSQ<(int)output.size(); iSQ++) {

				// get the shell quartet
				const ShellQuartet& sq = output[iSQ];
				if(status[iSQ] != FUNC_INOUT_SQ) continue;

				// let's see whether it has been already declared
				vector<ShellQuartet>::const_iterator it = find(outputSQList.begin(),outputSQList.end(),sq);
				if (it != outputSQList.end()) continue;

				// now print
				int nInts        = subFilesList[iSubFile].getLHSSQIntNum(sq);
				string name      = sq.formArrayName(section);
				string arrayType = getArrayType();
				string declare   = getArrayDeclare(lexical_cast<string>(nInts));
				string line      = arrayType + name + declare;
				printLine(nSpace,line,varfile);
			}
		}
	}

	// now close the result statment part
	varfile << endl;
	varfile.close();
}

void HRRInfor::formSubFiles(const SQIntsInfor& infor, const RR& hrr) 
{
	// set up a working copy of sub file record
	SubFileRecord record(section);
	record.init();
	int nLHS = 0;

	//
	// here we will form the sub file according to the printing order of rrsq list
	// basically, it's reverse order of rrsqlist
	//
	// here we can not consider the function parameter number for each sub file
	// this is because the sub file input/output can be only formed when we have
	// all of sub files formed
	//
	const list<RRSQ>& rrsqList = hrr.getRRSQList();
	for(list<RRSQ>::const_reverse_iterator it=rrsqList.rbegin(); it!=rrsqList.rend(); ++it) {

		// add this RRSQ
		record.updateRRSQ(*it);

		// count the LHS 
		const list<int>& LHS = it->getLHSIndexArray();
		const ShellQuartet& sq = it->getLHSSQ();
		nLHS += LHS.size();

		// do we reach the sub file limit?
		// if so we add it to the record list,
		// and everything restarts
		if (nLHS>nLHSForHRRSplit) {
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
	// we need to know which one is the module input/output
	for(int iSub=0; iSub<(int)subFilesList.size(); iSub++) {
		SubFileRecord& subFile = subFilesList[iSub];
		subFile.updateModuleOutput(infor,outputSQList);
		subFile.updateModuleInput(inputSQList);
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
	int statementFile = HRR1_FUNC_STATEMENT;
	if (section == HRR2) {
		statementFile = HRR2_FUNC_STATEMENT;
	}
	for(int iSub=0; iSub<(int)subFilesList.size(); iSub++) {

		// get the function parameter
		const SubFileRecord& record = subFilesList[iSub];

		// reserve space for the output
		string arg;
		arg.reserve(1000);

		// add in AB or CD
		if (side == BRA1 || side == BRA2 || side == BRA) {
			arg = arg + "const Double* A, const Double* B, "
		}else{	
			arg = arg + "const Double* C, const Double* D, "
		}

		// now let's consider the input shell quartets
		// all of module input, as well as function input
		// are marked here
		// we note that the bottom sq must be from VRR,
		// therefore it's a module input
		const vector<ShellQuartet>& rhs = record.getRHSSQList();
		const vector<int>&    rhsStatus = record.getRHSSQStatus();
		for(int iSQ=0; iSQ<(int)rhs.size(); iSQ++) {
			const ShellQuartet& sq = rhs[iSQ];
			if (rhsStatus[iSQ] == BOOTTOM_SQ) {
				Integral I(sq,0);
				string line = "const Double* " + I.getName() + ", ";
				arg = arg + line;
			}else{
				if (rhsStatus[iSQ] != FUNC_INOUT_SQ) continue;
				string line = "const Double* " + sq.formArrayName(section) + ", ";
				arg = arg + line;
			}
		}

		// now it's the output shell quartets
		// firstly we need to consider the output as function output
		const vector<ShellQuartet>& lhs = record.getLHSSQList();
		const vector<int>&    lhsStatus = record.getLHSSQStatus();
		bool hasABCD = false;
		for(int iSQ=0; iSQ<(int)lhs.size(); iSQ++) {
			if (lhsStatus[iSQ] == GLOBAL_RESULT_SQ) {
				hasABCD = true;
				continue;
			}
			if (lhsStatus[iSQ] != FUNC_INOUT_SQ) continue;
			const ShellQuartet& sq = lhs[iSQ];
			string line = "Double* " + sq.formArrayName(section) + ", ";
			arg = arg + line;
		}

		// now let's see whether we have global result
		if (hasABCD) {
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
		int fileIndex = iSub + 1;
		string name = infor.getWorkFuncName(true,section,fileIndex);
		string func = "void " + name + "(" + arg + ");";

		// open file stream, write the result
		string f    = infor.getWorkFuncName(false,statementFile);
		ofstream file;
		file.open(f.c_str(),std::ofstream::app);

		// add some comments
		file << endl;
		string line = "// HRR working function: " + 
			boost::lexical_cast<string>(fileIndex);
		printLine(0,line,file);
		printLine(0,func,file);
		file << endl;

		// now close file
		file.close();
	}

	// finally update the status of input and output
	// if in file split mode, all of inout sq must be 
	// in array form
	for(int iSQ=0; iSQ<(int)outputSQStatus.size(); iSQ++) {
		if (outputSQStatus[iSQ] == VARIABLE_SQ) {
			outputSQStatus[iSQ] = ARRAY_SQ;
		}
	}
	for(int iSQ=0; iSQ<(int)inputSQStatus.size(); iSQ++) {
		if (inputSQStatus[iSQ] == VARIABLE_SQ) {
			inputSQStatus[iSQ] = ARRAY_SQ;
		}
	}
}

HRRInfor::HRRInfor(const SQIntsInfor& infor, const RR& hrr):Infor(infor),hrrFileSplit(false),
	section(hrr.getSection()),nextSection(infor.nextSection(section)),oper(infor.getOper()),
	side(hrr.getSide()),inputSQList(hrr.getHRRBottomSQList()),outputSQList(hrr.getRRResultSQList()),
	inputSQStatus(inputSQList.size(),VARIABLE_SQ),outputSQStatus(outputSQList.size(),VARIABLE_SQ),
	outputSQIntNumList(outputSQList.size(),0)
{
	// let's count how many LHS for HRR
	int nLHS = hrr.countLHSIntNumbers();

	// now let's determine
	hrrFileSplit = false;
	if (nLHS>nHRRFileSplit) {
		hrrFileSplit = true;
	}

	// update the output shell quartet integral number
	const vector<set<int> >& intList = hrr.getRRUnsolvedIntList();
	for(int iSQ=0; iSQ<(int)outputSQList.size(); iSQ++) {
		const set<int>& intNumSet = intList[iSQ];
		outputSQIntNumList[iSQ] = intNumSet.size();
	}

	// update the output sq list status
	for(int iSQ=0; iSQ<(int)outputSQList.size(); iSQ++) {
		const ShellQuartet& sq = outputSQList[iSQ];
		if (sq.isSTypeSQ()) {
			outputSQStatus[iSQ] = BOTTOM_SQ;
		}else if (infor.isResult(sq)) {
			outputSQStatus[iSQ] = GLOBAL_RESULT_SQ;
		}
	}

	// update the input sq list status
	// in case we have bottom integral
	for(int iSQ=0; iSQ<(int)inputSQList.size(); iSQ++) {
		const ShellQuartet& sq = inputSQList[iSQ];
		if (sq.isSTypeSQ()) {
			inputSQStatus[iSQ] = BOTTOM_SQ;
		}
	}
}

void HRRInfor::updateOutputSQInArray(const vector<ShellQuartet>& moduleInput)
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

void HRRInfor::updateInputSQInArray(const vector<ShellQuartet>& moduleOutput)
{
	for(int iSQ=0; iSQ<(int)inputSQList.size(); iSQ++) {
		const ShellQuartet& sq = inputSQList[iSQ];
		if (inputSQStatus[iSQ] != VARIABLE_SQ) continue;
		vector<ShellQuartet>::const_iterator it = find(moduleOutput.begin(),moduleOutput.end(),sq);
		if (it != moduleOutput.end()) {
			inputSQStatus[iSQ] = ARRAY_SQ;
		}
	}
}
