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
#include "nonrr.h"
#include "boost/lexical_cast.hpp"
#include <boost/algorithm/string.hpp>   // string handling
#include "nonrrinfor.h"
using namespace printing;
using namespace inttype;
using namespace integral;
using namespace sqintsinfor;
using namespace nonrr;
using boost::lexical_cast;
using namespace boost;
using namespace nonrrinfor;

void NONRRInfor::declareArray(const SQIntsInfor& infor) const 
{
	//
	// If we do not have last section, we have nothing to declare
	// so the declare only happens when it's non-RR following by
	// derivatives calculation
	//
	if (nextSection == NULL_POS) {
		return;
	}
	int nSpace = 2;

	// now open the file
	string varFileName = infor.getWorkFuncName(false,section);
	ofstream varfile;
	varfile.open (varFileName.c_str(),std::ofstream::app);
	varfile << endl;

	string line;
	line = "/************************************************************";
	printLine(nSpace,line,varfile);
	line = " * declare the NONRR result shell quartets in array form";
	printLine(nSpace,line,varfile);
	line = " ************************************************************/";
	printLine(nSpace,line,varfile);

	// now print out all of result shell quartets
	// which are in array type
	string arrayType = getArrayType();
	for(int iSQ=0; iSQ<(int)outputSQList.size(); iSQ++) {
		const ShellQuartet& sq = outputSQList[iSQ];
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

	// now close the result statment part
	varfile << endl;
	varfile.close();
}

void HRRInfor::formSubFiles(const SQIntsInfor& infor, const NONRR& nonrr) 
{
	// set up a working copy of sub file record
	SubFileRecord record(section);
	record.init();
	int nLHS = 0;

	// set the limit value
	int nLHSLimit = nLHSForDerivSplit;
	if (section != DERIV) {
		nLHSLimit = nLHSForNonRRSplit;
	}

	//
	// here we will form the sub file according to the printing order of rrsq list
	// basically, it's reverse order of rrsqlist
	//
	// here we can not consider the function parameter number for each sub file
	// this is because the sub file input/output can be only formed when we have
	// all of sub files formed
	//
	const list<RRSQ>& rrsqList = nonrr.getRRSQList();
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
		if (nLHS>nLHSLimit) {
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

	// if this is the last section, we need to 
	// tell all of sub file records that this is 
	// the end
	if (isLastSection()) {
		for(int iSub=0; iSub<(int)subFilesList.size(); iSub++) {
			subFilesList[iSub].updateLHSSQStatus(GLOBAL_RESULT_SQ);
		}
	}

	// now let's print out the function statement
	int statementFile = DERIV_FUNC_STATEMENT;
	if (section == NON_RR) {
		statementFile = NON_RR_FUNC_STATEMENT;
	}
	for(int iSub=0; iSub<(int)subFilesList.size(); iSub++) {

		// get the function parameter
		const SubFileRecord& record = subFilesList[iSub];

		// reserve space for the output
		string arg;
		arg.reserve(1000);

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
		string line = "// NON-RR working function: " + 
			boost::lexical_cast<string>(fileIndex);
		printLine(0,line,file);
		printLine(0,func,file);
		file << endl;

		// now close file
		file.close();
	}
}

NONRRInfor::NONRRInfor(const SQIntsInfor& infor, const NONRR& nonrr):Infor(infor),
	nonrrFileSplit(false),section(nonrr.getSection()),nextSection(infor.nextSection(section)),
	inputSQList(nonrr.getBottomSQList()),outputSQList(nonrr.getResultSQList()),
	inputSQStatus(inputSQList.size(),VARIABLE_SQ),outputSQStatus(outputSQList.size(),VARIABLE_SQ),
	outputSQIntNumList(outputSQList.size(),0)
{
	// let's count how many LHS for non-RR
	int nLHS = nonrr.countLHSIntNumbers();

	// set the limit value
	int nLHSLimit = nLHSForDerivSplit;
	if (section != DERIV) {
		nLHSLimit = nLHSForNonRRSplit;
	}

	// now let's determine
	nonrrFileSplit = false;
	if (nLHS>nLHSLimit) {
		nonrrFileSplit = true;
	}

	// update the output shell quartet integral number
	const vector<set<int> >& intList = nonrr.getResultIntList();
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

void NONRRInfor::updateOutputSQInArray(const vector<ShellQuartet>& moduleInput)
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

void NONRRInfor::updateInputSQInArray(const vector<ShellQuartet>& moduleOutput)
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
