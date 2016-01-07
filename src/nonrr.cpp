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
#include "inttype.h"
#include "sqintsinfor.h"
#include "printing.h"
#include "rr.h"
#include "expfacinfor.h"
#include "nonrr.h"
#include "boost/lexical_cast.hpp"
using namespace inttype;
using namespace sqintsinfor;
using namespace printing;
using namespace rr;
using namespace expfacinfor;
using namespace nonrr;

NONRR::NONRR(const vector<ShellQuartet>& sqlist, 
		const vector<set<int> >& intList):codeSec(DERIV),derivOrder(sqlist[0].getDerivJobOrder()),
	oper(sqlist[0].getOper()),directions(1,NO_DERIV),
	inputSQList(sqlist),inputIntList(intList)
{
	// let's check that whether all of input shell quartets with 
	// derivatives information
	if (derivOrder > 0) {
		for(int iSQ=0; iSQ<(int)inputSQList.size(); iSQ++) {
			if (inputSQList[iSQ].getDerivJobOrder() != derivOrder) {
				cout << "the shell quartet is: " << inputSQList[iSQ].getName() << endl;
				crash(true, "one of shell quartet in IntDeriv constructor does not have deriv information");
			}	
		}
	}

	// now let's whether it's a valid non-RR non-Deriv work
	if (derivOrder == 0) {
		if (! isNONRROper(oper)) {
			crash(true, "the given operator can not do non-RR work in NONRR constructor");
		}
	}

	// reset the code section
	if (derivOrder == 0) {
		codeSec = NON_RR;
	}
	
	// initialize the direction information
	if (derivOrder == 0 && oper == THREEBODYKI) {
		directions.assign(3,NO_DERIV);
		directions[0] = DERIV_X;
		directions[1] = DERIV_Y;
		directions[2] = DERIV_Z;
	}
}

void NONRR::buildRRSQList()
{
	// the top loop is over the directions
	// for deriv job, the directions is with one NULL element
	for(int iDir=0; iDir<(int)directions.size(); iDir++) {
		int dir  = directions[iDir];

		//
		// now let's go to the build process
		// in this process, we convert the input 
		// sq and unsolved integral list into rrsq
		// push them into the rrsqlist
		//
		for(int iSQ=0; iSQ<(int)inputSQList.size(); iSQ++) {

			// for given sq, build rrsq
			const ShellQuartet& sq  = inputSQList[iSQ];
			const set<int>& intList = inputIntList[iSQ];

			// build the rrsq and insert the new one into the result list
			RRSQ rrsq(NULL_POS,NULL_POS,sq,intList,dir);
			rrsqList.push_back(rrsq);

			// now create new shell quartet list for the output 
			for(int item=0; item<rrsq.getNItems(); item++) {

				// get the corresponding sq
				const ShellQuartet& rhsSQ = rrsq.getRHSSQ(item);

				// obtain unsolved integral list
				set<int> newUnsolvedList;
				rrsq.getUnsolvedIntList(item,newUnsolvedList);

				// whether this sq is already contained in the tmp sq list?
				int pos = -1;
				for(int i=0; i<(int)resultSQList.size(); i++) {
					if (resultSQList[i] == rhsSQ) {
						pos = i;
						break;
					}
				}

				// if it's totally a new one, we just do simple update
				// else we need to merge the new unsolved list with the 
				// one already existing in the tmpUnsolvedIntList
				if (pos == -1) {
					resultSQList.push_back(rhsSQ);
					unsolvedIntList.push_back(newUnsolvedList);
				}else{

					// now we fetch the old unsolved list from tmpUnsolvedIntList
					set<int>& oldUnsolvedList = unsolvedIntList[pos];

					// we need to merge the new and old 
					for(set<int>::iterator it=newUnsolvedList.begin(); 
							it != newUnsolvedList.end(); ++it) {
						int val = *it;
						set<int>::iterator it2=oldUnsolvedList.find(val);
						if (it2 == oldUnsolvedList.end()) {
							oldUnsolvedList.insert(val);
						}
					}
				}
			}
		}
	}

	// finally, we need to sort the rrsq list
	rrsqList.sort();
}

int NONRR::countLHSIntNumbers() const
{
	int totalNInts = 0;
	for(list<RRSQ>::const_iterator it=rrsqList.begin(); it!=rrsqList.end(); ++it) {
		const list<int>& LHS = it->getLHSIndexArray();
		totalNInts += LHS.size();
	}
	return totalNInts;
}

size_t NONRR::evalDerivIntProcess(const SQIntsInfor& infor)
{
	// make sure this is a derivatives job
	if (derivOrder == 0) {
		crash(true, "it should be a derivative job in function of evalDerivIntProcess");
	}

	// this is the recorder of the result
	size_t nTotalLHSInts = 0;

	// firstly, let's build the rrsq list for the derivatives 
	buildRRSQList();
	nTotalLHSInts += (size_t)countLHSIntNumbers();

	// get the output shell quartets and it's integral list
	vector<ShellQuartet> outputSQList(resultSQList); 
	vector<set<int> > unsolvedList(unsolvedIntList); 

	// do we have additional non-RR work?
	if (isNONRROper(oper)) {

		// build the rrsq
		NONRR nonRRJob(outputSQList,unsolvedList);
		nonRRJob.buildRRSQList();
		nTotalLHSInts += (size_t)nonRRJob.countLHSIntNumbers();

		// now rewrite the shell quartet list
		outputSQList.clear();
		unsolvedList.clear();
		outputSQList = nonRRJob.getResultSQList();
		unsolvedList = nonRRJob.getResultIntSQList();
	}

	// now let's do RR work, consider whether it has HRR
	// this could be determined by 
	if (infor.hasHRR()) {

		// determine the first and second side
		RR hrr(HRR2,HRR,outputSQList,unsolvedList);
		int firstSide  = NULL_POS;
		int secondSide = NULL_POS;
		hrr.sideDeterminationInHRR(firstSide,secondSide);

		// now do the second side HRR
		// for generated code, second side is performed after
		// first side; so in generation of code second side 
		// is prior to the first side
		if (secondSide != NULL_POS) {
			hrr.generateRRSQList(secondSide);
			nTotalLHSInts += (size_t)hrr.countLHSIntNumbers();

			// now rewrite the shell quartet list
			outputSQList.clear();
			unsolvedList.clear();
			outputSQList = hrr.getHRRBottomSQList();
			unsolvedList = hrr.getHRRBottomIntList();
		}

		// do you have work on the first side?
		if (firstSide != NULL_POS) {

			// now do the first side
			RR hrr1(HRR1,HRR,outputSQList,unsolvedList);
			hrr1.generateRRSQList(firstSide);
			nTotalLHSInts += (size_t)hrr1.countLHSIntNumbers();

			// now rewrite the shell quartet list
			outputSQList.clear();
			unsolvedList.clear();
			outputSQList = hrr1.getHRRBottomSQList();
			unsolvedList = hrr1.getHRRBottomIntList();
		}
	}

	// finally let's do the VRR work
	// this is necessary for all of works
	bool doVRRWork = true;
	if (doVRRWork) {
		RR vrr(VRR,infor.getVRRMethod(),outputSQList,unsolvedList);
		vrr.generateRRSQList(NULL_POS);
		nTotalLHSInts += (size_t)vrr.countLHSIntNumbers();
	}

	// now basically everything is down
	return nTotalLHSInts;
}

void NONRR::nonRRPrint(const SQIntsInfor& infor) const
{
	// this is only for nonrr part
	if (codeSec != NON_RR) {
		crash(true, "fatal error in NONRR::nonRRPrint, only non_rr can call nonRRPrint");
	}

	// determine that how many space should be given for each line printing
	int nSpace = 2;
	if (resultIntegralHasAdditionalOffset(oper) && ! infor.fileSplit(codeSec)) {
		nSpace += 2;
	}

	// do we have declaration?
	// that's only possible when this is not the last section
	if (! infor.isLastSection(codeSec)) {

		// now get some vector to store the results
		vector<ShellQuartet> lhsSQList;
		lhsSQList.reserve(200);
		vector<int> lhsNumList;
		lhsNumList.reserve(200);

		// now grasp all of LHS shell quartets
		// we say that all of LHS must be defined, they are not null value
		for(list<RRSQ>::const_reverse_iterator it=rrsqList.rbegin(); it!=rrsqList.rend(); ++it) {
			const ShellQuartet& sq = it->getLHSSQ();
			if (infor.isResult(sq)) continue;
			lhsSQList.push_back(sq);
			const list<int>& LHS = it->getLHSIndexArray();
			int nLHS = LHS.size();
			lhsNumList.push_back(nLHS);
		}

		// let's print the heading part
		string varFileName = infor.getWorkFuncName(false,codeSec);
		ofstream varfile;
		varfile.open (varFileName.c_str(),std::ofstream::app);
		varfile << endl;
		string line;
		line = "/************************************************************";
		printLine(nSpace,line,varfile);
		line = " * declare the NON_RR result shell quartets";
		printLine(nSpace,line,varfile);
		line = " ************************************************************/";
		printLine(nSpace,line,varfile);

		// now print out all of shell quartets
		string arrayType = infor.getArrayType();
		for(int iSQ=0; iSQ<(int)lhsSQList.size(); iSQ++) {
			const ShellQuartet& sq = lhsSQList[iSQ];
			int nInts = lhsNumList[iSQ];
			string arrayName = sq.getName();
			string nLHSInts  = boost::lexical_cast<string>(nInts);
			string declare   = infor.getArrayDeclare(nLHSInts);
			line = arrayType + arrayName + declare;
			printLine(nSpace,line,varfile);
		}

		// now close the file
		varfile << endl;
		varfile.close();
	}

	// let's prepare something in constructing the code files
	// firstly, if the HRR part of codes is placed in a lot of functions
	// we need the function parameters
	vector<ShellQuartet> inputList;
	inputList.reserve(100);
	vector<ShellQuartet> outputList;
	outputList.reserve(100);

	// prepare the file index 
	// if we do file split, it starts from 1
	int index = -1;
	if (infor.fileSplit(codeSec)) index = 1;

	// shall we do restart for a new file?
	bool restart = false;

	// set the function parameter numbers
	int nFuncPar = 0;

	// set the lines of code
	int nLHS = 0;

	// do we have global results?
	bool hasABCD = false;

	// we need to know that totally how much rrsq we have
	int nRRSQ = rrsqList.size();
	int rrsqIndex = 1;

	// set the exponential information
	list<RRSQ>::const_reverse_iterator it=rrsqList.rbegin();
	const ShellQuartet& lhsSQ = it->getLHSSQ();
	ExpFacInfor expInfor(lhsSQ.getExpFacList(),lhsSQ.getExpFacListLen());

	// loop over the rr sq list
	for(it=rrsqList.rbegin(); it!=rrsqList.rend(); ++it) {

		// open the file
		string fileName = infor.getWorkFuncName(false,codeSec,index);
		ofstream myfile;
		myfile.open (fileName.c_str(),std::ofstream::app);
		it->nonRRPrint(codeSec,nSpace,infor,myfile);
		myfile.close();

		// if we are in file split mode
		if (infor.fileSplit(codeSec)) {

			// update input and output
			const ShellQuartet& sq = it->getLHSSQ();
			if (infor.isResult(sq)) {
				hasABCD = true;
			}else{
				outputList.push_back(sq);
				nFuncPar += 1;
			}
			for(int item=0; item<it->getNItems(); item++) {
				const ShellQuartet& rhsSQ = it->getRHSSQ(item);
				vector<ShellQuartet>::const_iterator it1 = std::find(inputList.begin(),inputList.end(),rhsSQ);
				if (it1 != inputList.end()) {
					continue;
				}
				vector<ShellQuartet>::const_iterator it2 = std::find(outputList.begin(),outputList.end(),rhsSQ);
				if (it2 != outputList.end()) {
					continue;
				}
				inputList.push_back(rhsSQ);
				nFuncPar += 1;
			}

			// update the nLHS
			const list<int>& LHS = it->getLHSIndexArray();
			nLHS += LHS.size();

			// shall we stop the file?
			if (nFuncPar>infor.maxParaForFunction) {
				restart = true;
			}else if (nLHS>infor.nLHSForNonRRSplit) {
				restart = true;
			} 
			
			// check exp infor
			if (infor.withExpFac()) {
				ExpFacInfor newExpInfor(sq.getExpFacList(),sq.getExpFacListLen());
				if (newExpInfor != expInfor) {
					restart = true;
				}
				expInfor = newExpInfor;
			}

			// check whether we reach the end
			if (rrsqIndex == nRRSQ) {
				restart = true;
			}

			// now let's deal with the situation we need to restart
			if (restart) {

				// firstly let's form the function name
				string funcName = infor.getWorkFuncName(true,codeSec,index);

				// form parameters
				string arg;
				for(int iSQ=0; iSQ<(int)inputList.size(); iSQ++) {
					const ShellQuartet& sq = inputList[iSQ];
					arg  = arg + "const Double* " + sq.getName() + ", ";
				}

				// now it's output
				for(int iSQ=0; iSQ<(int)outputList.size(); iSQ++) {
					const ShellQuartet& sq = outputList[iSQ];
					if (iSQ == (int)outputList.size()-1) {
						arg  = arg + "Double* " + sq.getName();
					}else{
						arg  = arg + "Double* " + sq.getName() + ", ";
					}
				}

				// do we add in abcd?
				if (hasABCD) {
					if (outputList.size()>0) {
						arg = arg + ", ";
					}
					arg = arg + "Double* abcd";
				}

				// now create the function prototype 
				string line = "void " + funcName + "( " + arg + " );";
				string prototype = infor.getWorkFuncName(false,NON_RR_FUNC_STATEMENT);
				ofstream pro;
				pro.open(prototype.c_str(),std::fstream::app);
				printLine(0,line,pro);
				pro << endl;
				pro.close();

				// clear input and output list
				inputList.clear();
				outputList.clear();

				// clear the counting
				nLHS = 0;
				nFuncPar = 0;

				// increase the index
				index += 1;

				// reset the hasABCD
				hasABCD = false;
			}

			// now finally reset the restart status
			restart = false;
		}

		// increase the rrsq index
		rrsqIndex += 1;
	}
}

void NONRR::derivPrint(const SQIntsInfor& infor) const
{
	// this is right now only applies for deriv section
	if (codeSec != DERIV) {
		crash(true, "print function for derivPrint in NONRR only for derivatives calculation");
	}

	// determine that how many space should be given for each line printing
	// for non-RR it's always has indent of 2
	int nSpace = 2;

	// additionally, for non-RR if it's inside additional loop like ESP etc.
	// we need to consider add more nSpace
	if (resultIntegralHasAdditionalOffset(oper) && ! infor.fileSplit(codeSec)) {
		nSpace += 2;
	}

	// let's prepare something in constructing the code files
	// firstly, if the HRR part of codes is placed in a lot of functions
	// we need the function parameters
	vector<ShellQuartet> inputList;
	inputList.reserve(100);

	// prepare the file index 
	// if we do file split, it starts from 1
	int index = -1;
	if (infor.fileSplit(codeSec)) index = 1;

	// shall we do restart for a new file?
	bool restart = false;

	// set the function parameter numbers
	// because the the output will always be "abcd"
	// so initialize the nFuncPar to be 1
	int nFuncPar = 1;

	// set the lines of code
	int nLHS = 0;

	// we need to know that totally how much rrsq we have
	int nRRSQ = rrsqList.size();
	int rrsqIndex = 1;

	// loop over the rr sq list
	for(list<RRSQ>::const_reverse_iterator it=rrsqList.rbegin(); it!=rrsqList.rend(); ++it) {

		// print out the codes
		string fileName = infor.getWorkFuncName(false,codeSec,index);
		ofstream myfile;
		myfile.open (fileName.c_str(),std::ofstream::app);
		it->nonRRPrint(codeSec,nSpace,infor,myfile);
		myfile.close();

		// if we are in file split mode
		if (infor.fileSplit(codeSec)) {

			// count the RHS
			for(int item=0; item<it->getNItems(); item++) {
				const ShellQuartet& rhsSQ = it->getRHSSQ(item);
				vector<ShellQuartet>::const_iterator it1 = std::find(inputList.begin(),inputList.end(),rhsSQ);
				if (it1 != inputList.end()) {
					continue;
				}
				inputList.push_back(rhsSQ);
				nFuncPar += 1;
			}

			// update the nLHS
			const list<int>& LHS = it->getLHSIndexArray();
			nLHS += LHS.size();

			// shall we stop the file?
			if (nFuncPar>infor.maxParaForFunction) {
				restart = true;
			}else if (nLHS>infor.nLHSForDerivSplit) {
				restart = true;
			}

			// now we reach the end of file
			if (rrsqIndex == nRRSQ) {
				restart = true;
			}

			// now let's deal with the situation we need to restart
			if (restart) {

				// firstly let's form the function name
				string funcName = infor.getWorkFuncName(true,codeSec,index);

				// form parameters
				string arg;

				// now add input shell quartets
				for(int iSQ=0; iSQ<(int)inputList.size(); iSQ++) {
					const ShellQuartet& sq = inputList[iSQ];
					arg  = arg + "const Double* " + sq.getName() + ", ";
				}
				arg = arg + "Double* abcd";

				// now create the function prototype 
				string line = "void " + funcName + "( " + arg + " );";
				string prototype = infor.getWorkFuncName(false,DERIV_FUNC_STATEMENT);
				ofstream pro;
				pro.open(prototype.c_str(),std::fstream::app);
				printLine(0,line,pro);
				pro << endl;
				pro.close();

				// clear input and output list
				inputList.clear();

				// clear the counting
				nLHS = 0;
				nFuncPar = 1;

				// increase the index
				index += 1;
			}

			// now finally reset the restart status
			restart = false;
		}

		// increase the rrsq index
		rrsqIndex += 1;
	}
}

void NONRR::arrayIndexTransformation(const SQIntsInfor& infor)
{
	// transform the integral index into the array index
	// this is performed to the local results
	// where the RHS is defined inside the module
	if (infor.withArrayIndex(codeSec)) {
		for(list<RRSQ>::iterator it=rrsqList.begin(); it!=rrsqList.end(); ++it) {
			for(list<RRSQ>::iterator it2=rrsqList.begin(); it2!=rrsqList.end(); ++it2) {
				it2->rhsArrayIndexTransform(*it);
			}
		}
	}

	// still the RHS
	// now let's see the RHS which are bottom integrals
	// these bottom integrals are defined in the previous section 
	if (infor.withArrayIndex(codeSec)) {
		for(list<RRSQ>::iterator it=rrsqList.begin(); it!=rrsqList.end(); ++it) {
			it->rhsArrayIndexTransform(resultSQList,unsolvedIntList);
		}
	}

	// now let's do the LHS
	// if this is the last section, the LHS will be just the 
	// final results, we do not need to do anything
	// if this is not the last section, we need to see what is 
	// the status for the next section
	// we note, both NON_RR and DERIV code section; they do not
	// have TMP_RESULTS
	// all of LHS results are module results or final results
	if (! infor.isLastSection(codeSec)) {

		// next section information
		int nextCode = infor.nextSection(codeSec);
		if (nextCode == NULL_POS) {
			crash(true, "this is meaningless that the next code section is null in NONRR::arrayIndexTransformation");
		}

		// whether it's with array
		if (infor.withArrayIndex(nextCode)) {
			for(list<RRSQ>::iterator it=rrsqList.begin(); it!=rrsqList.end(); ++it) {
				it->lhsArrayIndexTransform();
			}
		}
	}
}

