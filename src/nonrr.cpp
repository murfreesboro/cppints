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
#include "nonrrinfor.h"
#include "nonrr.h"
#include "boost/lexical_cast.hpp"
using namespace inttype;
using namespace sqintsinfor;
using namespace printing;
using namespace nonrrinfor;
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

int NONRR::countRHSIntNumbers() const 
{
	int totalNInts = 0;
	for(list<RRSQ>::const_iterator it=rrsqList.begin(); it!=rrsqList.end(); ++it) {
		totalNInts += it->countRHSIntegralNum();
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
	size_t nTotalRHSInts = 0;

	// firstly, let's build the rrsq list for the derivatives 
	buildRRSQList();
	nTotalRHSInts += (size_t)countRHSIntNumbers();

	// get the output shell quartets and it's integral list
	vector<ShellQuartet> outputSQList(resultSQList); 
	vector<set<int> > unsolvedList(unsolvedIntList); 

	// do we have additional non-RR work?
	if (isNONRROper(oper)) {

		// build the rrsq
		NONRR nonRRJob(outputSQList,unsolvedList);
		nonRRJob.buildRRSQList();
		nTotalRHSInts += (size_t)nonRRJob.countRHSIntNumbers();

		// now rewrite the shell quartet list
		outputSQList.clear();
		unsolvedList.clear();
		outputSQList = nonRRJob.getBottomSQList();
		unsolvedList = nonRRJob.getBottomIntSQList();
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
			nTotalRHSInts += (size_t)hrr.countRHSIntNumbers();

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
			nTotalRHSInts += (size_t)hrr1.countRHSIntNumbers();

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

		// for VRR we have simulated contraction coefficients
		// infor class knows about it
		int contractionDegree = infor.getVRRContDegree();
		size_t nVRRRHSInts = vrr.countRHSIntNumbers()*contractionDegree;
		nTotalRHSInts += nVRRRHSInts;
	}

	// now basically everything is down
	return nTotalRHSInts;
}

void NONRR::print(const SQIntsInfor& infor, const NONRRInfor& nonrrinfor) 
{
	// now let's conver the LHS part
	const vector<ShellQuartet>& outputList = nonrrinfor.getOutputSQList();
	const vector<int>& outputStatusList    = nonrrinfor.getOutputSQStatus();

	// we note, that all of LHS are module output
	// because we do not have temp results, so LHS never
	// to be RHS in the module, we just change the LHS part
	for(int iSQ=0; iSQ<(int)outputList.size(); iSQ++) {
		int status = outputStatusList[iSQ];
		if (! inArrayStatus(status)) continue;
		const ShellQuartet& sq = outputList[iSQ];
		for(list<RRSQ>::iterator it=rrsqList.begin(); it!=rrsqList.end(); ++it) {
			const ShellQuartet& lhsSQ = it->getLHSSQ();
			if (sq == lhsSQ) {
				it->updateLHSSQStatus(status);
				it->lhsArrayIndexTransform();
			}
		}
	}

	// let's see whether we convert the rhs into array form?
	bool handleRHSArrayForm = false;
	if (nonrrinfor.fileSplit()) {

		// get the total integral number
		// basically all of RHS are refering unsolved integral list
		// so we just count the unsolved integral list 
		int nInts = 0;
		for(int iSQ=0; iSQ<(int)unsolvedIntList.size(); iSQ++) {
			const set<int>& intList = unsolvedIntList[iSQ];
			nInts += intList.size();
		}

		// get the nLHS limit
		int nLHSLimit = nonrrinfor.nLHSForDerivSplit;
		if (codeSec == NON_RR) {
			nLHSLimit = nonrrinfor.nLHSForNonRRSplit;
		}

		// let's see whether we just convert them into var?
		Double nonRRCoef = nonrrinfor.getNonRRCoefs();
		if (nInts>nLHSLimit*nonRRCoef) {
			handleRHSArrayForm = true;
		}
	}else{

		// for the case of non-file split case, we should do
		// handling RHS array. Because the input for non-RR
		// may be in array form
		handleRHSArrayForm = true;
	}

	// now let's see whether we transform the RHS part
	if (handleRHSArrayForm) {
		const vector<int>& resultSQStatusList = infor.getInputSQStatus();
		for(list<RRSQ>::iterator it=rrsqList.begin(); it!=rrsqList.end(); ++it) {
			it->rhsArrayIndexTransform(resultSQList,resultSQStatusList,unsolvedIntList);
		}
	}

	// now let's print out the code
	if (nonrrinfor.fileSplit()) {

		// set the space 
		int nSpace = 2;

		// loop over sub files
		for(int iSub=0; iSub<nonrrinfor.getNSubFiles(); iSub++) {
			const SubFileRecord& record = nonrrinfor.getSubFileRecord(iSub);

			// now set up the file
			int fileIndex = iSub + 1;
			string filename = infor.getWorkFuncName(false,codeSec,fileIndex);
			ofstream myfile;
			myfile.open (filename.c_str(),std::ofstream::app);

			// let's convert the RHS from array to var
			if (! handleRHSArrayForm) {
				const vector<ShellQuartet>& rhs = record.getRHSSQList();
				const vector<int>&    rhsStatus = record.getRHSSQStatus();
				for(int iSQ=0; iSQ<(int)rhs.size(); iSQ++) {
					if (rhsStatus[iSQ] != FUNC_INOUT_SQ) continue;
					const ShellQuartet& sq = rhs[iSQ];
					for(int iSQ2=0; iSQ2<(int)resultSQList.size(); iSQ2++) {
						if (resultSQList[iSQ2] == sq) {
							const set<int>& intList = unsolvedIntList[iSQ2];

							// print out title
							string line;
							line = "/************************************************************";
							printLine(nSpace,line,myfile);
							line = " * convert from array into variable form: " + sq.getName();
							printLine(nSpace,line,myfile);
							line = " ************************************************************/";
							printLine(nSpace,line,myfile);

							// now do the convert
							int offset = 0;
							string arrayName = oriSQ.formArrayName(codeSec);
							for(set<int>::const_iterator it=intList.begin(); it!=intList.end(); ++it) {
								string rhs = arrayName + "[" + lexical_cast<string>(offset) + "]";
								int index = *it;
								Integral I(sq,index);
								string varName = I.formVarName(codeSec);
								string expression = "Double " + varName + " = " + rhs + ";";
								printLine(nSpace,expression,myfile);
								offset++;
							}
						}
					}
				}
			}

			// now let's print out the stuff here
			const vector<ShellQuartet>& lhs = record.getLHSSQList();
			for(list<RRSQ>::const_reverse_iterator it=rrsqList.rbegin(); it!=rrsqList.rend(); ++it) {
				const ShellQuartet& lhsSQ = it->getLHSSQ();
				vector<ShellQuartet>::const_iterator it = find(lhs.begin(),lhs.end(),lhsSQ);
				if (it != lhs.end()) {
					it->print(nSpace,infor,myfile);
				}
			}

			// now close the file
			myfile.close();
		}
	}else{

		// determine the space
		int oper   = workSQList[0].getOper();
		int nSpace = 2;
		if (resultIntegralHasAdditionalOffset(oper)) {
			nSpace += 2;
		}

		// create the RR file
		// this is already in a mod of a+
		string filename = infor.getWorkFuncName(false,codeSec);
		ofstream myfile;
		myfile.open (filename.c_str(),std::ofstream::app);

		// now do the printing work
		// we will print each rrsq in reverse order
		for(list<RRSQ>::const_reverse_iterator it=rrsqList.rbegin(); it!=rrsqList.rend(); ++it) {
			it->print(nSpace,infor,myfile);
		}

		// now close the whole file
		myfile.close();
	}
}

