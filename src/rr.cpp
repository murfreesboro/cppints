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
#include "rrsqsearch.h"
#include "sqintsinfor.h"
#include "rrints.h"
#include "printing.h"
#include "boost/lexical_cast.hpp"
#include "integral.h"
#include "inttype.h"
#include "intderiv.h"
#include "derivinfor.h"
#include "rr.h"
using namespace rrsqsearch;
using namespace sqintsinfor;
using namespace rrints;
using namespace printing;
using boost::lexical_cast;
using namespace integral;
using namespace inttype;
using namespace intderiv;
using namespace derivinfor;
using namespace rr;

///////////////////////////////////////////////////////////////////////////
//       ####     constructors related functions                         //
///////////////////////////////////////////////////////////////////////////
RR::RR(const int& rrType0, const vector<ShellQuartet>& inputSQList0,
	const SQIntsInfor& infor, int side0):rrType(rrType0),
	jobOrder(infor.getJobOrder()),side(side0),inputSQList(inputSQList0)
{
	// firstly, let's go to create the initsqlist for RR work
	vector<ShellQuartet> initSQList(inputSQList);
	int oper = inputSQList[0].getOper();
	if (needIntDeriv(oper)) {

		// build the int deriv 
		IntDeriv intDeriv(rrType,jobOrder,inputSQList);

		// set up the initial sq list for position search
		initSQList.clear();
		vector<set<int> > unsolvedIntList;
		intDeriv.getUnsolvedSQIntList(initSQList,unsolvedIntList);
	}

	// build the optimum RR path in advance 
	// and collecting the result path
	if (rrType == HRR) {
		RRSQSearch sqOptSearch(side,initSQList);
		optRRList = sqOptSearch.getSolvedSQList();
		posList   = sqOptSearch.getPosList(); 
	}else{
		RRSQSearch sqOptSearch(initSQList,rrType);
		optRRList = sqOptSearch.getSolvedSQList();
		posList   = sqOptSearch.getPosList(); 
	}

	// now ready to build the rrsqlist
	formRRSQList();
	bool doArrayIndexTransform = infor.withArrayIndex(rrType);
	if (doArrayIndexTransform) arrayIndexTransformation();
}

int RR::searchPos(const ShellQuartet& sq) const
{
	// for VRR, if it's the only non-S position
	// then we will just use it
	// no search at all
	if (rrType != HRR && sq.canDODirectRRPosSearch()) {
		return sq.getRRPos(); 
	}

	// search the matching shell quartet
	int pos = -1;
	for(int iSQ=0; iSQ<(int)optRRList.size(); iSQ++) {
		if (sq == optRRList[iSQ]) {
			pos = iSQ;
			break;
		}
	}

	// now get the position
	if (pos == -1) {
		cout << "shell quartet name " << sq.getName() << endl;
		crash(true,"We did not get position information in RR for this shell quartet");
		return -1;
	}else{
		return posList[pos];
	}
}

void RR::lhsIntegralCheck(const ShellQuartet& sq, set<int>& intList) const
{
	bool findSQ = false;
	set<int> workIntList(intList);
	for(list<RRSQ>::const_iterator it=rrsqList.begin(); it!=rrsqList.end(); ++it) {
		const ShellQuartet& lhsSQ = it->getLHSSQ();
		if (lhsSQ == sq) {
			findSQ = true;
			const list<int>& lhs = it->getLHSIndexArray();
			for(set<int>::const_iterator it2 = workIntList.begin(); 
					it2 != workIntList.end(); ++it2) {
				int val = *it2;
				list<int>::const_iterator it3 = find(lhs.begin(),lhs.end(),val);
				if (it3 != lhs.end()) intList.erase(val);
			}
		}
	}

	// finally, debug check
	if (! findSQ) {
		cout << "sq's name: " << sq.getName() << endl;
		crash(true,"missing the sq in RR::lhsIntegralCheck");
	}
}

bool RR::buildRRSQList(vector<ShellQuartet>& sqlist, vector<set<int> >& unsolvedIntList)
{
	//
	// first step, check the existing rrsqlist
	// whether the shell quartets in the sqlist
	// already been exisiting in the results?
	// in this case, we just need to update 
	// the result with the new unsolved integral 
	// list
	// the newSQIndex is used to hold the sq index
	// which is going into build process (-1 is old,
	// 1 is new)
	//
	bool finishBuildingProcess = true;
	vector<int> newSQIndex(sqlist.size(),-1);
	for(int iSQ=0; iSQ<(int)sqlist.size(); iSQ++) {

		// get the new sq going to be updated
		const ShellQuartet& sq = sqlist[iSQ];

		// we test that whether this sq is bottom type?
		// bottom sq should be the input for the given RR module
		// for example, VRR -> S integrals (directly calculated)
		// HRR -> bottom shell quartets (from VRR part)
		// if this is bottom sq, then we do nothing
		if (rrType == HRR) {
			if (! sq.canDoHRR(side)) continue;
		}else{
			if (sq.isSTypeSQ()) continue;
		}

		// now test whether we have it in the rrsqlist
		// if we have, do updating work
		// we note that sq in rrsqlist are different 
		// from each other, therefore we only have one
		// matched with the given sq
		bool matchOldSQ = false;
		const set<int>& intList = unsolvedIntList[iSQ];
		for(list<RRSQ>::iterator it=rrsqList.begin(); it!=rrsqList.end(); ++it) {
			const ShellQuartet& lhsSQ = it->getLHSSQ();
			if (lhsSQ == sq) {
				it->updateLHS(intList);
				matchOldSQ = true;
				break;
			}
		}

		// now update the index
		// and we know at least we have one sq
		// need to go into building process
		if (! matchOldSQ) {
			newSQIndex[iSQ] = 1;
			finishBuildingProcess = false;
		}
	}

	// before we go to building process, let's 
	// check whether we can finish the whole 
	// building process?
	if (finishBuildingProcess) return true;

	//
	// now let's go to the build process
	// in this process, we convert the input 
	// sq and unsolved integral list into rrsq
	// push them into the rrsqlist
	// then going to generate the new input 
	// for the next round of building process
	//

	// set up the vector to hold the input for 
	// next around of building process
	vector<ShellQuartet> tmpSQList;
	vector<set<int> > tmpUnsolvedIntList;

	// working for each sq - to make it into rrsq
	for(int iSQ=0; iSQ<(int)sqlist.size(); iSQ++) {

		// whether this sq is new one?
		if (newSQIndex[iSQ] < 0) {
			continue;
		}

		// for given sq, build rrsq
		const ShellQuartet& sq = sqlist[iSQ];
		const set<int>& intList = unsolvedIntList[iSQ];
		int pos = searchPos(sq);
		RRSQ rrsq(rrType,pos,sq,intList);

		// debug code
#ifdef RR_DEBUG
		cout << "now do expansion on position " << pos << " for sq " << sq.getName() << endl;
#endif

		// insert the new one into the result list
		// we do not need to check it anymore
		// this must be a new rrsq, and it's not 
		// the bottom type one
		rrsqList.push_back(rrsq);

		// now create new shell quartet list
		// and corresponding integral list so to do the next
		// round of building
		for(int item=0; item<rrsq.getNItems(); item++) {

			// get the corresponding sq
			const ShellQuartet& rhsSQ = rrsq.getRHSSQ(item);

			// obtain unsolved integral list
			set<int> newUnsolvedList;
			rrsq.getUnsolvedIntList(item,newUnsolvedList);

			// whether this sq is already contained in the tmp sq list?
			int pos = -1;
			for(int i=0; i<(int)tmpSQList.size(); i++) {
				if (tmpSQList[i] == rhsSQ) {
					pos = i;
					break;
				}
			}

			// if it's totally a new one, we just do simple update
			// else we need to merge the new unsolved list with the 
			// one already existing in the tmpUnsolvedIntList
			if (pos == -1) {
				tmpSQList.push_back(rhsSQ);
				tmpUnsolvedIntList.push_back(newUnsolvedList);
			}else{

				// now we fetch the old unsolved list from tmpUnsolvedIntList
				set<int>& oldUnsolvedList = tmpUnsolvedIntList[pos];

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

	// now based on th tmp LHS, we need to create the new ones
	sqlist.clear();
	unsolvedIntList.clear();
	for(int iSQ=0; iSQ<(int)tmpSQList.size(); iSQ++) {
		sqlist.push_back(tmpSQList[iSQ]);
		unsolvedIntList.push_back(tmpUnsolvedIntList[iSQ]);
	}

	// finally, we have to return false
	return false;
}

void RR::rrUpdating(vector<ShellQuartet>& sqlist, vector<set<int> >& unsolvedIntList)
{
	//
	// We note, that basically this function's structure 
	// is very similar to the function above 
	// buildRRSQList
	//

	while(true) {

		// do updating the the given integral list
		for(int iSQ=0; iSQ<(int)sqlist.size(); iSQ++) {
			const ShellQuartet& sq  = sqlist[iSQ];
			const set<int>& intList = unsolvedIntList[iSQ];
			for(list<RRSQ>::iterator it=rrsqList.begin(); it!=rrsqList.end(); ++it) {
				const ShellQuartet& lhsSQ = it->getLHSSQ();
				if (lhsSQ == sq) {
					it->updateLHS(intList);
				}
			}
		}

		// set up the vector to hold the input for 
		// next around of update process
		vector<ShellQuartet> tmpSQList;
		vector<set<int> > tmpUnsolvedIntList;

		// now it's time to prepare the integral for the 
		// next round
		for(int iSQ=0; iSQ<(int)sqlist.size(); iSQ++) {
			const ShellQuartet& sq  = sqlist[iSQ];
			const set<int>& intList = unsolvedIntList[iSQ];
			for(list<RRSQ>::const_iterator it=rrsqList.begin(); it!=rrsqList.end(); ++it) {
				const ShellQuartet& lhsSQ = it->getLHSSQ();
				if (lhsSQ == sq) {
					//it->debug_print();

					// now let's go to see the RHS part
					for(int item=0; item<it->getNItems(); item++) {

						// get the corresponding sq
						const ShellQuartet& rhsSQ = it->getRHSSQ(item);

						// is it a bottom SQ?
						if (rrType == HRR) {
							if (! rhsSQ.canDoHRR(side)) continue;
						}else{
							if (rhsSQ.isSTypeSQ()) continue;
						}

						// obtain unsolved integral list
						set<int> newIntList;
						it->getUnsolvedIntList(item,intList,newIntList);

						// now we are ready to push the raw results
						// get the position of rhsSQ in the result list
						int pos = -1;
						for(int i=0; i<(int)tmpSQList.size(); i++) {
							if (tmpSQList[i] == rhsSQ) {
								pos = i;
								break;
							}
						}

						// if it's totally a new one, we just add it
						// else we will insert the new unsolved 
						// integrals into the existing list by merging 
						// them together
						if (pos == -1) {
							tmpSQList.push_back(rhsSQ);
							tmpUnsolvedIntList.push_back(newIntList);
						}else{
							set<int>& oldUnsolvedList = tmpUnsolvedIntList[pos];
							for(set<int>::iterator it2=newIntList.begin(); 
									it2 != newIntList.end(); ++it2) {
								int val = *it2;
								set<int>::iterator it3=oldUnsolvedList.find(val);
								if (it3 == oldUnsolvedList.end()) {
									oldUnsolvedList.insert(val);
								}
							}
						}
					}
				}
			}
		}

		//
		// now check the raw results
		// it's mostly possible that the newly get
		// rhs integrals is already exsiting as LHS
		// we will abadon these ones
		// 
		for(int iSQ=0; iSQ<(int)tmpSQList.size(); iSQ++) {
			lhsIntegralCheck(tmpSQList[iSQ],tmpUnsolvedIntList[iSQ]);
		}

		// now we need to conver the tmp result into the real one
		// that is to say, convert RHS into the new RHS
		// so step into another loop
		sqlist.clear();
		unsolvedIntList.clear();
		for(int iSQ=0; iSQ<(int)tmpSQList.size(); iSQ++) {
			const ShellQuartet& sq  = tmpSQList[iSQ];
			const set<int>& intList = tmpUnsolvedIntList[iSQ];
			if (intList.size() > 0) {
				sqlist.push_back(sq);
				unsolvedIntList.push_back(intList);
			}
		}

		// see if it's the right time to step out?
		if (sqlist.size() == 0) break;
	}
}

void RR::completenessCheck() 
{
	// set up results to hold the missing integrals
	// in the completeness check
	vector<ShellQuartet> sqlist;
	vector<set<int> > unsolvedIntList;

	// now it's time to check the completeness
	// the first run
	bool inComplete = true;
	for(list<RRSQ>::const_iterator it=rrsqList.begin(); it!=rrsqList.end(); ++it) {

		// get the shell quartet as well as LHS index
		set<int> lhsIndex;
		const ShellQuartet& lhsSQ = it->getLHSSQ();
		it->formLHSIndexSet(lhsIndex);

		// now compare with the rest of rrsq
		for(list<RRSQ>::const_iterator it2=rrsqList.begin(); it2!=rrsqList.end(); ++it2) {
			set<int> intList;
			bool isComplete = it2->checkCompleteness(lhsSQ,lhsIndex,intList);
			if (! isComplete) {

				// update status
				inComplete = false;

				// now we are ready to push the results
				int pos = -1;
				for(int i=0; i<(int)sqlist.size(); i++) {
					if (sqlist[i] == lhsSQ) {
						pos = i;
						break;
					}
				}

				// if it's totally a new one, we just add it
				// else we will insert the new unsolved 
				// integrals into the existing list 
				if (pos == -1) {
					sqlist.push_back(lhsSQ);
					unsolvedIntList.push_back(intList);
				}else{
					set<int>& oldUnsolvedList = unsolvedIntList[pos];
					for(set<int>::iterator it3=intList.begin(); it3 != intList.end(); ++it3) {
						int val = *it3;
						oldUnsolvedList.insert(val);
					}
				}
			}
		}
	}

	// now let's see whether we do updating work?
	if (! inComplete) {
		rrUpdating(sqlist,unsolvedIntList);
	}

	// finally, for the incomplete case we need to 
	// check the completeness for second run
	// this time the rrsq should be complete, else
	// we just report a crash
	if (! inComplete) {
		for(list<RRSQ>::const_iterator it=rrsqList.begin(); it!=rrsqList.end(); ++it) {

			// get the shell quartet as well as LHS index
			set<int> lhsIndex;
			const ShellQuartet& lhsSQ = it->getLHSSQ();
			it->formLHSIndexSet(lhsIndex);

			// now compare with the rest of rrsq
			for(list<RRSQ>::const_iterator it2=rrsqList.begin(); it2!=rrsqList.end(); ++it2) {
				set<int> intList;
				bool isComplete = it2->checkCompleteness(lhsSQ,lhsIndex,intList);
				if (! isComplete) {
					crash(true, "Dead in the completeness check function!!!");
				}
			}
		}
	}
}

void RR::formRRSQList()
{
	// first step, build the rrsq list
	// we intentionally put everythig into a braket
	// so that to make the building process local
	bool doBuildingWork = true;
	if (doBuildingWork) {

		// set up the initial sq list and it's unsolved int list
		vector<ShellQuartet> sqlist(inputSQList);
		vector<set<int> > unsolvedIntList;

		// now form the unsolved integral list
		int oper = inputSQList[0].getOper();
		if (needIntDeriv(oper)) {
			IntDeriv intDeriv(rrType,jobOrder,inputSQList);
			sqlist.clear();
			intDeriv.getUnsolvedSQIntList(sqlist,unsolvedIntList);
			intDeriv.updateRRSQList(rrsqList);
		}else{
			int length = sqlist.size();
			unsolvedIntList.reserve(length);
			for(int iSQ=0; iSQ<length; iSQ++) {
				const ShellQuartet& sq = sqlist[iSQ];
				set<int> intList;
				sq.getIntegralList(intList);
				unsolvedIntList.push_back(intList);
			}
		}

		// now do the real building work
		// this is stepping into RR cycles
		while(true) {
			bool workFinish = buildRRSQList(sqlist,unsolvedIntList);
			if (workFinish) break;
		}
	}

	// now we should have all of stuff inside 
	// do sorting work
	// see the operator < for more information
	// about how we sort the rrsq
	rrsqList.sort();
#ifdef RR_DEBUG
	cout << "print out the result RRSQ list " << endl;
	for(list<RRSQ>::const_iterator it=rrsqList.begin(); it!=rrsqList.end(); ++it) {
		const ShellQuartet& lhsSQ = it->getLHSSQ();
		cout << "RRSQ LHS name " << lhsSQ.getName() << endl;
	}
#endif

	// now it's time to check the completeness
	completenessCheck(); 
}

void RR::arrayIndexTransformation()
{
	// transform the integral index into the array index
	// first step, transform the RHS into array index from the given LHS
	for(list<RRSQ>::iterator it=rrsqList.begin(); it!=rrsqList.end(); ++it) {
		for(list<RRSQ>::iterator it2=rrsqList.begin(); it2!=rrsqList.end(); ++it2) {
			it2->rhsArrayIndexTransform(*it);
		}
	}

	// now let's do the LHS
	for(list<RRSQ>::iterator it=rrsqList.begin(); it!=rrsqList.end(); ++it) {
		it->lhsArrayIndexTransform();
	}
}

///////////////////////////////////////////////////////////////////////////
//                  ####     utility functions                           //
///////////////////////////////////////////////////////////////////////////
int RR::sqStatusCheck(const SQIntsInfor& infor, const ShellQuartet& sq) const 
{
	//
	// for more information about the meaning of tmp results, module results etc.
	// please refer to constant integers defined in the head of the rrints.cpp
	//

	// by default, it's tmp result
	int status = TMP_RESULT;

	// is it the module result sq?
	bool isModuleResult = false;
	vector<ShellQuartet>::const_iterator it2 = find(inputSQList.begin(),inputSQList.end(),sq);
	if (it2 != inputSQList.end()) {
		status = MODULE_RESULT;
		isModuleResult = true;
	}

	// is it the final result?
	// since all of VRR result are done in contraction part
	// therefore all of VRR part of results can not be 
	// final results
	if (isModuleResult && rrType == HRR) {
		bool isFinalResult = infor.isResult(sq);
		if (isFinalResult) {
			status = FINAL_RESULT;
		}
	}

	// now return the result sq status
	return status;
}

/*
void RR::printHead(const int& nSpace, const int& status, const SQIntsInfor& infor, 
		const RRSQ& rrsq, ofstream& file) const 
{
	// obtain the information for this rrsq
	const ShellQuartet& sq = rrsq.getLHSSQ();
	const list<int>& lhs = rrsq.getLHSIndexArray();
	string name = sq.getName();
	int nInts = lhs.size();
	int nTotalInts = sq.getNInts();
	int diff = nTotalInts - nInts;

	// now print comment section to file
	file << endl;
	string line;
	line = "************************************************************";
	printLine(nSpace,line,file);

	// print out shell quartet infor
	// here if this is result sq, we will make a special mark
	if (status != TMP_RESULT) {
		line = " * result shell quartet name: " + name;
		printLine(nSpace,line,file);
	}else{
		line = " * shell quartet name: " + name;
		printLine(nSpace,line,file);
	}
	line = " * totally " + lexical_cast<string>(diff) + " integrals are omitted ";
	printLine(nSpace,line,file);

	// now finalize the comment printing
	line = " ************************************************************";
	printLine(nSpace,line,file);

	// check that wether we need to do additional declare for this section
	// of sq printing?
	if (infor.withArrayIndex(rrType)) {

		// for HRR/VRR, we print head for both tmp results and module results
		// this is just in normal way
		// for VRR, we also print it for tmp results and module results
		// however, for module results, things is a bit of different
		// that is, for non-composite shell quartets, we need to add _vrr modifier
		// so that to distingurish the final result one
		// also, if the position is derivatives, we also add it as modifier, too
		//
		// we note, that this should be in conststence with the print function
		// in rrsq
		if (status != FINAL_RESULT) {

			// get the vector name
			string vecName = "vector<Double> " + name;

			// do we need vrr modifier?
			if (rrType != HRR && status == MODULE_RESULT && ! infor.isComSQ()) {
				vecName = vecName + "_vrr";
			}

			// do we need to add in deriv modifier?
			int pos = rrsq.getPosition();
			if (isDerivInfor(pos)) {
				string derivInfor = symTransform(pos);
				vecName = vecName + "_" + derivInfor;
			}

			// now do array declare
			line = vecName + "(" + lexical_cast<string>(nInts) + ",0.0E0);";
			printLine(nSpace,line,file);
		}
	}
}
*/

void RR::print(const SQIntsInfor& infor, const string& filename) 
{
	// determine that how many space should be given for each line printing
	// for HRR it's always has indent of 2
	// for VRR we need to consider the contraction loop
	// however, if it's in split file choice, then nspace should be 2
	int nSpace = 2;
	int oper   = inputSQList[0].getOper();
	if (rrType != HRR && ! infor.splitCPPFile()) {
		nSpace     = getNSpaceByOper(oper);
	}

	// additionally, for HRR if it's inside additional loop like ESP etc.
	// we need to consider add more nSpace
	if (resultIntegralHasAdditionalOffset(oper) && rrType == HRR) {
		nSpace += 2;
	}

	// create the RR file
	// this is already in a mod of a+
	ofstream myfile;
	myfile.open (filename.c_str(),std::ofstream::app);

	// now do the printing work
	// we will print each rrsq in reverse order
	for(list<RRSQ>::reverse_iterator it=rrsqList.rbegin(); it!=rrsqList.rend(); ++it) {

#ifdef RR_DEBUG
		const ShellQuartet& lhssq = it->getLHSSQ();
		cout << "In printing function, current sq: " << lhssq.getName() << endl;
#endif

		// firstly, checking the status of the result sq
		const ShellQuartet& sq = it->getLHSSQ();
		int status = sqStatusCheck(infor,sq);

		// print code for this section of rrsq
		it->print(nSpace,status,infor,inputSQList,myfile);
	}

	// now close the whole file
	myfile.close();
}

void RR::getBottomSQList(vector<ShellQuartet>& sqlist) const
{
	//
	// for VRR part, bottom integrals is clear; it's just
	// the SSSS type of integrals which will be directly
	// calculated etc.
	// therefore, here the bottom sq list is only for 
	// HRR part.
	// We do not search rrsqlist, since it does not 
	// contain the bottom shell quartets. However,
	// all of shell quartet in rrsqlist must be same 
	// with the optRRList. Therefore we do searching
	// in the optRRList
	//
#ifdef RR_DEBUG
	if (rrType != HRR) {
		crash(true, "getBottomSQList is useless for VRR process!");
	}
#endif

	// now the work
	for(int iSQ=0; iSQ<(int)optRRList.size(); iSQ++) {
		const ShellQuartet& sq = optRRList[iSQ];
		if (! sq.canDoHRR(side)) {

			// additionally, we need to check whether this sq is 
			// already in the result list
			bool hasThisSQ = false;
			for(int i=0; i<(int)sqlist.size(); i++) {
				if (sq == sqlist[i]) hasThisSQ = true;
				break;
			}

			// now push into result
			if (! hasThisSQ) {
				sqlist.push_back(sq);
			}
		}
	}
}

/*
 
void RR::buildRRSQList()
{
	// firstly, set up the vector to hold the working 
	// shell quartets and unsolved integral list
	vector<ShellQuartet> sqlist;
	vector<set<int> > unsolvedIntList;

	// initilize the work
	sqlist = initSQList;
	unsolvedIntList.reserve(sqlist.size());
	for(int iSQ=0; iSQ<(int)sqlist.size(); iSQ++) {
		set<int> intList;
		sqlist[iSQ].getIntegralList(intList);
		unsolvedIntList.push_back(intList);
	}

	// now step into the loop until all of sq are bottom ones
	while(true) {

#ifdef DEBUG
		cout << "result RRSQ list currently " << endl;
		for(list<RRSQ>::const_iterator it=rrsqList.begin(); it!=rrsqList.end(); ++it) {
			const ShellQuartet& lhsSQ = it->getLHSSQ();
			cout << lhsSQ.getName() << endl;
		}
#endif

		// set up the tmp list to store this round of RHS result
		vector<ShellQuartet> tmpSQList;
		vector<set<int> > tmpUnsolvedIntList;

		// working for each sq - to make it into rrsq
		for(int iSQ=0; iSQ<(int)sqlist.size(); iSQ++) {

			// for each given sq list, build rrsq
			const ShellQuartet& sq = sqlist[iSQ];
			const set<int>& intList = unsolvedIntList[iSQ];
			int pos = searchPos(sq);
			RRSQ rrsq(rrType,pos,sq,intList);

#ifdef DEBUG
			const ShellQuartet& lhssq = rrsq.getLHSSQ();
			cout << "now do expansion on position " << pos << " for sq " 
				<< lhssq.getName() << endl;
#endif

			// insert the new one into the result list
			rrsqList.push_back(rrsq);

			// now for this rrsq, we have to do two things
			// firstly, for each RHS term of rrsq, we need to update
			// the existing LHS in rrsq list
			// second, we need to create new shell quartet list
			// and corresponding integral list so to do the next
			// round of creation
			for(int item=0; item<rrsq.getNItems(); item++) {

				// whether this is bottom sq?
				// bottom sq should be the input for the given RR
				// module
				// for example, VRR -> S integrals (directly calculated)
				// HRR -> bottom shell quartets (from VRR part)
				const ShellQuartet& rhsSQ = rrsq.getRHSSQ(item);
				if (rrType == HRR) {
					if (! rhsSQ.canDoHRR(side)) continue;
				}else{
					if (rhsSQ.isSTypeSQ()) continue;
				}

				// obtain unsolved integral list
				set<int> unsolvedList;
				rrsq.getUnsolvedIntList(item,unsolvedList);

				// updating work
				updateRRSQList(rhsSQ,unsolvedList);

				// whether this sq is already contained in next round of lhs sq list?
				if (tmpSQList.size() > 0) {
					vector<ShellQuartet>::iterator it = find(tmpSQList.begin(),tmpSQList.end(),rhsSQ);
					if (it != tmpSQList.end()) continue;
				}

				// now create tmp LHS for next round
				tmpSQList.push_back(rhsSQ);
				tmpUnsolvedIntList.push_back(unsolvedList);
			}
		}

		// now based on th tmp LHS, we need to create the 
		// real LHS
		sqlist.clear();
		unsolvedIntList.clear();
		for(int iSQ=0; iSQ<(int)tmpSQList.size(); iSQ++) {
			if (isNewRRSQ(tmpSQList[iSQ])) {
				sqlist.push_back(tmpSQList[iSQ]);
				unsolvedIntList.push_back(tmpUnsolvedIntList[iSQ]);
			}
		}

		// shall we break out?
		if (sqlist.size() == 0) break;
	}
}



void RR::buildRRSQList(const ShellQuartet& sq, const set<int>& intList)
{
	// build the rrsq
	int pos = searchPos(sq);
	RRSQ rrsq(rrType,pos,sq,intList);

#ifdef DEBUG
	cout << "now form rrsq on position " << pos << " for sq " << sq.getName() << endl;
#endif

	// update the rrsq
	updateRRSQList(rrsq);

	// now for this rrsq, we will recursively call the function
	// to build other rrsq
	for(int item=0; item<rrsq.getNItems(); item++) {

		// whether this is bottom sq?
		// bottom sq should be the input for the given RR
		// module
		// for example, VRR -> S integrals (directly calculated)
		// HRR -> bottom shell quartets (from VRR part)
		const ShellQuartet& rhsSQ = rrsq.getRHSSQ(item);
		if (rrType == HRR) {
			if (! rhsSQ.canDoHRR(side)) continue;
		}else{
			if (rhsSQ.isSTypeSQ()) continue;
		}

		// obtain unsolved integral list
		set<int> unsolvedList;
		rrsq.getUnsolvedIntList(item,unsolvedList);

		// finally, call the function itself
		buildRRSQList(rhsSQ,unsolvedList);
	}
}

*/
