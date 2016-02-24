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
#include "rrsqsearch.h"
#include "sqintsinfor.h"
#include "rrints.h"
#include "printing.h"
#include "boost/lexical_cast.hpp"
#include "integral.h"
#include "inttype.h"
#include "vrrinfor.h"
#include "hrrinfor.h"
#include "derivinfor.h"
#include "expfacinfor.h"
#include "rr.h"
using boost::lexical_cast;
using namespace rrsqsearch;
using namespace sqintsinfor;
using namespace rrints;
using namespace printing;
using namespace integral;
using namespace vrrinfor;
using namespace hrrinfor;
using namespace inttype;
using namespace derivinfor;
using namespace expfacinfor;
using namespace rr;

///////////////////////////////////////////////////////////////////////////
//       ####     RR forming related functions                           //
///////////////////////////////////////////////////////////////////////////
RR::RR(const int& codeSec0, const int& rrType0, const vector<ShellQuartet>& inputSQList0,
		const vector<set<int> >& inputUnsolvedIntList):codeSec(codeSec0),rrType(rrType0),side(NULL_POS),
	inputSQList(inputSQList0),workSQList(inputSQList0),initUnsolvedIntList(inputUnsolvedIntList)
{
	// firstly, make sure that this is only for RR work
	if (codeSec != HRR1 && codeSec != HRR2 && codeSec != VRR) {
		crash(true, "the input code section name in RR constructor is invalid");
	}
	if (! isValidVRRJob(rrType) && ! isValidHRRJob(rrType)) {
		crash(true, "the input rrType in RR constructor is invalid");
	}

	// also we check the length of workSQList comparing with initUnsolvedIntList
	if (workSQList.size() != initUnsolvedIntList.size()) {
		crash(true, "the input sq list in RR constructor has different length with inputUnsolvedIntList");
	}

	// remove exp factors etc. information
	// for VRR for workSQList
	if (codeSec == VRR) {
		for(int iSQ=0; iSQ<(int)workSQList.size(); iSQ++) {
			workSQList[iSQ].destroyMultipliers();
		}

		// save copies before work
		vector<ShellQuartet> sqlistCopy(workSQList);
		workSQList.clear();

		// now check the duplications
		// we need to do it for both shell quartets and 
		// corresponding integral list
		for(int iSQ=0; iSQ<(int)sqlistCopy.size(); iSQ++) {
			const ShellQuartet& sq = sqlistCopy[iSQ];
			vector<ShellQuartet>::const_iterator it = find(workSQList.begin(),workSQList.end(),sq);
			if (it == workSQList.end()) {
				workSQList.push_back(sq);
			}
		}
	}
}

void RR::updateSQIntListForVRR(vector<ShellQuartet>& sqlist, vector<set<int> >& unsolvedIntList) const
{
	// remove exp factors etc. information
	for(int iSQ=0; iSQ<(int)sqlist.size(); iSQ++) {
		sqlist[iSQ].destroyMultipliers();
	}

	// save copies before work
	vector<ShellQuartet> sqlistCopy(sqlist);
	vector<set<int> > unsolvedIntListCopy(unsolvedIntList);
	sqlist.clear();
	unsolvedIntList.clear();

	// now check the duplications
	// we need to do it for both shell quartets and 
	// corresponding integral list
	for(int iSQ=0; iSQ<(int)sqlistCopy.size(); iSQ++) {

		// find that whether we already have this sq 
		int pos = -1;
		const ShellQuartet& sq = sqlistCopy[iSQ];
		for(int i=0; i<(int)sqlist.size(); i++) {
			if (sqlist[i] == sq) {
				pos = i;
				break;
			}
		}

		// if this is totally new, we just push them back
		if (pos == -1) {
			sqlist.push_back(sq);
			const set<int>& intList = unsolvedIntListCopy[iSQ];
			unsolvedIntList.push_back(intList);
		}else{

			// for shell quartet, we do not need to do anymore
			// however if the shell quartet already exists, then we need to 
			// merge the duplicate integral list with the existing one

			// now let's get the new list and old list
			// old list: the integral list in the copy
			// new list: the integral list already contained 
			const set<int>& oldList = unsolvedIntListCopy[iSQ];
			set<int>& newList = unsolvedIntList[pos];
			for(set<int>::const_iterator it2 = oldList.begin(); it2 != oldList.end(); ++it2) {
				int val = *it2;
				set<int>::const_iterator it3 = newList.find(val);
				if (it3 == newList.end()) newList.insert(val);
			}
		}
	}
}

void RR::updateBottomSQIntsList(const ShellQuartet& sq, const set<int>& intList) 
{
	// let's check that whether the given sq is already in the bottom sq list
	int pos = -1;
	for(int iSQ=0; iSQ<(int)bottomSQList.size(); iSQ++) {
		if (bottomSQList[iSQ] == sq) {
			pos = iSQ;
			break;
		}
	}

	// now let's see whether we just push in the fresh new data?
	if (pos == -1) {
		bottomSQList.push_back(sq);
		bottomIntList.push_back(intList);
		return;
	}

	// now we need a merge work
	set<int>& oldList = bottomIntList[pos];
	for(set<int>::const_iterator it = intList.begin(); it != intList.end(); ++it) {
		int val = *it;
		set<int>::const_iterator it2 = oldList.find(val);
		if (it2 == oldList.end()) oldList.insert(val);
	}
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
		cout << "module name " << codeSec << endl;
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
	// the newSQIndex is used to hold the sq status
	// -1 is old, 1 is new
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
			if (! sq.canDoHRR(side)) {

				// for HRR we also need to update the bottom sq list
				const set<int>& intList = unsolvedIntList[iSQ];
				updateBottomSQIntsList(sq,intList); 
				continue;
			}
		}else{
			if (sq.isSTypeSQInRR()) continue;
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
		// this is all related to the RHS shell quartets
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

						// obtain unsolved integral list from the intList
						// the newIntlist contains correponding RHS integral 
						// index in forming the integrals in the intList
						set<int> newIntList;
						it->getUnsolvedIntList(item,intList,newIntList);

						// is it a bottom SQ?
						// for VRR case the bottom integral list
						// is in complete list
						// however, for HRR it may not
						// we update the bottom shell quartet 
						// as well as integral list here
						if (rrType != HRR) {
							if (rhsSQ.isSTypeSQInRR()) continue;
						}else{
							if (! rhsSQ.canDoHRR(side) && newIntList.size()>0) {
								updateBottomSQIntsList(rhsSQ,newIntList); 
								continue;
							}
						}

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
		// that is to say, convert RHS into the new LHS
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
		// basically, sqlist is the work SQ list, and 
		// we will make a copy of input unsolved int list here
		vector<ShellQuartet> sqlist(inputSQList);
		vector<set<int> > unsolvedIntList(initUnsolvedIntList);

		// remove all of modifier information for the work sq list
		// for VRR job, also we may need to revise unsolved int list
		if (isValidVRRJob(rrType)) {
			updateSQIntListForVRR(sqlist,unsolvedIntList);
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

void RR::generateRRSQList(const int& side0)
{
	// overwrite the side information
	// only used for HRR
	side = side0;
	if (rrType == HRR) {
		if (side != BRA && side != KET) {
			crash(true, "in the function of generateRRSQList in RR the input side information is not correct");
		}
	}

	// make a copy of work SQ List, we do not want to change it
	vector<ShellQuartet> initSQList(workSQList);

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
}

///////////////////////////////////////////////////////////////////////////
//                  ####     RR print related functions                  //
///////////////////////////////////////////////////////////////////////////
void RR::updateHRRInfor(const HRRInfor& infor)
{
	// check the rr type
	if (rrType != HRR) {
		crash(true, "RR::updateHRRInfor only applies for HRR");
	}

	// let's see whether we do file split
	if (infor.fileSplit()) {

		// for file split case, all of module output, and function output sq
		// has been recorded in the sub file list
		// so we just go over the sub file list to see the information
		for(int iSub=0; iSub<infor.getNSubFiles(); iSub++) {
			const SubFileRecord& record = infor.getSubFileRecord(iSub);
			const vector<ShellQuartet>& lhs = record.getLHSSQList();
			const vector<int>&    lhsStatus = record.getLHSSQStatus();
			for(int iSQ=0; iSQ<(int)lhs.size(); iSQ++) {
				int status = lhsStatus[iSQ];
				if (! inArrayStatus(status)) continue;
				const ShellQuartet& sq = lhs[iSQ];

				// now let's update all of LHS, and corresponding RHS
				for(list<RRSQ>::iterator it=rrsqList.begin(); it!=rrsqList.end(); ++it) {
					const ShellQuartet& lhsSQ = it->getLHSSQ();
					if (lhsSQ == sq) {
						it->updateLHSSQStatus(status);

						// it's possible that this LHS is used as RHS later in the same code section
						for(list<RRSQ>::iterator it2=rrsqList.begin(); it2!=rrsqList.end(); ++it2) {
							it2->rhsArrayIndexTransform(*it);
						}

						// also transform the LHS array index
						it->lhsArrayIndexTransform();
					}
				}
			}
		}

	}else{

		// we will use the HRR output sq list from hrr infor
		// however, it should be same with the input list 
		const vector<ShellQuartet>& outputList = infor.getOutputSQList();
		const vector<int>& outputStatusList    = infor.getOutputSQStatus();

		// in this case we need to update rrsq will the module output
		// if the module output sq in array form, we will update the status
		for(int iSQ=0; iSQ<(int)outputList.size(); iSQ++) {
			int status = outputStatusList[iSQ];
			const ShellQuartet& sq = outputList[iSQ];
			for(list<RRSQ>::iterator it=rrsqList.begin(); it!=rrsqList.end(); ++it) {
				const ShellQuartet& lhsSQ = it->getLHSSQ();
				if (sq == lhsSQ) {
					it->updateLHSSQStatus(status);

					// it's possible that this LHS is used as RHS later in the same code section
					// if this is not in array status we just continue
					if (! inArrayStatus(status)) continue;
					for(list<RRSQ>::iterator it2=rrsqList.begin(); it2!=rrsqList.end(); ++it2) {
						it2->rhsArrayIndexTransform(*it);
					}

					// also transform the LHS array index
					it->lhsArrayIndexTransform();
				}
			}
		}
	}

	// also update the rrsq with bottom sq list 
	const vector<int>& bottomSQStatusList = infor.getInputSQStatus();
	for(list<RRSQ>::iterator it=rrsqList.begin(); it!=rrsqList.end(); ++it) {
		it->rhsArrayIndexTransform(bottomSQList,bottomSQStatusList,bottomIntList);
	}
}

void RR::vrrPrint(const SQIntsInfor& infor, const VRRInfor& vrrinfor) const
{
	// print out VRR results
	vrrinfor.printVRRHead(infor);

	// the printing will be divided into two cases
	// one is with file split, and the other is not
	if (vrrinfor.fileSplit()) {

		// set the space 
		int nSpace = 2;

		// loop over sub files
		for(int iSub=0; iSub<vrrinfor.getNSubFiles(); iSub++) {
			const SubFileRecord& record = vrrinfor.getSubFileRecord(iSub);

			// now set up the file
			int fileIndex = iSub  + 1;
			string filename = infor.getWorkFuncName(false,VRR,fileIndex);
			ofstream myfile;
			myfile.open (filename.c_str(),std::ofstream::app);

			// we need to convert the input array shell quartet into variables
			if (iSub>0) {
				const vector<ShellQuartet>& rhs = record.getRHSSQList();
				const vector<int>&    rhsStatus = record.getRHSSQStatus();
				for(int iSQ=0; iSQ<(int)rhs.size(); iSQ++) {
					if (rhsStatus[iSQ] != FUNC_INOUT_SQ) continue;
					const ShellQuartet& sq = rhs[iSQ];
					for(list<RRSQ>::const_reverse_iterator it=rrsqList.rbegin(); it!=rrsqList.rend(); ++it) {
						const ShellQuartet& lhsSQ = it->getLHSSQ();
						if (sq == lhsSQ) {
							it->printArrayToVar(nSpace,myfile);
						}
					}
				}
			}

			// now let's print out the stuff here
			const vector<ShellQuartet>& lhs = record.getLHSSQList();
			for(list<RRSQ>::const_reverse_iterator it=rrsqList.rbegin(); it!=rrsqList.rend(); ++it) {
				const ShellQuartet& lhsSQ = it->getLHSSQ();
				vector<ShellQuartet>::const_iterator it2 = find(lhs.begin(),lhs.end(),lhsSQ);
				if (it2 != lhs.end()) {
					it->print(nSpace,infor,myfile);
				}
			}

			// finally we need to convert the output into array form, too
			const vector<int>&    lhsStatus = record.getLHSSQStatus();
			for(int iSQ=0; iSQ<(int)lhs.size(); iSQ++) {
				if (lhsStatus[iSQ] != FUNC_INOUT_SQ) continue;
				const ShellQuartet& sq = lhs[iSQ];
				for(list<RRSQ>::const_reverse_iterator it=rrsqList.rbegin(); it!=rrsqList.rend(); ++it) {
					const ShellQuartet& lhsSQ = it->getLHSSQ();
					if (sq == lhsSQ) {
						it->printVarToArray(nSpace,myfile);
					}
				}
			}

			// now close the file
			myfile.close();
		}
	}else{
	
		// determine the space
		int oper   = workSQList[0].getOper();
		int nSpace  = getNSpaceByOper(oper);
		if (resultIntegralHasAdditionalOffset(oper)) {
			nSpace += 2;
		}

		// create the RR file
		// this is already in a mod of a+
		string filename = infor.getWorkFuncName(false,VRR);
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

	// finally let's do contraction
	vrrinfor.vrrContraction(infor);
}

void RR::hrrPrint(const SQIntsInfor& infor, const HRRInfor& hrrinfor) 
{
	// this is only for HRR part
	if (codeSec != HRR1 && codeSec != HRR2) {
		crash(true, "fatal error in RR::hrrPrint, only HRR1/HRR2 can call hrrPrint");
	}

	// also check whether the section from hrr infor is same with this one
	if (codeSec != hrrinfor.getSection()) {
		crash(true, "fatal error in RR::hrrPrint, the section from hrrinfor is not same with the section in rr class");
	}

	// convert the rrsq in array index
	updateHRRInfor(hrrinfor);

	// let's do array declare here
	hrrinfor.declareArray(infor);

	// now let's print out the code
	if (hrrinfor.fileSplit()) {

		// set the space 
		int nSpace = 2;

		// loop over sub files
		for(int iSub=0; iSub<hrrinfor.getNSubFiles(); iSub++) {
			const SubFileRecord& record = hrrinfor.getSubFileRecord(iSub);

			// now set up the file
			int fileIndex = iSub + 1;
			string filename = infor.getWorkFuncName(false,codeSec,fileIndex);
			ofstream myfile;
			myfile.open (filename.c_str(),std::ofstream::app);

			// now print out the head
			string line;
			line = "/************************************************************";
			printLine(nSpace,line,myfile);
			line = " * initilize the HRR steps : build the AB/CD variables";
			printLine(nSpace,line,myfile);
			line = " ************************************************************/";
			printLine(nSpace,line,myfile);

			// is it AB or CD side?
			if (side == BRA) {
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

			// now let's print out the stuff here
			const vector<ShellQuartet>& lhs = record.getLHSSQList();
			for(list<RRSQ>::const_reverse_iterator it=rrsqList.rbegin(); it!=rrsqList.rend(); ++it) {
				const ShellQuartet& lhsSQ = it->getLHSSQ();
				vector<ShellQuartet>::const_iterator it2 = find(lhs.begin(),lhs.end(),lhsSQ);
				if (it2 != lhs.end()) {
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

		// now print out the head
		string line;
		line = "/************************************************************";
		printLine(nSpace,line,myfile);
		line = " * initilize the HRR steps : build the AB/CD variables";
		printLine(nSpace,line,myfile);
		line = " ************************************************************/";
		printLine(nSpace,line,myfile);

		// is it AB or CD side?
		if (side == BRA) {
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

		// now do the printing work
		// we will print each rrsq in reverse order
		for(list<RRSQ>::const_reverse_iterator it=rrsqList.rbegin(); it!=rrsqList.rend(); ++it) {
			it->print(nSpace,infor,myfile);
		}

		// now close the whole file
		myfile.close();
	}
}

void RR::sideDeterminationInHRR(int& firstSide, int& secondSide) const
{
	// firstly, check whether this is HRR
	if (rrType != HRR) {
		crash(true, "sideDeterminationInHRR is useless for non-HRR process!");
	}

	// we use the operator from the input sq list
	// they should have same operator
	int oper = inputSQList[0].getOper();

	// for some operators, we can not do HRR
	// now let's test it
	bool canHRR = canDOHRR(oper);
	if (! canHRR) {
		firstSide  = NULL_POS;
		secondSide = NULL_POS;
		return;
	}

	// determine the side from the operator
	int nBody = getOperOrder(oper);
	if (nBody == 1) {

		// no need to do HRR for one body integrals
		firstSide  = NULL_POS;
		secondSide = NULL_POS;
		return;

	}else if (nBody == 2 || nBody == 3) {

		// for two body and three body integrals,
		// ket side HRR is not necessary. Therefore
		// we only need to set the first side
		// and the first side is always BRA
		secondSide = NULL_POS;
		firstSide  = BRA;

		// now let's check whether for the 
		// first side we can do HRR
		bool canDOHRR = false;
		for(int iSQ=0; iSQ<(int)inputSQList.size(); iSQ++) {
			const ShellQuartet& sq = inputSQList[iSQ];
			if (sq.canDoHRR(firstSide)) {
				canDOHRR = true;
				break;
			}
		}
		if (! canDOHRR) {
			firstSide  = NULL_POS;
		}
		return;

	}else{

		//
		// for four body integrals, we need to see whether
		// ket side is first or the bra side is first?
		// 

		// firstly let's count the number of integrals
		// for each side 
		int nIntsBraSide = 0;
		int nIntsKetSide = 0;
		for(int iSide=0; iSide<2; iSide++) {
			int side = BRA;
			if (iSide == 1) side = KET;
			int num = 0;
			for(int iSQ=0; iSQ<(int)inputSQList.size(); iSQ++) {
				const ShellQuartet& sq = inputSQList[iSQ];
				num += sq.getNInts(side);
			}
			if (iSide == 0) {
				nIntsBraSide += num;
			}else{
				nIntsKetSide += num;
			}
		}

		// secondly let's whether HRR is not available 
		// for the side
		bool braCanDOHRR = false;
		bool ketCanDOHRR = false;
		for(int iSide=0; iSide<2; iSide++) {
			int side = BRA;
			if (iSide == 1) side = KET;
			for(int iSQ=0; iSQ<(int)inputSQList.size(); iSQ++) {
				const ShellQuartet& sq = inputSQList[iSQ];
				if (sq.canDoHRR(side)) {
					if (iSide == 0) {
						braCanDOHRR = true;
					}else{
						ketCanDOHRR = true;
					}
					break;
				}
			}
		}

		// now let's see whether both side are available?
		// we always choose the side that have less number
		// of integrals as the first side, since we need 
		// to carry the whole number of integral for the 
		// first side to the second side
		if (braCanDOHRR && ketCanDOHRR) {
			if (nIntsBraSide<=nIntsKetSide) {
				firstSide  = BRA;
				secondSide = KET;
			}else{
				firstSide  = KET;
				secondSide = BRA;
			}
			return;
		}

		// now there must be one side is not available 
		// for HRR
		// so let's see the situation
		if (!braCanDOHRR && ketCanDOHRR) {
			firstSide  = KET;
			secondSide = NULL_POS;
		}else if (!ketCanDOHRR && braCanDOHRR) {
			firstSide  = BRA;
			secondSide = NULL_POS;
		}else {
			firstSide  = NULL_POS;
			secondSide = NULL_POS;
		}
	}
}	

int RR::countLHSIntNumbers() const 
{
	// we note that the counting may contain some NULL LHS
	// however, this is only for approximation
	int totalNInts = 0;
	for(list<RRSQ>::const_iterator it=rrsqList.begin(); it!=rrsqList.end(); ++it) {
		const list<int>& LHS = it->getLHSIndexArray();
		totalNInts += LHS.size();
	}
	return totalNInts;
}

int RR::countRHSIntNumbers() const 
{
	int totalNInts = 0;
	for(list<RRSQ>::const_iterator it=rrsqList.begin(); it!=rrsqList.end(); ++it) {
		totalNInts += it->countRHSIntegralNum();
	}
	return totalNInts;
}

