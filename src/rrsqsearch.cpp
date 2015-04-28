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
#include "rrbuild.h"
#include "inttype.h"
#include "rrsqsearch.h"
using namespace inttype;
using namespace rrbuild;
using namespace rrsqsearch;

////////////////////////////////////////////////////////////////////////////////
//       ####          class related to RRShellQuart                          //
////////////////////////////////////////////////////////////////////////////////
RRShellQuart::RRShellQuart(const int& rrType, const ShellQuartet& sq):oriSQ(sq)
{
	// operator information
	int oper  = oriSQ.getOper();
	int order = getOperOrder(oper);

	// we try to get the length of items in the RR - just use BRA1 postition
	nPos = oriSQ.getNNonSShell();
	crash(nPos == 0,"The number of Non-S shells in RRSQ is zero");
	RRBuild testRR(rrType,oper,BRA1);
	nRRItems = testRR.getNItems();
	int n = nRRItems*nPos;   
	sqlist.reserve(n);
	posList.reserve(nPos);

	// now step into the real work code 
	for(int iOrder=0; iOrder<order; iOrder++) {

		// get the possible position
		int pos = BRA1;
		if (iOrder == 1) {
			pos = BRA2;
		} else if (iOrder == 2) {
			pos = KET1;
		} else if (iOrder == 3) {
			pos = KET2;
		}

		// testing whether we could form expansion on this pos
		// for S type of position or null position, we do nothing
		if(oriSQ.isNonSShell(pos)) {
			posList.push_back(pos);
		}else{
			continue;
		}	

		// form the general RR used in this RRSQ
		RRBuild generalRR(rrType,oper,pos);

		// form the possible shell quartet list in the RHS
		// push the sq into the sqlist
		generalRR.buildRRSQ(oriSQ,sqlist);
	}

	// finally, process the state and integral number for the sqlist
	n = sqlist.size();
	sqStateList.assign(n,NEW_SQ);  
	intNumList.assign(n,0);
	for(int i=0; i<n; i++) {
		if (sqlist[i].isnull()) {
			sqStateList[i] = NULL_SQ;
		}else{
			intNumList[i] = sqlist[i].getNInts();
		}
	}
}

/*
 * this is for the HRR use
 * finally we abandon it
RRShellQuart::RRShellQuart(const ShellQuartet& sq, const int& side):oriSQ(sq)
{
	// checking the side information
#ifdef RRSEARCH_DEBUG
	crash(side != BRA && side != KET, "Only BRA/KET are allowed in RRShellQuartet");
#endif

	// setting the position information etc.
	nPos = 2;
	nRRItems = 2; // this is constant for HRR
	int n = nRRItems*nPos;   
	sqlist.reserve(n);
	posList.reserve(nPos);
	int oper  = oriSQ.getOper();

	// now step into the real work code 
	for(int i=0; i<nPos; i++) {

		// get the possible position
		int pos = BRA1;
		if (side == KET) pos = KET1;
		if (i == 1) {
			pos = BRA2;
			if (side == KET) pos = KET2;
		}

		// testing whether we could form expansion on this pos
		if(oriSQ.isNonSShell(pos)) {
			posList.push_back(pos);
		}else{
			continue;
		}	

		// form the general RR used in this RRSQ
		RRBuild generalRR(HRR,oper,pos);

		// form the possible shell quartet list in the RHS
		// push the sq into the sqlist
		generalRR.buildRRSQ(oriSQ,sqlist);
	}

	// finally, process the state and integral number for the sqlist
	n = sqlist.size();
	sqStateList.assign(n,NEW_SQ);  
	intNumList.assign(n,0);
	for(int i=0; i<n; i++) {
		if (sqlist[i].isnull()) {
			sqStateList[i] = NULL_SQ;
		}else{
			intNumList[i] = sqlist[i].getNInts();
		}
	}
}
*/

void RRShellQuart::updateStateSQList(const list<ShellQuartet>& sqset)
{
	int n = sqlist.size();
	list<ShellQuartet>::const_iterator it;
	for(int i=0; i<n; i++) {
		if (sqlist[i].isnull()) continue;
		it = find(sqset.begin(),sqset.end(),sqlist[i]);
		if (it != sqset.end()) {
			sqStateList[i] = OLD_SQ;
		}
	}
}

void RRShellQuart::updateStateSQList(const vector<ShellQuartet>& sqset)
{
	int n = sqlist.size();
	vector<ShellQuartet>::const_iterator it;
	for(int i=0; i<n; i++) {
		if (sqlist[i].isnull()) continue;
		it = find(sqset.begin(),sqset.end(),sqlist[i]);
		if (it != sqset.end()) {
			sqStateList[i] = OLD_SQ;
		}
	}
}

void RRShellQuart::getMinIntPos(int& minPos, int& minIntNumber) const
{
	minPos       = 0;
	minIntNumber = -1;
	for(int iPos=0; iPos<nPos; iPos++) {

		// estimate integral number for the given position
		int minInt = 0;
		for(int iSQ=0; iSQ<nRRItems; iSQ++) {
			int index = iSQ+iPos*nRRItems;
			if (sqStateList[index] == NEW_SQ) {
				int intNum = intNumList[index];
				minInt += intNum;
			}
		}

		// do we get a smaller position?
		if (minIntNumber<0) {
			minIntNumber = minInt;
		}else{
			if (minIntNumber > minInt) {
				minPos = iPos;
				minIntNumber = minInt;
			}
		}
	}
}

void RRShellQuart::updateRHSSQList(list<ShellQuartet>& rhsList, const int& pos) const
{
	int p = pos*nRRItems;
	for(int i=0; i<nRRItems; i++) {

		// test the given position sq, they should be new one 
		// and non-S shell quartet
		// the null sq should be tested first
		const ShellQuartet& sq = sqlist[p+i];
		if (sq.isnull()) continue;
		bool isSTypeSQ = sq.isSTypeSQ();
		if (sqStateList[p+i] != NEW_SQ || isSTypeSQ) continue;

		// now let's go to see whether this shell quartet is 
		// already contained in rhsList
		list<ShellQuartet>::const_iterator it;
		it = find(rhsList.begin(),rhsList.end(),sq);
		if (it == rhsList.end()) {
			rhsList.push_back(sq);
		}
	}
}

void RRShellQuart::print() const
{
	cout << "LHS Shell quartet is " << oriSQ.getName() << endl;
	cout << "total number of possible non-S positions are " << nPos << endl;
	for(int iPos=0; iPos<nPos; iPos++) {
		cout << "position is " << posList[iPos] << endl;
		cout << "For this position, the possible RR RHS shell quartets:" << endl;
		for(int i=0; i<nRRItems; i++) {
			int index = iPos*nRRItems + i;
			const ShellQuartet& sq = sqlist[index];
			int intnumber = intNumList[index];
			int state     = sqStateList[index];
			if (state == NULL_SQ) {
				cout << "NULL SQ" << endl;
			}else if (state == OLD_SQ) {
				cout << sq.getName() << " OLD " << " integral number " << intnumber << endl;
			}else{
				cout << sq.getName() << " NEW " << " integral number " << intnumber << endl;
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
//       ####          class related to RRSQSearch (FOR VRR)                  //
////////////////////////////////////////////////////////////////////////////////
void RRSQSearch::searchOptPos(const list<ShellQuartet>& unsolvedMainSQList, 
		const list<ShellQuartet>& unsolvedAppendSQList)  
{

	//
	// This function is the key function to search the optimum position
	// for a given bunch of shell quartets. The basic idea is just loop
	// over all of possible position arragements.
	// The procedure is like below:
	// 0  combine the main list and append list together to form 
	//    the complete list;
	// 1  convert the original shell quartets into the rrSQ list
	//    so to obtain the RR expansion information;
	// 2  search all of possible position arrangements. Here
	//		we use an array to identify all of possible arrangement. 
	//		It's called loop_identifier.
	//	3  update the solvedSQlist and the unsolvedSQArch
	//	   so that to make the unsolvedSQList only contains the 
	//	   newly created sqs
	//
	
	// form the complete sq list for searching
	list<ShellQuartet> unsolvedSQList(unsolvedMainSQList);
	list<ShellQuartet>::const_iterator it;
	for(it=unsolvedAppendSQList.begin(); it != unsolvedAppendSQList.end(); ++it) {
		unsolvedSQList.push_back(*it);
	}

	// convert the original shell quartet into the rrsq
	// here we note that the RHS of rrsq would be compared
	// with both solved shell quartets and unsolved shell quartet archive.
	// In this step, we want to make sure that in the following 
	// opt search step, all of RHS shell quartets in the RR step
	// are all new shell quartets and they are not in the history
	//
	vector<RRShellQuart> rrSQList;
	int nSQ = unsolvedSQList.size();
	rrSQList.reserve(nSQ);
	for(it=unsolvedSQList.begin(); it != unsolvedSQList.end(); ++it) {

		// generate RR shell quartet
		RRShellQuart rrSQ(rrType,*it);

		// before we use, we need to upate the state of sqs inside rrSQ
		// this is "history(vertical) comparison"
		if (solvedSQList.size() > 0) {
			rrSQ.updateStateSQList(solvedSQList); // update the state with solved sq
		}
		rrSQ.updateStateSQList(unsolvedSQArch);  // update the state with unsolved sq archive

		// now push in the result rrsq
		//rrSQ.print();
		rrSQList.push_back(rrSQ);
	}

	// 
	// loop_max_pos records the maximum position number
	// for each shell quartet. For example, for shell
	// quartet 1 its non-S shell positions are
	// BRA1  BRA2  KET1 
	// then the getNPos() returns 3
	// and the maximum position number is 3-1 = 2
	// that means, the position could be 0 1 2 for 
	// shell quartet 1
	//
	vector<int> loop_max_pos(nSQ,0);
	for(int iSQ=0; iSQ<nSQ; iSQ++) {
		loop_max_pos[iSQ] = rrSQList[iSQ].getNPos()-1;
	}

	// initilize the search result 
	// we try to pick up the minimum int num for 
	// all possible position in each rrsq
	vector<int> optPosList(nSQ,0);
	int oriIntCount = 0;
	for(int iSQ=0; iSQ<nSQ; iSQ++) {
		int pos;
		int n;
		rrSQList[iSQ].getMinIntPos(pos,n);
		optPosList[iSQ] = pos;
		oriIntCount += n;
	}

	//
	// now let's search the optimum solution
	// the basic search process is that we set up a mapping between
	// all possible RR expansion arrangements with the vector of 
	// "loop_identifier". For example, if we have 2 shell quartets
	// and the possible RR expansion for them is given as:
	// sq1 :  BRA1(pos 0), BRA2(pos 1), and KET1(pos 2)
	// sq2 :  BRA1(pos 0), BRA2(pos 1)
	// then the loop_identifier is a vector containg 2 elements,
	// the first element in loop_identifier describes the sq1's 
	// possible pos, and the second element describes the sq2's 
	// possible pos. Therefore, 
	// loop_identifier[0][1] means the RR expansion for sq1 is 
	// BRA1, for sq2 is BRA2. If we loop over all possible 
	// elements in loop_identifier, then we go through all 
	// possible arrangement combinations for the given shell
	// quartet list.
	//
	vector<int> loop_identifier(nSQ, 0);

	//
	// on the other hand, we also need to perform the "horizontal 
	// comparison" in the opt search step. Basically, in the opt
	// search for each step, we have a given position list for 
	// the rrSQList; each rrSQ then has a list SQs in the RHS 
	// of RR expansion. We evaluate the integral number of these
	// "NEW RHS SQ"(see the vertical comparison), and sum over all 
	// of the integral numbers for each rrSQ and for all of rrSQList. 
	// We will search for the minimum integral number case.
	//
	// however, these RHS sqs (each rrSQ will have nRRItems and we 
	// have nSQ rrSQ therefore the total RHS sq number is nRRItems*nSQ).
	// may be same with each other. In this case, we only need to pick up
	// one of them to count into the integral number. This is called 
	// "horizontal comparison" comparing with the vertiaal one. In this sense,
	// in estimation we only pick up these RHS shell quartets only once.
	// if they appear again, we just omit them. For performing this job,
	// we set up this list to mornitor the integal number calcualtion.
	//
	list<ShellQuartet> hsqlist; // contains the picked up sq in horizontal comparison

	// now we start the opt position search
	int curPos = 0;             // current working position in loop_identifier
	while(true) {

		// counting integral numbers
		// push the first rrSQ's content into hsqlist
		// also count in the first rrsq's integral number
		int intNum = 0;
		hsqlist.clear();
		int pos0 = loop_identifier[0];
		const RRShellQuart& rrsq0 = rrSQList[0];
		for(int iRHSSQ=0; iRHSSQ<rrsq0.getNRRItems(); iRHSSQ++) {
			if (rrsq0.isNewSQ(pos0,iRHSSQ)) {
				hsqlist.push_back(rrsq0.getRHSSQ(pos0,iRHSSQ));
				intNum += rrsq0.intNumCount(pos0,iRHSSQ);
			}
		}
		
		// now counting the rest of rrsq
		// we will perform horizontal comparison
		list<ShellQuartet>::const_iterator it;
		for(int iSQ=1; iSQ<nSQ; iSQ++) {
			int pos = loop_identifier[iSQ];
			const RRShellQuart& rrsq = rrSQList[iSQ];
			for(int iRHSSQ=0; iRHSSQ<rrsq.getNRRItems(); iRHSSQ++) {
				if (rrsq.isNewSQ(pos,iRHSSQ)) {

					// now test that whether this sq is already in hsqlist
					const ShellQuartet& sq = rrsq.getRHSSQ(pos,iRHSSQ);
					it = find(hsqlist.begin(),hsqlist.end(),sq);
					if (it == hsqlist.end()) {
						hsqlist.push_back(sq);
						intNum += rrsq.intNumCount(pos,iRHSSQ);
					}
				}
			}
		}

		// whether we are able to get the better arrangement?
		if (intNum<oriIntCount) {
			optPosList = loop_identifier;
			oriIntCount = intNum;
		}

		// increase position
		if (loop_identifier[curPos] < loop_max_pos[curPos]) {
			loop_identifier[curPos] += 1;
		}else{

			// let's go to see whether we have loop over all possible arragements
			if (loop_identifier[nSQ-1] == loop_max_pos[nSQ-1]) {
				bool finish = true;
				for(int i=nSQ-2; i>=0; i--) {
					if (loop_identifier[i] != loop_max_pos[i]) {
						finish = false;
						break;
					}
				}
				if (finish) break;
			}

			// carry the position in the loop to the next possible one
			int iPos = curPos+1;
			while(iPos<nSQ) {
				if (loop_identifier[iPos] < loop_max_pos[iPos]) {
					loop_identifier[iPos] += 1;
					break;
				}
				iPos++;
			}

			// clear all of the position in front of iPos
			// so that to start a new round of searching
			for(int i=0; i<iPos; i++) {
				loop_identifier[i] = 0;
			}

			// set current position to the original
			curPos = 0;
		}
	}

	// debug print of the integral number
#ifdef RRSEARCH_DEBUG
	cout << "the minimum integral numbers in searchOptPos " << oriIntCount << endl;
#endif

	// now update the solvedSQList as well as the posList
	// only for the main ones
	for(int iSQ=0; iSQ<nSQ; iSQ++) {
		const RRShellQuart& rrsq = rrSQList[iSQ];
		const ShellQuartet& sq = rrsq.getLHSSQ();

		// is it the main one?
		list<ShellQuartet>::const_iterator it = find(unsolvedMainSQList.begin(),
				unsolvedMainSQList.end(),sq);
		if (it != unsolvedMainSQList.end()) {
			solvedSQList.push_back(sq);
			int var = rrsq.getVar(optPosList[iSQ]);
			posList.push_back(var);
		}
	}

	// now let's update the unsolvedSQArchive
	// we note that since each shell quartet in the RHS
	// is already been compared vertically, therefore 
	// they are new ones in the unsolvedSQArchive
	// we push them directly
	for(int iSQ=0; iSQ<nSQ; iSQ++) {
		const RRShellQuart& rrsq = rrSQList[iSQ];
		const ShellQuartet& sq = rrsq.getLHSSQ();

		// is it the main one?
		list<ShellQuartet>::const_iterator it = find(unsolvedMainSQList.begin(),
				unsolvedMainSQList.end(),sq);
		if (it != unsolvedMainSQList.end()) {
			int finalPos = optPosList[iSQ];
			rrsq.updateRHSSQList(unsolvedSQArch,finalPos);
		}
	}
}

bool RRSQSearch::canDoDirectParse(const ShellQuartet& sq)
{
	// if the given sq is null one or the pure S type
	// shell quartet, then we do not do anything
	// just ignore this shell quartet by returning true
	if (sq.isnull() || sq.isSTypeSQ()) {
		return true;
	}

	// determine the non-S positions
	bool canDO = sq.canDODirectRRPosSearch();
	if (! canDO) return false;

	// if we can do it, than do it here
	int pos = sq.getRRPos();
	posList.push_back(pos);
	solvedSQList.push_back(sq);
	return true;
}

bool RRSQSearch::canDoOptSearch(const list<ShellQuartet>& sqList,
		const list<ShellQuartet>& appendSQList) const 
{
	// for the given sqlist, we try to estimate the number of total
	// possible RR combinations
	// we note, that the number of 1073741824 is 4^15
	long long totalN = 1;
	list<ShellQuartet>::const_iterator it;
	for(it = sqList.begin(); it != sqList.end(); ++it) {
		int nPos = it->getNNonSShell();
		if (nPos>=1) totalN *= nPos;
		if (totalN > 1073741824) return false;
	}
	for(it = appendSQList.begin(); it != appendSQList.end(); ++it) {
		int nPos = it->getNNonSShell();
		if (nPos>=1) totalN *= nPos;
		if (totalN > 1073741824) return false;
	}
	return true;
}

void RRSQSearch::updateUnsolvedSQArch() 
{
	list<ShellQuartet>::iterator it;
	for(int iSQ=0; iSQ<(int)solvedSQList.size(); iSQ++) {
		it = find(unsolvedSQArch.begin(),unsolvedSQArch.end(),solvedSQList[iSQ]);
		if (it != unsolvedSQArch.end()) {
			unsolvedSQArch.erase(it);
		}
	}
}

RRSQSearch::RRSQSearch(const vector<ShellQuartet>& inputSQList, 
		const int& rrType0):rrType(rrType0)
{
	// reserve space for the result solved shell quartets
	// 10000 should be big enough for to this
	int bigNumber = 10000;
	solvedSQList.reserve(bigNumber);
	posList.reserve(bigNumber);

	// now initilize the unsolved sq archive
	// it's used to hold the unsolved SQs for the whole process
	int nInitSQ = inputSQList.size();
	for(int iSQ=0; iSQ<nInitSQ; iSQ++) {
		unsolvedSQArch.push_back(inputSQList[iSQ]);
	}

	// also we check the operator of input sq
	// we note, that in the VRR process, the input
	// sq should have the same kind of operators
	int oper = inputSQList[0].getOper();
	for(int iSQ=0; iSQ<nInitSQ; iSQ++) {
		int o = inputSQList[iSQ].getOper();
		if (o != oper) {
			cout << "first operator " << oper << endl;
			cout << "second operator " << o << endl;
			crash(true, "In the RRSQSearch class the input sq have different operators");
		}
	}

	// now perform the search
	// we note that we have to select the search method
	int method = rrProp(oper);
	RRSearchBasedOnTWOProperty(method,inputSQList);
	//print(method,inputSQList);
}

void RRSQSearch::pickupUnsolvedSQByLM(const int& P1,
		const int& P2, list<ShellQuartet>& unsolvedMainSQList, 
		list<ShellQuartet>& unsolvedAppendSQList) 
{
	// assign the L and M value
	int L = P1;
	int M = P2;

	// firstly, form the main shell quartets
	list<ShellQuartet>::const_iterator it;
	for(it = unsolvedSQArch.begin(); it != unsolvedSQArch.end(); ++it) {
		if (it->matchLM(L,M)) {
			// can we do direct parse?
			// if so, the solvedSQList would be updated here
			bool canDo = canDoDirectParse(*it);
			if (! canDo) {
				unsolvedMainSQList.push_back(*it);
			}
		}
	}

	// did we find anything?
	if (unsolvedMainSQList.size() == 0) return;

	//
	// now let's form the appended shell quartets
	// for L and M critiria, the appended shell quartets
	// has L between L and L-1
	// We note that sq with L-2 does not have RHS overlap
	// with the main shell quartets
	// M in RHS has M and M+1, they both has possibility
	// in overlaping with main shell quartets
	//
	//
	int nL = 2;
	int nM = 2;
	vector<int> Llist(nL,0); 
	vector<int> Mlist(nM,0);
	Llist[0] = L;
	if (L-1>=0) {
		Llist[1] = L-1;
	}else{
		nL = nL - 1;
	}
	Mlist[0] = M;
	Mlist[1] = M+1;

	// now let's process the LM combinations
	for(int iL=0; iL<nL; iL++) {
		for(int iM=0; iM<nM; iM++) {

			// we do not counting the main ones
			int l = Llist[iL];
			int m = Mlist[iM];
			if (l == L && m == M) continue;

			// search
			for(it = unsolvedSQArch.begin(); it != unsolvedSQArch.end(); ++it) {
				if (it->matchLM(l,m)) {
					bool canDo = canDoDirectParse(*it);
					if (! canDo) {
						unsolvedAppendSQList.push_back(*it);
					}
				}
			}
		}
	}
}

void RRSQSearch::pickupUnsolvedSQByL(const int& P, list<ShellQuartet>& unsolvedMainSQList, 
		list<ShellQuartet>& unsolvedAppendSQList) 
{
	// assign the L value
	int L = P;

	// firstly, form the main shell quartets
	list<ShellQuartet>::const_iterator it;
	for(it = unsolvedSQArch.begin(); it != unsolvedSQArch.end(); ++it) {
		if (it->matchL(L)) {
			// can we do direct parse?
			// if so, the solvedSQList would be updated here
			bool canDo = canDoDirectParse(*it);
			if (! canDo) {
				unsolvedMainSQList.push_back(*it);
			}
		}
	}

	// did we find anything?
	if (unsolvedMainSQList.size() == 0) return;

	//
	// now let's form the appended shell quartets
	// for L critiria, the appended shell quartets
	// has L between L and L-1
	// We note that sq with L-2 does not have RHS overlap
	// with the main shell quartets
	//
	//
	int nL = 2;
	vector<int> Llist(nL,0); 
	Llist[0] = L;
	if (L-1>=0) {
		Llist[1] = L-1;
	}else{
		nL = nL - 1;
	}

	// now let's process the LM combinations
	for(int iL=0; iL<nL; iL++) {

		// we do not counting the main ones
		int l = Llist[iL];
		if (l == L) continue;

		// search
		for(it = unsolvedSQArch.begin(); it != unsolvedSQArch.end(); ++it) {
			if (it->matchL(l)) {
				bool canDo = canDoDirectParse(*it);
				if (! canDo) {
					unsolvedAppendSQList.push_back(*it);
				}
			}
		}
	}
}

void RRSQSearch::pickupUnsolvedSQByLOper(const int& P1,
		const int& P2, list<ShellQuartet>& unsolvedMainSQList, 
		list<ShellQuartet>& unsolvedAppendSQList) 
{
	// assign the L and Operator value
	int L = P1;
	int O = P2;

	// for the RR expression, kinetic integrals will generate
	// the kinetic integral and two body overlap integral on 
	// the RHS of RR. If we step into the two body overlap
	// integrals, things becomes simple therefore we will
	// call the function of pickupUnsolvedSQByL. Since 
	// overlap only generate overlap integrals on LHS
	int rrprop = rrProp(O);
	if (rrprop == SQ_ON_L) {
		pickupUnsolvedSQByL(P1,unsolvedMainSQList,unsolvedAppendSQList);
		return;
	}

	// now let's concentrate on operators which generates 
	// several kind of integrals on RHS
	// firstly, form the main shell quartets
	list<ShellQuartet>::const_iterator it;
	for(it = unsolvedSQArch.begin(); it != unsolvedSQArch.end(); ++it) {
		if (it->matchLOper(L,O)) {
			// can we do direct parse?
			// if so, the solvedSQList would be updated here
			bool canDo = canDoDirectParse(*it);
			if (! canDo) {
				unsolvedMainSQList.push_back(*it);
			}
		}
	}

	// did we find anything?
	if (unsolvedMainSQList.size() == 0) return;

	//
	// we note, that for kinetic RR expansion since the RHS and RR has L, 
	// L-1, L-2 three kind of terms so that in considering the RHS overlap
	// case we must consider the friend shell quartets of L-1 and L-2
	// since they may overlap with main shell quartet on RHS of RR
	//
	vector<int> Llist; 
	vector<int> Olist; 
	if (O == KINETIC) {
		Llist.reserve(3);
		Llist.push_back(L);
		if (L-1>=0) {
			Llist.push_back(L-1);
		}
		if (L-2>=0) {
			Llist.push_back(L-2);
		}
	}else{
		crash(true,"Operator is not valid in RRSQSearch of pickupUnsolvedSQ");
	}
	selectOperListInRR(O, Olist);

	// now let's process the LM combinations
	for(int iL=0; iL<(int)Llist.size(); iL++) {
		for(int iO=0; iO<(int)Olist.size(); iO++) {
			int l    = Llist[iL];
			int oper = Olist[iO];

			// omit the main ones
			if (l == L && oper == O) continue;

			// search
			for(it = unsolvedSQArch.begin(); it != unsolvedSQArch.end(); ++it) {
				if (it->matchLOper(l,oper)) {
					bool canDo = canDoDirectParse(*it);
					if (! canDo) {
						unsolvedAppendSQList.push_back(*it);
					}
				}
			}
		}
	}
}

void RRSQSearch::pickupUnsolvedSQ(const int& pickupMethod, const int& P1,
		const int& P2, list<ShellQuartet>& unsolvedMainSQList, 
		list<ShellQuartet>& unsolvedAppendSQList) 
{
	//
	// Generally, this function is used to pick up the unsolved shell quartets
	// for the following searchOptPos function. It trys to answer this question,
	// that for some given shell quartet properties(for example, fixed L and M
	// value), what kind of shell quartets could be picked up together?
	//
	// In picking up process, we divide the shell quartets into two groups:
	// one group is called main shell quartet, they math the given property
	// 1 and 2;
	// the other group is called appended shell quartet. They did not match
	// the given property, however, such shell quartet are potentially has
	// the same RHS shell quartets with the main ones. Therefore, they 
	// are overlapped with the main shell quartets in this sense.
	//
	//

	if (pickupMethod == SQ_ON_L_M) {
		pickupUnsolvedSQByLM(P1,P2,unsolvedMainSQList,unsolvedAppendSQList);
	}else if (pickupMethod == SQ_ON_L_OPER) {
		pickupUnsolvedSQByLOper(P1,P2,unsolvedMainSQList,unsolvedAppendSQList);
	}else if (pickupMethod == SQ_ON_L) {
		pickupUnsolvedSQByL(P1,unsolvedMainSQList,unsolvedAppendSQList);
	}else{
		crash(true, "Incorrect pickup method in pickupUnsolvedSQ");
	}
}

void RRSQSearch::RRSearchBasedOnTWOProperty(const int& method,
		const vector<ShellQuartet>& initSQList)
{
	//
	// This routine is doing the RR position search search for 
	// the RR expansion that only two shell quartet properties 
	// involved.
	//
	// The algorithm for the search is like below:
	// 1  shell properties pick up, usually it's L and M combination
	//    or the L and operator combination. We note, this case
	//    is also used for the single shell property search process.
	//    For example, the RR expansion for overlap integrals.
	// 2  form the loop based on the two shell quartet properties;
	// 3  form the unsolvedSQList based on the fixed shell quartet properties;
	// 4  do the core function - searchOptPos
	// 5  update results - we need to make sure that unsolvedSQArch
	//    and solveSQList are not overlap with each other
	//

	// get shell quartet properties
	// one shell property is always L
	int maxL    = 0;
	int nInitSQ = initSQList.size(); 
	for(int iSQ=0; iSQ<nInitSQ; iSQ++) {
		int L = initSQList[iSQ].getLSum(); 
		if (L > maxL) maxL = L;
	}

	// the other property could be M or operator type
	vector<int> Plist;
	if (method == SQ_ON_L_M) {

		// the possible M cases is from 0 to maxM(maxL)
		// so total number is maxM+1
		int maxM   = maxL;
		int totalM = maxM+1;
		Plist.assign(totalM,0);
		for(int iM=1; iM<=maxM; iM++) {
			Plist[iM] = iM;
		}

	}else if (method == SQ_ON_L_OPER) {
		
		int O = initSQList[0].getOper();
		selectOperListInRR(O, Plist);
	
	}else if (method == SQ_ON_L) {
		
		//in this case, we actually do not have second
		//property defined. We just set it to be 0 
		Plist.assign(1,0);
	}

	// real work begins...
	for(int iP=0; iP<(int)Plist.size(); iP++) {
		for(int iL=maxL; iL>=0; iL--) {

			// fetch the second shell quartet property
			int property = Plist[iP];

			// step into the real loop - this is for each M value
			while(true) {

#ifdef RRSEARCH_DEBUG
				if (method == SQ_ON_L_OPER) {
					string operatorName = getOperStringName(property);
					cout << "each round of search, for L " << iL << " and operator " 
						<< operatorName << endl;
				}else{
					cout << "each round of search, for L " << iL << " and M value " 
						<< property << endl;
				}
#endif

				// for each searching process, what we use is the 
				// unsolved main sq list and unsolved append sq list 
				// to start the search
				// the main list is the target one
				// the append list is used to make the target one more
				// accurate
				list<ShellQuartet> unsolvedMainSQList; 
				list<ShellQuartet> unsolvedAppendSQList; 
				pickupUnsolvedSQ(method,iL,property,unsolvedMainSQList,unsolvedAppendSQList);
#ifdef RRSEARCH_DEBUG
				if (unsolvedMainSQList.size() > 0) {
					cout << "unsolved main sq list before search opt pos " << endl;
					for(list<ShellQuartet>::const_iterator it = unsolvedMainSQList.begin(); 
							it != unsolvedMainSQList.end(); ++it) {
						cout << it->getName() << endl;
					}
					cout << "unsolved append sq list before search opt pos " << endl;
					for(list<ShellQuartet>::const_iterator it = unsolvedAppendSQList.begin(); 
							it != unsolvedAppendSQList.end(); ++it) {
						cout << it->getName() << endl;
					}
					cout << endl;
				}
#endif

				// shall we just stop here for this loop?
				int nUnsolvedSQ = unsolvedMainSQList.size();
				if (nUnsolvedSQ == 0) {
					updateUnsolvedSQArch();
					break;
				}

				// before we do opt search
				// we need to keep an eye on the actual
				// number of shell quartets
				if (! canDoOptSearch(unsolvedMainSQList,unsolvedAppendSQList)) {
					crash(true, "the number of possible RR expansions is too large in RR research");
				}

				// now let's do opt search
				// we note, that the unsolved sq list should be 
				// updated inside this function so that unsolved
				// shell quartet list only contain the new created
				// shell quartets
				searchOptPos(unsolvedMainSQList,unsolvedAppendSQList);

				// update the sq archive list by comparing with the 
				// solvedSQList. sq archive should only contain these
				// unsolved sqs
				updateUnsolvedSQArch();

				// debug printing
#ifdef RRSEARCH_DEBUG
				cout << "number of solved shell quartets " << solvedSQList.size() << endl;
				for(int iSQ=0; iSQ<(int)solvedSQList.size(); iSQ++) {
					cout << solvedSQList[iSQ].getName() << " pos " << posList[iSQ] << endl;
				}

				cout << "number of unsolved shell quartets " << unsolvedSQArch.size() << endl;
				for(list<ShellQuartet>::const_iterator it=unsolvedSQArch.begin(); 
						it != unsolvedSQArch.end(); ++it) {
					cout << it->getName() << endl;
				}
				cout << endl;
#endif
			}
		}
	}
}

void RRSQSearch::print(const int& method, const vector<ShellQuartet>& initSQList) const
{

	// get shell quartet properties
	// one shell property is always L
	int maxL    = 0;
	int nInitSQ = initSQList.size(); 
	for(int iSQ=0; iSQ<nInitSQ; iSQ++) {
		int L = initSQList[iSQ].getLSum(); 
		if (L > maxL) maxL = L;
	}

	// the other property could be M or operator type
	// the code is copied from RR process based on 
	// two property value
	vector<int> Plist;
	if (method == SQ_ON_L_M) {
		// the possible M cases is from 0 to maxM(maxL)
		// so total number is maxM+1
		int maxM   = maxL;
		int totalM = maxM+1;
		Plist.assign(totalM,0);
		for(int iM=1; iM<=maxM; iM++) {
			Plist[iM] = iM;
		}
	}else if (method == SQ_ON_L_OPER) {
		int O = initSQList[0].getOper();
		selectOperListInRR(O, Plist);
	}else if (method == SQ_ON_L) {
		Plist.assign(1,0);
	}

	// printing work
	cout << "Result shell quartet printing" << endl;
	int totalInts = 0;
	for(int iP=0; iP<(int)Plist.size(); iP++) {
		for(int iL=maxL; iL>=0; iL--) {

			// property and L
			int L = iL;
			int P = Plist[iP];

			// check the solved list
			int nSQ = solvedSQList.size();
			for(int i=0; i<nSQ; i++) {
				const ShellQuartet& sq = solvedSQList[i];
				int pos = posList[i];

				// see whether we get the match sq
				bool match = false;
				if (method == SQ_ON_L_M) {
					if (sq.matchLM(L,P)) match = true;
				}else if (method == SQ_ON_L_OPER) {
					if (sq.matchLOper(L,P)) match = true;
				}else if (method == SQ_ON_L) {
					if (sq.matchL(L)) match = true;
				}else{
					crash(true,"Incorrect method passed in print function of RRSQSearch of VRR");
				}

				// print
				if (match) {
					totalInts += sq.getNInts();
					cout << sq.getName() << " RR expansion pos " << pos << " " 
						<< sq.getNInts() << endl;
				}
			}
		}
	}
	cout << "total integral number is " << totalInts << endl;
}

////////////////////////////////////////////////////////////////////////////////
//       ####          class related to RRSQSearch (FOR HRR)                  //
////////////////////////////////////////////////////////////////////////////////
bool RRSQSearch::isBottomForHRR(const ShellQuartet& sq, const int& side)
{
	// firstly, let's see whether the side is already
	// some given position?
	if (side == BRA1 || side == BRA2 || side == KET1 || side == KET2) {
		int pos = side;
		const Shell& s = sq.getShell(pos);
		if (s.getL() > 0) return false;
		return true;
	}

	// now the side is only bra/ket
	// how many non-S position it has?
	int nNonSPos = 0;
	if (side == BRA) {
		const Shell& b1 = sq.getShell(BRA1);
		const Shell& b2 = sq.getShell(BRA2);
		if (b1.getL() > 0) nNonSPos++;
		if (b2.getL() > 0) nNonSPos++;
	}else if (side == KET) {
		const Shell& k1 = sq.getShell(KET1);
		const Shell& k2 = sq.getShell(KET2);
		if (k1.getL() > 0) nNonSPos++;
		if (k2.getL() > 0) nNonSPos++;
	}else{
		cout << "side is " << side << endl;
		crash(true, "Illegal side pass in the isBottomforHRR function");
	}

	// now let's say whether it's the bottom?
	if (nNonSPos < 2) {
		return true;
	}	
	return false;
}

RRSQSearch::RRSQSearch(const int& side, const vector<ShellQuartet>& inputSQList):rrType(HRR)
{
	// reserve space for the result solved shell quartets
	// 10000 should be big enough for to this
	int bigNumber = 10000;
	solvedSQList.reserve(bigNumber);
	posList.reserve(bigNumber);

	// now perform the search
	HRRSearch(side,inputSQList);
	//print();
	
}

void RRSQSearch::HRRSearch(const int& side, const vector<ShellQuartet>& inputList) 
{
	//
	// This is the top routine for doing HRR search
	// we do HRR search for the given side
	//
	
	// firstly, let's decide whether we do direct RR search
	// or indirect HRR search
	bool doDirectWork = true;
	if (side == BRA || side == KET) doDirectWork = false;
	
	// now do the search
	if (doDirectWork) {

		// search process
		int pos = side;
		list<ShellQuartet> outputList;
		directHRRSearch(pos,inputList,outputList);

		// now add in the output result into solved list
		for(list<ShellQuartet>::const_iterator it = outputList.begin(); 
				it != outputList.end(); ++it) {
			solvedSQList.push_back(*it);
			posList.push_back(pos);
		}

	}else{

		// search process
		list<ShellQuartet> outputList;
		int resultPos = -1;
		resultPos = indirectHRRSearch(side,inputList,outputList);

		// now add in the output result into solved list
		for(list<ShellQuartet>::const_iterator it = outputList.begin(); 
				it != outputList.end(); ++it) {
			solvedSQList.push_back(*it);
			posList.push_back(resultPos);
		}
	}
}

void RRSQSearch::directHRRSearch(const int& pos, const vector<ShellQuartet>& inputSQList, 
		list<ShellQuartet>& outputList)
{
	// get the oper
	// for HRR, operator should be same for the given group of 
	// shell quartets. Therefore we just use the first shell
	// quartet
	int oper = inputSQList[0].getOper();
	RRBuild generalRR(HRR,oper,pos);
#ifdef RRSEARCH_DEBUG
	cout << "In Current HRRSearch variable is " << pos << endl; 
#endif

	// set the input LHS sq list
	// we assume here the input list should not have the 
	// repeating shell quartets
	list<ShellQuartet> lhsSQList;
	for(int iSQ=0; iSQ<(int)inputSQList.size(); iSQ++) {
		lhsSQList.push_back(inputSQList[iSQ]);
	}

	// now the real work
	while(true) {

		// debug test
		// for each round of search, let's take a look of 
		// input shell quartets for LHS, and current output
		// shell quartet of output
#ifdef RRSEARCH_DEBUG
		cout << "the initilial LHS shell quartet" << endl;
		for(list<ShellQuartet>::const_iterator it = lhsSQList.begin(); 
				it != lhsSQList.end(); ++it) {
			cout << it->getName() << endl;
		}
		cout << endl << endl;
		cout << "the current result output shell quartet" << endl;
		for(list<ShellQuartet>::const_iterator it = outputList.begin(); 
				it != outputList.end(); ++it) {
			cout << it->getName() << endl;
		}
#endif

		// let's push the LHS SQ into the output result
		// we note that these input sq are already unique ones
		// since we have do the search work below to make sure
		// its uniqueness
		for(list<ShellQuartet>::const_iterator it = lhsSQList.begin(); 
				it != lhsSQList.end(); ++it) {
			outputList.push_back(*it);
		}

		// set up a list to hold RHS shell quartets in this round
		vector<ShellQuartet> sqlist;
		int num = lhsSQList.size()*2;
		if (num > 0) sqlist.reserve(num);

		// form the RHS for each LHS sq
		for(list<ShellQuartet>::const_iterator it = lhsSQList.begin(); 
				it != lhsSQList.end(); ++it) {
			if (isBottomForHRR(*it,pos)) {
				continue;
			}else{
				generalRR.buildRRSQ(*it,sqlist);
			}
		}

#ifdef RRSEARCH_DEBUG
		cout << "RHS shell quartet list " << sqlist.size() << endl;
		for(int iSQ=0; iSQ<(int)sqlist.size(); iSQ++) {
			cout << sqlist[iSQ].getName() << endl;
		}
#endif

		// whether it's time to step out?
		// which means all of sq in the lhsSQList
		// are bottom sq
		if (sqlist.size() == 0) break;

		// now clear the lhs SQ list
		// and add in the new RHS sq into it to create new lhs sq for next round
		// the new LHS sq should be unique in a new batch of lhssqlist
		// and it should be new to the output result sq list
		lhsSQList.clear();
		for(int iSQ=0; iSQ<(int)sqlist.size(); iSQ++) {
			const ShellQuartet& sq = sqlist[iSQ];
			list<ShellQuartet>::iterator it  = find(lhsSQList.begin(),lhsSQList.end(),sq);
			list<ShellQuartet>::iterator it2 = find(outputList.begin(),outputList.end(),sq);
			if (it == lhsSQList.end() && it2 == outputList.end()) {
				lhsSQList.push_back(sq);
			}
		}

#ifdef RRSEARCH_DEBUG
		cout << "shell quartet for next round of search " << endl;
		for(list<ShellQuartet>::const_iterator it = lhsSQList.begin(); 
				it != lhsSQList.end(); ++it) {
			cout << it->getName() << endl;
		}
#endif
	}
}

int RRSQSearch::indirectHRRSearch(const int& side, const vector<ShellQuartet>& inputSQList, 
		list<ShellQuartet>& outputList)
{
	// get the oper
	// for HRR, operator should be same for the given group of 
	// shell quartets. Therefore we just use the first shell
	// quartet
	int oper = inputSQList[0].getOper();

	// collecting the possible expansion position
	int nPos = 2;
	vector<int> potentialPosList(nPos);
	potentialPosList[0] = BRA1;
	potentialPosList[1] = BRA2;
	if (side == KET) {
		potentialPosList[0] = KET1;
		potentialPosList[1] = KET2;
	}
#ifdef RRSEARCH_DEBUG
	cout << "In Current HRRSearch side is " << side << endl; 
#endif

	// set up two result List
	// one is for first variable
	// the other is for second variable
	list<ShellQuartet> firstVarResultList;
	list<ShellQuartet> secondVarResultList;
	int firstIntNum = 0;
	int secondIntNum = 0;

	// real work begins...
	for(int iPos=0; iPos<nPos; iPos++) {

		// initilize the loop
		int var = potentialPosList[iPos];
		RRBuild generalRR(HRR,oper,var);
#ifdef RRSEARCH_DEBUG
		cout << "In Current HRRSearch variable is " << var << endl; 
#endif

		// set up the output list to hold this position's result
		list<ShellQuartet> output;

		// set up a list to hold each round LHS shell quartets
		// initilize for the first round of search
		// we assume here the input list should not have the 
		// repeating shell quartets
		list<ShellQuartet> lhsSQList;
		for(int iSQ=0; iSQ<(int)inputSQList.size(); iSQ++) {
			lhsSQList.push_back(inputSQList[iSQ]);
		}

		while(true) {

			// debug test
			// for each round of search, let's take a look of 
			// input shell quartets for LHS, and current output
			// shell quartet of output
#ifdef RRSEARCH_DEBUG
			cout << "the initilial LHS shell quartet" << endl;
			for(list<ShellQuartet>::const_iterator it = lhsSQList.begin(); 
					it != lhsSQList.end(); ++it) {
				cout << it->getName() << endl;
			}
			cout << endl << endl;
			cout << "the current result output shell quartet" << endl;
			for(list<ShellQuartet>::const_iterator it = output.begin(); 
					it != output.end(); ++it) {
				cout << it->getName() << endl;
			}
#endif

			// let's push the LHS SQ into the output result
			// we note that these input sq are already unique ones
			// since we have do the search work below to make sure
			// its uniqueness
			for(list<ShellQuartet>::const_iterator it = lhsSQList.begin(); 
					it != lhsSQList.end(); ++it) {
				output.push_back(*it);
			}

			// set up a list to hold RHS shell quartets in this round
			vector<ShellQuartet> sqlist;
			int num = lhsSQList.size()*2;
			if (num > 0) sqlist.reserve(num);

			// form the RHS for each LHS sq
			for(list<ShellQuartet>::const_iterator it = lhsSQList.begin(); 
					it != lhsSQList.end(); ++it) {
				if (isBottomForHRR(*it,side)) {
					continue;
				}else{
					generalRR.buildRRSQ(*it,sqlist);
				}
			}

#ifdef RRSEARCH_DEBUG
			cout << "RHS shell quartet list " << sqlist.size() << endl;
			for(int iSQ=0; iSQ<(int)sqlist.size(); iSQ++) {
				cout << sqlist[iSQ].getName() << endl;
			}
#endif

			// whether it's time to step out?
			// which means all of sq in the lhsSQList
			// are bottom sq
			if (sqlist.size() == 0) break;

			// now clear the lhs SQ list
			// and add in the new RHS sq into it to create new lhs sq for next round
			// the new LHS sq should be unique in a new batch of lhssqlist
			// and it should be new to the output result sq list
			lhsSQList.clear();
			for(int iSQ=0; iSQ<(int)sqlist.size(); iSQ++) {
				const ShellQuartet& sq = sqlist[iSQ];
				list<ShellQuartet>::iterator it  = find(lhsSQList.begin(),lhsSQList.end(),sq);
				list<ShellQuartet>::iterator it2 = find(output.begin(),output.end(),sq);
				if (it == lhsSQList.end() && it2 == output.end()) {
					lhsSQList.push_back(sq);
				}
			}

#ifdef RRSEARCH_DEBUG
			cout << "shell quartet for next round of search " << endl;
			for(list<ShellQuartet>::const_iterator it = lhsSQList.begin(); 
					it != lhsSQList.end(); ++it) {
				cout << it->getName() << endl;
			}
#endif
		}

		// add in output into the result list
		for(list<ShellQuartet>::const_iterator it = output.begin(); it != output.end(); ++it) {
			if (iPos == 0) {
				firstVarResultList.push_back(*it);
				firstIntNum += it->getNInts();
			}else{
				secondVarResultList.push_back(*it);
				secondIntNum += it->getNInts();
			}
		}

#ifdef RRSEARCH_DEBUG
		cout << "the result output shell quartet after HRR work is done" << endl;
		for(list<ShellQuartet>::const_iterator it = output.begin(); 
				it != output.end(); ++it) {
			cout << it->getName() << endl;
		}
#endif

	}

	// now it's time to determine which one is best for HRR expansion
	int resultPos = -1;
	if (firstIntNum < secondIntNum) {
		resultPos = potentialPosList[0];
		for(list<ShellQuartet>::const_iterator it = firstVarResultList.begin(); 
				it != firstVarResultList.end(); ++it) {
			outputList.push_back(*it);
		}
	}else{
		resultPos = potentialPosList[1];
		for(list<ShellQuartet>::const_iterator it = secondVarResultList.begin(); 
				it != secondVarResultList.end(); ++it) {
			outputList.push_back(*it);
		}
	}
	return resultPos;
}

void RRSQSearch::print() const
{
	cout << "Result shell quartet printing" << endl;
	int nSQ = solvedSQList.size();
	int totalInts = 0;
	for(int i=0; i<nSQ; i++) {
		const ShellQuartet& sq = solvedSQList[i];
		int pos  = posList[i];
		int ints = sq.getNInts();
		totalInts += ints;
		cout << sq.getName() << " RR expansion pos " << pos << " " << ints << endl;
	}
	cout << "total integral number is " << totalInts << endl;
}
