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
#include "rrbuild.h"
#include "integral.h"
#include "inttype.h"
#include "boost/lexical_cast.hpp"
#include "printing.h"
#include "sqintsinfor.h"
#include "derivinfor.h"
#include "rrints.h"
using boost::lexical_cast;
using namespace rrbuild;
using namespace integral;
using namespace inttype;
using namespace printing;
using namespace sqintsinfor;
using namespace derivinfor;
using namespace rrints;

RRSQ::RRSQ(const int& rrType0, const int& pos, const ShellQuartet& sq,
		const set<int>& unsolvedIntegralList, int dir):rrType(rrType0),oper(sq.getOper()), 
	position(pos),direction(dir),lhsSQStatus(VARIABLE_SQ),oriSQ(sq)
{
	if (isDerivWork()) {
		RRBuild generalRR(oper,oriSQ.get1stDerivPos(),oriSQ.get2edDerivPos(),
				oriSQ.get1stDerivDir(),oriSQ.get2edDerivDir());
		setupExpression(unsolvedIntegralList,generalRR);
	}else if (isRRWork()) {	
		RRBuild generalRR(rrType,oper,position);
		setupExpression(unsolvedIntegralList,generalRR);
	}else {
		RRBuild generalRR(oper);
		setupExpression(unsolvedIntegralList,generalRR);
	}
}

void RRSQ::setupExpression(const set<int>& unsolvedIntegralList, const RRBuild& generalRR)
{
	// form the shell quartet list in the RHS
	vector<ShellQuartet> tmpSQList;
	int nPossibleItems = generalRR.getNItems();
	tmpSQList.reserve(nPossibleItems);
	generalRR.buildRRSQ(oriSQ,tmpSQList);
	
	// we drop these NULL shell quartets
	// but we need to remember the non-null sq position
	// so to use it in integral generation
	sqlist.reserve(nPossibleItems);
	sqPosList.reserve(nPossibleItems);
	int count = 0; // this corresponds to the index of sqlist
	for(int iSQ=0; iSQ<nPossibleItems; iSQ++) {
		const ShellQuartet& sq = tmpSQList[iSQ];
		if (sq.isnull()) {
			sqPosList.push_back(-1);
			continue;
		}else{
			sqPosList.push_back(count);
			sqlist.push_back(sq);
			count++;
		}
	}

	// the default rhs shell quartet status always local variable
	rhsSQStatus.assign(getNItems(),VARIABLE_SQ); 

	// initialize the RHS and coefficients
	RHS.reserve(getNItems());
	coe.reserve(getNItems());
	for(int item=0; item<getNItems(); item++) {
		list<int> tmpRHS;
		list<string> tmpCoe;
		RHS.push_back(tmpRHS);
		coe.push_back(tmpCoe);
	}

	// form each integral according to the RR expansion
	for(set<int>::const_iterator it = unsolvedIntegralList.begin(); 
			it != unsolvedIntegralList.end(); ++it) {

		// set up tmp result vector
		vector<string> c;
		vector<int> rhs;
		c.reserve(nPossibleItems);
		rhs.reserve(nPossibleItems);

		// create new integrals and form its RR expansion
		int index = *it; 
		LHS.push_back(index);
		Integral I(oriSQ,index);
		if (generalRR.isNonRRNonDerivWork() && oper == THREEBODYKI) {
			generalRR.build3BodyKIInt(I,direction,c,rhs); 
		}else{
			generalRR.buildRRInt(I,c,rhs);
		}

		// now push it into the result
		for(int item=0; item<nPossibleItems; item++) {

			// is it this term corresponding to null sq?
			if (sqPosList[item]<0) continue;

			// get the coefficients and rhs integral index for the given item
			int rhsIntIndex = rhs[item];
			const string& coefficients = c[item];

			// push them back into the result for this term
			int term = sqPosList[item];
			list<int>& rhsArray = RHS[term];
			list<string>& coeArray = coe[term];
			rhsArray.push_back(rhsIntIndex);
			coeArray.push_back(coefficients);
		}
	}
}

bool RRSQ::operator<(const RRSQ& rrsq) const 
{
	// for HRR, we need to determine it's side
	if (isRRWork() && rrType == HRR) {

		// determine the side
		int thisSide = BRA;
		if (position == KET1 || position == KET2) thisSide = KET;
		int otherSide = BRA;
		if (rrsq.position == KET1 || rrsq.position == KET2) otherSide = KET;
		crash(thisSide != otherSide, "operator < in RRSQ should have same side on HRR");

		// now do HRR comparing
		return oriSQ.hrrLessThan(rrsq.oriSQ,position);
	}

	// consider the non-RR non-Deriv case
	// whether the direction information is defined
	if (oper == THREEBODYKI && isNonRRNonDerivWork()) {
		if (direction == NO_DERIV || rrsq.direction == NO_DERIV) {
			crash(true,"illegal direction on non-RR non-derivative expression of RRSQ in operator <");
		}
		return direction > rrsq.direction;
	}

	// now the only possible case are deriv job
	// or VRR
	// they can be done by using the sq's operator <
	return (oriSQ < rrsq.oriSQ);
}

void RRSQ::getUnsolvedIntList(const int& sqindex, set<int>& unsolvedList) const
{
	// we note that it's only real integral index we care about
	const list<int>& rhs = RHS[sqindex];
	for(list<int>::const_iterator it=rhs.begin(); it!=rhs.end(); ++it) {
		int value = *it;
		if (value >=0) {
			unsolvedList.insert(value);
		}
	}
}

void RRSQ::getUnsolvedIntList(const int& sqindex, const set<int>& lhsIndex, 
		set<int>& unsolvedList) const
{
	// we note that it's only real integral index we care about
	const list<int>& rhs = RHS[sqindex];
	for(set<int>::const_iterator it=lhsIndex.begin(); it!=lhsIndex.end(); ++it) {
		int value = *it;

		// if the LHS is meaningless
		// we will have an error in pos = distance(itBegin,it2)
		// to get the correct position for rhs
		if (value<0) {
			crash(true, "the LHS integral index is meaningless in RRSQ::getUnsolvedIntList");
		};

		// get the position of this LHS integral 
		// we note, that this lhs must be appearing on LHS
		list<int>::const_iterator itBegin = LHS.begin();
		list<int>::const_iterator itEnd   = LHS.end();
		list<int>::const_iterator it2 = find(itBegin,itEnd,value);
		if (it2 == itEnd) {
			Integral I(oriSQ,value);
			cout << I.getName() << endl;
			crash(true, "why we do not have this LHS integral in getUnsolvedIntList?");
		}
		int pos = distance(itBegin,it2);

		// now let's add in the corresponding RHS integral
		list<int>::const_iterator it3 = rhs.begin();
		advance(it3,pos);
		int val = *it3;
		if (val>=0) {
			unsolvedList.insert(val);
		}
	}
}

void RRSQ::formLHSIndexSet(set<int>& intList) const
{
	// we note that it's only real integral index we care about
	for(list<int>::const_iterator it=LHS.begin(); it!=LHS.end(); ++it) {
		int value = *it;
		if (value >=0) {
			intList.insert(value);
		}
	}
}

void RRSQ::updateLHS(const set<int>& unsolvedIntList) 
{
	// do we really need to update it?
	// if this is full integrals, then we do not need to do it
	bool needUpdating = true;
	int totalIntNum = oriSQ.getNInts();
	if (totalIntNum == (int)LHS.size()) needUpdating = false;
	if (! needUpdating) return;

	// now we build rrbuild according to the information
	// set up in the constructor
	if (isDerivWork()) {
		RRBuild generalRR(oper,oriSQ.get1stDerivPos(),oriSQ.get2edDerivPos(),
				oriSQ.get1stDerivDir(),oriSQ.get2edDerivDir());
		updateLHSExpression(unsolvedIntList,generalRR); 
	}else if (isRRWork()) {
		RRBuild generalRR(rrType,oper,position);
		updateLHSExpression(unsolvedIntList,generalRR); 
	}else {
		if (oper == THREEBODYKI) {
			RRBuild generalRR(oper);
			updateLHSExpression(unsolvedIntList,generalRR); 
		}else{
			crash(true,"illegal situation in updateLHS, we do not know how to build RRBuild?");
		}
	}
}

void RRSQ::updateLHSExpression(const set<int>& unsolvedIntList, const RRBuild& generalRR) 
{
	// now do some prepare work in advance
	int nPossibleItems = generalRR.getNItems();

	// now let's go to compare this rrsq's LHS and the input new
	// unsolvedIntList. If they are totally same, we do not need 
	// to do anything. else we just append the new item to a 
	// proper position
	for(set<int>::const_iterator it = unsolvedIntList.begin(); 
			it!= unsolvedIntList.end(); ++it) {

		// firstly, get the integral index
		int index = *it;
		list<int>::iterator it2 = find(LHS.begin(),LHS.end(),index);
		if (it2 != LHS.end()) continue;

		// now we need to do real updating work
		// set up tmp result vector
		vector<string> c;
		vector<int> rhs;
		c.reserve(nPossibleItems);
		rhs.reserve(nPossibleItems);

		// create new integrals and form its RR expansion
		Integral I(oriSQ,index);
		if (generalRR.isNonRRNonDerivWork() && oper == THREEBODYKI) {
			generalRR.build3BodyKIInt(I,direction,c,rhs); 
		}else{
			generalRR.buildRRInt(I,c,rhs);
		}

		// now let's find the position to insert the new integrals
		// this is used to insert elements for RHS
		// if the intIndex is larger than all of old elements in LHS
		// then we just use -1, indicating we need to push back intIndex
		// else we will use insert
		//
		// we do this is because the advance function for iterator 
		// can not go beyond the last element. Therefore, if the given
		// intIndex should be push back, then we just do push back
		//
		int pos = -1;
		if (index>LHS.back()) {
			LHS.push_back(index);
		}else{
			int inc = 0;
			list<int>::iterator it3;
			for(it3=LHS.begin(); it3!=LHS.end(); ++it3) {
				pos = inc;
				int val = *it3;
				if (val > index) {
					break;
				}else{
					inc++;
				}
			}
			LHS.insert(it3,index);
		}

		// now push it into the result
		for(int item=0; item<nPossibleItems; item++) {

			// is it this term corresponding to null sq?
			if (sqPosList[item]<0) continue;

			// get the coefficients and rhs integral index for the given item
			int rhsIntIndex = rhs[item];
			const string& coefficients = c[item];

			// here we need to know, that the LHS data sequence should be 
			// maintained
			int term = sqPosList[item];
			list<int>& rhsArray = RHS[term];
			if (pos<0) {
				rhsArray.push_back(rhsIntIndex);
			}else{
				list<int>::iterator it1 = rhsArray.begin();
				advance(it1,pos);
				rhsArray.insert(it1,rhsIntIndex);
			}

			// update the coefficients
			list<string>& coeArray = coe[term];
			if (pos<0) {
				coeArray.push_back(coefficients);
			}else{
				list<string>::iterator it1 = coeArray.begin();
				advance(it1,pos);
				coeArray.insert(it1,coefficients);
			}
		}
	}
}

int RRSQ::countRHSIntegralNum() const 
{
	int count = 0;
	for(int item=0; item<getNItems(); item++) {
		const list<int>& rhs = RHS[item];
		for(list<int>::const_iterator it=rhs.begin(); it!=rhs.end(); ++it) {
			int value = *it;
			if (value>=0) count++;
		}
	}
	return count;
}

bool RRSQ::checkCompleteness(const ShellQuartet& sq, const set<int>& sqLHS, 
		set<int>& missingLHS) const
{
	// firstly, let's check whether the lhs of rrsq appears in this rhs
	// we note, that each sq could only appear in RHS once
	// since in RR all of terms are different
	int rhsItem = -1;
	for(int item=0; item<getNItems(); item++) {
		const ShellQuartet& rhsSQ = getRHSSQ(item);
		if (rhsSQ == sq) {
			rhsItem = item; 
			break;
		}
	}

	// shall we return now?
	if (rhsItem<0) return true;

	// create rhs set
	// we note that it's only real integral index we care about
	const list<int>& rhs = RHS[rhsItem];
	bool pass = true;
	for(list<int>::const_iterator it=rhs.begin(); it!=rhs.end(); ++it) {
		int value = *it;
		if (value>=0) {

			// now we try to find the same term in sq lhs
			set<int>::const_iterator it2 = sqLHS.find(value);

			// now here is the trouble...
			if (it2 == sqLHS.end()) {
				missingLHS.insert(value);
				pass = false;
			}
		}
	}
	return pass;
}

void RRSQ::rhsArrayIndexTransform(const RRSQ& rrsq)
{
	// let's check whether the lhs of rrsq appears in this rhs
	// we note, that each sq could only appear in RHS once
	// since in RR all of terms are different
	int rhsItem = -1;
	const ShellQuartet& sq = rrsq.getLHSSQ();
	for(int item=0; item<getNItems(); item++) {
		const ShellQuartet& rhsSQ = getRHSSQ(item);
		if (rhsSQ == sq) {
			rhsItem = item; 
			break;
		}
	}

	// shall we return now?
	if (rhsItem<0) return;

	// now let's refresh the rhs status
	int status = rrsq.getLHSSQStatus();
	if (status == GLOBAL_RESULT_SQ) {
		crash(true, "the global result should not appear on the RHS formula, wrong in RRSQ::rhsArrayIndexTransform");
	}
	rhsSQStatus[rhsItem] = status;

	// whether the sq status is not in array?
	// if so we do not do it
	if (! inArrayStatus(status)) {
		return;
	}

	// now replace the integral index with its position 
	list<int>& rhs = RHS[rhsItem];
	const list<int>& lhs = rrsq.getLHSIndexArray();
	for(list<int>::const_iterator it=lhs.begin(); it!=lhs.end(); ++it) {
		int val = *it;
		if (val>=0) {
			int pos = distance(lhs.begin(),it);
			replace(rhs.begin(),rhs.end(),val,pos);
		}
	}
}

void RRSQ::rhsArrayIndexTransform(const vector<ShellQuartet>& bottomSQList,
		const vector<int>& bottomSQStatus, const vector<set<int> >& unsolvedIntList)
{
	// let's check whether the RHS of this rrsq has any bottom integrals
	// and this RHS must be in array form
	// if not, we just simply bypass
	bool hasBottomSQ = false;
	vector<int> bottomSQPosList(getNItems(),-1);
	for(int item=0; item<getNItems(); item++) {
		const ShellQuartet& rhsSQ = getRHSSQ(item);
		for(int iSQ=0; iSQ<(int)bottomSQList.size(); iSQ++) {
			const ShellQuartet& sq = bottomSQList[iSQ];
			int status = bottomSQStatus[iSQ];
			if (rhsSQ == sq && inArrayStatus(status)) {
				hasBottomSQ = true;
				bottomSQPosList[item] = 1;
				break;
			}
		}
	}

	// shall we return now?
	if (! hasBottomSQ) return;

	// now replace the integral index with its position 
	for(int rhsItem=0; rhsItem<getNItems(); rhsItem++) {

		// let's see whether this term corresponding to bottom sq?
		if (bottomSQPosList[rhsItem] < 0) continue;

		// find out which position it corresponding to bottomSQ?
		int pos = -1;
		const ShellQuartet& rhsSQ = getRHSSQ(rhsItem);
		for(int iSQ=0; iSQ<(int)bottomSQList.size(); iSQ++) {
			const ShellQuartet& sq = bottomSQList[iSQ];
			if (rhsSQ == sq) {
				pos = iSQ;
				break;
			}
		}

		// double check whether we find it
		if (pos == -1) {
			cout << "rhsSQ " << rhsSQ.getName() << endl;
			cout << "lhsSQ " << oriSQ.getName() << endl;
			crash(true,"in rhsArrayIndexTransform with bottom sq list we did not find the proper bottom sq");
		}

		// update the RHS status
		int status = bottomSQStatus[pos];
		rhsSQStatus[rhsItem] = status;

		// now let's work
		list<int>& rhs = RHS[rhsItem];
		const set<int>& lhs = unsolvedIntList[pos];
		for(set<int>::const_iterator it=lhs.begin(); it!=lhs.end(); ++it) {
			int val = *it;
			if (val>=0) {
				int pos = distance(lhs.begin(),it);
				replace(rhs.begin(),rhs.end(),val,pos);
			}
		}
	}
}

void RRSQ::lhsArrayIndexTransform()
{
	int pos=0;
	for(list<int>::iterator it=LHS.begin(); it!=LHS.end(); ++it) {
		*it = pos;
		pos++;
	}
}

void RRSQ::print(const int& nSpace, const SQIntsInfor& infor, ofstream& file) const 
{
	// obtain the information for this rrsq
	string lhsName = oriSQ.getName();
	int nInts = LHS.size();
	int nTotalInts = oriSQ.getNInts();
	int diff = nTotalInts - nInts;

	// the position name
	string positionName = "NULL";
	if (position == BRA1) {
		positionName = "BRA1";
	}else if (position == BRA2) {
		positionName = "BRA2";
	}else if (position == KET1) {
		positionName = "KET1";
	}else if (position == KET2) {
		positionName = "KET2";
	}

	// the section name
	string sectionName = "VRR";
	if (isRRWork() && rrType == HRR) {
		sectionName = "HRR";
	}
	if (isDerivWork()) {
		sectionName = "DERIV";
	}
	if(isNonRRNonDerivWork()) {
		sectionName = "NON-RR";
	}

	// now print comment section to file
	file << endl;
	string line;
	line = "/************************************************************";
	printLine(nSpace,line,file);
	line = " * shell quartet name: " + lhsName;
	printLine(nSpace,line,file);
	if (isRRWork()) {
		line = " * expanding position: " + positionName;
		printLine(nSpace,line,file);
	}
	line = " * code section is: " + sectionName;
	printLine(nSpace,line,file);
	line = " * totally " + lexical_cast<string>(diff) + " integrals are omitted ";
	printLine(nSpace,line,file);
	if (direction == DERIV_X) {
		line = " * direction is on x";
		printLine(nSpace,line,file);
	}else if (direction == DERIV_Y) {
		line = " * direction is on y";
		printLine(nSpace,line,file);
	}else if (direction == DERIV_Z) {
		line = " * direction is on z";
		printLine(nSpace,line,file);
	}
	for(int item=0; item<getNItems(); item++) {
		const ShellQuartet& sq = sqlist[item];
		line = " * RHS shell quartet name: " + sq.getName();
		printLine(nSpace,line,file);
	}
	line = " ************************************************************/";
	printLine(nSpace,line,file);

	// now set up the declare of the vector 
	//
	// we do it when the rrsq is in array form, but it's not function in/out
	// because function in/out declare will be handled in the main file
	// so, it's only when the lhs sq status in array_sq, we do declare
	// this and all of codes are not in file split form.
	// 
	// take care of three body ki case, too
	//
	bool doDeclare = true;
	if (lhsSQStatus != ARRAY_SQ) doDeclare = false;
	if (isNonRRNonDerivWork() && oper == THREEBODYKI) {
		if (direction != DERIV_X) doDeclare = false;
	}
	if (doDeclare) {
		string arrayType = infor.getArrayType();
		string arrayName = oriSQ.formArrayName(rrType);
		string nLHSInts  = lexical_cast<string>(nInts);
		string declare   = infor.getArrayDeclare(nLHSInts);
		line = arrayType + arrayName + declare;
		printLine(nSpace,line,file);
	}

	// check that whether we need additional offset for the final result
	// here the nInts we use the total number of integrals
	// for the result shell quartes
	// for example, in ESP per grid points it includes all of integal results
	bool withAdditionalOffset = resultIntegralHasAdditionalOffset(oper);
	string additionalOffset;
	if (lhsSQStatus == GLOBAL_RESULT_SQ && withAdditionalOffset) {
		int nTolInts = infor.nInts();
		additionalOffset = determineAdditionalOffset(oper,nTolInts);
	}

	// we may need the += rather than = for the LHS
	bool usePlus =false;
	if (isNonRRNonDerivWork() && oper == THREEBODYKI) {
		if (direction != DERIV_X) usePlus = true;
	}

	// here we constructor vectors for printing
	// firstly it's the right hand side
	vector<vector<int> > rhs;
	rhs.reserve(getNItems());
	for(int i=0; i<getNItems(); i++) {

		// reserve space
		vector<int> tmp;
		tmp.reserve(LHS.size());

		// push into data
		const list<int>& rhsItems = RHS[i];
		for(list<int>::const_iterator it=rhsItems.begin(); it!=rhsItems.end(); ++it) {
			tmp.push_back(*it);
		}

		// now finish this term
		rhs.push_back(tmp);
	}

	// now it's the coefficients
	vector<vector<string> > coef;
	coef.reserve(getNItems());
	for(int i=0; i<getNItems(); i++) {

		// reserve space
		vector<string> tmp;
		tmp.reserve(LHS.size());

		// push into data
		const list<string>& coeItems = coe[i];
		for(list<string>::const_iterator it=coeItems.begin(); it!=coeItems.end(); ++it) {
			tmp.push_back(*it);
		}

		// now finish this term
		coef.push_back(tmp);
	}

	// print each integrals
	for(list<int>::const_iterator it=LHS.begin(); it!=LHS.end(); ++it) {

		// set up the expression
		int index = *it;
		string expression;
		expression.reserve(200);

		// firstly, it's the left hand side
		// if this is array or it's final result, we go here
		if (lhsSQStatus == GLOBAL_RESULT_SQ || inArrayStatus(lhsSQStatus)) {

			// compute the offset(array index) for the given integral
			//
			// we note, that if the shell quartet is the final result;
			// then it must contains all of integrals; therefore there's 
			// no neglecting integrals for this shell quartet
			//
			// because of this, we can use the index to re-create the 
			// integral (as you can see in the getOffset), even though
			// that this rrsq has been transformed into array index
			//
			// this explains why it's valid to pass index into 
			// the getOffset function
			int offset = index;
			if (lhsSQStatus == GLOBAL_RESULT_SQ) {
				offset = infor.getOffset(oriSQ,index);
			}

			// name
			string arrayName = oriSQ.formArrayName(rrType);
			if (lhsSQStatus == GLOBAL_RESULT_SQ) {
				arrayName = "abcd";
			}

			// now form the LHS
			if (withAdditionalOffset && lhsSQStatus == GLOBAL_RESULT_SQ) {
				arrayName = arrayName + "[" + additionalOffset + "+" + lexical_cast<string>(offset) + "]";
			}else{
				arrayName = arrayName + "[" + lexical_cast<string>(offset) + "]";
			}

			// consider the cases that we may need to use "+="
			if (usePlus) {
				expression = arrayName + " += ";
			}else{
				expression = arrayName + " = ";
			}
		}else{
			// now form the var type LHS
			// get the variable name first
			Integral I(oriSQ,index);
			string varName = I.formVarName(rrType);
			expression = "Double " + varName + " = ";
		}

		// determine the RHS expression
		int pos = distance(LHS.begin(),it);
		for(int item=0; item<getNItems(); item++) {

			// obtain the corresponding terms
			const vector<int>& rhsArray  = rhs[item];
			const vector<string>& c      = coef[item];

			// expression
			int rhsIndex = rhsArray[pos];
			const string& coefficients = c[pos];
			const ShellQuartet& sq = getRHSSQ(item);

			// whether this is null integral?
			if (rhsIndex==NULL_POS) continue;

			// considering the coefficients
			// we drop the multiplier of 1
			string k;
			if (coefficients[0] == '+' && coefficients[1] == '1' && coefficients[2] == '*') {
				unsigned pos = 3;
				k = coefficients.substr(pos);
				k = "+" + k;
			}else if (coefficients[0] == '-' && coefficients[1] == '1' && coefficients[2] == '*') {
				unsigned pos = 3;
				k = coefficients.substr(pos);
				k = "-" + k;
			}else{
				k = coefficients;
			}

			// determine the array status
			int rhsStatus = rhsSQStatus[item];
			bool rhsWithArray = false;
			if (inArrayStatus(rhsStatus)) {
				rhsWithArray = true;
			}

			// now generate the RHS
			if (rhsWithArray && ! sq.isSTypeSQ()) {

				// sq name
				string sqName = sq.formArrayName(rrType);
				if (k == "1") {
					expression += sqName + "[" + lexical_cast<string>(rhsIndex) + "]";
				}else{
					expression += k + "*";
					expression += sqName + "[" + lexical_cast<string>(rhsIndex) + "]";
				}

			}else{

				// now integral name
				Integral I(sqlist[item],rhsIndex);
				string intName = I.formVarName(rrType);
				if (k == "1") {
					expression += intName;
				}else{
					expression += k + "*";
					expression += intName;
				}
			}
		}

		// finally add simicolon
		expression += ";";
		//cout << expression << endl;

		// now print it to file
		printLine(nSpace,expression,file);
	}
}

void RRSQ::printArrayToVar(const int& nSpace, ofstream& file) const 
{
	// now print comment section to file
	file << endl;
	string line;
	line = "/************************************************************";
	printLine(nSpace,line,file);
	line = " * convert from array into variable form: " + oriSQ.getName();
	printLine(nSpace,line,file);
	line = " ************************************************************/";
	printLine(nSpace,line,file);

	// for using this function, the LHS must not be in array status
	if (inArrayStatus(lhsSQStatus)) {
		crash(true, "something wrong for RRSQ::printArrayToVar, LHS already been changed into array status");
	}
	if (lhsSQStatus == GLOBAL_RESULT_SQ) {
		crash(true, "something wrong for RRSQ::printArrayToVar, LHS is final result of abcd");
	}

	// print each integrals
	int offset = 0;
	string arrayName = oriSQ.formArrayName(rrType);
	for(list<int>::const_iterator it=LHS.begin(); it!=LHS.end(); ++it) {
		string rhs = arrayName + "[" + lexical_cast<string>(offset) + "]";
		int index = *it;
		Integral I(oriSQ,index);
		string varName = I.formVarName(rrType);
		string expression = "Double " + varName + " = " + rhs + ";";
		printLine(nSpace,expression,file);
		offset++;
	}
}

void RRSQ::printVarToArray(const int& nSpace, ofstream& file) const 
{
	// now print comment section to file
	file << endl;
	string line;
	line = "/************************************************************";
	printLine(nSpace,line,file);
	line = " * convert from variable into array form: " + oriSQ.getName();
	printLine(nSpace,line,file);
	line = " ************************************************************/";
	printLine(nSpace,line,file);

	// for using this function, the LHS must not be in array status
	if (inArrayStatus(lhsSQStatus)) {
		crash(true, "something wrong for RRSQ::printVarToArray, LHS already been changed into array status");
	}
	if (lhsSQStatus == GLOBAL_RESULT_SQ) {
		crash(true, "something wrong for RRSQ::printVarToArray, LHS is final result of abcd");
	}

	// print each integrals
	int offset = 0;
	string arrayName = oriSQ.formArrayName(rrType);
	for(list<int>::const_iterator it=LHS.begin(); it!=LHS.end(); ++it) {
		string lhs = arrayName + "[" + lexical_cast<string>(offset) + "]";
		int index = *it;
		Integral I(oriSQ,index);
		string varName = I.formVarName(rrType);
		string expression = lhs + " = " + varName + ";";
		printLine(nSpace,expression,file);
		offset++;
	}
}

void RRSQ::debug_print() const 
{
	// here we constructor vectors for printing
	// firstly it's the right hand side
	vector<vector<int> > rhs;
	rhs.reserve(getNItems());
	for(int i=0; i<getNItems(); i++) {

		// reserve space
		vector<int> tmp;
		tmp.reserve(LHS.size());

		// push into data
		const list<int>& rhsItems = RHS[i];
		for(list<int>::const_iterator it=rhsItems.begin(); it!=rhsItems.end(); ++it) {
			tmp.push_back(*it);
		}

		// now finish this term
		rhs.push_back(tmp);
	}

	// now it's the coefficients
	vector<vector<string> > coef;
	coef.reserve(getNItems());
	for(int i=0; i<getNItems(); i++) {

		// reserve space
		vector<string> tmp;
		tmp.reserve(LHS.size());

		// push into data
		const list<string>& coeItems = coe[i];
		for(list<string>::const_iterator it=coeItems.begin(); it!=coeItems.end(); ++it) {
			tmp.push_back(*it);
		}

		// now finish this term
		coef.push_back(tmp);
	}

	// print each integrals
	for(list<int>::const_iterator it=LHS.begin(); it!=LHS.end(); ++it) {

		// set up the expression
		int index = *it;
		string expression;

		// now form the var type LHS
		Integral I(oriSQ,index);
		expression = "Double " + I.getName() + " = ";

		// determine the RHS expression
		int pos = distance(LHS.begin(),it);
		for(int item=0; item<getNItems(); item++) {

			// obtain the corresponding terms
			const vector<int>& rhsArray  = rhs[item];
			const vector<string>& c      = coef[item];

			// expression
			int rhsIndex = rhsArray[pos];
			const string& coefficients = c[pos];

			// whether this is null integral?
			if (rhsIndex==NULL_POS) continue;

			// considering the coefficients
			// we drop the multiplier of 1
			string k;
			if (coefficients[0] == '+' && coefficients[1] == '1' && coefficients[2] == '*') {
				unsigned pos = 3;
				k = coefficients.substr(pos);
				k = "+" + k;
			}else if (coefficients[0] == '-' && coefficients[1] == '1' && coefficients[2] == '*') {
				unsigned pos = 3;
				k = coefficients.substr(pos);
				k = "-" + k;
			}else{
				k = coefficients;
			}

			// now create RHS integral term
			// we note, that RHS term in general are actually
			// all tmp results.
			Integral I(sqlist[item],rhsIndex);
			if (k == "1") {
				expression += I.getName();
			}else{
				expression += k + "*";
				expression += I.getName();
			}
		}

		// finally add simicolon
		expression += ";";

		// now print it to file
		cout << expression << endl;
	}
}

