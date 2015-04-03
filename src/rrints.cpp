//
// CPPINTS: A C++ Program to Generate Analytical Integrals Based on Gaussian
// Form Primitive Functions
//
// Copyright (C) 2012-2015 Fenglai Liu
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
		const set<int>& unsolvedIntegralList, int jobOrder0): rrType(rrType0), position(pos), 
	jobOrder(jobOrder0),isIntegralIndex(true),oriSQ(sq)
{
	// form the general RR used in this RRSQ
	int oper = oriSQ.getOper();
	RRBuild generalRR(rrType,oper,position,jobOrder);

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
		Integral I(sq,index);
		generalRR.buildRRInt(I,c,rhs);

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
	// check rr type
#ifdef DEBUG
	crash(rrType != rrsq.rrType, "operator < in RRSQ should be on same rr type");
#endif

	// for HRR, we need to determine it's side
	if (rrType == HRR) {

		// determine the side
		int thisSide = BRA;
		if (position == KET1 || position == KET2) thisSide = KET;
		int otherSide = BRA;
		if (rrsq.position == KET1 || rrsq.position == KET2) otherSide = KET;
		crash(thisSide != otherSide, "operator < in RRSQ should have same side on HRR");

		// now do HRR comparing
		return oriSQ.hrrLessThan(rrsq.oriSQ,thisSide);

	}else{

		// let's go to see the position
		// if position is simply the direvative information
		// then we simply compare the position infor
		if (isDerivInfor(position) && isDerivInfor(rrsq.position)) {
			return derivCompare(position,rrsq.position);
		}

		// now compare the shell quartet
		return (oriSQ < rrsq.oriSQ);
	}
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

		// omit these meaningless integrals
		if (value<0) continue;

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

	// now do some prepare work in advance
	int oper = oriSQ.getOper();
	RRBuild generalRR(rrType,oper,position,jobOrder);
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
		generalRR.buildRRInt(I,c,rhs);

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
	// firstly, we need to change the label of whether it's integral index
	isIntegralIndex = false;

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

void RRSQ::lhsArrayIndexTransform()
{
	int pos=0;
	for(list<int>::iterator it=LHS.begin(); it!=LHS.end(); ++it) {
		*it = pos;
		pos++;
	}
}

void RRSQ::print(const int& nSpace, const int& status, const SQIntsInfor& infor, 
		const vector<ShellQuartet>& moduleResultList, ofstream& file) const 
{
	// check whether it's integral index 
	// for VRR, we note that since the contraction part will
	// generate the final VRR result, and in the contraction
	// part we will need to add in VRR results into the HRR input
	// however, here in the RR section both LHS and RHS shell quartets
	// therefore will become internal to VRR
	// for HRR, LHS and RHS are all internal for HRR
	// Thus, we only need to determine the array form depending on the 
	// rrtype
	bool withArray = infor.withArrayIndex(rrType);
	if (isIntegralIndex && withArray) {
		crash(true, "print in RRSQ confilicts with input infor print choice");
	}

	// obtain the information for this rrsq
	string lhsName = oriSQ.getName();
	int nInts = LHS.size();
	int nTotalInts = oriSQ.getNInts();
	int diff = nTotalInts - nInts;

	// now print comment section to file
	file << endl;
	string line;
	line = "/************************************************************";
	printLine(nSpace,line,file);
	if (status != TMP_RESULT) {
		line = " * result shell quartet name: " + lhsName;
		printLine(nSpace,line,file);
	}else{
		line = " * shell quartet name: " + lhsName;
		printLine(nSpace,line,file);
	}
	line = " * totally " + lexical_cast<string>(diff) + " integrals are omitted ";
	printLine(nSpace,line,file);
	line = " ************************************************************/";
	printLine(nSpace,line,file);

	// now set up the declare of the vector
	// for the LHS shell quartet
	if(withArray && status != FINAL_RESULT) {
		string arrayType = infor.getArrayType();
		string arrayName = oriSQ.formArrayName(rrType,status,position,infor.isComSQ());
		string nLHSInts  = lexical_cast<string>(nInts);
		string declare   = infor.getArrayDeclare(nLHSInts);
		line = arrayType + arrayName + declare;
		printLine(nSpace,line,file);
	}

	// check that whether we need additional offset for the final result
	int oper = infor.getOper();
	bool withAdditionalOffset = resultIntegralHasAdditionalOffset(oper);
	string additionalOffset;
	if (status == FINAL_RESULT && withAdditionalOffset) {
		// here the nInts we use the total number of integrals
		// for the result shell quartes, which is final results
		int nTolInts = infor.nInts();
		additionalOffset = determineAdditionalOffset(oper,nTolInts);
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

		// firstly, it's the left hand side
		if (status == FINAL_RESULT) {

			// compute the offset for the given integral
			int offset = infor.getOffset(oriSQ,index);

			// we note, that any RRSQ in VRR is not final result
			// therefore we do not need to consider the case of 
			// "+="
			// they will be considered in VRR contraction part
			string arrayName = "abcd";
			if (withAdditionalOffset) {
				arrayName = arrayName + "[" + additionalOffset + "+" + lexical_cast<string>(offset) + "]";
			}else{
				arrayName = arrayName + "[" + lexical_cast<string>(offset) + "]";
			}
			expression = arrayName + " = ";
		}else{
			if (withArray) {
				string arrayName = oriSQ.formArrayName(rrType,status,position,infor.isComSQ());
				expression = arrayName + "[" + lexical_cast<string>(index) + "] = ";
			}else{

				// now form the var type LHS
				// get the variable name first
				Integral I(oriSQ,index);
				string varName = I.formVarName(rrType,status,position,infor.isComSQ());
				expression = "Double " + varName + " = ";
			}
		}
		//expression = oriSQ.getName() + "[" + lexical_cast<string>(index) + "]" + " = ";

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
			// we note, that RHS integral could be tmp results
			// or module results 
			// here it can not be the final results, since the RHS
			// must be defined above in the content
			// also keep on eye on all S type of integrals
			// this should be expressed only with integral
			//
			// on the other hand, we note that if the RHS integral
			// is final results, so how to deal with it?
			//
			// For VRR, this could happens. However, since we make 
			// everything related to the VRR result to be processed
			// in contraction step - then everything computed during 
			// the VRR process will be only TMP or MODULE type results.
			// Therefore, for VRR this is not a problem.
			//
			// for HRR, this does not happen. For each "pure" shell quartet,
			// it's HRR process only has one final result, the shell quartet
			// itself. For shell quartets with modifiers (such as multiple
			// exponents in derivatives calculation, or composite shell case),
			// each sub-shell quartet will has its own HRR path. Different 
			// sub-shell quartet HRR pathes are independent with each other.
			//
			const ShellQuartet& sq = sqlist[item];

			// whether the given sq is a module result?
			int rhsStatus = TMP_RESULT;
			for(int iSQ=0; iSQ<(int)moduleResultList.size(); iSQ++) {
				const ShellQuartet& resultSQ = moduleResultList[iSQ];
				if (resultSQ == sq) {
					rhsStatus = MODULE_RESULT;
					break;
				}
			}

			// for the rhs sq, we can not know it's position information
			// set it to be null
			int rhsPos = NULL_POS;

			// now generate the RHS
			// we note, that this RHS term will never to be the final results
			// therefore we will check the module result as well as temp result
			if (withArray && ! sq.isSTypeSQ()) {

				// sq name
				string sqName = sq.formArrayName(rrType,rhsStatus,rhsPos,infor.isComSQ());
				if (k == "1") {
					expression += sqName + "[" + lexical_cast<string>(rhsIndex) + "]";
				}else{
					expression += k + "*";
					expression += sqName + "[" + lexical_cast<string>(rhsIndex) + "]";
				}

			}else{

				// now integral name
				Integral I(sqlist[item],rhsIndex);
				string intName = I.formVarName(rrType,rhsStatus,rhsPos,infor.isComSQ());
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

		// now print it to file
		printLine(nSpace,expression,file);
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

