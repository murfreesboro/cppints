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
#include "printing.h"
#include "derivinfor.h"
#include "inttype.h"
using namespace inttype;
using namespace printing;
using namespace derivinfor;

void FirstOrderDerivInfor::print(ofstream& file) const
{
	// get the position information
	string pos = "BRA1";
	if (derivPos == BRA2) {
		pos = "BRA2";
	}else if (derivPos == KET1) {
		pos = "KET1";
	}else if (derivPos == KET2) {
		pos = "KET2";
	}

	// now print it
	string line = "// " + pos;
	printLine(0,line,file);

	// now it's directions
	for(int i=0; i<getDerivDirLen(); i++) {
		int dir = getDerivDirection(i);

		// transform them into string format
		string d = "X";
		if(dir == DERIV_Y) {
			d = "Y";
		}else if (dir == DERIV_Z) {
			d = "Z";
		}
		string line = "// " + d;
		printLine(0,line,file);
	}
}

void SecondOrderDerivInfor::print(ofstream& file) const
{
	// get the position information
	string pos1 = "BRA1";
	string pos2 = "BRA1";
	if (firstDerivPos == BRA2) {
		pos1 = "BRA2";
	}else if (firstDerivPos == KET1) {
		pos1 = "KET1";
	}else if (firstDerivPos == KET2) {
		pos1 = "KET2";
	}
	if (secondDerivPos == BRA2) {
		pos2 = "BRA2";
	}else if (secondDerivPos == KET1) {
		pos2 = "KET1";
	}else if (secondDerivPos == KET2) {
		pos2 = "KET2";
	}

	// now print it
	string line = "// " + pos1 + "  " + pos2;
	printLine(0,line,file);

	// now it's directions
	for(int i=0; i<getDerivDirLen(); i++) {
		int dir1 = -1; 
		int dir2 = -1;
		getDerivDirection(i,dir1,dir2);

		// transform them into string format
		string d1 = "X";
		string d2 = "X";
		if(dir1 == DERIV_Y) {
			d1 = "Y";
		}else if (dir1 == DERIV_Z) {
			d1 = "Z";
		}
		if(dir2 == DERIV_Y) {
			d2 = "Y";
		}else if (dir2 == DERIV_Z) {
			d2 = "Z";
		}
		string line = "// " + d1 + "  " + d2;
		printLine(0,line,file);
	}
}

int DerivInfor::getLen1stDerivInforArray() const 
{
	int nBody = getOperOrder(oper);
	if (redundantPos == NULL_POS) return nBody;
	return nBody-1;
}

int DerivInfor::getLen2edDerivInforArray() const 
{
	int nBody = getOperOrder(oper);
	if (redundantPos == NULL_POS) {
		int n = nBody*(nBody+1)/2;
		return n;
	}
	int n = nBody - 1;
	int num = n*(n+1)/2;
	return num;
}

bool DerivInfor::checkInfor()
{

	// firstly, right now the highest job order
	// we support is only 2
	// we only defines the second deriv information array
	if (jobOrder<=0 || jobOrder>2) {
		crash(true, "Invalid jobOrder passed into the DerivInfor class, checkInfor failed");
		return false;
	}

	// now let's see the situation when redudant 
	// pos is needed
	// right now only two situation that redundant position is 
	// not needed:
	// 1  one body integral
	// 2  NAI
	// why NAI do not require redundant position please see
	// the manual, we have explanation over there
	int nBody = getOperOrder(oper);
	if (nBody != 1 && hasDerivRedundantPos(oper)) {
		if (redundantPos != BRA1 && redundantPos != BRA2 
				&& redundantPos != KET1 && redundantPos != KET2) {
			crash(true, "Invalid redundant position in DerivInfor class, checkInfor failed");
			return false;
		}
	}	

	// everything good
	return true;
}

DerivInfor::DerivInfor(const int& oper0, const int& jobOrder0):oper(oper0),
	jobOrder(jobOrder0),redundantPos(NULL_POS)
{
	// checking the input information
	if (! checkInfor()) {
		crash(true, "check infor failed in DerivInfor constructor(without redudantPos)");
	}

	// initilize
	if (jobOrder == 1) {
		int len = getLen1stDerivInforArray(); 
		firstDerivInfor.reserve(len);
	}else if (jobOrder == 2) {
		int len = getLen2edDerivInforArray(); 
		secondDerivInfor.reserve(len);
	}

	// whether this is one body integral?
	int nBody = getOperOrder(oper);
	if (nBody == 1) {
		int pos = BRA1;
		if (jobOrder == 1) {
			FirstOrderDerivInfor deriv(pos);
			firstDerivInfor.push_back(deriv);
		}else if (jobOrder == 2) {
			SecondOrderDerivInfor deriv(pos,pos);
			secondDerivInfor.push_back(deriv);
		}
	}else if (nBody == 2 && ! hasDerivRedundantPos(oper)) {
		int pos1 = BRA1;
		int pos2 = BRA2;
		if (jobOrder == 1) {
			FirstOrderDerivInfor deriv1(pos1);
			firstDerivInfor.push_back(deriv1);
			FirstOrderDerivInfor deriv2(pos2);
			firstDerivInfor.push_back(deriv2);
		}else if (jobOrder == 2) {
			SecondOrderDerivInfor deriv11(pos1,pos1);
			secondDerivInfor.push_back(deriv11);
			SecondOrderDerivInfor deriv12(pos1,pos2);
			secondDerivInfor.push_back(deriv12);
			SecondOrderDerivInfor deriv22(pos2,pos2);
			secondDerivInfor.push_back(deriv22);
		}
	}else{
		crash(true, "Invalid situation in DerivInfor class without redundant pos, constructor failed");
	}
}

DerivInfor::DerivInfor(const int& oper0, const int& jobOrder0,
		const int& redundantPos0):oper(oper0),
	jobOrder(jobOrder0),redundantPos(redundantPos0)
{
	// checking the input information
	if (! checkInfor()) {
		crash(true, "check infor failed in DerivInfor constructor(wit redudantPos)");
	}

	// now do initilization
	if (jobOrder == 1) {
		int len = getLen1stDerivInforArray(); 
		firstDerivInfor.reserve(len);
	}else if (jobOrder == 2) {
		int len = getLen2edDerivInforArray(); 
		secondDerivInfor.reserve(len);
	}

	// now let's stuff the deriv infor
	int nBody = getOperOrder(oper);
	if (jobOrder == 1) {

		// initilize the position vector
		vector<int> posList;
		posList.reserve(4);
		for(int i=0; i<nBody; i++) {

			// set the original pos
			int pos = BRA1;
			if (i==1) {
				pos = BRA2;
			}else if (i==2) {
				pos = KET1;
			}else if (i==3) {
				pos = KET2;
			}

			// is it the redundant pos?
			if (pos == redundantPos) continue;

			// now the position is available
			posList.push_back(pos);
		}

		// now let's form the content
		for(int i=0; i<(int)posList.size(); i++) {
			int pos = posList[i];
			FirstOrderDerivInfor deriv(pos);
			firstDerivInfor.push_back(deriv);
		}

	}else if (jobOrder == 2) {

		// form the pos list
		// the first element is for derivPos1
		// the second element is fro derivPos2
		// here we just list every possible situation
		int num = getLen2edDerivInforArray();
		vector<int> posList(2*num);
		if (nBody == 2) {
			if (redundantPos == BRA1) {

				// bb
				posList[0] = BRA2;
				posList[1] = BRA2;
			}else{

				// aa
				posList[0] = BRA1;
				posList[1] = BRA1;
			}
		}else if (nBody == 3) {
			if (redundantPos == BRA1) {

				// bb
				posList[0] = BRA2;
				posList[1] = BRA2;

				// bc
				posList[2] = BRA2;
				posList[3] = KET1;

				// cc
				posList[4] = KET1;
				posList[5] = KET1;
			}else if (redundantPos == BRA2) {

				// aa
				posList[0] = BRA1;
				posList[1] = BRA1;

				// ac
				posList[2] = BRA1;
				posList[3] = KET1;

				// cc
				posList[4] = KET1;
				posList[5] = KET1;
			}else{

				// aa
				posList[0] = BRA1;
				posList[1] = BRA1;

				// ab
				posList[2] = BRA1;
				posList[3] = BRA2;

				// bb
				posList[4] = BRA2;
				posList[5] = BRA2;

			}
		}else if (nBody == 4) {
			if (redundantPos == BRA1) {

				// bb
				posList[0] = BRA2;
				posList[1] = BRA2;

				// bc
				posList[2] = BRA2;
				posList[3] = KET1;

				// bd
				posList[4] = BRA2;
				posList[5] = KET2;

				// cc
				posList[6] = KET1;
				posList[7] = KET1;

				// cd
				posList[8] = KET1;
				posList[9] = KET2;

				// dd
				posList[10] = KET2;
				posList[11] = KET2;
			}else if (redundantPos == BRA2) {

				// aa
				posList[0] = BRA1;
				posList[1] = BRA1;

				// ac
				posList[2] = BRA1;
				posList[3] = KET1;

				// ad
				posList[4] = BRA1;
				posList[5] = KET2;

				// cc
				posList[6] = KET1;
				posList[7] = KET1;

				// cd
				posList[8] = KET1;
				posList[9] = KET2;

				// dd
				posList[10] = KET2;
				posList[11] = KET2;
			}else if (redundantPos == KET1) {

				// aa
				posList[0] = BRA1;
				posList[1] = BRA1;

				// ab
				posList[2] = BRA1;
				posList[3] = BRA2;

				// ad
				posList[4] = BRA1;
				posList[5] = KET2;

				// bb
				posList[6] = BRA2;
				posList[7] = BRA2;

				// bd
				posList[8] = BRA2;
				posList[9] = KET2;

				// dd
				posList[10] = KET2;
				posList[11] = KET2;

			}else{

				// aa
				posList[0] = BRA1;
				posList[1] = BRA1;

				// ab
				posList[2] = BRA1;
				posList[3] = BRA2;

				// ac
				posList[4] = BRA1;
				posList[5] = KET1;

				// bb
				posList[6] = BRA2;
				posList[7] = BRA2;

				// bc
				posList[8] = BRA2;
				posList[9] = KET1;

				// cc
				posList[10] = KET1;
				posList[11] = KET1;

			}
		}else {
			crash(true, "Invalid situation in DerivInfor class with redundant pos, constructor failed");
		}

		// now construct the content
		for(int i=0; i<num; i++) {
			int pos1 = posList[2*i];
			int pos2 = posList[2*i+1];
			SecondOrderDerivInfor deriv(pos1,pos2);
			secondDerivInfor.push_back(deriv);
		}
	}
}

int DerivInfor::getOffset(const int& derivPos, const int& derivDir) const {
	int n=0;
	int len = getLen1stDerivInforArray(); 
	for(int i=0; i<len; i++) {
		const FirstOrderDerivInfor& infor = get1stDerivInfor(i);

		// now we either get to the same deriv infor
		// else we just add up all of directions 
		if (derivPos == infor.getDerivPos()) { 
			int l = infor.getDerivDirLen();
			for(int j=0; j<l; j++) {
				if (derivDir == infor.getDerivDirection(j)) {
					return n;
				}
				n++;
			}
		}else{
			n = n + infor.getDerivDirLen();
		}
	}

	// we should not get to here
	// just release a crash
	crash(true,"fail to get offset for the input dirv pos/dir(1st order) in DerivInfor class");
	return -1;
}

int DerivInfor::getOffset(const int& derivPos1, const int& derivPos2, 
		const int& derivDir1, const int& derivDir2) const {
	int n=0;
	int len = getLen2edDerivInforArray(); 
	for(int i=0; i<len; i++) {
		const SecondOrderDerivInfor& infor = get2edDerivInfor(i);
		int pos1, pos2;
		infor.getDerivPos(pos1,pos2);
		if (derivPos1 == pos1 && derivPos2 == pos2) {
			int l = infor.getDerivDirLen();
			for(int j=0; j<l; j++) {
				int dir1, dir2;
				infor.getDerivDirection(j,dir1,dir2);
				if (derivDir1 == dir1 && derivDir2 == dir2) {
					return n;
				}
				n++;
			}
		}else{
			n = n + infor.getDerivDirLen();
		}
	}

	// we should not get to here
	// just release a crash
	crash(true,"fail to get offset for the input dirv pos/dir(2ed order) in DerivInfor class");
	return -1;
}


int DerivInfor::getTotalNumDeriv() const
{
	int n=0;
	if(jobOrder == 1) { 
		int len = getLen1stDerivInforArray(); 
		for(int i=0; i<len; i++) {
			const FirstOrderDerivInfor& infor = get1stDerivInfor(i);
			n = n + infor.getDerivDirLen();
		}
	}else if (jobOrder == 2) {
		int len = getLen2edDerivInforArray(); 
		for(int i=0; i<len; i++) {
			const SecondOrderDerivInfor& infor = get2edDerivInfor(i);
			n = n + infor.getDerivDirLen();
		}
	}else{
		crash(true, "incorrect job order in getTotalNumDeriv of DerivInfor class");
	}
	return n;
}

