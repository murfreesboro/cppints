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
#include "general.h"
#include "inttype.h"
#include "sqints.h"
#include "infor.h"
using namespace inttype;
using namespace infor;
using namespace sqints;

extern void codeGen(const Infor& infor, const int& oper, const int& derivOrder);

int main(int argc, char* argv[]) {

	// we should have infor file passed in
	string inforFile = argv[1];
	Infor infor(inforFile);

	// now do jobs
	const vector<int>& joblist = infor.getJobList();
	int nJobs = joblist.size();
	int order = infor.derivOrder;
	for(int iJob=0; iJob<nJobs; iJob++) {
		int job = joblist[iJob];
		codeGen(infor,job,order);
	}

	/*
		SQInts sqints(infor,2,1,KINETIC,0);
		sqints.codeGeneration();
	 */

	/*
		Basis test(0,2,0);
		cout << test.getName() << endl;
		int oper = ERI;
		Integral I1(test,test,test,test,oper,100);
		cout << I1.getName() << endl;

		Shell s(3);
		vector<Basis> b;
		s.getBasis(b);
		for(int i=0; i<(int)b.size(); i++)
		cout << b[i].getName() << endl;

		if (s.hasThisBasisSet(test)) {
		cout << "yes!!" << endl;
		}

		ShellQuartet sq1(2,2,2,2,ERI,0);
		cout << sq1.getName() << endl;
		ShellQuartet sq2(2,1,2,1,ERI,0);
		cout << sq2.getName() << endl;
		if (sq1<sq2) {
		cout << "less than!!!!" << endl;
		}else{
		cout << "larger than!!!" << endl;
		}

		Integral I2(sq1,6);
		cout << I2.getName() << endl;

		vector<Integral> l;
		sq1.getIntegralList(l);
		for(int i=0; i<l.size(); i++) {
		cout << l[i].getName() << endl;
		}

		Integral I(test,test,test,test,ERI);
		if (sq1.hasThisIntegral(I)){
		cout << "Has this intergal!" << endl;
		}else{
		cout << "do not Has this intergal!" << endl;
		}

		RRBuild rr(OS,ERI,BRA1);
		ShellQuartet sq(2,1,2,3,ERI,0);
		cout << "original name " << sq.getName() << endl;
		vector<ShellQuartet> sqlist;
		rr.buildRRSQ(sq,sqlist); 
		for(int i=0; i<sqlist.size();i++) {
		cout << "name " << sqlist[i].getName() << endl;
		}

		RRBuild rr2(OS,KINETIC,BRA2);
		ShellQuartet sq2(2,2,-1,-1,KINETIC,0);
		cout << "original name " << sq2.getName() << endl;
		vector<ShellQuartet> sqlist2;
		rr2.buildRRSQ(sq2,sqlist2); 
		for(int i=0; i<sqlist2.size();i++) {
		cout << "name " << sqlist2[i].getName() << endl;
		}

		RRBuild rr3(OS,KINETIC,BRA1);
		ShellQuartet sq3(2,2,-1,-1,KINETIC,0);
		cout << "original name " << sq3.getName() << endl;
		vector<ShellQuartet> sqlist3;
		rr3.buildRRSQ(sq3,sqlist3); 
		for(int i=0; i<sqlist3.size();i++) {
		cout << "name " << sqlist3[i].getName() << endl;
		}
		*/

	/*
	// testing the (II|II)
	ShellQuartet sq0 (6, 0,6,0,ERI,0);
	ShellQuartet sq1 (7, 0,6,0,ERI,0);
	ShellQuartet sq2 (8, 0,6,0,ERI,0);
	ShellQuartet sq3 (9, 0,6,0,ERI,0);
	ShellQuartet sq4 (10,0,6,0,ERI,0);
	ShellQuartet sq5 (11,0,6,0,ERI,0);
	ShellQuartet sq6 (12,0,6,0,ERI,0);

	ShellQuartet sq7 (6, 0,7, 0,ERI,0);
	ShellQuartet sq8 (6, 0,8, 0,ERI,0);
	ShellQuartet sq9 (6, 0,9, 0,ERI,0);
	ShellQuartet sq10(6, 0,10,0,ERI,0);
	ShellQuartet sq11(6, 0,11,0,ERI,0);
	ShellQuartet sq12(6, 0,12,0,ERI,0);

	ShellQuartet sq13(7, 0,7,0,ERI,0);
	ShellQuartet sq14(8, 0,7,0,ERI,0);
	ShellQuartet sq15(9, 0,7,0,ERI,0);
	ShellQuartet sq16(10,0,7,0,ERI,0);
	ShellQuartet sq17(11,0,7,0,ERI,0);
	ShellQuartet sq18(12,0,7,0,ERI,0);

	ShellQuartet sq19(7, 0,8, 0,ERI,0);
	ShellQuartet sq20(7, 0,9, 0,ERI,0);
	ShellQuartet sq21(7, 0,10,0,ERI,0);
	ShellQuartet sq22(7, 0,11,0,ERI,0);
	ShellQuartet sq23(7, 0,12,0,ERI,0);

	ShellQuartet sq24(8, 0,8,0,ERI,0);
	ShellQuartet sq25(9, 0,8,0,ERI,0);
	ShellQuartet sq26(10,0,8,0,ERI,0);
	ShellQuartet sq27(11,0,8,0,ERI,0);
	ShellQuartet sq28(12,0,8,0,ERI,0);

	ShellQuartet sq29(8, 0,9, 0,ERI,0);
	ShellQuartet sq30(8, 0,10,0,ERI,0);
	ShellQuartet sq31(8, 0,11,0,ERI,0);
	ShellQuartet sq32(8, 0,12,0,ERI,0);

	ShellQuartet sq33(9, 0,9,0,ERI,0);
	ShellQuartet sq34(10,0,9,0,ERI,0);
	ShellQuartet sq35(11,0,9,0,ERI,0);
	ShellQuartet sq36(12,0,9,0,ERI,0);

	ShellQuartet sq37(9,0,10,0,ERI,0);
	ShellQuartet sq38(9,0,11,0,ERI,0);
	ShellQuartet sq39(9,0,12,0,ERI,0);

	ShellQuartet sq40(10,0,10,0,ERI,0);
	ShellQuartet sq41(11,0,10,0,ERI,0);
	ShellQuartet sq42(12,0,10,0,ERI,0);
	ShellQuartet sq43(10,0,11,0,ERI,0);
	ShellQuartet sq44(10,0,12,0,ERI,0);

	ShellQuartet sq45(11,0,11,0,ERI,0);
	ShellQuartet sq46(12,0,11,0,ERI,0);
	ShellQuartet sq47(11,0,12,0,ERI,0);
	ShellQuartet sq48(12,0,12,0,ERI,0);

	vector<ShellQuartet> sqlist;
	sqlist.push_back(sq0);
	sqlist.push_back(sq1);
	sqlist.push_back(sq2);
	sqlist.push_back(sq3);
	sqlist.push_back(sq4);
	sqlist.push_back(sq5);
	sqlist.push_back(sq6);
	sqlist.push_back(sq7);
	sqlist.push_back(sq8);
	sqlist.push_back(sq9);
	sqlist.push_back(sq10);
	sqlist.push_back(sq11);
	sqlist.push_back(sq12);
	sqlist.push_back(sq13);
	sqlist.push_back(sq14);
	sqlist.push_back(sq15);
	sqlist.push_back(sq16);
	sqlist.push_back(sq17);
	sqlist.push_back(sq18);
	sqlist.push_back(sq19);
	sqlist.push_back(sq20);
	sqlist.push_back(sq21);
	sqlist.push_back(sq22);
	sqlist.push_back(sq23);
	sqlist.push_back(sq24);
	sqlist.push_back(sq25);
	sqlist.push_back(sq26);
	sqlist.push_back(sq27);
	sqlist.push_back(sq28);
	sqlist.push_back(sq29);
	sqlist.push_back(sq30);
	sqlist.push_back(sq31);
	sqlist.push_back(sq32);
	sqlist.push_back(sq33);
	sqlist.push_back(sq34);
	sqlist.push_back(sq35);
	sqlist.push_back(sq36);
	sqlist.push_back(sq37);
	sqlist.push_back(sq38);
	sqlist.push_back(sq39);
	sqlist.push_back(sq40);
	sqlist.push_back(sq41);
	sqlist.push_back(sq42);
	sqlist.push_back(sq43);
	sqlist.push_back(sq44);
	sqlist.push_back(sq45);
	sqlist.push_back(sq46);
	sqlist.push_back(sq47);
	sqlist.push_back(sq48);
	*/

		/*
		// testing the (DD|DD)
		ShellQuartet sq0(2,0,2,0,ERI,0);
		ShellQuartet sq1(3,0,2,0,ERI,0);
		ShellQuartet sq2(4,0,2,0,ERI,0);
		ShellQuartet sq3(2,0,3,0,ERI,0);
		ShellQuartet sq4(2,0,4,0,ERI,0);
		ShellQuartet sq5(3,0,3,0,ERI,0);
		ShellQuartet sq6(4,0,3,0,ERI,0);
		ShellQuartet sq7(3,0,4,0,ERI,0);
		ShellQuartet sq8(4,0,4,0,ERI,0);
		vector<ShellQuartet> sqlist;
		sqlist.push_back(sq0);
		sqlist.push_back(sq1);
		sqlist.push_back(sq2);
		sqlist.push_back(sq3);
		sqlist.push_back(sq4);
		sqlist.push_back(sq5);
		sqlist.push_back(sq6);
		sqlist.push_back(sq7);
		sqlist.push_back(sq8);
		*/

		// testing the (PP|PP)
		/*
			ShellQuartet sq0(1,0,1,0,ERI,0);
			ShellQuartet sq1(1,0,2,0,ERI,0);
			ShellQuartet sq2(2,0,1,0,ERI,0);
			ShellQuartet sq3(2,0,2,0,ERI,0);
			vector<ShellQuartet> sqlist;
			sqlist.push_back(sq0);
			sqlist.push_back(sq1);
			sqlist.push_back(sq2);
			sqlist.push_back(sq3);

			ShellQuartet sq(2,2,2,2,ERI,0);
			vector<ShellQuartet> sqlist;
			sqlist.push_back(sq);
			RRSQSearch rrsearch(sqlist,OS);
			*/

		/*
			ShellQuartet sq(3,3,2,2,ERI,0);
			vector<ShellQuartet> sqlist;
			sqlist.push_back(sq);
			RRSQSearch rrsearch(sqlist,HRR);
			*/
}
