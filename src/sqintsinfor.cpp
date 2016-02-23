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
#include "inttype.h"
#include "shell.h"
#include "shellsymbol.h"
#include "basis.h"
#include "integral.h"
#include "sqintsinfor.h"
#include "nonrr.h"
#include "vrrinfor.h"
#include <boost/algorithm/string.hpp>   // string handling
using namespace printing;
using namespace inttype;
using namespace shell;
using namespace basis;
using namespace integral;
using namespace nonrr;
using namespace vrrinfor;
using namespace sqintsinfor;

void SQIntsInfor::formSQInfor()
{
	//
	// firstly, let's get the number of bodies 
	// from the input shell code
	//
	int nBody = inputShellCodes.size();

	//
	// one body integrals
	//
	if (nBody == 1) {

		// decompose the composite shell quartets into single ones
		// get the number of single sqs
		int bra1  = inputShellCodes[0];
		int nBra1 = getNShells(bra1);
		inputSQList.reserve(nBra1);
		braCoeOffset.reserve(nBra1);

		// is it a composite shell quartet?
		bool isCompositeSQ = false;
		if (nBra1>1) {
			isCompositeSQ = true;
		}

		// get the angular momentum 
		int LminBra1 = -1;
		int LmaxBra1 = -1;
		decodeL(bra1,LminBra1,LmaxBra1);

		// now create the shell quartets
		int cOffset = 0;
		for(int LBra1 = LminBra1; LBra1<= LmaxBra1; LBra1++) {

			// create shells
			Shell sbra1(LBra1);
			Shell sbra2;
			Shell sket1;
			Shell sket2;

			// now create division
			long long div = -1;
			if (isCompositeSQ) {
				div = codeSQ(LBra1);
			}

			// now create shell quartets
			int M = 0;
			ShellQuartet sq(sbra1,sbra2,sket1,sket2,oper,M,div);
			inputSQList.push_back(sq);

			// now offset
			braCoeOffset.push_back(cOffset);
			cOffset++;
		}
	}

	//
	// two body integrals
	//
	if (nBody == 2) {

		// decompose the composite shell quartets into single ones
		// get the number of single sqs
		int bra1  = inputShellCodes[0];
		int bra2  = inputShellCodes[1];
		int nBra1 = getNShells(bra1);
		int nBra2 = getNShells(bra2);
		inputSQList.reserve(nBra1*nBra2);
		braCoeOffset.reserve(nBra1*nBra2);

		// is it a composite shell quartet?
		bool isCompositeSQ = false;
		if (nBra1>1 || nBra2>1) {
			isCompositeSQ = true;
		}

		// get the angular momentum 
		int LminBra1 = -1;
		int LmaxBra1 = -1;
		decodeL(bra1,LminBra1,LmaxBra1);
		int LminBra2 = -1;
		int LmaxBra2 = -1;
		decodeL(bra2,LminBra2,LmaxBra2);

		// now create the shell quartets
		int cOffset = 0;
		for(int LBra2 = LminBra2; LBra2<= LmaxBra2; LBra2++) {
			for(int LBra1 = LminBra1; LBra1<= LmaxBra1; LBra1++) {

				// create shells
				Shell sbra1(LBra1);
				Shell sbra2(LBra2);
				Shell sket1;
				Shell sket2;

				// now create division
				long long div = -1;
				if (isCompositeSQ) {
					div = codeSQ(LBra1,LBra2);
				}

				// now create shell quartets
				int M = 0;
				ShellQuartet sq(sbra1,sbra2,sket1,sket2,oper,M,div);
				inputSQList.push_back(sq);

				// now offset
				braCoeOffset.push_back(cOffset);
				cOffset++;
			}
		}
	}

	//
	// three body integrals
	//
	if (nBody == 3) {

		// decompose the composite shell quartets into single ones
		// get the number of single sqs
		int bra1  = inputShellCodes[0];
		int bra2  = inputShellCodes[1];
		int ket1  = inputShellCodes[2];
		int nBra1 = getNShells(bra1);
		int nBra2 = getNShells(bra2);
		int nKet1 = getNShells(ket1);
		inputSQList.reserve(nBra1*nBra2*nKet1);
		braCoeOffset.reserve(nBra1*nBra2*nKet1);   
		ketCoeOffset.reserve(nBra1*nBra2*nKet1);   

		// is it a composite shell quartet?
		bool isCompositeSQ = false;
		if (nBra1>1 || nBra2>1 || nKet1>1) {
			isCompositeSQ = true;
		}

		// get the angular momentum 
		int LminBra1 = -1;
		int LmaxBra1 = -1;
		decodeL(bra1,LminBra1,LmaxBra1);
		int LminBra2 = -1;
		int LmaxBra2 = -1;
		decodeL(bra2,LminBra2,LmaxBra2);
		int LminKet1 = -1;
		int LmaxKet1 = -1;
		decodeL(ket1,LminKet1,LmaxKet1);

		// now create the shell quartets
		int ketCOffset = 0;
		for(int LKet1 = LminKet1; LKet1<= LmaxKet1; LKet1++) {
			int braCOffset = 0;
			for(int LBra2 = LminBra2; LBra2<= LmaxBra2; LBra2++) {
				for(int LBra1 = LminBra1; LBra1<= LmaxBra1; LBra1++) {

					// create shells
					Shell sbra1(LBra1);
					Shell sbra2(LBra2);
					Shell sket1(LKet1);
					Shell sket2;

					// now create division
					long long div = -1;
					if (isCompositeSQ) {
						div = codeSQ(LBra1,LBra2,LKet1);
					}

					// now create shell quartets
					int M = 0;
					ShellQuartet sq(sbra1,sbra2,sket1,sket2,oper,M,div);
					inputSQList.push_back(sq);

					// now offset
					braCoeOffset.push_back(braCOffset);
					ketCoeOffset.push_back(ketCOffset);
					braCOffset++;
				}
			}
			ketCOffset++;
		}
	}

	//
	// four body integrals
	//
	if (nBody == 4) {

		// decompose the composite shell quartets into single ones
		// get the number of single sqs
		int bra1  = inputShellCodes[0];
		int bra2  = inputShellCodes[1];
		int ket1  = inputShellCodes[2];
		int ket2  = inputShellCodes[3];
		int nBra1 = getNShells(bra1);
		int nBra2 = getNShells(bra2);
		int nKet1 = getNShells(ket1);
		int nKet2 = getNShells(ket2);
		inputSQList.reserve(nBra1*nBra2*nKet1*nKet2);
		braCoeOffset.reserve(nBra1*nBra2*nKet1*nKet2);
		ketCoeOffset.reserve(nBra1*nBra2*nKet1*nKet2);

		// is it a composite shell quartet?
		bool isCompositeSQ = false;
		if (nBra1>1 || nBra2>1 || nKet1>1 || nKet2>1) {
			isCompositeSQ = true;
		}

		// get the angular momentum 
		int LminBra1 = -1;
		int LmaxBra1 = -1;
		decodeL(bra1,LminBra1,LmaxBra1);
		int LminBra2 = -1;
		int LmaxBra2 = -1;
		decodeL(bra2,LminBra2,LmaxBra2);
		int LminKet1 = -1;
		int LmaxKet1 = -1;
		decodeL(ket1,LminKet1,LmaxKet1);
		int LminKet2 = -1;
		int LmaxKet2 = -1;
		decodeL(ket2,LminKet2,LmaxKet2);

		// now create the shell quartets
		int ketCOffset = 0;
		for(int LKet2 = LminKet2; LKet2<= LmaxKet2; LKet2++) {
			for(int LKet1 = LminKet1; LKet1<= LmaxKet1; LKet1++) {
				int braCOffset = 0;
				for(int LBra2 = LminBra2; LBra2<= LmaxBra2; LBra2++) {
					for(int LBra1 = LminBra1; LBra1<= LmaxBra1; LBra1++) {

						// create shells
						Shell sbra1(LBra1);
						Shell sbra2(LBra2);
						Shell sket1(LKet1);
						Shell sket2(LKet2);

						// now create division
						long long div = -1;
						if (isCompositeSQ) {
							div = codeSQ(LBra1,LBra2,LKet1,LKet2);
						}

						// now create shell quartets
						int M = 0;
						ShellQuartet sq(sbra1,sbra2,sket1,sket2,oper,M,div);
						inputSQList.push_back(sq);

						// now offset
						braCoeOffset.push_back(braCOffset);
						ketCoeOffset.push_back(ketCOffset);
						braCOffset++;
					}
				}
				ketCOffset++;
			}
		}
	}
}

void SQIntsInfor::formDerivSQList(const DerivInfor& derInfor, vector<ShellQuartet>& sqlist)
{
	//
	// there's an important note for forming the deriv sq list
	// the deriv sq list will be in dimension of 
	// derivSQList(nSQ,nDeriv)
	//
	// such arrangement is closely to the way we find the index
	// of given integral in the final results
	//
	// see the function of getOffset of SQIntsInfor for more 
	// information
	//


	// clear the input sqlist
	sqlist.clear();

	// compute the total number of deriv shell quartets
	int num = derInfor.getTotalNumDeriv();
	sqlist.reserve(inputSQList.size()*num);

	// now let's initilize the sq list
	if (derivOrder == 1) {
		int len = derInfor.getLen1stDerivInforArray(); 
		for(int i=0; i<len; i++) {
			const FirstOrderDerivInfor& infor = derInfor.get1stDerivInfor(i); 
			int l   = infor.getDerivDirLen();
			int pos = infor.getDerivPos();
			for(int j=0; j<l; j++) {
				int dir = infor.getDerivDirection(j);
				for(int iSQ=0; iSQ<(int)inputSQList.size(); iSQ++) {
					ShellQuartet sq(inputSQList[iSQ]);
					sq.addDerivInfor(pos,dir);
					sqlist.push_back(sq);
				}
			}
		}
	}else if (derivOrder == 2) {
		int len = derInfor.getLen2edDerivInforArray(); 
		for(int i=0; i<len; i++) {
			const SecondOrderDerivInfor& infor = derInfor.get2edDerivInfor(i); 
			int l   = infor.getDerivDirLen();
			int pos1= NULL_POS;
			int pos2= NULL_POS;
			infor.getDerivPos(pos1,pos2);
			for(int j=0; j<l; j++) {
				int dir1 = NO_DERIV;
				int dir2 = NO_DERIV;
				infor.getDerivDirection(j,dir1,dir2);
				for(int iSQ=0; iSQ<(int)inputSQList.size(); iSQ++) {
					ShellQuartet sq(inputSQList[iSQ]);
					sq.addDerivInfor(pos1,pos2,dir1,dir2);
					sqlist.push_back(sq);
				}
			}
		}
	}
}

void SQIntsInfor::formDerivInfor()
{
	//
	// derivOrder is already checked in the infor class,
	// so we here do not check it again
	//
	
	//
	// initialize the records
	// this is used when we evaluate the redundant position
	//
	derivRecords[0] = 0;
	derivRecords[1] = 0;
	derivRecords[2] = 0;
	derivRecords[3] = 0;

	// 
	// let's check whether we need to find 
	// which is the redundant position for derivatives?
	// here we use the operator order to determine the 
	// one body integral case
	//
	int nBody = getOperOrder(oper);
	if (nBody == 1 || ! hasDerivRedundantPos(oper)) {

		// now create the deriv infor
		// and form the derivSQList
		DerivInfor deriv(oper,derivOrder);
		derivInfor = deriv;
		formDerivSQList(derivInfor,derivSQList);
		return;
	}

	//
	// for the symmetrical shell quartets, which all
	// of shell components equal to each other; we 
	// also do not need to do evaluations
	// just pick up the last possible position
	//
	int redundantPosition = NULL_POS;
	int nCodes = inputShellCodes.size();
	if (nCodes==2) {
		if (inputShellCodes[0] == inputShellCodes[1]) {
			redundantPosition = BRA2;
		}
	}else if (nCodes == 3) {	
		if (inputShellCodes[0] == inputShellCodes[1] && inputShellCodes[0] == inputShellCodes[2]) {
			redundantPosition = KET1;
		}
	}else if (nCodes == 4) {	
		if (inputShellCodes[0] == inputShellCodes[1] && 
				inputShellCodes[0] == inputShellCodes[2] && 
				inputShellCodes[0] == inputShellCodes[3]) {
			redundantPosition = KET2;
		}
	}

	// now create the deriv infor
	// and form the derivSQList
	// for symmetrical shell quartet case
	if (redundantPosition != NULL_POS) {
		DerivInfor deriv(oper,derivOrder,redundantPosition);
		derivInfor = deriv;
		formDerivSQList(derivInfor,derivSQList);
		return;
	}

	//
	// now we need to guess which one is best to be the 
	// redundant position in the derivatives information
	// therefore, we will try every possible derivative 
	// position, and update the correspoding record value
	//
	redundantPosition = NULL_POS;
	size_t nTotalInts = 0;
	for(int iBody=0; iBody<nBody; iBody++) {

		// select the potential redundant position
		int pos = BRA1;
		if (iBody == 1) pos = BRA2;
		if (iBody == 2) pos = KET1;
		if (iBody == 3) pos = KET2;
		DerivInfor deriv(oper,derivOrder,pos);

		// form the sq list
		vector<ShellQuartet> sqlist;
		formDerivSQList(deriv,sqlist);

		// form the integral list
		// basically, this is full integral list
		vector<set<int> > intList;
		intList.reserve(sqlist.size());
		for(int iSQ=0; iSQ<(int)sqlist.size(); iSQ++) {
			set<int> list;
			sqlist[iSQ].getIntegralList(list);
			intList.push_back(list);
		}

		// now use NONRR to evaluate the integral number
		NONRR nonRR(sqlist,intList);
		size_t totalNum = nonRR.evalDerivIntProcess(*this);
		derivRecords[iBody] = totalNum;

		// shall we assign the first value?
		if (iBody == 0) {
			redundantPosition = BRA1;
			nTotalInts = totalNum;
		}else{
			if (totalNum<=nTotalInts) {
				redundantPosition = pos;
				nTotalInts = totalNum;
			}
		}
	}

	// set the minimum derivative integral number
	minDerivInts = nTotalInts;

	// now let's generate the data
	DerivInfor deriv(oper,derivOrder,redundantPosition);
	derivInfor = deriv;
	formDerivSQList(derivInfor,derivSQList);
}

SQIntsInfor::SQIntsInfor(const int& oper0, 
		const Infor& infor0, const int& codeBra1, const int& codeBra2, 
		const int& codeKet1,const int& codeKet2):Infor(infor0),withArray(false),
	doHRRWork(true),sectionInfor(6,NULL_POS),oper(oper0),minDerivInts(0)
{
	//
	// firstly, from the input shell code let's form the data
	//
	inputShellCodes.reserve(4);
	inputShellCodes.push_back(codeBra1);
	if (codeBra2>=0) inputShellCodes.push_back(codeBra2);
	if (codeKet1>=0) inputShellCodes.push_back(codeKet1); 
	if (codeKet2>=0) inputShellCodes.push_back(codeKet2);
	formSQInfor();

	//
	// now let's see whether we form the derivatives information
	//
	if (derivOrder>0) {
		formDerivInfor();
	}

	// let's see whether we do HRR work
	doHRRWork = weDOHRRWork();

	// clear the section information
	// prepare to add in running time
	sectionInfor.clear();
}

int SQIntsInfor::getVRRContDegree() const
{
	int con = 1;
	for(int iCode=0; iCode<(int)inputShellCodes.size();iCode++) {
		int code = inputShellCodes[iCode];
		int degree = getContractionDegree(code);
		con *= degree;
	}
	return con;
}

bool SQIntsInfor::isResult(const ShellQuartet& sq) const 
{
	if (derivOrder == 0) {
		for(int iSQ=0; iSQ<(int)inputSQList.size(); iSQ++) {
			if (inputSQList[iSQ] == sq) return true;
		}
	}else{
		for(int iSQ=0; iSQ<(int)derivSQList.size(); iSQ++) {
			if (derivSQList[iSQ] == sq) return true;
		}
	}
	return false;
}

int SQIntsInfor::nInts() const
{
	int nTolInts = 0;
	for(int iSQ=0; iSQ<(int)inputSQList.size(); iSQ++) {
		nTolInts += inputSQList[iSQ].getNInts(); 
	}

	// now let's the job order
	// for order = 0, this is just the result
	// for the derivatives case, we need to 
	// multiply the derivatives number
	if (derivOrder == 0) {
		return nTolInts;
	}
	int num = derivInfor.getTotalNumDeriv();
	return nTolInts*num;
}

int SQIntsInfor::getOffset(const ShellQuartet& sq, const int& index) const
{
	// now let's see the derivatives situation
	// because the derivSQList is in dimension of:
	// derivSQList(nSQ,nDeriv)
	// thus we will need to see where is the deriv section
	// in the nDeriv
	int derivOffset = 0;
	if (derivOrder>0) {

		// let's see where is the deriv information related 
		// to the sq during the whole deriv infor
		int derivPos = -1;
		if (derivOrder == 1) {

			// fetch the deriv infor
			int pos = sq.get1stDerivPos();
			int dir = sq.get1stDerivDir();

			// if they are null information, this is wrong
			if (pos == NULL_POS || dir == NO_DERIV) {
				crash(true, "in getOffset of SQIntsInfor the deriv position and direction are NULL, for 1st order");
			}

			// now let's see what's the position of the deriv infor
			derivPos = derivInfor.getOffset(pos, dir); 
		}else if (derivOrder == 2) {

			// fetch the deriv infor
			int pos1 = sq.get1stDerivPos();
			int dir1 = sq.get1stDerivDir();
			int pos2 = sq.get2edDerivPos();
			int dir2 = sq.get2edDerivDir();

			// if they are null information, this is wrong
			if (pos1 == NULL_POS || dir1 == NO_DERIV || pos2 == NULL_POS || dir2 == NO_DERIV) {
				crash(true, "in getOffset of SQIntsInfor the deriv position and direction are NULL, for 2ed order");
			}

			// now let's see what's the position of the deriv infor
			derivPos = derivInfor.getOffset(pos1,pos2,dir1,dir2); 
		}

		// the return pos should not be -1
		if (derivPos == -1) {
			crash(true, "fail to get the deriv position in getOffset function of SQIntsInfor");
		}

		// now let's count how many integrals for each deriv section
		int nTolInts = 0;
		for(int iSQ=0; iSQ<(int)inputSQList.size(); iSQ++) {
			nTolInts += inputSQList[iSQ].getNInts(); 
		}

		// now set the derivOffset
		derivOffset = nTolInts*derivPos;
	}

	// now to compute the offset within the composite shell quartet
	// we need to retrieve the basis shell information from 
	// the index and sq
	Integral I(sq,index);

	// now let's count the dimension for the whole shell quartet
	// which may be composite one
	int n1 = 0;
	int n2 = 0;
	int n3 = 0;

	// bra1 is not null always
	int code  = inputShellCodes[0];
	int lmin1 = -1;
	int lmax1 = -1;
	decodeL(code,lmin1,lmax1);
	n1 = getCartBas(lmin1,lmax1);

	// bra2
	int lmin2 = -1;
	int lmax2 = -1;
	if (inputShellCodes.size() >= 2) {
		code = inputShellCodes[1]; 
		decodeL(code,lmin2,lmax2);
		n2 = getCartBas(lmin2,lmax2);
	}

	// ket1
	int lmin3 = -1;
	int lmax3 = -1;
	if (inputShellCodes.size() >= 3) {
		code = inputShellCodes[2]; 
		decodeL(code,lmin3,lmax3);
		n3 = getCartBas(lmin3,lmax3);
	}

	// ket4
	int lmin4 = -1;
	int lmax4 = -1;
	if (inputShellCodes.size() >= 4) {
		code = inputShellCodes[3]; 
		decodeL(code,lmin4,lmax4);
	}

	// now compute the bra1 shell position for the first shell
	int pos1 = 0;
	const Basis& b1 = I.getBasis(BRA1);
	int L = b1.getL();
	int offset = getShellOffsetInCompositeShell(lmin1,L);
	pos1 = offset + b1.getLocalIndex();

	// bra2
	int pos2 = 0;
	if (lmin2>=0) {
		const Basis& b2 = I.getBasis(BRA2);
		L = b2.getL();
		offset = getShellOffsetInCompositeShell(lmin2,L);
		pos2 = offset + b2.getLocalIndex();
	}

	// ket1
	int pos3 = 0;
	if (lmin3>=0) {
		const Basis& k1 = I.getBasis(KET1);
		L = k1.getL();
		offset = getShellOffsetInCompositeShell(lmin3,L);
		pos3 = offset + k1.getLocalIndex();
	}

	// ket2
	int pos4 = 0;
	if (lmin4>=0) {
		const Basis& k2 = I.getBasis(KET2);
		L = k2.getL();
		offset = getShellOffsetInCompositeShell(lmin4,L);
		pos4 = offset + k2.getLocalIndex();
	}


	// now it's the global index
	int globalIndex = derivOffset + pos1 + pos2*n1 + pos3*n1*n2 + pos4*n1*n2*n3;
	return globalIndex;
}

void SQIntsInfor::getCoeOffset(const ShellQuartet& sq, 
		int& ic2Offset, int& jc2Offset) const
{
	// we compare the division to see which sq in the 
	// input sq list is matching the given sq
	// we only compare the division, that is to find 
	// the division information so that we can calculate
	// the coefficient position
	//
	// here you must use the inputSQList to check the position
	// the derivSQList contains repeat shell quartets in different
	// derivatives order
	//
	long long division = sq.getDivision();
	int pos = -1;
	for(int iSQ=0; iSQ<(int)inputSQList.size(); iSQ++) {
		long long division2 = inputSQList[iSQ].getDivision();
		if (division2 == division) {
			pos = iSQ;
			break;
		}
	}

	// now let's get the c2 coefficients
	// bra side is always existing
	ic2Offset = braCoeOffset[pos];

	// however, we may not have ket side
	// for one or two body integrals
	// here we need to use getOperOrder
	// oper like MOM it's treated like 
	// three body integrals, however;
	// it's RR expansion only on BRA side
	// therefore in terms of RR it's two
	// body integrals
	// in this sense we use oper order
	int nBody = getOperOrder(oper);
	if (nBody <= 2) {
		jc2Offset = -1;
		return;
	}else{
		jc2Offset = ketCoeOffset[pos];
	}
}

string SQIntsInfor::getWorkFuncName(bool inFuncName, int moduleName, 
		int iFile, bool finalFile) const
{
	// get the function name
	string file = getFuncName();

	// the above file name is the main body, now according to the 
	// module name we need revise
	string append;
	if (moduleName > 0) {
		if (moduleName == VRR) {
			append = "_vrr";
		}else if (moduleName == VRR_FUNC_STATEMENT) {
			append = "_vrrfunc";
		}else if (moduleName == VRR_HEAD) {
			append = "_vrrhead";
		}else if (moduleName == VRR_CONT) {
			append = "_vrrcont";
		}else if (moduleName == VRR_CONT_STATEMENT) {
			append = "_vrrcontfunc";
		}else if (moduleName == HRR1) {
			append = "_hrr1";
		}else if (moduleName == HRR1_FUNC_STATEMENT) {
			append = "_hrr1func";
		}else if (moduleName == HRR2) {
			append = "_hrr2";
		}else if (moduleName == HRR2_FUNC_STATEMENT) {
			append = "_hrr2func";
		}else if (moduleName == NON_RR) {
			append = "_nonrr";
		}else if (moduleName == NON_RR_FUNC_STATEMENT) {
			append = "_nonrrfunc";
		}else if (moduleName == DERIV) {
			append = "_deriv";
		}else if (moduleName == DERIV_FUNC_STATEMENT) {
			append = "_derivfunc";
		}else{
			cout << "module name: " << moduleName << endl;
			crash(true, "invalid module name given in SQIntsInfor::getWorkFuncName");
		}
	}

	// do we have the iFile defined?
	if (iFile >= 0) {
		append = append + "_" + boost::lexical_cast<string>(iFile);
	}

	// form the function name
	if (append.size() > 0) {
		file = file+ append;
	}

	// do we just need to the function name?
	if (inFuncName) {
		return file;
	}

	// finally make the file with suffix of cpp
	file = file + ".cpp";

	// now let's add in the dir information
	// if the module name is NULL, this must be the main cpp file
	// we do not place it in the tmp work dir
	bool withTmpWorkDir = true;
	if (finalFile) withTmpWorkDir = false;
	string f = getProjectFileDir(oper,file,derivOrder,withTmpWorkDir);
	return f;
}

string SQIntsInfor::getFuncName() const
{
	// first part of the file name is the hrr, vrr and the operator
	// then it's the shell symbol
	string file  = getProjectName();
	string oname = getOperStringName(oper);
	boost::algorithm::to_lower(oname);
	file = file + "_" + oname;
	for(int i=0; i<(int)inputShellCodes.size(); i++) {
		string shellSym = getShellName(inputShellCodes[i]);
		boost::algorithm::to_lower(shellSym);
		file = file + "_" + shellSym;
	}

	// finally, it's the derivative order
	if (derivOrder == 1) {
		string order = "d1";
		file = file + "_" + order;
	}else if (derivOrder == 2) {
		string order = "d2";
		file = file + "_" + order;
	}
	return file;
}

int SQIntsInfor::getCoeArrayLength(const int& side) const
{
	//
	// the length of coe array could be determined
	// from the input shell codes 
	// we can see that how many shell pairs for Bra/Ket side
	//
	
	//
	// some special case will be handled first
	// we note, that for integrals BRA side is always
	// exist, however; ket side may not
	// for MOM, even though it has KET1 position and 
	// used in RR process, it's two body operator
	// therefore if the side is ket, we need to have 
	// an eye to see that whether the key part exists
	// so we may return 0
	//
	int nBody = getOperOrder(oper);
	if (side == KET) {
		if (nBody<=2) return 0;
	}

	// angular codes
	int c1 = -1;
	int c2 = -1;
	if (side == BRA || side == BRA1 || side == BRA2) {
		if (nBody >= 1) c1 = inputShellCodes[0];  
		if (nBody >= 2) c2 = inputShellCodes[1];  
	}else if (side == KET || side == KET1 || side == KET2) {
		if (nBody >= 3) c1 = inputShellCodes[2];  
		if (nBody >= 4) c2 = inputShellCodes[3];  
	}else{
		cout << "side is " << side << endl;
		crash(true, "Input side is wrong, should be only BRA/KET");
	}

	// now get the number of shell
	int n1 = 0;
	int n2 = 0;
	if (c1 >= 0) n1 = getNShells(c1); 
	if (c2 >= 0) n2 = getNShells(c2); 

	// now return the result
	// we also do additional check
	if (n1 == 0) {
		crash(true, "Something is wrong in getCoeArrayLength");
	}
	if (n2 == 0) return n1;
	return n1*n2;
}

bool SQIntsInfor::areAllBottomSQ() const 
{
	if (derivOrder == 0 && ! isNONRROper(oper)) {
		for(int iSQ=0; iSQ<(int)inputSQList.size(); iSQ++) {
			const ShellQuartet& sq = inputSQList[iSQ];
			if (! sq.isSTypeSQ()) return false;
		}
		return true;
	}

	// here it's derivatives job
	// all shell quartets are not bottom SQ
	// even the S type of SQ
	return false;
}

bool SQIntsInfor::useBoostGamma() const
{
	// firstly, check the operator
	// if the operator does not involve incomplete gamma function
	// for bottom integral, then we do not use it
	if (useFmt(oper)) {

		// for the one with total L <=10, we do not use incomplete
		// gamma function. See the function of fmtIntegralsGeneration
		// in sqintsprint
		// get the L sum
		int maxL = 0;
		for(int iSQ=0; iSQ<(int)inputSQList.size(); iSQ++) {
			int LSum = inputSQList[iSQ].getLSum();
			if (LSum > maxL) maxL = LSum;
		}
		maxL += derivOrder;

		// maxL == 10 is the limiit
		if (maxL<=10) return false;
		return true;
	}	
	return false;
}

bool SQIntsInfor::weDOHRRWork() const
{
	// firstly, let's see whether HRR method is really
	// defined for use
	if (! defineHRR()) return false;

	// see whether we have operator can not do HRR?
	if (! canDOHRR(oper)) return false;

	// if this is one body integral, we do not need it too
	int nBody = inputShellCodes.size();
	if (nBody == 1) return false;

	// now let's see whether we have any non-RR section
	// whether the non-RR section just follows the VRR
	// section? this is only true when the input of 
	// non-RR section can be fully solved by VRR
	// consider the derivatives first
	if (getJobOrder() > 0) {

		// if derivatives is on the non-RR operator,
		// then as far as what we deal with now (three body KI);
		// we must do HRR
		if (isNONRROper(oper)) return true;

		// now is for all RR operators
		// for the deriv order == 1, it's only the this type of 
		// integrals (00|a), or (00|00) do not need to do HRR
		// for deriv order >=2, all integrals will involve
		// HRR step
		// we note that here this may not be fully accurate
		if (getJobOrder() == 1) {
			if (nBody == 3) {

				// for three body integral, test the LSum to see
				// whether it's (00|a)
				int LSum = 0;
				for(int iSQ=0; iSQ<(int)inputSQList.size(); iSQ++) {
					LSum += inputSQList[iSQ].getLSum(BRA);
				}
				if (LSum == 0) return false;
			}else{

				// for the other cases, namely two body and four body
				// integrals; let's test whether it's all bottom
				// integrals
				int LSum = 0;
				for(int iSQ=0; iSQ<(int)inputSQList.size(); iSQ++) {
					LSum += inputSQList[iSQ].getLSum();
				}
				if (LSum == 0) return false;
			}
		}

		// now for all of other cases, HRR is necessary
		return true;
	}

	// now let's consider the non-RR operator for deriv order = 0
	// in this case, we do not need to do it only when it's (00|a)
	// type of integral (a could > 0)
	if (isNONRROper(oper)) {
		if (oper == THREEBODYKI) {
			int LSum = 0;
			for(int iSQ=0; iSQ<(int)inputSQList.size(); iSQ++) {
				LSum += inputSQList[iSQ].getLSum(BRA);
			}
			if (LSum == 0) return false;
		}else{
			crash(true, "right now we only know non-RR operator is three body KI, hasHrr in sqintsinfor");
		}
	}

	// finally for RR case, deriv order is zero
	// in this case we just test the input shell quartet
	bool noHRR = true;
	for(int iSQ=0; iSQ<(int)inputSQList.size(); iSQ++) {
		const ShellQuartet& sq = inputSQList[iSQ];
		if(sq.canDoHRR(BRA)){
			noHRR = false;
			break;
		}
		if(sq.canDoHRR(KET)){
			noHRR = false;
			break;
		}
	}
	if (noHRR) return false;

	// now we just return true
	return true;
}

int SQIntsInfor::nextSection(const int& sec) const 
{
	// double check the section infor size
	if (sectionInfor.size() == 0) {
		crash(true, "illegal sectionInfor array size in nextSection of SQIntsInfor class, it's empty");
	}

	// this is the last section, no codes after this section
	if (sectionInfor[0] == sec) return NULL_POS;

	// now search the position
	int pos = -1;
	for(int i=0; i<(int)sectionInfor.size(); i++) {
		if (sectionInfor[i] == sec) {
			pos = i;
			break;
		}
	}

	// check
	if (pos < 0) {
		crash(true, "we did not find the input code section in nextSection of SQIntsInfor class");
	}

	// now let's return it
	// this is always exsiting
	return sectionInfor[pos-1];
}

string SQIntsInfor::getRedundantPos() const
{
	// firstly let's see whether we have a redundant position 
	// defined or not
	if (derivRecords[0] == 0 && derivRecords[1] == 0 && 
			derivRecords[2] == 0 && derivRecords[3] == 0) {
		return "NOT AVIALABLE";
	}

	string pos = "BRA1";
	size_t minNum = derivRecords[0];
	for(int i=1; i<4; i++) {
		if (derivRecords[i] > 0 && derivRecords[i]<minNum) {
			if (i == 1) pos = "BRA2";
			if (i == 2) pos = "KET1";
			if (i == 3) pos = "KET2";
		}
	}
	return pos;
}

void SQIntsInfor::headPrinting(ofstream& file) const 
{
	// first part, the comment of software license
	string line = "//";
	printLine(0,line,file);
	line = "// ";
	printLine(0,line,file);
	line = "// This code is generated from CPPINTS, a C++ program to generate the ";
	printLine(0,line,file);
	line = "// analytical integrals based on Gaussian form primitive functions. ";
	printLine(0,line,file);
	line = "// Copyright (C) 2015 The State University of New York at Buffalo ";
	printLine(0,line,file);
	line = "// This software uses the MIT license as below: ";
	printLine(0,line,file);
	line = "// ";
	printLine(0,line,file);
	line = "// Permission is hereby granted, free of charge, to any person obtaining ";
	printLine(0,line,file);
	line = "// a copy of this software and associated documentation files (the \"Software\"), ";
	printLine(0,line,file);
	line = "// to deal in the Software without restriction, including without limitation ";
	printLine(0,line,file);
	line = "// the rights to use, copy, modify, merge, publish, distribute, sublicense, ";
	printLine(0,line,file);
	line = "// and/or sell copies of the Software, and to permit persons to whom the Software ";
	printLine(0,line,file);
	line = "// is furnished to do so, subject to the following conditions: ";
	printLine(0,line,file);
	line = "// ";
	printLine(0,line,file);
	line = "// The above copyright notice and this permission notice shall be included in all ";
	printLine(0,line,file);
	line = "// copies or substantial portions of the Software. ";
	printLine(0,line,file);
	line = "// ";					    
	printLine(0,line,file);
	line = "// THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, ";
	printLine(0,line,file);
	line = "// INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR ";
	printLine(0,line,file);
	line = "// PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE ";
	printLine(0,line,file);
	line = "// FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR ";
	printLine(0,line,file);
	line = "// OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER  ";
	printLine(0,line,file);
	line = "// DEALINGS IN THE SOFTWARE. ";
	printLine(0,line,file);
	line = "// ";
	printLine(0,line,file);
	line = "// ";
	printLine(0,line,file);
	file << endl;

	// do we use vector or scr etc.?
	bool usescr = withSCRVec();
	bool useDoubleVec = withDoubleVec();
	bool withTBBVec = useTBBVec();
	bool withBoostGamma = useBoostGamma();

	// now second part, it's the include file
	line = "#include \"constants.h\""; 
	printLine(0,line,file);
	line = "#include <cstddef>"; 
	printLine(0,line,file);
	line = "#include <math.h>"; 
	printLine(0,line,file);
	if (withBoostGamma) {
		line = "#include <boost/math/special_functions/gamma.hpp>";
		printLine(0,line,file);
	}
	if (usescr) {
		line = "#include \"localmemscr.h\""; 
		printLine(0,line,file);
		line = "using namespace localmemscr;";
		printLine(0,line,file);
	}else {
		if (useDoubleVec){
			line = "#include <vector>"; 
			printLine(0,line,file);
			if (withTBBVec) {
				line = "#include \"tbb/scalable_allocator.h\"";
				printLine(0,line,file);
			}
		}
	}
	file << endl;

	// now do the typedef work
	// however, the information here is defined in the localmemscr
	// therefore, we do not repeat it if we use LocalMemScr
	if (! usescr) {
		line = "typedef int             Int;";
		printLine(0,line,file);
		line = "typedef size_t          UInt;";
		printLine(0,line,file);
		line = "#ifdef WITH_SINGLE_PRECISION";
		printLine(0,line,file);
		line = "typedef float           Double;";
		printLine(0,line,file);
		line = "#define THRESHOLD_MATH  0.0000001";
		printLine(0,line,file);
		line = "#else";
		printLine(0,line,file);
		line = "typedef double          Double;";
		printLine(0,line,file);
		line = "#define THRESHOLD_MATH  0.00000000000001";
		printLine(0,line,file);
		line = "#endif";
		printLine(0,line,file);
		file << endl;
	}

	// now if we use vector, we also do typedef here
	if (useDoubleVec) {
		line = "// typedef the vector type so that they all have the same type of DoubleVec";
		if (withTBBVec) {
			line = "typedef std::vector<Double,tbb::scalable_allocator<Double> >   DoubleVec; ";
		}else{
			line = "typedef std::vector<Double>   DoubleVec; ";
		}
		file << endl;
	}

	// print out variable comments
	line = "//";
	printLine(0,line,file);
	line = "//  here below is a list of variables used in the program";
	printLine(0,line,file);
	line = "//";
	printLine(0,line,file);
	line = "//  alpha is the bra1's exponent";
	printLine(0,line,file);
	line = "//  beta  is the bra2's exponent";
	printLine(0,line,file);
	line = "//  gamma is the ket1's exponent";
	printLine(0,line,file);
	line = "//  delta is the ket2's exponent";
	printLine(0,line,file);
	line = "//  A is the nuclear center for bra1";
	printLine(0,line,file);
	line = "//  B is the nuclear center for bra2";
	printLine(0,line,file);
	line = "//  C is the nuclear center for ket1";
	printLine(0,line,file);
	line = "//  D is the nuclear center for ket2";
	printLine(0,line,file);
	line = "//  P is the new center after bra1 combined with bra2";
	printLine(0,line,file);
	line = "//  Q is the new center after ket1 combined with ket2";
	printLine(0,line,file);
	line = "//  W is the new center after P combined with Q";
	printLine(0,line,file);
	line = "//  thresh value is threshold to perform significance check on primitive integrals";
	printLine(0,line,file);
	line = "//  pMax is maximum value of corresponding density matrix block(or value pair), used for ERI";
	printLine(0,line,file);
	line = "//  omega is the exponential factor used for operator in form of erf(omega*r12)/r12";
	printLine(0,line,file);
	line = "//  also omega could be the exponential factor used for operator in form of e^(-omega*r12^2)";
	printLine(0,line,file);
	line = "//";
	printLine(0,line,file);
	line = "//  variables:";
	printLine(0,line,file);
	line = "//";
	printLine(0,line,file);
	line = "//  zeta      = alpha + beta";
	printLine(0,line,file);
	line = "//  eta       = gamma + delta";
	printLine(0,line,file);
	line = "//  oned2z    = 1/(2*zeta)";
	printLine(0,line,file);
	line = "//  oned2e    = 1/(2*eta)";
	printLine(0,line,file);
	line = "//  onedz     = 1/zeta";
	printLine(0,line,file);
	line = "//  onede     = 1/eta";
	printLine(0,line,file);
	line = "//  kappa     = zeta + eta";
	printLine(0,line,file);
	line = "//  onedk     = 1/kappa";
	printLine(0,line,file);
	line = "//  oned2zeta = 1/(2*(alpha+beta+gamma))";
	printLine(0,line,file);
	line = "//  xi        = alpha*beta*onedz";
	printLine(0,line,file);
	line = "//  twoxi     = 2*alpha*beta*onedz";
	printLine(0,line,file);
	line = "//  rho       = zeta*eta*onedk";
	printLine(0,line,file);
	line = "//  rhod2zsq  = rho/(2*zeta*zeta)";
	printLine(0,line,file);
	line = "//  rhod2esq  = rho/(2*eta*eta)";
	printLine(0,line,file);
	line = "//  odorho    = omega/(rho+omega)";
	printLine(0,line,file);
	line = "//  rhodorho  = rho/(rho+omega)";
	printLine(0,line,file);
	line = "//  orhod2z2  = (rho/(2*zeta*zeta))*(omega/(rho+omega))";
	printLine(0,line,file);
	line = "//  orhod2e2  = (rho/(2*eta*eta))*(omega/(rho+omega))";
	printLine(0,line,file);
	line = "//  od2k      = (1/(2*kappa))*(omega/(rho+omega))";
	printLine(0,line,file);
	line = "//  adz       = alpha*onedz";
	printLine(0,line,file);
	line = "//  bdz       = beta*onedz";
	printLine(0,line,file);
	line = "//  gde       = gamma*onede";
	printLine(0,line,file);
	line = "//  gde       = delta*onede";
	printLine(0,line,file);
	line = "//";
	printLine(0,line,file);
	line = "//  input parameters based on primitive functions pair:";
	printLine(0,line,file);
	line = "//";
	printLine(0,line,file);
	line = "//  bra side shell pair is index as i";
	printLine(0,line,file);
	line = "//  inp2  is the number of primitive pairs";
	printLine(0,line,file);
	line = "//  iexp  is the array of 1/(alpha+beta)";
	printLine(0,line,file);
	line = "//  icoe  is the array of ic_bra1*ic_bra2";
	printLine(0,line,file);
	line = "//  ifac  is the array of pre-factor on bra side ";
	printLine(0,line,file);
	line = "//  for (SS|SS)^{m} etc. type of integrals";
	printLine(0,line,file);
	line = "//  ket side shell pair is index as j";
	printLine(0,line,file);
	line = "//  jnp2  is the number of primitive pairs";
	printLine(0,line,file);
	line = "//  jexp  is the array of 1/(gamma+delta)";
	printLine(0,line,file);
	line = "//  jcoe  is the array of jc_ket1*jc_ket2";
	printLine(0,line,file);
	line = "//  jfac  is the array of pre-factor on ket side ";
	printLine(0,line,file);
	line = "//  for (SS|SS)^{m} etc. type of integrals";
	printLine(0,line,file);
	line = "//";
	printLine(0,line,file);
	file << endl;

	// derivatives information
	if (getJobOrder() > 0) {
		line = "//";
		printLine(0,line,file);
		line = "// print out the information regarding of derivatives ";
		printLine(0,line,file);
		line = "// here below we count on all of RHS integrals(including the repeat ones)";
		printLine(0,line,file);
		line = "// this is used to simulate the FLOPS counting";
		printLine(0,line,file);
		line = "// BRA1 as redundant position, total RHS integrals evaluated as: " + getRedundantIntEvalNumber(0);
		printLine(0,line,file);
		line = "// BRA2 as redundant position, total RHS integrals evaluated as: " + getRedundantIntEvalNumber(1);
		printLine(0,line,file);
		line = "// KET1 as redundant position, total RHS integrals evaluated as: " + getRedundantIntEvalNumber(2);
		printLine(0,line,file);
		line = "// KET2 as redundant position, total RHS integrals evaluated as: " + getRedundantIntEvalNumber(3);
		printLine(0,line,file);
		line = "// the redundant position is: " + getRedundantPos();
		printLine(0,line,file);
		line = "//";
		printLine(0,line,file);
		file << endl;

		// now print out the derivative position and dir information
		line = "//";
		printLine(0,line,file);
		line = "// @@@@ derivative position-direction information";
		printLine(0,line,file);
		const DerivInfor& deriv = getDerivInfor();
		if (getJobOrder() == 1) {
			int len = deriv.getLen1stDerivInforArray();
			for(int i=0; i<len; i++) {
				const FirstOrderDerivInfor& firstDeriv = deriv.get1stDerivInfor(i);
				firstDeriv.print(file);
				line = "// ####";
				printLine(0,line,file);
			}
		}else if (getJobOrder() == 2) {
			int len = deriv.getLen2edDerivInforArray();
			for(int i=0; i<len; i++) {
				const SecondOrderDerivInfor& secondDeriv = deriv.get2edDerivInfor(i);
				secondDeriv.print(file);
				line = "// ####";
				printLine(0,line,file);
			}
		}
		file << endl;
	}
}

string SQIntsInfor::getArgList() const
{
	string arg;
	int intOperator = getOper();

	// this is for input without exponential factors
	// since we need iexpdiff etc. to create alpha, beta etc.
	// variables for case with exponential factors
	if (withExpFac()) {
		switch(intOperator) {
			case TWOBODYOVERLAP:
				arg = "const UInt& inp2, const Double* icoe, "
					"const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, " 
					"const Double* A, const Double* B, Double* abcd";
				break;
			case MOM:
				arg = "const UInt& inp2, const Double* icoe, "
					"const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, " 
					"const Double* A, const Double* B, const Double* C, Double* abcd";
				break;
			case THREEBODYOVERLAP:
				arg = "const UInt& inp2, const UInt& jnp2, const Double* icoe, "
					"const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, "
					"const Double* A, const Double* B, const Double* jcoe, "
					"const Double* jexp, const Double* C, Double* abcd";
				break;
			case THREEBODYKI:
				arg = "const UInt& inp2, const UInt& jnp2, const Double* icoe, "
					"const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, "
					"const Double* A, const Double* B, const Double* jcoe, "
					"const Double* jexp, const Double* C, Double* abcd";
				break;
			case NAI:
				arg = "const UInt& inp2, const UInt& nAtoms, const Double* icoe, " 
					"const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, "
					"const Double* A, const Double* B, const Double* N, const UInt* Z, " 
					"Double* abcd";
				break;
			case ESP:
				arg = "const UInt& inp2, const UInt& nGrids, const Double* icoe, " 
					"const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, "
					"const Double* A, const Double* B, const Double* R, " 
					"Double* abcd";
				break;
			case ERI:
				arg = "const UInt& inp2, const UInt& jnp2, const Double& thresh, const Double& pMax, const Double& omega, "
					"const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, "
					"const Double* A, const Double* B, const Double* jcoe, "
					"const Double* jexp, const Double* jexpdiff, const Double* jfac, const Double* Q, "
					"const Double* C, const Double* D, Double* abcd";
				break;
			case EXPR12:
				arg = "const UInt& inp2, const UInt& jnp2, const Double& thresh, const Double& omega, "
					"const Double* icoe, const Double* iexp, const Double* iexpdiff, const Double* ifac, const Double* P, "
					"const Double* A, const Double* B, const Double* jcoe, "
					"const Double* jexp, const Double* jexpdiff, const Double* jfac, const Double* Q, "
					"const Double* C, const Double* D, Double* abcd";
				break;
			case KINETIC:
				arg = "const UInt& inp2, const Double* icoe, "
					"const Double* iexp, const Double* iexpdiff, const Double* ifac, " 
					"const Double* P, const Double* A, const Double* B, Double* abcd";
				break;
			default:
				crash(true, "Invalid operator passed in getArgList");
				break;
		}
	}else{	
		switch(intOperator) {
			case TWOBODYOVERLAP:
				arg = "const UInt& inp2, const Double* icoe, "
					"const Double* iexp, const Double* ifac, const Double* P, " 
					"const Double* A, const Double* B, Double* abcd";
				break;
			case MOM:
				arg = "const UInt& inp2, const Double* icoe, "
					"const Double* iexp, const Double* ifac, const Double* P, " 
					"const Double* A, const Double* B, const Double* C, Double* abcd";
				break;
			case THREEBODYOVERLAP:
				arg = "const UInt& inp2, const UInt& jnp2, const Double* icoe, "
					"const Double* iexp, const Double* ifac, const Double* P, "
					"const Double* A, const Double* B, const Double* jcoe, "
					"const Double* jexp, const Double* C, Double* abcd";
				break;
			case THREEBODYKI:
				arg = "const UInt& inp2, const UInt& jnp2, const Double* icoe, "
					"const Double* iexp, const Double* iexpdiff, "
					"const Double* ifac, const Double* P, "
					"const Double* A, const Double* B, const Double* jcoe, "
					"const Double* jexp, const Double* C, Double* abcd";
				break;
			case NAI:
				arg = "const UInt& inp2, const UInt& nAtoms, const Double* icoe, " 
					"const Double* iexp, const Double* ifac, const Double* P, "
					"const Double* A, const Double* B, const Double* N, const UInt* Z, " 
					"Double* abcd";
				break;
			case ESP:
				arg = "const UInt& inp2, const UInt& nGrids, const Double* icoe, " 
					"const Double* iexp, const Double* ifac, const Double* P, "
					"const Double* A, const Double* B, const Double* R, " 
					"Double* abcd";
				break;
			case ERI:
				arg = "const UInt& inp2, const UInt& jnp2, const Double& thresh, const Double& pMax, const Double& omega, "
					"const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, "
					"const Double* A, const Double* B, const Double* jcoe, "
					"const Double* jexp, const Double* jfac, const Double* Q, "
					"const Double* C, const Double* D, Double* abcd";
				break;
			case EXPR12:
				arg = "const UInt& inp2, const UInt& jnp2, const Double& thresh, const Double& omega, "
					"const Double* icoe, const Double* iexp, const Double* ifac, const Double* P, "
					"const Double* A, const Double* B, const Double* jcoe, "
					"const Double* jexp, const Double* jfac, const Double* Q, "
					"const Double* C, const Double* D, Double* abcd";
				break;
			case KINETIC:
				arg = "const UInt& inp2, const Double* icoe, "
					"const Double* iexp, const Double* iexpdiff, const Double* ifac, " 
					"const Double* P, const Double* A, const Double* B, Double* abcd";
				break;
			default:
				crash(true, "Invalid operator passed in getArgList");
				break;
		}
	}

	// finally, consider that whether we have the scr class add in?
	if (withSCRVec()) {
		arg = arg + ", LocalMemScr& scr";
	}

	return arg;
}

