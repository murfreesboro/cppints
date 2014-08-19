//
// CPPINTS: A C++ Program to Generate Analytical Integrals Based on Gaussian
// Form Primitive Functions
//
// Copyright (C) 2012-2014 Fenglai Liu
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
#include "shell.h"
#include "shellsymbol.h"
#include "sqintsinfor.h"
#include <boost/algorithm/string.hpp>   // string handling
using namespace inttype;
using namespace shell;
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
		sqOffsetList.reserve(nBra1);
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
		int offset  = 0;
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
			sqOffsetList.push_back(offset);
			braCoeOffset.push_back(cOffset);
			offset += sq.getNInts();
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
		sqOffsetList.reserve(nBra1*nBra2);
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
		int offset  = 0;
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
				sqOffsetList.push_back(offset);
				braCoeOffset.push_back(cOffset);
				offset += sq.getNInts();
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
		sqOffsetList.reserve(nBra1*nBra2*nKet1);
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
		int offset     = 0;
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
					sqOffsetList.push_back(offset);
					braCoeOffset.push_back(braCOffset);
					ketCoeOffset.push_back(ketCOffset);
					offset += sq.getNInts();
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
		sqOffsetList.reserve(nBra1*nBra2*nKet1*nKet2);
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
		int offset     = 0;
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
						sqOffsetList.push_back(offset);
						braCoeOffset.push_back(braCOffset);
						ketCoeOffset.push_back(ketCOffset);
						offset += sq.getNInts();
						braCOffset++;
					}
				}
				ketCOffset++;
			}
		}
	}
}

int SQIntsInfor::determineHRRPos(const int& side, 
		const int& LCode1, const int LCode2) const
{
	// decode the L
	int Lmin1 = -1;
	int Lmax1 = -1;
	decodeL(LCode1,Lmin1,Lmax1);
	int Lmin2 = -1;
	int Lmax2 = -1;
	decodeL(LCode2,Lmin2,Lmax2);

	// now let's go to see the number of basis sets
	int pos = -1;
	int nBas1 = getCartBas(Lmin1,Lmax1);
	int nBas2 = getCartBas(Lmin2,Lmax2);
	if (side == BRA) {
		pos = BRA2;
		if (nBas1<nBas2) pos = BRA1;
	}else{
		pos = KET2;
		if (nBas1<nBas2) pos = KET1;
	}
	return pos;
}

bool SQIntsInfor::hasSTypeSQ(int side) const
{
	// consider that whether the the given 
	// side is S? For this case, HRR is 
	// actually not meaningful
	if (side == BRA1 || side == BRA2 || side == KET1 || side == KET2) {

		// get the LCode
		int pos = 0;
		if (side == BRA2) pos = 1;
		if (side == KET1) pos = 2;
		if (side == KET2) pos = 3;
		int LCode = inputShellCodes[pos];

		// get the angular momentum
		int Lmin = -1;
		int Lmax = -1;
		decodeL(LCode,Lmin,Lmax);

		// is it S shell?
		if (Lmin == S && Lmax == S) {
			return true;
		}

	}else{

		// for general side, we just check whether they have 
		// S shell
		int LCode1 = -1;
		int LCode2 = -1;
		if (side == BRA) {
			LCode1 = inputShellCodes[0];
			LCode2 = inputShellCodes[1];
		}else{
			LCode1 = inputShellCodes[2];
			LCode2 = inputShellCodes[3];
		}

		// decode the L
		int Lmin1 = -1;
		int Lmax1 = -1;
		decodeL(LCode1,Lmin1,Lmax1);
		int Lmin2 = -1;
		int Lmax2 = -1;
		decodeL(LCode2,Lmin2,Lmax2);

		// is it S shell?
		if (Lmin1 == S && Lmax1 == S) {
			return true;
		}
		if (Lmin2 == S && Lmax2 == S) {
			return true;
		}
	}

	// finally it's worhty to do HRR
	// in terms of this case
	return false;
}

void SQIntsInfor::sideDeterminationInHRR() 
{
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

	}else if (nBody == 2 || nBody == 3) {

		// for two body and three body integrals,
		// ket side is not necessary. Therefore
		// we only need to set the first side
		secondSide = NULL_POS;

		// we set the position for the case of 
		// composite shell, else we will search
		// the specific position later in rrsqsearch
		// class
		int bra1 = inputShellCodes[0];
		int bra2 = inputShellCodes[1];
		if (isCompositeShell(bra1) || isCompositeShell(bra2)) {
			firstSide = determineHRRPos(BRA,bra1,bra2); 
		}else{
			firstSide = BRA;
		}

		// check the s type of integral
		bool hasSSQ = hasSTypeSQ(firstSide);
		if (hasSSQ) firstSide = NULL_POS;

	}else{

		//
		// for four body integrals, two things need to be 
		// solved:
		// 1  shall we do bra side first or ket side first?
		// 2  for each bra/ket side, which position we are 
		// going to do with HRR expansion if the given shell
		// quartet is composite one?
		// 

		// now test each side
		int nBraInts = -1;
		int nKetInts = -1;
		int braSide  = -1;
		int ketSide  = -1;
		for(int iSide=0; iSide<2; iSide++) {

			// set L code
			int L1 = inputShellCodes[0];
			int L2 = inputShellCodes[1];
			if (iSide == 1) {
				L1 = inputShellCodes[2];
				L2 = inputShellCodes[3];
			}

			// we set the position for the case of 
			// composite shell, else we will search
			// the specific position later in rrsqsearch
			if (isCompositeShell(L1) || isCompositeShell(L2)) {
				if (iSide == 0) {
					braSide = determineHRRPos(BRA,L1,L2); 
				}else{
					ketSide = determineHRRPos(KET,L1,L2); 
				}
			}else{
				if (iSide == 0) {
					braSide = BRA;
				}else{
					ketSide = KET;
				}
			}

			// now we need to know how many integrals here
			int Lmin1 = -1;
			int Lmax1 = -1;
			decodeL(L1,Lmin1,Lmax1);
			int Lmin2 = -1;
			int Lmax2 = -1;
			decodeL(L2,Lmin2,Lmax2);

			// now let's go to see the number of basis sets
			int nBas1 = getCartBas(Lmin1,Lmax1);
			int nBas2 = getCartBas(Lmin2,Lmax2);

			if (iSide == 0) {
				nBraInts = nBas1*nBas2;
			}else{
				nKetInts = nBas1*nBas2;
			}
		}

		// now we have to dicide whether bra or ket should 
		// be taken first. Simply compare the number of 
		// basis set pairs
		if (nBraInts<nKetInts) {
			firstSide  = braSide;
			secondSide = ketSide;
		}else{
			firstSide  = ketSide;
			secondSide = braSide;
		}

		// finally, check each side, whether it's not 
		// worthy for HRR
		bool hasSSQ = hasSTypeSQ(firstSide);
		if (hasSSQ) firstSide = NULL_POS;
		hasSSQ = hasSTypeSQ(secondSide);
		if (hasSSQ) secondSide = NULL_POS;
	}
}	

SQIntsInfor::SQIntsInfor(const int& job, const int& oper0, 
		const Infor& infor0, const int& codeBra1, const int& codeBra2, 
		const int& codeKet1,const int& codeKet2):Infor(infor0),
	vrrWithArrayIndex(false),hrrWithArrayIndex(false),
	splitFile(false),firstSide(NULL_POS),
	secondSide(NULL_POS),jobOrder(job),oper(oper0)
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
	// now let's get the maxL to proceed to other settings
	// determine the max angular momentum sum for each single sq
	//
	int maxL = 0;
	for(int iSQ=0; iSQ<(int)inputSQList.size(); iSQ++) {
		int LSum = inputSQList[iSQ].getLSum();
		if (LSum > maxL) maxL = LSum;
	}
	maxL += jobOrder;

	// get the operator inforrmation
	// determine the nBody size
	int nBody = getRROrder(oper);
	bool candoHRR = canDOHRR(oper);
	
	// now let's determine the rr printing option
	if (enforceHRR == NO_ENFORCE_ON_RR) {

		// this is natrual determining step
		if (maxL <= maxL_hrrPrinting) {
			hrrWithArrayIndex = false;
		}else{
			hrrWithArrayIndex = true;
		}

		// then let's go to see the nBody of integrals
		// if nBody <=2; we do not use array in HRR
		// part
		if(nBody<=2) hrrWithArrayIndex = false;

	}else if (enforceHRR ==  ENFORCE_RR_WITH_VAR) {
		hrrWithArrayIndex = false;
	}else{
		hrrWithArrayIndex = true;
	}

	// VRR
	if (enforceVRR == NO_ENFORCE_ON_RR) {

		// this is natrual determining step
		if (maxL <= maxL_vrrPrinting) {
			vrrWithArrayIndex = false;
		}else{
			vrrWithArrayIndex = true;
		}

		//
		// then let's go to see the nBody of integrals
		// for VRR part
		// 1  if HRR presents, then the VRR's part only do
		//    2 body integrals (the rest of position will be S shell).
		//    Therefore, if HRR presents we only use variable
		// 2  if HRR does not present, then we see the nBody size
		//    similar with HRR, var form applies to the case
		//    that nBody<=2
		//
		if (candoHRR) {
			vrrWithArrayIndex = false;
		}else{
			if(nBody<=2) vrrWithArrayIndex = false;
		}

	}else if (enforceVRR ==  ENFORCE_RR_WITH_VAR) {
		vrrWithArrayIndex = false;
	}else{
		vrrWithArrayIndex = true;
	}

	// do we split file according to the maxl?
	// we note, that since the splitFile option will
	// change the array usage, therefore the parsing 
	// of single file option should follow the array form
	// selection
	if (maxL>maxL_singleFile) splitFile = true;

	// additionally, for the split file situation
	// we may have some exceptions
	
	// case 1: all of input shell quartet are S type integrals
	// this is trivial to perform file split, so not do it 
	bool  allSSQ = true;
	if (splitFile) {

		// exception 1 if all of shell quartet are 
		// pure s integrals, then we actually do not
		// need to do split file (VRR and HRR is nothing)
		for(int iSQ=0; iSQ<(int)inputSQList.size(); iSQ++) {
			if (! inputSQList[iSQ].isSTypeSQ()) {
				allSSQ = false;
				break;
			}
		}
	}
	if (allSSQ) splitFile = false;

	// case 2: if it's already in the split file mode,
	// then the HRR part must use the array form
	// since we can not pass in too many variables
	// to mass up the code
	if (splitFile) {
		hrrWithArrayIndex = true;
	}

	// final setting for HRR
	// if it's actually no HRR part in RR process
	// then we set it to be default choice: false
	if (! canDOHRR(oper)) {
		hrrWithArrayIndex = false;
	}

	// finally, let's determine HRR side information
	sideDeterminationInHRR();
}

bool SQIntsInfor::isResult(const ShellQuartet& sq) const 
{
	for(int iSQ=0; iSQ<(int)inputSQList.size(); iSQ++) {
		if (inputSQList[iSQ] == sq) return true;
	}
	return false;
}

int SQIntsInfor::getOffset(const ShellQuartet& sq) const
{
	// get the position
	int pos = -1;
	for(int iSQ=0; iSQ<(int)inputSQList.size(); iSQ++) {
		if (inputSQList[iSQ] == sq) {
			pos = iSQ;
			break;
		}
	}
	if (pos == -1) return pos;
	return sqOffsetList[pos];
}

void SQIntsInfor::getCoeOffset(const ShellQuartet& sq, 
		int& ic2Offset, int& jc2Offset) const
{
	// we compare the division to see which sq in the 
	// input sq list is matching the given sq
	// we only compare the division, that is to find 
	// the division information so that we can calculate
	// the coefficient position
	long long division = sq.getDivision();
	int pos = -1;
	for(int iSQ=0; iSQ<(int)inputSQList.size(); iSQ++) {
		long long division2 = inputSQList[iSQ].getDivision();
		if (division2 == division) {
			pos = iSQ;
			break;
		}
	}
#ifdef SQINTS_DEBUG
	cout << "pos is: " << pos << endl;
	crash(pos == -1, "Fail to get the correct sq position in getCoeOffset");
#endif

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

string SQIntsInfor::getFileName(bool withTmpWorkDir) const
{

	// get the function name
	string file = getFuncName();

	// now let's add in the dir information
	string f = getProjectFileDir(oper,file,jobOrder,withTmpWorkDir);
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
	if (jobOrder == 1) {
		string order = "d1";
		file = file + "_" + order;
	}else if (jobOrder == 2) {
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
	for(int iSQ=0; iSQ<(int)inputSQList.size(); iSQ++) {
		const ShellQuartet& sq = inputSQList[iSQ];
		if (! sq.isSTypeSQ()) return false;
	}
	return true;
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
		maxL += jobOrder;

		// maxL == 10 is the limiit
		if (maxL<=10) return false;
		return true;
	}	
	return false;


}

