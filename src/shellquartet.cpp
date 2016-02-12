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
#include "basis.h"
#include "shell.h"
#include "inttype.h"
#include "boost/lexical_cast.hpp"
#include "integral.h"
#include "derivinfor.h"
#include "shellquartet.h"
using namespace basis;
using namespace shell;
using namespace inttype;
using namespace integral;
using namespace derivinfor;
using boost::lexical_cast;
using namespace shellquartet;

int ShellQuartet::hrrCompare(const ShellQuartet& sq, const int& expandingPos) const {

	// checking the operator first
	if (O != sq.getOper()) {
		crash(true,"Only same type of integral could be compared in hrr less than function");
	}

	// in the HRR part, the deriv information should be null
	// let me also check it
	int jobOrder = getDerivJobOrder();
	if (jobOrder>0) {
		crash(true,"In hrr less than function we should not have any deriv information defined");
	}

	//
	// in sorting the shell quartets in HRR section, we may want to arrange the 
	// shell quartets which shares the same exponential factors information together
	// therefore, we firstly judge with the exps. information first
	//
	// secondly, within the same deriv section, we sort with the division information,
	// which is a kind of arbitrary
	//

	// if this shell quartet has different exp factor length 
	// with input one, then we will compare it here
	// where the exp factor list length is larger, then
	// the shell quartet is with more exponential factors
	// we will try to do it later in the code context
	// therefore as length is smaller, it's greater 
	// when we comparing with input sq
	if (expFacListLen != sq.expFacListLen) {
		if (expFacListLen < sq.expFacListLen) {
			return GT;
		}else{
			return LT;
		}
	}

	// now we need to compare the expotential factors
	// if the shell quartets share the sme expoential factor length
	// but their arrangement is different, in HRR they will never 
	// mix with each other. Therefore we only provide an arbitray
	// way to compare it
	if (expFacListLen>0) {
		const int* expList = sq.getExpFacList();
		for(int i=0; i<expFacListLen; i++) {
			int pos1 = expFacList[i];
			int pos2 = expList[i];
			if (pos1 != pos2) {
				if (pos1<pos2) {
					return GT;
				}else{
					return LT;
				}
			}
		}
	}
	
	// different division of shell quartet can be compared in arbitrary way
	// if the division is smaller, then it's bigger
	// let's do it in this way
	if (division != sq.getDivision()) {
		if (division < sq.getDivision()) return GT;
		return LT;
	}

	// now below it's the core algorithm for HRR recursive comparison.
	// We note, that the shell quartet who is smaller should actually
	// appear in the LHS, and the shell quartet who is larger should 
	// appear in the RHS of HRR recursive expression. Therefore, 
	// we are able to generate correct shell quartet list in RR that 
	// sq0 -> sq1 -> sq2 ...
	//
	// Generally, the rule is like this(for the given bra/ket side):
	//
	// 1 we compare the L sum for the given side between the two 
	// shell quartets. Because in HRR, the L Sum is decreasing from
	// LHS to RHS
	//
	// 2 if we pass 1 then we need to refer to the help of input 
	// expanding position. We note, that because for all of HRR process 
	// the expanding position is same(please refer to the doc or the 
	// rrsqsearch.cpp for more information), hence for the input shell 
	// quartet and this shell quartet, their expanding position must 
	// be same, too.
	//
	// Based on this point, we will compare the L on the expanding 
	// position. When L is larger, it must be on LHS; L is lower,
	// it must be on RHS. This is because the HRR is just trying to
	// raise up the given position angular momentum.
	//
	// 3 if we go here, then it means the given side is totally same.
	// then we return to the above function for futher process by
	// returning "EQ(==)"
	//
	int side = BRA;
	if (expandingPos == KET1 || expandingPos == KET2) side = KET;

	// get L Sum
	// smaller one is on RHS
	// larger one is on LHS
	int LSum1 = getLSum(side);
	int LSum2 = sq.getLSum(side);
	if (LSum1<LSum2) {
		return GT;
	}else if (LSum1>LSum2) {
		return LT;
	}

	// now L sum is same get the L for the expanding position
	int L     = -1;
	int LSQ   = -1;
	if (side == BRA) {

		// first step, checking the bra side shell
		if (bra1.isnull() || bra2.isnull()) {
			crash(true,"If bra2/bra1 is none, why we go to the hrr compare in shell quartet?");
		}

		// get shell
		const Shell& b1 = sq.getShell(BRA1);
		const Shell& b2 = sq.getShell(BRA2);
		if (expandingPos == BRA1) {
			L     = bra1.getL();
			LSQ   = b1.getL();
		}else if (expandingPos == BRA2) {
			L     = bra2.getL();
			LSQ   = b2.getL();
		}else{
			crash(true, "in function of hrrCompare wrong expanding position provided on bra side");
		}

	}else{

		// first step, checking the ket side shell
		if (ket1.isnull() || ket2.isnull()) {
			cout << "the other sq name " << sq.getName() << endl;
			cout << "this sq name " << getName() << endl;
			crash(true,"If ket2/ket1 is none, why we go to the hrr compare in shell quartet?");
		}

		// get shell
		const Shell& k1 = sq.getShell(KET1);
		const Shell& k2 = sq.getShell(KET2);
		if (expandingPos == KET1) {
			L     = ket1.getL();
			LSQ   = k1.getL();
		}else if (expandingPos == KET2) {
			L     = ket2.getL();
			LSQ   = k2.getL();
		}else{
			crash(true, "in function of hrrCompare wrong expanding position provided on ket side");
		}
	}

	// compare the L and LSQ for the given expanding position
	// larger one is on LHS
	// smaller one is on RHS
	if (L > LSQ) {
		return LT; 
	}else if (L < LSQ) {
		return GT; 
	} else {
		return EQ;
	}
}

bool ShellQuartet::hrrLessThan(const ShellQuartet& sq, const int& position) const {

	int compareResult = hrrCompare(sq, position); 
	if (compareResult == GT) {
		return false;
	}else if (compareResult == LT) {
		return true;
	}else{

		// now we need to think about how we can do
		// it's obvious that in the first compare,
		// the given side is same between this sq
		// and the input sq
		//
		// if the other side is suitable for HRR
		// compare, that means; both of the two
		// positions exist, then we compare
		// the other side 
		// We will do it with an arbitrary expanding 
		// pisition of KET1
		//
		// if the other side can not do HRR compare,
		// for example; three body integrals;
		// then we just compare the L of the single
		// shell in that side. This is arbitrary

		// now let's see whether we can do the other side
		// we can not do, because for the other side
		// we have a null shell therefore HRR does not apply
		// to this situation
		int otherSide = KET;
		if (position == KET1 || position == KET2) otherSide = BRA;
		bool canDoHRRCompare = true;
		if (otherSide == BRA) {
			const Shell& b1SQ = sq.getShell(BRA1);
			const Shell& b2SQ = sq.getShell(BRA2);
			const Shell& b1   = getShell(BRA1);
			const Shell& b2   = getShell(BRA2);
			if (b1SQ.getL()<0 || b2SQ.getL()<0 || b1.getL()<0 || b2.getL()<0) {
				canDoHRRCompare = false;
			}
		}else{
			const Shell& b1SQ = sq.getShell(KET1);
			const Shell& b2SQ = sq.getShell(KET2);
			const Shell& b1   = getShell(KET1);
			const Shell& b2   = getShell(KET2);
			if (b1SQ.getL()<0 || b2SQ.getL()<0 || b1.getL()<0 || b2.getL()<0) {
				canDoHRRCompare = false;
			}
		}

		// now let's compare
		if (canDoHRRCompare) {
			int newPos = KET2;
			if (otherSide == BRA) newPos = BRA2;
			int compareResult2 = hrrCompare(sq,newPos);
			if (compareResult2 == GT) {
				return false;
			}else if (compareResult2 == LT) {
				return true;
			}else{
				cout << "SQ1 name " << getName() << endl;
				cout << "SQ2 name " << sq.getName() << endl;
				crash(true,"Something wrong with the two shell quartets in operator hrrlessthan");
				return false;
			}
		}else{
			// in this case we have a null shell
			// this must be the KET2 shell and let's examine it
			if (otherSide != KET) {
				crash(true, "soemthing wrong for hrrlessthan, the other side must be KET");
			}
			const Shell& b2SQ = sq.getShell(KET2);
			const Shell& b2   = getShell(KET2);
			if (b2SQ.getL()>=0 || b2.getL()>=0) {
				crash(true, "soemthing wrong for hrrlessthan, the KET2 side shell must be NULL??");
			}

			// OK, let's compare the KET1 
			const Shell& b1SQ = sq.getShell(KET1);
			const Shell& b1   = getShell(KET1);
			int L1SQ = b1SQ.getL();
			int L1   = b1.getL();
			if (L1<L1SQ) {
				return false;
			}else if (L1>L1SQ) {
				return true;
			}else{
				cout << "SQ1 name " << getName() << endl;
				cout << "SQ2 name " << sq.getName() << endl;
				crash(true,"Something wrong with hrrlessthan, we fail to compare two sq when KET2 shell is null");
				return false;
			}
		}
	}
}

bool ShellQuartet::operator<(const ShellQuartet& sq) const {

	//
	// important note: the operator < means the section of codes 
	// appearing later then the comparing code section for sq
	//

	//
	// firstly, let take care of NON-VRR situation
	// the derivatives is evaluated first
	//
	int jobOrder = getDerivJobOrder();
	if (jobOrder > 0) {

		// let's see the derivative order
		// the order more higher, then it should be 
		// appearing later in the final codes
		int thisDerivOrder = 0;
		if (firstDerivPos >0) thisDerivOrder++;
		if (secondDerivPos >0) thisDerivOrder++;
		int otherDerivOrder = 0;
		if (sq.firstDerivPos >0) otherDerivOrder++;
		if (sq.secondDerivPos >0) otherDerivOrder++;
		if (thisDerivOrder != otherDerivOrder) {
			return thisDerivOrder > otherDerivOrder;
		}

		// compare the first deriv pos first
		// we follow the order x->y->z
		if (firstDerivPos != sq.firstDerivPos) {
			return firstDerivPos > sq.firstDerivPos;
		}

		// now the first deriv position same
		// we compare the  second deriv position
		// we follow the order x->y->z
		if (secondDerivPos != sq.secondDerivPos) {
			return secondDerivPos > sq.secondDerivPos;
		}

		// now all of deriv position are same
		// compare the direction
		if (firstDerivDir != sq.firstDerivDir) {
			return firstDerivDir > sq.firstDerivDir;
		}

		// compare the second direction
		if (secondDerivDir != sq.secondDerivDir) {
			return secondDerivDir > sq.secondDerivDir;
		}

		// now all of deriative information are same
		// the only possilility here is that shell
		// quartets are composite
		// let's check it here
		if (division>=0 && sq.division>=0 && sq.division != division) {
			if (division < sq.getDivision()) return false;
			return true;
		}
	}

	//
	// The arrangement is used to sort the shell quartets in the VRR generation.
	// In fact, the shell quartet who is smaller is on the LHS of RR, and the 
	// shell quartet who is larger is on the RHS of RR. Through such arrangement
	// we will make sure that in the whole VRR process we can generate the reiable
	// shell quartet sequence:
	// sq0 -> sq1 -> sq2 ....
	// and sq0 is the result shell quartet, the last shell quartet is always 
	// the (SS|O|SS)^{m} type, which is the largest one according to our discussion
	// here
	//
	// additionally, we note that for VRR process(where the operator< is applied),
	// we do not care about the information of division and exponetial factor.
	// They should be destroyed in the VRR section
	//
	
	//comparison in terms of different operator type
	//small operator means this type of shell quartets
	//generate the other type of shell quartets in RR (like kinetic
	//integral generates overlap integral in OS RR), but the other 
	//type of shell quartets is self dependent
	if (O != sq.getOper()) {
		int state = compareOper(O,sq.getOper());
		if (state == LT) {
			return true;
		}else if (state == GT) {
			return false;
		}else{
			cout << "in operator < two shell quartets' operators invalid in compare" << endl;
			crash(true, "We should not have them here");
			return false;
		}
	}

	//comparison in terms of M value
	//smaller M value term is less than the 
	//larger M value term
	//since larger M value term is always
	//after the smaller M value term in RR generation
	if (getM() < sq.getM()) {
		return true;
	}else if (getM() > sq.getM()) {
		return false;
	}

	// now comparing their total L
	// one thing is obvious that if total L
	// is smaller, then it's must on RHS of OS
	// rather than LHS
	// therefore, if L sum is larger, then the 
	// shell quartet is smaller (in all of RR
	// process, it's above the larger one)
	int LSum1 = getLSum();
	int LSum2 = sq.getLSum();
	if (LSum1 > LSum2) {
		return true;
	}else if (LSum1 < LSum2) {
		return false;
	}
	
	//
	//	From the recursive expression, we arrange the VRR integrals in the 
	//	way below(take the S,P and D as example):
	// SSSS 
	// PSSS > PSPS > PPPS > PPPP
	// DSSS 
	// DSPS > DPSS > DPPS > DSPP > DPPP
	// DSDS > DDSS
	// DPDS > DDPS > DDPP > DPDP
	// DDDS > DDDP > DDDD 
	// the upper layer is larger than the lower layer. The lower layer
	// is always in the LHS, and the upper layer is in the RHS of RR.
	// The rule here to sort the shell quartets is like this:
	// 1 consider the largest angular momentum maxL1, if maxL1(sq1) > maxL1(sq2)
	// then sq1 < sq2, if maxL1(sq1) < maxL1(sq2) then sq1 > sq2; else we continue
	// to go to 2;
	// 2 consider the secondly largest angular momentum maxL2, if 
	// maxL2(sq1) > maxL2(sq2) then sq1 < sq2, if maxL2(sq1) < maxL2(sq2) then 
	// sq1 > sq2; else we continue to go to 3;
	// 3 consider the thirdly largest angular momentum maxL3, if 
	// maxL3(sq1) > maxL3(sq2) then sq1 < sq2, if maxL3(sq1) < maxL3(sq2) then 
	// sq1 > sq2; else we continue to go to 4;
	// 4 consider the smallest angular momentum maxL4, if 
	// maxL4(sq1) > maxL4(sq2) then sq1 < sq2, if maxL4(sq1) < maxL4(sq2) then 
	// sq1 > sq2; else we continue;
	// 5 now we have the shell quartet that their angular momentum are same,
	// but with different arrangment orders. For example, (44|22) and (42|42).
	// now we need to tell them from each other.
	// 5.1 if L(sq1's bra1) < L(sq2's bra1), sq1 > sq2, L(sq1's bra1) > L(sq2's bra1), 
	// sq1 < sq2; else we go to 5.2;
	// 5.2 if L(sq1's bra2) < L(sq2's bra2), sq1 > sq2, L(sq1's bra2) > L(sq2's bra2),
	// sq1 < sq2; else we go to 5.3;
	// 5.3 if L(sq1's ket1) < L(sq2's ket1), sq1 > sq2, L(sq1's ket1) > L(sq2's ket1),
	// sq1 < sq2; now we should get a result.
	// 5.3 should be the final end. Since in the above we have compare the situation
	// of different operators and different M values. Therefore in this section,
	// M and operator should be same. In VRR, we can not have two same shell quartets
	// in comparison. Therefore, 5.3 should be the end.
	//
	// Here we note, that the 5.1-5.3 comparison actually has nothing to do with 
	// RR. We do not need it into the fine level. For example, (44|22) and (42|42)
	// should be in the same level of RR - this is because they can not have 
	// overlapped RHS terms in RR generation. therefore, the comparison rule above
	// is defined in arbtrary way.
	//

	// get shell information for sq
	const Shell& b1 = sq.getShell(BRA1);
	const Shell& b2 = sq.getShell(BRA2);
	const Shell& k1 = sq.getShell(KET1);
	const Shell& k2 = sq.getShell(KET2);

	// now get angular momentum values:self
	int selfL[4];
	selfL[0] = bra1.getL();
	selfL[1] = bra2.getL();
	selfL[2] = ket1.getL();
	selfL[3] = ket2.getL();
	sort(selfL, selfL+4);

	// now get angular momentum values:sq
	int sqL[4];
	sqL[0] = b1.getL();
	sqL[1] = b2.getL();
	sqL[2] = k1.getL();
	sqL[3] = k2.getL();
	sort(sqL, sqL+4);

	// compare each position of L
	// starting from biggest to smallest shell
	for(int i=3; i>=0; i--) {

		// are we reaching the non-shell which 
		// has the L < 0? In this case we go to
		// fine comparison
		if (selfL[i] < 0 || sqL[i] < 0) break;

		// now compare the L
		int L1 = selfL[i];
		int L2 = sqL[i];
		if (L1 < L2 ) {
			return false;
		}else if (L1 > L2) {
			return true;
		}
	}

	// finally, go to the fine comparison
	// comparing from bra1 until ket1
	// omit these null shell case
	for(int i=0; i<4; i++) {

		// get position
		int pos = BRA1;
		if (i == 1) {
			pos = BRA2;
		} else if (i == 2) {
			pos = KET1;
		} else if (i == 3) {
			pos = KET2;
		}

		// get the corresponding shell
		const Shell& s1 = getShell(pos);
		const Shell& s2 = sq.getShell(pos);

		// is this null shell?
		if(s1.isnull() || s2.isnull()) continue;

		// now compare
		int L1 = s1.getL();
		int L2 = s2.getL();
		if (L1 < L2 ) {
			return false;
		}else if (L1 > L2) {
			return true;
		}
	}

	// we should not get here
	// if we arrive here, then something must be wrong
	cout << getName() << endl;
	cout << sq.getName() << endl;
	crash(true, "Something wrong in operator < of shell quartet class!" );
	return false;
}

const Shell& ShellQuartet::getShell(const int& pos) const {
	if (pos == BRA1) {
		return bra1;
	}else if (pos == BRA2) {
		return bra2;
	}else if (pos == KET1) {
		return ket1;
	}else {
		crash(pos != KET2, "Wrong of shell passed in getShell in shell quartet");
		return ket2;
	}
}

bool ShellQuartet::operator==(const ShellQuartet& sq) const {
	if(O != sq.O || bra1 != sq.bra1 || bra2 != sq.bra2 ||
			ket1 != sq.ket1 || ket2 != sq.ket2 || mvalue != sq.mvalue 
			|| division != sq.division) { 
		return false;
	}

	// compare deriv infor
	if (firstDerivPos != sq.firstDerivPos || secondDerivPos != sq.secondDerivPos 
			|| firstDerivDir != sq.firstDerivDir || secondDerivDir != sq.secondDerivDir) {
		return false;
	}

	// if this shell quartet has different exp factor length 
	// with input one, of course they are not same
	if (expFacListLen != sq.expFacListLen) return false;

	// now we need to compare the expotential factors
	if (expFacListLen>0) {
		const int* expList = sq.getExpFacList();
		for(int i=0; i<expFacListLen; i++) {
			if (expFacList[i] != expList[i]) return false;
		}
	}
	return true;
}

string ShellQuartet::getName() const {
	string name = "SQ";
	string Oname = getOperStringName(O);
	name = name + "_" + Oname;
	name = name + "_" + bra1.getName();
	if (!bra2.isnull())
		name = name + "_" + bra2.getName();
	if (!ket1.isnull())
		name = name + "_" + ket1.getName();
	if (!ket2.isnull())
		name = name + "_" + ket2.getName();
	if (mvalue > 0)
		name = name + "_" + "M" + lexical_cast<string>(mvalue);
	if (division >= 0)
		name = name + "_" + "C" + lexical_cast<string>(division);

	// now let's add in deriv information
	int jobOrder = getDerivJobOrder();
	for(int i=1; i<=jobOrder; i++) {
		int pos = firstDerivPos; 
		if (i == 2) pos = secondDerivPos;
		if (pos > 0) {
			name = name + "_d";
			if (pos == BRA1) {
				name = name + "a";
			}else if (pos == BRA2) {
				name = name + "b";
			}else if (pos == KET1) {
				name = name + "c";
			}else if (pos == KET2) {
				name = name + "d";
			}else{
				crash(true,"in the getName of shellquartet class, when pos is >0 but value is invalid?");
			}
		}
		int dir = firstDerivDir;
		if (i == 2) dir = secondDerivDir;
		if(dir>0) {
			if (dir == DERIV_X) {
				name = name + "x";
			}else if (dir == DERIV_Y) {
				name = name + "y";
			}else if (dir == DERIV_Z) {
				name = name + "z";
			}else{
				crash(true,"in the getName of shellquartet class, when dir is > 0 but value is invalid?");
			}
		}
	}

	// exponential factors
	if (expFacListLen>0) {
		for(int i=0; i<expFacListLen; i++) {
			if (i == 0) {
				name = name + "_";
			}
			int pos = expFacList[i]; 
			if (pos > 0) {
				if (pos == BRA1) {
					name = name + "a";
				}else if (pos == BRA2) {
					name = name + "b";
				}else if (pos == KET1) {
					name = name + "c";
				}else if (pos == KET2) {
					name = name + "d";
				}
			}
		}
	}
	return name;
}

string ShellQuartet::formArrayName(const int& rrType) const
{
	string vecName = getName();

	// here we double check, for the S type of shell quartet;
	// we should not call this function
	if (isSTypeSQ()) {
		crash(true,"something wrong in ShellQuartet::formArrayName, it's a bottom shell quartet!!");
	}

	// do we need to add in the modifier of _vrr?
	// secondly, it's only for module result as well as without 
	// division modifier or exponential factor modifier
	if (isValidVRRJob(rrType)) {
			vecName = vecName + "_vrr_array";
	}

	return vecName;
}

void ShellQuartet::getIntegralList(set<int>& list) const {
	int nTotalInts = getNInts();
	for(int i=0; i<nTotalInts; i++) {
		list.insert(i);
	}
}

bool ShellQuartet::isSTypeSQ() const {
	//
	// the S Type of sq is related to the RR processing.
	// This is used to determine whether this is bottom
	// sq. We note, that as long as the position is not 
	// a RR expansion position(see the comments in 
	// getNNonSShell function), the shell in this position
	// should not be considered in judging the S type 
	// of SQ. Therefore, we introduce the operator 
	// order to further consider the case here 
	//
	// additionally, for the derivative case; any shell quartet
	// even though it's S; it's not bottom integral
	// thereofore, for derivatives and three body ki case 
	// it's just no S type SQ
	//
	
	// consider it's derivative order
	if (getDerivJobOrder() > 0) return false;

	// consider the three body KI case
	if (O == THREEBODYKI) return false;

	// now let's consider the RR case
	// here we need to consider the RR order
	// for the mom bottom integrals, we can
	// consider the (00|x) type as derived from (00|0)
	// therefore the (00|0) type integrals is not 
	// bottom integrals
	//
	// so we use RR order instead
	int L = 0;
	int nBody = getRROrder(O);
	if (! bra1.isnull()) L += bra1.getL();
	if (! bra2.isnull() && nBody>=2) L += bra2.getL();
	if (! ket1.isnull() && nBody>=3) L += ket1.getL();
	if (! ket2.isnull() && nBody==4) L += ket2.getL();
	if (L == 0) return true;
	return false;
}

/*
 * comment out this function, because it's not used
bool ShellQuartet::hasThisIntegral(const Integral& I) const {
	if (I.getMValue() != getM()) return false;
	const Basis& b1 = I.getBasis(BRA1);
	const Basis& b2 = I.getBasis(BRA2);
	const Basis& k1 = I.getBasis(KET1);
	const Basis& k2 = I.getBasis(KET2);
	if (!bra1.isnull() && !bra1.hasThisBasisSet(b1)) return false;
	if (!bra2.isnull() && !bra2.hasThisBasisSet(b2)) return false;
	if (!ket1.isnull() && !ket1.hasThisBasisSet(k1)) return false;
	if (!ket2.isnull() && !ket2.hasThisBasisSet(k2)) return false;
	if (I.getOper()!= getOper()) return false;
	return true;
}
*/

bool ShellQuartet::canDoHRR(const int& side) const {

	// if the passing parameter is a given specific postion
	// we will still turn it into the general side for check
	int rrSide = side;
	if (side == BRA1 || side == BRA2) rrSide = BRA;
	if (side == KET1 || side == KET2) rrSide = KET;

	if (rrSide == BRA) {
		int L1 = bra1.getL();
		int L2 = bra2.getL();
		
		// consider one body integral case
		if (L2 == NULL_POS) {
			return false;
		}else if (L1 == 0 || L2 == 0) {
			return false;
		}
	}else if (rrSide == KET) {
		int L1 = ket1.getL();
		int L2 = ket2.getL();

		// consider the case whether we do not
		// have ket side, or ket side only one body
		if (L1 == NULL_POS && L2 == NULL_POS) {
			return false;
		}else if (L1 >= 0 && L2 == NULL_POS) {
			return false;
		}else if (L1 == 0 || L2 == 0) {
			return false;
		}
	}else{
		cout << "given side " << rrSide << endl;
		crash(true,"Illegal side information given in canDoHRR in shell quartet class");
	}
	return true;
}

/*
 * currently this function is not used
int ShellQuartet::getMaxL() const
{
	int L1 = bra1.getL();
	int L2 = bra2.getL();
	int L3 = ket1.getL();
	int L4 = ket2.getL();

	// get the maxL
	int maxL = L1;
	if (maxL < L2) maxL = L2;
	if (maxL < L3) maxL = L3;
	if (maxL < L4) maxL = L4;
	return maxL;

}
*/

int ShellQuartet::getNInts() const
{
	int n = 1;
	if (! bra1.isnull()) n *= bra1.getBasisSetNumber(); 
	if (! bra2.isnull()) n *= bra2.getBasisSetNumber();
	if (! ket1.isnull()) n *= ket1.getBasisSetNumber();
	if (! ket2.isnull()) n *= ket2.getBasisSetNumber();
	return n;
}

int ShellQuartet::getNInts(const int& side) const
{
	int n = 1;
	if (side == BRA) {
		if (! bra1.isnull()) n *= bra1.getBasisSetNumber(); 
		if (! bra2.isnull()) n *= bra2.getBasisSetNumber();
	}else{
		if (! ket1.isnull()) n *= ket1.getBasisSetNumber();
		if (! ket2.isnull()) n *= ket2.getBasisSetNumber();
	}
	return n;
}

/*
int ShellQuartet::getNShell() const
{
	int n = 0;
	if (! bra1.isnull()) n++; 
	if (! bra2.isnull()) n++;
	if (! ket1.isnull()) n++;
	if (! ket2.isnull()) n++;
	return n;
}
*/

int ShellQuartet::getNNonSShell() const 
{
	//
	// Here we note, that the non-S shell
	// is closely related to where we 
	// can expand the RR. 
	// for some cases, the moment integrals; 
	// the ket1 appears in RR expansion
	// but it's never to be a position
	// for RR expansion. Therefore,
	// we do not count in ket1 for MOM
	// integrals. Therefore we use 
	// operator order to express it
	// clearly
	//
	int n = 0; 
	int nBody = getOperOrder(O);

	// bra1 is always considered
	int L1 = bra1.getL();
	if (L1>0) n++;

	// consider the bra2
	if (nBody>=2) {
		int L2 = bra2.getL();
		if (L2>0) n++; 
	}

	// ket1
	if (nBody>=3) {
		int L3 = ket1.getL();
		if (L3>0) n++; 
	}

	// finally it's ket2
	if (nBody==4) {
		int L4 = ket2.getL();
		if (L4>0) n++; 
	}
	return n;
}

int ShellQuartet::getNonSShellPos() const 
{
	//
	// this function is used to return 
	// the possible RR expansion position.
	// Since the return sequence is 
	// from bra1 to ket2, therefore
	// it's currently safe.
	//
	int L1  = bra1.getL();
	int L2  = bra2.getL();
	int L3  = ket1.getL();
	int L4  = ket2.getL();
	if (L1>0) return BRA1;
	if (L2>0) return BRA2; 
	if (L3>0) return KET1; 
	if (L4>0) return KET2; 
	return NULL_POS;
}

bool ShellQuartet::isNonSShell(const int& pos) const
{
	const Shell& s = getShell(pos);
	int L = s.getL();
	if (L>0) return true;
	return false;
}

int ShellQuartet::getLSum() const
{
	int n = 0; 
	int L1 = bra1.getL();
	int L2 = bra2.getL();
	int L3 = ket1.getL();
	int L4 = ket2.getL();
	if (L1>0) n += L1;
	if (L2>0) n += L2; 
	if (L3>0) n += L3; 
	if (L4>0) n += L4; 
	return n;
}

int ShellQuartet::getLSum(const int& side) const
{
	int n = 0; 
	int L1 = bra1.getL();
	int L2 = bra2.getL();
	int L3 = ket1.getL();
	int L4 = ket2.getL();
	if (side == BRA) {
		if (L1>0) n += L1;
		if (L2>0) n += L2; 
	}else{
		if (L3>0) n += L3; 
		if (L4>0) n += L4; 
	}
	return n;
}

bool ShellQuartet::matchLM(const int& l, const int& m) const 
{
	int totalL = getLSum();
	if (totalL == l && mvalue == m) return true;
	return false;
}

bool ShellQuartet::matchL(const int& l) const 
{
	int totalL = getLSum();
	if (totalL == l) return true;
	return false;
}

bool ShellQuartet::matchLOper(const int& l, const int& oper) const 
{
	int totalL = getLSum();
	if (totalL == l && O == oper) return true;
	return false;
}

bool ShellQuartet::isnull() const {
	if (isNullOper(O)) return true;
	int nBody = getOperOrder(O);
	if (nBody == 1 && bra1.isnull()) return true;
	if (nBody == 2 && (bra1.isnull() || bra2.isnull())) return true;
	if (nBody == 3 && (bra1.isnull() || bra2.isnull() || ket1.isnull())) return true;
	if (nBody == 4 && (bra1.isnull() || bra2.isnull() || ket1.isnull() || ket2.isnull())) {
		return true;
	}
	return false;
}

bool ShellQuartet::canDODirectRRPosSearch() const {

	// 
	// for VRR, some times we can determine the RR posistion
	// by the shell quartet itself. In this function we will
	// try it.
	//
	// it's only the non-S shell are 1 or 2 than we 
	// do direct parse
	//
	int nNonSShell = getNNonSShell();
	if (nNonSShell > 2) return false;

	// if non-S shell is 1, centainly no problem
	// is non-S shell is 2, whether the two positions are both P?
	// we note that the position should be only RR possible position
	// therefore we use the operator order to count the number
	bool canDOIt = false;
	if (nNonSShell == 1) canDOIt = true;
	if (nNonSShell == 2) {
		int nPPos = 0;
		int order = getOperOrder(O);
		const Shell& b1 = getShell(BRA1);
		if (b1.getL() == 1 && order>=1) nPPos++;
		const Shell& b2 = getShell(BRA2);
		if (b2.getL() == 1 && order>=2) nPPos++;
		const Shell& k1 = getShell(KET1);
		if (k1.getL() == 1 && order>=3) nPPos++;
		const Shell& k2 = getShell(KET2);
		if (k2.getL() == 1 && order>=4) nPPos++;
		if (nPPos == 2) canDOIt = true;
	}
	return canDOIt;
}

int ShellQuartet::getRRPos() const {

	// firstly test that we can do direct parse
	bool cando = canDODirectRRPosSearch();
	if (! cando) {
		crash(true, "In getRRPos we can not do direct parse for RR position in shell quartet");
	}

	// now let's see each situation
	// for non-S shell is 1
	// we simply return this shell
	int nNonSShell = getNNonSShell();
	if (nNonSShell == 1) return getNonSShellPos(); 

	// now the non-S shell is 2 and both of them are P shell
	// get the operator
	// here we need the operator order rather than rr order
	// since we are searching for the RR position
	int order = getOperOrder(O);
	if (order == 2) {

		// in this case, the shell quartet must be <P|O|P>
		// we expand the BRA2 position
		return BRA2;

	} else if (order == 3) {

		// in this case, the shell quartet must be <PP|O|S> or <PS|O|P> or <SP|O|P>
		const Shell& b1 = getShell(BRA2);
		const Shell& b2 = getShell(KET1);
		if (b1.getL() == 1) {
			return BRA2;
		}else if (b2.getL() == 1) {
			return KET1;
		}
		return BRA1;

	} else if (order == 4) {

		// in this case, we just try to find the shell who gives the P
		for(int i=0; i<4; i++) {

			// get the possible position
			int pos = BRA1;
			if (i == 1) {
				pos = BRA2;
			} else if (i == 2) {
				pos = KET1;
			} else if (i == 3) {
				pos = KET2;
			}
			const Shell& b = getShell(pos);
			if (b.getL() == 1) {
				return pos;
			}
		}
	}

	// in default, we return the NULL_POS
	return NULL_POS;
}

void ShellQuartet::addExpFac(int pos) 
{
	// first step, we need double check
	crash(expFacListLen>=MAX_EXP_FAC_LIST || expFacListLen<0, 
			"invalid exp fac length in ShellQuartet class, already >=MAX_EXP_FAC_LIST or < 0");
	if (pos != BRA1 && pos != BRA2 && pos != KET1 && pos != KET2) {
		crash(true, "invalid position value pass in addExpFac in ShellQuartet class");
	}

	// now add in value
	expFacList[expFacListLen] = pos; 
	expFacListLen++;

	// we need to re-shuffle the value so that to make it 
	// in order
	int list1[MAX_EXP_FAC_LIST];  
	for(int i=0; i<MAX_EXP_FAC_LIST; i++) list1[i] = expFacList[i];
	std::sort(list1, list1+MAX_EXP_FAC_LIST);

	// reset the expFacList
	for(int i=0; i<MAX_EXP_FAC_LIST; i++) expFacList[i] = NULL_POS;

	// now let's copy it back
	// just omit the NULL values
	int j=0;
	for(int i=0; i<MAX_EXP_FAC_LIST; i++) {
		if (list1[i] == NULL_POS) continue;
		expFacList[j] = list1[i];
		j++;
	}
}

string ShellQuartet::getExpFacMultiplers() const
{
	string multiplier;
	for(int i=0; i<MAX_EXP_FAC_LIST; i++) {

		// omit all of NULL values
		if (expFacList[i] == NULL_POS) continue;

		// add in * symbol for all variables except the first one
		if (multiplier.size() != 0) multiplier = multiplier + "*"; 

		// now add in the variables
		int pos = expFacList[i]; 
		if (pos == BRA1) {
			multiplier = multiplier + "alpha";
		}else if (pos == BRA2) {
			multiplier = multiplier + "beta";
		}else if (pos == KET1) {
			multiplier = multiplier + "gamma";
		}else if (pos == KET2) {
			multiplier = multiplier + "delta";
		}
	}
	return multiplier;
}
