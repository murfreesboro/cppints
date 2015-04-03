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

int ShellQuartet::hrrCompare(const ShellQuartet& sq, bool isBraSide) const {

	// checking the operator first
	if (O != sq.getOper()) {
		crash(true,"Only same type of integral could be compared in hrr less than function");
	}

	// different division of shell quartet has no comparing relationship
	if (division != sq.getDivision()) {
		crash(true,"Only same class of integral could be compared in hrr less than function");
	}

	// now below it's the core algorithm for HRR recursive comparison.
	// We note, that the shell quartet who is smaller should actually
	// appear in the LHS, and the shell quartet who is larger should 
	// appear in the RHS of HRR recursive expression. Therefore, 
	// we are able to generate correct shell quartet list in RR that 
	// sq0 -> sq1 -> sq2 ...
	//
	// For the HRR, the recursive relation has a feature that the total
	// L is keeping same between the LHS and RHS term, therefore it's 
	// not easy to use L sum to compare. However, here we note a 
	// feature, that for the RHS term, one term whose maxL higher than
	// LHS term, and both of two terms on RHS whose minL lower than LHS
	// term. So we will use this feature for comparison.
	// 
	// Generally, the rule is like this(for the given bra/ket side):
	// 1 we compare minimum of bra1 and bra2(or ket1 and ket2 if on
	// ket side), if minL1 < minL2, then sq1 > sq2; if minL1 > minL2
	// then sq1 < sq2; else we go to 2;
	// 2 we compare the max angular momentum to see who is larger,
	// if maxL1 < maxL2, then sq1 < sq2; if maxL1 > maxL2,
	// then sq1 > sq2; else we go to 3;
	// 3 if we pass 1 and 2 then it's sure that we have same angular
	// momentum arrangment on bra/ket side. then we compare the BRA1
	// or KET1 side value. If the value is larger, then it's smaller;
	// else it's larger. We have a detailed explanation below.
	// 4 if we go here, then it means the bra side is totally same.
	// then we repeat the above process for the ket side
	//
	// We note, that the above algorithm is not perfect. For some 
	// given shell quartets, it can not give the correct results.
	// For example, the ERI of (DP|**), according to the HRR equation;
	// for expansion on BRA1 position it gives:
	// (DP|**) = (PD|**) - ABi*(PP|**)
	// for this case, we can give the correct result because of 
	// step 3.
	//
	// However, the similar case is that we expands (PD|**)
	// on BRA2 position. Since HRR is convertible, so it 
	// gives:
	// (PD|**) = (DP|**) + ABi*(PP|**)
	// obviously this case we will give the wrong order between
	// (PD|**) and (DP|**). Therefore, here we have a potential
	// bug for the algorithm.
	//
	// If we follow the rule that L(BRA1)>= L(BRA2), and 
	// L(KET1) >= L(KET2); usually we are safe on using this 
	// algorithm. The exception is for SPD composite shell, 
	// where our algorithm will give wrong results. Therefore,
	// our code can not be used for SPD shell case.
	//
	int minL1 = -1;
	int maxL1 = -1;
	int minL2 = -1;
	int maxL2 = -1;
	int L     = -1;
	int LSQ   = -1;

	if (isBraSide) {

		// first step, checking the bra side shell
		if (bra1.isnull() || bra2.isnull()) {
			crash(true,"If bra2/bra1 is none, why we go to the hrr compare in shell quartet?");
		}

		// get shell
		const Shell& b1 = sq.getShell(BRA1);
		const Shell& b2 = sq.getShell(BRA2);

		// now take the min and max angular momentum value
		minL1 = min(bra1.getL(),bra2.getL());
		maxL1 = max(bra1.getL(),bra2.getL());
		minL2 = min(b1.getL(),b2.getL());
		maxL2 = max(b1.getL(),b2.getL());
		L     = bra1.getL();
		LSQ   = b1.getL();

	}else{

		// first step, checking the ket side shell
		if (ket1.isnull() || ket2.isnull()) {
			crash(true,"If ket2/ket1 is none, why we go to the hrr compare in shell quartet?");
		}

		// get shell
		const Shell& k1 = sq.getShell(KET1);
		const Shell& k2 = sq.getShell(KET2);

		// now take the min and max angular momentum value
		minL1 = min(ket1.getL(),ket2.getL());
		maxL1 = max(ket1.getL(),ket2.getL());
		minL2 = min(k1.getL(),k2.getL());
		maxL2 = max(k1.getL(),k2.getL());
		L     = ket1.getL();
		LSQ   = k1.getL();
	}

	// now consider the relation
	if (minL1 < minL2) {
		return GT;
	}else if (minL1 > minL2) {
		return LT;
	}else{

		// compare the maximum L
		if (maxL1 < maxL2) {
			return LT;
		}else if (maxL1 > maxL2) {
			return GT;
		}else{

			// compare the L and LSQ
			// we do this is because we have 
			// an unwritten rule that L(BRA1) >= L(BRA2)
			// L(KET1) >= L(KET2) for the input shell quartet
			// and HRR expansion is perfered on BRA2 and KET2
			if (L > LSQ) {
				return LT;
			}else if (L < LSQ) {
				return GT;
			} else {
				return EQ;
			}
		}
	}
}

bool ShellQuartet::hrrLessThan(const ShellQuartet& sq, const int& side) const {

	bool isBraSide = true;
	if (side == KET) isBraSide = false;
	int compareResult = hrrCompare(sq, isBraSide); 
	if (compareResult == GT) {
		return false;
	}else if (compareResult == LT) {
		return true;
	}else{
		int compareResult2 = hrrCompare(sq, ! isBraSide);
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
	}
}

bool ShellQuartet::operator<(const ShellQuartet& sq) const {

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
	}else{
		return true;
	}
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
	return name;
}

string ShellQuartet::formArrayName(const int& rrType, const int& status, 
		const int& pos, bool withModifier) const
{
	string vecName = getName();

	// do we need to add in the modifier of _vrr?
	if (rrType != HRR && status == MODULE_RESULT && ! withModifier) {
		vecName = vecName + "_vrr";
	}

	// do we need to add in deriv modifier?
	if (isDerivInfor(pos)) {
		string derivInfor = symTransform(pos);
		vecName = vecName + "_" + derivInfor;
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
	int L = 0;
	int nBody = getOperOrder(O);
	if (! bra1.isnull()) L += bra1.getL();
	if (! bra2.isnull() && nBody>=2) L += bra2.getL();
	if (! ket1.isnull() && nBody>=3) L += ket1.getL();
	if (! ket2.isnull() && nBody==4) L += ket2.getL();
	if (L == 0) return true;
	return false;
}

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

bool ShellQuartet::canDoHRR(const int& side) const {

	// if the passing parameter is a given specific postion
	// we will still turn it into the general side for check
	int rrSide = side;
	if (side == BRA1 || side == BRA2) rrSide = BRA;
	if (side == KET1 || side == KET2) rrSide = KET;

	if (rrSide == BRA) {
		int L1 = bra1.getL();
		int L2 = bra2.getL();
		if (L1 == 0 || L2 == 0) return false;
	}else if (rrSide == KET) {
		int L1 = ket1.getL();
		int L2 = ket2.getL();
		if (L1 == 0 || L2 == 0) return false;
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

