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
#include<algorithm>
#include "general.h"
#include "rrints.h"
#include "subfilerecord.h"
using namespace std;
using namespace rrints;
using namespace subfilerecord;

void SubFileRecord::updateFromRRSQ(const RRSQ& rrsq)
{
	// let's update the LHS shell quartet
	//
	// for derivatives all of lhs is final result
	// 
	// for non-RR, if sub file usd then all of LHS are array form
	// however, they are may be the final results, and will be changed later
	//
	// for VRR and HRR, the default LHS sq is variable form
	const ShellQuartet& lhsSQ = rrsq.getLHSSQ();
	vector<ShellQuartet>::const_iterator it = find(LHSSQList.begin(),LHSSQList.end(),lhsSQ);
	if (it == LHSSQList.end()) {
		LHSSQList.push_back(lhsSQ);
		if (moduleName == DERIV) {
			LHSSQStatus.push_back(GLOBAL_RESULT_SQ);
			hasABCD = true;
		}else if (moduleName == NON_RR) {
			LHSSQStatus.push_back(ARRAY_SQ);
		}else{
			LHSSQStatus.push_back(VARIABLE_SQ);
		}
		const list<int>& lhsIndex = rrsq.getLHSIndexArray();
		LHSSQIntNum.push_back(lhsIndex.size()); 
	}

	// now let's check the RHS
	// for non-RR and derivatives, all of input shell quartet are module input
	// in the case of sub files used, they must be in array
	// the default LHS shell quartet status is also variable form
	// we will update that later
	for(int item=0; item<it->getNItems(); item++) {
		const ShellQuartet& rhsSQ = rrsq.getRHSSQ(item);
		vector<ShellQuartet>::const_iterator it2 = find(RHSSQList.begin(),RHSSQList.end(),rhsSQ);
		if (it2 == RHSSQList.end()) {
			RHSSQList.push_back(rhsSQ);
			if (moduleName == DERIV || moduleName == NON_RR) {
				RHSSQStatus.push_back(ARRAY_SQ);
			}else{
				RHSSQStatus.push_back(VARIABLE_SQ);
			}
		}
	}
}

int SubFileRecord::getLHSIntNum(const ShellQuartet& outputSQ) const
{
	// for VRR case, the output shell quartet pass in
	// may have exponential information etc. contained
	// we need to destroy them before searching
	ShellQuartet lhsSQ(outputSQ);
	if (moduleName == VRR) {
		lhsSQ.destroyMultipliers();
	}

	// get the sq position
	int pos = -1;
	for(int iSQ=0; iSQ<(int)LHSSQList.size(); iSQ++) {
		const ShellQuartet& sq = LHSSQList[iSQ];
		if (lhsSQ == sq) {
			pos = iSQ;
			break;
		}
	}

	// check pos
	if (pos < 0) {
		crash(true, "we did not find the input lhs sq in SubFileRecord::getLHSIntNum");
	}

	// now retrun result
	return LHSSQIntNum[pos];
}

void SubFileRecord::updateModuleOutput(const SQIntsInfor& infor,
		const vector<ShellQuartet>& sqlist)
{
	for(int iSQ=0; iSQ<(int)LHSSQList.size(); iSQ++) {

		// let's see whether this is the module result?
		const ShellQuartet& sq = LHSSQList[iSQ];
		vector<ShellQuartet>::const_iterator it = find(sqlist.begin(),sqlist.end(),sq);
		if (it == sqlist.end()) continue;

		// now form the status
		LHSSQStatus[pos] = ARRAY_SQ;
		if (infor.isResult(sq)) {
			LHSSQStatus[pos] = GLOBAL_RESULT_SQ;
			hasABCD = true;
		}else{
			outputSQList.push_back(sq); 
		}
	}
}

void SubFileRecord::updateVRROutput(bool destroyMultiplerInfor,
		const SQIntsInfor& infor, const vector<ShellQuartet>& sqlist) 
{
	// we only do for VRR module
	if (moduleName != VRR) {
		crash(true, "module is not VRR for SubFileRecord::updateVRROutput");
	}

	// now let's see the result
	for(int iSQ=0; iSQ<(int)sqlist.size(); iSQ++) {

		// form the shell quartet
		ShellQuartet newSQ(sqlist[iSQ]);
		if (destroyMultiplerInfor) {
			newSQ.destroyMultipliers();
		}

		// get the sq position
		int pos = -1;
		for(int iSQ2=0; iSQ2<(int)LHSSQList.size(); iSQ2++) {
			const ShellQuartet& sq = LHSSQList[iSQ2];
			if (newSQ == sq) {
				pos = iSQ2;
				break;
			}
		}

		// check pos
		// then we know that whether we have this shell quartet
		// as module result in LHS
		if (pos < 0) continue;

		// determine that VRR and contraction is split or not
		// if multiper information is destroyed, it means that
		// input shell quartet is VRR module result
		// so VRR and contraction are done tother
		bool vrrContSplit = true;
		if (destroyMultiplerInfor) vrrContSplit = false;

		// for VRR module, to determine the module output shell
		// quartet is a bit of complicated
		if (vrrContSplit) {

			// in this case contraction always form the final
			// result, so the module output should be in array form
			LHSSQStatus[pos] = ARRAY_SQ;

			// here it always included into outputSQList
			outputSQList.push_back(newSQ); 
		}else{

			// for contraction is contained into VRR,
			// there's a possible case that the module
			// output is global result
			LHSSQStatus[pos] = ARRAY_SQ;
			if (infor.isResult(sqlist[iSQ])) {
				LHSSQStatus[pos] = GLOBAL_RESULT_SQ;
				hasABCD = true;
			}else{
				outputSQList.push_back(sqlist[iSQ]); 
			}
		}
	}
}
