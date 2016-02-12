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
	// for non-RR, if sub file used then all of LHS are function input/output
	// however, the LHS are may be the final results, and will be changed later
	//
	// for VRR and HRR, the default LHS sq is local to function
	const ShellQuartet& lhsSQ = rrsq.getLHSSQ();
	vector<ShellQuartet>::const_iterator it = find(LHSSQList.begin(),LHSSQList.end(),lhsSQ);
	if (it == LHSSQList.end()) {
		LHSSQList.push_back(lhsSQ);
		if (moduleName == DERIV) {
			LHSSQStatus.push_back(GLOBAL_RESULT_SQ);
		}else if (moduleName == NON_RR) {
			LHSSQStatus.push_back(FUNC_INOUT_SQ);
		}else{
			LHSSQStatus.push_back(FUNC_LOCAL_SQ);
		}
		const list<int>& lhsIndex = rrsq.getLHSIndexArray();
		LHSSQIntNum.push_back(lhsIndex.size()); 
	}

	// now let's check the RHS
	// for non-RR and derivatives, all of input shell quartet are function input
	// the default LHS shell quartet status is local to function
	// we will update that later
	for(int item=0; item<it->getNItems(); item++) {
		const ShellQuartet& rhsSQ = rrsq.getRHSSQ(item);
		vector<ShellQuartet>::const_iterator it2 = find(RHSSQList.begin(),RHSSQList.end(),rhsSQ);
		if (it2 == RHSSQList.end()) {
			RHSSQList.push_back(rhsSQ);
			if (moduleName == DERIV || moduleName == NON_RR) {
				RHSSQStatus.push_back(FUNC_INOUT_SQ);
			}else{
				if (moduleName == VRR) {
					if (rhsSQ.isSTypeSQ()) {
						RHSSQStatus.push_back(BOTTOM_SQ);
					}else{
						RHSSQStatus.push_back(FUNC_LOCAL_SQ);
					}
				}else{
					RHSSQStatus.push_back(FUNC_LOCAL_SQ);
				}
			}
		}
	}
}

int SubFileRecord::getLHSIntNum(const ShellQuartet& sq) const
{
	// get the sq position
	int pos = -1;
	for(int iSQ=0; iSQ<(int)LHSSQList.size(); iSQ++) {
		const ShellQuartet& sq2 = LHSSQList[iSQ];
		if (sq2 == sq) {
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

void SubFileRecord::getMValueLimit(int& lowerM, int& upperM) const
{
	int m1 = 1000; // lower one
	int m2 = 0;
	for(int iSQ=0; iSQ<(int)LHSSQList.size(); iSQ++) {
		const ShellQuartet& sq = LHSSQList[iSQ];
		int M = sq.getM();
		if (M<m1) m1 = M;
		if (M>m2) m2 = M;
	}

	// finally, since the RHS will requires the upper M value + 1
	// we do it here
	m2 += 1;
	lowerM = m1;
	upperM = m2;
}

void SubFileRecord::updateModuleOutput(const SQIntsInfor& infor,
		const vector<ShellQuartet>& sqlist)
{
	for(int iSQ=0; iSQ<(int)LHSSQList.size(); iSQ++) {

		// we only change the status of local sq
		// if it's not local, then nothing need to be changed
		if (LHSSQStatus[iSQ] != FUNC_LOCAL_SQ) continue;

		// let's see whether this is the module result?
		const ShellQuartet& sq = LHSSQList[iSQ];
		vector<ShellQuartet>::const_iterator it = find(sqlist.begin(),sqlist.end(),sq);
		if (it == sqlist.end()) continue;

		// now form the status
		// if file in split, the module output must be output for function, too
		LHSSQStatus[iSQ] = FUNC_INOUT_SQ;
		if (infor.isResult(sq)) {
			LHSSQStatus[iSQ] = GLOBAL_RESULT_SQ;
		}
	}
}

void SubFileRecord::updateModuleInput(const vector<ShellQuartet>& sqlist)
{
	for(int iSQ=0; iSQ<(int)RHSSQList.size(); iSQ++) {

		// we only change the status of local sq
		// if it's not local, then nothing need to be changed
		if (RHSSQStatus[iSQ] != FUNC_LOCAL_SQ) continue;

		// let's see whether this is the module result?
		const ShellQuartet& sq = RHSSQList[iSQ];
		vector<ShellQuartet>::const_iterator it = find(sqlist.begin(),sqlist.end(),sq);
		if (it == sqlist.end()) continue;

		// now form the status
		// if file in split, the module input must be input for function, too
		RHSSQStatus[iSQ] = FUNC_INOUT_SQ;
	}
}

void SubFileRecord::updateOutput(const SubFileRecord& record)
{
	const vector<ShellQuartet>& rhs = record.getRHSSQList();
	for(int iSQ=0; iSQ<(int)LHSSQList.size(); iSQ++) {

		// if this shell quartet status is already determined;
		// we move on
		const ShellQuartet& sq = LHSSQList[iSQ];
		if (LHSSQStatus[iSQ] != FUNC_LOCAL_SQ) continue;

		// this is the sq need to pass to next sub file record
		vector<ShellQuartet>::const_iterator it = find(rhs.begin(),rhs.end(),sq);
		if (it != rhs.end()) {
			LHSSQStatus[iSQ] = FUNC_INOUT_SQ;
		}
	}
}

void SubFileRecord::updateInput(const SubFileRecord& record)
{
	const vector<ShellQuartet>& lhs = record.getLHSSQList();
	for(int iSQ=0; iSQ<(int)RHSSQList.size(); iSQ++) {

		// if this shell quartet status is already determined;
		// we move on
		const ShellQuartet& sq = RHSSQList[iSQ];
		if (RHSSQStatus[iSQ] != FUNC_LOCAL_SQ) continue;

		// this is the sq generated from previous sub file
		vector<ShellQuartet>::const_iterator it = find(lhs.begin(),lhs.end(),sq);
		if (it != lhs.end()) {
			RHSSQStatus[iSQ] = FUNC_INOUT_SQ;
		}
	}
}
