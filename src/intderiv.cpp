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
#include "derivinfor.h"
#include "intderiv.h"
using namespace derivinfor;
using namespace intderiv;

IntDeriv::IntDeriv(const int& rrType0, const int& jobOrder0, 
		const vector<ShellQuartet>& sqlist):rrType(rrType0),
	jobOrder(jobOrder0),inputSQList(sqlist)
{
	// doing check first
	// this should be only for VRR
	if (rrType == HRR) {
		crash(true,"Currently IntDeriv class does not support HRR process");
	}

	// get the derivatives order
	// in default we stick to the deriv for x y and z
	// we may need to add in more information in the future
	int nDeriv = 3;
	int derivOrder = 1;

	// fetch the derivatives order array
	vector<int> derivOrderArray(nDeriv);
	formDerivOrderArray(derivOrder,derivOrderArray);

	// now create rrints array
	for(int iDeriv=0; iDeriv<nDeriv; iDeriv++) {

		// get the position
		int pos = derivOrderArray[iDeriv];

		// now loop over the input sq
		// we note, since the integral derivatives should be the 
		// bottom(result) part of VRR/HRR
		// therefore the unsolved integral list should be 
		// the full integral list of the corresponding sq
		for(int iSQ=0; iSQ<(int)inputSQList.size(); iSQ++) {

			// sq infor.
			const ShellQuartet& sq = inputSQList[iSQ];
			set<int> intList;
			sq.getIntegralList(intList);
			RRSQ rrsq(rrType,pos,sq,intList,jobOrder);

			// we note, that this rrsq is unique 
			// so directly push in
			rrsqList.push_back(rrsq);
		}
	}
}

void IntDeriv::updateRRSQList(list<RRSQ>& otherRRSQList) const
{
	for(list<RRSQ>::const_iterator it=rrsqList.begin(); it!=rrsqList.end(); ++it) {

		// firstly we should have a check
		// simply here we should perform a combination
		// rather than a merge work
		// therefore, any rrsq inside this class should not
		// appear in the otherRRSQList
		if (otherRRSQList.size() > 0) {
			for(list<RRSQ>::const_iterator it2=otherRRSQList.begin(); 
					it2!=otherRRSQList.end(); ++it2) {
				if (*it == *it2) {
					cout << "There are overlapping rrsq between input and this rrsqlist" << endl;
					crash(true,"Something wrong in the updateRRSQList in intderiv class");
				}
			}
		}

		// now let's do a simple adding 
		otherRRSQList.push_back(*it);
	}
}

void IntDeriv::getUnsolvedSQIntList(vector<ShellQuartet>& sqlist,
		vector<set<int> >& unsolvedIntList) const
{
	for(list<RRSQ>::const_iterator it=rrsqList.begin(); it!=rrsqList.end(); ++it) {

		// now let's go to see the RHS part
		for(int item=0; item<it->getNItems(); item++) {

			// get the corresponding sq
			const ShellQuartet& rhsSQ = it->getRHSSQ(item);

			// is it a bottom SQ?
			// we note that HRR does not apply here
			if (rhsSQ.isSTypeSQ()) continue;

			// obtain unsolved integral list
			set<int> newIntList;
			it->getUnsolvedIntList(item,newIntList);

			// now we are ready to push the raw results
			// get the position of rhsSQ in the result list
			int pos = -1;
			for(int i=0; i<(int)sqlist.size(); i++) {
				if (sqlist[i] == rhsSQ) {
					pos = i;
					break;
				}
			}

			// if it's totally a new one, we just add it
			// else we will insert the new unsolved 
			// integrals into the existing list by merging 
			// them together
			if (pos == -1) {
				sqlist.push_back(rhsSQ);
				unsolvedIntList.push_back(newIntList);
			}else{
				set<int>& oldUnsolvedList = unsolvedIntList[pos];
				for(set<int>::iterator it2=newIntList.begin(); 
						it2 != newIntList.end(); ++it2) {
					int val = *it2;
					set<int>::iterator it3=oldUnsolvedList.find(val);
					if (it3 == oldUnsolvedList.end()) {
						oldUnsolvedList.insert(val);
					}
				}
			}
		}
	}
}

