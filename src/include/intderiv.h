/**
 *
 * CPPINTS: A C++ Program to Generate Analytical Integrals Based on Gaussian
 * Form Primitive Functions
 *
 * Copyright (C) 2012-2014 Fenglai Liu
 * This softare uses the MIT license as below:
 *
 *	Permission is hereby granted, free of charge, to any person obtaining 
 *	a copy of this software and associated documentation files (the "Software"), 
 *	to deal in the Software without restriction, including without limitation 
 *	the rights to use, copy, modify, merge, publish, distribute, sublicense, 
 *	and/or sell copies of the Software, and to permit persons to whom the Software 
 *	is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all 
 * copies or substantial portions of the Software.
 *						    
 *	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
 *	INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
 *	PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
 *	FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
 *	ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * \file    intderiv.h
 * \brief   calculate the integral derivatives explicitly
 * \author  Fenglai Liu
 */
#ifndef INTDERIV_H
#define INTDERIV_H
#include "general.h"
#include "rrints.h"
#include "shellquartet.h"
using namespace rrints;
using namespace shellquartet;

namespace intderiv {

	/**
	 * RR is responsible to form vertical RR for the given group 
	 * of shell quartets
	 */
	class IntDeriv {

		private:

			int rrType;                        ///< rr type
			int jobOrder;                      ///< deriv order
			vector<ShellQuartet> inputSQList;  ///< the input shell quartets
			list<RRSQ> rrsqList;               ///< collecting the rrsq result

		public:

			///
			/// contruct the integral derivatives class
			/// \param rrType    VRR's RR algorithm name
			/// \param jobOrder  derivatives order
			/// \param sqlist    input shell quartet list
			///
			IntDeriv(const int& rrType0, const int& jobOrder0, 
					const vector<ShellQuartet>& sqlist);

			///
			/// destructor
			///
			~IntDeriv() { };

			///
			/// update this RRSQ results to another rrsqlist
			///
			void updateRRSQList(list<RRSQ>& otherRRSQList) const;

			///
			/// from the result rrsq list, get the unsolved shell quartets (RHS)
			/// as well as it's corresponding unsolved integral list
			///
			void getUnsolvedSQIntList(vector<ShellQuartet>& sqlist,
					vector<set<int> >& unsolveIntList) const;

	};
}


#endif

