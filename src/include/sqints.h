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
 * \file    sqints.h
 * \brief   for the given shell quartet, generate it's integrals
 * \author  Fenglai Liu
 */
#ifndef SQINTS_H
#define SQINTS_H
#include "general.h"
#include "sqintsinfor.h"
using namespace sqintsinfor;

namespace sqints {

	/**
	 * \class input shell quartet for RR
	 * \brief parsing the input L code so that to generate input shell quartets for this
	 *        recursive relation process
	 *
	 */
	class SQInts {

		private:

			//
			// infor class is contained here
			// just to make the reference to infor class more 
			// convenient
			//
			SQIntsInfor infor;                 ///< information of the project

			///
			/// printing the head part of cpp file
			/// comments, license etc.
			///
			void headPrinting(ofstream& file) const;

			///
			/// return the cpp function argument list
			///
			string getArgList() const;

			///
			/// working function to perform the HRR work for the given 
			/// side
			/// \param iSide  0 is for second side, 1 is for first side
			/// \param inputList input shell quartet for HRR work
			/// \param bottomSQList  the result SQ list after HRR work
			///
			void doCoreHRR(const int& iSide, const vector<ShellQuartet>& inputList, 
					vector<ShellQuartet>& bottomSQList) const; 

			///
			/// function to perform the entire HRR work 
			/// \return hrrResultList output hrr result shell quartets
			///
			void doHRR(vector<ShellQuartet>& hrrResultSQList) const; 

			///
			/// function to perform the entire VRR work 
			/// \input hrrResultList input hrr result shell quartets
			///
			void doVRR(const vector<ShellQuartet>& hrrResultSQList) const; 

			///
			/// main body to perform the RR work for the given shell quartet
			///
			void doRR() const;

			///
			/// if the file split requirement exists, we will do the 
			/// file assembling here so that to generate the working 
			/// cpp files
			///
			void assembleWorkingCPPFile(bool isHRR) const;

			///
			/// assemble the top CPP file here
			///
			void assembleTopCPPFile() const;

		public:

			///////////////////////////////////////////////////////////////////////
			//                                                                   //
			//                    constructor and destructor                     //
			//                                                                   //
			///////////////////////////////////////////////////////////////////////

			/**
			 * constructor for the 1 body shell quartet
			 * only bra1 exists
			 */
			SQInts(const Infor& infor0, const int& bra1, 
					const int& oper, 
					const int& order):infor(order,oper,infor0,bra1,NULL_POS,NULL_POS,NULL_POS) { };

			/**
			 * constructor for the 2 body shell quartet
			 * only bra1 and bra2 exist
			 */
			SQInts(const Infor& infor0, const int& bra1, 
					const int& bra2, const int& oper, 
					const int& order):infor(order,oper,infor0,bra1,bra2,NULL_POS,NULL_POS) { };

			/**
			 * constructor for the 3 body shell quartet
			 * bra1, bra2 and ket1 exist
			 */
			SQInts(const Infor& infor0, const int& bra1, 
					const int& bra2, const int& ket1, const int& oper, 
					const int& order):infor(order,oper,infor0,bra1,bra2,ket1,NULL_POS) { };

			/**
			 * constructor for the 4 body shell quartet
			 */
			SQInts(const Infor& infor0, const int& bra1, 
					const int& bra2, const int& ket1, const int& ket2, 
					const int& oper, 
					const int& order):infor(order,oper,infor0,bra1,bra2,ket1,ket2) { };

			/**
			 * destructor
			 */
			~SQInts() { };

			/**
			 * whether the cpp file corresponding to this 
			 * sqints already exist?
			 */
			bool isFileExist() const; 

			/**
			 * code generation
			 */
			void codeGeneration() const;

	};

}


#endif

