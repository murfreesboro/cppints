/**
 *
 * CPPINTS: A C++ Program to Generate Analytical Integrals Based on Gaussian
 * Form Primitive Functions
 *
 * Copyright (C) 2015 The State University of New York at Buffalo
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
			/// this is the working function to generate the integral codes
			/// this function will update infor in running time
			/// therefore it's not constant
			///
			void intCodeGeneration();

			///
			/// this function is to append the function prototype file to the
			/// main cpp file head
			///
			void appendFuncPrototype(int moduleName, ofstream& CPP) const;

			///
			/// form the function call in the main cpp file
			///
			void formFunctionCall(int moduleName, ofstream& CPP) const;

			///
			/// this function is used to form the "work" cpp file which is called
			/// by the main cpp file. For example, the VRR, HRR1 etc. because
			/// it's too large so it has its own cpp file to perform the work.
			///
			/// this function is to form this kind of work file
			///
			void formWorkFile(int moduleName) const;

			///
			/// this is to append the file to the main body file (represented by CPP)
			/// 
			void appendFile(int moduleName, ofstream& CPP, bool checkFile = true) const;

			///
			/// assemble the CPP file from all of pieces already 
			/// in the temp working dir
			///
			void assembleCPPFiles() const;

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
					const int& oper):infor(oper,infor0,bra1,NULL_POS,NULL_POS,NULL_POS) { };

			/**
			 * constructor for the 2 body shell quartet
			 * only bra1 and bra2 exist
			 */
			SQInts(const Infor& infor0, const int& bra1, 
					const int& bra2, const int& oper):infor(oper,infor0,bra1,bra2,NULL_POS,NULL_POS) { };

			/**
			 * constructor for the 3 body shell quartet
			 * bra1, bra2 and ket1 exist
			 */
			SQInts(const Infor& infor0, const int& bra1, 
					const int& bra2, const int& ket1, 
					const int& oper):infor(oper,infor0,bra1,bra2,ket1,NULL_POS) { };

			/**
			 * constructor for the 4 body shell quartet
			 */
			SQInts(const Infor& infor0, const int& bra1, 
					const int& bra2, const int& ket1, const int& ket2, 
					const int& oper):infor(oper,infor0,bra1,bra2,ket1,ket2) { };

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
			void codeGeneration();

	};

}


#endif

