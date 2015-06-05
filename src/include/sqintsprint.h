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
 * \file    sqintsprint.h
 * \brief   for the given SQInts project, we do printing work here
 * \author  Fenglai Liu
 */
#ifndef SQINTSPRINT_H
#define SQINTSPRINT_H
#include "general.h"
#include "sqintsinfor.h"
#include "shellquartet.h"
using namespace sqintsinfor;
using namespace shellquartet;

namespace sqintsprint {

	/**
	 * print out the codes directed by sqints
	 * We do not print out RR codes here. RR codes
	 * generation will be directed by SQInts itself,
	 * here we only print out other codes
	 */
	class SQIntsPrint {

		private:

			int rrType;                        ///< RR type for this printing VRR/HRR?
			SQIntsInfor infor;                 ///< information of the sqints project
			vector<ShellQuartet> rrSQList;     ///< convert the input sq into RR input sq list

			//////////////////////////////////////////////////////////
			//           #### internal member function:             //
			//////////////////////////////////////////////////////////

			///
			/// print out the bottom MOM integrals
			///
			void printMOMBottomIntegrals(const string& file) const;

			//////////////////////////////////////////////////////////
			//  #### internal member function: VRR Head             //
			//  we note, that if input sq are all bottom ones,      //
			//  we will print them here too                         //
			//////////////////////////////////////////////////////////

			///
			/// compute the total number of integrals for the 
			/// rrSQList
			///
			/// Please remember, the shell quartets inside 
			/// the rrSQList is depending on the RR. For example,
			/// the VRR of input shell quartet (D|SP) will generate
			/// rrSQList as (D|S)s, (F|S)p, (D|S)p. The integrals
			/// here is different from the (D|SP)
			///
			int nTotalInts() const;

			///
			/// whether you have VRR expansion on the given 
			/// variable position determined by input rr sq list? 
			/// \param var    the testing variable position
			///
			bool hasVRROnVar(const int& var) const;

			///
			/// from the input rr sq list, we will see what is 
			/// the max L sum
			///
			int getMaxLSum() const;

			///
			/// generate the fmt function code for bottom SSSS integrals
			///
			void fmtIntegralsGeneration(const int& maxLSum, 
					const int& oper, const int& nSpace, ofstream& file) const;

			///
			/// perform significant integral testing for integrals with fmt function
			///
			void fmtIntegralsTest(const int& maxLSum, 
					const int& oper, const int& nSpace, ofstream& file) const;

			///
			/// get the bottom integral name
			/// these bottom integrals use the fmt function
			///
			string getBottomIntName(const int& m, const int& oper) const;

			///
			/// print out VRR result->could be in vector format or 
			/// pure variable format
			///
			void vrrResultStatement(ofstream& myfile) const; 

			///
			/// print two body overlap integral's head
			///
			void printTwoBodyOverlapHead(ofstream& file) const;

			///
			/// print three body overlap integral's head
			///
			void printThreeBodyOverlapHead(ofstream& file) const;

			///
			/// print three body kinetic integral's head
			///
			void printThreeBodyKIHead(ofstream& file) const;

			///
			/// print the moment integral's head part
			///
			void printMOMHead(ofstream& file) const;

			///
			/// print ERI integral's head
			///
			void printERIHead(ofstream& file) const;

			///
			/// print NAI integral's head
			///
			void printNAIHead(ofstream& file) const;

			///
			/// print ESP integral's head
			///
			void printESPHead(ofstream& file) const;

			///
			/// print kinetic integral's head
			///
			void printKineticHead(ofstream& file) const;

		public:

			/**
			 * constructor - just make an copy of input
			 */
			SQIntsPrint(const int& rrType0, const SQIntsInfor& infor0, 
					const vector<ShellQuartet>& inputRRSQList):rrType(rrType0),
			infor(infor0),rrSQList(inputRRSQList) { };

			/**
			 * destructor
			 */
			~SQIntsPrint() { };

			///
			/// some times the Bottom integrals are complicated
			/// here we may need to generate them explicitly to
			/// the given file
			///
			void printBottomIntegrals(const string& name) const;

			///
			/// print out the ABX CDY etc. basic variables for the 
			/// HRR process
			///
			void printHRRSideVar(const int& side, const string& fileName) const;

			///
			/// generate the contraction part of codes for 
			/// VRR process. Only used for composite shell quartets
			/// for pure shell quartet, contraction should be 
			/// automatically done
			///
			void vrrContraction(const string& myfile) const;

			///
			/// print the braket closure for VRR
			///
			void printVRREnd(ofstream& file) const;

			///
			/// print the head for VRR part
			///
			void printVRRHead(const string& name) const;

			///
			/// return the argument list for the VRR part cpp file
			/// \param sqlist: VRR part of output shell quartet list
			///
			string getVRRArgList() const;

			///
			/// return the argument list for the HRR part cpp file
			///
			string getHRRArgList() const;

			///
			/// transform the function arg list get from getVRRArgList/getHRRArgList
			/// into the arg list that we can directly use in calling the function
			/// in the code
			///
			string transformArgList(string arg) const;

	};

}


#endif

