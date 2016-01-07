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
 * \file    vrrinfor.h
 * \brief   utility class to help vrr do contraction, declare work etc.
 * \author  Fenglai Liu
 */
#ifndef VRRINFOR_H
#define VRRINFOR_H
#include "general.h"
#include "infor.h"
#include "shellquartet.h"
#include "inttype.h"
#include "boost/lexical_cast.hpp"
using namespace infor;
using namespace shellquartet;
using namespace inttype;

namespace rr {
	class RR;
}

namespace sqintsinfor {
	class SQIntsInfor;
}

namespace vrrinfor {

	using namespace rr;
	using namespace sqintsinfor;


	/**
	 * \class VRRInfor
	 */
	class VRRInfor : public Infor {

		private:

			/////////////////////////////////////////////
			//             data members                //
			/////////////////////////////////////////////

			//
			// file split etc. information
			//
			bool vrrInFileSplit;               ///< does the vrr part of code in a single file?
			bool vrrContSplit;                 ///< does the vrr and contraction part of codes split also?
			int lastSection;                   ///< which is the previous section of VRR

			// 
			// general information for RR
			//
			int oper;                          ///< operator information
			vector<ShellQuartet> vrrSQList;    ///< vrr result shell quartet list (no exponential infor etc.)
			vector<set<int> > solvedIntList;   ///< the integral list corresponding result shell quartet
			vector<ShellQuartet> outputSQList; ///< output shell quartet list (with exponential infor etc.)
			vector<set<int> > outputIntList;   ///< the integral list corresponding output shell quartet

			/////////////////////////////////////////////
			// !!! inline functions                    //
			/////////////////////////////////////////////
         bool isLastSection() const {
				if (lastSection == NULL_POS) return true;
				return false;
			};

			/////////////////////////////////////////////
			// !!!        utility functions            //
			/////////////////////////////////////////////
			
			///
			/// whether you have VRR expansion on the given 
			/// variable position for the VRR results? 
			/// \param var    the testing variable position
			///
			bool hasVRROnVar(const int& var) const;

			///
			/// determine max L Sum for the VRR results
			///
			int getMaxLSum() const;

			///
			/// print out the bottom MOM integrals
			/// \param name: the file name for printing
			/// \param inputSQList: the result shell quartet list for this cpp file, got from sqintsinfor
			///
			void printMOMBottomIntegrals(const string& name, const vector<ShellQuartet>& inputSQList) const;

			///
			/// get the bottom integral name
			/// these bottom integrals use the fmt function
			///
			string getBottomIntName(const int& m, const int& oper) const;

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

			/////////////////////////////////////////////
			// !!!   head printing functions           //
			/////////////////////////////////////////////

			///
			/// print out the VRR result statement
			///
			void printResultStatement(ofstream& myfile, const SQIntsInfor& infor) const; 

			///
			/// print two body overlap integral's head
			///
			void printTwoBodyOverlapHead(ofstream& file, const SQIntsInfor& infor) const;

			///
			/// print three body overlap integral's head
			///
			void printThreeBodyOverlapHead(ofstream& file, const SQIntsInfor& infor) const;

			///
			/// print the moment integral's head part
			///
			void printMOMHead(ofstream& file, const SQIntsInfor& infor) const;

			///
			/// print ERI integral's head
			///
			void printERIHead(ofstream& file, const SQIntsInfor& infor) const;

			///
			/// print NAI integral's head
			///
			void printNAIHead(ofstream& file, const SQIntsInfor& infor) const;

			///
			/// print ESP integral's head
			///
			void printESPHead(ofstream& file, const SQIntsInfor& infor) const;

			///
			/// print kinetic integral's head
			///
			void printKineticHead(ofstream& file, const SQIntsInfor& infor) const;

			/////////////////////////////////////////////
			// !!!   contraction functions             //
			/////////////////////////////////////////////

			///
			/// this is to do the normal VRR contraction work
			/// that is to say, VRR and contraction are performed in a same file,
			/// either with file split or non-file-split module
			///
			void normalVRRContraction(const string& filename, const SQIntsInfor& infor) const;

			///
			/// if contraction is performed in different place with the VRR code,
			/// then in VRR part we need to do "preliminary digestion" work so
			/// that to digest all of VRR results into array form. So that in 
			/// functtion vrrContractionInSplit we can form the real contraction work
			///
			void preliminaryVRRContraction(const string& filename) const;

			///
			/// perform the VRR contraction work if the contraction is in different
			/// file comparing with the VRR code
			///
			void vrrContractionInSplit(const string& filename, const SQIntsInfor& infor, 
					const vector<ShellQuartet>& sqlist) const;

			/////////////////////////////////////////////
			// !!!   argument printing                 //
			/////////////////////////////////////////////

			///
			/// this function returns the function prototype for the VRR contraction
			/// function if the contraction code is performed in several parts
			/// \param infor      the information center
			/// \param sqlist     the output shell quartets of contraction work
			/// \param fileIndex  the file index for doing contraction work, starts from 1
			///
			string getVrrContractionFuncName(const SQIntsInfor& infor, 
					const vector<ShellQuartet>& sqlist, int fileIndex) const; 

		public:

			///
			/// the deriver function to print the head of VRR
			///
			void printVRRHead(const SQIntsInfor& infor) const;

			///
			/// the deriver function to perform contraction work
			///
			void vrrContraction(const SQIntsInfor& infor) const;

			///
			/// derive the argument list if VRR part of code is in a single file
			///
			string getVRRArgList(const SQIntsInfor& infor) const;

			///
			/// constructor to form vrr information
			///
			VRRInfor(const SQIntsInfor& infor, const RR& vrr);

			///
			/// destructor
			///
			~VRRInfor() { };

			///
			/// whether VRR is in file split mode?
			///
			bool fileSplit() const { return vrrInFileSplit; };

			///
			/// whether the contraction and VRR are split
			///
			bool vrrContractionSplit() const { return vrrContSplit; };
	};

}


#endif

