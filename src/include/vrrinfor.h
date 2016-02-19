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
#include "subfilerecord.h"
using namespace infor;
using namespace shellquartet;
using namespace subfilerecord;

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
			bool vrrInFileSplit;               ///< does the vrr part of code in sub files?
			bool vrrContSplit;                 ///< does the vrr and contraction split in different parts?
			int nextSection;                   ///< which is the next section of VRR
			int oper;                          ///< operator information

			// 
			// general information for RR
			//
			vector<ShellQuartet> vrrSQList;    ///< vrr result shell quartet list (no exponential infor etc.)
			vector<set<int> > solvedIntList;   ///< the integral list corresponding result shell quartet
			vector<ShellQuartet> outputSQList; ///< output shell quartet list (with exponential infor etc.)
			vector<set<int> > outputIntList;   ///< the integral list corresponding output shell quartet
			vector<int> outputSQStatus;        ///< the output shell quartet for VRR, in array or variable form?

			//
			// sub files record
			//
			vector<SubFileRecord> subFilesList;  ///< divide the whole VRR into different sub functions

			/////////////////////////////////////////////
			// !!! inline functions                    //
			/////////////////////////////////////////////
         bool isLastSection() const {
				if (nextSection == NULL_POS) return true;
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
			/// for the operator erf(r12)/r12, set up the scaling for bottom integrals
			///
			void setupErfPrefactors(const int& maxLSum, 
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
			void printResultStatement(ofstream& myfile) const; 

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
			/// print the exp(-omega*r12^2) head (EXPR12 operator)
			///
			void printEXPR12Head(ofstream& file, const SQIntsInfor& infor) const;

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
			/// this is the working function to do VRR contraction
			///
			/// if in VRR contraction split way, this function performs contraction
			/// only with the VRR result sq array; which should be done in VRR part.
			///
			/// if VRR and contraction not split, we just simply append the contraction
			/// to the VRR file
			///
			/// file index is the index for sub files. If it's set to -1, then VRR
			/// is not split and contraction appends to the VRR bottom.
			///
			void contraction(const SQIntsInfor& infor, 
					const vector<ShellQuartet>& sqlist, const int& fileIndex) const;

			/////////////////////////////////////////////
			// !!!    forming sub files                //
			/////////////////////////////////////////////

			///
			/// for the given shell quartet, whether it's the VRR output
			/// and how many contraction lines corresponding to this sq
			///
			int contrationCount(const ShellQuartet& sq) const;

			///
			/// this is the working function to form sub file records 
			///
			void formSubFiles(bool onlyOneSubFile, const SQIntsInfor& infor, const RR& vrr);

			///
			/// this is the function to derive the function arguments for the 
			/// given sub file record
			///
			/// the input subFileIndex starting from 0
			///
			string getVRRArgList(const int& subFileIndex, 
					const SQIntsInfor& infor, const SubFileRecord& record) const;

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
			/// this is the driver function for sub files forming
			/// we note, that this is the function to form the 
			/// choice of VRRSplit, also VRRContSplit
			///
			void subFilesForming(const SQIntsInfor& infor, const RR& vrr);

			///
			/// constructor to form vrr information, initialize the 
			/// content
			///
			VRRInfor(const SQIntsInfor& infor, const RR& vrr);

			///
			/// update the VRR output shell quartet list against
			/// other modules input, the input sq in the vector 
			/// are all in array form 
			///
			void updateOutputSQInArray(const vector<ShellQuartet>& moduleInput);

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

			/// 
			/// return the output sq list
			///
			const vector<ShellQuartet>& getOutputSQList() const { return outputSQList; };

			/// 
			/// return the output sq status list
			///
			const vector<int>& getOutputSQStatus() const { return outputSQStatus; };

			///
			/// get the number of sub file records
			///
			int getNSubFiles() const { return subFilesList.size(); };

			///
			/// return the given sub file record
			///
			const SubFileRecord& getSubFileRecord(int i) const { return subFilesList[i]; };
	};

}


#endif

