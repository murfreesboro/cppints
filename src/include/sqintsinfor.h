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
 * \file    sqintsinfor.h
 * \brief   processing the input information for the whole program
 * \author  Fenglai Liu
 */
#ifndef SQINTSINFOR_H
#define SQINTSINFOR_H
#include "general.h"
#include "infor.h"
#include "shellquartet.h"
#include "inttype.h"
#include "derivinfor.h"
#include "boost/lexical_cast.hpp"
using namespace infor;
using namespace shellquartet;
using namespace inttype;
using namespace derivinfor;

namespace vrrinfor {
	class VRRInfor;
}

namespace sqintsinfor {

	using namespace vrrinfor;


	/**
	 * \class SQIntsInfor
	 *
	 * This class is the messager for SQInts. In general, SQInts is the commander
	 * and all of data in the SQInts are actually folded here so that SQInts
	 * could send out this infor class to direct the movement of working 
	 * modules.
	 *
	 * For example, SQInts use this infor class to direct the printing in
	 * all kinds of modules, like VRR, HRR etc.
	 * 
	 * note for the derivatives shell quartets:
	 *
	 * the derivatives shell quartet is formed in dimension of (nSQ,nDeriv).
	 * nSQ is in same dimension of inputSQList, and nDeriv represents 
	 * the dimension for the derivatives.
	 *
	 * such arrangement is important. In function getOffset, we need such 
	 * arrangement so that we can figure out the index for the final abcd[..]
	 * array. On the other hand, in the function for printing derivatives 
	 * in nonrr.cpp, we also rely on such arrangement. In the function
	 * getSQIndexForDerivFile and getNumDerivFiles, the number of derivatives
	 * files are also formed based on this.
	 *
	 */
	class SQIntsInfor : public Infor {

		private:

			//
			// general information for the code sections
			// we keep a record of the with-array/file-split status for HRR, nonRR etc.
			// so that we can keep a track during the code generation
			//
			// we do not have with array setting for VRR. Because it's already in variable
			// form. Also we do not have file split for VRR too, because VRR is always 
			// the last module that we do not need to keep a copy of it
			//
			bool vrrInSingleFile;              ///< does the VRR code is in a separate file?
			bool vrrContractionSplit;          ///< whether we split the vrr and contraction
			bool hrr1InSingleFile;             ///< does the HRR1 code is in a separate file?
			bool hrr2InSingleFile;             ///< does the HRR2 code is in a separate file?
			bool nonRRInSingleFile;            ///< does the non-RR code is in a separate file?
			bool derivInSingleFile;            ///< does the deriv code is in a separate file?
			bool doHRRWork;                    ///< determine that whether we do HRR work
			vector<int> sectionInfor;          ///< section sequence information

			// 
			// general information for sqints
			//
			int oper;                          ///< operator for the input shell quartet
			size_t minDerivInts;               ///< this is the smallest number in derivRecords[4]
			size_t derivRecords[4];            ///< this is used to record evaluation results in formDerivInfor
			vector<int> inputShellCodes;       ///< keep a copy of input shell codes
			vector<int> braCoeOffset;          ///< for BRA side composite shell quartet 
			                                   ///< record coe offset
			vector<int> ketCoeOffset;          ///< for KET side composite shell quartet 
			                                   ///< record coe offset
			vector<ShellQuartet> inputSQList;  ///< decompose the input sq into single ones 
			vector<ShellQuartet> derivSQList;  ///< the shell quartet list for deriv job
			DerivInfor  derivInfor;            ///< the derivative information corresponding to derivSQList

			//
			// here we also set up a vector to take care of the VRR 
			// result shell quartet list
			// some times the VRR result can be passed into the nonRR
			// or the derviatives section. Suggesting that the nonRR 
			// or derivatives section is in file split, but the HRR1
			// section (just after VRR) is in variable form; then
			// the result of VRR will be fully in variable form, too.
			// This array is to correct this problem
			// we note, that this only happen to the VRR
			//
			vector<ShellQuartet> vrrSQInArray; ///< the shell quartet list must be in array form

			//
			// here we also set up a vector to take care of the HRR 
			// result shell quartet list
			// some times the HRR result can be passed the derviatives 
			// section just like the above situation. 
			//
			vector<ShellQuartet> hrrSQInArray; ///< the shell quartet list must be in array form

			///
			/// for the input shell codes, we form the rest of other data here
			/// 
			void formSQInfor();

			///
			/// according to the input derivInfor, we will add derivative information on the 
			/// input shell quartet list so that to form the output sqlist
			///
			void formDerivSQList(const DerivInfor& derInfor, vector<ShellQuartet>& sqlist);

			///
			/// form the derivInfor, as well as derivSQList for the derivatives job
			/// because the derivInfor and derivSQList must be formed together,
			/// therefore we do them together
			///
			void formDerivInfor();

			///
			/// determine that whether we do HRR work?
			/// we note, that this is only an estimation. We may have 
			/// cases that here we say it do HRR work, however practically
			/// it does not. However, this will not affect the real situation.
			///
			bool weDOHRRWork() const;

		public:

			///
			/// constructor for infor class
			///
			SQIntsInfor(const int& oper0, const Infor& infor, const int& codeBra1, 
					const int& codeBra2, const int& codeKet1,
					const int& codeKet2);

			///
			/// default detructor
			///
			~SQIntsInfor() { };

			///
			/// append the given module name to the section information
			///
			void appendCodeSection(const int& section) { sectionInfor.push_back(section); };

			///
			/// let's see what's the next code section for the given module?
			///
			int nextSection(const int& sec) const; 

			///
			/// whether the given section is the last section
			///
			bool isLastSection(const int& sec) const {
				if (sectionInfor[0] == sec) return true;
				return false;
			};

			///
			/// whether the section information array has this section?
			///
			bool hasSection(const int& sec) const {
				bool has = false;
				for(int i=0; i<(int)sectionInfor.size(); i++) {
					if (sectionInfor[i] == sec) {
						has = true;
						break;
					}
				}
				return has;
			};

			///
			/// return the section information array
			///
			const vector<int>& getSectionInfor() const { return sectionInfor; };

			///
			/// whether the given sq is in the result list?
			/// for derivOrder > 0, we search derivSQList, else
			/// for derivORder = 0, we search inputSQList
			///
			bool isResult(const ShellQuartet& sq) const;

			///
			/// whether is the composite sq?
			///
			bool isComSQ() const {
				if (inputSQList.size() > 1) return true;
				return false;
			};

			///
			/// whether the input SQ are all bottom ones?
			///
			bool areAllBottomSQ() const;

			///
			/// wether the calculation will go with exponential
			/// factors modification to VRR?
			///
			bool withExpFac() const {
				// for derivatives job, we have exponential 
				// factors adding in 
				if (getJobOrder() > 0) return true;

				// finally let's see whether it involves 
				// non-RR work like three body KI
				// this case it also needs modifier
				if (oper == THREEBODYKI) return true;

				// in default we return false
				return false;
			};

			///
			/// return the total number of integrals considering both 
			/// of the input shell quartet as well as the derivatives
			/// situation. For derivOrder = 0, the nInts() returns the 
			/// number of integrals determined by shell quartet, else
			/// it returns the number of integrals times number of derivatives
			///
			int nInts() const;

			///
			/// for the result shell quartet and it's index, we get its offset
			///
			/// the offset is computed in this way:
			/// - if the derivatives order is 0, then we will compute that
			/// what's the index for this given integral formed by sq and index
			/// in the whole composite shell quartet
			///
			/// - if the derivatives order is > 0, we will compute the index
			/// for the given integral not only within the composite shell
			/// quartet, but also among all of derivatives
			///
			int getOffset(const ShellQuartet& sq, const int& index) const;

			///
			/// for a given shell quartet, get it's coefficient 
			/// array offset
			/// we note that the shell quartet may not need to 
			/// be the input shell quartet. What we do here is 
			/// to use division information to get the offset 
			/// for coefficients, so for the derivative job, too
			///
			void getCoeOffset(const ShellQuartet& sq, 
					int& ic2Offset, int& jc2Offset) const;

			///
			/// copy the given file split status here as a record
			///
			void setFileSplitMode(const int& moduleName, bool doFileSplit);

			///
			/// return whether we use array for the given code section
			/// this is totally determined by the file split status
			///
			bool withArrayIndex(const int& codeSec) const { 
				if (codeSec == DERIV) {
					return derivInSingleFile;
				}else if (codeSec == NON_RR) {
					return nonRRInSingleFile;
				}else if (codeSec == HRR2) {
					return hrr2InSingleFile;
				}else if (codeSec == HRR1) {
					return hrr1InSingleFile;
				}else if (codeSec == VRR) {
					return vrrInSingleFile;
				}else if (codeSec == VRR_CONT) {
					return vrrContractionSplit;
				}else{
					cout << "codeSec " << codeSec << endl;
					crash(true,"the input code section name is invalid in withArrayIndex of SQIntsInfor");
				}
				return false;
			};

			///
			/// whether the given module is in file split mode?
			///
			bool fileSplit(const int& codeSec) const { 
				if (codeSec == DERIV) {
					return derivInSingleFile;
				}else if (codeSec == NON_RR) {
					return nonRRInSingleFile;
				}else if (codeSec == HRR2) {
					return hrr2InSingleFile;
				}else if (codeSec == HRR1) {
					return hrr1InSingleFile;
				}else if (codeSec == VRR) {
					return vrrInSingleFile;
				}else if (codeSec == VRR_CONT) {
					return vrrContractionSplit;
				}else{
					crash(true,"the input code section name is invalid in fileSplit of SQIntsInfor");
				}
				return false;
			};

			///
			/// whether the code is with array index
			/// if we have any file split mode, then this must be 
			/// with array index
			///
			bool withArray() const {
				if (vrrInSingleFile || vrrContractionSplit || hrr1InSingleFile || hrr2InSingleFile || 
						nonRRInSingleFile || derivInSingleFile) return true;
				return false;
			};

			///
			/// get the working file name(with .cpp and path), or the pure function name
			/// according to the user's choice
			/// \param inFuncName  whether the function returns the function name
			/// \param moduleName  the code section name, the allowed names are VRR,VRR_HEAD,HRR1 etc.
			///                    please see the general.h see the definitions
			/// \param iFile       if the module need to be decomposed into several files, this is the 
			///                    file index, we note that iFile starts from 1. If the iFile is set 
			///                    to 0, this will contain the function prototype information
			/// \param finalFile   whether this is the final cpp file?
			///
			string getWorkFuncName(bool inFuncName, int moduleName, int iFile = -1, 
					bool finalFile = false) const;

			///
			/// get the function name for sqints
			///
			string getFuncName() const;

			///
			/// get the operator
			///
			int getOper() const { return oper; };

			///
			/// get the job order
			///
			int getJobOrder() const { return derivOrder; };

			///
			/// can we really do HRR?
			///
			bool hasHRR() const { return doHRRWork; }; 

			/**
			 * return the input sq list
			 */
			const vector<ShellQuartet>& getInputSQList() const {
				return inputSQList;
			};

			/**
			 * return the deriv sq list
			 */
			const vector<ShellQuartet>& getDerivSQList() const {
				return derivSQList;
			};

			/**
			 * determine the total length of coefficient array for bra/ket side
			 */
			int getCoeArrayLength(const int& side) const;

			/**
			 * whether the code will be with scr form of vector?
			 */
			bool withSCRVec() const {
				if (withArray()) {
					if (useSCRVec()) return true;
				}
				return false;
			};

			/**
			 * whether the code will be with TBB or STD form of vector?
			 */
			bool withDoubleVec() const {
				if (withArray()) {
					if (! useSCRVec()) return true;
				}
				return false;
			};

			/**
			 * whether the code will use the boost library for imcomplete gamma function?
			 */
			bool useBoostGamma() const;

			///
			/// get the module name in string format
			///
			string getModuleName(int name) const {
				if (name == VRR) {
					return "VRR";
				}else if (name == HRR1) {
					return "HRR1";
				}else if (name == HRR2) {
					return "HRR2";
				}else if (name == NON_RR) {
					return "NON_RR";
				}else if (name == DERIV) {
					return "DERIV";
				}else{
					crash(true,"the input module name is invalid in SQIntsInfor::getModuleName");
				}
				return "NONE";
			};

			///
			/// get the integral number of records in string format
			///
			string getRedundantIntEvalNumber(int pos) const {
				return boost::lexical_cast<string>(derivRecords[pos]);
			};

			///
			/// get the number of shell quartets per each derivative information
			///
			int getNSQPerDeriv() const;

			///
			/// return the redundant position in string format
			///
			string getRedundantPos() const;

			///
			/// return the deriv infor
			///
			const DerivInfor&  getDerivInfor() const { return derivInfor; };

			///
			/// return the input shell code
			///
			const vector<int>& getShellCodeArray() const {
				return inputShellCodes;  
			};

			///
			/// from the output SQ list, we will pick up those who can not 
			/// do HRR work for the both sides. Which is to say, these sq
			/// will be passed to the VRR module for solving. However, because
			/// they will be used by the section of nonRR, deriv etc. and 
			/// these sections uses array form of input, therefore we need 
			/// them to be in array form, too
			///
			void updateVRRSQListInArray(const vector<ShellQuartet>& outputSQList); 

			///
			/// return the vector of vrrSQInArray
			///
			const vector<ShellQuartet>& getVRRSQInArray() const { return vrrSQInArray; };

			///
			/// from the output SQ list, we will pick up those who can do HRR work
			/// for either side
			///
			void updateHRRSQListInArray(const vector<ShellQuartet>& outputSQList); 

			///
			/// return the vector of vrrSQInArray
			///
			const vector<ShellQuartet>& getHRRSQInArray() const { return hrrSQInArray; };

			///
			/// printing out the head part for the cpp files
			/// it could be either work file(like vrr, hrr etc.)
			/// or the top cpp function file
			///
			void headPrinting(ofstream& file) const;

			///
			/// according to the job information, return the argument list
			/// for the top cpp function
			///
			string getArgList() const;

			///
			/// update the information from vrrinfor
			///
			void updateVRRInfor(const VRRInfor& vrrinfor);
	};

}


#endif

