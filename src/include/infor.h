/**
 *
 * CPPINTS: A C++ Program to Generate Analytical Integrals Based on Gaussian
 * Form Primitive Functions
 *
 * Copyright (C) 2012-2015 Fenglai Liu
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
 * \file    infor.h
 * \brief   processing the input information for the whole program
 * \author  Fenglai Liu
 */
#ifndef INFOR_H
#define INFOR_H
#include "general.h"
#define NO_ENFORCE_ON_RR       0   // no enforce for var/array form on RR, choose by program itself
#define ENFORCE_RR_WITH_ARRAY  1   // enforce to use array on RR
#define ENFORCE_RR_WITH_VAR    2   // enforce to use variable on RR
#define TBB_VEC                3   // we use std vector with TBB allocator
#define STD_VEC                4   // we use std vector with STD allocator
#define USE_SCR_VEC            5   // the vector is provided by scratch memory management class

namespace infor {

	/**
	 * \class LineParse
	 * \brief LineParse is used to parsing a line of information from the given text file.
	 */
	class LineParse {

		private:

			string         line;                ///< the line read in from text file
			vector<string> pieces;              ///< pieces getting from parsing the line
			int            number;              ///< number of pieces getting from the line
			bool           isEmpty;             ///< whether this is empty line
			bool           isComment;           ///< whether this is comment line

		public:
			LineParse(string input); 
			~LineParse() { };

			string findValue(int i)  const {
				crash(i>=number,"Pieces array is overflowed");
				return pieces.at(i);
			};

			int getNPieces () const {
				return number;
			};

			bool isCom() const {
				return isComment;
			};

			bool isEmp() const {
				return isEmpty;
			};
	};

	/**
	 * \class  WordConvert
	 * \brief  WordConvert is used to convert the word into different forms. For example,
	 * to convert a string into double or int. This is a utility class
	 */
	class WordConvert {

		public:
			WordConvert(){ };
			~WordConvert(){ };
			bool toInt(string s, int& x);
			bool compare(string s1, string s2);
			void capitalize(string& s);
			bool isInt(string s);
	};

	class Infor {

		public:

			//
			// additional note for the variable/array form usage in VRR/HRR:
			//
			// to use var/array form of variable, it's determined by several factors:
			// 1 if the var/array form is enforced by the user (this is determined in
			// enforceVRR and enforceHRR), then we will use the user specified form;
			// 2 if the user does not provide an enforced form, then var/aaray form
			// usage will be determined by several factors:
			//   2.1  maxL in terms of printing (maxL_vrrPrinting and maxL_hrrPrinting);
			//   2.2  nBody determined by the operator;
			//   2.3  whether we do HRR
			// please see the sqintsinfor.cpp (constructor) for more information.
			//
			// additional note for fmt function:
			// for f_{m}(t) calculation, you are able to choose M limit value from
			// 8, 9 or 10. The corresponding fmt error is 1.0e-13, 1.0e-13 and 
			// 1.0e-12. For more details, please see the doc file; where we provide
			// a hybrid scheme for computing the fmt function. 
			// fmt_error is an integer. For example, error is 1.0e-13 then it's 13.
			// default we use M_limit = 10
			//
			//
			int enforceVRR;        ///< enforced VRR part to use variable/array
			int enforceHRR;        ///< enforced HRR part to use variable/array
			int maxL;              ///< maximum angular momentum in integral generation
			int auxMaxL;           ///< maximum angular momentum for integral with aux shell
			int derivOrder;        ///< derivatives order for the integral code
			int maxL_vrrPrinting;  ///< maximum angular momentum used in setting VRR printing
			int maxL_hrrPrinting;  ///< maximum angular momentum used in setting HRR printing
			int maxL_singleFile;   ///< maximum angular momentum to keep single cpp file
			int vec_form;          ///< what kind of vector form we use? See above definition for TBB_VEC etc.
			int M_limit;           ///< in calculating fmt function, what is the M limit?
			int fmt_error;         ///< the fmt function error associating with fmt function
			string hrr_method;     ///< HRR method
			string vrr_method;     ///< VRR method
			vector<int> joblist;   ///< jobs going to be performed, namely integral operator

			///
			/// constructor for infor class
			///
			Infor(const string& input);

			///
			/// default destructor
			///
			~Infor() { };

			///
			/// whehter the user defined the HRR method?
			///
			bool hasHRR() const { 
				if (hrr_method == "none") return false;
				return true;
			};

			///
			/// whether the program use scr form to declare vectors?
			///
			bool useSCRVec() const { return (vec_form == USE_SCR_VEC); };

			///
			/// whether the program use the vector with TBB allocator?
			///
			bool useTBBVec() const { return (vec_form == TBB_VEC); };

			///
			/// whether the program use the standard vector with STD allocator?
			///
			bool useSTDVec() const { return (vec_form == STD_VEC); };

			///
			/// according to the vector usage type, return the vector type
			///
			string getArrayType() const {
				if (useSCRVec()) {
					return "Double* ";
				}else{
					return "DoubleVec ";
				}
			};

			///
			/// according to the vector usage type, return the vector definition part
			/// getNewMemPos is the member function for oject of "scr" in LocalMemScr class
			/// it will return a pointer for the given length of memory
			///
			string getArrayDeclare(const string& nInts) const {
				if (useSCRVec()) {
					return " = scr.getNewMemPos(" + nInts +");";
				}else{
					return "(" + nInts + ",0.0E0);";
				}
			};

			///
			/// get the VRR method
			///
			int getVRRMethod() const;

			///
			/// member function to return the job list
			///
			const vector<int>& getJobList() const { return joblist; };

			/**
			 * get the general project name (based on HRR and VRR)
			 */
			string getProjectName() const; 

			/**
			 * based on the operator, as well as integral file name
			 * and order, it returns where it's located (Dir)
			 * \param oper      operator
			 * \param fileName  file name for the project
			 * \param order     which derivatives order it's in?
			 * \param withTmpWorkDir whether it's involved with work dir?
			 */
			string getProjectFileDir(const int& oper, const string& fileName, 
					int order, bool withTmpWorkDir = true) const; 

			/**
			 * based on the operator, and order, it returns the tmp dir location
			 * \param oper      operator
			 * \param order     which derivatives order it's in?
			 */
			string getProjectTmpFileDir(const int& oper, int order) const; 

	};

}


#endif

