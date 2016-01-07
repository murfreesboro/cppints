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
 * \file    derivinfor.h
 * \brief   parsing the derivatives information 
 * \author  Fenglai Liu
 */
#ifndef DERIVINFOR_H
#define DERIVINFOR_H
#include "general.h"

namespace derivinfor {

	///
	/// this is a foundamental class to define the 1st order
	/// derivative information
	///
	class FirstOrderDerivInfor {

		private:

			int derivPos;            ///< bra1, bra2, ket1, ket2?
			int derivDirection[3];   ///< derivatives on x, y and z

		public:

			///
			/// constructor
			///
			FirstOrderDerivInfor(int pos):derivPos(pos) {
				derivDirection[0] = DERIV_X;
				derivDirection[1] = DERIV_Y;
				derivDirection[2] = DERIV_Z;
			};

			///
			/// default constructor
			///
			FirstOrderDerivInfor():derivPos(NULL_POS) {
				derivDirection[0] = NO_DERIV;
				derivDirection[1] = NO_DERIV;
				derivDirection[2] = NO_DERIV;
			};

			///
			/// destructor
			///
			~FirstOrderDerivInfor() { };

			///
			/// get the deriv pos
			///
			int getDerivPos() const { return derivPos; };

			///
			/// get the deriv direction for given index
			///
			int getDerivDirection(int i) const { return derivDirection[i]; };

			///
			/// return the length of first derivatives
			/// it's a constant, 3
			///
			int getDerivDirLen() const { return 3; }; 

			///
			/// print out the dirivative position and direction 
			/// information
			///
			void print(ofstream& file) const; 
	};

	///
	/// this is a foundamental class to define the 2ed order
	/// derivative information
	/// 
	/// for the direction array, it's a list of x, y or z
	/// we express the direction like this. For example,
	/// if the 1st deriv and 2ed deriv are on same position,
	/// for example, both are BRA1; then the direction array 
	/// will be like this:
	/// 1st direction:   x  x  x  y  y  z  
	/// 2ed direction:   x  y  z  y  z  z
	/// this representing the xx, xy, xz, yy, yz and zz
	/// in this case the rest of elements are NULL(-1)
	///
	/// if they are not at on the same position, for example; pos1 is BRA1
	/// and pos2 is BRA2; then the direction array will be like:
	/// 1st direction:   x  x  x  y  y  y  z  z  z  
	/// 2ed direction:   x  y  z  x  y  z  x  y  z
	/// this representing the xx, xy, xz, yx, yy, yz and zx, zy and zz
	///
	/// for the given shell quartet, one of the thing we assume here;
	/// is that (ab|cd) the shells a,b,c and d are different from
	/// each other. There's possibility that a may same as b; however
	/// in general we treat them differently
	///
	class SecondOrderDerivInfor {

		private:

			int firstDerivPos;            ///< bra1, bra2, ket1 or ket2 on the first derivatives position?
			int secondDerivPos;           ///< bra1, bra2, ket1 or ket2 on the second derivatives position?
			int firstDerivDirection[9];   ///< x, y or z information for first derivatives position?
			int secondDerivDirection[9];  ///< x, y or z information for second derivatives position?

		public:

			///
			/// default constructor
			///
			SecondOrderDerivInfor():firstDerivPos(NULL_POS),secondDerivPos(NULL_POS) {
				firstDerivDirection[0]  = NO_DERIV;
				firstDerivDirection[1]  = NO_DERIV;
				firstDerivDirection[2]  = NO_DERIV;
				firstDerivDirection[3]  = NO_DERIV;
				firstDerivDirection[4]  = NO_DERIV;
				firstDerivDirection[5]  = NO_DERIV;
				firstDerivDirection[6]  = NO_DERIV;
				firstDerivDirection[7]  = NO_DERIV;
				firstDerivDirection[8]  = NO_DERIV;
				secondDerivDirection[0] = NO_DERIV;
				secondDerivDirection[1] = NO_DERIV;
				secondDerivDirection[2] = NO_DERIV;
				secondDerivDirection[3] = NO_DERIV;
				secondDerivDirection[4] = NO_DERIV;
				secondDerivDirection[5] = NO_DERIV;
				secondDerivDirection[6] = NO_DERIV;
				secondDerivDirection[7] = NO_DERIV;
				secondDerivDirection[8] = NO_DERIV;
			};

			///
			/// constructor
			///
			SecondOrderDerivInfor(int pos1, int pos2):firstDerivPos(pos1),secondDerivPos(pos2) {

				// only two possibility case
				// one is two position are identical
				// the other one is that they are different
				if (firstDerivPos == secondDerivPos) {

					// xx
					firstDerivDirection[0]  = DERIV_X;  
					secondDerivDirection[0] = DERIV_X;

					// xy
					firstDerivDirection[1]  = DERIV_X;  
					secondDerivDirection[1] = DERIV_Y;

					// xz
					firstDerivDirection[2]  = DERIV_X;  
					secondDerivDirection[2] = DERIV_Z;

					// yy
					firstDerivDirection[3]  = DERIV_Y;  
					secondDerivDirection[3] = DERIV_Y;

					// yz
					firstDerivDirection[4]  = DERIV_Y;  
					secondDerivDirection[4] = DERIV_Z;
					
					// zz
					firstDerivDirection[5]  = DERIV_Z;  
					secondDerivDirection[5] = DERIV_Z;

					// the rest are null
					firstDerivDirection[6]  = NULL_POS;  
					secondDerivDirection[6] = NULL_POS;
					firstDerivDirection[7]  = NULL_POS;  
					secondDerivDirection[7] = NULL_POS;
					firstDerivDirection[8]  = NULL_POS;  
					secondDerivDirection[8] = NULL_POS;
				}else{

					// xx
					firstDerivDirection[0]  = DERIV_X;  
					secondDerivDirection[0] = DERIV_X;

					// xy
					firstDerivDirection[1]  = DERIV_X;  
					secondDerivDirection[1] = DERIV_Y;

					// xz
					firstDerivDirection[2]  = DERIV_X;  
					secondDerivDirection[2] = DERIV_Z;

					// yx
					firstDerivDirection[3]  = DERIV_Y;  
					secondDerivDirection[3] = DERIV_X;

					// yy
					firstDerivDirection[4]  = DERIV_Y;  
					secondDerivDirection[4] = DERIV_Y;
					
					// yz
					firstDerivDirection[5]  = DERIV_Y;  
					secondDerivDirection[5] = DERIV_Z;

					// zx
					firstDerivDirection[6]  = DERIV_Z;
					secondDerivDirection[6] = DERIV_X;
					
					// zy
					firstDerivDirection[7]  = DERIV_Z;
					secondDerivDirection[7] = DERIV_Y;
					
					// zz
					firstDerivDirection[8]  = DERIV_Z;
					secondDerivDirection[8] = DERIV_Z;
				}
			};

			///
			/// destructor
			///
			~SecondOrderDerivInfor() { };

			///
			/// get the deriv pos
			///
			void getDerivPos(int& pos1, int& pos2) const { 
				pos1 = firstDerivPos;  
				pos2 = secondDerivPos;     
			};

			///
			/// get the deriv direction for given index
			///
			void getDerivDirection(int i, int& dir1, int& dir2) const { 
					dir1 = firstDerivDirection[i];
					dir2 = secondDerivDirection[i];
			};

			///
			/// return the length of first derivatives
			/// these information also constant
			///
			int getDerivDirLen() const { 
				if (firstDerivPos == secondDerivPos) return 6;
				return 9;
			}; 

			///
			/// print out the dirivative position and direction 
			/// information
			///
			void print(ofstream& file) const; 
	};

	///
	/// based on the constant classes above (FirstOrderDerivInfor, SecondOrderDerivInfor etc.)
	/// we are able to form the derivative information for the given type of integrals.
	///
	/// Because of translation invariance, it's not all of derivatives that needs to be
	/// calculated. The derivatives needs to be computed called "unique derivatives",
	/// these derivatives could be derived from the unique ones because of translation 
	/// invariance are called "redundant derivatives". For more information about this,
	/// the user should refer to the manual.
	///
	/// the redundant derivatives have a feature, that this is all related to one given
	/// center. Therefore, such center is called redundant position. By given the 
	/// redundant position, we are able to figure out what are the unique derivatives
	/// and what are the redudant ones. 
	///
	/// we note that here in the infor array only unique ones are stored.
	///
	///
	class DerivInfor {

		private:

			int oper;             ///< the copy of integral operator
			int jobOrder;         ///< derivative order, 1, 2, 3 or even larger
			int redundantPos;     ///< redundant position
			vector<FirstOrderDerivInfor>  firstDerivInfor;    ///< first order deriv information
			vector<SecondOrderDerivInfor> secondDerivInfor;   ///< second order deriv information

		public:

			///
			/// constructor for the integral derivatives without the redundant position
			/// information. That means, the  redundant  position is set to NULL
			/// number of body integral will be specified by the operator
			///
			DerivInfor(const int& oper0, const int& jobOrder0);

			///
			/// constructor for the integral derivatives with redundant position information
			/// number of body integral will be specified by the operator
			///
			DerivInfor(const int& oper0, const int& jobOrder0, const int& redundantPos0);

			///
			/// default constructor - everything is meaningless
			///
			DerivInfor():oper(NULL_POS),jobOrder(NULL_POS),redundantPos(NULL_POS) { };

			///
			/// check the input information
			/// to see whether we can proceed
			///
			bool checkInfor();

			///
			/// destructor
			///
			~DerivInfor() { };

			///
			/// get the length of first deriv infor array
			///
			int getLen1stDerivInforArray() const; 

			///
			/// get the length of second deriv infor array
			///
			int getLen2edDerivInforArray() const; 

			///
			/// return the element in the first order infor array
			///
			const FirstOrderDerivInfor& get1stDerivInfor(int i) const {
			 return firstDerivInfor[i];
			};

			///
			/// return the element in the second order infor array
			///
			const SecondOrderDerivInfor& get2edDerivInfor(int i) const {
			 return secondDerivInfor[i];
			};

			///
			/// get the redundant pos
			///
			int getRedundantPos() const { return redundantPos; };

			///
			/// for the first order derivatives information (given position 
			/// and associated direction), we calculate it's offset which
			/// refers to it's position in the whole information array
			///
			int getOffset(const int& derivPos, const int& derivDir) const; 

			///
			/// for the second order derivatives information (given positions 
			/// and associated directions), we calculate it's offset which
			/// refers to it's position in the whole information array
			///
			int getOffset(const int& p1, const int& p2, 
					const int& dir1, const int& dir2) const; 

			///
			/// get the total number of derivatives
			///
			int getTotalNumDeriv() const;
	};
}

#endif
