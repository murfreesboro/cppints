#ifndef DERIVINFOR_H
#define DERIVINFOR_H
#include "libgen.h"

namespace derivinfor {

	///
	/// function to convert the position into string
	///
	inline string stringPos(Int pos) {
		if (pos == BRA1) return "BRA1";
		if (pos == BRA2) return "BRA2";
		if (pos == KET1) return "KET1";
		if (pos == KET2) return "KET2";
		crash(true, "invalid position in stringPos function");
		return "NONE";
	};

	///
	/// function to convert the direction into string
	///
	inline string stringDir(Int dir) {
		if (dir == DERIV_X) return "X";
		if (dir == DERIV_Y) return "Y";
		if (dir == DERIV_Z) return "Z";
		crash(true, "invalid direction in stringDir function");
		return "NONE";
	};

	class DerivInfor {

		private:

			LInt code;                    ///< the code file name
			Int order;                    ///< the order of the derivatives
			Int firstDerivPos;            ///< bra1, bra2, ket1 or ket2 on the first derivatives position?
			Int secondDerivPos;           ///< bra1, bra2, ket1 or ket2 on the second derivatives position?
			Int derivLength;              ///< 1st order deriv is always 3, 2ed deriv could be 6/9
			Int firstDerivDirection[9];   ///< x, y or z information for first derivatives position?
			Int secondDerivDirection[9];  ///< x, y or z information for second derivatives position?

		public:

			///
			/// default constructor
			///
			DerivInfor(LInt code0, Int order0):code(code0),order(order0),
			firstDerivPos(-1),secondDerivPos(-1),derivLength(0) {
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
			/// add in the pos information
			///
			void updatePos(Int pos1, Int pos2) {
				firstDerivPos = pos1;
				secondDerivPos = pos2;
			};

			///
			/// add in the dir information
			///
			void updateDir(Int dir1, Int dir2) {
				if (derivLength>=9) {
					crash(true, "illegal deriv length we have in updateDir");
				}
				firstDerivDirection[derivLength] = dir1;
				secondDerivDirection[derivLength] = dir2;
				derivLength++;
			};

			///
			/// destructor
			///
			~DerivInfor() { };

			///
			/// get the deriv pos
			///
			void getDerivPos(Int& pos1, Int& pos2) const { 
				pos1 = firstDerivPos;  
				pos2 = secondDerivPos;     
			};

			///
			/// get the deriv pos
			///
			Int getDerivPos() const { 
				return firstDerivPos;  
			};

			///
			/// get the deriv direction for given index
			///
			void getDerivDirection(Int i, Int& dir1, Int& dir2) const { 
				dir1 = firstDerivDirection[i];
				dir2 = secondDerivDirection[i];
			};

			///
			/// get the deriv direction for given index
			///
			Int getDerivDirection(Int i) const { 
				return firstDerivDirection[i];
			};

			///
			/// return the length of first derivatives
			/// these information also constant
			///
			Int getDerivDirLen() const { 
				return derivLength;
			}; 
	};

	class DerivInforArray {

		private:

			Int jobOrder;         ///< derivative order, 1, 2, 3 or even larger
			LInt code;            ///< integral code
			vector<DerivInfor>  derivInforVec;    

		public:

			DerivInforArray(Int jobOrder0, LInt code0):jobOrder(jobOrder0),code(code0) { };

			///
			/// destructor
			///
			~DerivInforArray() { };

			///
			/// read in the file so that to form the derivatives information 
			/// vector
			/// 
			void readInformation(const string& fileName);

			///
			/// get the length of first deriv infor array
			///
			Int getLenDerivInforArray() const { return derivInforVec.size(); }; 

			///
			/// return the element in the first order infor array
			///
			const DerivInfor& getDerivInfor(Int i) const {
				return derivInforVec[i];
			};

			///
			/// get the total number of derivatives
			///
			Int getTotalNumDeriv() const;

			///
			/// for debugging purpose
			///
			void print() const; 
	};
}

#endif
