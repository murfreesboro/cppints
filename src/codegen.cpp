//
// CPPINTS: A C++ Program to Generate Analytical Integrals Based on Gaussian
// Form Primitive Functions
//
// Copyright (C) 2012-2014 Fenglai Liu
// This softare uses the MIT license as below:
//
//	Permission is hereby granted, free of charge, to any person obtaining 
//	a copy of this software and associated documentation files (the "Software"), 
//	to deal in the Software without restriction, including without limitation 
//	the rights to use, copy, modify, merge, publish, distribute, sublicense, 
//	and/or sell copies of the Software, and to permit persons to whom the Software 
//	is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all 
// copies or substantial portions of the Software.
//						    
//	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
//	INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
//	PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
//	FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
//	ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
//
#include "infor.h"
#include "inttype.h"
#include "shellsymbol.h"
#include "sqints.h"
using namespace infor;
using namespace inttype;
using namespace sqints;

void codeGen(const Infor& infor, const int& oper, const int& derivOrder) 
{
	// get the max L from infor
	int max_L = infor.maxL;
	int aux_max_L = infor.auxMaxL;

	// before we regulaly recreated the codes
	// we may create codes not according to 
	// the order of operator
	
	//
	// additional note for two/three body ERI
	// and normal 4 body ERI:
	// in 4 body ERI we use shell pair order array
	// to generate LCode
	// and in two/three body ERI we still use 
	// the shell code. However, we note that 
	// all of the three are consistent with each
	// other in terms of LCode generation
	//

	
	// two body ERI
	// we note, that both of the two shell dimensions
	// should corresponding to aux shell, therefore
	// should be higher than the ordinary max_l
	if (oper == TWOBODYERI) {
		for(int n2=0; n2<MAX_INPUT_SHELL_TYPES; n2++) {
			for(int n1=n2; n1<MAX_INPUT_SHELL_TYPES; n1++) {

				// get the input shell code
				int L1    = INPUT_SHELL_ANG_MOM_CODE[n1]; // bra1
				int L2    = INPUT_SHELL_ANG_MOM_CODE[n2]; // ket1
				int LMin1 = -1;
				int LMax1 = -1;
				decodeL(L1,LMin1,LMax1);
				if (LMax1>aux_max_L) continue;
				int LMin2 = -1;
				int LMax2 = -1;
				decodeL(L2,LMin2,LMax2);
				if (LMax2>aux_max_L) continue;

				// now generate codes
				SQInts sqints(infor,L1,S,L2,S,ERI,derivOrder);
				sqints.codeGeneration();
			}
		}
		return;
	}

	// three body ERI
	// we note, that ket1 side shell dimension
	// should corresponding to aux shell, therefore
	// should be higher than the ordinary max_l
	if (oper == THREEBODYERI) {
		for(int n3=0; n3<MAX_INPUT_SHELL_TYPES; n3++) {
			for(int n2=0; n2<MAX_INPUT_SHELL_TYPES; n2++) {
				for(int n1=n2; n1<MAX_INPUT_SHELL_TYPES; n1++) {

					// get the input shell code
					int L1    = INPUT_SHELL_ANG_MOM_CODE[n1]; // bra1
					int L2    = INPUT_SHELL_ANG_MOM_CODE[n2]; // bra2
					int L3    = INPUT_SHELL_ANG_MOM_CODE[n3]; // ket1
					int LMin1 = -1;
					int LMax1 = -1;
					decodeL(L1,LMin1,LMax1);
					if (LMax1>max_L) continue;
					int LMin2 = -1;
					int LMax2 = -1;
					decodeL(L2,LMin2,LMax2);
					if (LMax2>max_L) continue;
					int LMin3 = -1;
					int LMax3 = -1;
					decodeL(L3,LMin3,LMax3);
					if (LMax3>aux_max_L) continue;

					// now generate codes
					SQInts sqints(infor,L1,L2,L3,S,ERI,derivOrder);
					sqints.codeGeneration();
				}
			}
		}
		return;
	}

	// check the number of body for the operator
	int nBody = getOperOrder(oper);

	// one body integral
	if (nBody == 1) {
		for(int n1=0; n1<MAX_INPUT_SHELL_TYPES; n1++) {

			// get the input shell code
			int L1   = INPUT_SHELL_ANG_MOM_CODE[n1];
			int lmin = -1;
			int lmax = -1;
			decodeL(L1,lmin,lmax);
			if (lmax>max_L) continue;

			// now generate codes
			SQInts sqints(infor,L1,oper,derivOrder);
			sqints.codeGeneration();
		}
	}

	// two body integral
	// n1 is for bra1, n2 for bra2
	// we have the requirement that L(bra1) >= L(bra2)
	// therefore n1>=n2
	if (nBody == 2) {
		if (oper == MOM) {

			//
			// for moment integrals, L indicates the moment order
			// Currently we have P, D and F provided (see the mom head
			// function in sqints.cpp, we have bottom integrals only
			// up to f so far). Therefore if you want to climb to 
			// higher order, please modify the code over there first
			//
			for(int L=P; L<=aux_max_L; L++) {
				for(int n2=0; n2<MAX_INPUT_SHELL_TYPES; n2++) {
					for(int n1=n2; n1<MAX_INPUT_SHELL_TYPES; n1++) {

						// get the input shell code
						int L1    = INPUT_SHELL_ANG_MOM_CODE[n1]; // bra1
						int L2    = INPUT_SHELL_ANG_MOM_CODE[n2]; // bra2
						int LMin1 = -1;
						int LMax1 = -1;
						decodeL(L1,LMin1,LMax1);
						if (LMax1>max_L) continue;
						int LMin2 = -1;
						int LMax2 = -1;
						decodeL(L2,LMin2,LMax2);
						if (LMax2>max_L) continue;

						// now generate codes
						SQInts sqints(infor,L1,L2,L,oper,derivOrder);
						sqints.codeGeneration();
					}
				}
			}
		}else{
			for(int n2=0; n2<MAX_INPUT_SHELL_TYPES; n2++) {
				for(int n1=n2; n1<MAX_INPUT_SHELL_TYPES; n1++) {

					// get the input shell code
					int L1    = INPUT_SHELL_ANG_MOM_CODE[n1]; // bra1
					int L2    = INPUT_SHELL_ANG_MOM_CODE[n2]; // bra2
					int LMin1 = -1;
					int LMax1 = -1;
					decodeL(L1,LMin1,LMax1);
					if (LMax1>max_L) continue;
					int LMin2 = -1;
					int LMax2 = -1;
					decodeL(L2,LMin2,LMax2);
					if (LMax2>max_L) continue;

					// now generate codes
					SQInts sqints(infor,L1,L2,oper,derivOrder);
					sqints.codeGeneration();
				}
			}
		}
	}

	// three body integral
	// n1 is for bra1, n2 for bra2, n3 for ket1
	// we have the requirement that L(bra1) >= L(bra2)
	// therefore n1>=n2
	if (nBody == 3) {
		for(int n3=0; n3<MAX_INPUT_SHELL_TYPES; n3++) {
			for(int n2=0; n2<MAX_INPUT_SHELL_TYPES; n2++) {
				for(int n1=n2; n1<MAX_INPUT_SHELL_TYPES; n1++) {

					// get the input shell code
					int L1    = INPUT_SHELL_ANG_MOM_CODE[n1]; // bra1
					int L2    = INPUT_SHELL_ANG_MOM_CODE[n2]; // bra2
					int L3    = INPUT_SHELL_ANG_MOM_CODE[n3]; // ket1
					int LMin1 = -1;
					int LMax1 = -1;
					decodeL(L1,LMin1,LMax1);
					if (LMax1>max_L) continue;
					int LMin2 = -1;
					int LMax2 = -1;
					decodeL(L2,LMin2,LMax2);
					if (LMax2>max_L) continue;
					int LMin3 = -1;
					int LMax3 = -1;
					decodeL(L3,LMin3,LMax3);
					if (LMax3>max_L) continue;

					// now generate codes
					SQInts sqints(infor,L1,L2,L3,oper,derivOrder);
					sqints.codeGeneration();
				}
			}
		}
	}

	// four body integral
	// we use the shell pair order array to form the LCode
	if (nBody == 4) {
		for(int n2=0; n2<MAX_SHELL_PAIR_NUMBER; n2++) {
			for(int n1=n2; n1<MAX_SHELL_PAIR_NUMBER; n1++) {

				// get the input shell code
				int LBra  = SHELL_PAIR_ORDER_ARRAY[n1];
				int LKet  = SHELL_PAIR_ORDER_ARRAY[n2];
				int LMin1 = -1;
				int LMax1 = -1;
				int LMin2 = -1;
				int LMax2 = -1;
				decodeSQ(LBra,LMin1,LMax1,LMin2,LMax2);
				if (LMax1>max_L) continue;
				if (LMax2>max_L) continue;
				int LMin3 = -1;
				int LMax3 = -1;
				int LMin4 = -1;
				int LMax4 = -1;
				decodeSQ(LKet,LMin3,LMax3,LMin4,LMax4);
				if (LMax3>max_L) continue;
				if (LMax4>max_L) continue;

				// next we have to form the LCode for single shell
				int L1 = codeL(LMin1,LMax1);
				int L2 = codeL(LMin2,LMax2);
				int L3 = codeL(LMin3,LMax3);
				int L4 = codeL(LMin4,LMax4);


				// now generate codes
				SQInts sqints(infor,L1,L2,L3,L4,oper,derivOrder);

				// test that whether this file already been exist?
				// this is because we may do two/thre body eri
				// ahead
				if (sqints.isFileExist()) continue;
				sqints.codeGeneration();
			}
		}
	}
}
