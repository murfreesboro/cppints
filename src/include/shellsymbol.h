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
 * \file    shellsymbol.h
 * \brief   describing the shell symbols
 * \author  Fenglai Liu
 *
 * currently we support the L up to L = 20
 * this is defined in the codes. If you want 
 * to increase this limit, please revise 
 * the shell pair order array first, then 
 * revise other shell data
 */
#ifndef SHELLSYMBOL_H
#define SHELLSYMBOL_H
#include "general.h"
#include "boost/lexical_cast.hpp"
#include <boost/algorithm/string.hpp>
using boost::lexical_cast;

/**
 * the known shell symbols 
 * for the L higher than 20, there's currently no
 * corresponding shell symbols. We use L21, L22 etc.
 * instead. See the getShellSymbol for more information
 *
 * we note, that these shell symbols does not contain
 * the composite shells. They are used "INSIDE" the 
 * program. SQInts class handles the composite shells
 * from input, and convert all of composite shell
 * into single shells before any practical work
 *
 * We note, if higher, we will use something like "L21"
 * to express the shell
 */
const string SHELL_NAME_LIST[ ] = {
	"S", "P", "D", "F", "G", "H", "I", "K", "L", "M", "N",
	"O", "Q", "R", "T", "U", "V", "W", "X", "Y", "Z"
};

/**
 * total number of shells appearing in the above array
 */
const int SHELL_NUMBER = 21; 

/**
 * the current highest L support in this program
 * we note, that it's determined by shell pair order
 * array
 */
const int MAX_L = 20; 

/////////////////////////////////////////////////////////////////////////
// the data below is related to the shells used to generate integrals  //
/////////////////////////////////////////////////////////////////////////

/**
 * this is possible shell angular momentums used as input for integrals
 */
const int S = 0;
const int P = 1;
const int D = 2;
const int F = 3;
const int G = 4;
const int H = 5;
const int I = 6;
const int K = 7;
const int L = 8;
const int M = 9;

/**
 * input shell names supported in the program
 */
const string INPUT_SHELL_NAME_LIST[] = {
	"S","SP","P","D","F","G","H","I", "K", "L", "M", "N",
	"O", "Q", "R", "T", "U", "V", "W", "X", "Y", "Z"
};

/**
 * the code for angular momentum numbers
 * associated with the input shell names
 */
const int INPUT_SHELL_ANG_MOM_CODE[] = {
	0, 100, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20
};

/**
 * this value indicates the maximum number of shell names from input
 */
const int MAX_INPUT_SHELL_TYPES = 22;

/**
 * this shell pair order array defines the 
 * shell pair angular momentum order
 *
 * for two arbitrary shell pairs, it answers 
 * the question that which one should be ahead,
 * and which one should be after.
 *
 * The array is generated through the program
 * in shell/util/shell_pair_order program
 * Please check the main src folder
 *
 * this feature is used in creating LCode
 * when combining two shell pairs (see 4 body
 * integral generation in codegen.cpp)
 *
 * below the comment part is the Lmin and Lmax
 * corresponding to the LCode. For example,
 * 100001  , //  1   1   0   1
 * the LCode is 100001, the bra1 shell is P 
 * (lmin=1 lmax=1); and the bra2 shell is 
 * SP (0 1).
 */
const int SHELL_PAIR_ORDER_ARRAY[ ]  =  { 
        0  , //  0   0   0   0 
      100  , //  0   1   0   0 
        1  , //  1   1   0   0 
        2  , //  2   2   0   0 
   100100  , //  0   1   0   1 
   100001  , //  1   1   0   1 
     1001  , //  1   1   1   1 
        3  , //  3   3   0   0 
   100002  , //  2   2   0   1 
     1002  , //  2   2   1   1 
        4  , //  4   4   0   0 
   100003  , //  3   3   0   1 
     1003  , //  3   3   1   1 
     2002  , //  2   2   2   2 
        5  , //  5   5   0   0 
   100004  , //  4   4   0   1 
     1004  , //  4   4   1   1 
     2003  , //  3   3   2   2 
        6  , //  6   6   0   0 
   100005  , //  5   5   0   1 
     1005  , //  5   5   1   1 
     2004  , //  4   4   2   2 
     3003  , //  3   3   3   3 
        7  , //  7   7   0   0 
   100006  , //  6   6   0   1 
     1006  , //  6   6   1   1 
     2005  , //  5   5   2   2 
     3004  , //  4   4   3   3 
        8  , //  8   8   0   0 
   100007  , //  7   7   0   1 
     1007  , //  7   7   1   1 
     2006  , //  6   6   2   2 
     3005  , //  5   5   3   3 
     4004  , //  4   4   4   4 
        9  , //  9   9   0   0 
   100008  , //  8   8   0   1 
     1008  , //  8   8   1   1 
     2007  , //  7   7   2   2 
     3006  , //  6   6   3   3 
     4005  , //  5   5   4   4 
       10  , //  10  10  0   0 
   100009  , //  9   9   0   1 
     1009  , //  9   9   1   1 
     2008  , //  8   8   2   2 
     3007  , //  7   7   3   3 
     4006  , //  6   6   4   4 
     5005  , //  5   5   5   5 
       11  , //  11  11  0   0 
   100010  , //  10  10  0   1 
     1010  , //  10  10  1   1 
     2009  , //  9   9   2   2 
     3008  , //  8   8   3   3 
     4007  , //  7   7   4   4 
     5006  , //  6   6   5   5 
       12  , //  12  12  0   0 
   100011  , //  11  11  0   1 
     1011  , //  11  11  1   1 
     2010  , //  10  10  2   2 
     3009  , //  9   9   3   3 
     4008  , //  8   8   4   4 
     5007  , //  7   7   5   5 
     6006  , //  6   6   6   6 
       13  , //  13  13  0   0 
   100012  , //  12  12  0   1 
     1012  , //  12  12  1   1 
     2011  , //  11  11  2   2 
     3010  , //  10  10  3   3 
     4009  , //  9   9   4   4 
     5008  , //  8   8   5   5 
     6007  , //  7   7   6   6 
       14  , //  14  14  0   0 
   100013  , //  13  13  0   1 
     1013  , //  13  13  1   1 
     2012  , //  12  12  2   2 
     3011  , //  11  11  3   3 
     4010  , //  10  10  4   4 
     5009  , //  9   9   5   5 
     6008  , //  8   8   6   6 
     7007  , //  7   7   7   7 
       15  , //  15  15  0   0 
   100014  , //  14  14  0   1 
     1014  , //  14  14  1   1 
     2013  , //  13  13  2   2 
     3012  , //  12  12  3   3 
     4011  , //  11  11  4   4 
     5010  , //  10  10  5   5 
     6009  , //  9   9   6   6 
     7008  , //  8   8   7   7 
       16  , //  16  16  0   0 
   100015  , //  15  15  0   1 
     1015  , //  15  15  1   1 
     2014  , //  14  14  2   2 
     3013  , //  13  13  3   3 
     4012  , //  12  12  4   4 
     5011  , //  11  11  5   5 
     6010  , //  10  10  6   6 
     7009  , //  9   9   7   7 
     8008  , //  8   8   8   8 
       17  , //  17  17  0   0 
   100016  , //  16  16  0   1 
     1016  , //  16  16  1   1 
     2015  , //  15  15  2   2 
     3014  , //  14  14  3   3 
     4013  , //  13  13  4   4 
     5012  , //  12  12  5   5 
     6011  , //  11  11  6   6 
     7010  , //  10  10  7   7 
     8009  , //  9   9   8   8 
       18  , //  18  18  0   0 
   100017  , //  17  17  0   1 
     1017  , //  17  17  1   1 
     2016  , //  16  16  2   2 
     3015  , //  15  15  3   3 
     4014  , //  14  14  4   4 
     5013  , //  13  13  5   5 
     6012  , //  12  12  6   6 
     7011  , //  11  11  7   7 
     8010  , //  10  10  8   8 
     9009  , //  9   9   9   9 
       19  , //  19  19  0   0 
   100018  , //  18  18  0   1 
     1018  , //  18  18  1   1 
     2017  , //  17  17  2   2 
     3016  , //  16  16  3   3 
     4015  , //  15  15  4   4 
     5014  , //  14  14  5   5 
     6013  , //  13  13  6   6 
     7012  , //  12  12  7   7 
     8011  , //  11  11  8   8 
     9010  , //  10  10  9   9 
       20  , //  20  20  0   0 
   100019  , //  19  19  0   1 
     1019  , //  19  19  1   1 
     2018  , //  18  18  2   2 
     3017  , //  17  17  3   3 
     4016  , //  16  16  4   4 
     5015  , //  15  15  5   5 
     6014  , //  14  14  6   6 
     7013  , //  13  13  7   7 
     8012  , //  12  12  8   8 
     9011  , //  11  11  9   9 
    10010  , //  10  10  10  10
   100020  , //  20  20  0   1 
     1020  , //  20  20  1   1 
     2019  , //  19  19  2   2 
     3018  , //  18  18  3   3 
     4017  , //  17  17  4   4 
     5016  , //  16  16  5   5 
     6015  , //  15  15  6   6 
     7014  , //  14  14  7   7 
     8013  , //  13  13  8   8 
     9012  , //  12  12  9   9 
    10011  , //  11  11  10  10
     2020  , //  20  20  2   2 
     3019  , //  19  19  3   3 
     4018  , //  18  18  4   4 
     5017  , //  17  17  5   5 
     6016  , //  16  16  6   6 
     7015  , //  15  15  7   7 
     8014  , //  14  14  8   8 
     9013  , //  13  13  9   9 
    10012  , //  12  12  10  10
    11011  , //  11  11  11  11
     3020  , //  20  20  3   3 
     4019  , //  19  19  4   4 
     5018  , //  18  18  5   5 
     6017  , //  17  17  6   6 
     7016  , //  16  16  7   7 
     8015  , //  15  15  8   8 
     9014  , //  14  14  9   9 
    10013  , //  13  13  10  10
    11012  , //  12  12  11  11
     4020  , //  20  20  4   4 
     5019  , //  19  19  5   5 
     6018  , //  18  18  6   6 
     7017  , //  17  17  7   7 
     8016  , //  16  16  8   8 
     9015  , //  15  15  9   9 
    10014  , //  14  14  10  10
    11013  , //  13  13  11  11
    12012  , //  12  12  12  12
     5020  , //  20  20  5   5 
     6019  , //  19  19  6   6 
     7018  , //  18  18  7   7 
     8017  , //  17  17  8   8 
     9016  , //  16  16  9   9 
    10015  , //  15  15  10  10
    11014  , //  14  14  11  11
    12013  , //  13  13  12  12
     6020  , //  20  20  6   6 
     7019  , //  19  19  7   7 
     8018  , //  18  18  8   8 
     9017  , //  17  17  9   9 
    10016  , //  16  16  10  10
    11015  , //  15  15  11  11
    12014  , //  14  14  12  12
    13013  , //  13  13  13  13
     7020  , //  20  20  7   7 
     8019  , //  19  19  8   8 
     9018  , //  18  18  9   9 
    10017  , //  17  17  10  10
    11016  , //  16  16  11  11
    12015  , //  15  15  12  12
    13014  , //  14  14  13  13
     8020  , //  20  20  8   8 
     9019  , //  19  19  9   9 
    10018  , //  18  18  10  10
    11017  , //  17  17  11  11
    12016  , //  16  16  12  12
    13015  , //  15  15  13  13
    14014  , //  14  14  14  14
     9020  , //  20  20  9   9 
    10019  , //  19  19  10  10
    11018  , //  18  18  11  11
    12017  , //  17  17  12  12
    13016  , //  16  16  13  13
    14015  , //  15  15  14  14
    10020  , //  20  20  10  10
    11019  , //  19  19  11  11
    12018  , //  18  18  12  12
    13017  , //  17  17  13  13
    14016  , //  16  16  14  14
    15015  , //  15  15  15  15
    11020  , //  20  20  11  11
    12019  , //  19  19  12  12
    13018  , //  18  18  13  13
    14017  , //  17  17  14  14
    15016  , //  16  16  15  15
    12020  , //  20  20  12  12
    13019  , //  19  19  13  13
    14018  , //  18  18  14  14
    15017  , //  17  17  15  15
    16016  , //  16  16  16  16
    13020  , //  20  20  13  13
    14019  , //  19  19  14  14
    15018  , //  18  18  15  15
    16017  , //  17  17  16  16
    14020  , //  20  20  14  14
    15019  , //  19  19  15  15
    16018  , //  18  18  16  16
    17017  , //  17  17  17  17
    15020  , //  20  20  15  15
    16019  , //  19  19  16  16
    17018  , //  18  18  17  17
    16020  , //  20  20  16  16
    17019  , //  19  19  17  17
    18018  , //  18  18  18  18
    17020  , //  20  20  17  17
    18019  , //  19  19  18  18
    18020  , //  20  20  18  18
    19019  , //  19  19  19  19
    19020  , //  20  20  19  19
    20020    //  20  20  20  20
};

/**
 * this is the total number of shell pairs
 * possibly appears in the program
 */
const int MAX_SHELL_PAIR_NUMBER = 253;

/**
 * here it's the data to define the unit to form
 * angular momentum code and LCode
 * such way is derived from the program to generate
 * shell pair data above
 */
const int ANG_UNIT = 100;
const int LCODE_UNIT_BRA2 = 1000;
const int LCODE_UNIT_KET1 = 1000000;
const int LCODE_UNIT_KET2 = 1000000000;

/////////////////////////////////////////////////////////////////////////
//           functions related to the above constant data              //
/////////////////////////////////////////////////////////////////////////

/**
 * get the shell symbol for giving L
 * we note that this function can not handle the 
 * composite shell
 * therefore do not use this function in sqints class!
 */
inline string getShellSymbol(int L) {
	crash(L<0, "Improper L given in getShellSymbol");
	if (L<SHELL_NUMBER) return SHELL_NAME_LIST[L];
	string n = lexical_cast<string>(L);
	string LSym = "L"+n;
	return LSym;
};

//
// here below is the functions used to deal with
// possible composite shell in SQInts class
//

/**
 * is it the composite shell?
 */
inline bool isCompositeShell(int code) {
	if (code < ANG_UNIT) {
		return false;
	}else{
		return true;
	}
};

inline string getShellName(const int& code) 
{
	int i;
	for(i=0; i<MAX_INPUT_SHELL_TYPES; i++) 
		if (code == INPUT_SHELL_ANG_MOM_CODE[i]) break;
	if(i == MAX_INPUT_SHELL_TYPES) crash(true, "Failed to find shell name in getShellName!!");
	return INPUT_SHELL_NAME_LIST[i];
}

inline int getShellCode(const string& name) 
{
	string name1 = name;
	boost::to_upper(name1);
	int i;
	for(i=0; i<MAX_INPUT_SHELL_TYPES; i++) 
		if (name1 == INPUT_SHELL_NAME_LIST[i]) break;
	if(i == MAX_INPUT_SHELL_TYPES) crash(true, "Failed to find shell code in getShellCode!!");
	return INPUT_SHELL_ANG_MOM_CODE[i];
}

/**
 * code angular momentum for an input shell 
 */
inline int codeL(int lmin, int lmax) {
	if (lmin == lmax) {
		return lmax;
	}else{
		return (lmin+lmax*ANG_UNIT);
	}
};

/**
 * decode L from the angular momentum code above
 */
inline void decodeL(const int& code, int& lmin, int& lmax) {
	if (code < ANG_UNIT) {
		lmin = code;
		lmax = code;
	}else{
		lmax = code/ANG_UNIT;
		lmin = code - lmax*ANG_UNIT;
	}
};

/**
 * get the number of sub-shells from the composite shells
 */
inline int getNShells(const int& code) {
	if (code < ANG_UNIT) {
		return 1;
	}else{
		int lmax = code/ANG_UNIT;
		int lmin = code - lmax*ANG_UNIT;
		return (lmax-lmin+1);
	}
};

/**
 * decode the angular momentum for the shell pairs (i,j|
 */
inline void decodeSQ(const int& code, int& iLmin, int& iLmax, int& jLmin, int& jLmax) {
	int jCode = code/LCODE_UNIT_BRA2;
	int iCode = code - jCode*LCODE_UNIT_BRA2;
	decodeL(iCode,iLmin,iLmax);
	decodeL(jCode,jLmin,jLmax);
};

/**
 * code the angular momentum for the one body integrals
 */
inline long long codeSQ(int iL) 
{
	return (long long)iL;
}

/**
 * code the angular momentum for the two body integrals
 *	i corresponds to row shell, j corresponds to col shell, integral will be (i|j)
 */
inline long long codeSQ(int iL, int jL) 
{
	long long iL0  = iL;
	long long jL0  = jL;
	long long c0   = jL0*LCODE_UNIT_BRA2;
	long long code = iL0 + c0;
	return code;
}

/**
 * code the angular momentum for three body integral
 *	the corresponding integral will be (ij|k)
 */
inline long long codeSQ(int iL, int jL, int kL) 
{
	long long iL0  = iL;
	long long jL0  = jL;
	long long kL0  = kL;
	long long c0   = jL0*LCODE_UNIT_BRA2;
	long long c1   = kL0*LCODE_UNIT_KET1;
	long long code = iL0 + c0 + c1;
	return code;
}

/**
 * code the angular momentum for four body integral
 *	the corresponding integral will be (ij|kl)
 */
inline long long codeSQ(int iL, int jL, int kL, int lL) 
{
	long long iL0  = iL;
	long long jL0  = jL;
	long long kL0  = kL;
	long long lL0  = lL;
	long long c0   = jL0*LCODE_UNIT_BRA2;
	long long c1   = kL0*LCODE_UNIT_KET1;
	long long c2   = lL0*LCODE_UNIT_KET2;
	long long code = iL0 + c0 + c1 + c2;
	return code;
}

/**
 * calculate number of Cartesian type of basis set functions
 */
inline int getCartBas(const int& lmin, const int& lmax) {
	return ((lmax+1)*(lmax+2)*(lmax+3)-lmin*(lmin+1)*(lmin+2))/6;
}

/**
 * calculate the offset for the given shell in a composite
 * shell quartet, for example; P shell in SPD shell
 */
inline int getShellOffsetInCompositeShell(const int& lmin, const int& L) {
	if (L == lmin) return 0;
	int l = L - 1;
	return ((l+1)*(l+2)*(l+3)-lmin*(lmin+1)*(lmin+2))/6;
};

#endif

