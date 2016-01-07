#ifndef SHELLPROP_H
#define SHELLPROP_H
#include "libgen.h"

namespace shellprop {

	/************************************************************************
	 *                          Data for Shell Structure                    *
	 ************************************************************************/
	/**
	 * shell types 
	 * 
	 * By given this shell types, we are able to calculate the number of 
	 * basis set functions for each shell(n is the given angular momentum
	 * number):
	 * 
	 * Cartesian:
	 * \f$\frac{(n+1)(n+2)(n+3)-n(n+1)(n+2)}{6}\f$
	 *
	 * Spherical:
	 * \f$(n+1)^{2}-n^{2}\f$
	 *
	 */
	const Int S = 0;
	const Int P = 1;
	const Int D = 2;
	const Int F = 3;
	const Int G = 4;
	const Int H = 5;
	const Int I = 6;
	const Int K = 7;
	const Int L = 8;
	const Int M = 9;

	/**
	 * maximum number of angular momentum allowed in the program
	 * for K,L,M we defined but we do not want user to go to
	 * such high L
	 */
	const Int MAX_L = 20;

	/**
	 * shell names supported in the program
	 * SPD would be split into SP and D in
	 * real calculation
	 */
	const string SHELL_NAME_LIST[] = {
		"S","SP","P","D","F","G","H","I", "K", "L", "M", "N",
		"O", "Q", "R", "T", "U", "V", "W", "X", "Y", "Z"
	};

	/**
	 * the code for angular momentum numbers
	 * associated with the shell names
	 *
	 * How to code the angular momentum for each shell?
	 * the unit digit is the Lmin,
	 * the tenth digit is the Lmax
	 */
	const Int SHELL_ANG_MOM_CODE[] = {
		0, 100, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20
	};

	/**
	 * this value indicates the maximum number of shell names
	 */
	const Int MAX_SHELL_TYPES = 22;

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
	const Int SHELL_PAIR_ORDER_ARRAY[ ]  =  { 
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
	const Int MAX_SHELL_PAIR_NUMBER = 253;

	/**
	 * here it's the data to define the unit to form
	 * angular momentum code and LCode
	 * such way is derived from the program to generate
	 * shell pair data above
	 */
	const Int ANG_UNIT = 100;
	const Int LCODE_UNIT_BRA2 = 1000;
	const Int LCODE_UNIT_KET1 = 1000000;
	const Int LCODE_UNIT_KET2 = 1000000000;

	/**
	 * code angular momentum for a shell 
	 * this should match the WORKING_SHELL_CODE array above
	 */
	inline Int codeL(Int lmin, Int lmax) {
		if (lmin == lmax) {
			return lmax;
		}else{
			return (lmin+lmax*ANG_UNIT);
		}
	};

	/**
	 * decode L from the angular momentum code above
	 */
	inline void decodeL(const Int& code, Int& lmin, Int& lmax) {
		if (code < ANG_UNIT) {
			lmin = code;
			lmax = code;
		}else{
			lmax = code/ANG_UNIT;
			lmin = code - lmax*ANG_UNIT;
		}
	};

	/**
	 * decode the angular momentum for the shell pairs (i,j|
	 */
	inline void decodeSQ(const Int& code, Int& iLmin, Int& iLmax, Int& jLmin, Int& jLmax) {
		Int jCode = code/LCODE_UNIT_BRA2;
		Int iCode = code - jCode*LCODE_UNIT_BRA2;
		decodeL(iCode,iLmin,iLmax);
		decodeL(jCode,jLmin,jLmax);
	};

	/**
	 * code the angular momentum for the one body integrals
	 */
	inline LInt codeSQ(Int iL) 
	{
		return (LInt)iL;
	}

	/**
	 * code the angular momentum for the two body integrals
	 *	i corresponds to row shell, j corresponds to col shell, integral will be (i|j)
	 */
	inline LInt codeSQ(Int iL, Int jL) 
	{
		LInt  iL0 = iL;
		LInt  jL0 = jL;
		LInt  c0  = jL0*LCODE_UNIT_BRA2;  
		LInt code = iL0 + c0;
		return code;
	}

	/**
	 * code the angular momentum for three body integral
	 *	the corresponding integral will be (ij|k)
	 */
	inline LInt codeSQ(Int iL, Int jL, Int kL) 
	{
		LInt  iL0 = iL;
		LInt  jL0 = jL;
		LInt  kL0 = kL;
		LInt  c0  = jL0*LCODE_UNIT_BRA2;  
		LInt  c1  = kL0*LCODE_UNIT_KET1;  
		LInt code = iL0 + c0 + c1;
		return code;
	}

	/**
	 * code the angular momentum for four body integral
	 *	the corresponding integral will be (ij|kl)
	 */
	inline LInt codeSQ(Int iL, Int jL, Int kL, Int lL) 
	{
		LInt  iL0 = iL;
		LInt  jL0 = jL;
		LInt  kL0 = kL;
		LInt  lL0 = lL;
		LInt  c0  = jL0*LCODE_UNIT_BRA2;  
		LInt  c1  = kL0*LCODE_UNIT_KET1;  
		LInt  c2  = lL0*LCODE_UNIT_KET2;  
		LInt code = iL0 + c0 + c1 + c2;
		return code;
	}

	/**
	 * number of Cartesian type of basis set functions
	 */
	inline Int getCartBas(const Int& lmin, const Int& lmax) {
		return ((lmax+1)*(lmax+2)*(lmax+3)-lmin*(lmin+1)*(lmin+2))/6;
	};

	/**
	 * return the offset for composite shell
	 * for example, we try to get the D2X basis set's global index in
	 * SPD shell. In this case, lmin is S, L is D and index is D2x index
	 * in D shell
	 */
	inline Int getBasOffset(const Int& lmin, const Int& L, const Int& index) {
		
		// if L is same with lmin, the index is just the global index
		if (L == lmin) {
			return index;
		}

		// now consider more
		Int lmax = L-1;
		return ((lmax+1)*(lmax+2)*(lmax+3)-lmin*(lmin+1)*(lmin+2))/6 + index;
	};
}

#endif

