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
 * \file    shell.h
 * \brief   describing the basis set in quantum chemistry
 * \author  Fenglai Liu
 * \note
 *
 * This module is used describe the class related to the "basis set functions". 
 * Originally, each basis set function is a combination of Gaussian primitive
 * functions:
 *
 * psi = sum_{mu}d_{mu}chi_{mu}
 * 
 * psi is the basis set function, and chi_{mu} is the primitive functions.
 * All of chi are on the same center as psi, and d_{mu} is some fixed 
 * coefficients. All of Gaussian primitive functions share the same angular
 * momentum with the basis set.
 * 
 * For each Gaussian primitive function, it has the form that:
 * 
 * chi = x^{l}y^{m}z^{n}e^{-alpha*r^{2}}
 * 
 * x^{l}y^{m}z^{n} is its angular momentum part, which is characterized by
 * three number of l, m, and n. The e^{-alpha*r^{2}} is its radial part,
 * so l,m,n combined with alpha and its prefactor of d_{mu}, then we know
 * all of information to get psi.
 *
 * In the recursive relation, we only care about the angular part information.
 * The radial information is directly folded in the S integral calculation
 * code. Therefore, we only use l,m,n to represent a shell in RR.
 *
 */
#ifndef BASIS_H
#define BASIS_H
#include "general.h"

namespace basis {

	class Basis {

		private:

			int l;    ///< angular momentum for X
			int m;    ///< angular momentum for Y
			int n;    ///< angular momentum for Z

		public:

			Basis(int l0, int m0, int n0): l(l0), m(m0), n(n0) {
				crash(l<NULL_POS||m<NULL_POS||n<NULL_POS,"Improper lmn given in the Basis");
			};
			Basis():l(NULL_POS), m(NULL_POS), n(NULL_POS) { };
			Basis(const Basis& b2):l(b2.l),m(b2.m),n(b2.n) { };
			~Basis() { };

			// eq 
			bool operator==(const Basis& b) const {
				if(l == b.l && m == b.m && n == b.n) {
					return true;
				}else{
					return false;
				}
			};

			// not equal
			bool operator!=(const Basis& b) const {
				if(l != b.l || m != b.m || n != b.n) {
					return true;
				}else{
					return false;
				}
			};

			/**
			 * operator for copy assign
			 */
			Basis& operator=(const Basis& b) {
				if (this == &b) {
					return *this;
				}else{
					l = b.l;             
					m = b.m;             
					n = b.n;             
					return *this;
				}
			};

			/**
			 * get the name of the basis set
			 */
			string getName() const;

			void getlmn(int& l0, int& m0, int& n0) const {
				l0 = l;
				m0 = m;
				n0 = n;
			};

			int getL() const {
				return l + m + n;
			};

			bool isnull() const {
				if (l == NULL_POS||m == NULL_POS||n == NULL_POS) return true;
				return false;
			};

			/**
			 * return the local index of this basis set in the shell
			 */
			int getLocalIndex() const;
	};
}


#endif

