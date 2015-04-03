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
 * \brief   describing the shells 
 * \author  Fenglai Liu
 * \note
 *
 * This module is used describe the class related to the "Shell". 
 * Shell actually is a group of basis set functions in the quantum chemistry,
 * all of these basis set functions share the same L, namely:
 * 
 * L = l+m+n
 * 
 * For example, Shell of L=1 has theree basis set functions, namely
 * Px  1,0,0(l,m,n)
 * Py  0,1,0
 * Pz  0,0,1
 *
 * In the recursive relation, similar to the basis set; we only care 
 * about the angular information of the shell. Furthermore, in RR we 
 * only use the non-composite shells. 
 */
#ifndef SHELL_H
#define SHELL_H
#include "general.h"

namespace basis {
	class Basis;
}

namespace shell {

	using namespace basis;

	class Shell {

		private:

			int L;     ///< the angular momentum for the shell

		public:

			// constructors and destructors
			Shell(int L0):L(L0) {
				crash(L< NULL_POS,"Imporper angular momentum passed in constructing shell");
				crash(L> 20,"The angular momentum you give is tooooo high! What's wrong?");
			};
			Shell():L(NULL_POS){ };
			Shell(const Shell& s0):L(s0.L) { };
			~Shell() { };

			// raising/lowering the angular momentum
			void changeL(int incL) {
				L = L + incL;
				if (L<0) L = NULL_POS;
			};

			// eq 
			bool operator==(const Shell& s) const {
				if(L == s.L) {
					return true;
				}else{
					return false;
				}
			};

			// not equal
			bool operator!=(const Shell& s) const {
				if(L != s.L) {
					return true;
				}else{
					return false;
				}
			};

			/**
			 * through the index of the basis set,
			 * we can get the l m and n so to get the corresponding basis set
			 */
			void getBasisSetFromIndex(const int& i, int& l, int& m, int& n) const;

			/**
			 * operator for copy assign
			 */
			Shell& operator=(const Shell& s) {
				if (this == &s) {
					return *this;
				}else{
					L = s.L;             
					return *this;
				}
			};

			/**
			 * whether the shell is null shell?
			 */
			bool isnull() const {
				if (L == NULL_POS) return true;
				return false;
			};

			int getL() const {
				return L;
			};

			/**
			 * get the number of Cartesian type of basis sets for
			 * this shell
			 */
			int getBasisSetNumber() const {
				return (L+1)*(L+2)/2;
			};

			/**
			 * get the basis sets for this shell
			 */
			void getBasis(vector<Basis>& basis) const;

			/**
			 * testing that whether we have this basis set
			 */
			bool hasThisBasisSet(const Basis& b) const;

			/**
			 * get the name of the shell
			 */
			string getName() const;
	};

}

#endif

