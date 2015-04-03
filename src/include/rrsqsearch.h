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
 * \file    rrsqsearch.h
 * \brief   initilize the shell quartet list on RR process
 * \author  Fenglai Liu
 */
#ifndef RRSQSEARCH_H
#define RRSQSEARCH_H
#include "general.h"
#include "shellquartet.h"
using namespace shellquartet;

namespace rrbuild {
	class RRBuild;
}

namespace rrsqsearch {

	using namespace rrbuild;

	//
	// here we will define some macro used to pick up
	// the shell quartets into the unsolvedSQList in
	// RRSQSearch process. See pickupUnsolvedSQ function
	//
	const int PICKUP_SQ_ON_L_M    = 0; // pick up shell quartets based on given L and M
	const int PICKUP_SQ_ON_L_OPER = 1; // pick up shell quartets based on given L and operator 

	///
	/// define the state of shell quartets used in RRShellQuart
	/// new shell quartets mean that these shell quartets
	/// would be counted in current RR searching;
	/// old shell quartets mean that these shell quartets 
	/// are already been considered in previous RR searching;
	/// null shell quartets means these shell quartet is not
	/// real ones.
	/// old shell quartets and null shell quartets are
	/// viewed as same in RR search
	///
	const int NEW_SQ   = 0;
	const int OLD_SQ   = 1;
	const int NULL_SQ  = 2;

	/**
	 * \class RRShellQuart
	 *
	 *	This class is used in the RR search process. In this class,
	 *	we generate all of possible RR expansions for the given shell
	 *	quartet. The steps are as follows:
	 *	-  counting all of non-S position for RR expansion;
	 *	-  for each non-S position, we do possible RR expansion and record
	 *	   the information (update sqlist, poslist etc.)
	 *	-  before the real use of RRShellQuart, we need to update the 
	 *	   state of the sqlist. Since sqlist may contain the shell 
	 *	   quartets already known, the integral number for these shell quartets
	 *	   are not counted in RR search process(this is called vertical comparison)
	 */
	class RRShellQuart {

		private:

			int nRRItems;                  ///< number of items used in RR
			int nPos;                      ///< number os possible non-s positions in oriSQ
			ShellQuartet oriSQ;            ///< original shell quartet, appearing on LHS
			vector<ShellQuartet> sqlist;   ///< records all shell quartets on each 
			                               ///< possible postion
			vector<int> intNumList;        ///< record the total integral number for each 
			                               ///< shell quartets
			vector<int> posList;           ///< RR position list in expansion
			vector<int> sqStateList;       ///< record the state of sq in sqlist.

		public:

			/**
			 * build the RRShellQuart in RR search process for VRR
			 * In RR search, the most important thing is to generate all of 
			 * possible positions for the given shell quartet. 
			 */
			RRShellQuart(const int& rrType, const ShellQuartet& sq);

			/**
			 * build the RRShellQuart in RR search process for VRR
			 * In RR search, the most important thing is to generate all of 
			 * possible positions for the given shell quartet. 
			 *
			 * abandoned
			 */
			//RRShellQuart(const ShellQuartet& sq, const int& side);

			/**
			 * destructor
			 */
			~RRShellQuart() { };

			/**
			 * In the RR search step, we need to set up the state for the new
			 * created sqlist. we want to make sure that every sq in the sqset
			 * will not be counted in the sqlist(that is to say, its state is "old")
			 * \param sqset the shell quartets used for comparison, in list form
			 */
			void updateStateSQList(const list<ShellQuartet>& sqset);

			/**
			 * In the RR search step, we need to set up the state for the new
			 * created sqlist. we want to make sure that every sq in the sqset
			 * will not be counted in the sqlist(that is to say, its state is "old")
			 * \param sqset the shell quartets used for comparison, in vector form
			 */
			void updateStateSQList(const vector<ShellQuartet>& sqset);

			/**
			 * for the RR search, return the integral number for the given sq
			 * \param pos  the variable position starting from 0 - it's not var itself!
			 * \param iSQ  which RHS shell quartet we are referring?
			 */
			int intNumCount(const int& pos, const int& iSQ) const {
				return intNumList[pos*nRRItems+iSQ];
			};

			/**
			 * for the RR search, return the RHS shell quartet state
			 * \param pos  the variable position starting from 0 - it's not var itself!
			 * \param iSQ  which RHS shell quartet we are referring?
			 */
			bool isNewSQ(const int& pos, const int& iSQ) const {
				int p = pos*nRRItems+iSQ;
				if (sqStateList[p] == NEW_SQ) return true;
				return false;
			};

			/**
			 * for the RR search, return the RHS shell quartet 
			 * \param pos  the variable position starting from 0 - it's not var itself!
			 * \param iSQ  which RHS shell quartet we are referring?
			 */
			const ShellQuartet& getRHSSQ(const int& pos, const int& iSQ) const {
				int p = pos*nRRItems+iSQ;
				return sqlist[p];
			};

			/*
			 * get the number of possible positions
			 */
			int getNPos() const {
				return nPos;
			};

			/*
			 * get the number of RR items
			 */
			int getNRRItems() const {
				return nRRItems;
			};

			/**
			 * we try to estimate the minimum integral number only
			 * from the rrSQ information
			 * \return minPos        estimate position for minimum integrals
			 * \return minIntNumber  estimate minimum integral number
			 */
			void getMinIntPos(int& minPos, int& minIntNumber) const;

			/**
			 * get the variable for the corresponding position
			 */
			int getVar(const int& pos) const {
				return posList[pos];
			};

			/**
			 * according to the given position, we retrieve the shell quartets
			 * on the RR RHS so that to start a new round of search 
			 * all of the null or old shell quartets are omitted
			 */
			void updateRHSSQList(list<ShellQuartet>& rhsList, const int& pos) const;

			/**
			 * get the LHS sq
			 */
			const ShellQuartet& getLHSSQ() const { return oriSQ; };

			/**
			 * debug printing
			 */
			void print() const;
	};

	/**
	 * This class is used to build the shell quartet list used in RR process
	 * furthermore, it's position in RR would also be determined
	 */
	class RRSQSearch {

		private:

			int rrType;                        ///< RR type
			vector<int> posList;               ///< for each sq, where the RR expansion did
			vector<ShellQuartet> solvedSQList; ///< the solved shell quartet list used in RR
			list<ShellQuartet> unsolvedSQArch; ///< keep a record of unsolved sq as archive

			/*********************************************************
			 *         the following functions are all for VRR       *
			 *********************************************************/

			///
			///  for the given sqlist, we try to estimate the total number
			///  of RR combinations. If it's too large then it's not 
			///  applicable. Therefore we return false
			///
			bool canDoOptSearch(const list<ShellQuartet>& sqlist1, 
					const list<ShellQuartet>& sqlist2) const;

			///
			/// this function parses the input SQ to see that 
			/// whether we can determine the expanding position
			/// of the given shell quartet directly
			/// this occurs in several cases:
			/// 1 only one position is non-S shell, and all of others are all S shell;
			/// 2 we have two positions are P shell, the rest are S shell;
			///
			bool canDoDirectParse(const ShellQuartet& sq);

			///
			/// update the given unsolved shell quartet archive by comparing 
			/// with the solvedSQList. If the unsolvedSQArch has any elements
			/// already in solved sq list, then we erase this element
			///
			void updateUnsolvedSQArch();

			///
			/// function for pick up connected LHS shell quartets for the following 
			/// searchOptPos function based on L and M combination
			/// \param P1, P2  shell quartet properties for picking up (P1-L,P2-M)
			/// \return unsolvedMainSQList  the main unsolvedSQList whose pos will be determined
			/// \return unsolvedAppendSQList  the other unsolvedSQList helps to determine main ones
			///
			void pickupUnsolvedSQByLM(const int& P1, const int& P2, 
					list<ShellQuartet>& unsolvedMainSQList,
					list<ShellQuartet>& unsolvedAppendSQList);

			///
			/// general function for pick up connected LHS shell quartets for the following 
			/// searchOptPos function based on L and operator combination
			/// \param P1, P2  shell quartet properties for picking up (P1-L,P2-operator)
			/// \return unsolvedMainSQList  the main unsolvedSQList whose pos will be determined
			/// \return unsolvedAppendSQList  the other unsolvedSQList helps to determine main ones
			///
			void pickupUnsolvedSQByLOper(const int& P1, const int& P2, 
					list<ShellQuartet>& unsolvedMainSQList,
					list<ShellQuartet>& unsolvedAppendSQList);

			///
			/// general function for pick up connected LHS shell quartets for the following 
			/// searchOptPos function based on L itself
			/// \param P1  shell quartet properties for picking up, it's LSum
			/// \return unsolvedMainSQList  the main unsolvedSQList whose pos will be determined
			/// \return unsolvedAppendSQList  the other unsolvedSQList helps to determine main ones
			///
			void pickupUnsolvedSQByL(const int& P1, list<ShellQuartet>& unsolvedMainSQList,
					list<ShellQuartet>& unsolvedAppendSQList);

			///
			/// general function for pick up connected LHS shell quartets for the following 
			/// searchOptPos function
			/// \param pickupMethod     what kind of pick up way we use? See above macro
			/// \param P1, P2           shell quartet properties for picking up
			/// \return unsolvedMainSQList  the main unsolvedSQList whose pos will be determined
			/// \return unsolvedAppendSQList  the other unsolvedSQList helps to determine main ones
			///
			void pickupUnsolvedSQ(const int& pickupMethod, const int& P1,
					const int& P2, list<ShellQuartet>& unsolvedMainSQList,
					list<ShellQuartet>& unsolvedAppendSQList);

			///
			/// for the given a group of shell quartet list, we search its optimum positions 
			/// \param mainsqlist  unsolved shell quartets whose position are going to be determined
			/// \param appendsqlist unsolved shell quartets who are here for getting accurate
			/// optimum position for the main shell quartets (their RHS could overlap with
			/// the main ones)
			///
			void searchOptPos(const list<ShellQuartet>& mainsqlist, 
					const list<ShellQuartet>& appendsqlist); 

			///
			/// searching the shell quartets based on two shell quartet properties 
			/// the detailed algorithm see the comment inside
			/// \param method:      pick up unsolved shell quartets into list
			/// \param inputSQList: the input shell quartet list for searching
			///
			void RRSearchBasedOnTWOProperty(const int& method, 
					const vector<ShellQuartet>& initSQList);

			/*********************************************************
			 *         the following functions are all for HRR       *
			 *********************************************************/

			///
			/// whether the input shell quartet is the bottom of HRR process?
			///
			bool isBottomForHRR(const ShellQuartet& sq, const int& side);

			///
			/// top routine for HRR search 
			///
			void HRRSearch(const int& side, const vector<ShellQuartet>& initSQList);

			///
			/// perform the indirect HRR search for the given side with respect to the 
			/// given list of shell quartets
			/// the indirect search means we only know which side the work is going to 
			/// be performed, we need to further investigate the perferable position
			/// for expansion
			/// we return result position for the given side and the output 
			/// result shell quartet list
			///
			int indirectHRRSearch(const int& side, const vector<ShellQuartet>& initSQList, 
					list<ShellQuartet>& resultList);

			///
			/// perform the direct HRR search for the given side with respect to the 
			/// given list of shell quartets
			/// the direct search means we know which side the work is going to 
			/// be performed, we just expanding the RR on the given position 
			///
			void directHRRSearch(const int& pos, const vector<ShellQuartet>& initSQList, 
					list<ShellQuartet>& resultList);

		public:

			///
			/// for the given group of input shell quartet list, we will give a list
			/// of shell quartets that is appearing in their VRR expansion. 
			//  What's more, the optimum position for each sq would be solved
			///
			RRSQSearch(const vector<ShellQuartet>& inputSQList, const int& rrType0);

			///
			/// for the given group of input shell quartet list, we will give a list
			/// of shell quartet for the HRR expansion on the given side
			///
			RRSQSearch(const int& side, const vector<ShellQuartet>& inputSQList);

			///
			/// destructor
			///
			~RRSQSearch() { };

			///
			/// return the pos list
			///
			const vector<int>& getPosList() const { return posList; };

			///
			/// return the solved sq list
			///
			const vector<ShellQuartet>& getSolvedSQList() const { return solvedSQList; };

			///
			/// debug printing for VRR
			/// 
			void print(const int& method, const vector<ShellQuartet>& sqlist) const;

			///
			/// debug printing for HRR
			/// 
			void print() const;
	};


}


#endif

