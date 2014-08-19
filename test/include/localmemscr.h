/**
 * \file    localmemscr.h
 * \author  Fenglai Liu and Jing Kong
 */
#ifndef LOCALMEMSCR_H
#define LOCALMEMSCR_H
#include "libgen.h"
#include<vector>
typedef std::vector<Double>   DoubleVec;
typedef std::vector<UInt>     UIntVec;
typedef std::vector<Int>      IntVec;

namespace localmemscr {


	/**
	 * \class   LocalMemScr
	 * \brief   manipulating local memory for single thread
	 *
	 * Currently the local memory scratch is only for double type of data
	 *
	 * The implementation of this class is very simple. However, it's not 
	 * flexible; for example; if you have a lot of places which need 
	 * memory, and these places are totally independent with each other;
	 * such class will cause memory waste since it will always hold the 
	 * memory until reset is called.
	 *
	 * Use it with your caution!
	 */
	class LocalMemScr {

		private:

			UInt lastMemUser;      ///< the last memory user position
			UInt lastMemLen;       ///< for the last memory user, it's required mem length
			DoubleVec memHolder;   ///< holding memory for the scratch

		public:

			/**
			 * initilize the local memory sratch
			 */
			LocalMemScr(const UInt& len):lastMemUser(0),lastMemLen(0),
			memHolder(len,ZERO) { };

			/**
			 * destructor
			 */
			~LocalMemScr() { };

			/**
			 * return a new memory position according to the length
			 * it required
			 */
			Double* getNewMemPos(const UInt& len) {

				// firstly check that whether there's no memory usr
				// if in this case, things is very simple
				if (lastMemLen == 0) {

					// check that whether the length is less than the 
					// whole memory length
					if (len>memHolder.size()) {
						crash(true,"localMemScr required memory is larger than what we have in total");
					}

					// now it's safe to do the work
					lastMemLen = len;
					return &memHolder[0];
				}

				// secondly, let's go to see that wether we can
				// append a new memory position to the memory user
				UInt lastPos    = lastMemUser;
				UInt lastLen    = lastMemLen;
				UInt newPos     = lastPos+lastLen;
				if (newPos>=memHolder.size() || newPos+len>=memHolder.size()) {
					crash(true, "in getNewMemPos vector overflow");
				}
				lastMemUser = newPos;
				lastMemLen  = len;
				return &memHolder[newPos];
			};

			/**
			 * reset the local memory scratch
			 * to erase all of memory users information
			 */
			void reset() {
				lastMemUser = 0;
				lastMemLen  = 0;
				memHolder.assign(memHolder.size(),ZERO);
			};
	};

}

#endif
