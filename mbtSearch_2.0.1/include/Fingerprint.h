// Fingerprint.h
//
// Copyright (c) 2016
// Universitaet Duisburg-Essen
// Campus Duisburg
// Institut fuer Soziologie
// Prof. Dr. Rainer Schnell
// Lotharstr. 65
// 47057 Duisburg 
//
// This file is part of the command line application "mbtSearch".
//
// "mbtSearch" is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// "mbtSearch" is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with "mbtSearch". If not, see <http://www.gnu.org/licenses/>.

#ifndef FINGERPRINT_H
#define FINGERPRINT_H

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Misc.h"

typedef unsigned int WORDTYPE;			// 32bit-words
#define WORD_LEN 32				// word-length in bits
#define BIT1 1u					// unsigned int literal 1

#define FOLDED_WORDS (128/WORD_LEN)		// array-length for hash-key

// Objects of class Fingerprint hold a single bit-vector of arbitrary length
// The class also provides functions for modifying and checking single bits of
// the bit-vector and to compute its cardinality.
// Abstract super class of specialized classes FingerprintTanimoto and
// FingerprintHamming
class Fingerprint {
	protected:

	char *mId;				// pointer to ID string
	WORDTYPE *mArray;			// array for stored bits
	WORDTYPE mHashArray[FOLDED_WORDS];	// 128 Bit folded Hash-Key
	int mLength;				// length of fingerprint in bits
	static int sCardinalityMap[0x10000];	// static 16-bit cardinality-map
	
	// calculate length of word-array
	inline int arrayLength() {
		return (mLength - 1) / WORD_LEN + 1;
	}
	
	// initialize word-array
	inline void allocate() {
		mLength = MAX(mLength, 128);
		mArray = new WORDTYPE[arrayLength()];
	}

	// calculate cardinality of 32bit-word
	inline int cardWord(WORDTYPE word) {
		return sCardinalityMap[word & 0xFFFF] + sCardinalityMap[(word >> 16) & 0xFFFF];
	}
	
	// constructor for empty fingerprint of given bit-length
	inline Fingerprint(int length) {
		mId = NULL;
		mLength = length;
		allocate();
		clear();
		fold();
	}
	
	// constructor for fingerprint based on given ascii-string
	inline Fingerprint(char * id, const char *str) {
		int l = 0;
		
		mId = id;
		
		while ((str[l] == '0') || (str[l] == '1')) {
			l++;
		}

		mLength = l;
		
		allocate();
		clear();
	
		for (int i = 0; i < l; i++) {
			if (str[i] != '0') {
				setBit(i);
			}
		}

		fold();
	}
	
	// copy-constructor for fingerprints
	inline Fingerprint(Fingerprint *print) {
	
		// copy id
		if (print->mId != NULL) {
			mId = new char[strlen(print->mId) + 1];
			strcpy(mId, print->mId);
		} else {
			mId = NULL;
		}
		
		// copy word-array 
		mLength = print->mLength;
		allocate();
	
		for (int i = 0; i < arrayLength(); i++) {
			mArray[i] = print->mArray[i];
		}

		// copy hash-key
		for (int i = 0; i < FOLDED_WORDS; i++) {
			mHashArray[i] = print->mHashArray[i];
		}
	}
	
	// destructor
	inline ~Fingerprint() {
		if (mId != NULL) {
			delete[] mId;
		}
		delete[] mArray;
	}

	public:

	// get id
	inline char *getId() {
		return mId;
	}

	// clear all bits
	inline void clear() {
		for (int i = 0; i < arrayLength(); i++) {
			mArray[i] = 0;
		}
	}
	
	// check if all bits are cleared
	inline int isEmpty() {
		for (int i = 0; i < arrayLength(); i++) {
			if (mArray[i] != 0) {
				return 0;
			}
		}
	
		return 1;
	}

	// checks if two prints are equal
	inline int isEqual(Fingerprint *print) {
		for (int i = 0; i < arrayLength(); i++) {
			if (mArray[i] != print->mArray[i]) {
				return 0;
			}
		}

		return 1;
	}
	
	// get length in bits
	inline int getLength() {
		return mLength;
	}
	
	// count set bits
	inline int cardinality() {
		int count = 0;
	
		for (int i = 0; i < arrayLength(); i++) {
			count += cardWord(mArray[i]);
		}
	
		return count;
	}

	// compute cardinality of intersection of two prints
	inline int intersectCard(Fingerprint *print) {
		int count = 0;

		for (int i = 0; i < arrayLength(); i++) {
			count += cardWord(mArray[i] & print->mArray[i]);
		}

		return count;
	}

	// get bit at position n
	inline WORDTYPE getBit(int n) {
		return (mArray[n / WORD_LEN] >> (n % WORD_LEN)) & BIT1;
	}
	
	// set bit at position n
	inline void setBit(int n) {
		mArray[n / WORD_LEN] |= (BIT1 << (n % WORD_LEN));
	}
	
	// unset bit at position n
	inline void unsetBit(int n) {
		mArray[n / WORD_LEN] &= 0xFFFF - (BIT1 << (n % WORD_LEN));
	}

	// join one bits of second print
	inline void join(Fingerprint *print) {
		for (int i = 0; i < arrayLength(); i++) {
			mArray[i] = mArray[i] | print->mArray[i];
		}
	}

	// compute hash-key
	inline void fold() {
		int len;
		
		len = arrayLength();
		
		for (int i = 0; i < FOLDED_WORDS; i++) {
			mHashArray[i] = mArray[i];
		}
		
		for (int i = FOLDED_WORDS; i < len; i++) {
			mHashArray[i % FOLDED_WORDS] ^= mArray[i];
		}
	}
	
	// static class function to initialise cardinality-map
	static void init();
	
	// convert fingerprint to ascii-string of size n
	void copyToString(char *str, int n);
};

// specialized Fingerprint class for computing Tanimoto similarity
class FingerprintTanimoto : public Fingerprint {
	public:

	inline FingerprintTanimoto(int length) : Fingerprint(length) {}
	inline FingerprintTanimoto(char * id, const char *str) : Fingerprint(id, str) {}
	inline FingerprintTanimoto(FingerprintTanimoto *print) : Fingerprint((Fingerprint*) print) {}
	inline int isEqual(Fingerprint *print) { return Fingerprint::isEqual((Fingerprint*) print); }
	inline void join(Fingerprint *print) { Fingerprint::join((Fingerprint*) print); }
	inline int intersectCard(Fingerprint *print) { return Fingerprint::intersectCard((Fingerprint*) print); }

	typedef float S;

	// compute tanimoto similarity
	inline S tanimoto(FingerprintTanimoto *print) {
		int count_and = 0;
		int count_or = 0;
		int min, len;

		len = arrayLength();
		min = print->arrayLength();
		
		// handle different bit-length
		if (min <= len) {
			for (int i = min; i < len; i++) {
				count_or += cardWord(mArray[i]);
			}
		} else {
			for (int i = len; i < min; i++) {
				count_or += cardWord(print->mArray[i]);
			}
			min = len;
		}
		
		for (int i = 0; i < min; i++) {
			int a = mArray[i] & print->mArray[i];
			int o = mArray[i] | print->mArray[i];
			count_and += cardWord(a);
			count_or  += cardWord(o);
		}
		
		return ((S) count_and) / count_or;
	}
	
	// compute tanimoto estimation on hash-keys
	inline S tanimotoXOR(FingerprintTanimoto *print, int AB) {

		int xorCount  = cardWord(mHashArray[0] ^ print->mHashArray[0])
			      + cardWord(mHashArray[1] ^ print->mHashArray[1])
			      + cardWord(mHashArray[2] ^ print->mHashArray[2])
			      + cardWord(mHashArray[3] ^ print->mHashArray[3]);
		
		return ((S) (AB - xorCount)) / (AB + xorCount);
	}

	// compute lower bound for 1D-Grid
	static inline int lowerBound(S minTanimotoSimilarity, int card) {
		return (int) ceil(minTanimotoSimilarity * card);
	}

	// compute upper bound for 1D-Grid
	static inline int upperBound(S minTanimotoSimilarity, int card) {
                return (int) (1.0 / minTanimotoSimilarity * card) + 1;
	}
};

// specialized Fingerprint class for computing Hamming distance
class FingerprintHamming : public Fingerprint {
	public:

	inline FingerprintHamming(int length) : Fingerprint(length) {}
	inline FingerprintHamming(char * id, const char *str) : Fingerprint(id, str) {}
	inline FingerprintHamming(FingerprintHamming *print) : Fingerprint((Fingerprint*) print) {}
	inline int isEqual(Fingerprint *print) { return Fingerprint::isEqual((Fingerprint*) print); }
	inline void join(Fingerprint *print) { Fingerprint::join((Fingerprint*) print); }
	inline int intersectCard(Fingerprint *print) { return Fingerprint::intersectCard((Fingerprint*) print); }

	typedef int S;

	// compute hamming distance
	inline S hamming(FingerprintHamming *print) {
		int count_xor = 0;
		int min, len;

		len = arrayLength();
		min = print->arrayLength();
		
		// handle different bit-length
		if (min <= len) {
			for (int i = min; i < len; i++) {
				count_xor += cardWord(mArray[i]);
			}
		} else {
			for (int i = len; i < min; i++) {
				count_xor += cardWord(print->mArray[i]);
			}
			min = len;
		}
		
		for (int i = 0; i < min; i++) {
			count_xor += cardWord(mArray[i] ^ print->mArray[i]);
		}
		
		return count_xor;
	}
	
	// compute hamming distance estimation on hash-keys
	inline S hammingXOR(FingerprintHamming *print) {

		return (  cardWord(mHashArray[0] ^ print->mHashArray[0])
		 	+ cardWord(mHashArray[1] ^ print->mHashArray[1])
			+ cardWord(mHashArray[2] ^ print->mHashArray[2])
			+ cardWord(mHashArray[3] ^ print->mHashArray[3])
		       );
	}

	// compute lower bound for 1D-Grid
	static inline int lowerBound(S maxHammingDistance, int card) {
		return card - maxHammingDistance;
	}

	// compute upper bound for 1D-Grid
	static inline int upperBound(S maxHammingDistance, int card) {
		return card + maxHammingDistance;
	}
};
#endif
