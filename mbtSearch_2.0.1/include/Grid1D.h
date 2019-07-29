// Grid1D.h
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

#ifndef GRID1D_H
#define GRID1D_H

#include <stdio.h>

#include "Fingerprint.h"
#include "FingerprintSorter.h"
#include "MultibitTree.h"
#include "ThreadPool.h"
#include "Tasks.h"
#include "FileSupport.h"

// Objects of template class Grid1D hold an array of instances of the template
// classes MultibitTree or UnionTree.
// In each tree all Fingerprints of the same cardinality are stored.
// Knowing the queries cardinality and the similarity threshold a search can
// be reduced on a relevant range of trees.
// The Grid1D also uses the template class ThreadPool for concurrently work
// on different different trees.
template<class TT>
class Grid1D {
	public:

	typedef TT T;
	typedef typename TT::P P;
	typedef typename TT::S S;
	
	private:

	P **mPrints;		// array of fingerprints
	T **mBuckets;		// array of search trees
	int mNBits;		// maximal size of fingerprints
	int64_t mSize;		// size of fingerprint-array used by mBuckets
	
	public:

	// constructor:
	//
	// read Fingerprints from file and sort them by cardinality and
	// create a search tree for each cardinality
	//
	// prints       : pointer on Fingerprint array
	// size         : size of <prints>
	// threads      : number of parallel threads passed to ThreadPool
	// leafLimit    : leaf limit parameter passed to all MultibitTrees
	Grid1D(const char *filename, int64_t size, int threads, int leafLimit) {
		FingerprintSorter<P> fpSorter;
		FileSupport<P> fileSupport;
		TaskCreateTree<T> task;

		ThreadPool< TaskCreateTree<T> > *workerPool = new ThreadPool< TaskCreateTree<T> >(threads);

		// read prints from file
		mSize = size;
		fileSupport.parseAllFingerprints(filename, &mPrints, &mSize, &mNBits);

		mBuckets = new T*[mNBits + 1];

		// sort prints into clusters by cardinality
		fpSorter.cluster(mPrints, mSize, mNBits);

		// create a MultibitTree for each cardinality cluster
		for (int i = 0; i < (mNBits + 1); i++) {
			if (fpSorter.getSize(i) > 0) {
				task.set(&mBuckets[i], mPrints, fpSorter.getStart(i), fpSorter.getEnd(i), mNBits, i, leafLimit);
				workerPool->dispatchTask(task);
			} else {
				mBuckets[i] = NULL;
			}
		}

		// wait for running threads
		workerPool->wait();

		delete workerPool;
	}

	// delete grid data
	~Grid1D() {
		for (int i = 0; i < mNBits; i++) {
			if (mBuckets[i]) {
				delete mBuckets[i];
			}
		}

		delete[] mBuckets;
	}

	// get size of print array
	inline int64_t getSize() {
		return mSize;
	}

	// get maximal print size
	inline int getNBits() {
		return mNBits;
	}

	// search query prints from input file in tree
	void searchFile(FILE *inputFile, S filter, QueryResult *queryResult, int threads) {
		P *queryPrint;
		int64_t i;
		int min, max, card;
		FileSupport<P> fileSupport;
		TaskSearchTreeRange<T> task;
		ThreadPool< TaskSearchTreeRange<T> > *workerPool = new ThreadPool< TaskSearchTreeRange<T> >(threads);

		if (inputFile != NULL) {
			i = 0;
			while (1) {
				// for each line parse fingerprint
				queryPrint = fileSupport.parseFingerprint(inputFile, i);

				if (queryPrint == NULL) {
					break;
				}

				// search only in buckets with suitable cardinality
				card = queryPrint->cardinality();
				min = MAX(P::lowerBound(filter, card), 0);
				max = MIN(P::upperBound(filter, card), mNBits + 1);

				// search print in grid
				task.set(mBuckets, min, max, queryResult, queryPrint, card, filter);
				workerPool->dispatchTask(task);

				i++;
			}

			queryResult->setSizeLastSearch(i);

			// wait for running threads
			workerPool->wait();
		}

		delete workerPool;
	}

	// search query prints from input file in tree with symdex pre-processing
	void searchFileSymdex(const char *filename, S filter, QueryResult *queryResult, int threads) {
		FingerprintSorter<P> fpSorter;
		P **queryPrints;
		int64_t sizePrints = 0;
		int nBits;
		int min, max;
		FileSupport<P> fileSupport;
		TaskSearchTreeSymdex<T> task;
		ThreadPool< TaskSearchTreeSymdex<T> > *workerPool = new ThreadPool< TaskSearchTreeSymdex<T> >(threads);

		// read query prints from file
		fileSupport.parseAllFingerprints(filename, &queryPrints, &sizePrints, &nBits);

		// sort prints into clusters by cardinality
		fpSorter.cluster(queryPrints, sizePrints, nBits);

		// for each cluster call asynchonous search-method
		for (int i = 0; i < (mNBits + 1); i++) {
			if (fpSorter.getSize(i) > 0) {
				// search only in buckets with suitable cardinality
				min = MAX(P::lowerBound(filter, i), 0);
				max = MIN(P::upperBound(filter, i), mNBits + 1);

				task.set(mBuckets, min, max, queryResult, queryPrints, fpSorter.getStart(i), fpSorter.getEnd(i), i, filter);
				workerPool->dispatchTask(task);
			}
		}

		queryResult->setSizeLastSearch(sizePrints);

		// wait for running threads
		workerPool->wait();

		delete workerPool;
	}
};
#endif
