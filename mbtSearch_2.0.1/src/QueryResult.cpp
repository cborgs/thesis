// QueryResult.cpp
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

#include "QueryResult.h"

// constructor
QueryResult::QueryResult(int sort, FILE *resultFile, const char *seperator, int64_t sizeTree) {
	mSize = 0;
	mRootNode = NULL;
	mSort = sort;
	mResultFile = resultFile;
	mSeperator = seperator;
	mCntXOR = 0;
	mCntSim = 0;
	mSizeTree = sizeTree;
	mSizeSearch = 0;

	pthread_mutex_init(&mAddMutex, NULL);
	pthread_mutex_init(&mCntXORMutex, NULL);
	pthread_mutex_init(&mCntSimMutex, NULL);
}

// destructor
QueryResult::~QueryResult() {
	QueryResultNode *node, *nextNode;
	if (mSort) {
		deleteNode(mRootNode);		// discard the whole result set recursivly
	} else {
		// discard the whole result iteratively for large scale results
		node = mRootNode;
		while (node != NULL) {
			if (node->mQueryId != NULL) {
				delete[] node->mQueryId;
			}
			nextNode = node->mRight;
			delete node;
			node = nextNode;
		}
	}
}

// delete a node and all its subnodes recursively
void QueryResult::deleteNode(QueryResultNode *node) {
	if (node != NULL) {
		if (node->mQueryId != NULL) {
			delete[] node->mQueryId;
		}
		deleteNode(node->mLeft);
		deleteNode(node->mRight);
		delete node;
	}
}

// add a Fingerprint and the corresponding similarity coefficient to the query result
void QueryResult::add(char *queryId, Fingerprint *print, float similarity) {
	// lock mutex
	pthread_mutex_lock(&mAddMutex);	

	// check if result file is specified
	if (mResultFile == NULL) {
		// if not collect results in memory
		QueryResultNode *newNode;
		QueryResultNode **nextNodePtr;
	
		// create new node
		newNode = new QueryResultNode;

		// copy query id, print pointer and similarity
		if (queryId != NULL) {
			newNode->mQueryId = new char[strlen(queryId) + 1];
			strcpy(newNode->mQueryId, queryId);
		} else {
			newNode->mQueryId = NULL;
		}

		newNode->mPrint = print;
		newNode->mSimilarity = similarity;
		newNode->mLeft = NULL;
		newNode->mRight = NULL;

		if (mSort) {
			// if query results have to be sorted, find the right location
			// to insert the new node
			nextNodePtr = &mRootNode;

			while (*nextNodePtr != NULL) {
				if ((*nextNodePtr)->mSimilarity >= similarity) {
					nextNodePtr = &((*nextNodePtr)->mRight);
				} else {
					nextNodePtr = &((*nextNodePtr)->mLeft);
				}
			}

			// insert new node
			*nextNodePtr = newNode;
		} else {
			// if query results need not to be sorted
			// insert at end of list
			if (mRootNode == NULL) {
				mRootNode = newNode;
			} else {
				mRootNode->mLeft->mRight = newNode;
			}

			mRootNode->mLeft = newNode; // use mLeft as pointer to end node
		}
	} else {
		// write results to file
		fprintf(mResultFile, "%s%s%s%s%.7f\n", queryId, mSeperator, print->getId(), mSeperator, similarity);
	}
	
	mSize++;

	// unlock mutex
	pthread_mutex_unlock(&mAddMutex);
}

// get statistics
void QueryResult::getStatistics(double *valuesPtr, double *percentsPtr) {
	valuesPtr[0] = (double)mCntXOR;
	valuesPtr[1] = (double)mCntSim;
	valuesPtr[2] = (double)(mSizeTree * mSizeSearch);

	percentsPtr[0] = (double)mCntXOR / (mSizeTree * mSizeSearch) * 100;
	percentsPtr[1] = (double)mCntSim / (mSizeTree * mSizeSearch) * 100;
	percentsPtr[2] = 100.0;
}

