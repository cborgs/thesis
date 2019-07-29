// FileSupport.h
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

#ifndef FILESUPPORT_H
#define FILESUPPORT_H

#include <stdio.h>
#include <string.h>

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include "Grid1D.h"

// maximal line size = length of ascii representation of fingerprint
#define STRSIZE 4000

// objects of template class FileSupport provide functionalities to
// read fingerprint files in csv like format

template<class T>
class FileSupport {
	private:

	// check, if a character is considered to be a white-space or seperator
	inline int isWS(char c) {
		return (
			(c == '"') ||
			(c == '\'') ||
			(c == ',') ||
			(c == ';') ||
			(c == ' ') ||
			(c == '\t')
		);
	}

	// check, if a character is considered to be end of line
	inline int isEOL(char c) {
		return (
			(c == 10) ||
			(c == 13) ||
			(c == 0)
		);
	}

	int parseLine(FILE *in, char *str, int64_t *idx1, int64_t *end1, int64_t *idx2, int64_t *end2);

	public:

	T* parseFingerprint(FILE *in, int64_t line);
	void parseAllFingerprints(const char *filename, T ***prints, int64_t *size, int *nBits);
};


// parse line from file, return number of parsed fields
template <class T>
int FileSupport<T>::parseLine(FILE *in, char *str, int64_t *idx1, int64_t *end1, int64_t *idx2, int64_t *end2) {
	// read line from file
	if (fgets(str, STRSIZE, in) == NULL) {
		return 0;
	}
	
	// parse line, find start of first field
	*idx1 = 0;
	while (isWS(str[*idx1]) && (*idx1 < STRSIZE-1)) {
		(*idx1)++;
	}

	// find end of first field
	*end1 = *idx1;
	while (!isWS(str[*end1]) && (*end1 < STRSIZE-1) && !isEOL(str[*end1])) {
		(*end1)++;
	}

	if (!isEOL(str[*end1])) {
		// find start of second field
		*idx2 = *end1;
		while (isWS(str[*idx2]) && (*idx2 < STRSIZE-1)) {
			(*idx2)++;
		}
		
		// find end of second field
		*end2 = *idx2;
		while (!isWS(str[*end2]) && (*end2 < STRSIZE-1) && !isEOL(str[*end2])) {
			(*end2)++;
		}
	} else {
		*idx2 = 0;
		*end2 = 0;
	}

	// cut field 1
	str[*end1] = 0;

	if (*idx2 != *end2) {
		// if there are two fields, cut field 2
		str[*end2] = 0;
		return 2;
	}

	return 1;
}

// read a fingerprint from file
template <class T>
T* FileSupport<T>::parseFingerprint(FILE *in, int64_t line) {
	int64_t fields;
	char *idStr;
	char str[STRSIZE];
	int64_t idx1, end1, idx2, end2;

	fields = parseLine(in, str, &idx1, &end1, &idx2, &end2);

	if (fields == 1) {
		// if there is only one field, use line as id
		idStr = new char[13];
		sprintf(idStr, "%012" PRId64, line + 1);
		// create fingerprint
		return new T(idStr, str+idx1);
	} else if (fields > 1) {
		// if there are two fields, use first string as id
		idStr = new char[end1-idx1+1];
		strcpy(idStr, str+idx1);
		// create fingerprint
		return new T(idStr, str+idx2);
	}

	return NULL;
}

// read fingerprint array from file
template <class T>
void FileSupport<T>::parseAllFingerprints(const char *filename, T ***prints, int64_t *size, int *nBits) {
	int64_t sizePrints;
	int maxLen;
	int len;
	char str[STRSIZE];
	FILE *in;
	T **readPrints;

	sizePrints = *size;

	// first pass: count number of prints in file
	if (sizePrints == 0) {
		in = fopen(filename, "r");

		if (in == NULL) {
			*size = 0;
			*nBits = 0;
			return;
		}

                while (fgets(str, STRSIZE, in)) {
			sizePrints++;
		}

		fclose(in);
	}

	// second pass: count number of prints in file
	maxLen = 0;

        in = fopen(filename, "r");

	if (in == NULL) {
		*size = 0;
		*nBits = 0;
		return;
	}

	// create array of Fingerprints
        readPrints = new T*[sizePrints];

	// read prints from file
        for (int64_t i = 0; i < sizePrints; i++) {
		// parse line from file
		readPrints[i] = parseFingerprint(in, i);

		if (readPrints[i] == NULL) {
			sizePrints = i;
			break;
                }
		
		// compute maximal length of Fingerprints
                len = readPrints[i]->getLength();

                if (len > maxLen) {
                        maxLen = len;
                }
        }

        fclose(in);

	*nBits = maxLen;
	*size = sizePrints;
	*prints = readPrints;
}
#endif
