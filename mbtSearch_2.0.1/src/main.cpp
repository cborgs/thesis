// main.cpp
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

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "Fingerprint.h"
#include "MultibitTree.h"
#include "UnionTree.h"
#include "Grid1D.h"
#include "Timer.h"

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

// types of required algorithms
typedef Grid1D< MultibitTree<FingerprintTanimoto> > GridMbtTan;
typedef Grid1D< MultibitTree<FingerprintHamming> > GridMbtHam;
typedef Grid1D< UnionTree<FingerprintTanimoto> > GridUnionTan;
typedef Grid1D< UnionTree<FingerprintHamming> > GridUnionHam;

// print usage information
void printUsage() {
	fprintf(stderr, "Usage: mbtSearch [ -i <query_prints.csv> ] -m <min-similarity>\n");
        fprintf(stderr, "  [ -o <results.csv> ] [ -d <csv-delimiter> ] [ -t <threads> ] [ -s <size> ]\n");
        fprintf(stderr, "  [ -l <leaf-limit> ] [ -p ] [ -a <algorithm> ] [ -v | -V ] [ -h ]\n");
        fprintf(stderr, "  <tree_prints.csv>\n");
}

// print user help
void printHelp() {
	printUsage();
	fprintf(stderr, "\n");
	fprintf(stderr, "  <tree_prints.csv>\n");
	fprintf(stderr, "          Mandatory parameter specifying the input csv-file containing\n");
	fprintf(stderr, "          the fingerprints to be loaded into the tree.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "  -i <query_prints.csv>\n");
	fprintf(stderr, "          Optional parameter specifying the input csv-file containing\n");
	fprintf(stderr, "          the fingerprints to search for. If this option is omitted the\n");
        fprintf(stderr, "          query fingerprints will be read from STDIN.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "  -m <min-similarity>\n");
	fprintf(stderr, "          Mandatory numeric value giving the lower bound of the similarity-\n");
        fprintf(stderr, "          coefficient to search for.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "  -o <results.csv>\n");
	fprintf(stderr, "          Optional parameter specifying the output csv-file containing\n");
	fprintf(stderr, "          the query results. If this option is omitted the results\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "  -d <delimiter>\n");
	fprintf(stderr, "          Optional parameter specifying the column delimiter for the\n");
        fprintf(stderr, "          output csv-file. If this option is omitted it defaults to ','.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "  -t <threads>\n");
	fprintf(stderr, "          Optional parameter specifying the number of parallel threads\n");
        fprintf(stderr, "          that shall be used to construct the search tree and perform\n");
        fprintf(stderr, "          the search within it. If this option is omitted it defaults to 1.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "  -s <size>\n");
	fprintf(stderr, "          Optional parameter specifying the maximal number of fingerprints\n");
        fprintf(stderr, "          that shall be loaded into the search tree\n");
        fprintf(stderr, "          If this option is omitted it defaults to 0 (unlimited).\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "  -l <leaf-limit>\n");
	fprintf(stderr, "          Optional parameter specifying the maximum number of fingerprints\n");
        fprintf(stderr, "          for which no further sub-tree shall be calculated\n");
        fprintf(stderr, "          If this option is omitted it defaults to 8.\n");
	fprintf(stderr, "\n");
        fprintf(stderr, "  -p      Optional switch to activate the symdex pre-processing.\n");
	fprintf(stderr, "\n");
        fprintf(stderr, "  -a <algorithm>\n");
        fprintf(stderr, "          Optional parameter specifying the search algorithm to be used.\n");
        fprintf(stderr, "          Possible values are:\n");
	fprintf(stderr, "\n");
        fprintf(stderr, "            mtan : multibit tree using tanimoto similarity (default)\n");
        fprintf(stderr, "            mham : multibit tree using hamming similarity\n");
        fprintf(stderr, "            utan : union tree using tanimoto similarity\n");
        fprintf(stderr, "            uham : union tree using hamming similarity\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "  -v      In verbose mode additional information will be printed out.\n");
	fprintf(stderr, "\n");
        fprintf(stderr, "  -V      In extended verbose mode additional statistic information will\n");
        fprintf(stderr, "          be collectd and printed. This slows down multi-threaded searches.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "  -h      Show this help.\n");
}

void printStatistics(QueryResult *query) {
	double values[3];
	double percents[3];

	query->getStatistics(values, percents);

	fprintf(stderr, "\nStatistics:\n\n");

	fprintf(stderr, "Checkpoint             Count   Percentage\n");
	fprintf(stderr, "XOR-Hash    %16.0f       %6.2f\n", values[0], percents[0]);
	fprintf(stderr, "Similarity  %16.0f       %6.2f\n", values[1], percents[1]);
	fprintf(stderr, "Total       %16.0f       %6.2f\n", values[2], percents[2]);

	fprintf(stderr, "\n");
}

// template function to load a search tree for a given algorithm
template <class T>
T* loadTree(const char * treeFileName, int64_t size, int threads, int leafLimit, int verbose) {
	Timer timer;
	T *grid;

	if (verbose >= 1) {
		fprintf(stderr, "Loading tree data from file '%s' ...\n", treeFileName);
		timer.start();
	}

	// load tree
	grid = new T(treeFileName, size, threads, leafLimit);
	size = grid->getSize();

	if (size <= 0) {
		fprintf(stderr, "Error: No data found while trying to load tree data!");
		exit(1);
	}

	// print performance information
	if (verbose >= 1) {
		fprintf(stderr, "%" PRId64 " tree loaded in %.3f seconds\n", size, (float)timer.stop() / 1000);
	}

	return grid;
}

// template function to start search with a given algorithm
template <class T>
void searchTree(T *grid, const char *inputFileName, int inputFileFlag, FILE *resultFile, int preproc, float minSimilarity, const char * delimiter, int verbose, int threads) {
	Timer timer;
	FILE *inputFile = stdin;
	int64_t treeSize = 0;

	if (verbose >= 1) {
		fprintf(stderr, "Searching prints ...\n");
		timer.start();
	}

	// initialize query result file
	if (verbose >= 2) {
		// set grid size for statistics
		treeSize = grid->getSize();
	}

        QueryResult queryResult(0, resultFile, delimiter, treeSize);

	// search tree for query prints
	if (preproc) {
		grid->searchFileSymdex(inputFileName, minSimilarity, &queryResult, threads);
	} else {
		if (inputFileFlag) {
			// if specified, open input file
			inputFile = fopen(inputFileName, "r");

			if (inputFile == NULL) {
				fprintf(stderr, "Error: Could not open input-file '%s'!\n", inputFileName);
				exit(1);
			}
		}

		grid->searchFile(inputFile, minSimilarity, &queryResult, threads);

		if (inputFileFlag) {
			fclose(inputFile);
		}
	}

	// print performance information
	if (verbose >= 1) {
		fprintf(stderr, "prints searched in %.3f seconds\n", (float)timer.stop() / 1000);
	}

	if (verbose >= 2) {
		printStatistics(&queryResult);
	}
	
	delete grid;
}

// command line main
int main (int argc, char **argv) {
	int c;
	const char *inputFileName = "STDIN";
	const char *treeFileName = NULL;
	const char *resultFileName = "STDOUT";
	FILE *resultFile = stdout;
	int resultFileFlag = 0;
	int inputFileFlag = 0;
	const char *delimiter = ",";
	int threads = 1;
	const char *sizeText = "unlimited";
	int64_t size = 0;
	int leafLimit = 8;
	int verbose = 0;
	int preproc = 0;
	float minSimilarity = 0.0;
	const char *algorithm = "mtan";

	// evaluate command line options
	while ((c = getopt(argc, argv, "i:m:o:d:t:s:l:a:pvVh")) != -1) {
		switch (c) {
			case 'i':
				inputFileName = optarg;
				inputFileFlag = 1;
				break;
			case 'm':
				sscanf(optarg, "%f", &minSimilarity);
				break;
			case 'o':
				resultFileName = optarg;
				resultFileFlag = 1;
				break;
			case 'd':
				delimiter = optarg;
				break;
			case 't':
				threads = atoi(optarg);
				break;
			case 's':
				sizeText = optarg;
				size = atoll(optarg);
				break;
			case 'l':
				leafLimit = atoi(optarg);
				break;
			case 'a':
				algorithm = optarg;
				break;
			case 'p':
				preproc = 1;
				break;
			case 'v':
				verbose = 1;
				break;
			case 'V':
				verbose = 2;
				break;
			case 'h':
				printHelp();
				exit(0);
			default:
				printUsage();
				exit(1);
		}
	}

	if (optind == argc - 1) {
		treeFileName = argv[optind];
	} else {
		printUsage();
		exit(1);
	}

	if (verbose >= 1) {
		fprintf(stderr, "input-file = '%s', tree-file = '%s', min-similarity = '%f'\n", inputFileName, treeFileName, minSimilarity);
		fprintf(stderr, "result-file = '%s', delimiter = '%s'\n", resultFileName, delimiter);
		fprintf(stderr, "worker-threads = %d, maximal number of prints in tree = '%s', leaf-limit = %d\n", threads, sizeText, leafLimit);
		fprintf(stderr, "symdex pre-processing = ");

		if (preproc) {
			fprintf(stderr, "on");
		} else {
			fprintf(stderr, "off");
		}

		fprintf(stderr, ", algorithm = %s\n\n", algorithm);
	}

        if (resultFileFlag) {
                // if specified, open result file
                resultFile = fopen(resultFileName, "w");

		if (resultFile == NULL) {
			fprintf(stderr, "Error: Could not create result-file '%s'!\n", resultFileName);
			exit(1);
		}

                // print column headers
                fprintf(resultFile, "query%sfingerprint%ssimilarity\n", delimiter, delimiter);
        }

	if ((minSimilarity <= 0.0) || (minSimilarity > 1.0)) {
		fprintf(stderr, "Error: min-similarity must be greater than 0.0 and less than or equal 1.0!\n");
		exit(1);
	}

	if (threads < 1) {
		fprintf(stderr, "Error: the number of worker-threads must be at least 1!\n");
		exit(1);
	}

	if (preproc && !inputFileFlag) {
		fprintf(stderr, "Error: An input file has to be provided if symdex pre-processing is switched on!\n");
		exit(1);
	}

        // initialize Fingerprint data structure (cardinality-map)
        Fingerprint::init();

	// load tree and search queries by the selected algorithm
	if (strcmp(algorithm, "mtan") == 0) {
		GridMbtTan *grid = loadTree<GridMbtTan>(treeFileName, size, threads, leafLimit, verbose);
		searchTree<GridMbtTan>(grid, inputFileName, inputFileFlag, resultFile, preproc, minSimilarity, delimiter, verbose, threads);
	} else if (strcmp(algorithm, "mham") == 0) {
		GridMbtHam *grid = loadTree<GridMbtHam>(treeFileName, size, threads, leafLimit, verbose);
		searchTree<GridMbtHam>(grid, inputFileName, inputFileFlag, resultFile, preproc, (1.0 - minSimilarity) * grid->getNBits(), delimiter, verbose, threads);
	} else if (strcmp(algorithm, "utan") == 0) {
		GridUnionTan *grid = loadTree<GridUnionTan>(treeFileName, size, threads, leafLimit, verbose);
		searchTree<GridUnionTan>(grid, inputFileName, inputFileFlag, resultFile, preproc, minSimilarity, delimiter, verbose, threads);
	} else if (strcmp(algorithm, "uham") == 0) {
		GridUnionHam *grid = loadTree<GridUnionHam>(treeFileName, size, threads, leafLimit, verbose);
		searchTree<GridUnionHam>(grid, inputFileName, inputFileFlag, resultFile, preproc, (1.0 - minSimilarity) * grid->getNBits(), delimiter, verbose, threads);
	} else {
		fprintf(stderr, "Error: Unknown algorithm '%s'!\n", algorithm);
		exit(1);
	}

	if (verbose >= 1) {
		fprintf(stderr, "Done.\n");
	}

	if (resultFileFlag) {
		fclose(resultFile);
	}

	exit(0);
}
