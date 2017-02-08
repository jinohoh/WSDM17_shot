#ifndef _OPTION_
#define _OPTION_

#include <getopt.h>

#include "Macros.hpp"

using namespace std;

enum optionEnum {
	/* model selection */
	TENSORFACT, SHOWMODEL, 

	/* parameter setting */
	SIZE_THREAD, SIZE_RANK, SIZE_MODE, SIZE_ITERATION, SIZE_EPSILON,

	/* method flag */
	M_SCAN, M_PLAIN,

	/* Loss flag */
	L_MODEL, L_CORE, 

	/* Termination condition */
	TC_ITER, TC_EPSILON
};

class Option {
private:
	Option();

public:

	static int flagMode;

	// File I/O variables
	static char *tensorFileName;
	static char *dbFileName;
	static char *modelFileName;

	// Dataset stats
	static uint szMode;
	static uint *szDim;
	static uint *szRank;

	// Method alternatives
	static uint method;
	static PRECISION initVar;
	static uint maxIteration;
	static double epsilon;

	// System variables
	static uint szThread;
	static uint szBufVec;
	static int flagVerbose;
	static bool flagFlip;
	static optionEnum flagMethod;
	static optionEnum flagLoss;
	static optionEnum flagCondTerminal;
	
	static bool parse(const int argc, char** argv);

	~Option();

	static void cleanup(){
		if (szDim != nullptr)
			delete[] szDim;
		if (szRank != nullptr)
			delete[] szRank;
	}
};

#endif // #ifndef _OPTION_
