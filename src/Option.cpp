#include "Option.hpp"

#include <cstring>
#include <iostream>
#include <vector>

int Option::flagMode;

char* Option::tensorFileName;
char* Option::dbFileName;
char* Option::modelFileName;
uint Option::szMode;
uint* Option::szDim;
uint* Option::szRank;

uint Option::method;
PRECISION Option::initVar;
uint Option::maxIteration;
double Option::epsilon;

uint Option::szThread;
uint Option::szBufVec;
int Option::flagVerbose;
bool Option::flagFlip;
optionEnum Option::flagMethod;
optionEnum Option::flagLoss;
optionEnum Option::flagCondTerminal;


const char short_options[] = "i:l:s:t:u:";

static struct option long_options[] = {
	{"showmodel", no_argument, &Option::flagMode, SHOWMODEL},
	{"verbose", no_argument, &Option::flagVerbose, 1},
	{"plain", no_argument, (int*)&Option::flagMethod, M_PLAIN},
	{"scan", no_argument, (int*)&Option::flagMethod, M_SCAN},
	{"use-modelloss", no_argument, (int*)&Option::flagLoss, L_MODEL},
	{"use-coreloss", no_argument, (int*)&Option::flagLoss, L_CORE},
	{"use-iteration", no_argument, (int*)&Option::flagCondTerminal, TC_ITER}, 
	{"use-epsilon", no_argument, (int*)&Option::flagCondTerminal, TC_EPSILON}, 
	{"size-mode",  required_argument, 0, SIZE_MODE},
	{"size-thread",  required_argument, 0, SIZE_THREAD},
	{"size-rank",  required_argument, 0, SIZE_RANK},
	{"iteration", required_argument, 0, SIZE_ITERATION}, 
	{"epsilon", required_argument, 0, SIZE_EPSILON},
	{0, 0, 0, 0}
};

bool Option::parse(const int argc, char** argv){

	int rankDefault = 4;
	int opt;

	/* Default parameter setting */
	Option::flagMode = TENSORFACT;
	
	Option::tensorFileName = nullptr;
	Option::dbFileName = nullptr;
	Option::modelFileName = nullptr;
	Option::szMode = 3;
	Option::szDim = nullptr;
	Option::szRank = nullptr;
	
	Option::method = 0;
	Option::initVar = 0.1;
	Option::maxIteration = 1;
	Option::epsilon = 0.01;
	
	Option::szThread = 4;
	Option::szBufVec = 500000;
	Option::flagVerbose = 0;
	Option::flagFlip = true;
	Option::flagMethod = M_SCAN;
	Option::flagLoss = L_MODEL;
	Option::flagCondTerminal = TC_ITER;

	int option_index =0;
	while (1)
	{
		opt = getopt_long(argc, argv, short_options, long_options, &option_index);
		if (-1 == opt)
			break;

		switch (opt)
		{
			case SIZE_MODE:
				Option::szMode = atoi(optarg);
				break;
			case SIZE_THREAD:
				Option::szThread = atoi(optarg);
				break;
			case SIZE_RANK: 
				{
					// Processing rank values
					vector<int> rankList;
					char* pch = strtok(optarg, ":");
					while (nullptr != pch){
						rankList.push_back(atoi(pch));
						pch = strtok(nullptr, ":");
					}
					Option::szRank = new uint[rankList.size()];
					for (uint i = 0; i < rankList.size(); i++)
						Option::szRank[i] = rankList[i];
					Option::szMode = rankList.size();
				}
				break;
			case SIZE_ITERATION:
				Option::maxIteration = atoi(optarg);
				break;
			case SIZE_EPSILON:
				Option::epsilon = atof(optarg);
				break;
		}
	}

	if (optind < argc)
	{
		switch (Option::flagMode){
			case TENSORFACT:
				if (argc - optind < 2)
					exit(0);
				
				Option::tensorFileName = argv[optind++];
				Option::dbFileName = new char[strlen(Option::tensorFileName) + 10];
				strcpy(Option::dbFileName, Option::tensorFileName);
				strcat(Option::dbFileName, ".db");
				Option::modelFileName = argv[optind++];
				break;
			case SHOWMODEL:
				Option::modelFileName = argv[optind++];
		}
	}

	if (nullptr == Option::szRank){
		Option::szRank = new uint[Option::szMode];
		for (uint i = 0; i < Option::szMode; i++)
			Option::szRank[i] = rankDefault;
	}

	Option::szDim = new uint[Option::szMode];
	for (uint i = 0; i < Option::szMode; i++)
		Option::szDim[i] = 0;

	return true;
}

Option::~Option(){
	if (nullptr != Option::dbFileName)
		delete[] Option::dbFileName;
	if (nullptr != Option::szRank)
		delete[] Option::szRank;
	if (nullptr != Option::szDim)
		delete[] Option::szDim;
}
