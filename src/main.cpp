#include <stdio.h>
#include "Macros.hpp"
#include "baseStruct.hpp"
#include "TensorFile.hpp"
#include "SolverTucker.hpp"
#include "RevEigSolver.hpp"
#include "Option.hpp"
#include <omp.h> // omp related functions.

int main(int argc, char** argv){

	// init.
	_srand(time(NULL));
	omp_set_dynamic(0);

	Option::parse(argc, argv);
	if (nullptr == Option::tensorFileName || nullptr == Option::modelFileName)
	{
		cout << "Not enough parameters" << endl;
		exit(0);
	} 
	else 
	{
		cout << "TensorFile: " << Option::tensorFileName << endl;
		cout << "ModelFile: " << Option::modelFileName << endl;
	}

	// build db
	TensorIndexWithMultipleFile<Record> index;
	index.buildDBFromTextFile(Option::tensorFileName, Option::dbFileName);

	// init SolverTransTucker, and run iterations.
	SolverTuckerFact<Record> solver(&index);
	solver.run();
	
	// Clean-up
	Option::cleanup();

	return 0;
}
