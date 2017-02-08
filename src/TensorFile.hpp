#ifndef _TENSORFILE_
#define _TENSORFILE_

#include "Macros.hpp"
#include "baseStruct.hpp"
#include "Option.hpp"
//#include "Matrix.hpp"

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <unistd.h> // open(), close(), lseek()


#include <fstream> // ifstream fp in TensorTextFileReader
#include <vector>  // std::vector
#include <algorithm> // std::make_heap, std::sort
#include <cstring> // strcpy

using namespace std;

template<typename Tuple>
class TensorTextFileReader {
private:
	ifstream fp;
	TensorTextFileReader() {
		fp = nullptr;
	};
public:
	TensorTextFileReader(const char *filename){
		fp.open(filename, std::ifstream::in);
		return;
	}
	~TensorTextFileReader(){
		if (fp.is_open())
			fp.close();
	}
	bool getRecord(RecordID rid, Tuple& oRecord);
	bool getNextLine(Tuple& oRecord){
		if (fp.eof())
			return false;

		for (uint i = 0; i < Tuple::szMode; i++)
			fp >> oRecord.dim[i];
		fp >> oRecord.value;

		if (fp.eof())
			return false;
		return true;
	}

	bool rewind(){ 
		fp.clear();
		fp.seekg(ios_base::beg);
		return true;
	}
};

template<typename Tuple>
class TensorBinaryFileHandler {
private:
	int fp;
	TensorBinaryFileHandler(){
		fp = -1;
	};
	char filename[256];
	char *buf;

public:
	static const bool HARD_CLEANUP = true; 	// cleanup with file removal.
	static const bool SOFT_CLEANUP = false;	// cleanup w.o. file removal.

	TensorBinaryFileHandler(const char *inFilename, int flag = 0){
		strcpy(filename, inFilename);
		fp = open(filename, O_RDWR | O_CREAT | flag, S_IRUSR | S_IWUSR);
		buf = new char[Tuple::getSize()];
	}
	~TensorBinaryFileHandler(){
		delete[] buf;
	}
	void getRecord(RecordID rid, Tuple& oRecord);
	void appendRecord(Tuple& iRecord){
		iRecord.encode(buf);
		int size = write(fp, buf, Tuple::getSize());
		assert(size == Tuple::getSize());
		return;
	}

	void showRecord(){
		int seekp = lseek(fp, 0, SEEK_CUR);
		lseek(fp, 0, SEEK_SET);

		Tuple temp;
		while (getNext(temp)){
			for (int i = 0; i < Tuple::szMode; i++)
				cout << temp.dim[i] << "\t";
			cout << temp.value << endl;
		}
		lseek(fp, seekp, SEEK_SET);
	}

	bool getNext(Tuple& oRecord){
		int size = read(fp, buf, Tuple::getSize());
		if (0 == size)
			return false;
		assert(Tuple::getSize() == size);
		oRecord.decode(buf);
		return true;
	}

	bool rewind(){
		return (0 == lseek(fp, 0, SEEK_SET));
	}

	void cleanup(bool hardness = SOFT_CLEANUP){
		close(fp);
		if (HARD_CLEANUP == hardness)
			remove(filename);
		fp = -1;
	}
};

template<typename Tuple>
class TensorIndex
{
private:
	char curMode;

public:
	class Scanner {
		static uint sCurMode;
	public:
		virtual bool getNext(Tuple&) = 0;
		static void setMode(char newMode){ sCurMode = newMode;}
		virtual ~Scanner(){}
	};
	virtual ~TensorIndex(){}
	inline void setMode(char iMode) {
		curMode  = iMode;
		//Scanner::setMode(iMode);
	}
	virtual void buildDBFromTextFile(const char *txtfilename, const char* dbfilename)=0;
};


template<typename Tuple>
class TensorMemHashIndex : public TensorIndex<Tuple> // NOTE. this is memory based structure.
{
private:
	PRECISION 	*records;
	long uint 	nnz;

public:
	class Scanner: public TensorIndex<Tuple>::Scanner {
		static char sCursor;
		char mPos;
		char mLength;
	public:
		Scanner(){
			mPos = 0;
			mLength = 0;
		}

		bool getNext (Tuple& oRecord){
			return false;
		}
	};

	void buildDBFromTextFile(const char *txtfilename, const char *dbfilename);
};

template<typename Tuple>
class TensorIndexWithMultipleFile : public TensorIndex<Tuple> 
{	// NOTE. this is a class for maintaining multiple files each of which is 
	// 1) sorted according to the different mode and 2) divided for parallelization
	uint *sortOrder;
	char dbname[256];
	TensorBinaryFileHandler<Tuple>** arrDBSeg;

	struct TuplefpPair { // This is structure for using heap.
		Tuple key;
		TensorBinaryFileHandler<Tuple> *hfp;
		TuplefpPair(){}
		TuplefpPair(Tuple inkey, TensorBinaryFileHandler<Tuple> *inhfp): key(inkey), hfp(inhfp){}
		TuplefpPair& operator= (const Tuple& rhs){
			key = rhs;
		}
		void init(const Tuple& inKey, TensorBinaryFileHandler<Tuple> *inhfp){
			key = inKey;
			hfp = inhfp;
		}
	};

public:
	class Scanner: public TensorIndex<Tuple>::Scanner {
	private:
		char 								id;
		static uint 						sCurMode;
		TensorIndexWithMultipleFile<Tuple> 	*parent;
		TensorBinaryFileHandler<Tuple> 		*fptr;
		static char 						snScanner; 
							   // This variable might be shared regardless of the value of T 
						      //  --> Assume that there is only one type of tensor during execution.
	public:
		Scanner(TensorIndex<Tuple> *inParentPtr){
			parent = (TensorIndexWithMultipleFile<Tuple>*) inParentPtr;
			id = __sync_fetch_and_add(&snScanner, 1);
			assert(id < Option::szThread);

			fptr = parent->arrDBSeg[sCurMode*Option::szThread + id];
		}
		~Scanner(){}

		void reinit(){
			assert(id < Option::szThread);
			fptr = parent->arrDBSeg[sCurMode*Option::szThread + id];
			fptr->rewind();
			return;
		}

		static void setMode(char newMode){
			sCurMode = newMode;
			snScanner = 0;
			system("echo 3 > /proc/sys/vm/drop_caches");
			return;
		}

		bool getNext(Tuple& oRecord){
			return fptr->getNext(oRecord);
		}

		bool rewind(){
			return fptr->rewind();
		}
	};

	TensorIndexWithMultipleFile(){
		sortOrder = new uint[Tuple::szMode];
		for (int i = 0; i < Tuple::szMode; i++)
			sortOrder[i] = i;
		if (0 == Tuple::szMode){
			sortOrder = nullptr;
			delete[] sortOrder;
		}
	}

	~TensorIndexWithMultipleFile(){
		cleanup();
	}

	void cleanup(){
		if (sortOrder != nullptr){
			delete[] sortOrder;
			sortOrder = nullptr;
		}

		if (arrDBSeg != nullptr){
			for (int i = 0; i < Tuple::szMode * Option::szThread; i++){
				arrDBSeg[i]->cleanup(TensorBinaryFileHandler<Tuple>::HARD_CLEANUP);
				delete arrDBSeg[i];
			}
			delete[] arrDBSeg;
			arrDBSeg = nullptr;
		}
	}

	static bool cmp(const Tuple& lhs, const Tuple& rhs, uint *sortOrder){
		for (int i = 0; i < Option::szMode; i++)
			if (lhs.dim[sortOrder[i]] < rhs.dim[sortOrder[i]])
				return true;
			else if (lhs.dim[sortOrder[i]] > rhs.dim[sortOrder[i]])
				return false;
		return false;
	}

	void buildDBFromTextFile(const char *txtfilename, const char *dbfilename);

	TensorBinaryFileHandler<Tuple>* dumpRecordIntoFile(const char*filename, std::vector<Tuple>& recList){
			TensorBinaryFileHandler<Tuple> *hSegment = new TensorBinaryFileHandler<Tuple>(filename, O_TRUNC);
			for (int i = 0; i < recList.size(); i++){
				hSegment->appendRecord(recList[i]);
			}
			hSegment->rewind();
			return hSegment;
	}
};

template<typename Tuple>
void TensorIndexWithMultipleFile<Tuple>::buildDBFromTextFile(const char *txtfilename, const char *dbfilename){

	strcpy(dbname, dbfilename);

	char tempFileName[256];
	char splitFileName[256];
	char buf[256];

	// BUILD a temporal binary files from the text file 
	TensorTextFileReader<Tuple> txtReader(txtfilename);
	arrDBSeg = new TensorBinaryFileHandler<Tuple>*[Tuple::szMode * Option::szThread];
	std::vector<TensorBinaryFileHandler<Tuple>*> listSortedTempFile;
	Tuple temp;

	std::vector<Tuple> tempVector;
	std::vector<TuplefpPair> heap;
	tempVector.reserve(Option::szBufVec);
	uint sizeSegmentTempFile = 0;

	for (int mode = 0; mode < Tuple::szMode; mode++){
		// Re-initialize the variables for an iterations.
		tempVector.resize(0);
		listSortedTempFile.resize(0);
		heap.resize(0);
		sizeSegmentTempFile = 0;
		txtReader.rewind();
		
		sortOrder[0] = mode;
		int j = 1;
		for (int i = 0; i < Tuple::szMode; i++)
			if (i != mode)
				sortOrder[j++] = i;

		/*
		 * Scanning records from txt file, and construct sorted temporal files for K-way merges.
		 */
		while (txtReader.getNextLine(temp)){
			tempVector.push_back(temp);
			assert(temp.dim[mode] >= 0);
			if (temp.dim[mode] >= Option::szDim[mode])
				Option::szDim[mode] = temp.dim[mode] + 1;
			if (tempVector.size() >= Option::szBufVec){
				// 1) Sort the record in memory, 2) dump it to temporal file, and 
				// 3) add the file desriptor for the new temporal fileto listSortedTempFile.
				std::sort(tempVector.begin(), tempVector.end(), 
						[this](const Tuple &a, const Tuple &b){ 
							return TensorIndexWithMultipleFile<Tuple>::cmp(a, b, this->sortOrder);});

				strcpy(tempFileName, txtfilename);
				sprintf(buf, ".seg%03d", sizeSegmentTempFile++);
				strcat(tempFileName, buf);
				listSortedTempFile.push_back(dumpRecordIntoFile(tempFileName, tempVector));

				// Cleanup the vector for next iteration
				tempVector.resize(0);
			}
		}

		// Do sorting for the remaining records in the tempVector.
		if (tempVector.size() > 0){
				std::sort(tempVector.begin(), tempVector.end(), 
						[this](const Tuple &a, const Tuple &b){ 
							return TensorIndexWithMultipleFile<Tuple>::cmp(a, b, this->sortOrder);});
			strcpy(tempFileName, txtfilename);
			sprintf(buf, ".seg%03d", sizeSegmentTempFile++);
			strcat(tempFileName, buf);
			listSortedTempFile.push_back(dumpRecordIntoFile(tempFileName, tempVector));
		}

		/*
		 * K-way merge sort, and distribute records to splitted files in binary format.
		 */
		// File initialization
		for (int tid = 0; tid < Option::szThread; tid++){ // File initialization
			sprintf(buf, ".m%d.t%d", mode, tid);
			strcpy(splitFileName, dbfilename);
			strcat(splitFileName, buf);
			TensorBinaryFileHandler<Tuple> *dbWriter = new TensorBinaryFileHandler<Tuple>(splitFileName, O_TRUNC);
			arrDBSeg[mode*Option::szThread + tid] = dbWriter;
		}

		// Initialize HEAP for efficient K-way merge sort.
		TuplefpPair tempPair;
		for (int i = 0; i < sizeSegmentTempFile; i++){
			if(listSortedTempFile[i]->getNext(temp)){
				tempPair.init(temp, listSortedTempFile[i]);
				heap.push_back(tempPair);
			}
		}
		std::make_heap(heap.begin(), heap.end(), 
					[this](const TuplefpPair &a, const TuplefpPair &b){ 
						return TensorIndexWithMultipleFile<Tuple>::cmp(a.key, b.key, this->sortOrder);});

		// Loop all records via HEAP, and distribute records into temporarly files.
		while (heap.size() > 0){
			tempPair = heap.front();
			std::pop_heap(heap.begin(), heap.end(), 
						[this](const TuplefpPair &a, const TuplefpPair &b){ 
							return TensorIndexWithMultipleFile<Tuple>::cmp(a.key, b.key, this->sortOrder);});
			heap.pop_back();

			arrDBSeg[mode*Option::szThread + (tempPair.key.dim[mode] % Option::szThread)]->appendRecord(tempPair.key);

			if (tempPair.hfp->getNext(tempPair.key)){
				heap.push_back(tempPair);
				std::push_heap(heap.begin(), heap.end(), 
						[this](const TuplefpPair &a, const TuplefpPair &b){ 
							return TensorIndexWithMultipleFile<Tuple>::cmp(a.key, b.key, this->sortOrder);});
			}
		}

		// clean up 
		for (int i = 0; i < sizeSegmentTempFile; i++){
			listSortedTempFile[i]->cleanup(TensorBinaryFileHandler<Tuple>::HARD_CLEANUP);
			delete listSortedTempFile[i];
		}
	}

	return;
}

#endif // #ifndef _TENSROFILE_
