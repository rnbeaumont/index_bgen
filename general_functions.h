// general_functions.h
#ifndef GENERAL_FUNCTIONS_H
#define GENERAL_FUNCTIONS_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <array>
#include <unordered_map>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/algorithm/string.hpp>
//#include "libs/iostreams/src/zlib.cpp"
using namespace std;
using namespace boost::iostreams;
using namespace boost;

struct bgen{
public:
	bgen();			//default constructor
	bgen(string file);	//constructor to open file when it's all written
	bool read_variant_id();	//read the next variant ID
	void read_variant_probabilities();	//read the probabilities from the previously read variant
	void skip_variant_probabilities();	//skip the probabilities and go to the next variant data block
	unsigned int get_m() const {return M;}		// get M
	unsigned int get_n() const {return N;}		// get N
	string get_rsid() const {return rsidS;}		//get rsidS
	bool get_rsid_compare() const {return rsidS.compare(".");}		//get rsidS
	unsigned int get_pos() const {return pos;}	//get pos
//	string get_chr() const {return chr;}	//get pos
	void seek(streampos varStart);
	unsigned int cSnp() const {return currentSnp;}	//get currentSnp
	string get_aA() const {return aA;}		//get aA
	string get_aB() const {return aB;}		//get aB
	vector<float> get_probs() {return probs;}	// get probs
	vector<float> probs;	// probabilities for last read variant
	int currentByte;
	ifstream input;	// ifstream containing opened bgen file
	int get_length();
	std::vector<std::string> alleleString;      //string containing allele labels (for v1.2 files with more than 2 alleles)

private:
	unsigned int offset;	//offset from start of file where genotype blocks start
	unsigned int LH;	//length in bytes of the header
	unsigned int M;		//number of SNPs
	unsigned int N;		//number of individuals
	bool compressed;	//whether the genotype data is zlib compressed
	int layout;		//layout (bgen version v1.x)
	bool sampleIds;		//sampleIDs - whether sample IDs are stored in the file
	unsigned int nv;	//number of individuals in current data block (must equal N but for some reason stored in each block too)
	unsigned int K;		//number of alleles (if layout==2, else undefined)
	unsigned int currentSnp;		//counter to keep track of which SNP has just been read
	bool readVarHead;	//bool to say whether we're part way through a variant or not
	vector<unsigned char> vid,rsid;	//variant identifier and rsid of current variant
	vector<unsigned char> chr;	//chromosome identifier (not always chromosome, especially if encoded using qctool!!!!!!!!!!!!!)
	unsigned int pos;	//position of current variant
	unsigned int nbytes;	// length of compressed data block (if compressed, else uninitiated!!)
	vector<unsigned char> alleleB;	//alleles
	vector<unsigned char> alleleA;
	string rsidS;	//string version of rsid
	string aA,aB;	// string versions of the alleles
};

#endif
