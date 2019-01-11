// index.h
#ifndef INDEX_H
#define INDEX_H
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
#include <boost/lexical_cast.hpp>
//#include <boost/iostreams/filtering_streambuf.hpp>
//#include <boost/iostreams/filtering_stream.hpp>
//#include <boost/iostreams/copy.hpp>
//#include <boost/iostreams/filter/gzip.hpp>
//#include <boost/algorithm/string.hpp>
//#include "libs/iostreams/src/zlib.cpp"
using namespace std;
using namespace boost::iostreams;
using namespace boost;

void index_snps(string inputBgenFile,string inputIndexFile);
void extract_snps(string inputBgenFile,string inputSampleFile,string inputIndexFile,string inputSNPList,string inputOutFile,bool inputGenB,string chr,bool inputRsidB);

#endif
