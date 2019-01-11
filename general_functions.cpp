#include "general_functions.h"
#include <zlib.h>

extern string chr;

bgen::bgen(){	//default constructor
	aA="";		//initialise string variables
	aB="";		//initialise string variables
	rsidS="";	//initialise string variables
	currentSnp=0;	//initialise currentSnp to 0 to be able to increment
	readVarHead=false;	//whether the head of thevariant has been read but not the probabilities
}

bgen::bgen(string file){	//constructor to open file
	using namespace std;
	unsigned char buf[4];	//buffer
	unsigned int magic;
	aA="";		//initialise string variables
	aB="";		//initialise string variables
	rsidS="";	//initialise string variables
	currentSnp=0;	//initialise currentSnp to 0 to be able to increment
	readVarHead=false;	//whether the head of thevariant has been read but not the probabilities
	input.open(file,ios::binary);	//open file
	if(!input.is_open()){
		cerr<<"ERROR: file "<<file<<" does not exists"<<endl;
		exit(-20);
	}
	input.read((char*)(&offset),4);	//read offset
	input.read((char*)(&LH),4);	//read LH
	input.read((char*)(&M),4);	//read M
	input.read((char*)(&N),4);	//read N
	input.read((char*)(&magic),4);	//read magic number
	if(!(magic==1852139362 || magic==0)){
		std::cerr<<"ERROR: Malformed bgen file. \'Magic number\' bytes not valid"<<std::endl;
		exit(-8);
	}
	// check whether there's anything in the free data area, and if so skip it
	if((LH-20)>0){
		input.seekg((LH-20),ios_base::cur);
	}
	//read the flags
	//Check which format the data is in
	input.read((char*)(&buf[0]),4);
	compressed=((buf[0] >> 0) & 1);	//CompressedSNPBlocks
	layout=((buf[0] >> 2) & 15);	//layout
	if(!(layout==1 || layout==2)){
		std::cerr<<"ERROR: file is in bgen format "<<layout<<". Only bgen files v1.1 and v1.2 supported"<<std::endl;
		exit(-8);
	}
	sampleIds=((buf[0]>>31) & 1);	//sampleIds		//TODO if this is 1 then I need to read the sample block!!!
	//Header block read so skip to the start of the genotypes
	input.seekg((offset+4),ios_base::beg);
	currentByte=4+offset;
}

int bgen::get_length(){
	int temp=input.tellg();
	input.seekg(0,input.end);
	int length=input.tellg();
	input.seekg(temp,input.beg);
	return length;
}

bool bgen::read_variant_id(){
	using namespace std;
	unsigned int lrsid;
	unsigned char buf[4];
//	unsigned long int uncompressedLength;
	vector<unsigned char> alleleTemp;
	//if v1.1 read number of individuals
	if(layout==1){
		input.read((char*)(&nv),4);
	}
	//read variant name
	//vid
	//set all bytes to zero as we'll be reading less than 4 bytes next!
	buf[0]=*(unsigned char*)"\0";
	buf[1]=*(unsigned char*)"\0";
	buf[2]=*(unsigned char*)"\0";
	buf[3]=*(unsigned char*)"\0";
	input.read((char*)(&buf),2);
	lrsid=(buf[1]<<8) | buf[0];
	vid.clear();
	vid.resize(lrsid);
	input.read((char*)(&vid[0]),lrsid);
//	input.read((char*)(&nv),4);//rsid
	input.read((char*)(&buf),2);
	lrsid=(buf[1]<<8) | buf[0];
	rsid.clear();
	rsid.resize(lrsid);
	input.read((char*)(&rsid[0]),lrsid);
	rsidS.clear();
	rsidS=string(rsid.begin(),rsid.end());
	//chr
	input.read((char*)(&buf),2);
	lrsid=(buf[1]<<8) | buf[0];
	chr.clear();
	chr.resize(lrsid);
	input.read((char*)(&chr[0]),lrsid);
	//position
	input.read((char*)(&pos),4);
	//if it's version 1.2 check how many alleles there are
	if(layout==2){
		input.read((char*)(&buf),2);
		K=(buf[1]<<8) | buf[0];
		alleleString.resize(K);
	}else{
		alleleString.resize(2);
	}
	//allele A
	input.read((char*)(&lrsid),4);
	alleleA.clear();
	alleleA.resize(lrsid);
	input.read((char*)(&alleleA[0]),lrsid);
	aA.clear();
	aA=string(alleleA.begin(),alleleA.end());
	alleleString[0]=string(alleleA.begin(),alleleA.end());
	//allele B
	input.read((char*)(&lrsid),4);
	alleleB.clear();
	alleleB.resize(lrsid);
	input.read((char*)(&alleleB[0]),lrsid);
	aB.clear();
	aB=string(alleleB.begin(),alleleB.end());
	alleleString[1]=string(alleleB.begin(),alleleB.end());
	//if v1.2 read in the other alleles
	if(layout==2){
		for(unsigned int i=2;i<K;++i){
			input.read((char*)(&lrsid),4);
			alleleTemp.clear();
			alleleTemp.resize(lrsid);
			input.read((char*)(&alleleTemp[0]),lrsid);		//read these in properly and store TODO
			alleleString[i]=string(alleleTemp.begin(),alleleTemp.end());
		}
	}
	//if v1.2 read in length of rest of data block
	if(layout==2){
		input.read((char*)(&nbytes),4);
//		if(compressed){
//			input.read((char*)(&uncompressedLength),4);
//			nbytes-=4;
//		}
	}else if(layout==1){	//v1.1
		//length of compressed data block
		if(compressed){
			input.read((char*)(&nbytes),4);
		}else{
			nbytes=6*N;
		}
	}else{
		cerr<<"ERROR: only bgen version 1.1 and 1.2 supported"<<endl;
		exit(-7);
	}
	//check number of individuals is the same as the number of individuals in the header block, else exit
	if(layout==1){
		if(nv!=N){
			cerr<<"ERROR: Snp ";
			for(vector<unsigned char>::const_iterator i=rsid.begin();i!=rsid.end();++i){
				cerr<<*i;
			}
			cerr<<" has fewer individuals than the bgen header block indicates. I don't know which individuals these are so exiting"<<endl;
			exit(10);
		}
	}
	currentSnp++;
	probs.clear();
	readVarHead=true;
	if(!input.eof()){
		return true;
	}else{
		return false;
	}
}

void bgen::read_variant_probabilities(){
	using namespace std;
	unsigned long int uncompressedLength=0;
	unsigned char *uncompressedBuf;
//	unsigned char probBuf[4];
	float probBuf,prob1,prob2;
	unsigned char *compressedData;
	int uncompressedValue;
//	unsigned char buf[4];
	unsigned int nalleles,minp,maxp,samplePloidy,phased,nbit,sampleMissing,mask,shift,fByte,lByte;
	unsigned int intBuf;

	if(!readVarHead){
		read_variant_id();
	}
	probs.clear();
	if(layout==1){
		uncompressedLength=6*N;
		uncompressedBuf=new unsigned char[uncompressedLength];
		if(compressed){
			compressedData=new unsigned char[nbytes];
			input.read((char*)(&compressedData[0]),nbytes);
			uncompressedValue=uncompress(uncompressedBuf,&uncompressedLength,(const Bytef*)compressedData,nbytes);
			delete [] compressedData;
			if(uncompressedValue!=0){
				cout<<"Error reading variant, exiting!"<<endl;
				exit(10);
			}
		}else{
			input.read((char*)(&uncompressedBuf[0]),uncompressedLength);
		}
		for(unsigned int i=0;i<uncompressedLength;i+=2){
			probBuf=(uncompressedBuf[i+1]<<8) | uncompressedBuf[i];
			probs.push_back(probBuf/32768.0);
		}
		delete [] uncompressedBuf;
	}else if(layout==2){
		if(compressed){
			input.read((char*)(&uncompressedLength),4);
			uncompressedBuf=new unsigned char[uncompressedLength];
			nbytes-=4;
			compressedData=new unsigned char[nbytes];
			input.read((char*)(&compressedData[0]),nbytes);
			uncompressedValue=uncompress(uncompressedBuf,&uncompressedLength,(const Bytef*)compressedData,nbytes);
			delete [] compressedData;
			if(uncompressedValue!=0){
				cout<<"Error reading variant, exiting!"<<endl;
				exit(10);
			}
		}else{
			uncompressedLength=nbytes;
			uncompressedBuf=new unsigned char[uncompressedLength];
			input.read((char*)(&uncompressedBuf[0]),uncompressedLength);
		}
		//number of individuals in data block (must equal N)
		nv=(uncompressedBuf[3]<<24) | (uncompressedBuf[2]<<16) | (uncompressedBuf[1]<<8) | uncompressedBuf[0];	//TODO check this equals header (N)
		if(nv != N){
			cerr<<"ERROR: number of individuals in data block differs from that indicated by the header block"<<endl;
			exit(-11);
		}
		//number of alleles (must equal K)
		nalleles=(uncompressedBuf[5]<<8) | uncompressedBuf[4];
		if(nalleles!=K){
			cerr<<"ERROR: number of alleles in data block differs from that indicated by the header block"<<endl;
			exit(-12);
		}
		//min ploidy
		minp=uncompressedBuf[6];
		//max ploidy
		maxp=uncompressedBuf[7];
		//is the data phased
		phased=uncompressedBuf[8+nv];
		//how many bits are used to store each probability
		nbit=uncompressedBuf[9+nv];
		mask=(pow(2,nbit)-1.0);
		if(minp==2 && maxp==2 && phased==0 && nalleles==2){
			for(unsigned int i=0;i<2*nv;i+=2){
				samplePloidy=0;
				sampleMissing=0;
				samplePloidy=((uncompressedBuf[8+i/2] >> 0) & 63);
				sampleMissing=((uncompressedBuf[8+i/2] >> 7) & 1);
				//read 2 (ploidy)
				fByte=10+nv+floor(i*nbit/8.0);  // work out which byte contains the start of my number
				lByte=10+nv+ceil((i+1)*nbit/8.0);       // work out which byte contains the end of my number
				intBuf=0;
				probBuf=0;
				shift=80+nv*8+(i*nbit)-fByte*8; // work out which bit in fByte the mask needs to start
				intBuf=(uncompressedBuf[fByte]>>shift);
				for(unsigned int j=1;j<=lByte-fByte-1;j++){     // loop from first->last
					intBuf=(uncompressedBuf[j+fByte] << (j*8-shift)) | static_cast<int>(intBuf);    //probBuf=probBuf | current byte
				}
				probBuf=intBuf & mask;  // shift relevant bytes to the right and mask
				probs.push_back(probBuf/mask);
				prob1=probBuf/mask;
				//do second probability
				fByte=10+nv+floor((1+i)*nbit/8.0);      // work out which byte contains the start of my number
				lByte=10+nv+ceil((2+i)*nbit/8.0);       // work out which byte contains the end of my number
				intBuf=0;
				probBuf=0;
				shift=80+nv*8+((1+i)*nbit)-fByte*8;     // work out which bit in fByte the mask needs to start
				intBuf=(uncompressedBuf[fByte]>>shift);
				for(unsigned int j=1;j<=lByte-fByte-1;j++){     // loop from first->last
					intBuf=(uncompressedBuf[j+fByte] << (j*8-shift)) | static_cast<int>(intBuf);    //probBuf=probBuf | current byte
				}
				probBuf=intBuf & mask;  // shift relevant bytes to the right and mask
				probs.push_back(probBuf/mask);
				prob2=probBuf/mask;
				if(sampleMissing==0){
					probs.push_back(1-prob1-prob2);
				}else{
					probs.push_back(0);
				}
			}
		}else{			// at the moment error, must fix this TODO
			cerr<<"ERROR: Not all samples are diploid - format not supported"<<endl;
			exit(11);
		}
		delete [] uncompressedBuf;
	}
	readVarHead=false;
}

void bgen::skip_variant_probabilities(){
	using namespace std;
	if(!readVarHead){
		read_variant_id();
	}
	input.seekg(nbytes,ios_base::cur);
	readVarHead=false;
}

void bgen::seek(streampos varStart){
	input.seekg(varStart,input.beg);
}
