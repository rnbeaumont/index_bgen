#include "general_functions.h"
#include "index.h"

void index_snps(string inputBgenFile,string inputIndexFile){
	//tell the user what we're doing
	unsigned int count=0;
	cout<<"\n\tIndexing Bgen file "<<inputBgenFile<<"..."<<endl;
	//open bgen
	bgen bGenFile(inputBgenFile);
	//open index file
	ofstream outFile;
	outFile.open(inputIndexFile);
	//read each block
	streampos varStart=bGenFile.currentByte;
	while(bGenFile.read_variant_id()){
		bGenFile.skip_variant_probabilities();
		string rsid=bGenFile.get_rsid();
		unsigned int pos=bGenFile.get_pos();
		string aA=bGenFile.get_aA();
		string aB=bGenFile.get_aB();
		outFile<<varStart<<"\t"<<rsid<<"\t"<<pos<<"\t"<<aA<<"\t"<<aB<<endl;
		varStart=bGenFile.input.tellg();
                count++;
                if(count==bGenFile.get_m()) break;
	}
}

void extract_snps(string inputBgenFile,string inputSampleFile,string inputIndexFile,string inputSNPList,string inputOutFile,bool inputGenB,string chr,bool inputRsidB){
	using boost::lexical_cast;
	ifstream SNPFile(inputSNPList);
	unordered_map<string, int> hash;
	vector <streampos> SNP_positions;
	vector <string> fields, outputVector;
	string tmp;
	unsigned int count=0, vsize=0;
	float dosage=0;
	//tell the user what we're doing
	cout<<"\n\tExtracting SNPs in file "<<inputSNPList<<" on chromosome "<<chr<<" from Bgen file "<<inputBgenFile<<"..."<<endl;
	//read in SNP list
	for(string str; getline(SNPFile,str); ){
		split(fields,str,is_any_of(" \t"));
		if(fields[0].compare(chr)==0){
			hash[fields[1]+"_"+fields[2]+"_"+fields[3]]=count;
			count++;
		}
	}
	//check where they are (and if the exist)!
	ifstream indexFile(inputIndexFile);
	for(string str; getline(indexFile,str); ){
		split(fields,str,is_any_of(" \t"));
		if(hash.find(fields[2]+"_"+fields[3]+"_"+fields[4])!=hash.end() || hash.find(fields[2]+"_"+fields[4]+"_"+fields[3])!=hash.end()){
			SNP_positions.push_back(lexical_cast<long long int>(fields[0].c_str()));	//add SNP positions to vector so I can loop over it to pull them out
		}
	}
	//open the bgen file
	bgen bGenFile(inputBgenFile);
	//open output file
	ofstream outFile;
	outFile.open(inputOutFile);
	//open the sample file and set up a vector from which to output
	if(!inputGenB){
		ifstream sampleFile(inputSampleFile);
		getline(sampleFile,tmp);
		split(fields,tmp,is_any_of(" \t"));
		outputVector.push_back(fields[0]);
		getline(sampleFile,tmp);
		for(string str; getline(sampleFile,str); ){
			split(fields,str,is_any_of(" \t"));
			outputVector.push_back(fields[0]);
		}
		vsize=outputVector.size();
	}
        if(vsize-1!=bGenFile.get_n()){
            cout<<"ERROR: Number of samples in sample file doesn't match that in the bgen file"<<endl;
            exit (EXIT_FAILURE);
        }
	//seek to each position and read the variant(s), store them in a vector for output. First do dosages, make it work then implement Gen
	for(vector<streampos>::iterator i=SNP_positions.begin();i!=SNP_positions.end();++i){
		bGenFile.seek(*i);	//seek to start of SNP
		bGenFile.read_variant_id();	//read variant id
		bGenFile.read_variant_probabilities();	//read probabilities
		if(!inputGenB){
			if(!inputRsidB){
				outputVector[0].append("\tchr"+chr+":"+to_string(bGenFile.get_pos())+"_"+bGenFile.get_aA()+"_"+bGenFile.get_aB());
			}else{
				if(bGenFile.get_rsid_compare()){
					outputVector[0].append("\t"+bGenFile.get_rsid()+"_"+bGenFile.get_aA()+"_"+bGenFile.get_aB());
				}else{
					outputVector[0].append("\t"+chr+":"+to_string(bGenFile.get_pos())+"_"+bGenFile.get_aA()+"_"+bGenFile.get_aB());
				}
			}
			for(unsigned int i=1;i<vsize;++i){
				dosage=bGenFile.probs[3*i-2]+2*bGenFile.probs[3*i-1];
				outputVector[i].append("\t"+to_string(dosage));
			}
		}else{
			outFile<<chr<<"\t"<<bGenFile.get_rsid()<<"\t"<<bGenFile.get_pos()<<"\t"<<bGenFile.get_aA()<<"\t"<<bGenFile.get_aB();
			for(vector<float>::iterator prob=bGenFile.probs.begin();prob!=bGenFile.probs.end();++prob){
				outFile<<"\t"<<*prob;
			}
			outFile<<endl;
		}
	}
	//output in the requested format
	if(!inputGenB){
		for(vector<string>::iterator i=outputVector.begin();i!=outputVector.end();++i){
			outFile<<*i<<endl;
		}
	}
}
