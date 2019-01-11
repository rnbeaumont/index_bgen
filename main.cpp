#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <unordered_map>
#include "general_functions.h"
#include "index.h"
#define STRLENPAR 500
using namespace std;
struct Parameters {
	string inputBgenFile;
	string inputSampleFile;
	string inputSNPList;
	string inputIndexFile;
	string inputOutFile;
	string inputChr;
	bool inputBgenFileB;
	bool inputSampleFileB;
	bool inputSNPListB;
	bool inputIndexFileB;
	bool inputOutFileB;
	bool inputGenB;
	bool inputChrB;
	bool inputRsidB;
};
Parameters read_options(int &argc,char *argv[]);
void print_usage();
void print_options(Parameters result);
void process_dosages();
void print_header();

int main( int argc, char *argv[] ){
	using namespace std;
//	print_header();
	Parameters inUse=read_options(argc,argv);
	if(!(inUse.inputBgenFileB && inUse.inputIndexFileB)){
		print_usage();
	}

	// If we're indexing the file
	if(!inUse.inputOutFileB && !inUse.inputGenB && !inUse.inputSNPListB && !inUse.inputSampleFileB){
		// index the file here
		index_snps(inUse.inputBgenFile,inUse.inputIndexFile);
	}else if(inUse.inputOutFileB && inUse.inputSNPListB && inUse.inputChrB && inUse.inputSampleFileB){
		// pull out the SNPs
		extract_snps(inUse.inputBgenFile,inUse.inputSampleFile,inUse.inputIndexFile,inUse.inputSNPList,inUse.inputOutFile,inUse.inputGenB,inUse.inputChr,inUse.inputRsidB);
	}else{
		print_usage();
	}
	puts( "\tDONE!\n\n" );
	return 0;
}

Parameters read_options(int &argc,char *argv[]){
	int m,n,	/* Loop counters */
	    l,		/* String length */
	    x;		/* Exit code. */
        char argcase[10];	/* List buffer */
	Parameters result;
//initialise result
	result.inputIndexFileB=false;
	result.inputSNPListB=false;
	result.inputSampleFileB=false;
	result.inputBgenFileB=false;
	result.inputOutFileB=false;
	result.inputGenB=false;
	result.inputChrB=false;

	for( n=1; n < argc; n++ )	/* Scan through args */
	{
		switch( (int)argv[n][0] )	/* Check for option character */
		{
		case '-': x = 0;	/* Bail out if 1. */
			  l = strlen( argv[n] );
			for( m = 1; m < l; ++m )	/* scan through options */
			{
				if( strcmp(&argv[n][1],"bgen")==0 ){
					if( (int)argv[n+1][0] == '-' ){
						puts( "Illegal syntax -- no string following option -bgen\n" );
						exit( 1 );
					}else{
						result.inputBgenFile=&argv[n+1][0];
						result.inputBgenFileB=true;
					}
					x = 1;
				}else if( strcmp(&argv[n][1],"sample")==0 ){
					if( (int)argv[n+1][0] == '-' ){
						puts( "Illegal syntax -- no string following option -sample\n" );
						exit( 1 );
					}else{
						result.inputSampleFile=&argv[n+1][0];
						result.inputSampleFileB=true;
/*						printf( "String = %s\n", s );*/
					}
					x = 1;
				}else if( strcmp(&argv[n][1],"snps")==0 ){
					if( (int)argv[n+1][0] == '-' ){
						puts( "Illegal syntax -- no string following option -snps\n" );
						exit( 1 );
					}else{
						result.inputSNPList=&argv[n+1][0];
						result.inputSNPListB=true;
					}
				}else if( strcmp(&argv[n][1],"index-file")==0 ){
					if( (int)argv[n+1][0] == '-' ){
						puts( "Illegal syntax -- no string following option -index-file\n" );
						exit( 1 );
					}else{
						result.inputIndexFile=&argv[n+1][0];
						result.inputIndexFileB=true;
					}
				}else if( strcmp(&argv[n][1],"out")==0 ){
					if( (int)argv[n+1][0] == '-' ){
						puts( "Illegal syntax -- no string following option -out\n" );
						exit( 1 );
					}else{
						result.inputOutFile=&argv[n+1][0];
						result.inputOutFileB=true;
					}
				}else if( strcmp(&argv[n][1],"gen")==0 ){
					result.inputGenB=true;
				}else if( strcmp(&argv[n][1],"chr")==0 ){
					if( (int)argv[n+1][0] == '-' ){
						puts( "Illegal syntax -- no string following option -chr\n" );
						exit( 1 );
					}else{
						result.inputChr=&argv[n+1][0];
						result.inputChrB=true;
					}
				}else if( strcmp(&argv[n][1],"rsid")==0 ){
					result.inputRsidB=true;
				}else{
					strcpy( argcase, &argv[n][1] );
					printf( "Unrecognised option : %s\n", &argv[n][0] ); /* THIS NEEDS FIXING TO OUTPUT THE WHOLE ARGUMENT!!!!!!! */
					print_usage() ;
				}
				if( x == 1 ) {
					break;
				}
			}
			break;
		default:  /*printf( "Unrecognised option : %s\n", &argv[n][0] );	Not option -- text. */
			break;
		}
	}
	return result;
}

void print_usage( ){
	printf( "\nUsage:\n" );
	printf( "\t-bgen            bgen file\n" );
	printf( "\t-sample          sample file\n" );
	printf( "\t-index-file      file containing SNP indexes\n" );
	printf( "\t-snps            SNPs to extract\n" );
	printf( "\t-out             output filename\n" );
	printf( "\t-chr             chr from which to extract SNPs\n" );
	printf( "\t-gen             (optional) specifies gen output rather than genotype dosages\n" );
	printf( "\t-rsid            (optional) include rsids in output file (where available) rather than chr:pos\n\n\n" );
	exit(0);
	return;
}

void print_options(Parameters result){
        cout<<"\nOptions in use:"<<endl;
        if(result.inputBgenFileB){
                cout<<"-bgen             "<<result.inputBgenFile<<endl;
        }
        if(result.inputSampleFileB){
                cout<<"-sample           "<<result.inputSampleFile<<endl;
        }
        if(result.inputSNPListB){
                cout<<"-snps             "<<result.inputSNPList<<endl;
        }
        if(result.inputIndexFileB){
                cout<<"-index-file       "<<result.inputIndexFile<<endl;
        }
        if(result.inputOutFileB){
                cout<<"-out              "<<result.inputOutFile<<endl;
        }
        if(result.inputGenB){
                cout<<"-gen              true"<<endl;
        }
        if(result.inputGenB){
                cout<<"-chr              "<<result.inputChr<<endl;
        }
        cout<<endl;
	return;
}

void print_header(){
	cout<<"\n\nindex_bgen\n"<<endl;
	return;
}
