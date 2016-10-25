#include <iostream>
#include <stdio.h>
#include "mstsegmenter.h"

using namespace std;

string fileM;
float K = 0.08; //default value
string fileOut;

//Command line instructions
void print_help(){

    printf("Print help: ");
    exit(0);
}

void parseArgs(int argc, char* argv[]){

    fileM="";
    cout<<"Parsing"<<endl;
    for(int ii=1; ii<argc; ii+=2){

        if(!strcmp(argv[ii], "-m")){
            QString input = argv[ii+1];
            fileM = input.toStdString();
        }
        else if(!strcmp(argv[ii], "-K")){
            double input = atof(argv[ii+1]);
            K = input;
            cout<<"K "<<K<<endl;
        }
        else if(!strcmp(argv[ii], "-out")){
            QString nf = argv[ii+1];
            fileOut = nf.toStdString();
        }
        else if(!strcmp(argv[ii], "-help")){
            print_help();
        }

        else
            cout<<"Invalid input, ignore ("<<argv[ii]<<")"<<endl;
    }
}

int main(int argc, char* argv[])
{
    cout << "Hello TriMST segmenter!" << endl;

    if(argc<3){
        cout<<"Not enough input arguments, type ./TriMST -help"<<endl;
        exit(-1);
    }

    parseArgs(argc, argv);
    cout<<"K "<<K<<endl;

    MSTsegmenter *MSTS = new MSTsegmenter(fileM, K);
    MSTS->callLoad();
    MSTS->callInit();
    MSTS->triggerSegmentation();

    MSTS->writeSegmentation(fileOut);
    return 0;

}

