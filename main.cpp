#include <iostream>
#include <stdio.h>
#include "mstsegmenter.h"

using namespace std;

string fileM;
float K = 0.08; //default value
float alpha = 200;
string fileOut;
string fileN;

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
        else if(!strcmp(argv[ii], "-a")){
            double input = atof(argv[ii+1]);
            alpha = input;
            cout<<"alpha "<<alpha<<endl;
        }
        else if(!strcmp(argv[ii], "-out")){
            QString nf = argv[ii+1];
            fileOut = nf.toStdString();
        }
        else if(!strcmp(argv[ii], "-ON")){
            QString nf = argv[ii+1];
            fileN = nf.toStdString();
            cout<<"ON "<<fileN.c_str()<<endl;
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
    MSTS->setAlpha(alpha);
    MSTS->callLoad();
    MSTS->callInit();
    MSTS->triggerSegmentation();

    MSTS->writeSegmentation(fileOut);
    if(strcmp(fileN.c_str(), "")){
        cout<<"Writing on "<<fileN.c_str()<<endl;
        MSTS->writeNumbers(fileN);
    }
    return 0;

}

