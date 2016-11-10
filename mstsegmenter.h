#ifndef MSTSEGMENTER_H
#define MSTSEGMENTER_H

//#include "LibTri/normals.h"
//#include "LibTri/Timer.h"
//#include "LibTri/Reader.h"

//#include <QString>
//#include <QStringList>
//#include <tr1/unordered_map>
//#include <queue>
//#include <tr1/unordered_set>
//#include <stdio.h>
#include "common.h"

//#define THR_FORMULA(cardinality, FK) (FK/cardinality)
//#define AREA_RATIO(Asmall, Abig) (Asmall/Abig)

namespace std { using namespace __gnu_cxx; }

typedef unsigned long long int edgekey;
typedef unsigned long int faceind;

//struct edgeWeight{

//    edgekey key;
//    double w;
//    int triangle1;
//    int triangle2;
//};

//struct compare{

//    bool operator()(edgeWeight ew1, edgeWeight ew2){
//        return ew1.w > ew2.w;
//    }
//};


////We need a structure for the regions (members, internal difference)
//struct segRegion{

//    segMembers members;
//    double regArea;
//    double internalDifference;
//    int regionIndex;
//    int size;
//};


class MSTsegmenter
{
public:
    MSTsegmenter();
    MSTsegmenter(string fileM, float FK);
    Mesh<Vertex3D, Triangle> mesh;

    vector<Vertex3D> facesBarycenters;
    vector<Normals> norms;

    inline void callLoad(){
        this->loadMesh();
        this->loadStructs();
    }

    inline edgekey getEdgeKey(faceind e1, faceind e2){

        if(e1<=e2)
            return (edgekey(e1) << 32 | edgekey(e2));
        else
            return (edgekey(e2) << 32 | edgekey(e1));
    }

    inline void callInit(){
        this->initSegmentation();
    }

    inline void triggerSegmentation(){
        this->mainIteration(10);
    }

    inline void setAlpha(double A){
        this->alpha = A;
    }

    void writeSegmentation(string fileOut);
    inline void writeNumbers(string fileN){

        FILE* f = fopen(fileN.c_str(), "a");

        fprintf(f, "%d %f --> %d\n", mesh.getTopSimplexesNum(), factorK, nRegs);

        fclose(f);
    }

private:
    void loadMesh();

    string filename;
    float factorK; //related to number of merges performed
    float BBDiagonal;
    int *clusterIndex;
    float *triAreas;
    Vertex3D meshBarycenter;
    int nRegs;
    float totMArea;
    float alpha;
    //internal differences
    //float *thresholds;

    //Distance related
    priority_queue<edgeWeight, vector<edgeWeight>, compare> sortedWeights;
    priority_queue<edgeWeight, vector<edgeWeight>, compare> ppSorted;
    std::tr1::unordered_map<edgekey, double> edgeDistances;

    //Regions
    vector<segRegion> Segmentation;

    /**
     * Set-up functions
     */
    void loadStructs();
    void calculateNormals();
    void retrieveTriBarycenters();
    void setTriangleAreas();
    void setMeshBarycenter();
    void getBoundingBoxDiagonal();

    /**
     * Distance calculators
     */
    double triangleDistanceFormula(int f1, int f2);
    void distanceBuilder();
    double mIntCalculator(int regInd1, int regInd2);
    bool checkMergePredicate(edgeWeight EW);
    Vertex3D halfPoint(Vertex3D v1, Vertex3D v2);

    inline std::tr1::unordered_map<edgekey, float> faceDistance(){
        std::tr1::unordered_map<edgekey, float> FD;

        for(unsigned int ii=0; ii<mesh.getTopSimplexesNum(); ii++){

            Triangle T = mesh.getTopSimplex(ii);

            for(int jj=0; jj<3; jj++){

                int f2 = T.TT(jj);
                edgekey ek = getEdgeKey(ii, f2);

                //If it is not already in the structure
                if(!FD.count(ek) && f2 >= 0){

                    float dist = triangleDistanceFormula(ii, f2);
                    assert(dist > 0);
                    FD[ek] = dist;
                }
            }
        }
        return FD;
    }

    inline std::tr1::unordered_map<edgekey, float> angleDistance(){
        std::tr1::unordered_map<edgekey, float> aDistances;

        for(int ii=0;ii<mesh.getTopSimplexesNum();ii++){
            Triangle T=mesh.getTopSimplex(ii);

            for(int jj=0;jj<3;jj++){
                int f2=T.TT(jj);
                edgekey ek=getEdgeKey(ii,f2);
                if(aDistances.count(ek)==0 && f2>=0){
                    float diffCenters[3];

                    //Vertices of the shared edge
                    faceind vi1=T.TV((jj+1)%3);
                    faceind vi2=T.TV((jj+2)%3);

                    Vertex3D v=mesh.getVertex(vi1);
                    Vertex3D w=mesh.getVertex(vi2);

                    double diffx=v.getX()-w.getX();
                    double diffy=v.getY()-w.getY();
                    double diffz=v.getZ()-w.getZ();

                    //Length of shared edge
                    double balanceF=sqrt(diffx*diffx + diffy*diffy + diffz*diffz);

                    float normF[3];
                    Normals N=norms.at(ii);
                    float eta;

                    //Normal vector for face ii
                    normF[0]=N.getNx();
                    normF[1]=N.getNy();
                    normF[2]=N.getNz();

                    Vertex3D C1=facesBarycenters.at(ii);
                    Vertex3D C2=facesBarycenters.at(f2);

                    diffCenters[0]=(C2.getX()-C1.getX());
                    diffCenters[1]=(C2.getY()-C1.getY());
                    diffCenters[2]=(C2.getZ()-C1.getZ());

                    Normals nn;
                    nn.setNX(diffCenters[0]);
                    nn.setNY(diffCenters[1]);
                    nn.setNZ(diffCenters[2]);

                    if(N.dotProd(nn) > 0){ //concave
                        eta=1.0;
                    }
                    else{ //convex
                        eta=0.2;
                    }

                    float dot=N.dotProd(norms.at(f2));
                    //Clamp dot product to be in [-1, 1]
                    if(dot > 1.0)
                        dot=1.0;
                    if(dot < -1.0)
                        dot=-1.0;

                    //Divide it by pi to have it normalized in [0,1]
                    float arcC=acos(dot)/M_PI;

                    float d_a=eta*(arcC)*balanceF;

                    aDistances[ek]=d_a;
                }
            }
        }
        return aDistances;
    }

    /**
     * Core functions
     */
    void initSegmentation();
    void mainIteration(int max);
    int lookForTriangle(int tIndex);
    int lookForTri(int tIndex);

    /**
     * Operations on regions
     */
    void mergeRegions(int R1, int R2, double difference);
    void updateInternalDifference(int regInd1, int regInd2, double difference);
    void reorderIndices();
    void postProcessMerge();
    double regionArea(int R);
};

#endif // MSTSEGMENTER_H
