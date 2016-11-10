#ifndef VERTEXMSTSEGMENTER_H
#define VERTEXMSTSEGMENTER_H

#include "common.h"

namespace std { using namespace __gnu_cxx; }

typedef unsigned long int vertexind;
typedef std::tr1::unordered_set<int> segMembers;

class VertexMSTSegmenter
{
public:
    VertexMSTSegmenter();
    VertexMSTSegmenter(string fileM, float kFactor);

    inline void callLoad(){
        this->loadMesh();
        this->loadStructs();
    }

    inline void setAlpha(float a){
        this->alpha = a;
    }
    inline void callReadF(string f){
        this->readFieldFile(f);
    }
    inline void callWriteF(string f){
        this->writeSegmOnFile(f);
    }

    inline void callInit(){
        this->initSegmentation();
    }
    inline void triggerSegmentation(){
        this->mainIteration();
    }


private:
    Mesh<Vertex3D, Triangle> mesh;
    string meshFile;
    float factorK;
    float BBDiagonal;
    int *clusterIndex;
    float *vertexFunction;
    float *triAreas;
    Vertex3D meshBarycenter;
    int nRegs;
    float totMArea;
    float alpha;
    vector<Normals> norms;

    //Distance related
    priority_queue<edgeWeight, vector<edgeWeight>, compare> sortedWeights;
    priority_queue<edgeWeight, vector<edgeWeight>, compare> ppSorted;
    std::tr1::unordered_map<edgekey, double> edgeDistances;

    //Regions
    vector<segRegion> Segmentation;

    /**
     * Set-up functions
     */
    void loadMesh();
    void loadStructs();
    void calculateNormals();
    void setTriangleAreas();
    void setMeshBarycenter();
    void getBoundingBoxDiagonal();

    /**
     * Distance calculators
     */
    //double triangleDistanceFormula(int f1, int f2);
    void distanceBuilder();
    std::tr1::unordered_map<edgekey, float> vertexDistance();
    std::tr1::unordered_map<edgekey, float> fieldDistance();
    double mIntCalculator(int regInd1, int regInd2);
    bool checkMergePredicate(edgeWeight EW);
    Vertex3D halfPoint(Vertex3D v1, Vertex3D v2);

    /**
     * Core functions
     */
    void initSegmentation();
    void mainIteration();
    int lookForVertex(int tIndex);
    int lookForVert(int tIndex);

    inline edgekey getEdgeKey(faceind e1, faceind e2){

        if(e1<=e2)
            return (edgekey(e1) << 32 | edgekey(e2));
        else
            return (edgekey(e2) << 32 | edgekey(e1));
    }

    /**
     * Operations on regions
     */
    void mergeRegions(int R1, int R2, double difference);
    void updateInternalDifference(int regInd1, int regInd2, double difference);
    void reorderIndices();
    void postProcessMerge();
    vector<double> createRegionAreas();

    inline double triangleArea(Vertex3D v0, Vertex3D v1, Vertex3D v2){

        double edge0 = v0.distance(v1);
        double edge1 = v1.distance(v2);
        double edge2 = v2.distance(v0);
        double semiperimeter = (edge0 + edge1 + edge2)/2.0;

        double area = sqrt(semiperimeter*(semiperimeter-edge0)*(semiperimeter-edge1)*(semiperimeter-edge2));

        return area;
    }

    /**
     * Input - Output
     */
    void readFieldFile(string fieldFile);
    void writeSegmOnFile(string fileOut);
};

#endif // VERTEXMSTSEGMENTER_H
