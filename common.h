#ifndef COMMON_H
#define COMMON_H

#include "LibTri/normals.h"
#include "LibTri/Timer.h"
#include "LibTri/Reader.h"
#include <tr1/unordered_map>
#include <queue>
#include <tr1/unordered_set>
#include <numeric>
#include <QString>
#include <QStringList>
#include <cfloat>
#include <fstream>

#define THR_FORMULA(cardinality, FK) (FK/cardinality)
#define AREA_RATIO(Asmall, Abig) (Asmall/Abig)

/// types for indices of triangles and for adjacencies between them
typedef unsigned long int faceind;
typedef unsigned long long int edgekey;

typedef std::tr1::unordered_set<int> segMembers;

/// types for grid initialization
//typedef signed char coordind;
//typedef int gridkey;


struct edgeWeight{

    edgekey key;
    double w;
    int triangle1;
    int triangle2;
};

struct compare{

    bool operator()(edgeWeight ew1, edgeWeight ew2){
        return ew1.w > ew2.w;
    }
};

//We need a structure for the regions (members, internal difference)
struct segRegion{

    segMembers members;
    double regArea;
    double internalDifference;
    int regionIndex;
    int size;
};

#endif // COMMON_H
