#include "mstsegmenter.h"


/***************************************************************************************************************************
 * TriMST
 * Constructors
 **************************************************************************************************************************/

MSTsegmenter::MSTsegmenter(){

}

MSTsegmenter::MSTsegmenter(string fileM, float FK)
{
    this->filename=fileM;
    this->factorK = FK;
}

/***************************************************************************************************************************
 * TriMST
 * Set-up functions
 **************************************************************************************************************************/

void MSTsegmenter::loadMesh(){

    mesh = Mesh<Vertex3D, Triangle>();

    QString qmn = QString::fromStdString(filename);
    QStringList qsl = qmn.split(".");

    QString suffix = qsl.back();
    if(!suffix.compare("tri"))
        Reader::readMeshFile(mesh, filename);
    else if(!suffix.compare("off"))
        Reader::readOFFMesh(mesh, filename);
    else{
        cout<<"Format not recognized: "<<suffix.toStdString().c_str()<<endl;
        exit(0);
    }

    mesh.build();

    cout<<"Loaded mesh with "<<mesh.getNumVertex()<<" vertices and "<<mesh.getTopSimplexesNum()<<" top simplexes"<<endl;
}

void MSTsegmenter::loadStructs(){

    clusterIndex = new int[mesh.getTopSimplexesNum()];

    for(int ii=0; ii<mesh.getTopSimplexesNum(); ii++){
        clusterIndex[ii] = -1;
        //thresholds[ii] = THR_FORMULA(1, factorK);
    }

    //Segmentation.reserve(mesh.getTopSimplexesNum());

    for(unsigned int ii=0; ii<mesh.getTopSimplexesNum(); ii++){
        segRegion SR;
        SR.internalDifference = THR_FORMULA(1, factorK);
        SR.members.insert(ii);
        SR.regionIndex = ii;
        SR.size = 1;
        Segmentation.push_back(SR);
    }

    calculateNormals();
    retrieveTriBarycenters();
    setMeshBarycenter();
    setTriangleAreas();
    getBoundingBoxDiagonal();
    cout<<"All set "<<endl;
}

void MSTsegmenter::calculateNormals(){

    for(int ii=0; ii<mesh.getTopSimplexesNum(); ii++){

        Triangle t = mesh.getTopSimplex(ii);
        Normals aux = Normals(mesh.getVertex(t.TV(0)), mesh.getVertex(t.TV(1)), mesh.getVertex(t.TV(2)));
        norms.push_back(aux);
    }

    cout<<"Computed normals for "<<norms.size()<<" triangles"<<endl;
}

void MSTsegmenter::retrieveTriBarycenters(){


    for(int ii=0; ii<mesh.getTopSimplexesNum(); ii++){
        Vertex3D aux;
        Triangle T = mesh.getTopSimplex(ii);
        int indices[3];
        for(int a=0;a<3;a++)
            indices[a]=T.TV(a);
        Vertex3D v3d[3];
        for(int a=0;a<3;a++)
            v3d[a]=mesh.getVertex(indices[a]);

        float valX=0.0, valY=0.0, valZ=0.0;
        for(int a=0;a<3;a++){
            valX+=v3d[a].getX();
            valY+=v3d[a].getY();
            valZ+=v3d[a].getZ();
        }
        aux.setX(valX/3.0);
        aux.setY(valY/3.0);
        aux.setZ(valZ/3.0);

        facesBarycenters.push_back(aux);
    }

    cout<<"Computed "<<facesBarycenters.size()<<" barycenters"<<endl;
}

void MSTsegmenter::setMeshBarycenter(){

    float sumX=0.0; float sumY=0.0; float sumZ=0.0;

    for(int ii=0;ii<mesh.getNumVertex();ii++){

        sumX += mesh.getVertex(ii).getX();
        sumY += mesh.getVertex(ii).getY();
        sumZ += mesh.getVertex(ii).getZ();
    }
    meshBarycenter.setX(sumX/(float)mesh.getNumVertex());
    meshBarycenter.setY(sumY/(float)mesh.getNumVertex());
    meshBarycenter.setZ(sumZ/(float)mesh.getNumVertex());

    cout<<"Set mesh barycenter"<<endl;

}

void MSTsegmenter::setTriangleAreas(){

    triAreas = new float[mesh.getTopSimplexesNum()];

    for(int ii=0; ii<mesh.getTopSimplexesNum(); ii++)
        triAreas[ii] = mesh.TArea(ii);
}

void MSTsegmenter::getBoundingBoxDiagonal(){

    Vertex3D minC, maxC;

    minC.setX(mesh.getVertex(0).getX());
    minC.setY(mesh.getVertex(0).getY());
    minC.setZ(mesh.getVertex(0).getZ());

    maxC.setX(mesh.getVertex(0).getX());
    maxC.setY(mesh.getVertex(0).getY());
    maxC.setZ(mesh.getVertex(0).getZ());

    for(int ii=1; ii<mesh.getNumVertex(); ii++){

        Vertex3D vv=mesh.getVertex(ii);

        //Update if necessary
        if(vv.getX()<minC.getX())
            minC.setX(vv.getX());
        if(vv.getY()<minC.getY())
            minC.setY(vv.getY());
        if(vv.getZ()<minC.getZ())
            minC.setZ(vv.getZ());

        if(vv.getX()>maxC.getX())
            maxC.setX(vv.getX());
        if(vv.getY()>maxC.getY())
            maxC.setY(vv.getY());
        if(vv.getZ()>maxC.getZ())
            maxC.setZ(vv.getZ());
    }

    //diagonal length
    BBDiagonal = maxC.distance(minC);
}


/*************************************************************************************************************************
 * TriMST
 * Distance related functions
 *************************************************************************************************************************/

double MSTsegmenter::triangleDistanceFormula(int f1, int f2){

    //first triangle
    Triangle t1 = mesh.getTopSimplex(f1);

    int indexT, indV1;
    Vertex3D v1, v2, halfP;
    indexT = -1;

    for(int ii=0;ii<3;ii++){

        if(t1.TT(ii)==f2){
            indexT=f2;
            break;
        }
    }
    //There's something wrong, f1 and f2 are not adjacent
    if(indexT < 0){
        cout<<"Not valid faces"<<endl;
        exit(-1);
    }

    // Common edge and vertices
    Edge* shared = t1.TE((indexT)%3);
    indV1 = shared->EV(0);
    v1 = mesh.getVertex(indV1);
    indV1 = shared->EV(1);
    v2 = mesh.getVertex(indV1);

    halfP = halfPoint(v1, v2);

    //Corresponding centroids
    Vertex3D C1 = facesBarycenters.at(f1);
    Vertex3D C2 = facesBarycenters.at(f2);

    return C1.distance(halfP) + halfP.distance(C2);
}

Vertex3D MSTsegmenter::halfPoint(Vertex3D v1, Vertex3D v2){

    double coordX = (v1.getX() + v2.getX())/2;
    double coordY = (v1.getY() + v2.getY())/2;
    double coordZ = (v1.getZ() + v2.getZ())/2;

    Vertex3D v =  Vertex3D(coordX, coordY, coordZ);
    return v;
}

void MSTsegmenter::distanceBuilder(){

    cout<<"Start distances"<<endl;

    std::tr1::unordered_map<edgekey, float> FD = faceDistance();
    std::tr1::unordered_map<edgekey, float> AD = angleDistance();
    for(int ii=0; ii<mesh.getTopSimplexesNum(); ii++){

        Triangle T = mesh.getTopSimplex(ii);
        for(int jj=0;jj<3;jj++){
            int indexAdj = T.TT(jj);
            if(indexAdj>mesh.getTopSimplexesNum() || indexAdj<0)
                continue;
            edgekey ek = getEdgeKey((faceind)ii, (faceind)indexAdj);
            if(edgeDistances.count(ek)) //key already inserted, ignore and move on
                continue;
            double TDist = (FD[ek] + 150*AD[ek])/BBDiagonal;//triangleDistanceFormula(ii, indexAdj);
            edgeDistances[ek] = TDist;
            edgeWeight ew;
            ew.key = ek;
            ew.w = TDist;
            if(ii<indexAdj){
                ew.triangle1=ii;
                ew.triangle2=indexAdj;
            }
            else{
                ew.triangle2=ii;
                ew.triangle1=indexAdj;
            }
            sortedWeights.push(ew);
            ppSorted.push(ew);
        }
    }
}

double MSTsegmenter::mIntCalculator(int regInd1, int regInd2){

    double D1, D2;
    int region1 = lookForTri(regInd1);
    int region2 = lookForTri(regInd2);

    if(region1<0)
        D1 = factorK;
    else
        D1 = Segmentation.at(region1).internalDifference + THR_FORMULA(Segmentation.at(region1).members.size(), factorK);

    if(region2<0)
        D2 = factorK;
    else
        D2 = Segmentation.at(region2).internalDifference + THR_FORMULA(Segmentation.at(region2).members.size(), factorK);

    if(D1 < D2)
        return D1;
    return D2;
}


bool MSTsegmenter::checkMergePredicate(edgeWeight EW){

    int index1 = EW.triangle1;
    int index2 = EW.triangle2;

    if(EW.w < mIntCalculator(index1, index2))
        return true;
    return false;
}

/**********************************************************************************************************************
 * TriMST
 * Core functions
 **********************************************************************************************************************/

void MSTsegmenter::initSegmentation(){

    //initialize the edge queue
    distanceBuilder();
    cout<<"Distance built"<<endl;
}

void MSTsegmenter::mainIteration(int max){

    int ni = 0;
    int progressive = 0;//mesh.getTopSimplexesNum();
    while(!sortedWeights.empty()){

        edgeWeight first = sortedWeights.top();
        sortedWeights.pop();

        //Pop elem from queue
        //Regions R1 and R2
        int find1 = lookForTri(first.triangle1);
        int find2 = lookForTri(first.triangle2);

        double D1 = Segmentation.at(find1).internalDifference;
        double D2 = Segmentation.at(find2).internalDifference;
        if(first.w < D1 && first.w < D2){
            mergeRegions(find1, find2, first.w);
            progressive++;
        }

    }
    cout<<"Merging"<<endl;
    postProcessMerge();
//    for(int ii=0; ii<Segmentation.size(); ii++){
//        cout<<"Size "<<Segmentation.at(ii).size<<endl;
//    }
    cout<<"Merged, now reordering"<<endl;
    reorderIndices();
}

/**
 * @brief MSTsegmenter::lookForTri
 * @param tIndex index of the triangle
 * @return The index of the region to which tIndex belongs, -1 if it belongs to no region
 */
int MSTsegmenter::lookForTriangle(int tIndex){

    for(unsigned int ii=0; ii<Segmentation.size(); ii++){
        segMembers sm = Segmentation.at(ii).members;
        if(sm.count(tIndex)){
            cout<<"found "<<ii<<endl;
            return ii;
        }
    }
    return -1;
}

int MSTsegmenter::lookForTri(int tIndex){

    int search = tIndex;//Segmentation.at(tIndex).regionIndex;
    //int rIndex = tIndex;
    int times = 0;

    //cout<<search<<" "<<Segmentation.at(search).regionIndex<<endl;
    while(search != Segmentation.at(search).regionIndex){
        search = Segmentation.at(search).regionIndex;
        times++;
    }
    Segmentation.at(tIndex).regionIndex = search;
    //cout<<"T "<<times<<endl;
    return search;
}

/**********************************************************************************************************************
 * TriMST
 * Operations on regions
 **********************************************************************************************************************/

void MSTsegmenter::mergeRegions(int R1, int R2, double difference){

    int winner, other;
    if(R1 < R2){
        winner=R1;
        other=R2;
    }
    else{
        winner=R2;
        other=R1;
    }
    segMembers::iterator it;
    segMembers toReassign = Segmentation.at(other).members;
    for(it=toReassign.begin(); it != toReassign.end(); ++it){
        Segmentation.at(winner).members.insert(*it);
    }
    Segmentation.at(other).regionIndex = Segmentation.at(winner).regionIndex;
    updateInternalDifference(R1, R2, difference);
    Segmentation.at(winner).size += Segmentation.at(other).size;
    Segmentation.at(other).size = 0;
}


void MSTsegmenter::updateInternalDifference(int regInd1, int regInd2, double difference){

    int lower;
    if(regInd1 < regInd2)
        lower = regInd1;
    else
        lower = regInd2;

    int card = Segmentation[lower].members.size();
    //cout<<"C "<<card<<endl;

    double newDiff = difference + THR_FORMULA(card, factorK);
    //cout<<"Old "<<difference<<" new "<<newDiff<<endl;

    Segmentation[regInd1].internalDifference = newDiff;
    Segmentation[regInd2].internalDifference = newDiff;
}

void MSTsegmenter::reorderIndices(){

    std::tr1::unordered_map<int, int> placed;
    int progressive=0;

    for(int ii=0; ii<Segmentation.size(); ii++){

        clusterIndex[ii] = lookForTri(ii);
    }

    //other step to give indices in [0, n-1]
    for(int ii=0; ii<Segmentation.size(); ii++){

        int index = clusterIndex[ii];
        if(placed.count(index))
            clusterIndex[ii] = placed[index];
        else{
            placed[index] = progressive;
            clusterIndex[ii] = progressive++;
        }
    }
    cout<<"Counted "<<progressive<<endl;
}

void MSTsegmenter::postProcessMerge(){

    while (!ppSorted.empty()) {
        edgeWeight ew = ppSorted.top();
        ppSorted.pop();

        int t1 = ew.triangle1;
        int t2 = ew.triangle2;
        int r1 = lookForTri(t1);
        int r2 = lookForTri(t2);

        if(r1==r2)
            continue;
        //cout<<"S1 "<<Segmentation.at(r1).size<<" S2 "<<Segmentation.at(r2).size<<endl;
        int thresholdMerge = mesh.getTopSimplexesNum()/1000.0;
        if(Segmentation.at(r1).size > thresholdMerge && Segmentation.at(r2).size > thresholdMerge)
            continue;

        mergeRegions(r1, r2, ew.w);
    }
}

/**********************************************************************************************************************
 * TriMST
 * Input - Output
 **********************************************************************************************************************/

void MSTsegmenter::writeSegmentation(string fileOut){

    FILE* f = fopen(fileOut.c_str(), "w");

    for(int ii=0; ii<mesh.getTopSimplexesNum(); ii++)
        fprintf(f, "%d\n", clusterIndex[ii]);

    fclose(f);
}
