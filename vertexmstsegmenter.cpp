#include "vertexmstsegmenter.h"

/***************************************************************************************************************************
 * TriMST
 * Constructors
 **************************************************************************************************************************/

VertexMSTSegmenter::VertexMSTSegmenter()
{
}

VertexMSTSegmenter::VertexMSTSegmenter(string fileM, float kFactor){

    this->meshFile = fileM;
    this->factorK = kFactor;
}


/***************************************************************************************************************************
 * TriMST
 * Set-up functions
 **************************************************************************************************************************/

void VertexMSTSegmenter::loadMesh(){

    mesh = Mesh<Vertex3D, Triangle>();

    QString qmn = QString::fromStdString(meshFile);
    QStringList qsl = qmn.split(".");
    QString suffix = qsl.back();

    if(!suffix.compare("tri"))
        Reader::readMeshFile(mesh, meshFile);
    else if(!suffix.compare("off"))
        Reader::readOFFMesh(mesh, meshFile);
    else{
        cout<<"Format not recognized: "<<suffix.toStdString().c_str()<<endl;
        exit(0);
    }

    mesh.build();

    cout<<"Loaded mesh with "<<mesh.getNumVertex()<<" vertices and "<<mesh.getTopSimplexesNum()<<" top simplexes"<<endl;
}

void VertexMSTSegmenter::loadStructs(){

    clusterIndex = new int[mesh.getNumVertex()];

    calculateNormals();
    setMeshBarycenter();
    setTriangleAreas();
    getBoundingBoxDiagonal();
    vector<double> auxAreas = createRegionAreas();

    for(int ii=0; ii<mesh.getNumVertex(); ii++)
        clusterIndex[ii] = -1;

    for(int ii=0; ii<mesh.getNumVertex(); ii++){

        segRegion SR;
        SR.members.insert(ii);
        SR.regArea = auxAreas.at(ii);
        double areaRatio = AREA_RATIO(SR.regArea, totMArea);
        SR.internalDifference = THR_FORMULA(areaRatio, factorK);
        SR.regionIndex=ii;
        SR.size=1;
        Segmentation.push_back(SR);
    }
}

void VertexMSTSegmenter::calculateNormals(){

    for(int ii=0; ii<mesh.getTopSimplexesNum(); ii++){

        Triangle t = mesh.getTopSimplex(ii);
        Normals aux = Normals(mesh.getVertex(t.TV(0)), mesh.getVertex(t.TV(1)), mesh.getVertex(t.TV(2)));
        norms.push_back(aux);
    }

    cout<<"Computed normals for "<<norms.size()<<" triangles"<<endl;
}

void VertexMSTSegmenter::setMeshBarycenter(){

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

void VertexMSTSegmenter::setTriangleAreas(){

    triAreas = new float[mesh.getTopSimplexesNum()];
    totMArea = 0.0;

    for(int ii=0; ii<mesh.getTopSimplexesNum(); ii++){
        triAreas[ii] = mesh.TArea(ii);
        totMArea += triAreas[ii];
    }
}

void VertexMSTSegmenter::getBoundingBoxDiagonal(){

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

vector<double> VertexMSTSegmenter::createRegionAreas(){

    vector<double> ret;
    for(int ii=0; ii<mesh.getNumVertex(); ii++)
        ret.push_back(0.0);

    for(int ii=0; ii<mesh.getTopSimplexesNum(); ii++){

        Triangle T = mesh.getTopSimplex(ii);
        Vertex3D V0 = mesh.getVertex(T.TV(0));
        Vertex3D V1 = mesh.getVertex(T.TV(1));
        Vertex3D V2 = mesh.getVertex(T.TV(2));
        Vertex3D triBaryc = Vertex3D((V0.getX()+V1.getX()+V2.getX())/3.0, (V0.getY()+V1.getY()+V2.getY())/3.0, (V0.getZ()+V1.getZ()+V2.getZ())/3.0);
        Vertex3D half01((V0.getX()+V1.getX())/2.0, (V0.getY()+V1.getY())/2.0, (V0.getZ()+V1.getZ())/2.0);
        Vertex3D half12((V1.getX()+V2.getX())/2.0, (V1.getY()+V2.getY())/2.0, (V1.getZ()+V2.getZ())/2.0);
        Vertex3D half20((V0.getX()+V2.getX())/2.0, (V0.getY()+V2.getY())/2.0, (V0.getZ()+V2.getZ())/2.0);

        int indV0 = T.TV(0);
        int indV1 = T.TV(1);
        int indV2 = T.TV(2);

        ret[indV0] += (triangleArea(V0, half01, triBaryc) + triangleArea(V0, triBaryc, half20));
        ret[indV1] += (triangleArea(V1, half12, triBaryc) + triangleArea(V1, triBaryc, half01));
        ret[indV2] += (triangleArea(V2, half12, triBaryc) + triangleArea(V2, triBaryc, half20));
    }
    return ret;
}

/*************************************************************************************************************************
 * TriMST
 * Distance related functions
 *************************************************************************************************************************/

std::tr1::unordered_map<edgekey, float> VertexMSTSegmenter::vertexDistance(){

    std::tr1::unordered_map<edgekey, float> FD;

    for(unsigned int ii=0; ii<mesh.getNumVertex(); ii++){

        Vertex3D T = mesh.getVertex(ii);
        vector<int> vneigh = mesh.VV(ii);

        for(int jj=0; jj<vneigh.size(); jj++){

            int f2 = vneigh.at(jj);//T.TT(jj);
            edgekey ek = getEdgeKey(ii, f2);

            //If it is not already in the structure
            if(!FD.count(ek) && f2 >= 0){

                float dist = T.distance(mesh.getVertex(f2));//triangleDistanceFormula(ii, f2);
                assert(dist > 0);
                FD[ek] = dist;
            }
        }
    }
    return FD;
}

std::tr1::unordered_map<edgekey, float> VertexMSTSegmenter::fieldDistance(){

    std::tr1::unordered_map<edgekey, float> FD;

    for(unsigned int ii=0; ii<mesh.getNumVertex(); ii++){

        Vertex3D T = mesh.getVertex(ii);
        vector<int> vneigh = mesh.VV(ii);

        for(int jj=0; jj<vneigh.size(); jj++){

            int f2 = vneigh.at(jj);//T.TT(jj);
            edgekey ek = getEdgeKey(ii, f2);

            //If it is not already in the structure
            if(!FD.count(ek) && f2 >= 0){

                float dist = fabs(vertexFunction[ii]-vertexFunction[f2]);
                assert(dist > 0);
                FD[ek] = dist;
            }
        }
    }
    return FD;
}

void VertexMSTSegmenter::distanceBuilder(){

    for(unsigned int ii=0; ii<mesh.getNumVertex(); ii++){

        Vertex3D T = mesh.getVertex(ii);
        vector<int> vneigh = mesh.VV(ii);

        for(int jj=0; jj<vneigh.size(); jj++){

            int f2 = vneigh.at(jj);//T.TT(jj);
            edgekey ek = getEdgeKey(ii, f2);

            //If it is not already in the structure
            if(edgeDistances.count(ek) || f2 >= mesh.getNumVertex()) //key already inserted, ignore and move on
                continue;
            float dist1 = fabs(vertexFunction[ii]-vertexFunction[f2]);
            float dist2 = T.distance(mesh.getVertex(f2));
            assert(dist1 >= 0 && dist2>=0);
            edgeDistances[ek] = dist2*alpha+dist1/BBDiagonal;

            edgeWeight ew;
            ew.key = ek;
            ew.w = edgeDistances[ek];
            if(ii<f2){
                ew.triangle1=ii;
                ew.triangle2=f2;
            }
            else{
                ew.triangle2=ii;
                ew.triangle1=f2;
            }
            sortedWeights.push(ew);
            ppSorted.push(ew);
        }
    }
}


/**********************************************************************************************************************
 * TriMST
 * Core functions
 **********************************************************************************************************************/

void VertexMSTSegmenter::initSegmentation(){

    //initialize the edge queue
    distanceBuilder();
    cout<<"Distance built"<<endl;
}

void VertexMSTSegmenter::mainIteration(){

    int progressive = 0;
    while (!sortedWeights.empty()) {

        edgeWeight first = sortedWeights.top();
        sortedWeights.pop();
        int find1 = lookForVertex(first.triangle1);
        int find2 = lookForVertex(first.triangle2);
        //cout<<"Find1 "<<find1<<" Find2 "<<find2<<" V1 "<<first.triangle1<<" V2 "<<first.triangle2<<endl;

        double D1 = Segmentation.at(find1).internalDifference;
        double D2 = Segmentation.at(find2).internalDifference;

        if(first.w < D1 && first.w < D2){
            mergeRegions(find1, find2, first.w);
            progressive++;
        }
    }
    cout<<"Post-processing..."<<endl;
    //for(int kk=0;kk<mesh.getNumVertex();kk++)
      //  cout<<"Vert "<<kk<<" "<<clusterIndex[kk]<<endl;
    postProcessMerge();
    cout<<"Done, now reordering..."<<endl;
    reorderIndices();
}

int VertexMSTSegmenter::lookForVertex(int tIndex){

    int search = tIndex;
    int times = 0;

    while(search != Segmentation.at(search).regionIndex){
        search = Segmentation.at(search).regionIndex;
        times++;
    }
    Segmentation.at(tIndex).regionIndex = search;
    return search;
}


/**********************************************************************************************************************
 * TriMST
 * Operations on regions
 **********************************************************************************************************************/

void VertexMSTSegmenter::mergeRegions(int R1, int R2, double difference){

    int winner, other;
    //cout<<"M"<<endl;
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
    Segmentation.at(winner).size += Segmentation.at(other).size;
    Segmentation.at(other).size = 0;
    Segmentation.at(winner).regArea += Segmentation.at(other).regArea;
    updateInternalDifference(R1, R2, difference);
}

void VertexMSTSegmenter::updateInternalDifference(int regInd1, int regInd2, double difference){

    int lower;
    if(regInd1 < regInd2)
        lower = regInd1;
    else
        lower = regInd2;

    double card1 = Segmentation[lower].regArea;
    double card = AREA_RATIO(card1, totMArea);

    double newDiff = difference + THR_FORMULA(card, factorK);

    Segmentation[regInd1].internalDifference = newDiff;
    Segmentation[regInd2].internalDifference = newDiff;
}

void VertexMSTSegmenter::reorderIndices(){

    std::tr1::unordered_map<int, int> placed;
    int progressive=0;

    for(int ii=0; ii<Segmentation.size(); ii++){

        clusterIndex[ii] = lookForVertex(ii);
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
    nRegs = progressive;
}

void VertexMSTSegmenter::postProcessMerge(){

    double thresholdMerge = mesh.MArea()/800.0;
    while (!ppSorted.empty()) {
        edgeWeight ew = ppSorted.top();
        ppSorted.pop();

        int t1 = ew.triangle1;
        int t2 = ew.triangle2;
        int r1 = lookForVertex(t1);
        int r2 = lookForVertex(t2);

        if(r1==r2)
            continue;
        if(Segmentation.at(r1).regArea > thresholdMerge && Segmentation.at(r2).regArea > thresholdMerge)
            continue;

        mergeRegions(r1, r2, ew.w);
    }
}


/**********************************************************************************************************************
 * TriMST
 * Input-Output Operations
 **********************************************************************************************************************/

void VertexMSTSegmenter::readFieldFile(string fieldFile){

    FILE* f = fopen(fieldFile.c_str(), "r");

    vertexFunction = new float[mesh.getNumVertex()];
    int num;
    float fv;
    fscanf(f, "%d", &num);
    cout<<"Num "<<num<<" V "<<mesh.getNumVertex()<<endl;

    for(int ii=0; ii<num; ii++){
        fscanf(f, "%f", &fv);
        vertexFunction[ii] = fv;//functionValue.push_back(fv);
    }

    cout<<"read "<<num<<"elements"<<endl;
}

void VertexMSTSegmenter::writeSegmOnFile(string fileOut){
    FILE* f = fopen(fileOut.c_str(), "w");

    for(int ii=0; ii<mesh.getNumVertex(); ii++)
        fprintf(f, "%d\n", clusterIndex[ii]);

    fclose(f);
}
