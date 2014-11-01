#ifndef DFEATURE_VECTOR_H
#define DFEATURE_VECTOR_H

#include <OpenMesh/Core/IO/MeshIO.hh>

#include <vector>
#include <Eigen/Eigen>

#include <fstream>
#include "Align.h"

class FeatureVector;


class RefMesh {
public :
    int vvs;
    DTriMesh *mesh;
    std::vector<AffineAlign*> align;
    std::vector<int> d;
    std::vector< std::vector< std::pair<int,int> > > rd;
    std::vector<double> w,c;
    std::vector< std::pair<int, std::pair<int,int> > > bfsq;
    bool usec;
    int fvNum;
    ~RefMesh();
    RefMesh(DTriMesh &ms);
    void getFeature(DTriMesh &ms, FeatureVector &fv);
    void printCij(std::ofstream &fout);
    double getWij(int e);
    static double normalScale;
    static int root;
};

class FeatureVector {
public :
    std::vector<Eigen::Matrix3d> s,r,dr,logdr,logr;
    std::vector<bool> isConst;
    std::vector<Eigen::Vector3d> constPoint;
    FeatureVector() {}
    FeatureVector(std::vector<double> weight, std::vector<FeatureVector> &fvs);
    void setConstPoint(int i, Eigen::Vector3d v);
    void loadConstPoint(std::istream& cin);
    void loadConstPoint(std::istream& cin, DTriMesh& mesh);
    void blendFrom(std::vector<double> weight, std::vector<FeatureVector> &fvs);
    void calcLogRFromR();
    void calcRFromLogR();
    void resize(const FeatureVector &other);
};

std::ostream& operator<<(std::ostream& cout, FeatureVector &fv);

FeatureVector operator +(const FeatureVector &a, const FeatureVector &b);
FeatureVector operator -(const FeatureVector &a, const FeatureVector &b);
double operator *(const FeatureVector &a, const FeatureVector &b);
FeatureVector operator *(const FeatureVector &a, const double &b);
#endif
