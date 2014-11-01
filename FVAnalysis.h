#ifndef FVANALYSIS_H
#define FVANALYSIS_H

#include "ARAPDeform.h"

class FVAnalysis {
public:
    FVAnalysis(ARAPDeform* deform);
    void ckRefFeatureVector();
    void work();
    double calcRs(FeatureVector &ax, double &l, double &r);
    void makeOrtho(FeatureVector &ax);
    void writeMesh(std::string name, int id);
    void solveFV(FeatureVector &fv);

    std::vector<FeatureVector> fvs, axis;
    std::vector< std::pair<double,double> > region;
    ARAPDeform *deform;
    DTriMesh outputMesh;
};

#endif // FVANALYSIS_H
