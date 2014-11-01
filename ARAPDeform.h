#ifndef ARAP_DEFORM_H
#define ARAP_DEFORM_H

#include "DMEngine.h"
#include "FeatureVector.h"
#include "LQSolver.h"

class ARAPDeform {
public :
    DMEngine *eng;
    RefMesh ref;
    std::vector<DTriMesh*> meshs;
    std::vector<FeatureVector> fvs;
    LQSolver lqsolver;
    bool needAlign;
    bool needCalcRs;
    int maxIterTime;
    int newtonIterTime;
    double iterEps;
    double newtonIterEps;
    std::pair<double,double> rlimit;
    double SRp;
    bool CalcLocalFirst;

    bool AtAChanged;

    ARAPDeform(DMEngine &eng, DTriMesh &refMesh, std::vector<DTriMesh*> ms);

    // given feature vector, solve mesh
    void solve(FeatureVector fv, DTriMesh &mesh, bool predec = true);

    // given feature vector weights, solve mesh
    void solve(std::vector<double> weight, DTriMesh &mesh);

    // given some const point, try to blend fvs to  it
    void solve(FeatureVector fv, std::vector<double> &weight, DTriMesh &mesh);
    virtual void solve2(FeatureVector fv, std::vector<double> &weight, DTriMesh &mesh);

    // given feature vector,  try to blend fvs to  it
    void solve(FeatureVector &fv, std::vector<double> &weight);

    // given mesh, try to blend fvs to it
    void solve(DTriMesh &mesh, std::vector<double> &weight);

    void preDecMatrix(FeatureVector &fv, DTriMesh &mesh);
    virtual double globalOptimize(FeatureVector &fv, DTriMesh &mesh);
    double globalOptimizeFast(FeatureVector &fv, DTriMesh &mesh);
    virtual double localOptimize(FeatureVector &fv, DTriMesh &mesh);
    void localOptimizeFast(FeatureVector &fv, DTriMesh &mesh);
    virtual double getWeightResidual(FeatureVector &fv, const std::vector<double> &weight, DTriMesh &mesh);
    std::vector<double> getGradient(FeatureVector &fv, const std::vector<double> &weight, DTriMesh &mesh, double &zeroResult);

    void dfsInit(int i, std::vector<int> &visit, FeatureVector &fv, DTriMesh &mesh);
    void bfsInit(int i, std::vector<int> &visit, FeatureVector &fv, DTriMesh &mesh);
    void bfsInit2(FeatureVector &fv, DTriMesh &mesh);
    void writeIterMesh(DTriMesh &mesh, std::string name, int id);
    void ckAlign(DTriMesh &mesh);
    void initWeight(std::vector<double> &weight);
    void initMesh(FeatureVector &fv, DTriMesh &mesh);
};

class SR_ARAPDeform : public ARAPDeform {
public:
    SR_ARAPDeform(DMEngine &eng, DTriMesh &refMesh, std::vector<DTriMesh*> ms)
        : ARAPDeform(eng, refMesh, ms) {}
    double globalOptimize(FeatureVector &fv, DTriMesh &mesh);
    double localOptimize(FeatureVector &fv, DTriMesh &mesh);
};

class T_ARAPDeform : public ARAPDeform {
public:
    ARAPDeform *org;
    DMEngine eng2;
    std::vector<bool> ev;


    T_ARAPDeform(DMEngine &eng, DTriMesh &refMesh, std::vector<DTriMesh*> ms);
    ~T_ARAPDeform() {delete org;}
    double getWeightResidual(FeatureVector &fv, const std::vector<double> &weight, DTriMesh &mesh);
    void solve2(FeatureVector fv, std::vector<double> &weight, DTriMesh &mesh);
};

class T2_ARAPDeform : public ARAPDeform {
public:
    ARAPDeform *org;
    DMEngine eng2;
    std::vector<bool> ev;


    T2_ARAPDeform(DMEngine &eng, DTriMesh &refMesh, std::vector<DTriMesh*> ms);
    ~T2_ARAPDeform() {delete org;}
    void solve2(FeatureVector fv, std::vector<double> &weight, DTriMesh &mesh);
};


void ARAPTest();
#endif
