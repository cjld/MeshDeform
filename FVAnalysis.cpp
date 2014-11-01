#include "FVAnalysis.h"

using namespace std;
using namespace Eigen;

int ckf(double x) {
    if (x>Eps) return 1; else
        if (x<-Eps) return -1;
    return 0;
}

FVAnalysis::FVAnalysis(ARAPDeform *deform) {
    this->fvs.resize(deform->fvs.size());
    for (int i=0; i<fvs.size(); i++) {
        cout << "Copy " << i << endl;
        fvs[i] = deform->fvs[i];
    }
    this->deform = deform;
    deform->maxIterTime = 20;
    deform->iterEps = 0;
    outputMesh = *(deform->ref.mesh);
}

void FVAnalysis::ckRefFeatureVector() {
    FeatureVector &fv = fvs[0];
    for (int i=0; i<fv.s.size(); i++) {
        double rs = (fv.s[i] - Matrix3d::Identity()).norm();
        if (ckf(rs) != 0) {
            cout << "!!Error fv0.s " << i << " not identity " << rs << endl;
        }
    }
    for (int i=0; i<fv.logdr.size(); i++) {
        double rs = fv.logdr[i].norm();
        if (ckf(rs) != 0) {
            cout << "!!Error fv0.logdr " << i << " not zero " << rs << endl;
        }
    }
}

void FVAnalysis::writeMesh(std::string name, int id) {
    id += 1000;
    stringstream ss;
    ss << name << id << ".obj";
    cout << "WriteMesh " << ss.str() << endl;
    OpenMesh::IO::write_mesh(outputMesh, ss.str().c_str());
}

void FVAnalysis::solveFV(FeatureVector &fv) {
    FeatureVector fv2 = fv + fvs[0];
    fv2.calcRFromLogR();
    deform->solve(fv2, outputMesh);
}

void FVAnalysis::work() {

    freopen("log.txt","w",stdout);
    this->ckRefFeatureVector();
    for (int i=0; i<fvs.size(); i++)
        fvs[i].calcLogRFromR();
    for (int i=1; i<fvs.size(); i++)
        fvs[i] = fvs[i] - fvs[0];
    for (int i=1; i<fvs.size(); i++) {
        double mi = 1e300, ml = 1e300, mr = -ml;
        int mid = 0;
        for (int j=1; j<fvs.size(); j++) {
            double l = 1e300 ,r = -1e300;
            double rs = calcRs(fvs[j], l, r);
            cout << "calc main axis : " << j << " rs : " << rs;
            cout << " l : " << l << " r : " << r << endl;
            if (rs < mi) {
                mi = rs, mid = j;
                ml = l, mr = r;
            }
        }
        cout << "find axis " << i << " id : "  << mid << " rs : " << mi;
        cout << " l : " << ml << " r : " << mr << endl;
        cerr << "find axis " << i << " id : "  << mid << " rs : " << mi;
        cerr << " l : " << ml << " r : " << mr << endl;
        axis.push_back(fvs[mid]);
        region.push_back(make_pair(ml, mr));
        FeatureVector &fv = axis[i-1];
        this->makeOrtho(fv);
        this->solveFV(fv);
        this->writeMesh("axis",i);
    }
}

double FVAnalysis::calcRs(FeatureVector &ax, double &l, double &r) {
    double mo = ax * ax;
    cout << "mo : " << mo << endl;
    mo = sqrt(mo);
    if (ckf(mo)==0) return 1e30;
    double ans=0;
    for (int i=1; i<fvs.size(); i++) {
        double dot = fvs[i] * ax / mo;
        l = min(l,dot / mo);
        r = max(r,dot / mo);
        double mo2 = fvs[i] * fvs[i];
        ans += mo2 - dot*dot;
    }
    return ans;
}

void FVAnalysis::makeOrtho(FeatureVector &ax) {
    double mo = ax * ax;
    double ckrs = 0;
    for (int i=1; i<fvs.size(); i++) {
        double dot = fvs[i] * ax;
        dot /= mo;
        fvs[i] = fvs[i] - (ax * dot);
        ckrs += fvs[i] * fvs[i];
        cout << "cut " << i << " " << dot << endl;
    }
    cout << "ckrs : " << ckrs << endl;
}
