#include "FeatureVector.h"
#include <unsupported/Eigen/MatrixFunctions>
#include <cmath>
#include <iostream>

#include <omp.h>

#ifdef DUSE_OPENMP
#define DOMP_END \
}
#else
#define DOMP_END ;
#endif
using namespace std;

double RefMesh::normalScale = 0.3;
int RefMesh::root = 0;

RefMesh::~RefMesh() {
    for (int i=0; i<align.size(); i++)
        delete align[i];
}

double cotan(Eigen::Vector3d a, Eigen::Vector3d b) {
    double na = a.norm(), nb = b.norm();
    if (na<Eps || nb<Eps) return 0;
    double cos = a.dot(b) / (na*nb);
    if (cos == 1) return 1;
    return cos / sqrt(1-cos*cos);
}

double tan2(Eigen::Vector3d a, Eigen::Vector3d b) {
    double na = a.norm(), nb = b.norm();
    if (na<Eps || nb<Eps) return 0;
    double cos = a.dot(b) / (na*nb);
    double theta = acos(cos)/2;
    double ans = tan(theta);
    if (ans>=0 && ans<=100) return ans / nb;
    return 1 / nb;
}

double sexp(double x) {
    if (x <= 0) return exp(x);
    return 1+x;
}

RefMesh::RefMesh(DTriMesh &ms) : mesh(&ms) {
    usec = false, fvNum = 0;
    ms.update_face_normals();
    ms.update_vertex_normals();
    align.resize(mesh->n_vertices());
    d.resize(mesh->n_vertices());
    rd.resize(mesh->n_vertices());
    vvs=0;
    for (int i=0; i<mesh->n_vertices(); i++) {
        DTriMesh::VertexHandle vi(i);
        DTriMesh::Point p = mesh->point(vi);
        std::vector<Eigen::Vector3d> vec;
        d[i] = vvs;
        double lens = 0;
        for (DTriMesh::VertexVertexIter vvi = mesh->vv_iter(vi); vvi; vvi++) {
            rd[vvi.handle().idx()].push_back(std::make_pair(i, vvs));
            Eigen::Vector3d q = OtoE(mesh->point(vvi) - p);
            vec.push_back(q);
            vvs++;
            lens += q.norm();
        }
        vvs = d[i];
        double ss = 0;
        for (DTriMesh::VertexVertexIter vvi = mesh->vv_iter(vi); vvi; vvi++) {
            int j = vvs - d[i];
            int pre = j==0 ? (int)vec.size()-1 : j-1;
            int next = j+1==vec.size() ? 0 : j+1;
/*
            w.push_back(sqrt( std::max(Eps, 0.5*(
                        cotan(vec[pre],vec[pre]-vec[j]) +
                        cotan(vec[next],vec[next]-vec[j])
            ))));
*/
#ifndef DUSE_TAN
            w.push_back( sqrt( sexp( 0.5*(
                        cotan(vec[pre],vec[pre]-vec[j]) +
                        cotan(vec[next],vec[next]-vec[j])
            ))));
            /*
            if (vec[pre].dot(vec[j]) == vec[pre].norm() * (vec[pre]-vec[j]).norm()) {
                cout << "pl i : " << i << " j : " << j << endl;
            }
            */
#else

            w.push_back(tan2(vec[pre],vec[j]) +
                        tan2(vec[next],vec[j]));
#endif
            //w.push_back(1);

            /*
            w.push_back( sqrt( abs (0.5*(
                        cotan(vec[pre],vec[pre]-vec[j]) +
                        cotan(vec[next],vec[next]-vec[j])
            ))));
            */

            ss += w[w.size()-1];
            vvs++;
        }
        for (int j=d[i]; j<w.size(); j++)
            w[j] = w[j] / ss;

        // add normal
        vec.push_back( OtoE(mesh->normal(vi) * (lens / vec.size() * RefMesh::normalScale)) );
        align[i] = new AffineAlign(vec);
    }
    c.resize(w.size(), 0);
    std::vector<bool> visit(mesh->n_vertices(), false);
    int qBeg = 0;
    for (int i=root; i<visit.size(); i++) if (!visit[i]) {
        visit[i] = true;
        bfsq.push_back(make_pair(i,make_pair(-1,-1) ));
        while (qBeg < bfsq.size()) {
            int i = bfsq[qBeg++].first;

            DTriMesh::VertexHandle vi(i);
            int vvs = d[i];

            for (DTriMesh::VertexVertexIter vvi = mesh->vv_iter(vi); vvi; vvi++, vvs++) {
                int j = vvi.handle().idx();
                if (!visit[j]) {
                    visit[j] = true;
                    bfsq.push_back(make_pair( j, make_pair(i,vvs) ));
                }
            }
        }
    }
}

void RefMesh::getFeature(DTriMesh &ms, FeatureVector &fv) {

    static int callTime = -1;
    callTime ++;


    ms.update_face_normals();
    ms.update_vertex_normals();

    fv.s.resize(align.size());
    fv.r.resize(align.size());
    fv.isConst.resize(align.size(), 0);
    fv.constPoint.resize(align.size());
/*
    ofstream fout;
    if (callTime == 1) {
        fout.open("color.txt");
        fout << ms.n_vertices() << endl;
    }
    */

    //ofstream cout((string("fvv") + (char)(callTime+'0')).c_str());

    for (int i=0; i<ms.n_vertices(); i++) {
        DTriMesh::VertexHandle vi(i);
        DTriMesh::Point p = ms.point(vi);
        std::vector<Eigen::Vector3d> vec;
        double lens = 0;
        for (DTriMesh::VertexVertexIter vvi = ms.vv_iter(vi); vvi; vvi++) {
            Eigen::Vector3d q = OtoE(ms.point(vvi) - p);
            vec.push_back(q);
            lens += q.norm();
        }
        // add normal
        vec.push_back( OtoE(ms.normal(vi) * (lens / vec.size() * RefMesh::normalScale)) );
        if (vec.size() <= 1) {
            fv.s[i] = fv.r[i] = Eigen::Matrix3d::Identity();
            continue;
        }
        Eigen::Matrix3d mat = align[i]->calc(vec);
        polarDec(mat, fv.r[i], fv.s[i]);
        /*
        if (i == 12049 && callTime < 2) {
            stringstream ss;
            ss << "tt" << callTime << ".txt";
            ofstream fout(ss.str().c_str());
            fout << "lens : " << lens << endl;
            for (int i=0; i<vec.size(); i++)
                fout << vec[i].transpose() << endl;
        }
        if (fv.s[i].determinant() < 0) {
            fout << 1 << endl;
            cout << "~Warning fv.s " << i << " < 0 id : " << callTime << "\n";
        } else fout << 0 << endl;
        */
    }

    fv.dr.resize(vvs);
    fv.logdr.resize(vvs);

    if (usec)
        this->fvNum++;

    for (int i=0; i<ms.n_vertices(); i++) {
        DTriMesh::VertexHandle vi(i);
        int vs = d[i];
        for (DTriMesh::VertexVertexIter vvi = ms.vv_iter(vi); vvi; vvi++) {
            fv.dr[vs] = fv.r[i].transpose() * fv.r[vvi.handle().idx()];
            fv.logdr[vs] = log(fv.dr[vs]);
            if (usec)
                c[vs] = (c[vs] * (this->fvNum-1) + fv.logdr[vs].norm() + (fv.s[i] - Eigen::Matrix3d::Identity()).norm() ) / this->fvNum;
            vs++;
        }
    }
}

FeatureVector::FeatureVector(std::vector<double> weight, std::vector<FeatureVector> &fvs) {
    s.resize(fvs[0].s.size());
    r.resize(fvs[0].s.size());
    dr.resize(fvs[0].dr.size());
    logdr.resize(fvs[0].dr.size());
    isConst.resize(s.size(), 0);
    constPoint.resize(s.size());
    this->blendFrom(weight, fvs);
}


void FeatureVector::blendFrom(std::vector<double> weight, std::vector<FeatureVector> &fvs) {


#ifdef DUSE_OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
    for (int i=0; i<s.size(); i++) {
        r[i] = Eigen::Matrix3d::Zero();
        Eigen::Matrix3d mat(Eigen::Matrix3d::Zero());
        for (int j=0; j<weight.size(); j++) {
            mat += weight[j] * fvs[j].s[i];
            r[i] += weight[j] * log(fvs[j].r[i]);
        }
        r[i] = exp(r[i]);
        s[i] = mat;
    }
DOMP_END;

#ifdef DUSE_OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
    for (int i=0; i<dr.size(); i++) {
        Eigen::Matrix3d mat(Eigen::Matrix3d::Zero());
        for (int j=0; j<weight.size(); j++)
            mat += weight[j] * fvs[j].logdr[i];
        logdr[i] = mat;
        dr[i] = exp(mat);
    }
DOMP_END;

}


std::ostream& operator<<(std::ostream& cout, FeatureVector &fv) {
    cout << "S : " << std::endl;
    for (int i=0; i<fv.s.size(); i++) {
        cout << i << std::endl << fv.s[i] << std::endl;
        cout << "det : " << fv.s[i].determinant() << endl;
    }
    cout << "dR : " << std::endl;
    for (int i=0; i<fv.dr.size(); i++) {
        cout << i << std::endl << fv.dr[i] << std::endl;
        cout << "det : " << fv.dr[i].determinant() << endl;
    }
    cout << "log(dR) : " << std::endl;
    for (int i=0; i<fv.logdr.size(); i++)
        cout << i << std::endl << fv.logdr[i] << std::endl;
    return cout;
}

void FeatureVector::setConstPoint(int i, Eigen::Vector3d v) {
    isConst[i] = true;
    constPoint[i] = v;
}

void FeatureVector::loadConstPoint(std::istream& cin) {
    int n;
    cin >> n;
    for (int i=0; i<n; i++) {
        int m;
        cin >> m;
        std::vector<int> ids(m);
        for (int i=0; i<m; i++)
            cin >> ids[i];
        for (int i=0; i<m; i++) {
            double x,y,z;
            cin >> x >> y >> z;
            this->setConstPoint(ids[i], Eigen::Vector3d(x,y,z));
        }
    }
}

void FeatureVector::loadConstPoint(std::istream& cin, DTriMesh& mesh) {
    int n;
    cin >> n;
    for (int i=0; i<n; i++) {
        int m;
        cin >> m;
        std::vector<int> ids(m);
        for (int i=0; i<m; i++)
            cin >> ids[i];
        for (int i=0; i<m; i++) {
            double x,y,z;
            cin >> x >> y >> z;
            this->setConstPoint(ids[i], OtoE(mesh.point(DTriMesh::VertexHandle(ids[i]))));
        }
    }
}

double RefMesh::getWij(int e) {
    //if (fvNum) return exp(-c[e]);// * w[e];
    return w[e];
}

void RefMesh::printCij(std::ofstream &fout) {
    std::vector<double> color;
    for (int i=0; i<mesh->n_vertices(); i++) {
        DTriMesh::VertexHandle vi(i);
        int vs = d[i];
        int s=0;
        double w=0;
        for (DTriMesh::VertexVertexIter vvi = mesh->vv_iter(vi); vvi; vvi++) {
            w += this->getWij(vs);
            s++;
            vs++;
        }
        color.push_back(w/s);
    }
    double mx = color[0], mi = mx;
    for (int i=1; i<color.size(); i++) {
        mx = std::max(mx, color[i]);
        mi = std::min(mi, color[i]);
    }
    fout << color.size() << endl;
    for (int i=0; i<color.size(); i++) {
        fout << (color[i]-mi) / (mx-mi) << '\n';
    }
    std::cout << "Color Max : " << mx << endl;
    std::cout << "Color Min : " << mi << endl;
    std::cout << "Color Max-Min : " << mx - mi << endl;
    fout << endl;
}

void FeatureVector::calcLogRFromR() {
    logr.resize(r.size());
    for (int i=0; i<logr.size(); i++)
        logr[i] = log(r[i]);
}

void FeatureVector::calcRFromLogR() {
    r.resize(logr.size());
    dr.resize(logdr.size());
    for (int i=0; i<logr.size(); i++)
        r[i] = exp(logr[i]);
    for (int i=0; i<logdr.size(); i++)
        dr[i] = exp(logdr[i]);
    isConst.resize(s.size(), false);
}

void FeatureVector::resize(const FeatureVector &other) {
    s.resize(other.s.size());
    logdr.resize(other.logdr.size());
    logr.resize(other.logr.size());
}

FeatureVector operator +(const FeatureVector &a, const FeatureVector &b) {
    FeatureVector ans;
    ans.resize(a);
    for (int i=0; i<ans.s.size(); i++) {
        ans.s[i] = a.s[i] + b.s[i];
        ans.logr[i] = a.logr[i] + b.logr[i];
    }
    for (int i=0; i<ans.logdr.size(); i++) {
        ans.logdr[i] = a.logdr[i] + b.logdr[i];
    }
    return ans;
}

FeatureVector operator -(const FeatureVector &a, const FeatureVector &b) {
    FeatureVector ans;
    ans.resize(a);
    for (int i=0; i<ans.s.size(); i++) {
        ans.s[i] = a.s[i] - b.s[i];
        ans.logr[i] = a.logr[i] - b.logr[i];
    }
    for (int i=0; i<ans.logdr.size(); i++) {
        ans.logdr[i] = a.logdr[i] - b.logdr[i];
    }
    return ans;
}


double dot(const Eigen::Matrix3d &a, const Eigen::Matrix3d &b) {
    return
            a(0,0) * b(0,0) + a(0,1) * b(0,1) + a(0,2) * b(0,2) +
            a(1,0) * b(1,0) + a(1,1) * b(1,1) + a(1,2) * b(1,2) +
            a(2,0) * b(2,0) + a(2,1) * b(2,1) + a(2,2) * b(2,2);
}

double operator *(const FeatureVector &a, const FeatureVector &b) {
    double ans = 0;
    for (int i=0; i<a.s.size(); i++) {
        ans += dot(a.s[i], b.s[i]);
    }
    for (int i=0; i<a.logdr.size(); i++) {
        ans += dot(a.logdr[i], b.logdr[i]);
    }
    return ans;
}

FeatureVector operator *(const FeatureVector &a, const double &b) {
    FeatureVector ans;
    ans.resize(a);
    for (int i=0; i<ans.s.size(); i++) {
        ans.s[i] = a.s[i] * b;
        ans.logr[i] = a.logr[i] * b;
    }
    for (int i=0; i<ans.logdr.size(); i++) {
        ans.logdr[i] = a.logdr[i] * b;
    }
    return ans;
}
