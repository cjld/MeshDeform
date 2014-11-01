#include "Align.h"
#include <Eigen/Geometry>
#include <iostream>
#include <fstream>
#include <cassert>
#include <Eigen/SVD>
#include <unsupported/Eigen/MatrixFunctions>
using namespace Eigen;
using namespace std;

ofstream fout("debug.txt");


Eigen::Matrix3d exp(Eigen::Matrix3d x) {
    //return x.exp();
    double theta = sqrt(x(0,1)*x(0,1) + x(0,2)*x(0,2) + x(1,2)*x(1,2));
    if (abs(theta) == 0) return Eigen::Matrix3d::Identity();
    x /= theta;
    return Eigen::Matrix3d::Identity() +
            x * sin(theta) +
            x * x * (1 - cos(theta));
}

Eigen::Matrix3d log(Eigen::Matrix3d x) {
    //return x.log();
    double theta = (x.trace() - 1) / 2;
    theta = acos(max(-1.0, min(1.0, theta)));
    if (abs(theta) == 0) return Eigen::Matrix3d::Zero();
    return (theta / (2 * sin(theta))) * (x - x.transpose());
}

Vector3d crossNorm(Vector3d va, Vector3d vb, double len) {
    va = va.cross(vb);
    if (va.norm() < Eps) return va;
    return va * (len / va.norm());
}
/*
Vector3d genNorm(Vector3d v) {
    Vector3d p = v;
    int i=0;
    while (i<2 && abs(p(i))<Eps) i++;
    p(i) = 0;
    return crossNorm(v, p);
}
*/
double polarDec(const Eigen::Matrix3d &a, Eigen::Matrix3d &r, Eigen::Matrix3d &s) {
    JacobiSVD<Eigen::MatrixXd> svd(a, ComputeThinU | ComputeThinV);
    r = svd.matrixU() * svd.matrixV().transpose();
    s = svd.matrixV() * svd.singularValues().asDiagonal() * svd.matrixV().transpose();
/*
    if (flag) {
    cout << "\nSVD dec : " << endl;
    cout << "A : " << endl << a << endl;
    cout << "U : " << endl;
    cout << svd.matrixU() << endl;
    cout << "D : " << endl;
    cout << svd.singularValues() << endl;
    cout << "V : " << endl;
    cout << svd.matrixV() << endl;
    }
*/
    if (r.determinant()<0) {
        Vector3d sv = svd.singularValues();
        int minsv = 0;
        r = svd.matrixU();
        if (sv(1) < sv(minsv)) minsv = 1;
        if (sv(2) < sv(minsv)) minsv = 2;
        if (sv(minsv) < -Eps) {
            cerr << "polar dec Error, min singular values <= 0 :" << endl;
            cerr << a << endl;
        }
        //cout << "Min SV " << sv(minsv) << " " << minsv << endl;
        r.col(minsv) *= -1;
        sv(minsv) *= -1;
        //cout << "R :\n" << r << endl;
        r = r * svd.matrixV().transpose();
        s = svd.matrixV() * sv.asDiagonal() * svd.matrixV().transpose();
        return sv.sum();
    }
    return svd.singularValues().sum();
}


double polarDec(const Eigen::Matrix3d &a, Eigen::Matrix3d &r) {
    JacobiSVD<Eigen::MatrixXd> svd(a, ComputeThinU | ComputeThinV);
    r = svd.matrixU() * svd.matrixV().transpose();
    if (r.determinant()<0) {
        //cout << "Determinant of R = -1" << endl;
        Vector3d sv = svd.singularValues();
        int minsv = 0;
        r = svd.matrixU();
        if (sv(1) < sv(minsv)) minsv = 1;
        if (sv(2) < sv(minsv)) minsv = 2;
        if (sv(minsv) < -Eps) {
            cerr << "polar dec Error, min singular values <= 0 :" << endl;
            cerr << a << endl;
        }
        //cout << "Min SV " << sv(minsv) << " " << minsv << endl;
        r.col(minsv) *= -1;
        sv(minsv) *= -1;
        //cout << "R :\n" << r << endl;
        r = r * svd.matrixV().transpose();
        return sv.sum();
    }
    return svd.singularValues().sum();
}


AffineAlign::AffineAlign(std::vector<Eigen::Vector3d> &v) : p(v), AtA(Matrix3d::Zero()) {
    if (flag) fout << p.size() << endl;
    for (int i=0; i<p.size(); i++) {
        AtA += p[i] * p[i].transpose();
        if (flag) fout << p[i] << endl << endl;
    }
    if (flag) {
        fout << "AtA" << endl;
        fout << AtA << endl;
        fout << AtA.determinant() << endl;
    }
    AtA = AtA.inverse().eval();
}

Matrix3d AffineAlign::calc(const std::vector<Vector3d> &v) {
    if (v.size() != p.size()) {
        cout << "!!Error v.size() != p.size()" << endl;
        fout << "!!Error v.size() != p.size()" << endl;
    }
    Vector3d vx(Vector3d::Zero()), vy(Vector3d::Zero()), vz(Vector3d::Zero());
    for (int i=0; i<p.size(); i++) {
        vx += p[i]*v[i](0);
        vy += p[i]*v[i](1);
        vz += p[i]*v[i](2);
    }
    vx = AtA * vx;
    vy = AtA * vy;
    vz = AtA * vz;
    Matrix3d res;
    res << vx, vy, vz;
    return res.transpose();
}

double AffineAlign::residual(Eigen::Matrix3d m, std::vector<Eigen::Vector3d> v) {
    double rs = 0;
    for (int i=0; i<v.size(); i++)
        rs += (v[i]-m*p[i]).squaredNorm();
    return rs;
}

bool AffineAlign::ckt(Eigen::Matrix3d m, std::vector<Eigen::Vector3d> v) {
    double rs = residual(m,v);
    double eps= 1e-3;
    for (int i=0; i<3; i++) for (int j=0; j<3; j++)
        for (int f=-1; f<=1; f+=2) {
            Matrix3d mat = m;
            mat(i,j) += f*eps;
            double ts = residual(mat,v);
            cout << i<< "," << j << " : " << ts-rs << endl;
            if (ts<rs-Eps)
                return false;
        }
    return true;
}


RotateAlign::RotateAlign(std::vector<Vector3d> &v1, std::vector<Vector3d> &v2) : p(&v1), q(&v2) {
    assert(p->size() == q->size());
}

Matrix3d RotateAlign::calc() {
    res = 0;
    Matrix3d s(Matrix3d::Zero());
    for (int i=0; i<q->size(); i++) {
        res += (*q)[i].squaredNorm();
        res += (*p)[i].squaredNorm();
        s += (*q)[i] * (*p)[i].transpose();
    }
    Matrix3d r,y;
    res -= 2*polarDec(s,r,y);
    return r.transpose();
}

double RotateAlign::residual(const Matrix3d &m) {
    double rs=0;
    for (int i=0; i<q->size(); i++)
        rs += ((*p)[i] - m*(*q)[i]).squaredNorm();
    return rs;
}

bool RotateAlign::ckt(Matrix3d m) {
    double rs = residual(m);
    double eps= 1e-3;
    cout << "residual : " << endl;
    cout << rs << endl;
    cout << this->res << endl;
    cout << this->res-rs << endl;
    for (int i=0; i<3; i++) for (int j=0; j<3; j++)
        for (int f=-1; f<=1; f+=2) {
            Matrix3d mat = m;
            mat(i,j) += f*eps;
            Matrix3d a,b;
            polarDec(mat,a,b);
            double ts = residual(a);
            cout << i<< "," << j << " : " << ts << " " << ts-rs << " " << a.determinant() << endl;
            if (ts<rs-Eps) {
                cout << "Error!" << endl;
                return false;
            }
        }
    return true;
}


void RotateAlign::AlignAtoB(DTriMesh &A, DTriMesh &B) {
    vector<Vector3d> v1, v2;
    Vector3d vs1(Vector3d::Zero()), vs2(Vector3d::Zero());
    for (int i=0; i<A.n_vertices(); i++) {
        DTriMesh::VertexHandle vi(i);
        DTriMesh::Point pa = A.point(vi), pb = B.point(vi);
        Vector3d va(pa[0], pa[1], pa[2]), vb(pb[0], pb[1], pb[2]);
        v1.push_back(va);
        v2.push_back(vb);
        vs1 += va, vs2 += vb;
    }
    vs1 /= A.n_vertices();
    vs2 /= A.n_vertices();

    for (int i=0; i<v1.size(); i++) {
        v1[i] -= vs1;
        v2[i] -= vs2;
    }

    RotateAlign align(v2,v1);
    Matrix3d rot = align.calc();

    for (int i=0; i<v1.size(); i++) {
        v1[i] = rot * v1[i] + vs2;
        DTriMesh::Point p(v1[i](0), v1[i](1), v1[i](2));
        A.point(DTriMesh::VertexHandle(i)) = p;
    }
}

void RotateAlign::AlignAtoBCenter(DTriMesh &A, DTriMesh &B) {
	vector<Vector3d> v1, v2;
	Vector3d vs1(Vector3d::Zero()), vs2(Vector3d::Zero());
	for (int i=0; i<A.n_vertices(); i++) {
		DTriMesh::VertexHandle vi(i);
		DTriMesh::Point pa = A.point(vi), pb = B.point(vi);
		Vector3d va(pa[0], pa[1], pa[2]), vb(pb[0], pb[1], pb[2]);
		v1.push_back(va);
		v2.push_back(vb);
		vs1 += va, vs2 += vb;
	}
	vs1 /= A.n_vertices();
	vs2 /= A.n_vertices();

	for (int i=0; i<v1.size(); i++) {
		v1[i] -= vs1;
		v2[i] -= vs2;
	}

	RotateAlign align(v2,v1);
	Matrix3d rot = align.calc();

	for (int i=0; i<v1.size(); i++) {
		v1[i] = rot * v1[i]+vs2;
		DTriMesh::Point p(v1[i](0), v1[i](1), v1[i](2));
		A.point(DTriMesh::VertexHandle(i)) = p;
	}
}


OpenMesh::Vec3d EtoO(const Eigen::Vector3d &v) {
    return OpenMesh::Vec3d(v(0),v(1),v(2));
}

Eigen::Vector3d OtoE(const OpenMesh::Vec3d &v) {
    return Eigen::Vector3d(v[0],v[1],v[2]);
}
