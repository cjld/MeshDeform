#ifndef DALIGN_H
#define DALIGN_H

#include <Eigen/Eigen>
#include <vector>

#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#define Eps 1e-10

extern int flag;

struct DTraits : public OpenMesh::DefaultTraits {
    typedef OpenMesh::Vec3d Point; // use double-values points
    typedef OpenMesh::Vec3d Normal; // use double-values points
    VertexAttributes( OpenMesh::Attributes::Normal );
    FaceAttributes( OpenMesh::Attributes::Normal );
};

typedef OpenMesh::TriMesh_ArrayKernelT<DTraits> DTriMesh;

class AffineAlign {
public :
    std::vector<Eigen::Vector3d> p;
    Eigen::Matrix3d AtA;
    AffineAlign(std::vector<Eigen::Vector3d> &v);
    Eigen::Matrix3d calc(const std::vector<Eigen::Vector3d> &v);
    double residual(Eigen::Matrix3d m, std::vector<Eigen::Vector3d> v);
    bool ckt(Eigen::Matrix3d m, std::vector<Eigen::Vector3d> v);
};

class RotateAlign {
public :
    std::vector<Eigen::Vector3d> *p, *q;
    double res;

    RotateAlign(std::vector<Eigen::Vector3d>& v1, std::vector<Eigen::Vector3d>& v2);
    Eigen::Matrix3d calc();
    double residual(const Eigen::Matrix3d &m);
    bool ckt(Eigen::Matrix3d m);

    static void AlignAtoB(DTriMesh &A, DTriMesh &B);
    static void AlignAtoBCenter(DTriMesh &A, DTriMesh &B);
};

double polarDec(const Eigen::Matrix3d &a, Eigen::Matrix3d &r, Eigen::Matrix3d &s);
double polarDec(const Eigen::Matrix3d &a, Eigen::Matrix3d &r);

OpenMesh::Vec3d EtoO(const Eigen::Vector3d &v);
Eigen::Vector3d OtoE(const OpenMesh::Vec3d &v);

Eigen::Matrix3d exp(Eigen::Matrix3d);
Eigen::Matrix3d log(Eigen::Matrix3d);
#endif
