#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <ctime>
#include <iostream>
#include <tuple>

#include "deformsf.h"

#include "ARAPDeform.h"
#include "FVAnalysis.h"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <unsupported/Eigen/MatrixFunctions>
#include "DeformCaller.h"

using namespace std;
using namespace Eigen;

int flag = 0;

#ifndef DBUILD_LIB

void ckAlign() {
    srand(time(0));
    Matrix3d t = Matrix3d::Random();
    Matrix3d r,s;
    polarDec(t,r,s);
    vector<Vector3d> v1, v2;
    for (int i=0; i<2; i++) {
        Vector3d v(Vector3d::Random());
        v1.push_back(v);
        //v2.push_back(r*v);
        v2.push_back(Vector3d::Random());
    }
    RotateAlign align(v1, v2);
    Matrix3d rr = align.calc();

    cout << rr << endl;
    cout << r << endl;
    cout << (r-rr.transpose()).norm() << endl;

    align.ckt(rr);
}

Matrix3d getRandRot() {
    srand(time(0));
    Matrix3d x(Matrix3d::Random()),t1,t2;
    polarDec(x,t1,t2);
    return t1;
}

void ckMatrixLogExp() {
    Matrix3d x = getRandRot();
//    /x << -1,0,0,0,-1,0,0,0,1;
    cout << x << endl << endl;
    cout << log(x) << endl << endl;
    cout << exp(log(x)) << endl << endl;
    cout << exp(log(x)) - x << endl << endl;


    cout << x.log() << endl << endl;
    cout << x.log().exp() << endl << endl;
    cout << x.log().exp() - x << endl << endl;

}

void testConvexOptimize() {
    TriSecSolver s1;
    NewtonIterSolver s2;
    TriNewtonIterSolver s3;

    s1.eps = s2.eps = s3.eps = 1e-6;
    s3.convexFunction = s2.convexFunction = s1.convexFunction = [](double x) -> double {x = abs(x - 10); return pow(x,2.123)+1;};

    double best, num;
    double l = 0, r = 100;
    s1.solve(l, r, best, num);
    cout << endl;
    cout << "TriSecSolver : " << endl;
    cout << "   iterTime : " << s1.iterTime;
    cout << "   best : " << best << "\n   num : " << num << endl;
    cout << endl;

    s2.solve(l, r, best, num);
    cout << endl;
    cout << "NewtonSolver : " << endl;
    cout << "   iterTime : " << s2.iterTime;
    cout << "   best : " << best << "\n   num : " << num << endl;
    cout << endl;

    s3.solve(l, r, best, num);
    cout << endl;
    cout << "TriNewtonSolver : " << endl;
    cout << "   iterTime : " << s3.iterTime;
    cout << "   best : " << best << "\n   num : " << num << endl;
    cout << endl;
}

void testFVA() {
    DMEngine eng;
    vector<DTriMesh*> ms;
    string path = "E:/project/SIGA2014/dataset/horse/";
    int num = 49;
    cout << "Input Path : " << endl;
    cin >> path;
    cout << "Input Number : " << endl;
    cin >> num;
    for (int i=1; i <= num; i++) {
        DTriMesh *mesh = new DTriMesh();
        stringstream ss;
        ss << path << "/" << i << ".obj";
        cout << "Read Mesh : " << ss.str() << endl;
        OpenMesh::IO::read_mesh(*mesh, ss.str().c_str());
        ms.push_back(mesh);
    }
    cout << "Ref Id : " << endl;
    int x;
    cin >> x;
    swap(ms[0], ms[x]);
    ARAPDeform deform(eng, *ms[0], ms);
    cout << "ARAP construct done" << endl;
    FVAnalysis fva(&deform);
    cout << "FVA construct done" << endl;
    fva.work();
}

#define FORE(i,v) for (decltype(v.begin()) i = v.begin(); i!=v.end(); i++)

void FuckFaceMesh() {
    string objname = "E:\\SIGA2014\\dataset\\faces\\1.obj";
    map<int,bool> has;
    ifstream cin(objname);
    char x;
    string ha;
    getline(cin,ha);
    cout << ha;
    vector< tuple<int,int,int> > fs;
    ofstream cout("result.obj");
    int t=1;
    map<int,Vector3d> pp;
    while (cin>>x) {
        if (x != 'f') {
            double a,b,c;
            cin >> a >> b >> c;
            pp[t++] = Vector3d(a,b,c);
        } else {
            int a,b,c;
            cin >> a >> x >> x >> b >> x >> x >> c >> x >> x;
            has[a] = has[b] = has[c] = true;
            fs.push_back(make_tuple(a,b,c));
        }
    }
    map<int,int> rid;
    t=1;
    FORE(x,has) {
        int i = x->first;
        rid[i] = t++;
    }
    FORE(x,has) {
        int i = x->first;
        cout << "v " << pp[i].transpose() << endl;
    }
    FORE(x,fs) {
        cout << "f " << rid[get<0>(*x)] << " " << rid[get<1>(*x)] << " " << rid[get<2>(*x)] << endl;
    }
}

int main() {
	DeformSF::main();
    //CallMain();
    //testFVA();
    //return 0;
    DMEngine eng;
    DTriMesh mesh1, mesh2, mesh3;
    cout << "load mesh" << endl;
    //OpenMesh::IO::read_mesh(mesh1, "C:/Users/CJLD/Desktop/Lab/bars/2.obj");
    //OpenMesh::IO::read_mesh(mesh2, "C:/Users/CJLD/Desktop/Lab/bars/4.obj");
    //OpenMesh::IO::read_mesh(mesh3, "C:/Users/CJLD/Desktop/Lab/bars/6.obj");

    //OpenMesh::IO::read_mesh(mesh1, "E:/project/SIGA2014/dataset/scape/49.obj");
    //OpenMesh::IO::read_mesh(mesh2, "E:/project/SIGA2014/dataset/scape/50.obj");
    //OpenMesh::IO::read_mesh(mesh3, "E:/project/SIGA2014/dataset/scape/19.obj");

    OpenMesh::IO::read_mesh(mesh1, "E:/project/scape/cat/alignmodels/3.obj");
    OpenMesh::IO::read_mesh(mesh2, "E:/project/scape/cat/alignmodels/6.obj");

    //OpenMesh::IO::read_mesh(mesh1, "C:/Users/CJLD/Desktop/Lab/bar2/bar.obj");
    //OpenMesh::IO::read_mesh(mesh2, "C:/Users/CJLD/Desktop/Lab/bar2/bar_deform2.obj");

    //OpenMesh::IO::read_mesh(mesh1, "C:/Users/CJLD/Desktop/Lab/box3D/box3D.obj");
    //OpenMesh::IO::read_mesh(mesh2, "C:/Users/CJLD/Desktop/Lab/box3D/box3D_90.obj");
    //OpenMesh::IO::read_mesh(mesh2, "C:/Users/CJLD/Desktop/Lab/ARAP/build-ARAP-VS2012x64-Debug/t4.obj");

    long long tt = clock();
#ifdef DRELEASE
    freopen("log.txt","w",stdout);
#endif

    vector<DTriMesh*> ms;
    ms.push_back(&mesh1);
    ms.push_back(&mesh2);
    cout << "calc feature vector" << endl;
    ARAPDeform deform(eng, mesh1, ms);
    deform.needAlign = true;
    deform.maxIterTime = 0;
    deform.iterEps = 0;
    vector<double> weight;
    double t = 0.671;
    weight.push_back(1-t);
    weight.push_back(t);
    DTriMesh result(mesh1);
    cout << "solve" << endl;

    //FeatureVector fv = deform.fvs[0];
    //fv.blendFrom(weight, deform.fvs);
    //fv.loadConstPoint(ifstream("C:/Users/CJLD/Desktop/Lab/scape/constPoint/c6.txt"));
    //deform.solve2(fv, weight, result);
    deform.solve(weight, result);
    //deform.solve(fv, result);

    cerr << "Total Time : " << (clock() - tt) * 1.0 / CLOCKS_PER_SEC << endl;

    cout << "write mesh" << endl;
    OpenMesh::IO::write_mesh(result, "result.obj");

    return 0;
}

#endif
