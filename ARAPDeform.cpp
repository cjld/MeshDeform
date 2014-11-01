#include "ARAPDeform.h"
#include <fstream>
#include <queue>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <sstream>


using namespace std;
using namespace Eigen;

#include <omp.h>

#ifdef DUSE_OPENMP
#define DOMP_END \
}
#else
#define DOMP_END ;
#endif

#define floatTime ((clock()-tt)*1.0 / CLOCKS_PER_SEC)

bool isNan(double fN) {
    return !(fN==fN);
}

bool isNan(Eigen::Matrix3d m) {
    for (int i=0; i<3; i++) for (int j=0; j<3; j++)
        if (isNan(m(i,j))) return true;
    return false;
}


bool isNan(Eigen::Vector3d v) {
    for (int i=0; i<3; i++)
        if (isNan(v(i))) return true;
    return false;
}

ARAPDeform::ARAPDeform(DMEngine &eng, DTriMesh &refMesh, std::vector<DTriMesh*> ms) : ref(refMesh), meshs(ms), lqsolver(&eng) {
    int cpunum = omp_get_num_procs();//cpunum=1;
    omp_set_num_threads(cpunum);
    SRp = 1;

    iterEps = 0.01;
    newtonIterEps = 1e-3;
    needAlign = false;
    maxIterTime = 100;
    newtonIterTime = 5;
    AtAChanged = true;
    needCalcRs = true;
    CalcLocalFirst = false;
    this->eng = &eng;
    rlimit = make_pair(0,1);
    fvs.resize(ms.size());
    //ofstream fout("fv.txt");

 /*
#ifdef DUSE_OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
*/

    for (int i=0; i<fvs.size(); i++) {
        printf("Calc %d feature vector \n", i);
        ref.getFeature(*meshs[i], fvs[i]);
        //fout << "F"<<i<<endl;
        //fout << fvs[i] << endl;
    }
//DOMP_END;
}

void ARAPDeform::initMesh(FeatureVector &fv, DTriMesh &mesh) {
    long long tt = clock();
    cout << "init feature vector ... ";
    bfsInit2(fv, mesh);
    cout << " time : " << floatTime << endl;
    return;
    vector<int> visit(fv.s.size(), 0);
    for (int i=0; i<visit.size(); i++)
        if (!visit[i]) {
            //fv.r[i] = fvs[0].r[i];
            bfsInit(i, visit, fv, mesh);
        }
}

void ARAPDeform::writeIterMesh(DTriMesh &mesh, string name, int id) {
#ifdef WRITE_ITER_MESH
    DTriMesh ms = mesh;
    RotateAlign::AlignAtoB(ms, *(ref.mesh));
    stringstream ss;
    ss << name << id+1000 << ".obj";
    cout << "Write mesh " + ss.str() << endl;
    OpenMesh::IO::write_mesh(ms, ss.str().c_str());
#endif
}

void ARAPDeform::solve(FeatureVector fv, DTriMesh &mesh, bool predec) {
    if (CalcLocalFirst) {
        double presr = this->SRp;
        this->SRp = 0;
        localOptimize(fv, mesh);
        this->SRp = presr;
    } else
        initMesh(fv, mesh);

    if (predec)
        preDecMatrix(fv, mesh);

    double preRes = 1e10;
    double res = globalOptimize(fv, mesh);
    int iter = 0;
    writeIterMesh(mesh, "iter", iter);

    while (abs(preRes - res) > iterEps && iter < maxIterTime) {
        iter ++;
        cout << "Iter time : " << iter << endl;
        preRes = res;
        localOptimize(fv, mesh);
        res = globalOptimize(fv, mesh);
        writeIterMesh(mesh, "iter", iter);
        cout << "DResidual : " << abs(preRes - res) << endl;
    }

    ckAlign(mesh);
}

void ARAPDeform::bfsInit(int i, std::vector<int> &visit, FeatureVector &fv, DTriMesh &mesh) {
    queue<int> q;
    q.push(i);
    visit[i] = true;
    while (!q.empty()) {
        i = q.front(); q.pop();
        Matrix3d &r = fv.r[i];
        DTriMesh::VertexHandle vi(i);
        int vs = ref.d[i];
        for (DTriMesh::VertexVertexIter vvi = mesh.vv_iter(vi); vvi; vvi++) {
            int j = vvi.handle().idx();
            if (!visit[j]) {
                fv.r[j] = r * fv.dr[vs];
                q.push(j);
                visit[j] = true;
            }
            vs++;
        }
    }
}

void ARAPDeform::bfsInit2(FeatureVector &fv, DTriMesh &mesh) {
    for (int i=0; i<ref.bfsq.size(); i++) {
        int j = ref.bfsq[i].first;
        int fa = ref.bfsq[i].second.first;
        int ei = ref.bfsq[i].second.second;
        if (fa != -1)
            fv.r[j] = fv.r[fa] * fv.dr[ei];
    }
}

void ARAPDeform::dfsInit(int i, std::vector<int> &visit, FeatureVector &fv, DTriMesh &mesh) {
    DTriMesh::VertexHandle vi(i);
    visit[i] = 1;
    Matrix3d &r = fv.r[i];
    int vs = ref.d[i];
    for (DTriMesh::VertexVertexIter vvi = mesh.vv_iter(vi); vvi; vvi++) {
        int j = vvi.handle().idx();
        if (!visit[j]) {
            fv.r[j] = r * fv.dr[vs];
            dfsInit(j, visit, fv, mesh);
        }
        vs++;
    }
}


void ARAPDeform::preDecMatrix(FeatureVector &fv, DTriMesh &mesh) {
    if (!AtAChanged) return;
    AtAChanged = false;
    lqsolver.saveA = true;
    globalOptimize(fv, mesh);
    lqsolver.saveA = false;
}

double ARAPDeform::globalOptimize(FeatureVector &fv, DTriMesh &mesh) {
    if (!lqsolver.saveA) return globalOptimizeFast(fv, mesh);

    cout << "Global Optimize ...... " << endl;
    long long tt = clock();
    vector< pair< pair<int,int>, double > > data;
    vector<double> vb;
    int n=0;
    for (int i=0; i<mesh.n_vertices(); i++) {
        DTriMesh::VertexHandle vi(i);
        Eigen::Matrix3d &ri = fv.r[i];
        int ei = ref.d[i];
        if (fv.isConst[i]) {
            for (int dim=0; dim<3; dim++) {
                data.push_back( make_pair( make_pair(n++, i*3+dim), 1 ));
                vb.push_back(0);
            }
        }
        for (DTriMesh::VertexVertexIter vvi = mesh.vv_iter(vi); vvi; vvi++, ei++) {
            int j = vvi.handle().idx();
            DTriMesh::VertexHandle vj(j);
            //double wij = ref.w[ei];
            double wij = ref.getWij(ei);
            Eigen::Matrix3d &rij = fv.dr[ei];
            Eigen::Matrix3d &sj = fv.s[j];
            DTriMesh::Point qj = ref.mesh->point(vj);
            int ej = ref.d[j];
            for (DTriMesh::VertexVertexIter vvj = mesh.vv_iter(vj); vvj; vvj++, ej++) {
                int k = vvj.handle().idx();
                double wjk = ref.w[ej];
                DTriMesh::Point tqjk = ref.mesh->point(DTriMesh::VertexHandle(k)) - qj;
                Eigen::Vector3d qjk(tqjk[0],tqjk[1],tqjk[2]);
                Eigen::Vector3d c = ri*(rij*(sj*qjk));
                if (fv.isConst[k]) c -= fv.constPoint[k];
                if (fv.isConst[j]) c += fv.constPoint[j];
                for (int dim=0; dim<3; dim++) {
                    if (!fv.isConst[k])
                        data.push_back(make_pair( make_pair(n, k*3+dim), wij * wjk ));
                    if (!fv.isConst[j])
                        data.push_back(make_pair( make_pair(n, j*3+dim), - wij * wjk ));
                    vb.push_back( wij * wjk * c(dim) );
                    n++;
                }
            }
        }
    }
    vector<double> result;
    cout << "matlab solve" << endl;
    lqsolver.needrs = true;
    double rs = lqsolver.solve(n, fv.s.size()*3, data, vb, result);
    if (lqsolver.saveA) return 0;
    lqsolver.needrs = false;
    cout << "copy ans" << endl;
    for (int i=0; i<mesh.n_vertices(); i++) {
        DTriMesh::VertexHandle vi(i);

        if (fv.isConst[i]) {
            mesh.point(vi) = EtoO(fv.constPoint[i]);
            continue;
        }

        for (int dim=0; dim<3; dim++)
            mesh.point(vi)[dim] = result[i*3+dim];
    }

    cout << "!!!Global Optimize Done , Time : " << (clock()-tt)*1.0 / CLOCKS_PER_SEC << " Rs : " << rs*rs << endl;
    return rs*rs;
}


double ARAPDeform::globalOptimizeFast(FeatureVector &fv, DTriMesh &mesh) {
    cout << "Global Optimize Fast ...... " << endl;
    long long tt = clock();
    vector<Vector3d> cv(fv.s.size()), rcv(fv.logdr.size(), Vector3d::Zero());
    vector<double> cs(fv.s.size());
    vector<Matrix3d> cm(fv.s.size());
    auto &rd = ref.rd;

#ifdef DUSE_OPENMP
#pragma omp parallel
{
#pragma omp for
#endif

    for (int j=0; j<mesh.n_vertices(); j++) {
        Eigen::Matrix3d &cmj = cm[j];
        double &csj = cs[j];

        cv[j] = Vector3d::Zero();
        cm[j] = Matrix3d::Zero();
        cs[j] = 0;

        for (int k=0; k<rd[j].size(); k++) {
            int i = rd[j][k].first;
            int ei = rd[j][k].second;

            double wij = ref.getWij(ei);
            Eigen::Matrix3d &rij = fv.dr[ei];
            Eigen::Matrix3d &sj = fv.s[j];

            cmj += (wij*wij) * fv.r[i] * rij * sj;
            csj += (wij*wij);
        }
    }
/*
    for (int i=0; i<mesh.n_vertices(); i++) {
        DTriMesh::VertexHandle vi(i);
        Eigen::Matrix3d &ri = fv.r[i];
        int ei = ref.d[i];
        for (DTriMesh::VertexVertexIter vvi = mesh.vv_iter(vi); vvi; vvi++, ei++) {
            int j = vvi.handle().idx();
            double wij = ref.getWij(ei);
            Eigen::Matrix3d &rij = fv.dr[ei];
            Eigen::Matrix3d &sj = fv.s[j];
            cm[j] += (wij*wij) * ri * rij * sj;
            cs[j] += (wij*wij);
        }
    }
*/
DOMP_END;

#ifdef DUSE_OPENMP
#pragma omp parallel
{
#pragma omp for
#endif

    for (int i=0; i<mesh.n_vertices(); i++) {
        DTriMesh::VertexHandle vi(i);
        int ei = ref.d[i];
        DTriMesh::Point qi = ref.mesh->point(vi);
        for (DTriMesh::VertexVertexIter vvi = mesh.vv_iter(vi); vvi; vvi++, ei++) {
            int j = vvi.handle().idx();
            DTriMesh::VertexHandle vj(j);
            Vector3d qij = OtoE(ref.mesh->point(vj) - qi);
            double wij = ref.w[ei];
            Eigen::Vector3d c = cm[i] * (qij * (wij*wij));
            if (fv.isConst[i]) c += fv.constPoint[i] * (cs[i] * (wij*wij));
            if (fv.isConst[j]) c -= fv.constPoint[j] * (cs[i] * (wij*wij));
            if (!fv.isConst[i]) cv[i] -= c;
            //if (!fv.isConst[j]) cv[j] += c;
            if (!fv.isConst[j]) rcv[ei] += c;
        }
    }

DOMP_END;

    vector<double> dataB(cv.size()*3), result;

#ifdef DUSE_OPENMP
#pragma omp parallel
{
#pragma omp for
#endif

    for (int j=0; j<mesh.n_vertices(); j++) {
        for (int k=0; k<rd[j].size(); k++) {
            int ei = rd[j][k].second;
            cv[j] += rcv[ei];
        }
        dataB[j*3] = cv[j](0);
        dataB[j*3+1] = cv[j](1);
        dataB[j*3+2] = cv[j](2);
    }

DOMP_END;


    lqsolver.solve(dataB.size(), dataB, result);

    for (int i=0; i<mesh.n_vertices(); i++) {
        DTriMesh::VertexHandle vi(i);

        if (fv.isConst[i]) {
            mesh.point(vi) = EtoO(fv.constPoint[i]);
            continue;
        }

        mesh.point(vi)[0] = result[i*3+0];
        mesh.point(vi)[1] = result[i*3+1];
        mesh.point(vi)[2] = result[i*3+2];
    }

    double residual = 0;
    if (needCalcRs) {
        vector<double> rss(omp_get_num_procs(), 0);

#ifdef DUSE_OPENMP
#pragma omp parallel
{
#pragma omp for
#endif

        for (int i=0; i<mesh.n_vertices(); i++) {
            double &residual = rss[omp_get_thread_num()];
            DTriMesh::VertexHandle vi(i);
            int ei = ref.d[i];
            DTriMesh::Point qi = ref.mesh->point(vi);
            DTriMesh::Point pi = mesh.point(vi);
            for (DTriMesh::VertexVertexIter vvi = mesh.vv_iter(vi); vvi; vvi++, ei++) {
                int j = vvi.handle().idx();
                DTriMesh::VertexHandle vj(j);
                Vector3d qij = OtoE(ref.mesh->point(vj) - qi);
                Vector3d pij = OtoE(mesh.point(vj) - pi);
                double wij = ref.w[ei];
                assert(abs(cs[i]) != 0);
                double ds = 0;
                ds += pij.squaredNorm() * cs[i];
                ds -= (cm[i] * qij).dot(pij) * 2;
                ds += (fv.s[i] * qij).squaredNorm() * cs[i];
                residual += ds * wij * wij;
            }
        }
DOMP_END;

        for (int i=0; i<rss.size(); i++) residual += rss[i];

    }

    cout << "Global Optimize Fast Done , Time : " << (clock()-tt)*1.0 / CLOCKS_PER_SEC << " Rs : " << residual << endl;
    return residual;
}

double ARAPDeform::localOptimize(FeatureVector &fv, DTriMesh &mesh) {
    localOptimizeFast(fv, mesh);
    return 0;

    cout << "Local Optimize ...... " << endl;
    long long tt = clock();
    double residual = 0;
    for (int i=0; i<mesh.n_vertices(); i++) {
        DTriMesh::VertexHandle vi(i);
        Eigen::Matrix3d &ri = fv.r[i];
        int ei = ref.d[i];
        vector<Eigen::Vector3d> vq, vp;
        for (DTriMesh::VertexVertexIter vvi = mesh.vv_iter(vi); vvi; vvi++, ei++) {
            int j = vvi.handle().idx();
            DTriMesh::VertexHandle vj(j);
            //double wij = ref.w[ei];
            double wij = ref.getWij(ei);
            Eigen::Matrix3d &rij = fv.dr[ei];
            Eigen::Matrix3d &sj = fv.s[j];
            DTriMesh::Point qj = ref.mesh->point(vj);
            DTriMesh::Point pj = mesh.point(vj);
            int ej = ref.d[j];
            for (DTriMesh::VertexVertexIter vvj = mesh.vv_iter(vj); vvj; vvj++, ej++) {
                int k = vvj.handle().idx();
                double wjk = ref.w[ej];
                DTriMesh::Point tqjk = ref.mesh->point(DTriMesh::VertexHandle(k)) - qj;
                DTriMesh::Point tpjk = mesh.point(DTriMesh::VertexHandle(k)) - pj;
                Eigen::Vector3d qjk(tqjk[0],tqjk[1],tqjk[2]);
                Eigen::Vector3d pjk(tpjk[0],tpjk[1],tpjk[2]);
                qjk = (rij*(sj*qjk));
                vq.push_back(qjk * (wij*wjk));
                vp.push_back(pjk * (wij*wjk));
            }
        }
        RotateAlign align(vp, vq);
        ri = align.calc();
        residual += align.res;
    }
    cout << "!!!Local Optimize Done , residual : " << residual << " Time : " << (clock()-tt)*1.0 / CLOCKS_PER_SEC << endl;
    return residual;
}

void ARAPDeform::localOptimizeFast(FeatureVector &fv, DTriMesh &mesh) {
    cout << "Local Optimize Fast...... " << endl;
    long long tt = clock();
    double residual = 0;
    vector<Eigen::Matrix3d> mats(fv.s.size());

#ifdef DUSE_OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
    for (int i=0; i<mesh.n_vertices(); i++) {
        DTriMesh::VertexHandle vi(i);
        int ei = ref.d[i];
        Eigen::Matrix3d &mat = mats[i];
        mat = Matrix3d::Zero();
        DTriMesh::Point qi = ref.mesh->point(vi);
        DTriMesh::Point pi = mesh.point(vi);
        for (DTriMesh::VertexVertexIter vvi = mesh.vv_iter(vi); vvi; vvi++, ei++) {
            int j = vvi.handle().idx();
            DTriMesh::VertexHandle vj(j);
            double wij = ref.w[ei];
            Vector3d tqij = OtoE(ref.mesh->point(vj) - qi);
            Vector3d tpij = OtoE(mesh.point(vj) - pi);
            mat += ((wij * wij) * tqij) * tpij.transpose();
        }
    }
DOMP_END;

#ifdef DUSE_OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
    for (int i=0; i<mesh.n_vertices(); i++) {
        DTriMesh::VertexHandle vi(i);
        Eigen::Matrix3d &ri = fv.r[i];
        int ei = ref.d[i];
        Matrix3d mat(Matrix3d::Zero());
        for (DTriMesh::VertexVertexIter vvi = mesh.vv_iter(vi); vvi; vvi++, ei++) {
            int j = vvi.handle().idx();
            double wij = ref.getWij(ei);
            Eigen::Matrix3d &rij = fv.dr[ei];
            Eigen::Matrix3d &sj = fv.s[j];
            mat += rij * sj * mats[j] * (wij * wij);
        }
        polarDec(mat, ri);
        ri.transposeInPlace();
    }
DOMP_END;

    cout << "Local Optimize Fast Done , " << " Time : " << (clock()-tt)*1.0 / CLOCKS_PER_SEC << endl;
}


double ARAPDeform::getWeightResidual(FeatureVector &fv, const std::vector<double> &weight, DTriMesh &mesh) {
    cout << "GetWR : ";
    for (int i=0; i<weight.size(); i++)
        cout << weight[i] << " ";
    cout << endl;
    fv.blendFrom(weight, this->fvs);
    initMesh(fv, mesh);
    return globalOptimize(fv, mesh);
}


std::vector<double> ARAPDeform::getGradient(FeatureVector &fv, const std::vector<double> &weight, DTriMesh &mesh, double &zeroResult) {
    vector<double> g, w=weight;
    double rs = getWeightResidual(fv, weight, mesh);
    zeroResult = rs;
    double eps = 1e-6;
    for (int i=1; i<weight.size(); i++) {
        w[0] -= eps;
        w[i] += eps;

        double rs2 = getWeightResidual(fv, w, mesh);
        g.push_back(- (rs2 - rs) / eps);

        w[0] += eps;
        w[i] -= eps;
    }
    return g;
}

double unifromGradient(std::vector<double> &weight) {
    double ans=0,s=0;
    for (int i=0; i<weight.size(); i++) {
        ans += weight[i] * weight[i];
        s += weight[i];
    }
    weight.insert(weight.begin(), -s);
    s = 1/sqrt(s*s + ans);
    for (int i=0; i<weight.size(); i++) weight[i] *= s;
    return sqrt(ans);
}

void ARAPDeform::solve(std::vector<double> weight, DTriMesh &mesh) {
    cout << "blend weight" <<endl;
    FeatureVector fv(weight, fvs);
    cout << "solve fv" << endl;
    solve(fv, mesh);
}

void ARAPDeform::ckAlign(DTriMesh &mesh) {
    if (needAlign)
        RotateAlign::AlignAtoB(mesh, *(ref.mesh));
}

void ARAPDeform::solve(FeatureVector fv, std::vector<double> &weight, DTriMesh &mesh) {
    initWeight(weight);
    fv.blendFrom(weight, fvs);
    initMesh(fv, mesh);

    globalOptimize(fv, mesh);

    double preRes = 1e10;
    double res = 1e9;
    int iter = 0;
    writeIterMesh(mesh, "iter", iter);

    while (abs(preRes - res) > iterEps && iter < maxIterTime) {
        iter ++;
        cout << "Iter time : " << iter << endl;
        preRes = res;

        solve(mesh, weight);
        fv.blendFrom(weight, fvs);

        res = localOptimize(fv, mesh);

        globalOptimize(fv, mesh);

        writeIterMesh(mesh, "iter", iter);
        cout << "DResidual : " << abs(preRes - res) << endl;
    }

    ckAlign(mesh);
}

std::vector<double> muladd(double x, std::vector<double> a, std::vector<double> b) {
    vector<double> ans(a.size());
    for (int i=0; i<a.size(); i++)
        ans[i] = x*a[i] + b[i];
    return ans;
}

void ARAPDeform::solve2(FeatureVector fv, std::vector<double> &weight, DTriMesh &mesh) {
    preDecMatrix(fv, mesh);
    int t=0;
    double teps = this->newtonIterEps;
    double best = sqrt(weight.size()*1.), num;
    TriNewtonIterSolver cvs;
    std::vector<double> g;
    cvs.eps = teps;
    cvs.convexFunction = [&](double x) -> double {
        return getWeightResidual(fv, muladd(x,g,weight), mesh);
    };
    if (weight.size() > 1)
    while (t++ < newtonIterTime) {
        double l=-abs(best)/sqrt(weight.size()*1.)/2, r=abs(best);
        g = getGradient(fv, weight, mesh, cvs.zeroResult);
        double wr = unifromGradient(g);
        cout << "Gradient Norm : " << wr << endl;
        cout << "Uni-Gradient Vector : ";
        for (int i=0; i<g.size(); i++)
            cout << g[i] << " ";
        cout << endl;
        best = 0;
        cvs.solve(l, r, best, num);
        cout << "Best : " << best << endl;
        weight = muladd(best,g,weight);
        cout << "Weight : ";
        for (int i=0; i<weight.size(); i++)
            cout << weight[i] << " ";
        cout << endl;
        writeIterMesh(mesh, "sf", t);
        if (abs(best) < teps) break;
    }
    fv.blendFrom(weight, fvs);
    solve(fv, mesh, false);
}


void ARAPDeform::solve(FeatureVector &fv, std::vector<double> &weight) {
    initWeight(weight);

    vector< pair< pair<int,int>, double > > data;
    vector<double> vb;
    int n=0;

    for (int i=0; i<fv.s.size(); i++) {
        for (int x=0; x<3; x++) for (int y=0; y<3; y++) {
            vb.push_back(fv.s[i](x,y));
            for (int j=0; j<fvs.size(); j++)
                data.push_back( make_pair( make_pair(n, j), fvs[j].s[i](x,y) ) );
            n++;
        }
    }

    for (int i=0; i<fv.logdr.size(); i++) {
        for (int x=0; x<3; x++) for (int y=0; y<3; y++) {
            vb.push_back(fv.logdr[i](x,y));
            for (int j=0; j<fvs.size(); j++)
                data.push_back( make_pair( make_pair(n, j), fvs[j].logdr[i](x,y) ) );
            n++;
        }
    }

    lqsolver.solve(n, weight.size(), data, vb, weight, true);

    cout << "Weights : ";
    double s = 0;
    for (int i=0; i<weight.size(); i++) {
        s += weight[i];
        cout << weight[i] << " ";
    }
    cout << endl << "S : " << s << endl;
}


void ARAPDeform::solve(DTriMesh &mesh, std::vector<double> &weight) {
    initWeight(weight);
    FeatureVector fv;
    ref.getFeature(mesh, fv);
    solve(fv, weight);
}

void ARAPDeform::initWeight(std::vector<double> &weight) {
    if (weight.size() != fvs.size()) {
        weight.resize(fvs.size(), 0);
        weight[0] = 1;
    }
}

double T_ARAPDeform::getWeightResidual(FeatureVector &fv, const std::vector<double> &weight, DTriMesh &mesh) {
    double prers = ARAPDeform::getWeightResidual(fv, weight, mesh);
    Vector3d s1(Vector3d::Zero()), s2(s1);
    int s=0;
    for (int i=0; i<ev.size(); i++)
        if (ev[i]) {
            s++;
            s1 += fv.constPoint[i];
            s2 += OtoE(mesh.point(DTriMesh::VertexHandle(i)));
        }
    s1 = (s1-s2) / s;

    double rs = 0;
    for (int i=0; i<ev.size(); i++)
        if (ev[i]) {
            cout << "fvCP : " << fv.constPoint[i].transpose();
            cout << " meshP : " << OtoE(mesh.point(DTriMesh::VertexHandle(i))).transpose();
            cout << " ds : " <<(fv.constPoint[i] - OtoE(mesh.point(DTriMesh::VertexHandle(i))) - s1).transpose();
            cout << endl;
            rs += (fv.constPoint[i] - OtoE(mesh.point(DTriMesh::VertexHandle(i))) - s1).squaredNorm();
        }
    cout << "VpRS : " << rs << endl;
    double ws = 0;
    for (int i=0; i<weight.size(); i++) {
        ws += abs(weight[i]);
    }
    cout << "ws : " << ws << endl;
    rs = rs * ws + std::max(0.0, rs - 1.2);
    cout << "WSRS : " << rs << endl;
    return rs;
}

void T_ARAPDeform::solve2(FeatureVector fv, std::vector<double> &weight, DTriMesh &mesh) {
    ev.resize(fv.isConst.size(), false);
    fv.isConst.swap(ev);

    org->AtAChanged = this->AtAChanged;
    org->maxIterTime = this->maxIterTime;
    org->newtonIterEps = this->newtonIterEps;
    org->newtonIterTime = this->newtonIterTime;
    org->iterEps = this->iterEps;

    cout << "TSolve2" << endl;

    preDecMatrix(fv, *(this->ref.mesh));
    int t=0;
    double teps = this->newtonIterEps;
    double best = sqrt(sqrt(weight.size()*1.)), num;
    TriNewtonIterSolver cvs;
    std::vector<double> g;
    cvs.eps = teps;
    cvs.convexFunction = [&](double x) -> double {
        return getWeightResidual(fv, muladd(x,g,weight), mesh);
    };
    while (t++ < newtonIterTime) {
        double l=-abs(best)/sqrt(weight.size()*1.)/2, r=abs(best);
        g = getGradient(fv, weight, mesh, cvs.zeroResult);
        double wr = unifromGradient(g);
        cout << "Gradient Norm : " << wr << endl;
        cout << "Uni-Gradient Vector : ";
        for (int i=0; i<g.size(); i++)
            cout << g[i] << " ";
        cout << endl;
        best = 0;
        cvs.solve(l, r, best, num);
        cout << "Best : " << best << endl;
        weight = muladd(best,g,weight);
        cout << "Weight : ";
        for (int i=0; i<weight.size(); i++)
            cout << weight[i] << " ";
        cout << endl;
        writeIterMesh(mesh, "sf", t);
        if (abs(best) < teps) break;
    }
    fv.isConst.swap(ev);
    org->solve(fv, mesh);
}

T_ARAPDeform::
T_ARAPDeform(DMEngine &eng, DTriMesh &refMesh, std::vector<DTriMesh*> ms)
    : ARAPDeform(eng, refMesh, ms) {
    org = new ARAPDeform(eng2, refMesh, ms);
}

// T2 code

void T2_ARAPDeform::solve2(FeatureVector fv, std::vector<double> &weight, DTriMesh &mesh) {
    ev.resize(fv.isConst.size(), false);
    fv.isConst.swap(ev);

    org->AtAChanged = this->AtAChanged;
    org->maxIterTime = this->maxIterTime;
    org->newtonIterEps = this->newtonIterEps;
    org->newtonIterTime = this->newtonIterTime;
    org->iterEps = this->iterEps;
    this->needCalcRs = false;

    cout << "T2Solve2" << endl;

    preDecMatrix(fv, *(this->ref.mesh));
    double teps = this->newtonIterEps;

    LQSolverWithLimit::newtonIterEps = this->newtonIterEps;
    LQSolverWithLimit::newtonMaxIterTime = this->newtonIterTime;
    //LQSolverWithLimit::iterEps = this->iterEps;

    vector<int> constList;
    vector<Vector3d> constPList;
    for (int i=0; i<ev.size(); i++)
        if (ev[i]) {
            constList.push_back(i);
            constPList.push_back(fv.constPoint[i]);
        }

    auto convert = [&](vector<Vector3d> &v) -> VectorXd {
        VectorXd res(v.size() * 3);
        Vector3d s(Vector3d::Zero());
        for (int i=0; i<v.size(); i++) s+=v[i];
        s /= v.size();
        for (int i=0; i<v.size(); i++) {
            auto x = v[i] - s;
            res(i*3) = x(0);
            res(i*3+1) = x(1);
            res(i*3+2) = x(2);
        };
        return res;
    };

    VectorXd targetB = convert(constPList);

    auto getX = [&](void) -> Eigen::VectorXd {
        ARAPDeform::getWeightResidual(fv, weight, mesh);
        Eigen::VectorXd res(constList.size() * 3);
        vector<Vector3d> v(constPList.size());
        for (int i=0; i<constList.size(); i++)
            v[i] = OtoE(mesh.point(DTriMesh::VertexHandle(constList[i])));
        return convert(v);
    };

    LQSolverWithLimit::GetJacobiFunction setJacobi =
            [&](
            std::vector<Eigen::VectorXd> &A,
            Eigen::VectorXd &b,
            std::vector< std::pair<double, double> > &limit
            ) {
        limit.resize(fvs.size()-1);
        for (int i=1; i<fvs.size(); i++)
            limit[i-1] = make_pair(rlimit.first - weight[i], rlimit.second - weight[i]);
        b = getX();
        A.resize(fvs.size()-1);
        double eps2 = 1e-6;
        for (int i=1; i<fvs.size(); i++) {
            weight[0] -= eps2;
            weight[i] += eps2;

            A[i-1] = (getX() - b) / eps2;

            weight[0] += eps2;
            weight[i] -= eps2;
        }
        b = targetB - b;
    };

    LQSolverWithLimit::ReturnRsFunction getRs =
            [&](std::vector<double> &x) {
        cout << "Get Delta : ";
        for (int i=0; i<x.size(); i++) {
            weight[i+1] += x[i];
            weight[0] -=  x[i];
            cout << weight[i+1] << "(" << x[i] << ") ";
        }
        cout << endl;
    };

    if (weight.size() > 1)
    LQSolverWithLimit::NewtonIter(setJacobi, getRs);

    fv.isConst.swap(ev);
    org->solve(fv, mesh);
}

T2_ARAPDeform::
T2_ARAPDeform(DMEngine &eng, DTriMesh &refMesh, std::vector<DTriMesh*> ms)
    : ARAPDeform(eng, refMesh, ms) {
    org = new ARAPDeform(eng2, refMesh, ms);
}

void ARAPTest() {
    DMEngine eng;
    DTriMesh mesh1, mesh2, mesh3, result;

    OpenMesh::IO::read_mesh(mesh1, "E:/project/SIGA2014/dataset/scape_smooth/49.obj");
    OpenMesh::IO::read_mesh(mesh2, "E:/project/SIGA2014/dataset/scape_smooth/50.obj");
    OpenMesh::IO::read_mesh(mesh3, "tmp.obj");

    vector<DTriMesh*> meshs;
    meshs.push_back(&mesh1);
    meshs.push_back(&mesh2);

    T_ARAPDeform deform(eng, mesh1, meshs);

    FeatureVector fv = deform.fvs[0];
    fv.loadConstPoint(ifstream("E:/project/SIGA2014/dataset/scape_smooth/handles/c1.txt"), mesh3);

    vector<double> weight;
    weight.resize(meshs.size(),0);
    weight[0] = 1;

    result = mesh1;

    deform.solve2(fv, weight, result);

    OpenMesh::IO::write_mesh(result, "result.obj");
    system("pause");
}
