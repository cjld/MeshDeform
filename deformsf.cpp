#include "deformsf.h"

using namespace std;

#define sqr(x) ((x) * (x))

struct EdgeInfo {
	double length, n_len[2];
	double angle;

	/**
	 *    3
	 *   /  \
	 *  0 -> 1
	 *	 \  /
	 *    2
	 */
	DTriMesh::VertexHandle vh[4];
	Eigen::Vector3d x[4], e, n[2];

	EdgeInfo(DTriMesh &mesh, DTriMesh::HalfedgeHandle it) 
	{
		vh[0] = mesh.to_vertex_handle(it);
		vh[1] = mesh.from_vertex_handle(it);
		vh[2] = mesh.to_vertex_handle(
			mesh.next_halfedge_handle(it));
		vh[3] = mesh.to_vertex_handle(
			mesh.next_halfedge_handle(
			mesh.opposite_halfedge_handle(it)));

		for (int i = 0; i < 4; i++)
			x[i] = OtoE(mesh.point(vh[3]));

		e = x[1] - x[0];
		length = e.norm();

		n[0] = (x[2] - x[0]).cross(e);
		n[1] = e.cross(x[3] - x[0]);
		n_len[0] = n[0].norm();
		n_len[1] = n[1].norm();

		angle = acos(n[0].dot(n[1]) / (n_len[1] * n_len[0]));
		if (n[0].cross(n[1]).dot(e) < 0)
			angle = M_PI - angle;
		else
			angle = M_PI + angle;
	}

	// Derivative of angle (x0~3) and length(x0~1)
	void getDerivative(Eigen::Vector3d d_angle[4], Eigen::Vector3d d_len[2])
	{
		d_len[0] = - e / length;
		d_len[1] = - d_len[0];

		d_angle[0] =
			(x[2] - x[1]).dot(e) / (length * sqr(n_len[0])) * n[0] +
			(x[3] - x[1]).dot(e) / (length * sqr(n_len[1])) * n[1];
		d_angle[1] = -(
			(x[2] - x[0]).dot(e) / (length * sqr(n_len[0])) * n[0] +
			(x[3] - x[0]).dot(e) / (length * sqr(n_len[1])) * n[1]);
		d_angle[2] = n[0] * length / sqr(n_len[0]);
		d_angle[3] = n[1] * length / sqr(n_len[1]);
	}

	// return the sum of the areas of the two triangles sharing the edge
	double getArea() 
	{
		return n_len[0] + n_len[1];
	}
};

void getTriFaceVertex(DTriMesh &mesh, DTriMesh::FaceHandle fi,
          /* output */DTriMesh::VertexHandle vh[3])
{
	DTriMesh::FaceVertexIter fvi = mesh.fv_iter(fi);
	vh[0] = *fvi;
	fvi++;
	vh[1] = *fvi;
	fvi++;
	vh[2] = *fvi;
	assert(vh[0] != vh[1] && vh[1] != vh[2] && vh[0] != vh[2]);
}

DeformSF::DeformSF(DMEngine &eng, std::vector<DTriMesh*> ms)
{
	this->eng = &eng;
	this->setTerm();

	meshs = ms;
	iterEps = 1e-10;

	// calc feature vector
	for (int i=0; i<meshs.size(); i++)
		fvs.push_back(getFeature(*meshs[i]));

	this->edgeSize = ((DTriMesh *)meshs[0])->n_edges();
	this->vertexSize = ((DTriMesh *)meshs[0])->n_vertices();
	this->featureSize = edgeSize * 2 + 1;

	this->w = this->getWeight(*ms[0]);
	this->fvs.resize(ms.size());
	for (int i = 0; i < ms.size(); i++)
		fvs[i] = getFeature(*ms[i]);
}

DeformSF::Feature DeformSF::getFeature(DTriMesh &mesh) {
	Feature fv(featureSize);
	int i = 0;
	DTriMesh::EdgeIter it = mesh.edges_begin();
	while (it != mesh.edges_end()) {
		EdgeInfo info(mesh, mesh.halfedge_handle(*it, 0));
		fv[i] = info.length;
		fv[i+1] = info.angle;

		it++, i+=2;
	}

	double v=0;
	for (DTriMesh::FaceIter fi = mesh.faces_begin(); fi != mesh.faces_end(); fi++) {
		DTriMesh::VertexHandle vh[3];
		getTriFaceVertex(mesh, *fi, vh);
		Eigen::Vector3d x1 = OtoE(mesh.point(vh[0]));
		Eigen::Vector3d x2 = OtoE(mesh.point(vh[1]));
		Eigen::Vector3d x3 = OtoE(mesh.point(vh[2]));
		v += x1.cross(x2).dot(x3);
	}
	fv[i] = v / 6;
	return fv;
}

void DeformSF::setTerm(double stretchTerm, double bendTerm, double volumeTerm)
{
	this->stretchTerm = stretchTerm;
	this->bendTerm = bendTerm;
	this->volumeTerm = volumeTerm;
}

DeformSF::Feature DeformSF::getWeight(DTriMesh &mesh)
{
	Feature w(featureSize);

	int i = 0;
	DTriMesh::EdgeIter it = mesh.edges_begin();
	while (it != mesh.edges_end()) {
		EdgeInfo info(mesh, mesh.halfedge_handle(*it, 0));
		w[i] = sqrt(this->stretchTerm) / info.length;
		w[i+1] = sqrt(this->bendTerm / info.getArea()) * info.length;

		it++, i+=2;
	}

	double v=0;
	for (DTriMesh::FaceIter fi = mesh.faces_begin(); fi != mesh.faces_end(); fi++) {
		DTriMesh::VertexHandle vh[3];
		getTriFaceVertex(mesh, *fi, vh);
		//int a = vh[0].idx();
		Eigen::Vector3d x1 = OtoE(mesh.point(vh[0]));
		Eigen::Vector3d x2 = OtoE(mesh.point(vh[1]));
		Eigen::Vector3d x3 = OtoE(mesh.point(vh[2]));
		v += x1.cross(x2).dot(x3);
	}
	w[i] = sqrt(this->volumeTerm) / v;
	return w;
}

#define MP3(x,y,z) make_pair( make_pair((x),(y)), (z))
#define push_back_V3d(data, x, y, z) \
if (!fv.isConst[j]) { \
	data.push_back(MP3((x), (y)*3, z(0))); \
	data.push_back(MP3((x), (y)*3+1, z(1))); \
	data.push_back(MP3((x), (y)*3+2, z(2))); \
}

DMSpMatrixData DeformSF::getJacob(FeatureVector &fv, DTriMesh &mesh) {
	std::vector< std::pair< std::pair<int,int>, double> > data;

	int i = 0;
	DTriMesh::EdgeIter it = mesh.edges_begin();
	Eigen::Vector3d d_angle[4], d_len[2];
	while (it != mesh.edges_end()) {
		EdgeInfo info(mesh, mesh.halfedge_handle(*it, 0));
		info.getDerivative(d_angle, d_len);
		for (int j = 0; j < 4; j++)
			push_back_V3d(data, i, info.vh[j].idx(), d_angle[j]);
		for (int j = 0; j < 2; j++)
			push_back_V3d(data, i + 1, info.vh[j].idx(), d_len[j]);
		it++, i+=2;
	}

	for (DTriMesh::FaceIter fi = mesh.faces_begin(); fi != mesh.faces_end(); fi++) {
		DTriMesh::VertexHandle vh[3];
		getTriFaceVertex(mesh, *fi, vh);
		Eigen::Vector3d x1 = OtoE(mesh.point(vh[0]));
		Eigen::Vector3d x2 = OtoE(mesh.point(vh[1]));
		Eigen::Vector3d x3 = OtoE(mesh.point(vh[2]));
		Eigen::Vector3d n = (x2 - x1).cross(x3 - x1);

		for (int j = 0; j < 3; j++)
			push_back_V3d(data, i, vh[j].idx(), n);
	}

	return data;
}

void DeformSF::solve(FeatureVector fv, DTriMesh &mesh) {

	for (int i = 0; i < fv.isConst.size(); i++)
		if (fv.isConst[i])
			mesh.point(DTriMesh::VertexHandle(i)) = EtoO(fv.constPoint[i]);
	int it = 20;
	DMSpMatrixData data;
	Feature &target = fvs[0], now;
	DMMatrix f_(*eng, "f_", featureSize, 1, &target[0]);
	DMMatrix w(*eng, "w", featureSize, 1, &this->w[0]);
	DMMatrix dx(*eng, "dx", xSize, 1);
	while (it--) {
		now = getFeature(mesh);
		data = getJacob(fv, mesh);
		DMSpMatrix J(*eng, "J", featureSize, xSize, data, true);
		DMMatrix f(*eng, "f", featureSize, 1, &now[0]);
		eng->Eval("J = J.*w");
		eng->Eval("f = w.*(f-f_)");
		eng->Eval("dx = (J'*J)\\-J'*f");

		dx.Update();
		for (int i = 0; i < this->vertexSize; i++)
			if (!fv.isConst[i]) {
				for (int dim = 0; dim < 3; dim++)
					mesh.point(DTriMesh::VertexHandle(i))[dim] = dx.GetNum(i*3+dim, 1);
			}
	}
}

void DeformSF::main() {
	string ws = "F:/Desktop/Lab/deform_misc/MeshDeform/build/workspace";
}