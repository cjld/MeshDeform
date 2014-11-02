#include "deformsf.h"

struct EdgeInfo {
	double length;
	double angle;
	DTriMesh::VertexHandle vh[4];
	// derivative of angle
	Eigen::Vector3d d_angle[4];
	// derivative of length
	Eigen::Vector3d d_len;
	// sum of 2 sides of edge
	double area;
};


EdgeInfo getEdgeInfo(DTriMesh &mesh, DTriMesh::HalfedgeHandle it) {
	EdgeInfo info;
	info.vh[0] = mesh.to_vertex_handle(it);
	info.vh[1] = mesh.from_vertex_handle(it);
	info.vh[2] = mesh.to_vertex_handle(
	                 mesh.next_halfedge_handle(it));
	info.vh[3] = mesh.to_vertex_handle(
	                 mesh.next_halfedge_handle(
                          mesh.opposite_halfedge_handle(it)));

	Eigen::Vector3d x1 = OtoE(mesh.point(info.vh[0]));
	Eigen::Vector3d x2 = OtoE(mesh.point(info.vh[1]));
	Eigen::Vector3d x3 = OtoE(mesh.point(info.vh[2]));
	Eigen::Vector3d x4 = OtoE(mesh.point(info.vh[3]));

	Eigen::Vector3d e = x2 - x1;
	Eigen::Vector3d n1 = (x3 - x2).cross(x3 - x1);
	Eigen::Vector3d n2 = (x4 - x1).cross(x4 - x2);

	double n1n = n1.norm();
	double n2n = n2.norm();
	info.area = n1n + n2n;

	info.angle = acos(n1.dot(n2) / (n1n * n2n));
	info.length = e.norm();

	info.d_len = - e / e.norm();
	// TODO: for all 4 vertices
	info.d_angle[0] = (x3 - x2).dot(e) * n1 / (e.norm() * n1n * n1n) +
	                  (x4 - x2).dot(e) * n2 / (e.norm() * n2n * n2n);
	//info.d_angle_x34 = e.length() /

	return info;
}

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
}

Eigen::VectorXd DeformSF::getFeature(DTriMesh &mesh) {
	Eigen::VectorXd fv(featureSize);
	int i = 0;
	DTriMesh::EdgeIter it = mesh.edges_begin();
	while (it != mesh.edges_end()) {
		EdgeInfo info = getEdgeInfo(mesh, mesh.halfedge_handle(*it, 0));
		fv(i) = info.length;
		fv(i+1) = info.angle;

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
	fv(i) = v / 6;
	return fv;
}

void DeformSF::setTerm(double stretchTerm, double bendTerm, double volumeTerm)
{
	this->stretchTerm = stretchTerm;
	this->bendTerm = bendTerm;
	this->volumeTerm = volumeTerm;
}

Eigen::VectorXd DeformSF::getWeight(DTriMesh &mesh)
{
	Eigen::VectorXd w(featureSize);

	int i = 0;
	DTriMesh::EdgeIter it = mesh.edges_begin();
	while (it != mesh.edges_end()) {
		EdgeInfo info = getEdgeInfo(mesh, mesh.halfedge_handle(*it, 0));
		w(i) = sqrt(this->stretchTerm) / info.length;
		w(i+1) = sqrt(this->bendTerm / info.area) * info.length;

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
	w(i) = sqrt(this->volumeTerm) / v;
	return w;
}

DMSpMatrix DeformSF::getJacob(DTriMesh &mesh) {
	std::vector< std::pair< std::pair<int,int>, double> > data;
	DMSpMatrix j(*eng, "J", xSize, featureSize, data, true);

	int i = 0;
	DTriMesh::EdgeIter it = mesh.edges_begin();
	while (it != mesh.edges_end()) {
		EdgeInfo info = getEdgeInfo(mesh, mesh.halfedge_handle(*it, 0));
		data.push_back( make_pair( make_pair() ));
		it++, i+=2;
	}

	return j;
}
