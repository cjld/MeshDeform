#include "deformsf.h"

struct EdgeInfo {
	double length;
	double angle;
	// derivative of angle
	Eigen::Vector3d d_angle[4];
	// derivative of length
	Eigen::Vector3d d_len;
};


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

EdgeInfo getEdgeInfo(DTriMesh &mesh, DTriMesh::HalfedgeHandle it) {
	EdgeInfo info;

	Eigen::Vector3d x1 = OtoE(mesh.point(mesh.to_vertex_handle(it)));
	Eigen::Vector3d x2 = OtoE(mesh.point(mesh.from_vertex_handle(it)));
	Eigen::Vector3d x3 = OtoE(mesh.point(mesh.to_vertex_handle(
					         mesh.next_halfedge_handle(it))));
	Eigen::Vector3d x4 = OtoE(mesh.point(mesh.to_vertex_handle(
	                         mesh.next_halfedge_handle(
	                             mesh.opposite_halfedge_handle(it)))));

	Eigen::Vector3d e = x2 - x1;
	Eigen::Vector3d n1 = (x3 - x2).cross(x3 - x1);
	Eigen::Vector3d n2 = (x4 - x1).cross(x4 - x2);

	double n1n = n1.norm();
	double n2n = n2.norm();

	info.angle = acos(n1.dot(n2) / (n1n * n2n));
	info.length = e.norm();

	info.d_len = - e / e.norm();
	info.d_angle[0] = (x3 - x2).dot(e) * n1 / (e.norm() * n1n * n1n) +
	                  (x4 - x2).dot(e) * n2 / (e.norm() * n2n * n2n);
	//info.d_angle_x34 = e.length() /

	return info;
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
		DTriMesh::FaceVertexIter fvi = mesh.fv_iter(*fi);
		DTriMesh::VertexHandle v1 = *fvi;
		fvi++;
		DTriMesh::VertexHandle v2 = *fvi;
		fvi++;
		DTriMesh::VertexHandle v3 = *fvi;
		assert(v1 != v2 && v1 != v3 && v2 != v3);
		Eigen::Vector3d x1 = OtoE(mesh.point(v1));
		Eigen::Vector3d x2 = OtoE(mesh.point(v2));
		Eigen::Vector3d x3 = OtoE(mesh.point(v3));
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
