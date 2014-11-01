#ifndef DEFORMSF_H
#define DEFORMSF_H

#include "ARAPDeform.h"

/*
 * A Implementation of
 * Example-Driven Deformations Based on Discrete Shells
 * Volume xx (200y), Number z, pp. 1–11
 * Stefan Fröhlich, Mario Botsch
 * */
class DeformSF
{
public:
	DMEngine *eng;
	std::vector<DTriMesh*> meshs;
	std::vector<Eigen::VectorXd> fvs;

	// weight for feature
	Eigen::VectorXd w;
	int featureSize;
	int vertexSize;
	int edgeSize;

	DeformSF(DMEngine &eng, std::vector<DTriMesh*> ms);

	// set these term to blend energy
	double stretchTerm, bendTerm, volumeTerm;
	// iterator therhold
	double iterEps;


	void setTerm(
	        double stretchTerm = 100,
	        double bendTerm = 1,
	        double volumeTerm = 1000);


	// given const point, solve mesh
	void solve(FeatureVector fv, DTriMesh &mesh);

	// given feature vector weights, solve mesh
    void solve(std::vector<double> weight, DTriMesh &mesh);

	Eigen::VectorXd getFeature(DTriMesh &mesh);
};

#endif // DEFORMSF_H
