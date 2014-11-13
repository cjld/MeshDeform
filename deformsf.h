#ifndef DEFORMSF_H
#define DEFORMSF_H

#include "ARAPDeform.h"

/**
 * A Implementation of
 * Example-Driven Deformations Based on Discrete Shells
 * Volume xx (200y), Number z, pp. 1–11
 * Stefan Fröhlich, Mario Botsch
 * */

class DeformSF
{
	typedef std::vector<double> Feature;
public:
	DMEngine *eng;
	std::vector<DTriMesh*> meshs;
	std::vector<Feature> fvs;

	// weight for feature
	Feature w;
	int featureSize;
	int vertexSize;
	int edgeSize;
	int xSize;

	DeformSF(DMEngine &eng, std::vector<DTriMesh*> ms);

	// set these term to blend energy
	double stretchTerm, bendTerm, volumeTerm;
	// iterator therhold
	double iterEps;

	// default value from the paper
	void setTerm(
	        double stretchTerm = 100,
	        double bendTerm = 1,
	        double volumeTerm = 1000);

	// given const point, solve mesh
	void solve(FeatureVector fv, DTriMesh &mesh);

	// given feature vector weights, solve mesh
    void solve(std::vector<double> weight, DTriMesh &mesh);

	// have not been squared
	Feature getFeature(DTriMesh &mesh);
	Feature getWeight(DTriMesh &mesh);
	DMSpMatrixData getJacob(FeatureVector &fv, DTriMesh &mesh);
	void static main();
};

#endif // DEFORMSF_H
