#ifndef DMENGINE_H
#define DMENGINE_H
#include "engine.h"
#include "mex.h"
#include <string>
#include <iostream>
#include <cstring>
#include <vector>
#include <algorithm>

class DMEngine {

    Engine *ep;

public:
    DMEngine() {
        int statusmat = 0;
        static bool firstOpen = false;
        if (firstOpen)
            ep = engOpen("");
        else
            ep = engOpenSingleUse(NULL, NULL,&statusmat);
        firstOpen = false;
        if (!ep) {
            std::cerr << "\nCan't start MATLAB engine\n" << std::endl;
        }
    }

    ~DMEngine() {
        engClose(ep);
    }

    int Eval(std::string s) {
        return engEvalString(ep, s.c_str());
    }

    int PutVariable(std::string var_name, const mxArray *ap) {
        return engPutVariable(ep, var_name.c_str(), ap);
    }

    mxArray* GetVariable(std::string var_name) {
        return engGetVariable(ep, var_name.c_str());
    }

    void Pause() {system("pause");}
};

class DMMatrix {

public:
    mxArray *p;
    DMEngine *eng;
    std::string name;
    int n,m;

    DMMatrix(DMEngine& eng, std::string name, int n, int m, double *data = NULL) {
        this->n = n, this->m = m;
        this->eng = &eng;
        this->name = name;
        if (!data) {
            Update();
            return;
        }
        p = mxCreateDoubleMatrix(n,m,mxREAL);
        memcpy((void *)(mxGetPr(p)), (void *)data, sizeof(double)*n*m);
        eng.PutVariable(name, p);
        p = eng.GetVariable(name);
    }

    double& GetNum(int i, int j) {return mxGetPr(p)[i+j*n];}

    ~DMMatrix() {
        mxDestroyArray(p);
    }

    friend std::ostream& operator<<(std::ostream& out, DMMatrix &mat) {
        out << mat.name << " : \n";
        for (int i=0; i<mat.n; i++) {
            for (int j=0; j<mat.m; j++)
                out << mat.GetNum(i,j) << '\t';
            out << '\n';
        }
        out << std::endl;
        return out;
    }

    void Update() {p = eng->GetVariable(name);}
};

typedef std::vector< std::pair< std::pair<int, int>, double> > DMSpMatrixData;

class DMSpMatrix {

public:
    mxArray *p;
    DMEngine *eng;
    std::string name;
    int n,m,l;

	static void maintain(DMSpMatrixData &data) {
		std::sort(data.begin(), data.end(),
			[](const std::pair< std::pair<int, int>, double>& a, const std::pair< std::pair<int, int>, double>& b) -> bool {
			return
				a.first.second < b.first.second ||
				a.first.second == b.first.second && a.first.first < b.first.first;
		});
		int j = 0;
		for (int i = 1; i < data.size(); i++)
			if (data[i].first != data[j].first) {
			j++;
			if (i != j) data[j] = data[i];
			}
			else
				data[j].second += data[i].second;
		data.erase(data.begin() + j, data.end());
	}

    DMSpMatrix(DMEngine& eng, std::string name, int n, int m, DMSpMatrixData data, bool NeedSort = false) {
        this->n = n, this->m = m;
        this->eng = &eng;
        this->name = name;
        this->l = (int)data.size();
        p = mxCreateSparse(n,m,data.size(),mxREAL);
        if (NeedSort)
			maintain(data);
        for (int i=0; i<(int)data.size(); i++) {
            mxGetIr(p)[i] = data[i].first.first;
            mxGetJc(p)[data[i].first.second+1] ++;
            mxGetPr(p)[i] = data[i].second;
        }
        for (int i=1; i<=m; i++) mxGetJc(p)[i] += mxGetJc(p)[i-1];
        eng.PutVariable(name, p);
        p = eng.GetVariable(name);
    }

    std::pair< std::pair<int,int>, double> GetNum(int i) {
        return std::make_pair( std::make_pair(mxGetIr(p)[i], mxGetJc(p)[i]), mxGetPr(p)[i]);
    }

    friend std::ostream& operator<<(std::ostream& out, DMSpMatrix &mat) {
        out << mat.name << " : \n";
        for (int i=0; i<mat.l; i++) {
            auto val = mat.GetNum(i);
            out << "{(" << val.first.first << ", " << val.first.second << ") = " << val.second << "}\n";
        }
        out << std::endl;
        return out;
    }

    ~DMSpMatrix() {
        mxDestroyArray(p);
    }

    void Update() {p = eng->GetVariable(name);}
};

#endif // DMENGINE_H
