#ifndef lightJetChiSquareMinimumSolver_h
#define lightJetChiSquareMinimumSolver_h

#include <vector>
#include <cmath>
#include <memory>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "TMath.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompBase.h"

#include "Math/GenVector/LorentzVector.h"

using namespace std;
using namespace ROOT::Math;
typedef LorentzVector<PxPyPzE4D<double>> XYZTLorentzVector;

/*
 * Chi-square minimizer for recoiling particles
 *
 */

struct recoil_minimizer_input {
    vector<XYZTLorentzVector> *ps; // recoiling particles
    vector<double> *jetPtWidths;
    vector<double> *jetPhiWidths;
    vector<double> *jetEtaWidths;

    double *dx, *dy, *dz;
};

struct recoil_minimizer_data {
    unsigned int n_ps;

    double dxCheck, dyCheck, dzCheck;

    vector<double> minDeltasX;
    vector<double> minDeltasY;
    vector<double> minDeltasZ;

    vector<TMatrixD> cov_rad;       // covariance matrix; radial (pT, phi, eta)
    vector<TMatrixD> cov;           // covariance matrix
    vector<TMatrixD> inv_sum_x_cov; // cov[i] * inv_sum_cov

    TMatrixD inv_sum_cov; // inverted sum of covariance matrices
};

class lightJetChiSquareMinimumSolver
{
  private:
    const bool do3D_;
    const bool debug = true;

    recoil_minimizer_input input_;
    recoil_minimizer_data data_;

    TDecompBase *inverter3D_; // matrix inverter

    double chi2_;

    void calcMin();
    void checkSize(recoil_minimizer_input &);
    void Eval_covariance(const recoil_minimizer_input &,
                         recoil_minimizer_data &);
    void Eval_cov_sum(recoil_minimizer_data &);

  public:
    lightJetChiSquareMinimumSolver(vector<XYZTLorentzVector> &,
                                   vector<double> &, vector<double> &,
                                   vector<double> &, double &, double &,
                                   double &, bool);
    lightJetChiSquareMinimumSolver(const unsigned int, double &, double &,
                                   double &, bool);
    lightJetChiSquareMinimumSolver(const lightJetChiSquareMinimumSolver &other);
    ~lightJetChiSquareMinimumSolver();

    inline void Init_data(const unsigned int &, recoil_minimizer_data &);

    void setupEquations(vector<XYZTLorentzVector> &, vector<double> &,
                        vector<double> &, vector<double> &);
    void printResults();

    int nJets() { return data_.n_ps; };

    // void setRecoil(double , double, double );

    vector<double> *getMinDeltasX() { return &data_.minDeltasX; };
    vector<double> *getMinDeltasY() { return &data_.minDeltasY; };
    vector<double> *getMinDeltasZ() { return &data_.minDeltasZ; };

    double getChiSquare();
};

#endif
