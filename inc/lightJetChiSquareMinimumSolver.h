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

class lightJetChiSquareMinimumSolver
{

  private:
    const bool do3D_;
    const bool debug = false;

    vector<double> jetPxWidths2_, jetPyWidths2_, jetPzWidths2_, jetPxPyWidths_,
        jetPxPzWidths_, jetPyPzWidths_;

    const double &dx_, &dy_, &dz_;
    double dxCheck_, dyCheck_, dzCheck_;

    vector<TMatrixD> cov2D_;           // covariance matrix for 2D case
    vector<TMatrixD> cov3D_;           // covariance matrix for 2D case
    vector<TMatrixD> inv_sum_x_cov2D_; // cov2D_[i] * inv_sum_cov2D_
    vector<TMatrixD> inv_sum_x_cov3D_; // cov3D_[i] * inv_sum_cov3D_

    TDecompBase *inverter2D_;
    TDecompBase *inverter3D_;

    TMatrixD inv_sum_cov2D_; // inverted sum of covariance matrices for 2D
    TMatrixD inv_sum_cov3D_; // inverted sum of covariance matrices for 3D

    vector<double> minDeltasX_;
    vector<double> minDeltasY_;
    vector<double> minDeltasZ_;

    double chi2_;

    unsigned int nJets_;

    void setCartesianWidths(vector<XYZTLorentzVector> &, vector<double> &,
                            vector<double> &, vector<double> &);

    void calcMin();
    inline void calcMin_2D();
    inline void calcMin_3D();

    void calcSigmas();
    inline void calcSigmas_2D();
    inline void calcSigmas_3D();

    void checkSize(vector<XYZTLorentzVector> &, vector<double> &,
                   vector<double> &, vector<double> &);

  public:
    lightJetChiSquareMinimumSolver(vector<XYZTLorentzVector> &,
                                   vector<double> &, vector<double> &,
                                   vector<double> &, double &, double &,
                                   double &, bool);
    lightJetChiSquareMinimumSolver(int, double &, double &, double &, bool);

    lightJetChiSquareMinimumSolver(const lightJetChiSquareMinimumSolver &other);

    void setupEquations(vector<XYZTLorentzVector> &, vector<double> &,
                        vector<double> &, vector<double> &);
    void printResults();

    int nJets() { return nJets_; };

    double jetPxWidth2(int i) { return jetPxWidths2_.at(i); };
    double jetPyWidth2(int i) { return jetPyWidths2_.at(i); };
    double jetPzWidth2(int i) { return jetPzWidths2_.at(i); };
    double jetPxPyWidth(int i) { return jetPxPyWidths_.at(i); };
    double jetPxPzWidth(int i) { return jetPxPzWidths_.at(i); };
    double jetPyPzWidth(int i) { return jetPyPzWidths_.at(i); };

    // void setRecoil(double , double, double );

    vector<double> *getMinDeltasX() { return &minDeltasX_; };
    vector<double> *getMinDeltasY() { return &minDeltasY_; };
    vector<double> *getMinDeltasZ() { return &minDeltasZ_; };

    ~lightJetChiSquareMinimumSolver();

    double getChiSquare();
};

#endif
