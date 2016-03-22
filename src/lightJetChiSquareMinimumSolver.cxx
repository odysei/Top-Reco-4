#ifndef lightJetChiSquareMinimumSolver_cxx
#define lightJetChiSquareMinimumSolver_cxx

#include "lightJetChiSquareMinimumSolver.h"
#include "TDecompSVD.h"
#include "TDecompLU.h"

lightJetChiSquareMinimumSolver::lightJetChiSquareMinimumSolver(
    vector<XYZTLorentzVector> &jets, vector<double> &jetPtWidths,
    vector<double> &jetPhiWidths, vector<double> &jetEtaWidths, double &dx,
    double &dy, double &dz, bool do3D)
    : do3D_(do3D), dx_(dx), dy_(dy), dz_(dz), inv_sum_cov2D_(2, 2),
      inv_sum_cov3D_(3, 3)
{
    const unsigned int n_jets = jets.size();
    nJets_ = n_jets;
    dxCheck_ = 0.;
    dyCheck_ = 0.;
    dzCheck_ = 0.;
    chi2_ = 0.;

    jetPxWidths2_.resize(n_jets, 0.);
    jetPyWidths2_.resize(n_jets, 0.);
    jetPzWidths2_.resize(n_jets, 0.);
    jetPxPyWidths_.resize(n_jets, 0.);
    jetPxPzWidths_.resize(n_jets, 0.);
    jetPyPzWidths_.resize(n_jets, 0.);
    minDeltasX_.resize(n_jets, 0.);
    minDeltasY_.resize(n_jets, 0.);
    minDeltasZ_.resize(n_jets, 0.);
    cov2D_.reserve(n_jets);
    cov3D_.reserve(n_jets);
    inv_sum_x_cov2D_.reserve(n_jets);
    inv_sum_x_cov3D_.reserve(n_jets);
    for (unsigned int i = 0; i < n_jets; ++i) {
        cov2D_.push_back(TMatrixD(2, 2));
        cov3D_.push_back(TMatrixD(3, 3));
        inv_sum_x_cov2D_.push_back(TMatrixD(2, 2));
        inv_sum_x_cov3D_.push_back(TMatrixD(3, 3));
    }

    inverter2D_ = new TDecompLU(2);
    inverter3D_ = new TDecompLU(3);

    checkSize(jets, jetPtWidths, jetPhiWidths, jetEtaWidths);
    setCartesianWidths(jets, jetPtWidths, jetPhiWidths, jetEtaWidths);
    calcSigmas();
}

lightJetChiSquareMinimumSolver::lightJetChiSquareMinimumSolver(
    int nObjects, double &dx, double &dy, double &dz, bool do3D)
    : do3D_(do3D), dx_(dx), dy_(dy), dz_(dz), inv_sum_cov2D_(2, 2),
      inv_sum_cov3D_(3, 3)
{
    const unsigned int n_jets = nObjects;
    nJets_ = n_jets;
    dxCheck_ = 0.;
    dyCheck_ = 0.;
    dzCheck_ = 0.;
    chi2_ = 0.;

    jetPxWidths2_.resize(n_jets, 0.);
    jetPyWidths2_.resize(n_jets, 0.);
    jetPzWidths2_.resize(n_jets, 0.);
    jetPxPyWidths_.resize(n_jets, 0.);
    jetPxPzWidths_.resize(n_jets, 0.);
    jetPyPzWidths_.resize(n_jets, 0.);
    minDeltasX_.resize(n_jets, 0.);
    minDeltasY_.resize(n_jets, 0.);
    minDeltasZ_.resize(n_jets, 0.);
    cov2D_.reserve(n_jets);
    cov3D_.reserve(n_jets);
    inv_sum_x_cov2D_.reserve(n_jets);
    inv_sum_x_cov3D_.reserve(n_jets);
    for (unsigned int i = 0; i < n_jets; ++i) {
        cov2D_.push_back(TMatrixD(2, 2));
        cov3D_.push_back(TMatrixD(3, 3));
        inv_sum_x_cov2D_.push_back(TMatrixD(2, 2));
        inv_sum_x_cov3D_.push_back(TMatrixD(3, 3));
    }

    inverter2D_ = new TDecompLU(2);
    inverter3D_ = new TDecompLU(3);
    // cout << "Light jet chi square constructor" << endl;
    // cout << "dx is " << dx_ << endl;
    // cout << "dy is " << dy_ << endl;
    // cout << "dz is " << dz_ << endl;
    // cout << "dxCheck is " << dxCheck_ << endl;
    // cout << "dyCheck is " << dyCheck_ << endl;
    // cout << "dzCheck is " << dzCheck_ << endl;
    // cout << "End light jet chi square constructor" << endl;
}

lightJetChiSquareMinimumSolver::lightJetChiSquareMinimumSolver(
    const lightJetChiSquareMinimumSolver &other)
    : do3D_(other.do3D_), dx_(other.dx_), dy_(other.dy_), dz_(other.dz_),
      inv_sum_cov2D_(other.inv_sum_cov2D_), inv_sum_cov3D_(other.inv_sum_cov3D_)
{
    nJets_ = other.nJets_;
    dxCheck_ = other.dxCheck_;
    dyCheck_ = other.dyCheck_;
    dzCheck_ = other.dzCheck_;
    chi2_ = other.chi2_;

    jetPxWidths2_ = other.jetPxWidths2_;
    jetPyWidths2_ = other.jetPyWidths2_;
    jetPzWidths2_ = other.jetPzWidths2_;
    jetPxPyWidths_ = other.jetPxPyWidths_;
    jetPxPzWidths_ = other.jetPxPzWidths_;
    jetPyPzWidths_ = other.jetPxPzWidths_;
    minDeltasX_ = other.minDeltasX_;
    minDeltasY_ = other.minDeltasY_;
    minDeltasZ_ = other.minDeltasZ_;
    cov2D_ = other.cov2D_;
    cov3D_ = other.cov3D_;
    inv_sum_x_cov2D_ = other.inv_sum_x_cov2D_;
    inv_sum_x_cov3D_ = other.inv_sum_x_cov3D_;

    inverter2D_ = new TDecompLU(2);
    *inverter2D_ = *(other.inverter2D_);
    inverter3D_ = new TDecompLU(3);
    *inverter3D_ = *(other.inverter3D_);
}

lightJetChiSquareMinimumSolver::~lightJetChiSquareMinimumSolver()
{
    delete inverter2D_;
    delete inverter3D_;
}

// void lightJetChiSquareMinimumSolver::setRecoil(double x, double y, double z)
//{
//  cout << "dx is " << dx_ << endl;
//  cout << "dy is " << dy_ << endl;
//  cout << "dz is " << dz_ << endl;
//  cout << "dxCheck is " << dxCheck_ << endl;
//  cout << "dyCheck is " << dyCheck_ << endl;
//  cout << "dzCheck is " << dzCheck_ << endl;
//
//  //cout << "Setting recoil" << endl;
////  dxCheck_=x;
////  dyCheck_=y;
////  if(do3D_) dzCheck_=z;
//
////  cout << "dx is " << dx_ << endl;
////  cout << "dy is " << dy_ << endl;
////  cout << "dz is " << dz_ << endl;
////  cout << "dxCheck is " << dxCheck_ << endl;
////  cout << "dyCheck is " << dyCheck_ << endl;
////  cout << "dzCheck is " << dzCheck_ << endl;
//}

void lightJetChiSquareMinimumSolver::setupEquations(
    vector<XYZTLorentzVector> &jets, vector<double> &jetPtWidths,
    vector<double> &jetPhiWidths, vector<double> &jetEtaWidths)
{
    checkSize(jets, jetPtWidths, jetPhiWidths, jetEtaWidths);
    if (jets.size() != jetPxWidths2_.size()) {
        cerr << "Unequal number of cartesian and radial jets!" << endl;
        return;
    }
    setCartesianWidths(jets, jetPtWidths, jetPhiWidths, jetEtaWidths);
    calcSigmas();
}

void lightJetChiSquareMinimumSolver::setCartesianWidths(
    vector<XYZTLorentzVector> &jets, vector<double> &jetPtWidths,
    vector<double> &jetPhiWidths, vector<double> &jetEtaWidths)
{
    if (debug)
        cout << "lightJetChiSquareMinimumSolver::setCartesianWidths\n";
    for (unsigned int i = 0; i < nJets_; i++) {
        if (debug) {
            cout << "Setting cartesian widths:" << endl;
            cout << "jetPtWidth = " << jetPtWidths[i] << endl;
            cout << "jetPhiWidth = " << jetPhiWidths[i] << endl;
            cout << "jetEtaWidth = " << jetEtaWidths[i] << endl;
        }

        /*      double halfPt2 = 0.5*jets[i].Pt()*jets[i].Pt();
        //      double sigmaPt2 = log(1+jetPtWidths[i]);//WHY?
              double sigmaPt2 = jetPtWidths[i];//try just setting it exactly
        equal to what we put in
              sigmaPt2*=sigmaPt2;
              double sigmaPhi2 = jetPhiWidths[i]*jetPhiWidths[i];
              double sigmaEta = jetEtaWidths[i];
              double sigmaEta2 = sigmaEta*sigmaEta;

              double expSigmaPt2 = exp(sigmaPt2);
              double expTwoSigmaPt2 = expSigmaPt2*expSigmaPt2;
              double expMinusSigmaPhi2 = exp(-sigmaPhi2);
              double expMinusTwoSigmaPhi2 = expMinusSigmaPhi2*expMinusSigmaPhi2;
              double expMinusSigmaPhi2Over2 = exp(-0.5*sigmaPhi2);
              double expSigmaEta2 = exp(sigmaEta2);
              double expTwoSigmaEta2 = expSigmaEta2*expSigmaEta2;
             double expSigmaEta2OverTwo = exp(-0.5*sigmaEta2);//Should really
        have a minus in the name
              double cosPhi = cos(jets[i].Phi());
              double sinPhi = sin(jets[i].Phi());
              double cosTwoPhi = cos(2.*jets[i].Phi());
              double sinTwoPhi = sin(2.*jets[i].Phi());
              double coshTwoEta = cosh(2.*jets[i].Eta());
              double sinhEta  = sinh(jets[i].Eta());
              double sinhEta2 = sinhEta*sinhEta;
              jetPxWidths2_[i]   =  halfPt2* ( expTwoSigmaPt2 * (1 +
        cosTwoPhi*expMinusTwoSigmaPhi2)
                              - expSigmaPt2 * (1 + cosTwoPhi) *
        expMinusSigmaPhi2) ;
              jetPyWidths2_[i]   =  halfPt2* ( expTwoSigmaPt2 * (1 -
        cosTwoPhi*expMinusTwoSigmaPhi2)
                              - expSigmaPt2 * (1 - cosTwoPhi) *
        expMinusSigmaPhi2) ;
              jetPzWidths2_[i]   =  halfPt2* ( expTwoSigmaPt2 *
        (expTwoSigmaEta2*coshTwoEta - 1) //FIXMEis there a typo here?
        expMinusTwoSigmaEta instead?
                              - 2.*expSigmaPt2*expSigmaEta2*sinhEta2 );
              jetPxPyWidths_[i] = halfPt2 * sinTwoPhi *( expTwoSigmaPt2
        *expMinusTwoSigmaPhi2
                                - expSigmaPt2 * expMinusSigmaPhi2);
              jetPxPzWidths_[i] =
        2.*halfPt2*cosPhi*sinhEta*expSigmaEta2OverTwo*expMinusSigmaPhi2Over2 * (
        expTwoSigmaPt2 - expSigmaPt2 );
              jetPyPzWidths_[i] =
        2.*halfPt2*sinPhi*sinhEta*expSigmaEta2OverTwo*expMinusSigmaPhi2Over2 * (
        expTwoSigmaPt2 - expSigmaPt2 );

              //Temp quick "fix" to try setting PxPy to square roots (or
        something like that) of the Px and Py
            double px_temp = jets[i].Pt()*cos(jets[i].Phi());
            double py_temp = jets[i].Pt()*sin(jets[i].Phi());
            jetPxWidths2_[i] = abs(px_temp);
            jetPyWidths2_[i] = abs(py_temp);
            jetPxPyWidths_[i] = sqrt(abs(px_temp*py_temp));
            cout<< "jetPxPyWidths is "<< jetPxPyWidths_[i];
            //jetPyPxWidths_[i] = jetPxPyWidths_[i];*/

        const double Pt2 = pow(jets[i].Pt(), 2);
        const double sigmaPt2 = pow(jetPtWidths[i], 2);
        const double sigmaPhi2 = pow(jetPhiWidths[i], 2);
        const double cosPhi = cos(jets[i].Phi());
        const double sinPhi = sin(jets[i].Phi());
        const double cosPhi2 = pow(cosPhi, 2);
        const double sinPhi2 = pow(sinPhi, 2);
        jetPxWidths2_[i] = sigmaPt2 * cosPhi2 + sigmaPhi2 * Pt2 * sinPhi2;
        jetPyWidths2_[i] = sigmaPt2 * sinPhi2 + sigmaPhi2 * Pt2 * cosPhi2;
        jetPxPyWidths_[i] =
            sigmaPt2 * sinPhi * cosPhi - sigmaPhi2 * Pt2 * sinPhi * cosPhi;
        jetPzWidths2_[i] = 0;
        jetPxPzWidths_[i] = 0;
        jetPyPzWidths_[i] = 0;

        if (debug) {
            cout << "Calculating widths:\n"
                 << "pt  is " << jets[i].Pt() << " with width of "
                 << log(1 + jetPtWidths[i]) << "\n"
                 << "phi is " << jets[i].Phi() << " with width of "
                 << log(1 + jetPhiWidths[i]) << "\n"
                 << "px  is " << jets[i].Pt() * cos(jets[i].Phi())
                 << " with width of " << sqrt(jetPxWidths2_[i]) << "\n"
                 << "py  is " << jets[i].Pt() * sin(jets[i].Phi())
                 << " with width of " << sqrt(jetPyWidths2_[i]) << "\n"
                 << "pxpy  width is " << sqrt(jetPxPyWidths_[i]) << "\n"
                 << "correlation coefficient is "
                 << jetPxPyWidths_[i] /
                        (sqrt(jetPxWidths2_[i]) * sqrt(jetPyWidths2_[i]))
                 << endl;
        }
    }
}

void lightJetChiSquareMinimumSolver::calcSigmas()
{
    inv_sum_cov2D_.Zero();
    inv_sum_cov3D_.Zero();

    if (do3D_)
        return calcSigmas_3D();
    return calcSigmas_2D();
}

inline void lightJetChiSquareMinimumSolver::calcSigmas_2D()
{
    for (unsigned int i = 0; i < nJets_; ++i) {
        const double px2 = jetPxWidths2_[i];
        const double py2 = jetPyWidths2_[i];
        const double pxpy = jetPxPyWidths_[i];

        double array2D[4] = {0.};
        array2D[0] = px2;
        array2D[1] = pxpy;
        array2D[2] = pxpy;
        array2D[3] = py2;
        cov2D_[i] = TMatrixD(2, 2, array2D);
        inv_sum_cov2D_ += cov2D_[i];
    }

    // bool checkDecomp = false;
    dynamic_cast<TDecompLU *>(inverter2D_)->SetMatrix(TMatrixD(inv_sum_cov2D_));
    // checkDecomp = inverter2D_->Decompose();
    dynamic_cast<TDecompLU *>(inverter2D_)->Invert(inv_sum_cov2D_);
    // inv_sum_cov2D_.Print();

    for (unsigned int i = 0; i < nJets_; ++i)
        inv_sum_x_cov2D_[i] = cov2D_[i] * inv_sum_cov2D_;
}

inline void lightJetChiSquareMinimumSolver::calcSigmas_3D()
{
    for (unsigned int i = 0; i < nJets_; ++i) {
        const double px2 = jetPxWidths2_[i];
        const double py2 = jetPyWidths2_[i];
        const double pz2 = jetPzWidths2_[i];
        const double pxpy = jetPxPyWidths_[i];
        const double pxpz = jetPxPzWidths_[i];
        const double pypz = jetPyPzWidths_[i];

        double array3D[9] = {0.};
        array3D[0] = px2;
        array3D[1] = pxpy;
        array3D[2] = pxpz;
        array3D[3] = pxpy;
        array3D[4] = py2;
        array3D[5] = pypz;
        array3D[6] = pxpz;
        array3D[7] = pypz;
        array3D[8] = pz2;
        cov3D_[i] = TMatrixD(3, 3, array3D);
        inv_sum_cov3D_ += cov3D_[i];
    }

    // bool checkDecomp = false;
    dynamic_cast<TDecompLU *>(inverter3D_)->SetMatrix(TMatrixD(inv_sum_cov3D_));
    // checkDecomp = inverter3D_->Decompose();
    dynamic_cast<TDecompLU *>(inverter3D_)->Invert(inv_sum_cov3D_);
    // inv_sum_cov3D_.Print();

    for (unsigned int i = 0; i < nJets_; ++i)
        inv_sum_x_cov3D_[i] = cov3D_[i] * inv_sum_cov3D_;
}

void lightJetChiSquareMinimumSolver::calcMin()
{
    // cout << "dx is " << dx_ << endl;
    // cout << "dy is " << dy_ << endl;
    // cout << "dz is " << dz_ << endl;
    // cout << "dxCheck is " << dxCheck_ << endl;
    // cout << "dyCheck is " << dyCheck_ << endl;
    // cout << "dzCheck is " << dzCheck_ << endl;

    if (do3D_)
        return calcMin_3D();
    return calcMin_2D();
}

inline void lightJetChiSquareMinimumSolver::calcMin_2D()
{
    if (dxCheck_ == dx_ && dyCheck_ == dy_)
        return;

    // cout << "Calculating minimum chi^2" << endl;

    dxCheck_ = dx_;
    dyCheck_ = dy_;

    const double dArray2D[2] = {dx_, dy_};
    const TVectorD dVec2D(2, dArray2D);

    for (unsigned int i = 0; i < nJets_; ++i) {
        const TVectorD jet_delta = inv_sum_x_cov2D_[i] * dVec2D;

        minDeltasX_[i] = jet_delta[0];
        minDeltasY_[i] = jet_delta[1];
    }

    chi2_ = dVec2D * (inv_sum_cov2D_ * dVec2D);
}

inline void lightJetChiSquareMinimumSolver::calcMin_3D()
{
    if (dxCheck_ == dx_ && dyCheck_ == dy_ && dzCheck_ == dz_)
        return;

    // cout << "Calculating minimum chi^2" << endl;

    dxCheck_ = dx_;
    dyCheck_ = dy_;
    dzCheck_ = dz_;

    const double dArray3D[3] = {dx_, dy_, dz_};
    const TVectorD dVec3D(3, dArray3D);

    for (unsigned int i = 0; i < nJets_; ++i) {
        const TVectorD jet_delta = inv_sum_x_cov3D_[i] * dVec3D;

        minDeltasX_[i] = jet_delta[0];
        minDeltasY_[i] = jet_delta[1];
        minDeltasZ_[i] = jet_delta[2];
    }

    chi2_ = dVec3D * (inv_sum_cov3D_ * dVec3D);
}

void lightJetChiSquareMinimumSolver::printResults()
{
    for (unsigned int i = 0; i < nJets_; i++) {
        cout << "delta px " << i + 1 << " = " << minDeltasX_[i] << endl;
        cout << "delta py " << i + 1 << " = " << minDeltasY_[i] << endl;
        // cout<< "jetptwidths = "<< jetPtWidths<<endl;
        if (do3D_) {
            cout << "delta pz " << i + 1 << " = " << minDeltasZ_[i] << endl;
        }
    }
}

double lightJetChiSquareMinimumSolver::getChiSquare()
{
    calcMin();
    //  vector<double>::iterator thisDeltaX = minDeltasX_.begin();
    //  double deltaXCheck(0.);
    //  double deltaYCheck(0.);
    //  for(vector<double>::iterator thisDeltaY = minDeltasY_.begin();
    //  thisDeltaY != minDeltasY_.end(); thisDeltaX++, thisDeltaY++)
    //    {
    //      deltaXCheck+=*thisDeltaX;
    //      deltaYCheck+=*thisDeltaY;
    //    }
    //  cout << "delta x = " << dx_ << " and delta x check = " << deltaXCheck <<
    //  endl;
    //  cout << "delta y = " << dy_ << " and delta y check = " << deltaYCheck <<
    //  endl;
    // cout<<"light jet get chi square"<<endl;
    // printResults();
    // cout<< "light jet chi2 = " << chi2_<<endl;
    return chi2_;
}

void lightJetChiSquareMinimumSolver::checkSize(vector<XYZTLorentzVector> &jets,
                                               vector<double> &jetPtWidths,
                                               vector<double> &jetPhiWidths,
                                               vector<double> &jetEtaWidths)
{
    if (jets.size() != jetPtWidths.size()) {
        cout << "Unequal number of jets and jet pT widths!" << endl;
        cout << "there are " << jets.size() << " jets and "
             << jetPtWidths.size() << " jet pt widths" << endl;
        return;
    }
    if (jets.size() != jetPhiWidths.size()) {
        cout << "Unequal number of jets and jet phi widths!" << endl;
        cout << "there are " << jets.size() << " jets and "
             << jetPhiWidths.size() << " jet phi widths" << endl;
        return;
    }
    if (jets.size() != jetEtaWidths.size()) {
        cout << "Unequal number of jets and jet eta widths!" << endl;
        cout << "there are " << jets.size() << " jets and "
             << jetEtaWidths.size() << " jet eta widths" << endl;
        return;
    }
    return;
}

#endif