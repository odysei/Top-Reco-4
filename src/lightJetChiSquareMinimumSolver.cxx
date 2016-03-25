#ifndef lightJetChiSquareMinimumSolver_cxx
#define lightJetChiSquareMinimumSolver_cxx

#include "lightJetChiSquareMinimumSolver.h"
#include "TDecompSVD.h"
#include "TDecompLU.h"

lightJetChiSquareMinimumSolver::lightJetChiSquareMinimumSolver(
    vector<XYZTLorentzVector> &input_ps, vector<double> &jetPtWidths,
    vector<double> &jetPhiWidths, vector<double> &jetEtaWidths, double &dx,
    double &dy, double &dz, bool do3D)
    : do3D_(do3D)
{
    input_.dx = &dx;
    input_.dy = &dy;
    input_.dz = &dz;
    input_.ps = &input_ps;
    input_.jetPtWidths = &jetPtWidths;
    input_.jetPhiWidths = &jetPhiWidths;
    input_.jetEtaWidths = &jetEtaWidths;
    chi2_ = 0.;

    Init_data(input_ps.size(), data_);

    inverter3D_ = new TDecompLU(3);

    checkSize(input_);
    Eval_covariance(input_, data_);
    Eval_cov_sum(data_);
}

lightJetChiSquareMinimumSolver::lightJetChiSquareMinimumSolver(
    const unsigned int nObjects, double &dx, double &dy, double &dz, bool do3D)
    : do3D_(do3D)
{
    input_.dx = &dx;
    input_.dy = &dy;
    input_.dz = &dz;
    chi2_ = 0.;

    Init_data(nObjects, data_);

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

inline void lightJetChiSquareMinimumSolver::Init_data(const unsigned int &sz,
                                                      recoil_minimizer_data &da)
{
    da.n_ps = sz;
    da.dxCheck = 0.;
    da.dyCheck = 0.;
    da.dzCheck = 0.;

    da.minDeltasX.resize(sz, 0.);
    da.minDeltasY.resize(sz, 0.);
    da.minDeltasZ.resize(sz, 0.);
    da.cov_rad.reserve(sz);
    da.cov.reserve(sz);
    da.inv_sum_x_cov.reserve(sz);
    for (unsigned int i = 0; i < sz; ++i) {
        da.cov_rad.push_back(TMatrixD(3, 3));
        da.cov.push_back(TMatrixD(3, 3));
        da.inv_sum_x_cov.push_back(TMatrixD(3, 3));
    }
}

lightJetChiSquareMinimumSolver::lightJetChiSquareMinimumSolver(
    const lightJetChiSquareMinimumSolver &other)
    : do3D_(other.do3D_)
{
    input_ = other.input_;
    data_ = other.data_;
    chi2_ = other.chi2_;

    inverter3D_ = new TDecompLU(3);
    *inverter3D_ = *(other.inverter3D_);
}

lightJetChiSquareMinimumSolver::~lightJetChiSquareMinimumSolver()
{
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
    vector<XYZTLorentzVector> &input_ps, vector<double> &jetPtWidths,
    vector<double> &jetPhiWidths, vector<double> &jetEtaWidths)
{
    input_.ps = &input_ps;
    input_.jetPtWidths = &jetPtWidths;
    input_.jetPhiWidths = &jetPhiWidths;
    input_.jetEtaWidths = &jetEtaWidths;

    checkSize(input_);
    Eval_covariance(input_, data_);
    Eval_cov_sum(data_);
}

void lightJetChiSquareMinimumSolver::checkSize(recoil_minimizer_input &in)
{
    if (in.ps->size() != in.jetPtWidths->size()) {
        cout << "Unequal number of jets and jet pT widths!" << endl;
        cout << "there are " << in.ps->size() << " jets and "
             << in.jetPtWidths->size() << " jet pt widths" << endl;
        return;
    }
    if (in.ps->size() != in.jetPhiWidths->size()) {
        cout << "Unequal number of jets and jet phi widths!" << endl;
        cout << "there are " << in.ps->size() << " jets and "
             << in.jetPhiWidths->size() << " jet phi widths" << endl;
        return;
    }
    if (in.ps->size() != in.jetEtaWidths->size()) {
        cout << "Unequal number of jets and jet eta widths!" << endl;
        cout << "there are " << in.ps->size() << " jets and "
             << in.jetEtaWidths->size() << " jet eta widths" << endl;
        return;
    }
    return;
}

// void lightJetChiSquareMinimumSolver::Eval_covariance(
//     const recoil_minimizer_input &in, recoil_minimizer_data &da)
// {
// //     if (do3D_)
// //         return Eval_covariance_3D(in, da);
// //     return Eval_covariance_2D(in, da);
//     return Eval_covariance_3D(in, da);
//
//     // old part
//     if (debug)
//         cout << "lightJetChiSquareMinimumSolver::Eval_covariance\n";
//     for (unsigned int i = 0; i < da.n_ps; i++) {
//         if (debug) {
//             cout << "Setting cartesian widths:" << endl;
//             cout << "jetPtWidth = " << (*in.jetPtWidths)[i] << endl;
//             cout << "jetPhiWidth = " << (*in.jetPhiWidths)[i] << endl;
//             cout << "jetEtaWidth = " << (*in.jetEtaWidths)[i] << endl;
//         }
//
//         /*      double halfPt2 = 0.5*jets[i].Pt()*jets[i].Pt();
//         //      double sigmaPt2 = log(1+jetPtWidths[i]);//WHY?
//               double sigmaPt2 = jetPtWidths[i];//try just setting it exactly
//         equal to what we put in
//               sigmaPt2*=sigmaPt2;
//               double sigmaPhi2 = jetPhiWidths[i]*jetPhiWidths[i];
//               double sigmaEta = jetEtaWidths[i];
//               double sigmaEta2 = sigmaEta*sigmaEta;
//
//               double expSigmaPt2 = exp(sigmaPt2);
//               double expTwoSigmaPt2 = expSigmaPt2*expSigmaPt2;
//               double expMinusSigmaPhi2 = exp(-sigmaPhi2);
//               double expMinusTwoSigmaPhi2 =
//               expMinusSigmaPhi2*expMinusSigmaPhi2;
//               double expMinusSigmaPhi2Over2 = exp(-0.5*sigmaPhi2);
//               double expSigmaEta2 = exp(sigmaEta2);
//               double expTwoSigmaEta2 = expSigmaEta2*expSigmaEta2;
//              double expSigmaEta2OverTwo = exp(-0.5*sigmaEta2);//Should really
//         have a minus in the name
//               double cosPhi = cos(jets[i].Phi());
//               double sinPhi = sin(jets[i].Phi());
//               double cosTwoPhi = cos(2.*jets[i].Phi());
//               double sinTwoPhi = sin(2.*jets[i].Phi());
//               double coshTwoEta = cosh(2.*jets[i].Eta());
//               double sinhEta  = sinh(jets[i].Eta());
//               double sinhEta2 = sinhEta*sinhEta;
//               jetPxWidths2_[i]   =  halfPt2* ( expTwoSigmaPt2 * (1 +
//         cosTwoPhi*expMinusTwoSigmaPhi2)
//                               - expSigmaPt2 * (1 + cosTwoPhi) *
//         expMinusSigmaPhi2) ;
//               jetPyWidths2_[i]   =  halfPt2* ( expTwoSigmaPt2 * (1 -
//         cosTwoPhi*expMinusTwoSigmaPhi2)
//                               - expSigmaPt2 * (1 - cosTwoPhi) *
//         expMinusSigmaPhi2) ;
//               jetPzWidths2_[i]   =  halfPt2* ( expTwoSigmaPt2 *
//         (expTwoSigmaEta2*coshTwoEta - 1) //FIXMEis there a typo here?
//         expMinusTwoSigmaEta instead?
//                               - 2.*expSigmaPt2*expSigmaEta2*sinhEta2 );
//               jetPxPyWidths_[i] = halfPt2 * sinTwoPhi *( expTwoSigmaPt2
//         *expMinusTwoSigmaPhi2
//                                 - expSigmaPt2 * expMinusSigmaPhi2);
//               jetPxPzWidths_[i] =
//         2.*halfPt2*cosPhi*sinhEta*expSigmaEta2OverTwo*expMinusSigmaPhi2Over2
//         * (
//         expTwoSigmaPt2 - expSigmaPt2 );
//               jetPyPzWidths_[i] =
//         2.*halfPt2*sinPhi*sinhEta*expSigmaEta2OverTwo*expMinusSigmaPhi2Over2
//         * (
//         expTwoSigmaPt2 - expSigmaPt2 );
//
//               //Temp quick "fix" to try setting PxPy to square roots (or
//         something like that) of the Px and Py
//             double px_temp = jets[i].Pt()*cos(jets[i].Phi());
//             double py_temp = jets[i].Pt()*sin(jets[i].Phi());
//             jetPxWidths2_[i] = abs(px_temp);
//             jetPyWidths2_[i] = abs(py_temp);
//             jetPxPyWidths_[i] = sqrt(abs(px_temp*py_temp));
//             cout<< "jetPxPyWidths is "<< jetPxPyWidths_[i];
//             //jetPyPxWidths_[i] = jetPxPyWidths_[i];*/
//
// //         // newer version
// //         const double Pt2 = pow((*in.ps)[i].Pt(), 2);
// //         const double sigmaPt2 = pow((*in.jetPtWidths)[i], 2);
// //         const double sigmaPhi2 = pow((*in.jetPhiWidths)[i], 2);
// //         const double cosPhi = cos((*in.ps)[i].Phi());
// //         const double sinPhi = sin((*in.ps)[i].Phi());
// //         const double cosPhi2 = pow(cosPhi, 2);
// //         const double sinPhi2 = pow(sinPhi, 2);
// //         jetPxWidths2_[i] = sigmaPt2 * cosPhi2 + sigmaPhi2 * Pt2 * sinPhi2;
// //         jetPyWidths2_[i] = sigmaPt2 * sinPhi2 + sigmaPhi2 * Pt2 * cosPhi2;
// //         jetPxPyWidths_[i] =
// //             sigmaPt2 * sinPhi * cosPhi - sigmaPhi2 * Pt2 * sinPhi *
// cosPhi;
// //         jetPzWidths2_[i] = 0;
// //         jetPxPzWidths_[i] = 0;
// //         jetPyPzWidths_[i] = 0;
// //
// //         if (debug) {
// //             cout << "Calculating widths:\n"
// //                  << "pt  is " << (*in.ps)[i].Pt() << " with width of "
// //                  << (*in.jetPtWidths)[i] << "\n"
// //                  << "phi is " << (*in.ps)[i].Phi() << " with width of "
// //                  << (*in.jetPhiWidths)[i] << "\n"
// //                  << "px  is " << (*in.ps)[i].Pt() * cos((*in.ps)[i].Phi())
// //                  << " with width of " << sqrt(jetPxWidths2_[i]) << "\n"
// //                  << "py  is " << (*in.ps)[i].Pt() * sin((*in.ps)[i].Phi())
// //                  << " with width of " << sqrt(jetPyWidths2_[i]) << "\n"
// //                  << "pxpy  width is " << sqrt(jetPxPyWidths_[i]) << "\n"
// //                  << "correlation coefficient is "
// //                  << jetPxPyWidths_[i] /
// //                         (sqrt(jetPxWidths2_[i]) * sqrt(jetPyWidths2_[i]))
// //                  << endl;
// //         }
//     }
// }

inline void lightJetChiSquareMinimumSolver::Eval_covariance(
    const recoil_minimizer_input &in, recoil_minimizer_data &da)
{
    if (debug)
        cout << "lightJetChiSquareMinimumSolver::Eval_covariance\n";
    for (unsigned int i = 0; i < da.n_ps; i++) {
        if (debug) {
            cout << "Setting cartesian widths:" << endl;
            cout << "jetPtWidth = " << (*in.jetPtWidths)[i] << endl;
            cout << "jetPhiWidth = " << (*in.jetPhiWidths)[i] << endl;
            cout << "jetEtaWidth = " << (*in.jetEtaWidths)[i] << endl;
        }
        const double pT = (*in.ps)[i].Pt();
        const double phi = (*in.ps)[i].Phi();
        const double eta = (*in.ps)[i].Eta();

        TMatrixD cov(3, 3);
        TMatrixD R(3, 3);
        TMatrixD RT(3, 3);

        cov[0][0] = pow((*in.jetPtWidths)[i], 2);
        cov[0][1] = 0;
        cov[0][2] = 0;
        cov[1][0] = 0;
        cov[1][1] = pow((*in.jetPhiWidths)[i], 2);
        cov[1][2] = 0;
        cov[2][0] = 0;
        cov[2][1] = 0;
        cov[2][2] = pow((*in.jetEtaWidths)[i], 2);
        R[0][0] = cos(phi);
        R[0][1] = -pT * sin(phi);
        R[0][2] = 0;
        R[1][0] = sin(phi);
        R[1][1] = pT * cos(phi);
        R[1][2] = 0;
        R[2][0] = tan(eta);
        R[2][1] = 0;
        R[2][2] = pT * (1 + pow(tan(eta), 2));
        RT.Transpose(R);

        // tranform to cartesian coordinate system
        da.cov[i] = R * cov * RT;
        da.cov_rad[i] = cov;

        if (debug) {
            cout << "Calculating widths:\n"
                 << "pT: " << pT << " +- " << (*in.jetPtWidths)[i] << endl
                 << "phi: " << phi << " +- " << (*in.jetPhiWidths)[i] << endl
                 << "eta: " << eta << " +- " << (*in.jetEtaWidths)[i] << endl
                 << "px: " << pT * cos(phi) << " +- " << sqrt(cov[0][0]) << endl
                 << "py: " << pT * sin(phi) << " +- " << sqrt(cov[1][1]) << endl
                 << "pz: " << pT * tan(eta) << " +- " << sqrt(cov[1][1])
                 << endl;

            cout << "Covariance matrix in cartesian coords:\n";
            da.cov[i].Print();
        }
    }
}

inline void
lightJetChiSquareMinimumSolver::Eval_cov_sum(recoil_minimizer_data &da)
{
    da.inv_sum_cov.ResizeTo(3, 3);
    da.inv_sum_cov.Zero();
    for (unsigned int i = 0; i < da.n_ps; ++i)
        da.inv_sum_cov += da.cov[i];

    dynamic_cast<TDecompLU *>(inverter3D_)->SetMatrix(TMatrixD(da.inv_sum_cov));
    // checkDecomp = inverter3D_->Decompose();
    dynamic_cast<TDecompLU *>(inverter3D_)->Invert(da.inv_sum_cov);
    // da.inv_sum_cov3D.Print();

    for (unsigned int i = 0; i < da.n_ps; ++i)
        da.inv_sum_x_cov[i] = da.cov[i] * da.inv_sum_cov;
}

void lightJetChiSquareMinimumSolver::calcMin()
{
    // cout << "dx is " << dx_ << endl;
    // cout << "dy is " << dy_ << endl;
    // cout << "dz is " << dz_ << endl;
    // cout << "dxCheck is " << dxCheck_ << endl;
    // cout << "dyCheck is " << dyCheck_ << endl;
    // cout << "dzCheck is " << dzCheck_ << endl;

    if (do3D_) {
        if (data_.dxCheck == *(input_.dx) && data_.dyCheck == *(input_.dy) &&
            data_.dzCheck == *(input_.dz))
            return;
    } else {
        if (data_.dxCheck == *(input_.dx) && data_.dyCheck == *(input_.dy))
            return;
    }

    // cout << "Calculating minimum chi^2" << endl;

    data_.dxCheck = *(input_.dx);
    data_.dyCheck = *(input_.dy);
    if (do3D_)
        data_.dzCheck = *(input_.dz);
    else
        data_.dzCheck = 0;

    const double dp_arr[3] = {data_.dxCheck, data_.dyCheck, data_.dzCheck};
    const TVectorD dp(3, dp_arr);

    for (unsigned int i = 0; i < data_.n_ps; ++i) {
        const TVectorD dp_i = data_.inv_sum_x_cov[i] * dp;

        data_.minDeltasX[i] = dp_i[0];
        data_.minDeltasY[i] = dp_i[1];
        data_.minDeltasZ[i] = dp_i[2];
    }

    chi2_ = dp * (data_.inv_sum_cov * dp);
}

void lightJetChiSquareMinimumSolver::printResults()
{
    for (unsigned int i = 0; i < data_.n_ps; ++i) {
        cout << "delta px " << i + 1 << " = " << data_.minDeltasX[i] << endl;
        cout << "delta py " << i + 1 << " = " << data_.minDeltasY[i] << endl;
        cout << "delta pz " << i + 1 << " = " << data_.minDeltasZ[i] << endl;
    }
}

double lightJetChiSquareMinimumSolver::getChiSquare()
{
    calcMin();
    //  vector<double>::iterator thisDeltaX = minDeltasX.begin();
    //  double deltaXCheck(0.);
    //  double deltaYCheck(0.);
    //  for(vector<double>::iterator thisDeltaY = minDeltasY.begin();
    //  thisDeltaY != minDeltasY.end(); thisDeltaX++, thisDeltaY++)
    //    {
    //      deltaXCheck+=*thisDeltaX;
    //      deltaYCheck+=*thisDeltaY;
    //    }
    //  cout << "delta x = " << dx_ << " and delta x check = " << deltaXCheck <<
    //  endl;
    //  cout << "delta y = " << dy_ << " and delta y check = " << deltaYCheck <<
    //  endl;
    // cout<<"light jet get chi square"<<endl;
    // cout<< "light jet chi2 = " << chi2_<<endl;
    return chi2_;
}

#endif