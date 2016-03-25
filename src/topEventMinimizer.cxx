#include "topEventMinimizer.h"

#include <numeric>

#include <TH1D.h>
#include <TGraph.h>
#include <TCanvas.h>

topEventMinimizer::topEventMinimizer(vector<XYZTLorentzVector> nonTopObjects,
                                     vector<double> nonTopObjectPtWidths,
                                     vector<double> nonTopObjectPhiWidths,
                                     vector<double> nonTopObjectEtaWidths,
                                     double mTop, double sigmaMTop, double mW,
                                     double sigmaMW, double totalTopPx,
                                     double totalTopPy, double totalTopPz)
    : dx_(0.), dy_(0.), dz_(0.),
      nonTopChiSquare_(lightJetChiSquareMinimumSolver(nonTopObjects.size(), dx_,
                                                      dy_, dz_, false))
{
    // cout << "Basic constructor" << endl;
    nTops_ = 0;
    totalTopPx_ = totalTopPx;
    totalTopPy_ = totalTopPy;
    totalTopPz_ = totalTopPz;
    nonTopObjects_ = nonTopObjects;
    nonTopObjectPtWidths_ = nonTopObjectPtWidths;
    nonTopObjectPhiWidths_ = nonTopObjectPhiWidths;
    nonTopObjectEtaWidths_ = nonTopObjectEtaWidths;
    mTop_ = mTop;
    sigmaMTop_ = sigmaMTop;
    mW_ = mW;
    sigmaMW_ = sigmaMW;
    maxConsideredChiSquareRoot_ = 30;

    setupMap();

    calcNonTopMomentum();

    nonTopChiSquare_.setupEquations(nonTopObjects_, nonTopObjectPtWidths_,
                                    nonTopObjectPhiWidths_,
                                    nonTopObjectEtaWidths_);

    nonTopObjects_PxDeltasBest_ = vector<double>(nonTopObjects.size(), 0.);
    nonTopObjects_PyDeltasBest_ = vector<double>(nonTopObjects.size(), 0.);
    nonTopObjects_PzDeltasBest_ = vector<double>(nonTopObjects.size(), 0.);

    initializeChiSquares();
    Initialize_minimizers(outerMin_, innerMin_);

    innerMinStatus = -1;
    outerMinStatus = -1;
    outerMin_Edm = -1;
}

topEventMinimizer::topEventMinimizer(
    vector<XYZTLorentzVector> allObjects, vector<double> allObjectPtWidths,
    vector<double> allObjectPhiWidths, vector<double> allObjectEtaWidths,
    vector<int> bJets, vector<int> firstWDaughters,
    vector<int> secondWDaughters, vector<bool> isLeptonicTopDecay, double mTop,
    double sigmaMTop, double mW, double sigmaMW, double totalTopPx,
    double totalTopPy, double totalTopPz)
    : dx_(0.), dy_(0.), dz_(0.),
      nonTopChiSquare_(lightJetChiSquareMinimumSolver(
          allObjects.size() - (int)accumulate(isLeptonicTopDecay.begin(),
                                              isLeptonicTopDecay.end(), 0),
          dx_, dy_, dz_, false))
{
    // cout << "constructor with input tops" << endl;
    nTops_ = 0;
    bJets_ = bJets;
    firstWDaughters_ = firstWDaughters;
    secondWDaughters_ = secondWDaughters;
    allObjects_ = allObjects;
    allObjectPtWidths_ = allObjectPtWidths;
    allObjectPhiWidths_ = allObjectPhiWidths;
    allObjectEtaWidths_ = allObjectEtaWidths;
    isLeptonicTopDecay_ = isLeptonicTopDecay;
    mTop_ = mTop;
    sigmaMTop_ = sigmaMTop;
    mW_ = mW;
    sigmaMW_ = sigmaMW;
    totalTopPx_ = totalTopPx;
    totalTopPy_ = totalTopPy;
    totalTopPz_ = totalTopPz;
    maxConsideredChiSquareRoot_ = 30;

    setupMap();

    if (!checkInputSizes())
        return;

    // setBJets();
    // setWDaughters();
    setNonTopObjectCollections();

    nonTopChiSquare_.setupEquations(nonTopObjects_, nonTopObjectPtWidths_,
                                    nonTopObjectPhiWidths_,
                                    nonTopObjectEtaWidths_);

    // make the tops
    for (int iTop = 0; iTop < (int)bJets.size(); iTop++) {
        if (isLeptonicTopDecay.at(iTop)) {
            addLeptonicTop(bJets.at(iTop), firstWDaughters.at(iTop));
        } else {
            addHadronicTop(bJets.at(iTop), firstWDaughters.at(iTop),
                           secondWDaughters.at(iTop));
        }
    }

    initializeChiSquares();
    initializeDeltas();
    Initialize_minimizers(outerMin_, innerMin_);

    innerMinStatus = -1;
    outerMinStatus = -1;
    outerMin_Edm = -1;
}

void topEventMinimizer::Initialize_minimizers(ROOT::Math::Minimizer *&outer,
                                              ROOT::Math::Minimizer *&inner)
{
    outer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    outer->SetMaxFunctionCalls(1000000);
    outer->SetTolerance(0.001);
    outer->SetPrintLevel(5);

    inner = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    inner->SetMaxFunctionCalls(1000000);
    inner->SetTolerance(0.001);
    inner->SetPrintLevel(0);
}

topEventMinimizer::~topEventMinimizer()
{
    // cout << "destructor" << endl;

    delete outerMin_;
    delete innerMin_;
}

void topEventMinimizer::setupMap()
{
    // typedef void *FnTestPtr();
    // map < string, FnTestPtr > testMap;
    // testMap["lah"] = topEventMinimizer::testfunc;
    // typedef void (*FnPtr)(int, double, double, double, double);
    // map < string, FnPtr > funcMap;
    funcMap["getTop"] = &topEventMinimizer::getTop;
    funcMap["getW"] = &topEventMinimizer::getW;
    funcMap["getNonTopObject4"] = &topEventMinimizer::getNonTopObject4;
    funcMap["getBJet"] = &topEventMinimizer::getBJet;
    funcMap["getWDaughter1"] = &topEventMinimizer::getWDaughter1;
    funcMap["getWDaughter2"] = &topEventMinimizer::getWDaughter2;
}

void topEventMinimizer::testfunc() {}

bool topEventMinimizer::checkInputSizes()
{
    if (bJets_.size() != firstWDaughters_.size() ||
        bJets_.size() != secondWDaughters_.size() ||
        firstWDaughters_.size() != secondWDaughters_.size()) {
        cout << "Check sizes of input b-jets and W daughter vectors" << endl;
        return false;
    }
    return true;
}

// void topEventMinimizer::setBJets()
//{
////  bJetPxs_.clear();
////  bJetPys_.clear();
////  bJetPzs_.clear();
////  bJetEs_.clear();
////
////  bJetPtWidths_.clear();
////  bJetEtaWidths_.clear();
////  bJetPhiWidths_.clear();
////
////  for(int i = 0; i < nTops_; i++ )
////    {
////      cout << "Top " << i+1 << " b-jet has object index: " << bJets_.at(i)
///<< endl;
////
////      bJetPxs_.push_back( allObjects_.at(bJets_.at(i)).Px() );
////      bJetPys_.push_back( allObjects_.at(bJets_.at(i)).Py() );
////      bJetPzs_.push_back( allObjects_.at(bJets_.at(i)).Pz() );
////      bJetEs_ .push_back( allObjects_.at(bJets_.at(i)).E()  );
////
////      bJetPtWidths_ .push_back( allObjectPtWidths_ .at(i) );
////      bJetEtaWidths_.push_back( allObjectEtaWidths_.at(i) );
////      bJetPhiWidths_.push_back( allObjectPhiWidths_.at(i) );
////    }
//
//}

void topEventMinimizer::setNonTopObjectCollections()
{
    for (int iObj = 0; iObj < int(allObjects_.size()); iObj++) {
        for (int i = 0; i < nTops_; i++) {
            if (iObj == bJets_.at(i))
                continue;
            if (iObj == firstWDaughters_.at(i))
                continue;
            if (iObj == secondWDaughters_.at(i))
                continue;
        }
        nonTopObjects_.push_back(allObjects_.at(iObj));
        nonTopObjectPts_.push_back(allObjects_.at(iObj).Pt());
        nonTopObjectPhis_.push_back(allObjects_.at(iObj).Phi());
        nonTopObjectPtWidths_.push_back(allObjectPtWidths_.at(iObj));
        nonTopObjectPhiWidths_.push_back(allObjectPhiWidths_.at(iObj));
        nonTopObjectEtaWidths_.push_back(allObjectEtaWidths_.at(iObj));

        nonTopObjects_PxDeltasBest_.push_back(0.);
        nonTopObjects_PyDeltasBest_.push_back(0.);
        nonTopObjects_PzDeltasBest_.push_back(0.);
    }

    calcNonTopMomentum();
}

// void topEventMinimizer::setWDaughters()
//{
////  WDaughterPxs_.clear();
////  WDaughterPys_.clear();
////  WDaughterPzs_.clear();
////  WDaughterEs_.clear();
////
////  WDaughterPtWidths_.clear();
////  WDaughterEtaWidths_.clear();
////  WDaughterPhiWidths_.clear();
////
////  for(int i = 0; i < nTops_; i++ )
////    {
////      cout << "Top " << i+1 << " 1st W daughter has object index: " <<
/// firstWDaughters_ .at(i) << endl;
////      cout << "Top " << i+1 << " 2nd W daughter has object index: " <<
/// secondWDaughters_.at(i) << endl;
////
////
/// WDaughterPxs_.push_back(make_pair(allObjects_.at(firstWDaughters_.at(i)).Px(),(secondWDaughters_.at(i)<0)?0.:allObjects_.at(secondWDaughters_.at(i)).Px()));
////
/// WDaughterPys_.push_back(make_pair(allObjects_.at(firstWDaughters_.at(i)).Py(),(secondWDaughters_.at(i)<0)?0.:allObjects_.at(secondWDaughters_.at(i)).Py()));
////
/// WDaughterPzs_.push_back(make_pair(allObjects_.at(firstWDaughters_.at(i)).Pz(),(secondWDaughters_.at(i)<0)?0.:allObjects_.at(secondWDaughters_.at(i)).Pz()));
////      WDaughterEs_
///.push_back(make_pair(allObjects_.at(firstWDaughters_.at(i)).E()
///,(secondWDaughters_.at(i)<0)?0.:allObjects_.at(secondWDaughters_.at(i)).E()
///));
////
////      WDaughterPtWidths_
///.push_back(make_pair(allObjectPtWidths_.at(firstWDaughters_.at(i)),(secondWDaughters_.at(i)<0)?0.:allObjectPtWidths_.at(secondWDaughters_.at(i))));
////
/// WDaughterEtaWidths_.push_back(make_pair(allObjectEtaWidths_.at(firstWDaughters_.at(i)),(secondWDaughters_.at(i)<0)?0.:allObjectEtaWidths_.at(secondWDaughters_.at(i))));
////
/// WDaughterPhiWidths_.push_back(make_pair(allObjectPhiWidths_.at(firstWDaughters_.at(i)),(secondWDaughters_.at(i)<0)?0.:allObjectPhiWidths_.at(secondWDaughters_.at(i))));
////    }
//}

void topEventMinimizer::initializeDeltas()
{
    if (nTops_ == 0) {
        cerr << "No tops have been added yet!" << endl;
        return;
    }

    bJets_PtDeltas_ = vector<double>(nTops_, 0.);
    bJets_PhiDeltas_ = vector<double>(nTops_, 0.);
    bJets_EtaDeltas_ = vector<double>(nTops_, 0.);
    firstWDaughters_PtDeltas_ = vector<double>(nTops_, 0.);
    firstWDaughters_PhiDeltas_ = vector<double>(nTops_, 0.);
    firstWDaughters_EtaDeltas_ = vector<double>(nTops_, 0.);
    topMassDeltas_ = vector<double>(nTops_, 0.);
    WMassDeltas_ = vector<double>(nTops_, 0.);

    bJets_PtDeltasBest_ = vector<double>(nTops_, 0.);
    bJets_PhiDeltasBest_ = vector<double>(nTops_, 0.);
    bJets_EtaDeltasBest_ = vector<double>(nTops_, 0.);
    firstWDaughters_PtDeltasBest_ = vector<double>(nTops_, 0.);
    firstWDaughters_PhiDeltasBest_ = vector<double>(nTops_, 0.);
    firstWDaughters_EtaDeltasBest_ = vector<double>(nTops_, 0.);
    secondWDaughters_PtDeltasBest_ = vector<double>(nTops_, 0.);
    secondWDaughters_PhiDeltasBest_ = vector<double>(nTops_, 0.);
    secondWDaughters_EtaDeltasBest_ = vector<double>(nTops_, 0.);
    topMassDeltasBest_ = vector<double>(nTops_, 0.);
    WMassDeltasBest_ = vector<double>(nTops_, 0.);

    ellipseAnglesBest_ = vector<double>(nTops_, 0.);

    ellipseAnglesCurrent_ = vector<double>(nTops_, 0.);
    topMassDeltasCurrent_ = vector<double>(nTops_, 0.);

    ellipseAnglesInnerBest_ = vector<double>(nTops_, 0.);
    topMassDeltasInnerBest_ = vector<double>(nTops_, 0.);
}

void topEventMinimizer::initializeChiSquares()
{
    chi2_ = 0.;
    nonTopChi2_ = 0.;
    hadChi2_ = 0.;
    topChi2_ = 0.;
    topMassChi2_ = 0.;

    chi2Best_ = 1.e99;
    innerChi2Best_ = 1.e99;
    topChi2Best_ = 0.;
    hadChi2Best_ = 0.;
    topMassChi2Best_ = 0.;
    nonTopChi2Best_ = 0.;

    thisInnerChi2Best_ = 1.e99;
    thisTopMassChi2Best_ = 0.;
    thisNonTopChi2Best_ = 0.;
    thisHadChi2Best_ = 0.;
}

topSystemChiSquare *topEventMinimizer::makeLeptonicTop(int ibJet,
                                                       int iWDaughter1)
{
    return dynamic_cast<topSystemChiSquare *>(new leptonicTopSystemChiSquare(
        allObjects_.at(ibJet).Px(), allObjects_.at(ibJet).Py(),
        allObjects_.at(ibJet).Pz(), allObjects_.at(ibJet).E(),
        allObjectPtWidths_.at(ibJet), allObjectEtaWidths_.at(ibJet),
        allObjectPhiWidths_.at(ibJet), allObjects_.at(iWDaughter1).Px(),
        allObjects_.at(iWDaughter1).Py(), allObjects_.at(iWDaughter1).Pz(),
        allObjects_.at(iWDaughter1).E(), allObjectPtWidths_.at(iWDaughter1),
        allObjectEtaWidths_.at(iWDaughter1),
        allObjectPhiWidths_.at(iWDaughter1), mTop_, sigmaMTop_, mW_, sigmaMW_));
}

topSystemChiSquare *topEventMinimizer::makeLeptonicTop(
    double bJetPx, double bJetPy, double bJetPz, double bJetE,
    double bJetPtWidth, double bJetEtaWidth, double bJetPhiWidth,
    double WDaughter1Px, double WDaughter1Py, double WDaughter1Pz,
    double WDaughter1E, double WDaughter1PtWidth, double WDaughter1EtaWidth,
    double WDaughter1PhiWidth, double mTop, double sigmaMTop, double mW,
    double sigmaMW)
{
    return dynamic_cast<topSystemChiSquare *>(new leptonicTopSystemChiSquare(
        bJetPx, bJetPy, bJetPz, bJetE, bJetPtWidth, bJetEtaWidth, bJetPhiWidth,
        WDaughter1Px, WDaughter1Py, WDaughter1Pz, WDaughter1E,
        WDaughter1PtWidth, WDaughter1EtaWidth, WDaughter1PhiWidth, mTop,
        sigmaMTop, mW, sigmaMW));
}

topSystemChiSquare *
topEventMinimizer::makeHadronicTop(int ibJet, int iWDaughter1, int iWDaughter2)
{
    return dynamic_cast<topSystemChiSquare *>(new hadronicTopSystemChiSquare(
        allObjects_.at(ibJet).Px(), allObjects_.at(ibJet).Py(),
        allObjects_.at(ibJet).Pz(), allObjects_.at(ibJet).E(),
        allObjectPtWidths_.at(ibJet), allObjectEtaWidths_.at(ibJet),
        allObjectPhiWidths_.at(ibJet), allObjects_.at(iWDaughter1).Px(),
        allObjects_.at(iWDaughter1).Py(), allObjects_.at(iWDaughter1).Pz(),
        allObjects_.at(iWDaughter1).E(), allObjectPtWidths_.at(iWDaughter1),
        allObjectEtaWidths_.at(iWDaughter1),
        allObjectPhiWidths_.at(iWDaughter1), allObjects_.at(iWDaughter2).Px(),
        allObjects_.at(iWDaughter2).Py(), allObjects_.at(iWDaughter2).Pz(),
        allObjects_.at(iWDaughter2).E(), allObjectPtWidths_.at(iWDaughter2),
        allObjectEtaWidths_.at(iWDaughter2),
        allObjectPhiWidths_.at(iWDaughter2), mTop_, sigmaMTop_, mW_, sigmaMW_));
}

topSystemChiSquare *topEventMinimizer::makeHadronicTop(
    double bJetPx, double bJetPy, double bJetPz, double bJetE,
    double bJetPtWidth, double bJetEtaWidth, double bJetPhiWidth,
    double WDaughter1Px, double WDaughter1Py, double WDaughter1Pz,
    double WDaughter1E, double WDaughter1PtWidth, double WDaughter1EtaWidth,
    double WDaughter1PhiWidth, double WDaughter2Px, double WDaughter2Py,
    double WDaughter2Pz, double WDaughter2E, double WDaughter2PtWidth,
    double WDaughter2EtaWidth, double WDaughter2PhiWidth, double mTop,
    double sigmaMTop, double mW, double sigmaMW)
{
    return dynamic_cast<topSystemChiSquare *>(new hadronicTopSystemChiSquare(
        bJetPx, bJetPy, bJetPz, bJetE, bJetPtWidth, bJetEtaWidth, bJetPhiWidth,
        WDaughter1Px, WDaughter1Py, WDaughter1Pz, WDaughter1E,
        WDaughter1PtWidth, WDaughter1EtaWidth, WDaughter1PhiWidth, WDaughter2Px,
        WDaughter2Py, WDaughter2Pz, WDaughter2E, WDaughter2PtWidth,
        WDaughter2EtaWidth, WDaughter2PhiWidth, mTop, sigmaMTop, mW, sigmaMW));
}

void topEventMinimizer::addLeptonicTop(int ibJet, int iWDaughter1)
{
    topSystemChiSquares_.push_back(
        make_pair(makeLeptonicTop(ibJet, iWDaughter1), true));
    nTops_ += 1;
}

void topEventMinimizer::addLeptonicTop(
    double bJetPx, double bJetPy, double bJetPz, double bJetE,
    double bJetPtWidth, double bJetEtaWidth, double bJetPhiWidth,
    double WDaughter1Px, double WDaughter1Py, double WDaughter1Pz,
    double WDaughter1E, double WDaughter1PtWidth, double WDaughter1EtaWidth,
    double WDaughter1PhiWidth, double mTop, double sigmaMTop, double mW,
    double sigmaMW)
{
    topSystemChiSquares_.push_back(make_pair(
        makeLeptonicTop(bJetPx, bJetPy, bJetPz, bJetE, bJetPtWidth,
                        bJetEtaWidth, bJetPhiWidth, WDaughter1Px, WDaughter1Py,
                        WDaughter1Pz, WDaughter1E, WDaughter1PtWidth,
                        WDaughter1EtaWidth, WDaughter1PhiWidth, mTop, sigmaMTop,
                        mW, sigmaMW),
        true));

    nTops_ += 1;
}

void topEventMinimizer::addHadronicTop(int ibJet, int iWDaughter1,
                                       int iWDaughter2)
{
    topSystemChiSquares_.push_back(
        make_pair(makeHadronicTop(ibJet, iWDaughter1, iWDaughter2), false));
    nTops_ += 1;
}

void topEventMinimizer::addHadronicTop(
    double bJetPx, double bJetPy, double bJetPz, double bJetE,
    double bJetPtWidth, double bJetEtaWidth, double bJetPhiWidth,
    double WDaughter1Px, double WDaughter1Py, double WDaughter1Pz,
    double WDaughter1E, double WDaughter1PtWidth, double WDaughter1EtaWidth,
    double WDaughter1PhiWidth, double WDaughter2Px, double WDaughter2Py,
    double WDaughter2Pz, double WDaughter2E, double WDaughter2PtWidth,
    double WDaughter2EtaWidth, double WDaughter2PhiWidth, double mTop,
    double sigmaMTop, double mW, double sigmaMW)
{
    topSystemChiSquares_.push_back(make_pair(
        makeHadronicTop(bJetPx, bJetPy, bJetPz, bJetE, bJetPtWidth,
                        bJetEtaWidth, bJetPhiWidth, WDaughter1Px, WDaughter1Py,
                        WDaughter1Pz, WDaughter1E, WDaughter1PtWidth,
                        WDaughter1EtaWidth, WDaughter1PhiWidth, WDaughter2Px,
                        WDaughter2Py, WDaughter2Pz, WDaughter2E,
                        WDaughter2PtWidth, WDaughter2EtaWidth,
                        WDaughter2PhiWidth, mTop, sigmaMTop, mW, sigmaMW),
        false));
    nTops_ += 1;
}

void topEventMinimizer::printTopConstituents()
{
    int iTop = 0;
    for (auto it = topSystemChiSquares_.begin();
         it != topSystemChiSquares_.end(); ++it, ++iTop) {
        cout << "This is top number " << iTop + 1 << endl;
        (*it).first->printTopConstituents();
    }
}

void topEventMinimizer::calcTopMassRanges()
{
    int iTop = 0;
    for (auto it = topSystemChiSquares_.begin();
         it != topSystemChiSquares_.end(); ++it, ++iTop) {
        // cout << "This is top number " << iTop+1 << endl;
        if (!((*it).first->hasTopMassRange()))
            (*it).first->calcTopMassRange();
    }
}

void topEventMinimizer::printNonTopObjects()
{
    for (int iObj = 0; iObj < int(nonTopObjects_.size()); iObj++) {
        cout << "Light jet " << iObj + 1
             << "\npx = " << nonTopObjects_.at(iObj).Px()
             << "\npy = " << nonTopObjects_.at(iObj).Py()
             << "\npz = " << nonTopObjects_.at(iObj).Pz() << endl;
    }
}

void topEventMinimizer::buildBestNonTopObjects()
{
    // nonTopChiSquare_.printResults();

    vector<double> *minJetDeltasX = nonTopChiSquare_.getMinDeltasX();
    vector<double> *minJetDeltasY = nonTopChiSquare_.getMinDeltasY();
    vector<double> *minJetDeltasZ = nonTopChiSquare_.getMinDeltasZ();
    if (minJetDeltasX->size() != nonTopObjects_.size()) {
        cout << "bad size in light jet modification check!" << endl;
        return;
    }
    auto itx = minJetDeltasX->begin();
    auto ity = minJetDeltasY->begin();
    auto itz = minJetDeltasZ->begin();
    int i = 0;
    for (auto it = nonTopObjects_.begin(); it != nonTopObjects_.end();
         ++it, ++itx, ++ity, ++i) {
        // cout << "This is jet number " << i+1 << endl;
        nonTopObjects_PxDeltasBest_[i] = *itx;
        nonTopObjects_PyDeltasBest_[i] = *ity;
        nonTopObjects_PzDeltasBest_[i] = *itz;
    }
}

void topEventMinimizer::calcWDaughterEllipses()
{
    int iTop = 0;
    for (vector<pair<topSystemChiSquare *, bool>>::const_iterator
             thisTopChiSquare = topSystemChiSquares_.begin();
         thisTopChiSquare != topSystemChiSquares_.end();
         thisTopChiSquare++, iTop++) {
        // cout << "Calculating second W daughter ellipse for top number " <<
        // iTop+1 << endl;
        (*thisTopChiSquare).first->calcWDaughter2Ellipse();
    }
}

void topEventMinimizer::calcNonTopMomentum()
{
    nonTopPx_ = 0;
    nonTopPy_ = 0;
    nonTopPz_ = 0;

    for (int iObj = 0; iObj < int(nonTopObjects_.size()); iObj++) {
        nonTopPx_ += nonTopObjects_.at(iObj).Px();
        nonTopPy_ += nonTopObjects_.at(iObj).Py();
        nonTopPz_ += nonTopObjects_.at(iObj).Pz();
    }
}

void topEventMinimizer::setupNonTopChiSquare()
{
    getDxDyFromEllipses();
    nonTopChi2_ = nonTopChiSquare_.getChiSquare();
}

void topEventMinimizer::setRecoil(double px, double py, double pz)
{
    dx_ = -1. * px;
    dy_ = -1. * py;
    dz_ = -1. * pz;

    // cout << "dx is " << dx_ << endl;
    // cout << "dy is " << dy_ << endl;
    // cout << "dz is " << dz_ << endl;
}

void topEventMinimizer::getDxDyFromEllipses()
{
    double sumTopPx = 0., topPx = 0.;
    double sumTopPy = 0., topPy = 0.;
    double sumTopPz = 0., topPz = 0.;
    double topE = 0.;

    int iTop = 0;

    for (vector<pair<topSystemChiSquare *, bool>>::const_iterator
             thisTopChiSquare = topSystemChiSquares_.begin();
         thisTopChiSquare != topSystemChiSquares_.end();
         thisTopChiSquare++, iTop++) {
        // cout << "This is top number " << iTop+1 << endl;

        // Add this top momentum to sum
        (*thisTopChiSquare).first->getTopMomentum(topPx, topPy, topPz, topE);
        sumTopPx += topPx;
        sumTopPy += topPy;
        sumTopPz += topPz;

        // If the top decays hadronically, calculate the second W daughter (jet)
        // deltas at this point on the ellipse
        if (!((*thisTopChiSquare).second)) {
            (*thisTopChiSquare).first->calcWDaughter2Deltas();
        }
    }

    double sumDeltaTopPx = sumTopPx - totalTopPx_;
    double sumDeltaTopPy = sumTopPy - totalTopPy_;
    double sumDeltaTopPz = sumTopPz - totalTopPz_;
    // setRecoil(sumTopPx, sumTopPy, sumTopPz);
    setRecoil(sumDeltaTopPx, sumDeltaTopPy, sumDeltaTopPz);
}

void topEventMinimizer::findStartingValues(int nPoints)
{
    // cout << "Determining starting values for the ellipse angles and top mass
    // deltas" << endl;

    double startingChi2 = 1.e9;
    double twoPiOverN = 2. * 3.14159265359 / (double)nPoints;
    int whichTop, step;
    double angle;
    vector<double> angles;
    // next 3 lines is my addition
    ellipseAngles_.clear();
    ellipseAngles_.push_back(0.);
    ellipseAngles_.push_back(0.);

    calcWDaughterEllipses();

    // cout << "Beginning loop over all possible angle combinations" << endl;

    for (int i = 0; i < pow(nPoints, nTops_); i++) {
        whichTop = 1;
        angles.clear();
        for (int j = 0; j < nTops_; j++) {
            step = (i / whichTop) % nPoints;
            angle = twoPiOverN * step;
            // cout << "Setting angle for top " << j+1 << " to " << angle <<
            // endl;
            (topSystemChiSquares_.at(j)).first->setEllipseAngle(angle);
            whichTop *= nPoints;
            angles.push_back(angle);
        }
        getDxDyFromEllipses();

        calcHadronicChiSquare(); // W daughter 2 deltas are calculated in
                                 // getDxDyFromEllipses
        double thisChi2 =
            nonTopChiSquare_.getChiSquare() + getHadronicChiSquare();
        // cout << "This chi2 is: " << thisChi2 << endl;
        if (thisChi2 < startingChi2) {
            startingChi2 = thisChi2;
            ellipseAngles_ = angles;
        }
    }

    // cout << "The starting value for the inner chi^2 is " << startingChi2 <<
    // endl;

    // set angles corresponding to minimum chi^2
    for (int iTop = 0; iTop < nTops_; iTop++) {
        // cout<<ellipseAngles_.at(iTop)<<endl;
        (topSystemChiSquares_.at(iTop))
            .first->setEllipseAngle(ellipseAngles_.at(iTop));
    }

    // cout<<"before setup"<<endl;
    setupNonTopChiSquare();
    // cout<<"before calc hadchi"<<endl;
    calcHadronicChiSquare();

    // cout << "Corrected non-top objects:" << endl;
    // nonTopChiSquare_.printResults();
}

double topEventMinimizer::innerMinimizationOperator(const double *inputDeltas)
{
    // cout<<"at innermin operator"<<endl;
    // printTopConstituents();
    vector<double> ellipseAnglesCurrent;
    vector<double> topMassDeltasCurrent;

    int i = 0, iTop = 0;

    for (vector<pair<topSystemChiSquare *, bool>>::const_iterator
             thisTopChiSquare = topSystemChiSquares_.begin();
         thisTopChiSquare != topSystemChiSquares_.end();
         thisTopChiSquare++, iTop++) {
        ellipseAnglesCurrent.push_back(inputDeltas[i]);
        (*thisTopChiSquare).first->setEllipseAngle(inputDeltas[i]);
        i += 1;

        topMassDeltasCurrent.push_back(inputDeltas[i]);
        (*thisTopChiSquare).first->setTopMassDelta(inputDeltas[i]);
        topMassDeltas_.at(iTop) = inputDeltas[i];
        i += 1;

        // if ((*thisTopChiSquare).second == false){
        //    (*thisTopChiSquare).first->printWDaughter2();
        //}
    }

    setupNonTopChiSquare();
    calcHadronicChiSquare();
    calcTopMassChiSquare();
    double innerChi2 = nonTopChi2_ + hadChi2_ + topMassChi2_;

    // cout<< "Innermin"<<endl;
    // cout<< "Innermin hadronic chi2 = " << hadChi2_ <<endl;
    // cout<< "Innermin nontopchi2 = " << nonTopChi2_ <<endl;

    if (innerChi2 < thisInnerChi2Best_) {
        // cout << "I found a new inner chi^2 minimum: " << innerChi2 << endl;

        thisInnerChi2Best_ = innerChi2;
        thisNonTopChi2Best_ = nonTopChi2_;
        thisHadChi2Best_ = hadChi2_;
        thisTopMassChi2Best_ = topMassChi2_;

        ellipseAnglesInnerBest_ = ellipseAnglesCurrent;
        topMassDeltasInnerBest_ = topMassDeltasCurrent;
        // FIXME also save these best non-top objects
    }

    return innerChi2;
}

void topEventMinimizer::minimizeNonTopChiSquare()
{
    // cout << "Doing inner minimization" << endl;

    // cout<<"before infunctor"<<endl;
    // Set up the functor
    ROOT::Math::Functor innerFunc(
        this, &topEventMinimizer::innerMinimizationOperator,
        2 * nTops_); // 1 angle + 1 top mass delta per top system

    // cout<<"before setfunc in"<<endl;
    // Set up the minimization piece:
    innerMin_->SetFunction(innerFunc);

    // cout<<"before find start val"<<endl;
    // Find starting values for the ellipse angles and top mass deltas
    findStartingValues(50);
    // cout<<"after find start val"<<endl;

    // Set up the parameters

    int iPar = 0;
    for (int iTop = 0; iTop < nTops_; iTop++) {
        // cout<<"inside itop loop"<<endl;
        // ellipse angle
        TString parName = "theta_";
        parName += iTop;
        innerMin_->SetVariable(iPar, string(parName), ellipseAngles_[iTop],
                               0.02 * 3.14159265359);
        iPar += 1;

        // top mass delta
        bool hasHighEdge;
        double deltaMTopRangeLow, deltaMTopRangeHigh;
        (topSystemChiSquares_.at(iTop))
            .first->getTopMassDeltaRange(hasHighEdge, deltaMTopRangeLow,
                                         deltaMTopRangeHigh);
        deltaMTopRangeHigh =
            min(deltaMTopRangeHigh, maxConsideredChiSquareRoot_);
        parName = "deltaMTop_";
        parName += iTop;
        if (hasHighEdge) {
            // cout << "Current top mass delta is " << topMassDeltas_.at(iTop)
            // << endl;
            // cout << "deltaMTop range is " << deltaMTopRangeLow << " to " <<
            // deltaMTopRangeHigh << endl;
            innerMin_->SetLimitedVariable(
                iPar, string(parName), topMassDeltas_.at(iTop), 0.1,
                deltaMTopRangeLow, deltaMTopRangeHigh);
        } else {
            // cout << "Current top mass delta is " << topMassDeltas_.at(iTop)
            // << endl;
            // cout << "deltaMTop lower edge is "<< deltaMTopRangeLow << endl;
            innerMin_->SetLowerLimitedVariable(iPar, string(parName),
                                               topMassDeltas_.at(iTop), 0.1,
                                               deltaMTopRangeLow);
        }
        iPar += 1;
    }

    // cout << "Starting the inner minimization" << endl;
    innerMin_->Minimize();

    // cout << "Ending value inner chi^2 reported by Minuit is " <<
    // innerMin_->MinValue() << endl;
    // cout << "Compare to best value I found: " << thisInnerChi2Best_ << endl;

    // Set the best values corresponding to the minimum of the chi square

    iPar = 0;
    for (int iTop = 0; iTop < nTops_; iTop++) {
        // ellipse angle
        ellipseAngles_[iTop] = innerMin_->X()[iPar];
        (topSystemChiSquares_.at(iTop))
            .first->setEllipseAngle(ellipseAnglesInnerBest_.at(iTop));
        iPar += 1;

        // top mass delta
        topMassDeltas_.at(iTop) = innerMin_->X()[iPar];
        (topSystemChiSquares_.at(iTop))
            .first->setTopMassDelta(topMassDeltasInnerBest_.at(iTop));
        iPar += 1;
    }

    // Recalculate chi^2 at the minimum from the best delta values
    setupNonTopChiSquare();
    calcHadronicChiSquare();
    calcTopMassChiSquare();

    // cout << "Inner chi^2 minimum is " << nonTopChi2_+hadChi2_+topMassChi2_ <<
    // endl;
    // cout << "Non-top chi^2 is " << nonTopChi2_ << endl;
    // cout << "Hadronic chi^2 is " << hadChi2_ << endl;
    // cout << "Top mass chi^2 is " << topMassChi2_ << endl;
    // cout << "Best inner chi^2 minimum is " << thisInnerChi2Best_ << endl;
}

double topEventMinimizer::outerMinimizationOperator(const double *inputDeltas)
{
    // std::cout << "at outermin"<<std::endl;
    // printTopConstituents();
    // reset the inner chi^2 minimum for this outer minimizer step
    thisInnerChi2Best_ = 1.e99;
    thisTopMassChi2Best_ = 0.;
    thisNonTopChi2Best_ = 0.;
    thisHadChi2Best_ = 0.;

    int i = 0, iTop = 0;
    for (auto it = topSystemChiSquares_.begin();
         it != topSystemChiSquares_.end(); ++it, ++iTop) {
        // cout << "This is top " << iTop+1 << endl;

        // std::cout<<"inside itop loop"<<std::endl;
        // Set top object deltas and recalculate momenta
        (*it).first->setDeltas(inputDeltas[i], inputDeltas[i + 1],
                               inputDeltas[i + 2], // b-jet deltas
                               inputDeltas[i + 3], inputDeltas[i + 4],
                               inputDeltas[i + 5], // first W daughter deltas
                               inputDeltas[i + 6]  // W mass delta
                               );
        bJets_PtDeltas_.at(iTop) = inputDeltas[i];
        bJets_PhiDeltas_.at(iTop) = inputDeltas[i + 1];
        bJets_EtaDeltas_.at(iTop) = inputDeltas[i + 2];
        firstWDaughters_PtDeltas_.at(iTop) = inputDeltas[i + 3];
        firstWDaughters_PhiDeltas_.at(iTop) = inputDeltas[i + 4];
        firstWDaughters_EtaDeltas_.at(iTop) = inputDeltas[i + 5];
        WMassDeltas_.at(iTop) = inputDeltas[i + 6];

        // std::cout<<"after set stuff"<<std::endl;
        // Calculate the top mass range:
        // done in getTopMassDeltaRange call in inner minimization routine

        // setup the new second W daughter ellipse
        (*it).first->setupWDaughter2Ellipse();

        // std::cout<<"after setupdaughter"<<std::endl;
        // increment i
        i += 7;
    }

    // Calculate the inner piece
    minimizeNonTopChiSquare();

    // std::cout<<"after minnontop"<<std::endl;
    // Calculate the outer piece
    calcTopChiSquare();
    chi2_ = topChi2_ + nonTopChi2_ + hadChi2_ + topMassChi2_;
    // cout<<"topChi2 = "<< topChi2_ <<endl;
    // std::cout<<"after calctopchi"<<std::endl;

    if (chi2_ < chi2Best_) {
        if (debug)
            cout << "New minimum total chi^2: " << chi2_ << endl;

        // update the chi^2 and its components
        chi2Best_ = chi2_;
        topChi2Best_ = topChi2_;
        innerChi2Best_ = thisInnerChi2Best_;
        hadChi2Best_ = thisHadChi2Best_;
        topMassChi2Best_ = thisTopMassChi2Best_;
        nonTopChi2Best_ = thisNonTopChi2Best_;

        // save the current delta values
        bJets_PtDeltasBest_ = bJets_PtDeltas_;
        bJets_PhiDeltasBest_ = bJets_PhiDeltas_;
        bJets_EtaDeltasBest_ = bJets_EtaDeltas_;
        firstWDaughters_PtDeltasBest_ = firstWDaughters_PtDeltas_;
        firstWDaughters_PhiDeltasBest_ = firstWDaughters_PhiDeltas_;
        firstWDaughters_EtaDeltasBest_ = firstWDaughters_EtaDeltas_;
        WMassDeltasBest_ = WMassDeltas_;
        topMassDeltasBest_ = topMassDeltasInnerBest_;
        ellipseAnglesBest_ = ellipseAnglesInnerBest_;
    }

    return chi2_;
}

void topEventMinimizer::minimizeTotalChiSquare()
{
    //     std::cout<<"at min"<<std::endl;
    const int nParameters = 7 * nTops_; // 3 per b-jet + 3 per W daughter 1 + 1
                                        // per W mass = 7 per top system

    //     std::cout<<"before set functor"<<std::endl;
    // Set up the functor
    ROOT::Math::Functor func(
        this, &topEventMinimizer::outerMinimizationOperator, nParameters);

    //     std::cout<<"before setfunc"<<std::endl;
    // Set up the minimization piece:
    outerMin_->SetFunction(func);

    //     std::cout<<"before set min param"<<std::endl;
    // Setup the minimizer parameters

    int iPar = 0, iTop = 0;
    const double max = maxConsideredChiSquareRoot_;
    for (auto it = topSystemChiSquares_.begin();
         it != topSystemChiSquares_.end(); ++it, ++iTop) {
        ostringstream convert; // stream used for the (int) conversion
        convert << iTop;
        const string iTop_str = convert.str();
        const string par1 = "bJetPtDelta_" + iTop_str;
        outerMin_->SetLimitedVariable(iPar, par1, bJets_PtDeltas_.at(iTop), 0.1,
                                      -max, max);
        ++iPar;
        const string par2 = "bJetPhiDelta_" + iTop_str;
        outerMin_->SetLimitedVariable(iPar, par2, bJets_PhiDeltas_.at(iTop),
                                      0.1, -max, max);
        ++iPar;
        const string par3 = "bJetEtaDelta_" + iTop_str;
        outerMin_->SetLimitedVariable(iPar, par3, bJets_EtaDeltas_.at(iTop),
                                      0.1, -max, max);
        ++iPar;
        const string par4 = "WDaughter1PtDelta_" + iTop_str;
        outerMin_->SetLimitedVariable(
            iPar, par4, firstWDaughters_PtDeltas_.at(iTop), 0.1, -max, max);
        ++iPar;
        const string par5 = "WDaughter1PhiDelta_" + iTop_str;
        outerMin_->SetLimitedVariable(
            iPar, par5, firstWDaughters_PhiDeltas_.at(iTop), 0.1, -max, max);
        ++iPar;
        const string par6 = "WDaughter1EtaDelta_" + iTop_str;
        outerMin_->SetLimitedVariable(
            iPar, par6, firstWDaughters_EtaDeltas_.at(iTop), 0.1, -max, max);
        ++iPar;
        const string par7 = "deltaMW_" + iTop_str;
        outerMin_->SetLimitedVariable(iPar, par7, WMassDeltas_.at(iTop), 0.1,
                                      -max, max);
        ++iPar;
    }

    // std::cout<<"before minimize"<<std::endl;
    // cout << "Starting outer minimization" << endl;
    outerMin_->Minimize();

    // std::cout<<"after minimise"<<std::endl;

    // cout << "Minimum chi^2 values reported by Minuit:" << endl;
    // cout << "Total chi^2 is " << outerMin_->MinValue() << endl;
    // cout << "Inner chi^2 is " << innerMin_->MinValue() << endl;

    // cout << "Minimum values I found:" << endl;
    // cout << "Best total chi^2 is " << chi2Best_ << endl;
    // cout << "Best inner chi^2 is " << innerChi2Best_ << endl;
    // cout << "Final outer chi^2 check: " <<
    // innerChi2Best_+topChi2Best_-chi2Best_ << endl;
    // cout << "Final inner chi^2 check: " <<
    // hadChi2Best_+topMassChi2Best_+nonTopChi2Best_-innerChi2Best_ << endl;
    // cout << "Best hadronic W daughter chi^2 is " << hadChi2Best_ << endl;
    // cout << "Best top mass chi^2 is " << topMassChi2Best_ << endl;
    // cout << "Best non-top chi^2 is " << nonTopChi2Best_ << endl;
    // cout << "Best top system chi^2 is " << topChi2Best_ << endl;

    // cout << "Printing outer min results" << endl;
    outerMin_->SetPrintLevel(1);
    outerMin_->PrintResults();
    cout << "Outer min status is " << outerMin_->Status() << endl;

    // cout << "Printing inner min results" << endl;
    innerMin_->SetPrintLevel(4);
    innerMin_->PrintResults();
    cout << "Inner min status is " << innerMin_->Status() << endl;

    outerMinStatus = outerMin_->Status();
    innerMinStatus = innerMin_->Status();
    outerMin_Edm = outerMin_->Edm();

    setBestValues();
}

void topEventMinimizer::setBestValues()
{
    // cout << "Setting delta values corresponding to the minimum total chi^2 I
    // found: " << chi2Best_ << endl;

    int iTop = 0;
    for (vector<pair<topSystemChiSquare *, bool>>::const_iterator
             thisTopChiSquare = topSystemChiSquares_.begin();
         thisTopChiSquare != topSystemChiSquares_.end();
         thisTopChiSquare++, iTop++) {
        // cout << "This is top number " << iTop+1 << endl;
        // cout << "b-jet pt delta  = " << bJets_PtDeltasBest_ .at(iTop) <<
        // endl;
        // cout << "b-jet phi delta = " << bJets_PhiDeltasBest_.at(iTop) <<
        // endl;
        // cout << "b-jet eta delta = " << bJets_EtaDeltasBest_.at(iTop) <<
        // endl;
        // cout << "first W daughter pt delta  = " <<
        // firstWDaughters_PtDeltasBest_ .at(iTop) << endl;
        // cout << "first W daughter phi delta = " <<
        // firstWDaughters_PhiDeltasBest_.at(iTop) << endl;
        // cout << "first W daughter eta delta = " <<
        // firstWDaughters_EtaDeltasBest_.at(iTop) << endl;
        // cout << "W   mass delta = " << WMassDeltasBest_  .at(iTop) << endl;
        // cout << "Top mass delta = " << topMassDeltasBest_.at(iTop) << endl;
        // cout << "Second W daughter ellipse angle = " <<
        // ellipseAnglesBest_.at(iTop) << endl;

        (*thisTopChiSquare).first->setTopMassDelta(topMassDeltasBest_.at(iTop));
        (*thisTopChiSquare)
            .first->setDeltas(bJets_PtDeltasBest_.at(iTop),
                              bJets_PhiDeltasBest_.at(iTop),
                              bJets_EtaDeltasBest_.at(iTop),
                              firstWDaughters_PtDeltasBest_.at(iTop),
                              firstWDaughters_PhiDeltasBest_.at(iTop),
                              firstWDaughters_EtaDeltasBest_.at(iTop),
                              WMassDeltasBest_.at(iTop));

        (*thisTopChiSquare).first->setupWDaughter2Ellipse();
        (*thisTopChiSquare).first->calcWDaughter2Ellipse();
        (*thisTopChiSquare).first->setEllipseAngle(ellipseAnglesBest_.at(iTop));

        double ptDelta, phiDelta, etaDelta;
        (*thisTopChiSquare)
            .first->getWDaughter2Deltas(ptDelta, phiDelta, etaDelta);
        secondWDaughters_PtDeltasBest_.at(iTop) = ptDelta;
        secondWDaughters_PhiDeltasBest_.at(iTop) = phiDelta;
        secondWDaughters_EtaDeltasBest_.at(iTop) = etaDelta;
        // cout << "second W daugther pt  delta is " << ptDelta  << endl;
        // cout << "second W daugther eta delta is " << etaDelta << endl;
        // cout << "second W daugther phi delta is " << phiDelta << endl;
    }

    setupNonTopChiSquare();
    calcHadronicChiSquare();
    calcTopMassChiSquare();
    calcTopChiSquare();
    calcTotalChiSquare();
    buildBestNonTopObjects();

    // cout << "Total chi^2 check: " << chi2_-chi2Best_ << endl;
    // cout << "Inner chi^2 check: " <<
    // topMassChi2_+nonTopChi2_+hadChi2_-innerChi2Best_ << endl;
    // cout << "Top system chi^2 check: " << topChi2_-topChi2Best_ << endl;
    // cout << "Hadronic W daughter chi^2 check: " << hadChi2_-hadChi2Best_ <<
    // endl;
    // cout << "Top mass chi^2 check: " << topMassChi2_-topMassChi2Best_ <<
    // endl;
    // cout << "Non-top chi^2 check: " << nonTopChi2_-nonTopChi2Best_ << endl;
}

void topEventMinimizer::getBestDeltas(
    vector<double> &bJetPtDeltas, vector<double> &bJetPhiDeltas,
    vector<double> &bJetEtaDeltas, vector<double> &firstWDaughterPtDeltas,
    vector<double> &firstWDaughterPhiDeltas,
    vector<double> &firstWDaughterEtaDeltas,
    vector<double> &secondWDaughterPtDeltas,
    vector<double> &secondWDaughterPhiDeltas,
    vector<double> &secondWDaughterEtaDeltas, vector<double> &WMassDeltas,
    vector<double> &topMassDeltas, vector<double> &nonTopObjectPxDeltas,
    vector<double> &nonTopObjectPyDeltas, vector<double> &nonTopObjectPzDeltas)
{
    bJetPtDeltas.clear();
    bJetPhiDeltas.clear();
    bJetEtaDeltas.clear();

    firstWDaughterPtDeltas.clear();
    firstWDaughterPhiDeltas.clear();
    firstWDaughterEtaDeltas.clear();

    secondWDaughterPtDeltas.clear();
    secondWDaughterPhiDeltas.clear();
    secondWDaughterEtaDeltas.clear();

    WMassDeltas.clear();

    topMassDeltas.clear();

    nonTopObjectPxDeltas.clear();
    nonTopObjectPyDeltas.clear();
    nonTopObjectPzDeltas.clear();

    // cout << "Filling top object deltas" << endl;
    for (int iTop = 0; iTop < nTops_; iTop++) {
        bJetPtDeltas.push_back(bJets_PtDeltasBest_.at(iTop));
        bJetPhiDeltas.push_back(bJets_PhiDeltasBest_.at(iTop));
        bJetEtaDeltas.push_back(bJets_EtaDeltasBest_.at(iTop));

        firstWDaughterPtDeltas.push_back(
            firstWDaughters_PtDeltasBest_.at(iTop));
        firstWDaughterPhiDeltas.push_back(
            firstWDaughters_PhiDeltasBest_.at(iTop));
        firstWDaughterEtaDeltas.push_back(
            firstWDaughters_EtaDeltasBest_.at(iTop));

        secondWDaughterPtDeltas.push_back(
            secondWDaughters_PtDeltasBest_.at(iTop));
        secondWDaughterPhiDeltas.push_back(
            secondWDaughters_PhiDeltasBest_.at(iTop));
        secondWDaughterEtaDeltas.push_back(
            secondWDaughters_EtaDeltasBest_.at(iTop));

        WMassDeltas.push_back(WMassDeltasBest_.at(iTop));

        topMassDeltas.push_back(topMassDeltasBest_.at(iTop));
    }

    const unsigned int sz = nonTopObjects_.size();
    nonTopObjectPxDeltas.reserve(sz);
    nonTopObjectPyDeltas.reserve(sz);
    nonTopObjectPzDeltas.reserve(sz);
    // cout << "Now filling non-top object deltas" << endl;
    for (unsigned int i = 0; i < sz; ++i) {
        nonTopObjectPxDeltas.push_back(nonTopObjects_PxDeltasBest_[i]);
        nonTopObjectPyDeltas.push_back(nonTopObjects_PyDeltasBest_[i]);
        nonTopObjectPzDeltas.push_back(nonTopObjects_PyDeltasBest_[i]);
    }
}

double topEventMinimizer::getOneTopMassChiSquare(int iTop)
{
    return (topSystemChiSquares_.at(iTop)).first->getTopMassChiSquare();
}

double topEventMinimizer::getOneBChiSquare(int iTop)
{
    // cout<<"inside getonebchiSquare"<<endl;
    double toreturn = (topSystemChiSquares_.at(iTop)).first->getBChiSquare();
    return toreturn;
}

double topEventMinimizer::getOneWDaughter1ChiSquare(int iTop)
{
    double toreturn =
        (topSystemChiSquares_.at(iTop)).first->getWDaughter1ChiSquare();
    return toreturn;
}

double topEventMinimizer::getOneWMassChiSquare(int iTop)
{
    double toreturn =
        (topSystemChiSquares_.at(iTop)).first->getWMassChiSquare();
    return toreturn;
}

void topEventMinimizer::calcTopChiSquare()
{
    topChi2_ = 0;
    for (vector<pair<topSystemChiSquare *, bool>>::const_iterator
             thisTopChiSquare = topSystemChiSquares_.begin();
         thisTopChiSquare != topSystemChiSquares_.end(); thisTopChiSquare++) {
        topChi2_ += (*thisTopChiSquare).first->getChiSquare();
    }
}

void topEventMinimizer::calcHadronicChiSquare()
{
    hadChi2_ = 0;
    for (vector<pair<topSystemChiSquare *, bool>>::const_iterator
             thisTopChiSquare = topSystemChiSquares_.begin();
         thisTopChiSquare != topSystemChiSquares_.end(); thisTopChiSquare++) {
        if (!((*thisTopChiSquare).second))
            hadChi2_ += (*thisTopChiSquare).first->getHadronicChiSquare();
    }
}

void topEventMinimizer::calcTopMassChiSquare()
{
    topMassChi2_ = 0;
    for (vector<pair<topSystemChiSquare *, bool>>::const_iterator
             thisTopChiSquare = topSystemChiSquares_.begin();
         thisTopChiSquare != topSystemChiSquares_.end(); thisTopChiSquare++) {
        topMassChi2_ += (*thisTopChiSquare).first->getTopMassChiSquare();
    }
}

void topEventMinimizer::calcTotalChiSquare()
{
    chi2_ = nonTopChiSquare_.getChiSquare() + getTopChiSquare() +
            getHadronicChiSquare() + getTopMassChiSquare();
}

double topEventMinimizer::getChiSquare()
{
    calcTotalChiSquare();
    // cout << "Current chi^2 is " << chi2_ << endl;
    return chi2_;
}

double topEventMinimizer::getTopChiSquare()
{
    calcTopChiSquare();
    return topChi2_;
}

double topEventMinimizer::getTopMassChiSquare()
{
    calcTopMassChiSquare();
    return topMassChi2_;
}

double topEventMinimizer::getHadronicChiSquare()
{
    calcHadronicChiSquare();
    return hadChi2_;
}

double topEventMinimizer::getNonTopChiSquare()
{
    nonTopChi2_ = nonTopChiSquare_.getChiSquare();
    return nonTopChi2_;
}

void topEventMinimizer::getBJet(int whichTop, double &px, double &py,
                                double &pz, double &e)
{
    if (whichTop >= nTops_) {
        px = 0;
        py = 0;
        pz = 0;
        e = 0;
        return;
    }

    (topSystemChiSquares_.at(whichTop)).first->getBJet(px, py, pz, e);
}

void topEventMinimizer::getWDaughter1(int whichTop, double &px, double &py,
                                      double &pz, double &e)
{
    if (whichTop >= nTops_) {
        px = 0;
        py = 0;
        pz = 0;
        e = 0;
        return;
    }

    (topSystemChiSquares_.at(whichTop)).first->getWDaughter1(px, py, pz, e);
}

void topEventMinimizer::getWDaughter2(int whichTop, double &px, double &py,
                                      double &pz, double &e)
{
    if (whichTop >= nTops_) {
        px = 0;
        py = 0;
        pz = 0;
        e = 0;
        return;
    }

    (topSystemChiSquares_.at(whichTop)).first->getWDaughter2(px, py, pz, e);
}

XYZTLorentzVector topEventMinimizer::getConverter(string whichFunc,
                                                  int whichTop)
{
    /* typedef void (*FnPtr)(int, double, double, double, double);
     map < string, FnPtr > funcMap;
     funcMap["getTop"] = getTop;
     funcMap["getW"] = getW;
     funcMap["getNonTopObject"] = getNonTopObject;
     funcMap["getBJet"] = getBJet;
     funcMap["getWDaughter1"] = getWDaughter1;
     funcMap["getWDaughter2"] = getWDaughter2;*/

    double px, py, pz, e;
    px = 0;
    py = 0;
    pz = 0;
    e = 0;

    // typedef map<string, TLorentzVector> tmap;
    // for (tmap::iterator h = funcMap.begin(); h != funcMap.end(); h++){
    //    std::size_t foundString = (h->first).find( whichFunc );
    //    if (foundString != string::npos){

    // FnPtr f = funcMap[whichFunc];
    //(*f)(whichTop, px, py, pz, e);

    // if (whichFunc == "getNonTopObject"){

    //}

    (this->*funcMap[whichFunc])(whichTop, px, py, pz, e);
    XYZTLorentzVector toreturn;
    toreturn.SetPxPyPzE(px, py, pz, e);

    return toreturn;
}

void topEventMinimizer::getTop(int whichTop, double &px, double &py, double &pz,
                               double &e)
{
    if (whichTop >= nTops_) {
        px = 0;
        py = 0;
        pz = 0;
        e = 0;
        return;
    }

    (topSystemChiSquares_.at(whichTop)).first->getTopMomentum(px, py, pz, e);
}

void topEventMinimizer::getW(int whichTop, double &px, double &py, double &pz,
                             double &e)
{
    if (whichTop >= nTops_) {
        px = 0;
        py = 0;
        pz = 0;
        e = 0;
        return;
    }

    double px1, py1, pz1, e1;
    double px2, py2, pz2, e2;

    getWDaughter1(whichTop, px1, py1, pz1, e1);
    getWDaughter2(whichTop, px2, py2, pz2, e2);

    px = px1 + px2;
    py = py1 + py2;
    pz = pz1 + pz2;
    e = e1 + e2;
}

void topEventMinimizer::getNonTopObject(int whichObject, double &px, double &py,
                                        double &pz)
{
    if (whichObject >= (int)nonTopObjects_.size()) {
        px = 0;
        py = 0;
        pz = 0;
        return;
    }

    px = nonTopObjects_[whichObject].Px() +
         nonTopObjects_PxDeltasBest_[whichObject];
    py = nonTopObjects_[whichObject].Py() +
         nonTopObjects_PyDeltasBest_[whichObject];
    pz = nonTopObjects_[whichObject].Pz() +
         nonTopObjects_PzDeltasBest_[whichObject];
}

void topEventMinimizer::getNonTopObject4(int whichObject, double &px,
                                         double &py, double &pz, double &e)
{
    if (whichObject >= (int)nonTopObjects_.size()) {
        px = 0;
        py = 0;
        pz = 0;
        e = 0;
        return;
    }

    px = nonTopObjects_[whichObject].Px() +
         nonTopObjects_PxDeltasBest_[whichObject];
    py = nonTopObjects_[whichObject].Py() +
         nonTopObjects_PyDeltasBest_[whichObject];
    pz = nonTopObjects_[whichObject].Pz() +
         nonTopObjects_PzDeltasBest_[whichObject];
    e = sqrt(pow(nonTopObjects_[whichObject].M(), 2) + pow(px, 2) + pow(py, 2) +
             pow(pz, 2));
}

void topEventMinimizer::plotEllipses(TString plotName)
{
    // cout << "Plotting the second W daughter ellipses" << endl;
    // calcWDaughterEllipses();

    int nPoints(2000);
    // int nPoints(5);
    double twoPiOverN = 2. * 3.14159265359 / (double)nPoints;

    vector<TGraph> ellipses = vector<TGraph>(nTops_, TGraph(nPoints + 1));
    TGraph points(nTops_);

    TMatrixD *thisEllipse;
    double thisTheta = 0;
    double thetaArray[3] = {};
    TVectorD thisWDaughter2Perp(3, thetaArray);
    double thisWDaughter2X, thisWDaughter2Y;
    double sumPx(nonTopPx_), sumPy(nonTopPy_), sumPz(nonTopPz_);

    double maxX(0), minX(1.e9);
    double maxY(0), minY(1.e9);

    int iTop = 0;

    for (vector<pair<topSystemChiSquare *, bool>>::const_iterator
             thisTopChiSquare = topSystemChiSquares_.begin();
         thisTopChiSquare != topSystemChiSquares_.end();
         thisTopChiSquare++, iTop++) {
        // cout << "Top number " << iTop+1 << endl;

        thisEllipse =
            (*thisTopChiSquare).first->getHomogeneousWDaughterEllipse();
        // cout << "Hperp is: " << endl;
        // thisEllipse->Print();

        // reset angles, starting with theta=0
        // cout << "resetting angles" << endl;
        thisTheta = 0;
        // cout << "Ellipse angle: " << thisTheta << endl;
        (*thisTopChiSquare).first->setEllipseAngle(thisTheta);
        thetaArray[0] = 1;
        thetaArray[1] = 0;
        thetaArray[2] = 1;

        // second W daughter
        thisWDaughter2Perp = TVectorD(3, thetaArray);
        // thisWDaughter2Perp.Print();
        thisWDaughter2Perp *= *thisEllipse;
        // thisWDaughter2Perp.Print();

        double thisWDaughter1Px, thisWDaughter1Py, thisWDaughter1Pz,
            thisWDaughter1E;
        (*thisTopChiSquare)
            .first->getWDaughter1(thisWDaughter1Px, thisWDaughter1Py,
                                  thisWDaughter1Pz, thisWDaughter1E);
        sumPx += thisWDaughter1Px;
        sumPy += thisWDaughter1Py;
        sumPz += thisWDaughter1Pz;

        double thisBJetPx, thisBJetPy, thisBJetPz, thisBJetE;
        (*thisTopChiSquare)
            .first->getBJet(thisBJetPx, thisBJetPy, thisBJetPz, thisBJetE);
        sumPx += thisBJetPx;
        sumPy += thisBJetPy;
        sumPz += thisBJetPz;

        // cout << "Current sum px = " << sumPx << endl;
        // cout << "Current sum py = " << sumPy << endl;

        if (iTop < nTops_ - 1) {
            // cout << "Not the last top" << endl;
            thisWDaughter2X = thisWDaughter2Perp[0];
            thisWDaughter2Y = thisWDaughter2Perp[1];
        } else {
            // cout << "Last top" << endl;
            // cout << "sum px = " << sumPx << endl;
            // cout << "sum py = " << sumPy << endl;
            thisWDaughter2X = -sumPx - thisWDaughter2Perp[0];
            thisWDaughter2Y = -sumPy - thisWDaughter2Perp[1];
        }
        // thisWDaughter2X=thisWDaughter2Perp[0];
        // thisWDaughter2Y=thisWDaughter2Perp[1];

        ellipses[iTop].SetPoint(0, thisWDaughter2X, thisWDaughter2Y);

        maxX = max(maxX, thisWDaughter2X);
        minX = min(minX, thisWDaughter2X);
        maxY = max(maxY, thisWDaughter2Y);
        minY = min(minY, thisWDaughter2Y);

        // loop over points around the ellipse

        for (int iPoint = 1; iPoint <= nPoints; iPoint++) {
            // cout << "At point " << iPoint << endl;
            thisTheta = (double)iPoint * twoPiOverN;
            // cout << "Ellipse angle: " << thisTheta << endl;
            (*thisTopChiSquare).first->setEllipseAngle(thisTheta);
            thisWDaughter2Perp[0] = cos(thisTheta);
            thisWDaughter2Perp[1] = sin(thisTheta);
            thisWDaughter2Perp[2] = 1;

            thisWDaughter2Perp *= *thisEllipse;
            // thisWDaughter2Perp.Print();

            if (iTop < nTops_ - 1) {
                thisWDaughter2X = thisWDaughter2Perp[0];
                thisWDaughter2Y = thisWDaughter2Perp[1];
            } else {
                thisWDaughter2X = -sumPx - thisWDaughter2Perp[0];
                thisWDaughter2Y = -sumPy - thisWDaughter2Perp[1];
            }
            // thisWDaughter2X=thisWDaughter2Perp[0];
            // thisWDaughter2Y=thisWDaughter2Perp[1];

            ellipses[iTop].SetPoint(iPoint, thisWDaughter2X, thisWDaughter2Y);

            maxX = max(maxX, thisWDaughter2X);
            minX = min(minX, thisWDaughter2X);
            maxY = max(maxY, thisWDaughter2Y);
            minY = min(minY, thisWDaughter2Y);

            // cout << "Current plot bounds: " << minX << "\t" << maxX << "\t"
            // << minY << "\t" << maxY << endl;

        } // end loop over points

        if (iTop == 0)
            ellipses[iTop].SetLineColor(kRed);
        else if (iTop == 1)
            ellipses[iTop].SetLineColor(kBlue);
        else
            ellipses[iTop].SetLineColor(kBlack);

        // draw second W daughter momentum at the minimum chi^2
        // cout << "Ellipse angle at the minimum chi^2 for top " << iTop << " is
        // " << ellipseAngles_.at(iTop) << endl;
        thisWDaughter2Perp[0] = cos(ellipseAngles_.at(iTop));
        thisWDaughter2Perp[1] = sin(ellipseAngles_.at(iTop));
        thisWDaughter2Perp[2] = 1;

        // thisEllipse->Print();

        thisWDaughter2Perp *= *thisEllipse;

        if (iTop < nTops_ - 1) {
            points.SetPoint(iTop, thisWDaughter2Perp[0], thisWDaughter2Perp[1]);
        } else {
            points.SetPoint(iTop, -sumPx - thisWDaughter2Perp[0],
                            -sumPy - thisWDaughter2Perp[1]);
        }
        // points.SetPoint(iTop,thisWDaughter2Perp[0],thisWDaughter2Perp[1]);

    } // end loop over tops

    // cout << "End loop over tops" << endl;
    // cout << "sum px = " << sumPx << endl;
    // cout << "sum py = " << sumPy << endl;

    points.SetMarkerStyle(24);
    points.SetMarkerSize(3);
    // points.Print();

    // now plot the ellipses
    TH1D drawer("drawer", "drawer", 1, minX - 0.05 * (maxX - minX),
                maxX + 0.05 * (maxX - minX));
    TCanvas canv(plotName, plotName, 800, 800);
    drawer.Draw("");
    drawer.SetAxisRange(minY - 0.05 * (maxY - minY),
                        maxY + 0.05 * (maxY - minY), "Y");
    drawer.Draw("AXIS");
    points.Draw("PSAME");
    for (int i = 0; i < nTops_; i++) {
        // cout << "Top number " << i+1 << endl;
        // cout << "Number of points is " << ellipses[i].GetN() << endl;
        ellipses.at(i).Draw("LSAME");
        ellipses.at(i).Clear();
    }
    canv.SaveAs(plotName + ".pdf");
    canv.Clear();
    points.Clear();
    drawer.Clear();
}
