#ifndef TOPEVENTMINIMIZER
#define TOPEVENTMINIMIZER

#include "neutrinoSolutions.h"
#include "topSystemChiSquare.h"
#include "hadronicTopSystemChiSquare.h"
#include "leptonicTopSystemChiSquare.h"
#include "lightJetChiSquareMinimumSolver.h"
#include "TLorentzVector.h"
#include <map>

class topEventMinimizer
{
  private:
    int nTops_;

    vector<int> bJets_;
    vector<int> firstWDaughters_;
    vector<int> secondWDaughters_;

    vector<XYZTLorentzVector> allObjects_;
    vector<double> allObjectPtWidths_;
    vector<double> allObjectPhiWidths_;
    vector<double> allObjectEtaWidths_;

    vector<bool> isLeptonicTopDecay_;

    // vector<XYZTLorentzVector> bJetLorentzVectors_;

    // vector<double> bJetPxs_, bJetPys_, bJetPzs_, bJetEs_;

    // vector<double> bJetPtWidths_, bJetEtaWidths_, bJetPhiWidths_;

    // vector<pair<XYZTLorentzVector, XYZTLorentzVector> >
    // WDaughterLorentzVectors_;

    // vector<pair<double, double> > WDaughterPxs_, WDaughterPys_,
    // WDaughterPzs_, WDaughterEs_;

    // vector<pair<double, double> > WDaughterPtWidths_, WDaughterEtaWidths_,
    // WDaughterPhiWidths_;

    vector<XYZTLorentzVector> nonTopObjects_;
    vector<double> nonTopObjectPts_;
    vector<double> nonTopObjectPhis_;
    vector<double> nonTopObjectPtWidths_;
    vector<double> nonTopObjectPhiWidths_;
    vector<double> nonTopObjectEtaWidths_;

    double nonTopPx_;
    double nonTopPy_;
    double nonTopPz_;

    double mTop_;
    double sigmaMTop_;

    double mW_;
    double sigmaMW_;

    double totalTopPx_;
    double totalTopPy_;
    double totalTopPz_;

    vector<double> bJets_PtDeltas_;
    vector<double> bJets_PhiDeltas_;
    vector<double> bJets_EtaDeltas_;
    vector<double> firstWDaughters_PtDeltas_;
    vector<double> firstWDaughters_PhiDeltas_;
    vector<double> firstWDaughters_EtaDeltas_;
    vector<double> topMassDeltas_;
    vector<double> WMassDeltas_;

    vector<double> bJets_PtDeltasBest_;
    vector<double> bJets_PhiDeltasBest_;
    vector<double> bJets_EtaDeltasBest_;
    vector<double> firstWDaughters_PtDeltasBest_;
    vector<double> firstWDaughters_PhiDeltasBest_;
    vector<double> firstWDaughters_EtaDeltasBest_;
    vector<double> topMassDeltasBest_;
    vector<double> WMassDeltasBest_;

    vector<double> secondWDaughters_PtDeltasBest_;
    vector<double> secondWDaughters_PhiDeltasBest_;
    vector<double> secondWDaughters_EtaDeltasBest_;

    vector<double> topMassDeltasCurrent_;
    vector<double> topMassDeltasInnerBest_;

    vector<double> nonTopObjects_PxDeltasBest_;
    vector<double> nonTopObjects_PyDeltasBest_;
    vector<double> nonTopObjects_PzDeltasBest_;

    // vector<pair<double, double> > WDaughterMasses_;

    double chi2_;

    vector<pair<topSystemChiSquare *, bool>> topSystemChiSquares_;

    double nonTopChi2_;
    double dx_, dy_, dz_;

    lightJetChiSquareMinimumSolver nonTopChiSquare_;

    vector<double> ellipseAngles_;
    vector<double> ellipseAnglesBest_;
    vector<double> ellipseAnglesCurrent_;
    vector<double> ellipseAnglesInnerBest_;

    double hadChi2_;
    double topChi2_;
    double topMassChi2_;

    double chi2Best_;
    double innerChi2Best_;
    double topChi2Best_;
    double hadChi2Best_;
    double topMassChi2Best_;
    double nonTopChi2Best_;

    double thisInnerChi2Best_;
    double thisTopMassChi2Best_;
    double thisNonTopChi2Best_;
    double thisHadChi2Best_;

    double maxConsideredChiSquareRoot_;

    // ROOT::Math::Minimizer* ellipseAngleMin_;
    // ROOT::Math::Minimizer* topMassMin_;
    ROOT::Math::Minimizer *innerMin_;
    ROOT::Math::Minimizer *outerMin_;

    bool checkInputSizes();

    // void setBJets();
    void setNonTopObjectCollections();
    // void setWDaughters();

    //  void initializeDeltas();
    void initializeChiSquares();
    void Initialize_minimizers(ROOT::Math::Minimizer *&outer,
                               ROOT::Math::Minimizer *&inner);

    void setRecoil(double, double, double);

    void calcWDaughterEllipses();

    void getDxDyFromEllipses();
    void calcNonTopMomentum();

    void buildBestNonTopObjects();

    void setupNonTopChiSquare();

    void calcHadronicChiSquare();
    void calcTopMassChiSquare();
    void calcTopChiSquare();

    /*  double getHadronicChiSquare();
      double getTopChiSquare();
      double getTopMassChiSquare();
      double getNonTopChiSquare();
    */
    void setBestValues();

    void setupMap();
    void testfunc();

    typedef void (topEventMinimizer::*FnPtr)(int, double &, double &, double &,
                                             double &);
    map<string, FnPtr> funcMap;

    // typedef void (topEventMinimizer::FnTestPtr)();
    // map < string, FnTestPtr* > testMap;

  public:
    // flags
    const bool debug = true;

    // status observables
    int innerMinStatus;
    int outerMinStatus;

    double outerMin_Edm;

    // functions
    topEventMinimizer(vector<XYZTLorentzVector>, vector<double>, vector<double>,
                      vector<double>, vector<int>, vector<int>, vector<int>,
                      vector<bool>, double, double, double, double, double,
                      double, double);

    topEventMinimizer(vector<XYZTLorentzVector>, vector<double>, vector<double>,
                      vector<double>, double, double, double, double, double,
                      double, double);

    ~topEventMinimizer();

    topSystemChiSquare *makeLeptonicTop(int, int);
    topSystemChiSquare *makeHadronicTop(int, int, int);

    topSystemChiSquare *makeLeptonicTop(double, double, double, double, double,
                                        double, double, double, double, double,
                                        double, double, double, double, double,
                                        double, double, double);
    topSystemChiSquare *makeHadronicTop(double, double, double, double, double,
                                        double, double, double, double, double,
                                        double, double, double, double, double,
                                        double, double, double, double, double,
                                        double, double, double, double, double);

    void addLeptonicTop(int, int);
    void addHadronicTop(int, int, int);

    void addLeptonicTop(double, double, double, double, double, double, double,
                        double, double, double, double, double, double, double,
                        double, double, double, double);
    void addHadronicTop(double, double, double, double, double, double, double,
                        double, double, double, double, double, double, double,
                        double, double, double, double, double, double, double,
                        double, double, double, double);

    void printTopConstituents();
    void printNonTopObjects();

    void calcTopMassRanges();

    // double ellipseAngleMinimizationOperator(const double* );
    double outerMinimizationOperator(const double *);
    double innerMinimizationOperator(const double *);

    void findStartingValues(int);
    void minimizeNonTopChiSquare();
    void minimizeTotalChiSquare();

    void calcTotalChiSquare();
    double getChiSquare();

    double getBestTotalChiSquare() { return chi2Best_; };
    double getBestTopSystemChiSquare() { return topChi2Best_; };
    double getBestHadronicChiSquare() { return hadChi2Best_; };
    double getBestTopMassChiSquare() { return topMassChi2Best_; };
    double getBestNonTopChiSquare() { return nonTopChi2Best_; };

    void getBestDeltas(vector<double> &, vector<double> &, vector<double> &,
                       vector<double> &, vector<double> &, vector<double> &,
                       vector<double> &, vector<double> &, vector<double> &,
                       vector<double> &, vector<double> &, vector<double> &,
                       vector<double> &, vector<double> &);

    void getBJet(int, double &, double &, double &, double &);
    void getWDaughter1(int, double &, double &, double &, double &);
    void getWDaughter2(int, double &, double &, double &, double &);
    void getTop(int, double &, double &, double &, double &);
    void getW(int, double &, double &, double &, double &);
    void getNonTopObject(int, double &, double &, double &);
    void getNonTopObject4(int, double &, double &, double &, double &);

    XYZTLorentzVector getConverter(string, int);

    void plotEllipses(TString);

    double getHadronicChiSquare();
    double getTopChiSquare();
    double getTopMassChiSquare();
    double getNonTopChiSquare();

    double getOneTopMassChiSquare(int);
    double getOneBChiSquare(int);
    double getOneWDaughter1ChiSquare(int);
    double getOneWMassChiSquare(int);

    void initializeDeltas();
};

#endif
