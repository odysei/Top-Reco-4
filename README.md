#Top quark reconstruction

This repository contains a C++ implementation of the methods described in "Object reconstruction in collider events containing top quarks".


## Instructions to run

The code can be run from the command line:

root -b -l -q loadFiles.C'(<output_dir>,whichLoop,maxLoops)'

where <output_dir> is a string containing the path to the directory where the output is written; whichLoop and maxLoops are integers used to split processing into multiple jobs. For each job to run on 100 events, maxLoops=100 and whichLoop will be the starting event number. For example, to run on the first 100 events and write the output to the directory output_files, call:

root -b -l -q loadFiles.C'("output_files",0,100)'



## Class and file descriptions


### WDaughterEllipseCalculator

A C++ implementation of the algorithms presented in [Betchart et al, arXiv:1305.1878](http://arxiv.org/abs/1305.1878) (in Python). 
The framework has been extended to consider either leptonic or hadronic top decays.  In practice this does not introduce any changes to the class.


#### Constructor

The following information must be input to the constructor:

- four-momentum of the b-jet
- four-momentum of one measured W daughter
- assumed top mass
- assumed W mass
- assumed second W daughter mass

Coordinates of the b-jet and first W daughter and the top and W masses are passed by reference in order to allow for corrections to the four-momenta, as needed in the top system chi square class. FIXME vary the second W daughter mass?

#### Methods of interest

- setupEllipse(double mTop, double mW, double WDaughter2Mass)
Set the assumed top, W, and second W daughter masses and recalculate the different variables defined in the Betchart paper.

- calcWDaughterEllipse()
Calculate the second W daughter ellipse in the transverse plane px, py (Hperp).

- calcExtendedWDaughterEllipse()
Calculate the second W daughter ellipse in the extended representation px, py, pz (Nperp).

- getExtendedWDaughterEllipse()
Return Hperp.

- getHomogeneousWDaughterEllipse()
Return Nperp.

- getWDaughterMomentum(double theta)
Given a point on the second W daughter ellipse, defined by an angle theta, calculate the corresponding momentum (px, py, pz).


### lightJetMinimumChiSquareSolver

In the Betchart paper, the addition of a MET constraint is equivalent to
requiring that the two neutrino ellipses intersect. This can be achieved by
varying within their resolutions all the objects present in the event that are
not associated with the top quarks. At a hadron collider these would be light
jets. 
Note: could in theory also apply to leptons, but since these are so well
measured at the LHC we simply focus on jets. However the code can deal with
non-zero lepton resolutions, see below the topEventMinimizer method
setNonTopObjectCollections which are input to the lightJetMinimumChiSquareSolver
class.
In the general case we want the non-top objects to exactly cancel the momentum
imbalance from the decay of the top quarks.

We can write and minimize a chi square variable representing the distance
between the measured and the corrected light jets which make the ellipses
intersect. Transforming from polar to cartesian coordinates yields an analytic
solution for the minimum by inverting a covariance matrix of the sum of all the
object resolutions. The purpose of this class, given an input collection of
objects and their resolutions, is to calculate the minimum and corresponding
corrections to the objects.

#### Constructor

There are two constructors, which require different inputs.

1. With list of input objects and resolutions:
- vector of four-vectors of light jets
- vector of pT resolutions (one for each jet)
- vector of phi resolutions
- vector of eta resolutions
- coordinates (2D or 3D) of the momentum imbalance to cancel

2. With number of objects only:
- number of lights jets
- coordinates (2D or 3D) of the momentum imbalance to cancel

In both cases the momentum imbalance vector is passed by reference, in order
toallow for changes (as the other objects in the event are allowed to vary, the
sum of top momenta will also vary).

#### Methods of interest

- Eval_covariance()
Produces covariance matrix in cartesian and polar/spherical coordinates
based on polar (pT, phi, eta) input.

- Eval_cov_sum()
Once the light jet object collections have been set, calculate the
globalcovariance matrix and invert it. This method only needs to be called once
per event!

- setupEquations()
In case the second constructor was used, pass in the vectors of non-top
objectfour-momenta and resolutions. Calls setCartesianWidths and calcSigmas.

- calcMin()
Calculate the vector of deltas (corrections to the light jets) defined in
Eq.C.7, i.e., calculate the B_i matrices and multiply by the displacement
vector. From the vector of deltas, calculate the chi square. This method needs
to be called only when the displacement vector changes.

- getChiSquare()
Calls calcMin and returns the minimum of the chi square.

- getMinDeltasX/Y/Z()
Returns a vector of px/py/pz deltas at the minimum of the chi square.

### topSystemChiSquare

Parent top system chi square class. The two different decay channels, leptonic and hadronic, are derived classes with decay-specific members and methods; the parent class contains members and methods common to both channels.

Construct a chi square for one top quark and its decay products, which is contructed from corrections (deltas) to the b-jet, one measured W daughter, top mass, and W mass. The chi square represents the degree of likelihood that a top quark candidate four-momentum can be reconstructed from the input decay product momenta. 

If the decay is leptonic, the neutrino must be reconstructed in order to obtain the top momentum, using the Betchart method. If the decay is hadronic, the b-jet and one of the measured W daughters can be used to construct the ellipse that the second W daughter must lie on, also with the Betchart method. The possible second W daughter momenta can then be compared to the measured  second W daughter. 

Each instance has a WDaughterEllipseCalculator instance. The b-jet and first W daughter momenta and top and W masses are expressed as functions of the measured or assumed values, the different object widths, and corrections to the various objects. The WDaughterEllipseCalculator is instantiated with these values.

One limitation of the Betchart paper is that solutions do not exist in the case where the variable Z^2, defined in Section 2.4.2, is negative. We introduce a method to get around this issue: Z^2 == 0 is equivalent to a quadratic equation in the top mass squared. We can therefore solve for the range(s) of top mass values where Z^2 is positive, and restrict corrections to the top mass to stay within the corresponding range.

#### Constructor

The topSystemChiSquare constructor takes as input the four-momentum of a b-jet and its pT, eta and phi widths, the four-momentum of a measured W daughter (lepton or light quark jet) and its pT, eta and phi widths, an assumed top mass and width, an assumed W mass and width, and the assumed mass of the second W daughter (neutrino or light quark jet). The coordinates are passed by reference in order to be able to easily update the momenta when each object is varied.

The corrected b-jet and W daughter four-momenta and the top and W masses are used to initialize an instance of the WDaughterEllipseCalculator class.

#### Methods of interest

- Various setDelta functions (one for each object) 
Called at each step in the minimizer to set new delta values and update the corresponding variables (object momenta, masses).

- setupWDaughter2Ellipse
Calls the setupEllipse method of the WDaughterEllipseCalculator class. Allows to update the ellipse with each new step in the minimization.

- calcWDaughter2Ellipse
Calculates the ellipse matrix in both representations (Hperp, homogeneous and Nperp, extended). Since the extended representation involves inverting the Hperp matrix, care must be taken to ensure that it is invertible. This is equivalent to requiring Z^2>0, and so the method calcTopMassRange (see below) is first called.

- calcTopMassRange
Calculates the roots mTopEdge1 and mTopEdge 2 of a quadratic equation in mTop^2, then finds the range in which Z^2 is positive. Three cases must be considered:
	1. Both roots are negative: only need to check whether Z^2 is positive in the range mTop>0.
	2. One positive and one negative root: two intervals to check, [0; sqrt(mTopEdge2)] and [sqrt(mTopEdge2); +inf[.
	3. Both roots are positive: three intervals to check, [0; sqrt(mTopEdge1)], [sqrt(mTopEdge1); sqrt(mTopEdge2)] and [sqrt(mTopEdge2);+inf[. 
In case more than one top mass interval in which Z^2>0 exists, pick the one either containing or closest to the input assumed top mass. 
Also calculate the corresponding delta range in which the top mass will be allowed to vary during the minimization.

- getTopMassChiSquare
Calculates and returns the component of the top system chi square dependent on the top mass delta. The top mass must be varied inside the inner minimizer (see topEventMinimizer below) because the top mass range where Z^2 is positive depends on the b-jet and first W daughter momenta (and W mass). The seven corresponding parameters belong to the outer minimizer. For each step in the outer minimization, the inner chi square is minimized. Inner minimizer parameters are the top mass deltas and angles corresponding to one point on the ellipse, and thus one second W daughter momentum (one each per top system). 


### leptonicTopSystemChiSquare

Class derived from topSystemChiSquare. Implements methods specific to leptonic top decays.

#### Constructor

The leptonicTopSystemChiSquare constructor takes inputs identical to the topSystemChiSquare constructor, which is called within. See above for documentation.

#### Methods of interest

- setEllipseAngle
Sets the position on the second W daughter ellipse; updates the second W daughter momentum; updates the top quark momentum.

- calcChiSquare
Calculates the top system chi square components for parameters other than the top mass: b-jet pT, eta, phi, first W daughter pT, eta, phi, and W mass. Since the second W daughter in this mode is a neutrino, whose true four-momentum we do not know, it does not contribute to the top system chi square.


### hadronicTopSystemChiSquare

Class derived from topSystemChiSquare. Implements methods specific to hadronic top decays, mainly having to do with the fact that the second W daughter is visible in this decay mode. It is integrated into the overall chi square.

#### Constructor

In addition to the inputs required by the topSystemChiSquare constructor, which is called inside, the hadronicTopSystemChiSquare constructor also takes as arguments the measured four-momentum of the second W daughter as well as its pT, eta and phi widths. The coordinates are again passed by reference to allow for easy propagation of corrections.

#### Methods of interest

- setEllipseAngle
See above method in leptonicTopSystemChiSquare description. In addition, after the second W daughter momentum is updated, the function calcWDaughter2Deltas is called (see below).

- calcWDaughter2Deltas
Given an angle on the second W daughter ellipse, and thus a momentum, calculate the variation with respect to the input measured momentum. In this sense the second W daughter deltas are different from the b-jet or first W daughter (lepton or quark) deltas: the inner minimizer parameters are the ellipse angles and top mass deltas; one ellipse angle _defines_ the second W daughter deltas, while for the other (outer minimizer) parameters each step in the minimization starts with new delta values.

- calcChiSquare
See above description in leptonicTopSystemChiSquare class.

- getHadronicChiSquare
Calculates and returns the second W daughter contribution to the top system chi square. However, since the corresponding variable, the ellipse angle, is a parameter in the inner minimizer, the chi square must be calculated separately.


### topEventMinimizer

Functions at the event level to minimize the per-event chi square defined in Eq. 1, which is constructed from the sum of top system chi squares (one per top assumed to be present in the event, each represented by one topSystemChiSquare object) and the non-top object chi square (combining all measured objects in the event not assumed to originate from a top decay, and represented by one lightJetMinimiumChiSquareSolver object).

As can be seen in Eq. 3, the minimum chi square can be expressed as a series of nested minimizations. We implement an outer and an inner minimizer.

Although in practice the code has been tested using two assumed tops, it can in theory accommodate any number of tops, decaying in either mode. This is achieved with a vector of pairs of (a pointer to) a topSystemChiSquare object, and a boolean holding the mode of decay (true for leptonic, false for hadronic). The user then adds the desired number of tops to the event using the methods addLeptonicTop and addHadronicTop (see below).

#### Constructor

Two constructors (get rid of one?) exist:

1. Input only non-top objects and later construct the tops
The only necessary inputs in this case are a vector of non-top object four-momenta, vectors of pT, eta, and phi widths, and the assumed top and W masses and widths. The user can then construct tops later using the addLeptonicTop and/or addHadronicTop methods.

2. Input all measured objects and construct the tops directly
Input a vector of four-momenta of all the measured objects in the event, the corresponding vectors of pT, eta and phi widths, the assumed top and W masses and widths, vectors of integers holding the indices of the b-jets, first and second W daughters in the vector of all objects, and a vector of booleans holding the type of decay (leptonic or hadronic). The different top systems are created inside the constructor.


#### Methods of interest

- addLeptonicTop
Add a leptonically decaying top. Calls the leptonicTopSystemChiSquare constructor; the inputs are identical.

- addHadronicTop
Add a hadronically decaying top. Calls the hadronicTopSystemChiSquare constructor; the inputs are identical.

- findStartingValues
Find starting values for the ellipse angles by looping around ellipses. At each point the sum of top momenta is recalculated and the non-top objects are set to recoil against it; the sum of the non-top (light jet) and hadronic chi squares is then computed. The values corresponding to the minimum are then used as a starting point for the inner minimizer.

- minimizeNonTopChiSquare
The inner minimizer parameters are the ellipse angles and top mass deltas (one each per top system). For each top system the top mass delta parameter is restricted in range by the edges computed in the topSystemChiSquare::calcTopMassRange method (see above).  The chi square to minimize is the sum of top mass, hadronic second W daughter, and non-top object chi squares.

- minimizeTotalChiSquare
The outer minimizer parameters are the b-jet pT, eta and phi deltas, the first W daughter pT, eta and phi deltas, and the W mass deltas (one each per top system). At each step in the minimization procedure, the second W daughter ellipses are recalculated using the updated b-jet and first W daughter momenta and W mass. Next the inner piece is calculated, ie, the inner minimization performed. The total chi square is defined as the sum of the minimum of the inner chi square at this outer step, and the sum of top system chi squares.

- getBestDeltas
Once the full minimization has been run, retrieve the delta values at the minimum.


### topReconstructionFromLHE

Class that tests the above code on the simulated ttbar sample (see below). It gives a direct example of how to set everything up, do the minimization, and retrieve the output. Also writes some values out to a TTree in a file "output.root" which can then be inspected.

Currently set up to do semi-leptonic ttbar decays.

### lhe.root

Simulation sample: top quark pair production in association with two jets. The sample was generated with Madgraph and the output LHE file is read into an instance of the topReconstructionFromLHE class.

### loadFiles

Compile and run the code on the test file lhe.root.



# Top-Reco-3
