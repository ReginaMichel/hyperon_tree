/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* AliAnaysisTaskHypertritonKFTreeLocal
 *
 */

#include <iostream>

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"

#include "AliMCEvent.h"
#include "AliEventCuts.h"
#include "AliTimeRangeCut.h"
#include "AliPIDResponse.h"
#include "AliMultSelection.h"
#include "AliESDtrackCuts.h"
#include "TChain.h"

// added for retro correction
#include "TGeoManager.h"
#include "AliCDBManager.h"
#include "AliGeomManager.h"
#include "TGeoGlobalMagField.h"
#include "AliGRPManager.h"
#include "AliPID.h"
#include "AliExternalTrackParam.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliAnalysisTaskSE.h"
#include "AliPIDResponse.h"
#include "AliAnalysisTaskWeakDecayVertexer.h"
#include "AliTrackerBase.h"
#include "AliESDv0.h"

// added for weighting
#include "AliPDG.h"
#include "AliPWGFunc.h"

#include "AliAnalysisTaskHypertritonKFTreeLocal.h"

/// your analysis class
class AliAnalysisTaskHypertritonKFTreeLocal;
/// classimp: necessary for root
ClassImp(AliAnalysisTaskHypertritonKFTreeLocal)
///_____________________________________________________________________________
AliAnalysisTaskHypertritonKFTreeLocal::AliAnalysisTaskHypertritonKFTreeLocal() : AliAnalysisTaskSE(),
fEventCuts(),
fTimeRangeCut(),
fPIDResponse(nullptr),
PrimVertex(),
He3CandTrackIdPos(10),
He3CandTrackIdNeg(10),
PionCandTrackIdPos(1500),
PionCandTrackIdNeg(1500),
HyperTritonCandidates(50),
kDoQA(false),
kIsMC(false),
kRunRetroCorrection(false),
fESDtrackCuts(nullptr),
fQAList(nullptr),
fOutputList(nullptr),
fMCEvent(nullptr),
hNumberOfEvents(nullptr),
histoEventCentrality(nullptr),
CandidateTree(nullptr),
CentralityPercentile(-999),
mass(-999),
ErrorMass(-999),
px(-999),
py(-999),
pz(-999),
p(-999),
pT(-999),
Rapidity(-999),
Charge(0),
Chi2PerNDF(-999),
CosPointingAngle(-999),
DistanceToPV(-999),
DeviationFromPV(-999),
DistanceToPVXY(-999),
DeviationFromPVXY(-999),
runN(0),
magB(-999),
massUncorr(-999),
massCorrDiff(-999),
ErrorMassUncorr(-999),
pxUncorr(-999),
pyUncorr(-999),
pzUncorr(-999),
pUncorr(-999),
pTUncorr(-999),
pxPionUncorr(-999),
pyPionUncorr(-999),
pzPionUncorr(-999),
pPionUncorr(-999),
pTPionUncorr(-999),
pxHeUncorr(-999),
pyHeUncorr(-999),
pzHeUncorr(-999),
pHeUncorr(-999),
pTHeUncorr(-999),
xSecVertexUncorr(-999),
ySecVertexUncorr(-999),
zSecVertexUncorr(-999),
xSecVertexCorrDiff(-999),
ySecVertexCorrDiff(-999),
zSecVertexCorrDiff(-999),
SecVertexCorrDiff(-999),
massSV(-999),
massSVDiff(-999),
ErrorMassSV(-999),
pxSV(-999),
pySV(-999),
pzSV(-999),
pSV(-999),
pTSV(-999),
xSecVertexSV(-999),
ySecVertexSV(-999),
zSecVertexSV(-999),
xSecVertexSVDiff(-999),
ySecVertexSVDiff(-999),
zSecVertexSVDiff(-999),
SecVertexSVDiff(-999),
xPionTrack(-999),
yPionTrack(-999),
xPionPV(-999),
yPionPV(-999),
xPionTPC(-999),
yPionTPC(-999),
xPionSV(-999),
yPionSV(-999),
xHeTrack(-999),
yHeTrack(-999),
xHePV(-999),
yHePV(-999),
xHeTPC(-999),
yHeTPC(-999),
xHeSV(-999),
yHeSV(-999),
lxPionTPC(-999),
lxPionSV(-999),
lxHeTPC(-999),
lxHeSV(-999),
massTopo(-999),
ErrorMassTopo(-999),
pxTopo(-999),
pyTopo(-999),
pzTopo(-999),
pTopo(-999),
pTTopo(-999),
pxVarianceTopo(-999),
pyVarianceTopo(-999),
pzVarianceTopo(-999),
RapidityTopo(-999),
Chi2PerNDFTopo(-999),
CosPointingAngleTopo(-999),
DistanceToPVTopo(-999),
DeviationFromPVTopo(-999),
DistanceToPVXYTopo(-999),
DeviationFromPVXYTopo(-999),
DecayLength(-999),
ErrorDecayLength(-999),
DecayLengthXY(-999),
ErrorDecayLengthXY(-999),
DistanceOfDaughters(-999),
DeviationOfDaughters(-999),
DistanceOfDaughtersXY(-999),
DeviationOfDaughtersXY(-999),
pxPion(-999),
pyPion(-999),
pzPion(-999),
pPion(-999),
pTPion(-999),
pxPionSV(-999),
pyPionSV(-999),
pzPionSV(-999),
pPionSV(-999),
pTPionSV(-999),
ChargePion(-999),
DCAPion(-999),
DeviationFromPVPion(-999),
DCAPionXY(-999),
DeviationFromPVPionXY(-999),
NCrossedRowsTPCPion(-999),
NPIDClusterTPCPion (-999),
TPCMomPion(-999),
TPCnSigmaPion(-999),
HasPointOnITSLayer0Pion(false),
HasPointOnITSLayer1Pion(false),
HasPointOnITSLayer2Pion(false),
HasPointOnITSLayer3Pion(false),
HasPointOnITSLayer4Pion(false),
HasPointOnITSLayer5Pion(false),
HasITSRefitPion(false),
PIDForTrackingPion(-999),
DistanceToSecVertPion(-999),
DeviationToSecVertPion(-999),
pxHe(-999),
pyHe(-999),
pzHe(-999),
pHe(-999),
pTHe(-999),
pxHeSV(-999),
pyHeSV(-999),
pzHeSV(-999),
pHeSV(-999),
pTHeSV(-999),
ChargeHe(-999),
DCA3He(-999),
DeviationFromPV3He(-999),
DCA3HeXY(-999),
DeviationFromPV3HeXY(-999),
NCrossedRowsTPC3He(-999),
NPIDClusterTPC3He(-999),
TPCMom3He(-999),
TPCnSigma3He(-999),
TPCnSigma3H(-999),
HasPointOnITSLayer0He3(false),
HasPointOnITSLayer1He3(false),
HasPointOnITSLayer2He3(false),
HasPointOnITSLayer3He3(false),
HasPointOnITSLayer4He3(false),
HasPointOnITSLayer5He3(false),
HasITSRefitHe3(false),
PIDForTrackingHe3(-999),
DistanceToSecVertHe(-999),
DeviationToSecVertHe(-999),
///Monte Carlo
weight(-999),
pxMC(-999),
pyMC(-999),
pzMC(-999),
pxHeMC(-999),
pyHeMC(-999),
pzHeMC(-999),
pHeMC(-999),
pTHeMC(-999),
pxPionMC(-999),
pyPionMC(-999),
pzPionMC(-999),
pPionMC(-999),
pTPionMC(-999),
pxVariance(-999),
pyVariance(-999),
pzVariance(-999),
pxHeVariance(-999),
pyHeVariance(-999),
pzHeVariance(-999),
pxPionVariance(-999),
pyPionVariance(-999),
pzPionVariance(-999),
pMC(-999),
pTMC(-999),
RapidityMC(-999),
DecayLengthMC(-999),
DecayLengthXYMC(-999),
xSecVertexMC(-999),
ySecVertexMC(-999),
zSecVertexMC(-999),
xSecVertex(-999),
ySecVertex(-999),
zSecVertex(-999),
rSecVertex(-999),
xSecVertexVariance(-999),
ySecVertexVariance(-999),
zSecVertexVariance(-999),
GeneratedTreeMC(nullptr),
fHistNsigmaTPCvsP3He(nullptr),
fHistNsigmaTPCvsPPion(nullptr),
fHistpxTrueRecHe3(nullptr),
fHistpyTrueRecHe3(nullptr),
fHistpzTrueRecHe3(nullptr),
fHistMomPion(nullptr),
fHistMomHe3(nullptr)
{
  /// default constructor, don't allocate memory here!
  /// this is used by root for IO purposes, it needs to remain empty
}
///_____________________________________________________________________________
AliAnalysisTaskHypertritonKFTreeLocal::AliAnalysisTaskHypertritonKFTreeLocal(const char* name) : AliAnalysisTaskSE(name),
fEventCuts(),
fTimeRangeCut(),
fPIDResponse(nullptr),
PrimVertex(),
He3CandTrackIdPos(10),
He3CandTrackIdNeg(10),
PionCandTrackIdPos(1500),
PionCandTrackIdNeg(1500),
HyperTritonCandidates(50),
kDoQA(false),
kIsMC(false),
kRunRetroCorrection(false),
fESDtrackCuts(nullptr),
fQAList(nullptr),
fOutputList(nullptr),
fMCEvent(nullptr),
hNumberOfEvents(nullptr),
histoEventCentrality(nullptr),
CandidateTree(nullptr),
CentralityPercentile(-999),
mass(-999),
ErrorMass(-999),
px(-999),
py(-999),
pz(-999),
p(-999),
pT(-999),
Rapidity(-999),
Charge(0),
Chi2PerNDF(-999),
CosPointingAngle(-999),
DistanceToPV(-999),
DeviationFromPV(-999),
DistanceToPVXY(-999),
DeviationFromPVXY(-999),
runN(0),
magB(-999),
massUncorr(-999),
massCorrDiff(-999),
ErrorMassUncorr(-999),
pxUncorr(-999),
pyUncorr(-999),
pzUncorr(-999),
pUncorr(-999),
pTUncorr(-999),
pxPionUncorr(-999),
pyPionUncorr(-999),
pzPionUncorr(-999),
pPionUncorr(-999),
pTPionUncorr(-999),
pxHeUncorr(-999),
pyHeUncorr(-999),
pzHeUncorr(-999),
pHeUncorr(-999),
pTHeUncorr(-999),
xSecVertexUncorr(-999),
ySecVertexUncorr(-999),
zSecVertexUncorr(-999),
xSecVertexCorrDiff(-999),
ySecVertexCorrDiff(-999),
zSecVertexCorrDiff(-999),
SecVertexCorrDiff(-999),
massSV(-999),
massSVDiff(-999),
ErrorMassSV(-999),
pxSV(-999),
pySV(-999),
pzSV(-999),
pSV(-999),
pTSV(-999),
xSecVertexSV(-999),
ySecVertexSV(-999),
zSecVertexSV(-999),
xSecVertexSVDiff(-999),
ySecVertexSVDiff(-999),
zSecVertexSVDiff(-999),
SecVertexSVDiff(-999),
xPionTrack(-999),
yPionTrack(-999),
xPionPV(-999),
yPionPV(-999),
xPionTPC(-999),
yPionTPC(-999),
xPionSV(-999),
yPionSV(-999),
xHeTrack(-999),
yHeTrack(-999),
xHePV(-999),
yHePV(-999),
xHeTPC(-999),
yHeTPC(-999),
xHeSV(-999),
yHeSV(-999),
lxPionTPC(-999),
lxPionSV(-999),
lxHeTPC(-999),
lxHeSV(-999),
massTopo(-999),
ErrorMassTopo(-999),
pxTopo(-999),
pyTopo(-999),
pzTopo(-999),
pTopo(-999),
pTTopo(-999),
pxVarianceTopo(-999),
pyVarianceTopo(-999),
pzVarianceTopo(-999),
RapidityTopo(-999),
Chi2PerNDFTopo(-999),
CosPointingAngleTopo(-999),
DistanceToPVTopo(-999),
DeviationFromPVTopo(-999),
DistanceToPVXYTopo(-999),
DeviationFromPVXYTopo(-999),
DecayLength(-999),
ErrorDecayLength(-999),
DecayLengthXY(-999),
ErrorDecayLengthXY(-999),
DistanceOfDaughters(-999),
DeviationOfDaughters(-999),
DistanceOfDaughtersXY(-999),
DeviationOfDaughtersXY(-999),
pxPion(-999),
pyPion(-999),
pzPion(-999),
pPion(-999),
pTPion(-999),
pxPionSV(-999),
pyPionSV(-999),
pzPionSV(-999),
pPionSV(-999),
pTPionSV(-999),
ChargePion(-999),
DCAPion(-999),
DeviationFromPVPion(-999),
DCAPionXY(-999),
DeviationFromPVPionXY(-999),
NCrossedRowsTPCPion(-999),
NPIDClusterTPCPion (-999),
TPCMomPion(-999),
TPCnSigmaPion(-999),
HasPointOnITSLayer0Pion(false),
HasPointOnITSLayer1Pion(false),
HasPointOnITSLayer2Pion(false),
HasPointOnITSLayer3Pion(false),
HasPointOnITSLayer4Pion(false),
HasPointOnITSLayer5Pion(false),
HasITSRefitPion(false),
PIDForTrackingPion(-999),
DistanceToSecVertPion(-999),
DeviationToSecVertPion(-999),
pxHe(-999),
pyHe(-999),
pzHe(-999),
pHe(-999),
pTHe(-999),
pxHeSV(-999),
pyHeSV(-999),
pzHeSV(-999),
pHeSV(-999),
pTHeSV(-999),
ChargeHe(-999),
DCA3He(-999),
DeviationFromPV3He(-999),
DCA3HeXY(-999),
DeviationFromPV3HeXY(-999),
NCrossedRowsTPC3He(-999),
NPIDClusterTPC3He(-999),
TPCMom3He(-999),
TPCnSigma3He(-999),
TPCnSigma3H(-999),
HasPointOnITSLayer0He3(false),
HasPointOnITSLayer1He3(false),
HasPointOnITSLayer2He3(false),
HasPointOnITSLayer3He3(false),
HasPointOnITSLayer4He3(false),
HasPointOnITSLayer5He3(false),
HasITSRefitHe3(false),
PIDForTrackingHe3(-999),
DistanceToSecVertHe(-999),
DeviationToSecVertHe(-999),
///Monte Carlo
weight(-999),
pxMC(-999),
pyMC(-999),
pzMC(-999),
pxHeMC(-999),
pyHeMC(-999),
pzHeMC(-999),
pHeMC(-999),
pTHeMC(-999),
pxPionMC(-999),
pyPionMC(-999),
pzPionMC(-999),
pPionMC(-999),
pTPionMC(-999),
pxVariance(-999),
pyVariance(-999),
pzVariance(-999),
pxHeVariance(-999),
pyHeVariance(-999),
pzHeVariance(-999),
pxPionVariance(-999),
pyPionVariance(-999),
pzPionVariance(-999),
pMC(-999),
pTMC(-999),
RapidityMC(-999),
DecayLengthMC(-999),
DecayLengthXYMC(-999),
xSecVertexMC(-999),
ySecVertexMC(-999),
zSecVertexMC(-999),
xSecVertex(-999),
ySecVertex(-999),
zSecVertex(-999),
rSecVertex(-999),
xSecVertexVariance(-999),
ySecVertexVariance(-999),
zSecVertexVariance(-999),
GeneratedTreeMC(nullptr),
fHistNsigmaTPCvsP3He(nullptr),
fHistNsigmaTPCvsPPion(nullptr),
fHistpxTrueRecHe3(nullptr),
fHistpyTrueRecHe3(nullptr),
fHistpzTrueRecHe3(nullptr),
fHistMomPion(nullptr),
fHistMomHe3(nullptr)
{
  /// constructor
  /// define the input of the analysis: in this case we take a 'chain' of events
  /// this chain is created by the analysis manager, so no need to worry about it,
  /// it does its work automatically
  
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  /// trees
  DefineOutput(3, TTree::Class());
  DefineOutput(4, TTree::Class());
}
///_____________________________________________________________________________
AliAnalysisTaskHypertritonKFTreeLocal::~AliAnalysisTaskHypertritonKFTreeLocal()
{
  /// destructor
  /// at the end of your task, it is deleted from memory by calling this function
  if(fESDtrackCuts) delete fESDtrackCuts;
  if(fQAList) delete fQAList;
  if(fOutputList) delete fOutputList;
  if(CandidateTree) delete CandidateTree;
  if(GeneratedTreeMC) delete GeneratedTreeMC;
}

///_____________________________________________________________________________
void AliAnalysisTaskHypertritonKFTreeLocal::UserCreateOutputObjects()
{
  /// create output objects
  ///
  /// this function is called ONCE at the start of your analysis (RUNTIME)
  /// here you ceate the histograms that you want to use
  ///
  /// the histograms are in this case added to a tlist, this list is in the end saved
  /// to an output file
  
  /// Runs which are bad for TPC PID due to a problem in the TPC gain in one or two sectors towards the end of the run
  //  fEventCuts.UseTimeRangeCut(); /// at the moment implemented by hand
  
  /// Create object for track selection
  fESDtrackCuts = new AliESDtrackCuts("fESDtrackCuts");
  fESDtrackCuts->SetRequireTPCRefit(kTRUE);
  fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
  fESDtrackCuts->SetMinNClustersTPC(60);
  fESDtrackCuts->SetMaxChi2PerClusterTPC(5.);
  fESDtrackCuts->SetEtaRange(-0.9,0.9);
  
  fQAList = new TList();          /// list which will contain all of QA histograms
  fQAList->SetOwner(kTRUE);
  
  fOutputList = new TList();     /// list which will contain all of your histograms at the end of the analysis, the contents of this list are written to the output file
  fOutputList->SetOwner(kTRUE);
  
  hNumberOfEvents = new TH1F("hNumberOfEvents","Events at different selection steps ",7,0,7);
  hNumberOfEvents->GetXaxis()->SetTitle("Selection step");
  hNumberOfEvents->Sumw2();
  fOutputList->Add(hNumberOfEvents);
  
  histoEventCentrality = new TH1I("histoEventCentrality","Events vs centrality percentile (V0M)",100,0,100);
  histoEventCentrality->GetXaxis()->SetTitle("Centrality percentile");
  histoEventCentrality->Sumw2();
  fOutputList->Add(histoEventCentrality);
  
  /// QA histograms
  if (kDoQA) {
    /// Add event selection QA plots
    fEventCuts.AddQAplotsToList(fQAList);
    /// Add own QA plots
    fHistNsigmaTPCvsPPion = new TH2F("fHistNsigmaTPCvsPPion","#it{N}#sigma^{TPC}_{#pi} vs #it{p} (GeV/#it{c});#it{p} (GeV/#it{c});#it{N}#sigma^{TPC}_{#pi}", 100,0,10,100,-5.,5.);
    fHistNsigmaTPCvsPPion->Sumw2();
    fQAList->Add(fHistNsigmaTPCvsPPion);
    
    fHistNsigmaTPCvsP3He = new TH2F("fHistNsigmaTPCvsP3He","#it{N}#sigma^{TPC}_{^{3}He} vs #it{p} (GeV/#it{c});#it{p} (GeV/#it{c});#it{N}#sigma^{TPC}_{^{3}He}", 100,0,10,100,-5.,5.);
    fHistNsigmaTPCvsP3He->Sumw2();
    fQAList->Add(fHistNsigmaTPCvsP3He);
  }
  
  /// Data
    ///  Tree with hypertriton candidates
    CandidateTree = new TTree("CandidateTree","CandidateTree"); /// Tree has to be created to ensure that all output files are created, otherwise the task fails output validation
    CandidateTree->Branch("CentralityPercentile",&CentralityPercentile,"CentralityPercentile/F");
    CandidateTree->Branch("weight",&weight,"weight/F");
    CandidateTree->Branch("mass",&mass,"mass/F");
    CandidateTree->Branch("ErrorMass",&ErrorMass,"ErrorMass/F");
    CandidateTree->Branch("px",&px,"px/F");
    CandidateTree->Branch("py",&py,"py/F");
    CandidateTree->Branch("pz",&pz,"pz/F");
    CandidateTree->Branch("p",&p,"p/F");
    CandidateTree->Branch("pT",&pT,"pT/F");
    CandidateTree->Branch("Rapidity",&Rapidity,"Rapidity/F");
    CandidateTree->Branch("Charge",&Charge,"Charge/I");
    CandidateTree->Branch("Chi2PerNDF",&Chi2PerNDF,"Chi2PerNDF/F");
    CandidateTree->Branch("CosPointingAngle",&CosPointingAngle,"CosPointingAngle/F");
    CandidateTree->Branch("DistanceToPV",&DistanceToPV,"DistanceToPV/F");
    CandidateTree->Branch("DistanceToPVXY",&DistanceToPVXY,"DistanceToPVXY/F");
    CandidateTree->Branch("DeviationFromPV",&DeviationFromPV,"DeviationFromPV/F");
    CandidateTree->Branch("DeviationFromPVXY",&DeviationFromPVXY,"DeviationFromPVXY/F");
    
    if (kRunRetroCorrection) {
        CandidateTree->Branch("runN",&runN,"runN/I");
        CandidateTree->Branch("magB",&magB,"magB/D");
        CandidateTree->Branch("massUncorr",&massUncorr,"massUncorr/F");
        CandidateTree->Branch("massCorrDiff",&massCorrDiff,"massCorrDiff/F");
        CandidateTree->Branch("ErrorMassUncorr",&ErrorMassUncorr,"ErrorMassUncorr/F");
        CandidateTree->Branch("pxUncorr",&pxUncorr,"pxUncorr/F");
        CandidateTree->Branch("pyUncorr",&pyUncorr,"pyUncorr/F");
        CandidateTree->Branch("pzUncorr",&pzUncorr,"pzUncorr/F");
        CandidateTree->Branch("pUncorr",&pUncorr,"pUncorr/F");
        CandidateTree->Branch("pTUncorr",&pTUncorr,"pTUncorr/F");
        CandidateTree->Branch("pxPionUncorr",&pxPionUncorr,"pxPionUncorr/F");
        CandidateTree->Branch("pyPionUncorr",&pyPionUncorr,"pyPionUncorr/F");
        CandidateTree->Branch("pzPionUncorr",&pzPionUncorr,"pzPionUncorr/F");
        CandidateTree->Branch("pPionUncorr",&pPionUncorr,"pPionUncorr/F");
        CandidateTree->Branch("pTPionUncorr",&pTPionUncorr,"pTPionUncorr/F");
        CandidateTree->Branch("pxHeUncorr",&pxHeUncorr,"pxHeUncorr/F");
        CandidateTree->Branch("pyHeUncorr",&pyHeUncorr,"pyHeUncorr/F");
        CandidateTree->Branch("pzHeUncorr",&pzHeUncorr,"pzHeUncorr/F");
        CandidateTree->Branch("pHeUncorr",&pHeUncorr,"pHeUncorr/F");
        CandidateTree->Branch("pTHeUncorr",&pTHeUncorr,"pTHeUncorr/F");
        CandidateTree->Branch("xSecVertexUncorr",&xSecVertexUncorr,"xSecVertexUncorr/F");
        CandidateTree->Branch("ySecVertexUncorr",&ySecVertexUncorr,"ySecVertexUncorr/F");
        CandidateTree->Branch("zSecVertexUncorr",&zSecVertexUncorr,"zSecVertexUncorr/F");
        CandidateTree->Branch("xSecVertexCorrDiff",&xSecVertexCorrDiff,"xSecVertexCorrDiff/F");
        CandidateTree->Branch("ySecVertexCorrDiff",&ySecVertexCorrDiff,"ySecVertexCorrDiff/F");
        CandidateTree->Branch("zSecVertexCorrDiff",&zSecVertexCorrDiff,"zSecVertexCorrDiff/F");
        CandidateTree->Branch("SecVertexCorrDiff",&SecVertexCorrDiff,"SecVertexCorrDiff/F");
        CandidateTree->Branch("massSV",&massSV,"massSV/F");
        CandidateTree->Branch("massSVDiff",&massSVDiff,"massSVDiff/F");
        CandidateTree->Branch("ErrorMassSV",&ErrorMassSV,"ErrorMassSV/F");
        CandidateTree->Branch("pxSV",&pxSV,"pxSV/F");
        CandidateTree->Branch("pySV",&pySV,"pySV/F");
        CandidateTree->Branch("pzSV",&pzSV,"pzSV/F");
        CandidateTree->Branch("pSV",&pSV,"pSV/F");
        CandidateTree->Branch("pTSV",&pTSV,"pTSV/F");
        CandidateTree->Branch("xSecVertexSV",&xSecVertexSV,"xSecVertexSV/F");
        CandidateTree->Branch("ySecVertexSV",&ySecVertexSV,"ySecVertexSV/F");
        CandidateTree->Branch("zSecVertexSV",&zSecVertexSV,"zSecVertexSV/F");
        CandidateTree->Branch("xSecVertexSVDiff",&xSecVertexSVDiff,"xSecVertexSVDiff/F");
        CandidateTree->Branch("ySecVertexSVDiff",&ySecVertexSVDiff,"ySecVertexSVDiff/F");
        CandidateTree->Branch("zSecVertexSVDiff",&zSecVertexSVDiff,"zSecVertexSVDiff/F");
        CandidateTree->Branch("SecVertexSVDiff",&SecVertexSVDiff,"SecVertexSVDiff/F");
        CandidateTree->Branch("xPionTrack",&xPionTrack,"xPionTrack/F");
        CandidateTree->Branch("yPionTrack",&yPionTrack,"yPionTrack/F");
        CandidateTree->Branch("xPionPV",&xPionPV,"xPionPV/F");
        CandidateTree->Branch("yPionPV",&yPionPV,"yPionPV/F");
        CandidateTree->Branch("xPionTPC",&xPionTPC,"xPionTPC/F");
        CandidateTree->Branch("yPionTPC",&yPionTPC,"yPionTPC/F");
        CandidateTree->Branch("xPionSV",&xPionSV,"xPionSV/F");
        CandidateTree->Branch("yPionSV",&yPionSV,"yPionSV/F");
        CandidateTree->Branch("xHeTrack",&xHeTrack,"xHeTrack/F");
        CandidateTree->Branch("yHeTrack",&yHeTrack,"yHeTrack/F");
        CandidateTree->Branch("xHePV",&xHePV,"xHePV/F");
        CandidateTree->Branch("yHePV",&yHePV,"yHePV/F");
        CandidateTree->Branch("xHeTPC",&xHeTPC,"xHeTPC/F");
        CandidateTree->Branch("yHeTPC",&yHeTPC,"yHeTPC/F");
        CandidateTree->Branch("xHeSV",&xHeSV,"xHeSV/F");
        CandidateTree->Branch("yHeSV",&yHeSV,"yHeSV/F");
        CandidateTree->Branch("lxPionTPC",&lxPionTPC,"lxPionTPC/D");
        CandidateTree->Branch("lxPionSV",&lxPionSV,"lxPionSV/D");
        CandidateTree->Branch("lxHeTPC",&lxHeTPC,"lxHeTPC/D");
        CandidateTree->Branch("lxHeSV",&lxHeSV,"lxHeSV/D");
        CandidateTree->Branch("HasITSRefitPion",&HasITSRefitPion,"HasITSRefitPion/O");
        CandidateTree->Branch("HasITSRefitHe3",&HasITSRefitHe3,"HasITSRefitHe3/O");
    }
    
    CandidateTree->Branch("massTopo",&massTopo,"massTopo/F");
    CandidateTree->Branch("ErrorMassTopo",&ErrorMassTopo,"ErrorMassTopo/F");
    CandidateTree->Branch("pxTopo",&pxTopo,"pxTopo/F");
    CandidateTree->Branch("pyTopo",&pyTopo,"pyTopo/F");
    CandidateTree->Branch("pzTopo",&pzTopo,"pzTopo/F");
    CandidateTree->Branch("pTopo",&pTopo,"pTopo/F");
    CandidateTree->Branch("pTTopo",&pTTopo,"pTTopo/F");
    CandidateTree->Branch("RapidityTopo",&RapidityTopo,"RapidityTopo/F");
    CandidateTree->Branch("Chi2PerNDFTopo",&Chi2PerNDFTopo,"Chi2PerNDFTopo/F");
    CandidateTree->Branch("CosPointingAngleTopo",&CosPointingAngleTopo,"CosPointingAngleTopo/F");
    CandidateTree->Branch("DistanceToPVTopo",&DistanceToPVTopo,"DistanceToPVTopo/F");
    CandidateTree->Branch("DistanceToPVXYTopo",&DistanceToPVXYTopo,"DistanceToPVXYTopo/F");
    CandidateTree->Branch("DeviationFromPVTopo",&DeviationFromPVTopo,"DeviationFromPVTopo/F");
    CandidateTree->Branch("DeviationFromPVXYTopo",&DeviationFromPVXYTopo,"DeviationFromPVXYTopo/F");
    
    CandidateTree->Branch("DecayLength",&DecayLength,"DecayLength/F");
    CandidateTree->Branch("ErrorDecayLength",&ErrorDecayLength,"ErrorDecayLength/F");
    CandidateTree->Branch("DecayLengthXY",&DecayLengthXY,"DecayLengthXY/F");
    CandidateTree->Branch("ErrorDecayLengthXY",&ErrorDecayLengthXY,"ErrorDecayLengthXY/F");
    
    /// Daughter variables
    CandidateTree->Branch("DistanceOfDaughters",&DistanceOfDaughters,"DistanceOfDaughters/F");
    CandidateTree->Branch("DistanceOfDaughtersXY",&DistanceOfDaughtersXY,"DistanceOfDaughtersXY/F");
    
    CandidateTree->Branch("DeviationOfDaughters",&DeviationOfDaughters,"DeviationOfDaughters/F");
    CandidateTree->Branch("DeviationOfDaughtersXY",&DeviationOfDaughtersXY,"DeviationOfDaughtersXY/F");
    
    CandidateTree->Branch("pxPion",&pxPion,"pxPion/F");
    CandidateTree->Branch("pyPion",&pyPion,"pyPion/F");
    CandidateTree->Branch("pzPion",&pzPion,"pzPion/F");
    CandidateTree->Branch("pPion",&pPion,"pPion/F");
    CandidateTree->Branch("pTPion",&pTPion,"pTPion/F");
    CandidateTree->Branch("pxPionSV",&pxPionSV,"pxPionSV/F");
    CandidateTree->Branch("pyPionSV",&pyPionSV,"pyPionSV/F");
    CandidateTree->Branch("pzPionSV",&pzPionSV,"pzPionSV/F");
    CandidateTree->Branch("pPionSV",&pPionSV,"pPionSV/F");
    CandidateTree->Branch("pTPionSV",&pTPionSV,"pTPionSV/F"); 
    CandidateTree->Branch("ChargePion",&ChargePion,"ChargePion/I");
    CandidateTree->Branch("DCAPion",&DCAPion,"DCAPion/F");
    CandidateTree->Branch("DeviationFromPVPion",&DeviationFromPVPion,"DeviationFromPVPion/F");
    CandidateTree->Branch("DCAPionXY",&DCAPionXY,"DCAPionXY/F");
    CandidateTree->Branch("DeviationFromPVPionXY",&DeviationFromPVPionXY,"DeviationFromPVPionXY/F");
    CandidateTree->Branch("NCrossedRowsTPCPion",&NCrossedRowsTPCPion,"NCrossedRowsTPCPion/I");
    CandidateTree->Branch("NPIDClusterTPCPion",&NPIDClusterTPCPion,"NPIDClusterTPCPion/I");
    CandidateTree->Branch("TPCMomPion",&TPCMomPion,"TPCMomPion/F");
    CandidateTree->Branch("TPCnSigmaPion",&TPCnSigmaPion,"TPCnSigmaPion/F");
    CandidateTree->Branch("HasPointOnITSLayer0Pion",&HasPointOnITSLayer0Pion,"HasPointOnITSLayer0Pion/O");
    CandidateTree->Branch("HasPointOnITSLayer1Pion",&HasPointOnITSLayer1Pion,"HasPointOnITSLayer1Pion/O");
    CandidateTree->Branch("HasPointOnITSLayer2Pion",&HasPointOnITSLayer2Pion,"HasPointOnITSLayer2Pion/O");
    CandidateTree->Branch("HasPointOnITSLayer3Pion",&HasPointOnITSLayer3Pion,"HasPointOnITSLayer3Pion/O");
    CandidateTree->Branch("HasPointOnITSLayer4Pion",&HasPointOnITSLayer4Pion,"HasPointOnITSLayer4Pion/O");
    CandidateTree->Branch("HasPointOnITSLayer5Pion",&HasPointOnITSLayer5Pion,"HasPointOnITSLayer5Pion/O");
    CandidateTree->Branch("PIDForTrackingPion",&PIDForTrackingPion,"PIDForTrackingPion/I");
    CandidateTree->Branch("DistanceToSecVertPion",&DistanceToSecVertPion,"DistanceToSecVertPion/F");
    CandidateTree->Branch("DeviationToSecVertPion",&DeviationToSecVertPion,"DeviationToSecVertPion/F");
    
    CandidateTree->Branch("pxHe",&pxHe,"pxHe/F");
    CandidateTree->Branch("pyHe",&pyHe,"pyHe/F");
    CandidateTree->Branch("pzHe",&pzHe,"pzHe/F");
    CandidateTree->Branch("pHe",&pHe,"pHe/F");
    CandidateTree->Branch("pTHe",&pTHe,"pTHe/F");
    CandidateTree->Branch("pxHeSV",&pxHeSV,"pxHeSV/F");
    CandidateTree->Branch("pyHeSV",&pyHeSV,"pyHeSV/F");
    CandidateTree->Branch("pzHeSV",&pzHeSV,"pzHeSV/F");
    CandidateTree->Branch("pHeSV",&pHeSV,"pHeSV/F");
    CandidateTree->Branch("pTHeSV",&pTHeSV,"pTHeSV/F"); 
    CandidateTree->Branch("ChargeHe",&ChargeHe,"ChargeHe/I");
    CandidateTree->Branch("DCA3He",&DCA3He,"DCA3He/F");
    CandidateTree->Branch("DeviationFromPV3He",&DeviationFromPV3He,"DeviationFromPV3He/F");
    CandidateTree->Branch("DCA3HeXY",&DCA3HeXY,"DCA3HeXY/F");
    CandidateTree->Branch("DeviationFromPV3HeXY",&DeviationFromPV3HeXY,"DeviationFromPV3HeXY/F");
    CandidateTree->Branch("NCrossedRowsTPC3He",&NCrossedRowsTPC3He,"NCrossedRowsTPC3He/I");
    CandidateTree->Branch("NPIDClusterTPC3He",&NPIDClusterTPC3He,"NPIDClusterTPC3He/I");
    CandidateTree->Branch("TPCMom3He",&TPCMom3He,"TPCMom3He/F");
    CandidateTree->Branch("TPCnSigma3He",&TPCnSigma3He,"TPCnSigma3He/F");
    CandidateTree->Branch("TPCnSigma3H",&TPCnSigma3H,"TPCnSigma3H/F");
    CandidateTree->Branch("HasPointOnITSLayer0He3",&HasPointOnITSLayer0He3,"HasPointOnITSLayer0He3/O");
    CandidateTree->Branch("HasPointOnITSLayer1He3",&HasPointOnITSLayer1He3,"HasPointOnITSLayer1He3/O");
    CandidateTree->Branch("HasPointOnITSLayer2He3",&HasPointOnITSLayer2He3,"HasPointOnITSLayer2He3/O");
    CandidateTree->Branch("HasPointOnITSLayer3He3",&HasPointOnITSLayer3He3,"HasPointOnITSLayer3He3/O");
    CandidateTree->Branch("HasPointOnITSLayer4He3",&HasPointOnITSLayer4He3,"HasPointOnITSLayer4He3/O");
    CandidateTree->Branch("HasPointOnITSLayer5He3",&HasPointOnITSLayer5He3,"HasPointOnITSLayer5He3/O");
    CandidateTree->Branch("PIDForTrackingHe3",&PIDForTrackingHe3,"PIDForTrackingHe3/I");
    CandidateTree->Branch("DistanceToSecVertHe",&DistanceToSecVertHe,"DistanceToSecVertHe/F");
    CandidateTree->Branch("DeviationToSecVertHe",&DeviationToSecVertHe,"DeviationToSecVertHe/F");
    
    CandidateTree->Branch("xSecVertex",&xSecVertex,"xSecVertex/F");
    CandidateTree->Branch("ySecVertex",&ySecVertex,"ySecVertex/F");
    CandidateTree->Branch("zSecVertex",&zSecVertex,"zSecVertex/F");
    CandidateTree->Branch("rSecVertex",&rSecVertex,"rSecVertex/F");
    
    CandidateTree->Branch("xSecVertexVariance",&xSecVertexVariance,"xSecVertexVariance/F");
    CandidateTree->Branch("ySecVertexVariance",&ySecVertexVariance,"ySecVertexVariance/F");
    CandidateTree->Branch("zSecVertexVariance",&zSecVertexVariance,"zSecVertexVariance/F");
    
    if (kIsMC){
      ///MC truth information
      /// true momentum of the hypertriton candidate
      CandidateTree->Branch("pxMC",&pxMC,"pxMC/F");
      CandidateTree->Branch("pyMC",&pyMC,"pyMC/F");
      CandidateTree->Branch("pzMC",&pzMC,"pzMC/F");
      CandidateTree->Branch("pMC",&pMC,"pMC/F");
      CandidateTree->Branch("pTMC",&pTMC,"pTMC/F");
      /// true momenta of the daughter tracks
      CandidateTree->Branch("pxHeMC",&pxHeMC,"pxHeMC/F");
      CandidateTree->Branch("pyHeMC",&pyHeMC,"pyHeMC/F");
      CandidateTree->Branch("pzHeMC",&pzHeMC,"pzHeMC/F");
      CandidateTree->Branch("pHeMC",&pHeMC,"pHeMC/F");
      CandidateTree->Branch("pTHeMC",&pTHeMC,"pTHeMC/F");
      CandidateTree->Branch("pxPionMC",&pxPionMC,"pxPionMC/F");
      CandidateTree->Branch("pyPionMC",&pyPionMC,"pyPionMC/F");
      CandidateTree->Branch("pzPionMC",&pzPionMC,"pzPionMC/F");
      CandidateTree->Branch("pPionMC",&pPionMC,"pPionMC/F");
      CandidateTree->Branch("pTPionMC",&pTPionMC,"pTPionMC/F");
      
      /// Variance of the momentun of the hypertriton candidate
      CandidateTree->Branch("pxVariance",&pxVariance,"pxVariance/F");
      CandidateTree->Branch("pyVariance",&pyVariance,"pyVariance/F");
      CandidateTree->Branch("pzVariance",&pzVariance,"pzVariance/F");
      CandidateTree->Branch("pxVarianceTopo",&pxVarianceTopo,"pxVarianceTopo/F");
      CandidateTree->Branch("pyVarianceTopo",&pyVarianceTopo,"pyVarianceTopo/F");
      CandidateTree->Branch("pzVarianceTopo",&pzVarianceTopo,"pzVarianceTopo/F");
      /// Variance of the momentum of the daughter tracks
      CandidateTree->Branch("pxHeVariance",&pxHeVariance,"pxHeVariance/F");
      CandidateTree->Branch("pyHeVariance",&pyHeVariance,"pyHeVariance/F");
      CandidateTree->Branch("pzHeVariance",&pzHeVariance,"pzHeVariance/F");
      CandidateTree->Branch("pxPionVariance",&pxPionVariance,"pxPionVariance/F");
      CandidateTree->Branch("pyPionVariance",&pyPionVariance,"pyPionVariance/F");
      CandidateTree->Branch("pzPionVariance",&pzPionVariance,"pzPionVariance/F");
      
      CandidateTree->Branch("RapidityMC",&RapidityMC,"RapidityMC/F");
      CandidateTree->Branch("DecayLengthMC",&DecayLengthMC,"DecayLengthMC/F");
      CandidateTree->Branch("DecayLengthXYMC",&DecayLengthXYMC,"DecayLengthXYMC/F");
      
      CandidateTree->Branch("xSecVertexMC",&xSecVertexMC,"xSecVertexMC/F");
      CandidateTree->Branch("ySecVertexMC",&ySecVertexMC,"ySecVertexMC/F");
      CandidateTree->Branch("zSecVertexMC",&zSecVertexMC,"zSecVertexMC/F");
    }
  
  if (kIsMC){
    /// MC output
      /// generated MC
      GeneratedTreeMC = new TTree("GeneratedTreeMC","GeneratedTreeMC"); /// Tree has to be created to ensure that all output files are created, otherwise the task fails output validation
      GeneratedTreeMC->Branch("CentralityPercentile",&CentralityPercentile,"CentralityPercentile/F");
      GeneratedTreeMC->Branch("weight",&weight,"weight/F");
      GeneratedTreeMC->Branch("pMC",&pMC,"pMC/F");
      GeneratedTreeMC->Branch("pTMC",&pTMC,"pTMC/F");
      GeneratedTreeMC->Branch("RapidityMC",&RapidityMC,"RapidityMC/F");
      GeneratedTreeMC->Branch("DecayLengthMC",&DecayLengthMC,"DecayLengthMC/F");
    
    if (kDoQA) {
      /// QA and corrections
        fHistpxTrueRecHe3 = new TH2F("fHistpxTrueRecHe3", "(#it{p}_{x}^{true} - #it{p}_{x}^{rec}) vs #it{p}_{x}^{rec}", 100, 0.1, 10.1, 120, -3, 3);
        fHistpxTrueRecHe3->GetXaxis()->SetTitle("#it{p}_{x}^{rec} (GeV/#it{c})");
        fHistpxTrueRecHe3->GetYaxis()->SetTitle("(#it{p}_{x}^{true} - #it{p}_{x}^{rec}) (GeV/#it{c})");
        fHistpxTrueRecHe3->Sumw2();
        fQAList->Add(fHistpxTrueRecHe3);
        
        fHistpyTrueRecHe3 = new TH2F("fHistpyTrueRecHe3", "(#it{p}_{y}^{true} - #it{p}_{y}^{rec}) vs #it{p}_{y}^{rec}", 100, 0.1, 10.1, 120, -3, 3);
        fHistpyTrueRecHe3->GetXaxis()->SetTitle("#it{p}_{y}^{rec} (GeV/#it{c})");
        fHistpyTrueRecHe3->GetYaxis()->SetTitle("(#it{p}_{y}^{true} - #it{p}_{y}^{rec}) (GeV/#it{c})");
        fHistpyTrueRecHe3->Sumw2();
        fQAList->Add(fHistpyTrueRecHe3);
        
        fHistpzTrueRecHe3 = new TH2F("fHistpzTrueRecHe3", "(#it{p}_{z}^{true} - #it{p}_{z}^{rec}) vs #it{p}_{z}^{rec}", 100, 0.1, 10.1, 120, -3, 3);
        fHistpzTrueRecHe3->GetXaxis()->SetTitle("#it{p}_{z}^{rec} (GeV/#it{c})");
        fHistpzTrueRecHe3->GetYaxis()->SetTitle("(#it{p}_{z}^{true} - #it{p}_{z}^{rec}) (GeV/#it{c})");
        fHistpzTrueRecHe3->Sumw2();
        fQAList->Add(fHistpzTrueRecHe3);
        
        fHistMomPion = new TH1F("fHistMomPion","Momentum distribution of #pi from ^{3}_{#Lambda}H;#it{p} (GeV/#it{c});Counts",100,0,10);
        fHistMomPion->Sumw2();
        fQAList->Add(fHistMomPion);
        
        fHistMomHe3 = new TH1F("fHistMomHe3","Momentum distribution of ^{3}He from ^{3}_{#Lambda}H;#it{p} (GeV/#it{c});Counts",100,0,10);
        fHistMomHe3->Sumw2();
        fQAList->Add(fHistMomHe3);
    }
  }
  
  PostData(1,fQAList);
  PostData(2,fOutputList);
  
  int NumberOfContainer = 3;
  
    PostData(NumberOfContainer++,CandidateTree);
    if (kIsMC) PostData(NumberOfContainer++,GeneratedTreeMC);
}
///_____________________________________________________________________________
void AliAnalysisTaskHypertritonKFTreeLocal::UserExec(Option_t *)
{
  /// user exec
  /// this function is called once for each event
  /// the manager will take care of reading the events from file, and with the static function InputEvent() you
  /// have access to the current event.
  /// once you return from the UserExec function, the manager will retrieve the next event from the chain
  
    // __ 2018 PbPb pass 3__ // Bethe-Bloch parameters obtained from Janik Ditzel, Benjamin Dönigus, and Esther Bartsch
    
    // Triton parameters
    fBetheParamsT[0] = 0.648689;
    fBetheParamsT[1] = 56.6706;
    fBetheParamsT[2] = -1.63243e-10;
    fBetheParamsT[3] = 2.46921;
    fBetheParamsT[4] = 16.8531;
    fBetheParamsT[5] = 0.06;
    
    // Helium parameters
    fBetheParamsHe[0] = 1.70184;
    fBetheParamsHe[1] = 28.4426;
    fBetheParamsHe[2] = 3.21871e-12;
    fBetheParamsHe[3] = 2.06952;
    fBetheParamsHe[4] = 2.77971;
    fBetheParamsHe[5] = 0.06;  
    
  if ( !PassedEventSelection() ) {
    PostData(1,fQAList);
    PostData(2,fOutputList);
    return;
  }
  
  /// set Magnetic field for KF
  KFParticle::SetField(fInputEvent->GetMagneticField());
  
  /// Create Primary Vertex with KF
  PrimVertex = CreateKFVertex(fInputEvent->GetPrimaryVertex());
  
  if (kIsMC){
    ProcessMC();
  } else {
    ProcessESD();
  }
  
  PostData(1,fQAList);
  PostData(2,fOutputList);
  
  int NumberOfContainer = 3;
    PostData(NumberOfContainer++,CandidateTree);
    if (kIsMC) PostData(NumberOfContainer++,GeneratedTreeMC);
}
///_____________________________________________________________________________
void AliAnalysisTaskHypertritonKFTreeLocal::Terminate(Option_t *)
{
  /// terminate
  /// called at the END of the analysis (when all events are processed)
}
///_____________________________________________________________________________
void AliAnalysisTaskHypertritonKFTreeLocal::ProcessESD()
{
  /// container for found 3He candidates
  He3CandTrackIdPos.clear();
  He3CandTrackIdNeg.clear();
  
  /// container for found pion candidates
  PionCandTrackIdPos.clear();
  PionCandTrackIdNeg.clear();
  
  /// check for helium-3 candidates and store their track number
  for(Int_t i=0; i < fInputEvent->GetNumberOfTracks(); i++) {
    AliESDtrack* track = static_cast<AliESDtrack*>(fInputEvent->GetTrack(i));         /// get a track (type AliESDtrack) from the event
    if(!track) continue;
    if(!PassedBasicTrackQualityCuts(track)) continue;              /// skip track is it fails basic track selection
    
    if (kDoQA) {
      /// Fill QA histogram
      Double_t dEdxSigmaPion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
      if(dEdxSigmaPion < 5) fHistNsigmaTPCvsPPion->Fill(track->P(),dEdxSigmaPion);
      
        Double_t dEdxSigmaHe3 = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kHe3);
        if(dEdxSigmaHe3 < 5) fHistNsigmaTPCvsP3He->Fill(track->P(), dEdxSigmaHe3);
    }
    
    Int_t ChargeTrack = (Int_t) track->Charge();
    
    /// Needed for both decay channels
    if(PionSelection(track)){
      if (ChargeTrack > 0) {
        PionCandTrackIdPos.push_back(i);
      } else{
        PionCandTrackIdNeg.push_back(i);
      }
    }
    
    if(Helium3Selection(track)){
      if (ChargeTrack > 0) {
        He3CandTrackIdPos.push_back(i);
      } else{
        He3CandTrackIdNeg.push_back(i);
      }
    }
  }

    /// Loop only over correct charge combinations -> unlike sign pairs
    /// Automatically removes like sign pairs
  TwoBodyDecay(He3CandTrackIdPos,PionCandTrackIdNeg);
  TwoBodyDecay(He3CandTrackIdNeg,PionCandTrackIdPos);
}
///____________________________________________________________
void AliAnalysisTaskHypertritonKFTreeLocal::TwoBodyDecay(const vector<Int_t>& He3CandTrackId, const vector<Int_t>& PionCandTrackId){
  /// next event if there are no 3He or pion candidates
  if ( (Int_t) He3CandTrackId.size() <= 0) return;
  if ( (Int_t) PionCandTrackId.size() <= 0) return;
  
  /// combine pion candidates with 3He candidates
  for(unsigned int iPion=0; iPion < (UInt_t) PionCandTrackId.size(); iPion++) {
    AliESDtrack* ESDtrackPion = static_cast<AliESDtrack*>(fInputEvent->GetTrack(PionCandTrackId[iPion]));         /// get a track (type AliESDtrack) from the event
    
    Int_t ChargePionTrack = (Int_t) ESDtrackPion->Charge();
    /// create KFParticle for pion
    Int_t pdgPion = 211;
    if (ChargePionTrack < 0) {
      pdgPion = -211;
    }
    KFParticle3LH kfpDaughter2   = CreateKFParticle(ESDtrackPion, 0.13957f, 1);
 //   KFParticle kfpDaughter2 = CreateKFTrack(ESDtrackPion, pdgPion);
    
    for (unsigned int jHe3=0; jHe3 < (UInt_t) He3CandTrackId.size(); jHe3++) {
      /// Make sure not to use a single track twice
      if ( PionCandTrackId[iPion]==He3CandTrackId[jHe3] ) continue;
      
      AliESDtrack* ESDtrackHe3 = static_cast<AliESDtrack*>(fInputEvent->GetTrack(He3CandTrackId[jHe3]));
      Int_t Charge3He = (Int_t) ESDtrackHe3->Charge();
      /// create KFParticle for 3He
      Int_t pdg3He = 1000020030;
      if (Charge3He < 0) {
        pdg3He = -1000020030;
      }
      
      KFParticle3LH kfpDaughter1 = CreateKFParticle(ESDtrackHe3, 2.80839f, 2);
  //    KFParticle kfpDaughter1 = CreateKFTrack(ESDtrackHe3, pdg3He);
      
      if (kfpDaughter1.GetDistanceFromParticle(kfpDaughter2) > 3) continue;
      
      /// create the kfp mother and histograms of transverse momenta and mass
      KFParticle3LH kfpMother;
      kfpMother.AddDaughters(kfpDaughter2, kfpDaughter1, 2); // (daughter1, daughter2, constructMethod (0 or 2)
      // KFParticle kfpMother(kfpDaughter1,kfpDaughter2);
      /// preselect hypertriton candidates before retro-correction
      if (kfpMother.GetNDF()<0) continue;
      if (kfpMother.GetChi2()<0) continue;
      if (kfpMother.GetChi2()>10000) continue; /// protection against infinite
      Float_t DecayLengthUncorr;
      Float_t ErrorDecayLengthUncorr;
      kfpMother.GetDecayLength(DecayLengthUncorr, ErrorDecayLengthUncorr);
      if (DecayLengthUncorr > 60) continue;
      
      /// Fill uncorrected daughter variables:
      pxPionUncorr = ESDtrackPion->Px();
      pyPionUncorr = ESDtrackPion->Py();
      pzPionUncorr = ESDtrackPion->Pz();
      pPionUncorr = TMath::Sqrt(pxPionUncorr * pxPionUncorr + pyPionUncorr * pyPionUncorr + pzPionUncorr * pzPionUncorr);
      pTPionUncorr = TMath::Sqrt(pxPionUncorr * pxPionUncorr + pyPionUncorr * pyPionUncorr);
      pxHeUncorr = 2* ESDtrackHe3->Px();
      pyHeUncorr = 2* ESDtrackHe3->Py();
      pzHeUncorr = 2* ESDtrackHe3->Pz();
      pHeUncorr = TMath::Sqrt(pxHeUncorr * pxHeUncorr + pyHeUncorr * pyHeUncorr + pzHeUncorr * pzHeUncorr);
      pTHeUncorr = TMath::Sqrt(pxHeUncorr * pxHeUncorr + pyHeUncorr * pyHeUncorr);
      kfpMother.GetMass(massUncorr, ErrorMassUncorr);
      if ( massUncorr < 2.92 || massUncorr > 3.07 ) continue;
      pxUncorr = kfpMother.GetPx();
      pyUncorr = kfpMother.GetPy();
      pzUncorr = kfpMother.GetPz();
      pUncorr = TMath::Sqrt(pxUncorr * pxUncorr + pyUncorr * pyUncorr + pzUncorr * pzUncorr);
      pTUncorr = TMath::Sqrt(pxUncorr * pxUncorr + pyUncorr * pyUncorr);
      
      if (kRunRetroCorrection) {
        if (!ProcessRetroCorrection(ESDtrackHe3, pdg3He, ESDtrackPion, pdgPion, kfpMother)) continue;
      }
      weight = 1;
      CandidateTree->Fill();
      
    }
  }
}
///_____________________________________________________________________________
void AliAnalysisTaskHypertritonKFTreeLocal::ProcessMC(){
  /// Monte Carlo only part of the task
  
  /// Protect against missing MC trees
  AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if(!mcH){
    AliError("No MC Event Handler available");
    return;
  }
  
  fMCEvent = mcH->MCEvent();
  
  if(!fMCEvent){
    AliError("No MC Event, but MC Data required");
    return;
  }
  
  ///----------------------------------------------------------------------------------------------
  /// create a translation table: ESDTrackId(mcTrackId) = ESDTrackIndex, or = -1 if there is no ESDTrack
  ///--------------------------------------------------------------------------------------------------
  vector<Int_t> ESDTrackId(fMCEvent->GetNumberOfTracks(), -1);
  
  /// loop over all ESD tracks and fill ESDTrackId
  for(Int_t i=0; i < fInputEvent->GetNumberOfTracks(); i++) {
    AliESDtrack* track = static_cast<AliESDtrack*>(fInputEvent->GetTrack(i)); /// get a track (type AliESDtrack) from the event
    if(!track) continue;
    if(!PassedBasicTrackQualityCuts(track)) continue;  /// skip track if it fails basic track selection
    
    ESDTrackId[abs(track->GetLabel())] = i;
  }
  
  /// loop over MC all particles
  for (Int_t igen = 0; igen < fMCEvent->GetNumberOfTracks(); igen++){
    AliVParticle *mcParticle = fMCEvent->GetTrack(igen);
    if(!mcParticle) continue;
    
    /// check if we have a hypertriton
    if( abs(mcParticle->PdgCode()) != 1010010030 ) continue;
    /// preselection on rapidity
    RapidityMC = mcParticle->Y();
    if (abs(RapidityMC) > 1.) continue;
    
    pMC = mcParticle->P();
    pTMC = mcParticle->Pt();
    
    //For MC candidate tree
    pxMC = mcParticle->Px();
    pyMC = mcParticle->Py();
    pzMC = mcParticle->Pz();
    
    // calculating weight for events in candidate and generated tree:
    //Blast-Wave Function
    AliPWGFunc *Functions = new AliPWGFunc();
    Functions->SetVarType(AliPWGFunc::kdNdpt);
    //Blast-Wave Fit Parameters
    Double_t massbw = 2.98919;
    Double_t beta = 0.779589;
    Double_t temp = 0.159845;
    Double_t n = 0.120707;
    Double_t norm = 1.12658e+09;
    Double_t massbw2 = 2.98911;
    Double_t beta2 = 0.728534;
    Double_t temp2 = 0.334698;
    Double_t n2 = 0.519438;
    Double_t norm2 = 30080.8;
    //Blast-Wave Distribution
    TF1 *bw_mult = Functions->GetBGBW (massbw, beta, temp, n, norm, "");
    TF1 *bw_mult2 = Functions->GetBGBW (massbw2, beta2, temp2, n2, norm2, "");
    if (CentralityPercentile < 10) {
        weight = bw_mult->Eval(pTMC);
    } else {
        weight = bw_mult2->Eval(pTMC);
    }
    
    Int_t Id3He = -1;
    Int_t IdPion = -1;
    
    /// due to delta electrons a loop over daughters is needed
    AliVParticle* mcDaughter;
    Int_t pdgDaughter = -1;
    Int_t n3He=0;
    Int_t nPion=0;
    
    for (int idaughter = mcParticle->GetDaughterFirst(); idaughter <= mcParticle->GetDaughterLast(); idaughter++) {
      mcDaughter = fMCEvent->GetTrack(idaughter);
      if(!mcDaughter) continue;
      pdgDaughter = mcDaughter->PdgCode();
      
      if (abs(pdgDaughter) == 1000020030) { /// 3He
        n3He++;
        Id3He = idaughter;
      } else if (abs(pdgDaughter) == 211){ /// pion
        nPion++;
        IdPion = idaughter;
      }
    }
    if (n3He > 1 || nPion > 1) continue;
    
    TwoBodyDecayMC(ESDTrackId, Id3He, IdPion);
  }
}
///____________________________________________________________
void AliAnalysisTaskHypertritonKFTreeLocal::TwoBodyDecayMC(const vector<Int_t>& ESDTrackId, Int_t Id3He, Int_t IdPion){
  
  AliVParticle* mcDaughter1 = fMCEvent->GetTrack(Id3He); /// 3He
  AliVParticle* mcDaughter2 = fMCEvent->GetTrack(IdPion); /// pion
  if (!mcDaughter1) return;
  if (!mcDaughter2) return;
  
  Int_t pdgDaughter1 = mcDaughter1->PdgCode();
  if (abs(pdgDaughter1)!= 1000020030) return;
  Int_t pdgDaughter2 = mcDaughter2->PdgCode();
  if (abs(pdgDaughter2)!= 211) return;
  
  /// fill variables for all generated MC hypertritons
  /// calculated decay length
  Double_t prodVertex[3];
  mcDaughter1->XvYvZv(prodVertex);
  
  xSecVertexMC = prodVertex[0];
  ySecVertexMC = prodVertex[1];
  zSecVertexMC = prodVertex[2];
  
  Double_t primVertex[3];
  fInputEvent->GetPrimaryVertex()->GetXYZ(primVertex);
  
  Double_t Difference[3];
  for (Int_t iCord=0; iCord < 3; iCord++) {
    Difference[iCord] = primVertex[iCord]-prodVertex[iCord];
  }
  
  DecayLengthMC = TMath::Sqrt((Difference[0]*Difference[0]) + (Difference[1]*Difference[1]) + (Difference[2]*Difference[2])); /// approximation of the helix curve
  DecayLengthXYMC = TMath::Sqrt((Difference[0]*Difference[0]) + (Difference[1]*Difference[1]));
  
  GeneratedTreeMC->Fill();
  
  /// look for the ESD trackId
  Int_t ESDId1 = ESDTrackId[Id3He];
  Int_t ESDId2 = ESDTrackId[IdPion];
  /// check that both daughters produced an ESD track (with basic track cut)
  if( ESDId1 == -1 || ESDId2 == -1 ) return;
  
  AliESDtrack* ESDTrack1 = static_cast<AliESDtrack*>(fInputEvent->GetTrack(ESDId1));
  AliESDtrack* ESDTrack2 = static_cast<AliESDtrack*>(fInputEvent->GetTrack(ESDId2));
  if (!ESDTrack1) return;
  if (!ESDTrack2) return;
  
  pxHeMC = mcDaughter1->Px();
  pyHeMC = mcDaughter1->Py();
  pzHeMC = mcDaughter1->Pz();
  pHeMC = TMath::Sqrt(pxHeMC * pxHeMC + pyHeMC * pyHeMC + pzHeMC * pzHeMC);
  pTHeMC = TMath::Sqrt(pxHeMC * pxHeMC + pyHeMC * pyHeMC);
  
  pxPionMC =  mcDaughter2->Px();
  pyPionMC =  mcDaughter2->Py();
  pzPionMC =  mcDaughter2->Pz();
  pPionMC = TMath::Sqrt(pxPionMC * pxPionMC + pyPionMC * pyPionMC + pzPionMC * pzPionMC);
  pTPionMC = TMath::Sqrt(pxPionMC * pxPionMC + pyPionMC * pyPionMC);
  
  if (kDoQA) {
    /// QA histograms
    /// Fill histograms for momentum range
    fHistMomHe3->Fill(mcDaughter1->P());
    
    fHistpxTrueRecHe3->Fill(2* ESDTrack1->Px(), pxHeMC-2* ESDTrack1->Px());
    fHistpyTrueRecHe3->Fill(2* ESDTrack1->Py(), pyHeMC-2* ESDTrack1->Py());
    fHistpzTrueRecHe3->Fill(2* ESDTrack1->Pz(), pzHeMC-2* ESDTrack1->Pz());
    
    fHistMomPion->Fill(mcDaughter2->P());
  }
  
  ///    Check identity of daughters
  if ( !Helium3Selection(ESDTrack1) ) return;
  if ( !PionSelection(ESDTrack2) ) return;
  
  KFParticle3LH kfpDaughter1 = CreateKFParticle(ESDTrack1, 2.80839f, 2);
  KFParticle3LH kfpDaughter2 = CreateKFParticle(ESDTrack2, 0.13957f, 1);
  
  // Distance of Daughters Cut
  if (kfpDaughter1.GetDistanceFromParticle(kfpDaughter2) > 3) return;
  
  /// create the kfp mother
  KFParticle3LH kfpMother;
  kfpMother.AddDaughters(kfpDaughter2, kfpDaughter1, 2);

  /// preselect hypertriton candidates before retro-correction
  if (kfpMother.GetNDF()<0) return;
  if (kfpMother.GetChi2()<0) return;
  if (kfpMother.GetChi2()>10000) return; /// protection against infinite
  Float_t DecayLengthUncorr;
  Float_t ErrorDecayLengthUncorr;
  kfpMother.GetDecayLength(DecayLengthUncorr, ErrorDecayLengthUncorr);
  if (DecayLengthUncorr > 60) return;
  
  /// Fill uncorrected variables:
  pxPionUncorr = ESDTrack2->Px();
  pyPionUncorr = ESDTrack2->Py();
  pzPionUncorr = ESDTrack2->Pz();
  pPionUncorr = TMath::Sqrt(pxPionUncorr * pxPionUncorr + pyPionUncorr * pyPionUncorr + pzPionUncorr * pzPionUncorr);
  pTPionUncorr = TMath::Sqrt(pxPionUncorr * pxPionUncorr + pyPionUncorr * pyPionUncorr);
  pxHeUncorr = 2* ESDTrack1->Px();
  pyHeUncorr = 2* ESDTrack1->Py();
  pzHeUncorr = 2* ESDTrack1->Pz();
  pHeUncorr = TMath::Sqrt(pxHeUncorr * pxHeUncorr + pyHeUncorr * pyHeUncorr + pzHeUncorr * pzHeUncorr);
  pTHeUncorr = TMath::Sqrt(pxHeUncorr * pxHeUncorr + pyHeUncorr * pyHeUncorr);
  kfpMother.GetMass(massUncorr, ErrorMassUncorr);
  if ( massUncorr < 2.92 || massUncorr > 3.07 ) return;
  pxUncorr = kfpMother.GetPx();
  pyUncorr = kfpMother.GetPy();
  pzUncorr = kfpMother.GetPz();
  pUncorr = TMath::Sqrt(pxUncorr * pxUncorr + pyUncorr * pyUncorr + pzUncorr * pzUncorr);
  pTUncorr = TMath::Sqrt(pxUncorr * pxUncorr + pyUncorr * pyUncorr);
  
  if (kRunRetroCorrection) {
      if (!ProcessRetroCorrection(ESDTrack1, pdgDaughter1, ESDTrack2, pdgDaughter2, kfpMother)) return;
  }
  
  CandidateTree->Fill();
  
}
///____________________________________________________________
Bool_t AliAnalysisTaskHypertritonKFTreeLocal::ProcessRetroCorrection(AliESDtrack* trackHelium, Int_t pdgHelium, AliESDtrack* trackPion, Int_t pdgPion, KFParticle3LH kfpMother) {
    
    /// 1. get position of secondary vertex
    kfpMother.TransportToDecayVertex(); 
    KFVertex SecondaryVertexUncorr(kfpMother);
    Double_t SV[3];
    SV[0] = SecondaryVertexUncorr.GetX();
    SV[1] = SecondaryVertexUncorr.GetY();
    SV[2] = SecondaryVertexUncorr.GetZ();
    xSecVertexUncorr = SV[0];
    ySecVertexUncorr = SV[1];
    zSecVertexUncorr = SV[2];
    Double_t covSV[6];
    for (Int_t j = 0; j < 6; j++){
        covSV[j] = SecondaryVertexUncorr.GetCovariance(j);
    }
    AliESDVertex* secvertex = new AliESDVertex(SV, covSV, 1, 1, "MySecVertex");
    
    /// 2. retro correction for esd information with sv
    /// get Magnetic field and run number
    magB = fInputEvent->GetMagneticField();
    runN = fInputEvent->GetRunNumber();
    if (!gGeoManager) {
		AliCDBManager::Instance()->SetRaw(1);
		AliCDBManager::Instance()->SetRun(runN);
		AliGeomManager::LoadGeometry();
		AliGeomManager::ApplyAlignObjsFromCDB("GRP ITS TPC TRD");
	}
	if (!TGeoGlobalMagField::Instance()->GetField()) {
		AliGRPManager gm;
		if(!gm.ReadGRPEntry()) {
			printf("Cannot get GRP entry");
		}
		if( !gm.SetMagField() ) {
			printf("Problem with magnetic field setup");
		}
	}
	
	// Position of the primary vertex
    Double_t primVertex[3];
    fInputEvent->GetPrimaryVertex()->GetXYZ(primVertex);
	// Covariance matrix of the vertex position
	// Set all element at 0, this matrix and the quantities related to this matrix are not useful to us
	Double_t *covmatrix = new Double_t[6];
	covmatrix[0] = 0.; covmatrix[1] = 0.; covmatrix[2] = 0.; 
	covmatrix[3] = 0.; covmatrix[4] = 0.; covmatrix[5] = 0.;
	AliESDVertex* primvertex = new AliESDVertex(primVertex, covmatrix, 1, 1, "MyPrimVertex");
    
    xPionTrack = trackPion->GetX();
    yPionTrack = trackPion->GetY();
    xHeTrack = trackHelium->GetX();
    yHeTrack = trackHelium->GetY();
    
    // Propagation to DCA to SV, to get lxHeSV and lxPionSV, but without material correction
    Double_t dztemp2[2], covartemp2[3];
    trackHelium ->PropagateToDCA(secvertex, magB, 250, dztemp2, covartemp2);
    trackPion   ->PropagateToDCA(secvertex, magB, 250, dztemp2, covartemp2);
    
    // x-position of the tracks in local coordinate systems
    lxHeSV = trackHelium->GetX();
    lxPionSV = trackPion->GetX();
    
    // factor 2 for helium, because the track "momenta" in the ESDs are stored as rigidities
    pxHeSV = 2* trackHelium->Px();
    pyHeSV = 2* trackHelium->Py();
    pzHeSV = 2* trackHelium->Pz();
    pHeSV = TMath::Sqrt(pxHeSV * pxHeSV + pyHeSV * pyHeSV + pzHeSV * pzHeSV);
    pTHeSV = TMath::Sqrt(pxHeSV * pxHeSV + pyHeSV * pyHeSV);
    pxPionSV = trackPion->Px();
    pyPionSV = trackPion->Py();
    pzPionSV = trackPion->Pz();
    pPionSV = TMath::Sqrt(pxPionSV * pxPionSV + pyPionSV * pyPionSV + pzPionSV * pzPionSV);
    pTPionSV = TMath::Sqrt(pxPionSV * pxPionSV + pyPionSV * pyPionSV);
    
    // create new KF mother as a test
    KFParticle3LH kfpHeliumSV = CreateKFParticle(trackHelium, 2.80839f, 2);
    KFParticle3LH kfpPionSV   = CreateKFParticle(trackPion, 0.13957f, 1);
 //   KFParticle kfpHeliumSV = CreateKFTrack(trackHelium, pdgHelium);
 //   KFParticle kfpPionSV = CreateKFTrack(trackPion, pdgPion);
    KFParticle3LH kfpMotherSV;
    kfpMotherSV.AddDaughters(kfpPionSV, kfpHeliumSV, 2); // (daughter1, daughter2, constructMethod (0 or 2)
 //   KFParticle kfpMotherSV(kfpHeliumSV,kfpPionSV);
    kfpMotherSV.GetMass(massSV, ErrorMassSV);
    massSVDiff = massSV - massUncorr;
    pxSV = kfpMotherSV.GetPx();
    pySV = kfpMotherSV.GetPy();
    pzSV = kfpMotherSV.GetPz();
    pSV = TMath::Sqrt(pxSV * pxSV + pySV * pySV + pzSV * pzSV);
    pTSV = TMath::Sqrt(pxSV * pxSV + pySV * pySV);
    kfpMotherSV.TransportToDecayVertex(); 
    KFVertex SecondaryVertexSV(kfpMotherSV);
    xSecVertexSV = SecondaryVertexSV.GetX();
    ySecVertexSV = SecondaryVertexSV.GetY();
    zSecVertexSV = SecondaryVertexSV.GetZ();
    xSecVertexSVDiff = SecondaryVertexSV.GetX() - xSecVertexUncorr;
    ySecVertexSVDiff = SecondaryVertexSV.GetY() - ySecVertexUncorr;
    zSecVertexSVDiff = SecondaryVertexSV.GetZ() - zSecVertexUncorr; 
    SecVertexSVDiff = TMath::Sqrt(xSecVertexSVDiff * xSecVertexSVDiff + ySecVertexSVDiff * ySecVertexSVDiff + zSecVertexSVDiff * zSecVertexSVDiff);
    
    // Propagation to DCA to PV, without material correction
    Double_t dztemp[2], covartemp[3];
    trackHelium ->PropagateToDCA(primvertex, magB, 250, dztemp, covartemp);
    trackPion   ->PropagateToDCA(primvertex, magB, 250, dztemp, covartemp);
    
    xPionPV = trackPion->GetX();
    yPionPV = trackPion->GetY();
    xHePV = trackHelium->GetX();
    yHePV = trackHelium->GetY();
    
    // Propagate track from the DCA to the primary vertex to the inner wall of the TPC (R = 85 cm)
    // !!!! WITH MATERIAL CORRECTION !!!!
    lxHeTPC = 0.;
    lxPionTPC = 0.;
    if( !trackHelium->GetXatLabR(85., lxHeTPC, magB, 1) ) return false;
    if( !trackPion->GetXatLabR(85., lxPionTPC, magB, 1) ) return false;
    if( !AliTrackerBase::PropagateTrackTo(trackHelium, lxHeTPC, GetMassForTracking(trackHelium->GetPIDForTracking()), 3, kFALSE, 0.75, kFALSE, kTRUE ) ) return false;
    if( !AliTrackerBase::PropagateTrackTo(trackPion, lxPionTPC, GetMassForTracking(trackPion->GetPIDForTracking()), 3, kFALSE, 0.75, kFALSE, kTRUE ) ) return false;
    
    xPionTPC = trackPion->GetX();
    yPionTPC = trackPion->GetY();
    xHeTPC = trackHelium->GetX();
    yHeTPC = trackHelium->GetY();
    
    // Propagate track from the inner wall of the TPC to the V0 decay point
    // !!!! WITH MATERIAL CORRECTION !!!! and with correct mass hypothesis for the energy-loss correction
    if( !AliTrackerBase::PropagateTrackTo(trackHelium, lxHeSV, GetMassForTracking(AliPID::kHe3), 3, kFALSE, 0.75, kFALSE, kTRUE ) ) return false;
    if( !AliTrackerBase::PropagateTrackTo(trackPion, lxPionSV, GetMassForTracking(AliPID::kPion), 3, kFALSE, 0.75, kFALSE, kTRUE ) ) return false;
    
    xPionSV = trackPion->GetX();
    yPionSV = trackPion->GetY();
    xHeSV = trackHelium->GetX();
    yHeSV = trackHelium->GetY(); 

    /// 3. generate new daughters
    // create new daughters
    KFParticle3LH kfpHeliumCorr = CreateKFParticle(trackHelium, 2.80839f, 2);
    KFParticle3LH kfpPionCorr   = CreateKFParticle(trackPion, 0.13957f, 1);
 //   KFParticle kfpHeliumCorr = CreateKFTrack(trackHelium, pdgHelium);
 //   KFParticle kfpPionCorr = CreateKFTrack(trackPion, pdgPion);
    // Distance of Daughters Cut
    if (kfpHeliumCorr.GetDistanceFromParticle(kfpPionCorr) > 3) return false;
    // Fill daughter variables
    FillDaughterVariables(kfpHeliumCorr,kfpPionCorr);
    // store if tracks had an ITS refit during the track reconstruction
    HasITSRefitHe3 = ((trackHelium->GetStatus() & AliESDtrack::kITSrefit) != 0 );
    HasITSRefitPion = ((trackPion->GetStatus() & AliESDtrack::kITSrefit) != 0 );
    
    /// 4. combine them to a new mother
    KFParticle3LH kfpMotherCorr;
    kfpMotherCorr.AddDaughters(kfpPionCorr, kfpHeliumCorr, 2); // (daughter1, daughter2, constructMethod (0 or 2)
 //   KFParticle kfpMotherCorr(kfpHeliumCorr,kfpPionCorr);
    // Fill corrected variables
    FillPionVariables(trackPion, kfpPionCorr, kfpMotherCorr);      
    FillHe3Variables(trackHelium, kfpHeliumCorr, kfpMotherCorr);
    if (kIsMC){
        TPCnSigma3He = fPIDResponse->NumberOfSigmasTPC(trackHelium, AliPID::kHe3);
        TPCnSigma3H = fPIDResponse->NumberOfSigmasTPC(trackHelium, AliPID::kTriton);
    } else {
        TPCnSigma3He = Bethe(*trackHelium, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe);
        TPCnSigma3H = Bethe(*trackHelium, AliPID::ParticleMass(AliPID::kTriton), 1, fBetheParamsT);
    }
    // Selection on the hypertriton properties
    if (!HypertritonCandidateSelection(kfpMotherCorr)) return false;
    massCorrDiff = mass - massUncorr;
    // careful: transports mother
    FillDistanceToSecondaryVertex(kfpHeliumCorr,kfpPionCorr,kfpMotherCorr);
    
    xSecVertexCorrDiff = xSecVertex - xSecVertexUncorr;
    ySecVertexCorrDiff = ySecVertex - ySecVertexUncorr;
    zSecVertexCorrDiff = zSecVertex - zSecVertexUncorr; 
    SecVertexCorrDiff = TMath::Sqrt(xSecVertexCorrDiff * xSecVertexCorrDiff + ySecVertexCorrDiff * ySecVertexCorrDiff + zSecVertexCorrDiff * zSecVertexCorrDiff);
    if (SecVertexCorrDiff > 5) return false;
    
    // careful! GetDCA includes a propagation to the PV
    DCA3He = GetDCA(trackHelium,"3D");
    DeviationFromPV3He = kfpHeliumCorr.GetDeviationFromVertex(PrimVertex);
    DCA3HeXY = GetDCA(trackHelium,"xy");
    DeviationFromPV3HeXY = kfpHeliumCorr.GetDeviationFromVertexXY(PrimVertex);
    DCAPion = GetDCA(trackPion,"3D");
    DeviationFromPVPion = kfpPionCorr.GetDeviationFromVertex(PrimVertex);
    DCAPionXY = GetDCA(trackPion,"xy");
    DeviationFromPVPionXY = kfpPionCorr.GetDeviationFromVertexXY(PrimVertex); 
    
    return true;
}
///____________________________________________________________
Double_t AliAnalysisTaskHypertritonKFTreeLocal::GetMassForTracking(Int_t fPIDForTracking) {
  Int_t pid = fPIDForTracking;
  if (pid<AliPID::kPion) pid = AliPID::kPion;
  Double_t m = AliPID::ParticleMass(pid);
  return m;
}
///____________________________________________________________
Bool_t AliAnalysisTaskHypertritonKFTreeLocal::PassedEventSelection() {
  if(!fInputEvent){
    AliError("Reconstructed Event not available");
    return false;
  }
  hNumberOfEvents->Fill(0.5);
  
  /// Event selection using the AliEventCuts class
  if (!fEventCuts.AcceptEvent(fInputEvent)) {
    return false;
  }
  hNumberOfEvents->Fill(1.5);
  
  ///Vertex Contributors (only for first checks), should be included in accept events as well
  if ( fInputEvent->GetPrimaryVertex()->GetNContributors() < 2 ) return false;
  if ( fInputEvent->GetPrimaryVertex()->GetZ() < -10 || fInputEvent->GetPrimaryVertex()->GetZ() > 10 ) return false;
  hNumberOfEvents->Fill(2.5);
  
  CentralityPercentile = 300;
  AliMultSelection *MultSelection = (AliMultSelection*) fInputEvent->FindListObject("MultSelection");
  if(!MultSelection) {
    /// If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
    AliWarning("AliMultSelection object not found!");
    PostData(1,fQAList);
    return false;
  }else{
    CentralityPercentile = MultSelection->GetMultiplicityPercentile("V0M");
  }
  hNumberOfEvents->Fill(3.5);
  /// centrality selection 0-90%
  if (CentralityPercentile < 0.0) return false;
  if (CentralityPercentile >= 90.0) return false;
  hNumberOfEvents->Fill(4.5);
  
  ///Time-Range Selection (for LHC18r) (for nor standalone but could be integrated into AliEventCuts: fEventCuts.UseTimeRangeCut(); )
  fTimeRangeCut.InitFromEvent(InputEvent());
  const Bool_t cutThisEvent = fTimeRangeCut.CutEvent(InputEvent());
  if (cutThisEvent) return false;
  hNumberOfEvents->Fill(5.5);
  
  /// Get PID response object
  fPIDResponse = fInputHandler->GetPIDResponse();
  if(!fPIDResponse) {
    AliError("No PID Response found");
    return false;
  }
  hNumberOfEvents->Fill(6.5);
  
  histoEventCentrality->Fill(CentralityPercentile);
  
  return true;
}
///____________________________________________________________
Bool_t AliAnalysisTaskHypertritonKFTreeLocal::PassedBasicTrackQualityCuts (AliESDtrack* track) {
  //Basic Track selection for the analysis
  
  Bool_t Passed = fESDtrackCuts->AcceptTrack(track);
  return Passed;
}
///____________________________________________________________
Bool_t AliAnalysisTaskHypertritonKFTreeLocal::Helium3Selection (AliESDtrack* track)  {
  
  /// PID selection
  if (kIsMC){
      if (abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kHe3)) > 5) return false;
  } else {
      if (abs(Bethe(*track, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe)) > 5) return false;
  }
  return true;
}
///____________________________________________________________
Bool_t AliAnalysisTaskHypertritonKFTreeLocal::PionSelection (AliESDtrack* track)  {
  
  /// Apply stronger preselection on number of crossed rows in the TPC for pions
  if ( (Int_t) track->GetTPCCrossedRows() < 70) return false;
  
  Double_t dEdxSigmaPion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
  /// Select pions with its specific energy loss in the TPC
  if (abs(dEdxSigmaPion) > 5) return false; 
  
  // preselection
  // typical momentum range of pions from hypertriton decays
  if ( track->P() < 0.1 ) return false;
  if ( track->P() > 1.2 ) return false;
  
  return true;
}
///____________________________________________________________
Double_t AliAnalysisTaskHypertritonKFTreeLocal::Bethe(const AliESDtrack& track, Double_t mass, Int_t charge, Double_t* params){
  /// Calculates number of sigma deviation from expected dE/dx in TPC
  /// param track particle track
  /// param mass mass hypothesis of particle
  /// param charge particle charge hypothesis
  /// param params Parameters of Aleph parametrization of Bethe Energy-loss
  Double_t expected = charge*charge*AliExternalTrackParam::BetheBlochAleph(charge*track.GetInnerParam()->GetP()/mass,params[0],params[1],params[2],params[3],params[4]);
  Double_t sigma = expected*params[5];
  if (TMath::IsNaN(expected)) return -999;
  return (track.GetTPCsignal() - expected) / sigma;
}
///____________________________________________________________
Bool_t AliAnalysisTaskHypertritonKFTreeLocal::HypertritonCandidateSelection (KFParticle3LH kfpMother){
  
  /// Select hypertriton candidates and fill variables
  if (kfpMother.GetNDF()<0) return false;
  if (kfpMother.GetChi2()<0) return false;
  if (kfpMother.GetChi2()>10000) return false; /// protection against infinite
  
  Charge = kfpMother.GetQ();
  /// remove reconstructed candidates with wrong charge
  if ( abs(Charge)!= 1 ) return false;  
  
  Bool_t MassCalculated = kfpMother.GetMass(mass, ErrorMass);
  if (MassCalculated !=0) return false;
  /// preselection
  if ( mass < 2.94 || mass > 3.05 ) return false;
  
  /// rapidity selection
  if ( ((kfpMother.E() - kfpMother.Pz()) > 0) && (kfpMother.E() + kfpMother.Pz()) >= 0 ) Rapidity = kfpMother.GetRapidity(); 

  Chi2PerNDF = kfpMother.GetChi2()/kfpMother.GetNDF();
  
  /// fill variables
  DistanceToPV =  kfpMother.GetDistanceFromVertex(PrimVertex);
  DistanceToPVXY =  kfpMother.GetDistanceFromVertexXY(PrimVertex);
  DeviationFromPV = kfpMother.GetDeviationFromVertex(PrimVertex);
  DeviationFromPVXY = kfpMother.GetDeviationFromVertexXY(PrimVertex);
  
  ///preselection on DCA of hypertriton to the primary vertex
  if (DistanceToPV > 5) return false;
  
  px = kfpMother.GetPx();
  py = kfpMother.GetPy();
  pz = kfpMother.GetPz();
  p = TMath::Sqrt(px * px + py * py + pz * pz);
  pT = TMath::Sqrt(px * px + py * py);
  
  if (kIsMC) {
    pxVariance = kfpMother.Covariance(9);
    pyVariance = kfpMother.Covariance(14);
    pzVariance = kfpMother.Covariance(20);
  }
  
  CosPointingAngle = CalculatePointingAngle(kfpMother,PrimVertex);
  
  KFParticle3LH kfpMotherTopo;
  kfpMotherTopo = kfpMother;
  
  kfpMotherTopo.SetProductionVertex(PrimVertex); /// if primary vertex is modified above use the modified version
  
  Bool_t MassCalculatedTopo = kfpMotherTopo.GetMass(massTopo, ErrorMassTopo);
  if (MassCalculatedTopo!=0) {
    massTopo = -999;
    ErrorMassTopo = -999;
  }
  
  Chi2PerNDFTopo = kfpMotherTopo.GetChi2()/kfpMotherTopo.GetNDF();
  
  DistanceToPVTopo =  kfpMotherTopo.GetDistanceFromVertex(PrimVertex);
  DistanceToPVXYTopo =  kfpMotherTopo.GetDistanceFromVertexXY(PrimVertex);
  DeviationFromPVTopo = kfpMotherTopo.GetDeviationFromVertex(PrimVertex);
  DeviationFromPVXYTopo = kfpMotherTopo.GetDeviationFromVertexXY(PrimVertex);
  
  pxTopo = kfpMotherTopo.GetPx();
  pyTopo = kfpMotherTopo.GetPy();
  pzTopo = kfpMotherTopo.GetPz();
  pTopo = TMath::Sqrt(pxTopo * pxTopo + pyTopo * pyTopo + pzTopo * pzTopo);
  pTTopo = TMath::Sqrt(pxTopo * pxTopo + pyTopo * pyTopo);
  
  if (kIsMC) {
    pxVarianceTopo = kfpMotherTopo.Covariance(9);
    pyVarianceTopo = kfpMotherTopo.Covariance(14);
    pzVarianceTopo = kfpMotherTopo.Covariance(20);
  }
  if ( ((kfpMotherTopo.E() - kfpMotherTopo.Pz()) > 0) && (kfpMotherTopo.E() + kfpMotherTopo.Pz()) >= 0 ) RapidityTopo = kfpMotherTopo.GetRapidity();
  
  bool DecayLengthCalculated = kfpMotherTopo.GetDecayLength(DecayLength, ErrorDecayLength); /// returns 0 is correctly calculated
  if (DecayLengthCalculated != 0) {
    DecayLength = -999;
    ErrorDecayLength = -999;
  }
  bool DecayLengthXYCalculated = kfpMotherTopo.GetDecayLengthXY(DecayLengthXY, ErrorDecayLengthXY); /// returns 0 is correctly calculated
  if (DecayLengthXYCalculated != 0){
    DecayLengthXY = -999;
    ErrorDecayLengthXY = -999;
  }  
  
  if (DecayLength > 50 || DecayLength < 0) return false;
  
  CosPointingAngleTopo = CalculatePointingAngle(kfpMotherTopo,PrimVertex);
  /// preselection
  if (CosPointingAngleTopo < 0.95) return false;
  
  return true;
  
}
///____________________________________________________________
Double_t AliAnalysisTaskHypertritonKFTreeLocal::GetDCA (AliESDtrack *track , TString type)  {
  
  /// Calculate DCA to primary vertex
  Double_t impactParameter[2];
  Double_t covarianceMatrix[3];
  if (!track->PropagateToDCA(fInputEvent->GetPrimaryVertex(),fInputEvent->GetMagneticField(),10000,impactParameter,covarianceMatrix)) return -9999;
  
  if(!strncmp(type,"xy",2)) {
    return impactParameter[0];
  } else if (!strncmp(type,"z",1)){
    return impactParameter[1];
  } else if (!strncmp(type,"3D",2)){
    return sqrt(impactParameter[0]*impactParameter[0] + impactParameter[1]*impactParameter[1]);
  }
  
  AliError("### Error ####: Wrong DCA type given \n");
  return -999;
}
///____________________________________________________________
KFParticle AliAnalysisTaskHypertritonKFTreeLocal::CreateKFTrack(AliESDtrack *track, int pdgCode){
  /// Input variables for KF should be float instead of double because not all functions are compatible with both
  
  /// GetTrack parameters
  Double_t trackParameter[6];
  Double_t covMatrix[21];
  
  track->GetXYZ(trackParameter);
  track->GetPxPyPz(&trackParameter[3]);
  track->GetCovarianceXYZPxPyPz(covMatrix);
  
  Int_t Charge = (Int_t) track->Charge();
  if (abs(pdgCode) == 1000020030){
    Charge = Charge*2; /// The exact value seems to be not relevant
    for (int i=3; i<6; i++) {
      trackParameter[i] = trackParameter[i]*2;
    }
    for (int i=6; i<21; i++) {
      covMatrix[i] = covMatrix[i]*2;  /// scale mom space entries of cov matrix by 2
      if (i==9 || i==13 || i==14 || i==18 || i==19 || i==20 ) {
        covMatrix[i] = covMatrix[i]*2;  /// scale mom mom entries of cov matrix by 4
      }
    }
  }
  
  /// Interface to KFParticle
  KFPTrack kfpTrk;
  /// Set the values
  kfpTrk.SetParameters((Float_t) trackParameter[0],(Float_t) trackParameter[1],(Float_t) trackParameter[2],(Float_t) trackParameter[3],(Float_t) trackParameter[4],(Float_t) trackParameter[5]);
  Float_t covF[21];
  for (Int_t i = 0; i<21;i++) { covF[i] = (Float_t) covMatrix[i]; }
  kfpTrk.SetCovarianceMatrix(covF);
  kfpTrk.SetCharge(Charge);
  kfpTrk.SetNDF(1); /// where is this coming from why 1  /// track should be 2?
  
  /// Get Chi2perNDF
  Float_t TrackChi2perNDF = 999;
  Int_t  nClustersTPC = track->GetTPCNcls();
  if ( nClustersTPC > 5) {
    TrackChi2perNDF = (Float_t) track->GetTPCchi2()/Float_t(nClustersTPC - 5);
  }
  kfpTrk.SetChi2(TrackChi2perNDF);
  
  /// Build KFParticle
  KFParticle KFTrk(kfpTrk, pdgCode);
  
  return KFTrk;
}
///____________________________________________________________
KFParticle3LH AliAnalysisTaskHypertritonKFTreeLocal::CreateKFParticle(AliESDtrack *track, float Mass, Int_t Charge) {

  Double_t fP[6];
  track -> GetXYZ(fP);
  track -> PxPyPz(fP+3);
  Int_t fQ = track -> Charge()*TMath::Abs(Charge);
  fP[3] *= TMath::Abs(Charge);
  fP[4] *= TMath::Abs(Charge);
  fP[5] *= TMath::Abs(Charge);

  Double_t pt=1./TMath::Abs(track -> GetParameter()[4]) * TMath::Abs(Charge);
  Double_t cs=TMath::Cos(track -> GetAlpha()), sn=TMath::Sin(track -> GetAlpha());
  Double_t r=TMath::Sqrt((1.-track -> GetParameter()[2])*(1.+track -> GetParameter()[2]));

  Double_t m00=-sn, m10=cs;
  Double_t m23=-pt*(sn + track -> GetParameter()[2]*cs/r), m43=-pt*pt*(r*cs - track -> GetParameter()[2]*sn);
  Double_t m24= pt*(cs - track -> GetParameter()[2]*sn/r), m44=-pt*pt*(r*sn + track -> GetParameter()[2]*cs);
  Double_t m35=pt, m45=-pt*pt*track -> GetParameter()[3];

  m43*=track -> GetSign();
  m44*=track -> GetSign();
  m45*=track -> GetSign();

  const Double_t *cTr = track -> GetCovariance();
  Double_t fC[21];
  fC[0 ] = cTr[0]*m00*m00;
  fC[1 ] = cTr[0]*m00*m10; 
  fC[2 ] = cTr[0]*m10*m10;
  fC[3 ] = cTr[1]*m00; 
  fC[4 ] = cTr[1]*m10; 
  fC[5 ] = cTr[2];
  fC[6 ] = m00*(cTr[3]*m23 + cTr[10]*m43); 
  fC[7 ] = m10*(cTr[3]*m23 + cTr[10]*m43); 
  fC[8 ] = cTr[4]*m23 + cTr[11]*m43; 
  fC[9 ] = m23*(cTr[5]*m23 + cTr[12]*m43)  +  m43*(cTr[12]*m23 + cTr[14]*m43);
  fC[10] = m00*(cTr[3]*m24 + cTr[10]*m44); 
  fC[11] = m10*(cTr[3]*m24 + cTr[10]*m44); 
  fC[12] = cTr[4]*m24 + cTr[11]*m44; 
  fC[13] = m23*(cTr[5]*m24 + cTr[12]*m44)  +  m43*(cTr[12]*m24 + cTr[14]*m44);
  fC[14] = m24*(cTr[5]*m24 + cTr[12]*m44)  +  m44*(cTr[12]*m24 + cTr[14]*m44);
  fC[15] = m00*(cTr[6]*m35 + cTr[10]*m45); 
  fC[16] = m10*(cTr[6]*m35 + cTr[10]*m45); 
  fC[17] = cTr[7]*m35 + cTr[11]*m45; 
  fC[18] = m23*(cTr[8]*m35 + cTr[12]*m45)  +  m43*(cTr[13]*m35 + cTr[14]*m45);
  fC[19] = m24*(cTr[8]*m35 + cTr[12]*m45)  +  m44*(cTr[13]*m35 + cTr[14]*m45); 
  fC[20] = m35*(cTr[9]*m35 + cTr[13]*m45)  +  m45*(cTr[13]*m35 + cTr[14]*m45);

  KFParticle3LH *part = new KFParticle3LH();
  part->Create(fP,fC,fQ,Mass);
  return *part;
}
///____________________________________________________________
KFVertex AliAnalysisTaskHypertritonKFTreeLocal::CreateKFVertex(const AliVVertex* vertex){
  
  /// GetTrack parameters
  Double_t param[6];
  Double_t cov[6];
  
  vertex->GetXYZ(param);
  vertex->GetCovarianceMatrix(cov);
  
  KFPVertex kfpVtx;
  /// Set the values
  Float_t paramF[3] = {(Float_t) param[0],(Float_t) param[1],(Float_t) param[2]};
  kfpVtx.SetXYZ(paramF);
  Float_t covF[6] = {(Float_t) cov[0],(Float_t) cov[1],(Float_t) cov[2],
    (Float_t) cov[3],(Float_t) cov[4],(Float_t) cov[5]};
  kfpVtx.SetCovarianceMatrix(covF);
  KFVertex KFVtx(kfpVtx);
  return KFVtx;
}
///____________________________________________________________
Float_t AliAnalysisTaskHypertritonKFTreeLocal::CalculatePointingAngle(KFParticle3LH KFPart, KFVertex KFVtx){
  
  KFPart.TransportToDecayVertex(); /// After SetProductionVertex the particle is stored at its production vertex but the information at the decay vertex is needed
  
  /// Store position of secondary vertex
  xSecVertex = KFPart.GetX();
  ySecVertex = KFPart.GetY();
  zSecVertex = KFPart.GetZ();
  rSecVertex = TMath::Sqrt(xSecVertex * xSecVertex + ySecVertex * ySecVertex);
  
  xSecVertexVariance = KFPart.Covariance(0);
  ySecVertexVariance = KFPart.Covariance(2);
  zSecVertexVariance = KFPart.Covariance(5);
  
  Double_t v[3];
  v[0] = KFPart.GetX() - KFVtx.GetX();
  v[1] = KFPart.GetY() - KFVtx.GetY();
  v[2] = KFPart.GetZ() - KFVtx.GetZ();
  
  Double_t p[3];
  p[0] = KFPart.GetPx();
  p[1] = KFPart.GetPy();
  p[2] = KFPart.GetPz();
  
  Float_t vnorm3 = TMath::Sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  Float_t pnorm3 = TMath::Sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
  
  Double_t pointAngle   = (v[0]*p[0]+v[1]*p[1]+v[2]*p[2])/(vnorm3*pnorm3); /// cos(pointing angle)
  
  return (Float_t) pointAngle;
  
}
///____________________________________________________________
void AliAnalysisTaskHypertritonKFTreeLocal::FillHe3Variables (AliESDtrack* track, KFParticle3LH KFPart, KFParticle3LH KFMoth){
  
  pxHe = 2* track->Px();
  pyHe = 2* track->Py();
  pzHe = 2* track->Pz();
  pHe = TMath::Sqrt(pxHe * pxHe + pyHe * pyHe + pzHe * pzHe);
  pTHe = TMath::Sqrt(pxHe * pxHe + pyHe * pyHe);
  
  ChargeHe = 2*track->Charge();
  
  TPCMom3He = track->GetTPCmomentum();
  
  HasPointOnITSLayer0He3 = track->HasPointOnITSLayer(0);
  HasPointOnITSLayer1He3 = track->HasPointOnITSLayer(1);
  HasPointOnITSLayer2He3 = track->HasPointOnITSLayer(2);
  HasPointOnITSLayer3He3 = track->HasPointOnITSLayer(3);
  HasPointOnITSLayer4He3 = track->HasPointOnITSLayer(4);
  HasPointOnITSLayer5He3 = track->HasPointOnITSLayer(5);
  
  NCrossedRowsTPC3He = (Int_t) track->GetTPCCrossedRows();
  NPIDClusterTPC3He = (Int_t) track->GetTPCsignalN();
  
  PIDForTrackingHe3 = track->GetPIDForTracking();
  
  if (kIsMC) {
    Double_t covMatrix[21];
    track->GetCovarianceXYZPxPyPz(covMatrix);
    pxHeVariance = covMatrix[9];
    pyHeVariance = covMatrix[14];
    pzHeVariance = covMatrix[20];
  }
}
///____________________________________________________________
void AliAnalysisTaskHypertritonKFTreeLocal::FillPionVariables (AliESDtrack* track, KFParticle3LH KFPart, KFParticle3LH KFMoth){
  
  pxPion = track->Px();
  pyPion = track->Py();
  pzPion = track->Pz();
  pPion = TMath::Sqrt(pxPion * pxPion + pyPion * pyPion + pzPion * pzPion);
  pTPion = TMath::Sqrt(pxPion * pxPion + pyPion * pyPion);

  ChargePion = track->Charge();
  
  TPCMomPion = track->GetTPCmomentum();
  TPCnSigmaPion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
  
  HasPointOnITSLayer0Pion = track->HasPointOnITSLayer(0);
  HasPointOnITSLayer1Pion = track->HasPointOnITSLayer(1);
  HasPointOnITSLayer2Pion = track->HasPointOnITSLayer(2);
  HasPointOnITSLayer3Pion = track->HasPointOnITSLayer(3);
  HasPointOnITSLayer4Pion = track->HasPointOnITSLayer(4);
  HasPointOnITSLayer5Pion = track->HasPointOnITSLayer(5);
  
  NCrossedRowsTPCPion = (Int_t) track->GetTPCCrossedRows();
  NPIDClusterTPCPion =  (Int_t) track->GetTPCsignalN();
  
  PIDForTrackingPion = track->GetPIDForTracking();
  
  if (kIsMC) {
    Double_t covMatrix[21];
    track->GetCovarianceXYZPxPyPz(covMatrix);
    pxPionVariance = covMatrix[9];
    pyPionVariance = covMatrix[14];
    pzPionVariance = covMatrix[20];
  }
}
///____________________________________________________________
void AliAnalysisTaskHypertritonKFTreeLocal::FillDaughterVariables (KFParticle3LH kfpDaughter1, KFParticle3LH kfpDaughter2){
  
  DistanceOfDaughters = kfpDaughter1.GetDistanceFromParticle(kfpDaughter2);
  DeviationOfDaughters = kfpDaughter1.GetDeviationFromParticle(kfpDaughter2);
  
  DistanceOfDaughtersXY = kfpDaughter1.GetDistanceFromParticleXY(kfpDaughter2);
  DeviationOfDaughtersXY = kfpDaughter1.GetDeviationFromParticleXY(kfpDaughter2);
}
///____________________________________________________________
void AliAnalysisTaskHypertritonKFTreeLocal::FillDistanceToSecondaryVertex(KFParticle3LH kfpHelium, KFParticle3LH kfpPion, KFParticle3LH kfpMother){
  
  kfpMother.TransportToDecayVertex(); /// After SetProductionVertex the particle is stored at its production vertex but the information at the decay vertex is needed
  KFVertex SecondaryVertex(kfpMother);
  
  DistanceToSecVertHe = kfpHelium.GetDistanceFromVertex(SecondaryVertex);
  DeviationToSecVertHe = kfpHelium.GetDeviationFromVertex(SecondaryVertex);
  
  DistanceToSecVertPion = kfpPion.GetDistanceFromVertex(SecondaryVertex);
  DeviationToSecVertPion = kfpPion.GetDeviationFromVertex(SecondaryVertex);
}
