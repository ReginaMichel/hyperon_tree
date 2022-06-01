#if !defined (__CINT__) || defined (__CLING__)
R__ADD_INCLUDE_PATH($ALICE_ROOT)
#include "ANALYSIS/macros/AddTaskPIDResponse.C"

R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
#include "OADB/macros/AddTaskPhysicsSelection.C"
#include "OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"

R__LOAD_LIBRARY(../AliAnalysisTaskHypertritonKFTreeLocal_cxx.so)

#include "AliAnalysisAlien.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "TStopwatch.h"

R__ADD_INCLUDE_PATH($KFPARTICLE_ROOT)
#include "../AliAnalysisTaskHypertritonKFTreeLocal.h"
#include "../AddHypertritonKFTreeLocal.C"

#endif

//=====================================================================================================================================================


//Definitions
#define ALIPHYSICS_VER  "vAN-20211105_ROOT6-1"
#define GRIDWorkingDir  "HYPERTRITON_TREE"
#define AnalysisMacro   "Analysis"
#define AnalysisTask    "AliAnalysisTaskHypertritonKFTreeLocal"
#define TTL             20000
#define nRunsPerMaster  10


//Functions
AliAnalysisGrid *CreateAlienHandler (Int_t iMC, const char *mode, Bool_t merge);
void LoadAnalysisTask (Int_t iMC, AliAnalysisManager *mgr);
void EventHandler (AliAnalysisManager *mgr);
void SetPath (Int_t iMC, AliAnalysisAlien *alien);
void SetInputRuns (AliAnalysisAlien *alien, const char *mode, Int_t iMC);
void SetAdditionalLibraries (AliAnalysisAlien *alien);
void LoadPhysicsSelection();
void LoadPIDResponse ();
void LoadCentrality ();
void LoadKF();


//______________________________________________________________________________________________________________________________________________________
void LoadAnalysisTask (Int_t iMC, AliAnalysisManager *mgr)  {
    
    UInt_t triggerMask = (AliVEvent::kINT7|AliVEvent::kCentral|AliVEvent::kSemiCentral);
        
    Bool_t Run2Body = true;   /// if false  2-body output still created but empty
    Bool_t Run3Body = false;   /// if false  3-body output still created but empty
    Bool_t OnlyBackground = false;
    Bool_t IsMC = true;     /// If true task is running using MC Truth information /// false means running as inout would be data
    Bool_t DoQA = true;
    
    gROOT->LoadMacro("AliAnalysisTaskHypertritonKFTreeLocal.cxx++g");
    gROOT->LoadMacro("AddHypertritonKFTreeLocal.C");
    AliAnalysisTaskHypertritonKFTreeLocal *task = AddHypertritonKFTreeLocal(triggerMask, Run2Body, Run3Body, IsMC, DoQA, OnlyBackground);

}
//______________________________________________________________________________________________________________________________________________________
void runAnalysis (Int_t iMC=0, const char *mode="termiante", Bool_t merge=kTRUE)  {
// void runAnalysis (Int_t iMC=0, const char *mode="full", Bool_t merge=kTRUE)  {
    
    //Grid Connection
    TGrid::Connect("alien://");
    
    //Alien Handler
    AliAnalysisGrid *alienHandler = CreateAlienHandler (iMC,mode,merge);
    if (!alienHandler) return;
    
    //Analysis Manager
    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisManager");
    mgr->SetGridHandler(alienHandler);
    
    //Event Handler, PID, Physics Selection & Centrality
    EventHandler(mgr);
    LoadPhysicsSelection ();
    LoadCentrality();
    LoadPIDResponse();
    LoadKF ();

    //Analysis Task
    LoadAnalysisTask (iMC,mgr);
    
    //Start Analysis
    if (!mgr->InitAnalysis()) return;
    mgr->PrintStatus();
    mgr->StartAnalysis("grid");
}
//______________________________________________________________________________________________________________________________________________________
AliAnalysisGrid *CreateAlienHandler (Int_t iMC, const char *mode, Bool_t merge)  {
    
    //Alien Handler
    AliAnalysisAlien *alien = new AliAnalysisAlien();
    SetAdditionalLibraries (alien);
    alien->SetOverwriteMode();
    alien->SetCheckCopy(kFALSE);
    alien->SetRunMode(mode);
    alien->SetNtestFiles(10);
    alien->SetAPIVersion("V1.1x");
    alien->SetAliPhysicsVersion(ALIPHYSICS_VER);
    alien->AddIncludePath("$ALICE_PHYSICS/include");
    alien->AddIncludePath("$ALICE_ROOT/include");
    alien->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include -I$KFPARTICLE_ROOT -I$KFPARTICLE_ROOT/include");

    SetPath (iMC,alien);
    alien->SetDataPattern("*AliESDs.root");
    SetInputRuns (alien,mode,iMC);
    alien->SetNrunsPerMaster(nRunsPerMaster);
    alien->SetGridWorkingDir (Form("%s_%d",GRIDWorkingDir,iMC));
    alien->SetGridOutputDir("OUTPUT");
    alien->SetAnalysisSource(Form("%s.cxx",AnalysisTask));
    alien->SetAdditionalLibs(Form("%s.cxx %s.h libpythia6_4_21.so",AnalysisTask,AnalysisTask));
    alien->SetMergeViaJDL(merge);
    alien->SetMaxMergeStages(2);
    alien->SetAnalysisMacro(Form("%s_%d.C",AnalysisMacro,iMC));
    alien->SetSplitMaxInputFileNumber(50);
    alien->SetMasterResubmitThreshold(90);
    alien->SetTTL(TTL);
    alien->SetExecutable(Form("%s_%d.sh",AnalysisMacro,iMC));
    alien->SetInputFormat("xml-single");
    alien->SetJDLName(Form("%s_%d.jdl",AnalysisMacro,iMC));
    alien->SetMergeExcludes("EventStat_temp.root");
    alien->SetPrice(1);
    alien->SetSplitMode("se");
    return alien;
}
//______________________________________________________________________________________________________________________________________________________________
void LoadKF()  {
    
    gInterpreter->ProcessLine(".include $ROOTSYS/include");
    gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
    gInterpreter->ProcessLine(".include $KFPARTICLE_ROOT/include");

}
//______________________________________________________________________________________________________________________________________________________
void EventHandler (AliAnalysisManager *mgr)  {
    
    //ESD Input Handler
    AliESDInputHandler *esdH = new AliESDInputHandler();
    mgr->SetInputEventHandler(esdH);

    //MC Event Handler
    AliMCEventHandler *mcHandler = new AliMCEventHandler();
    mcHandler->SetReadTR(false);
    mgr->SetMCtruthEventHandler(mcHandler);
}
//______________________________________________________________________________________________________________________________________________________
void LoadPIDResponse ()  {
    
    Bool_t isMC=kTRUE;
    AliAnalysisTaskPIDResponse *pidTask = AddTaskPIDResponse(isMC);
}
//______________________________________________________________________________________________________________________________________________________
void LoadPhysicsSelection()  {
    
    Bool_t isMC = kTRUE;
    AliPhysicsSelectionTask *PhySel = AddTaskPhysicsSelection(isMC);
}
//______________________________________________________________________________________________________________________________________________________
void LoadCentrality ()  {

    Bool_t lCalibration = kFALSE;
    AliMultSelectionTask *MultSelTask = AddTaskMultSelection(lCalibration);
}
//______________________________________________________________________________________________________________________________________________________
void SetPath (Int_t iMC, AliAnalysisAlien *alien)  {
    
    const char *path_mc_production;
   // if (iMC==0) path_mc_production=("/alice/sim/2020/LHC20d2a/"); pass1 
    if (iMC==1) path_mc_production=("/alice/sim/2020/LHC20d2a/"); // 0-10%
    if (iMC==2) path_mc_production=("/alice/sim/2020/LHC20d2b/"); // 10-50%
    if (iMC==3) path_mc_production=("/alice/sim/2020/LHC20d2c/"); // 50%-90%
    
                                                                  // pass3
    if (iMC==4) path_mc_production=("/alice/sim/2020/LHC20g7a/"); // 0-10%
    if (iMC==5) path_mc_production=("/alice/sim/2020/LHC20g7b/"); // 10-50%
    if (iMC==6) path_mc_production=("/alice/sim/2020/LHC20g7c/"); // 50%-90%

    
    alien->SetGridDataDir(path_mc_production);
}
//______________________________________________________________________________________________________________________________________________________
void SetInputRuns (AliAnalysisAlien *alien, const char *mode, Int_t iMC)  {
    
//     //Run List
     Int_t run[] = {
            297595, 297590, 297588, 297558, 297544, 297542, 297541, 297540, 297537, 297512, 297483, 297481, 297479, 297452, 297451, 297450, 297446, 297442, 297441, 297415, 297414, 297413, 297406, 297405, 297380, 297379, 297372, 297367, 297366, 297363, 297336, 297335, 297333, 297332, 297317, 297315, 297312, 297311, 297310, 297278, 297222, 297221, 297219, 297218, 297196, 297195, 297194, 297193, 297133, 297132, 297129, 297128, 297124, 297123, 297119, 297118, 297117, 297085, 297035, 297031, 297029, 296966, 296941, 296938, 296935, 296934, 296932, 296931, 296930, 296903, 296900, 296899, 296894, 296890, 296852, 296851, 296850, 296849, 296848, 296839, 296838, 296836, 296835, 296799, 296794, 296793, 296790, 296787, 296786, 296785, 296784, 296781, 296752, 296750, 296749, 296694, 296693, 296691, 296690, 296623, 296622, 296621, 296619, 296618, 296616, 296615, 296594, 296553, 296552, 296551, 296550, 296549, 296548, 296547, 296516, 296512, 296511, 296510, 296509, 296472, 296433, 296424, 296423, 296420, 296419, 296415, 296414, 296383, 296381, 296380, 296379, 296378, 296377, 296376, 296375, 296312, 296309, 296307, 296304, 296303, 296280, 296279, 296273, 296270, 296269, 296247, 296246, 296244, 296243, 296242, 296241, 296240, 296198, 296197, 296196, 296195, 296194, 296192, 296191, 296143, 296142, 296135, 296134, 296133, 296132, 296123, 296074, 296068, 296066, 296065, 296063, 296062, 296060, 296016, 295947, 295945, 295943, 295942, 295941, 295937, 295936, 295913, 295910, 295909, 295908, 295881, 295861, 295860, 295859, 295856, 295855, 295854, 295853, 295831, 295829, 295826, 295825, 295822, 295819, 295818, 295816, 295791, 295788, 295786, 295763, 295762, 295759, 295758, 295755, 295754, 295725, 295723, 295721, 295719, 295718, 295717, 295714, 295712, 295677, 295676, 295675, 295673, 295671, 295668, 295667, 295666, 295615, 295612, 295611, 295610, 295589, 295588, 295587, 295586, 295585 };
           


       Int_t nRuns = sizeof(run)/sizeof(Int_t);
       for ( Int_t i=0; i<nRuns; i++ )
           alien->AddRunNumber(run[i]);
}
//______________________________________________________________________________________________________________________________________________________
void SetAdditionalLibraries (AliAnalysisAlien *alien)  {
    
    alien->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS -I$ALICE_PHYSICS/include -I$ALICE_ROOT/PYTHIA6 -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS  -I$ALICE_PHYSICS/OADB -I$ALICE_PHYSICS/PWGHF -I$ALICE_PHYSICS/PWGHF/base -I$ALICE_PHYSICS/PWGHF/vertexingHF -I$ALICE_PHYSICS/PWG/FLOW/Base -I$ALICE_PHYSICS/PWG/FLOW/Tasks -I$ALICE_PHYSICS/JETAN -I$ALICE_PHYSICS/PWG/Tools -g");

    alien->SetAdditionalLibs("libCDB.so libSTEER.so  libCORRFW.so libPWGflowBase.so libPWGflowTasks.so libGui.so libProof.so libMinuit.so libXMLParser.so libRAWDatabase.so libRAWDatarec.so libCDB.so libSTEERBase.so libSTEER.so libTPCbase.so libTOFbase.so libTOFrec.so libTRDbase.so libVZERObase.so libVZEROrec.so libT0base.so libT0rec.so libpythia6.so libEGPythia6.so libAliPythia6.so");
    
}
//______________________________________________________________________________________________________________________________________________________
