#ifndef __TRANSVERSITY_VARIABLE_OBJECTS_SEEN__
#define __TRANSVERSITY_VARIABLE_OBJECTS_SEEN__
#include "TObject.h"
#include "TLorentzVector.h"

#include "TransversityUtils.hxx"

static const Int_t kNThreshMax = 10;

struct PartStruct {
  PartStruct(){
    Reset();
  }
  Int_t PDG;
  Double_t Momentum;
  TLorentzVector FourMomentum;
  UInt_t StdHepPosition;
  void Reset(){
    PDG = 0;
    Momentum = 0;
    FourMomentum = TLorentzVector(0,0,0,0);
    StdHepPosition = 0;
  }
};

struct TransversityVars : public TObject {
private:

  bool IsInGev; //!

  Int_t SmearingMag_PhiDeg; //!
  Int_t SmearingMag_MuMmtmMeV; //!
  Int_t SmearingMag_PMmtmMeV; //!
  bool DoSmear = false; //!

public:

  TransversityVars(bool InGeV=true, Double_t TargetBE_MeV=0xdeadbeef,
    Int_t NThresh=0, Int_t* Threshs_MeV=0,
    Int_t SmearingMag_PhiDeg=0, Int_t SmearingMag_MuMmtmMeV=0,
    Int_t SmearingMag_PMmtmMeV=0, bool ds=false,
    TString generatorName="");

//******************************************************************************
//                     Pertinent Particle Properties
//******************************************************************************

  Int_t NThresh;
  Int_t* Threshs_MeV; //[NThresh]

  TString GeneratorName;

//******************************************************************************
//                     Event Properties
//******************************************************************************

//Generator reaction code
  Int_t NeutConventionReactionCode;
  Double_t ReconNuEnergy;
  Double_t ReconTargetMass;

  Double_t ERec_CCQE;
  Double_t ERec_E_PMu;

//******************************************************************************
//                     Pertinent Particle Properties
//******************************************************************************

//Neutrino
  Int_t IncNeutrino_PDG;
  TLorentzVector IncNeutrino_4Mom_MeV;

//Struck Nucleon
  Int_t StruckNucleonPDG;
  TLorentzVector StruckNucleon_4Mom_MeV;

//Muon
  Int_t Muon_PDG;
  TLorentzVector Muon_4Mom_MeV;
  TVector3 Muon_Pt_MeV;
  Double_t Muon_MomMag_MeV_GSmear;

//First Proton
  Int_t FirstProton_PDG;
  TLorentzVector FirstProton_4Mom_MeV;
  Int_t FirstProton_StdHepPosition;

//Highest Momentum Proton
  Int_t HMProton_PDG;
  TLorentzVector HMProton_4Mom_MeV;
  Double_t HMProton_MomMag_MeV_GSmear;
  Int_t HMProton_StdHepPosition;

//Highest Momentum Trackable
  Int_t HMTrackable_PDG;
  TLorentzVector HMTrackable_4Mom_MeV;
  Double_t HMTrackable_MomMag_MeV_GSmear;

//******************************************************************************
//                       'Verse Variable Values
//******************************************************************************

//DeltaPhiT
  Double_t DeltaPhiT_HMProton_deg;
  Double_t DeltaPhiT_FirstProton_deg;
  Double_t DeltaPhiT_HMTrackable_deg;
  Double_t DeltaPhiT_HMTrackable_GSmear_deg;

//DeltaPt
  TVector3 DeltaPT_HMProton_MeV;
  TVector3 DeltaPT_FirstProton_MeV;

//DeltaAlphat
  Double_t DeltaAlphaT_HMProton_deg;
  Double_t DeltaAlphaT_FirstProton_deg;
//DeltaP_TT
  Int_t DeltaP_TT_PionPDG;
  Double_t DeltaP_TT;

//******************************************************************************
//                       Smeared Sample States
//******************************************************************************

  Double_t Muon_phi_GSmear_deg;
  Double_t HMTrackable_phi_GSmear_deg;

//******************************************************************************
//                       Final State Particles
//******************************************************************************

  Int_t NFinalStateParticles;

  Int_t NProtons;
  Int_t NGammas;
  Int_t NNeutrons;
  Int_t NPiPlus;
  Int_t NPiZero;
  Int_t NPiMinus;
  Int_t NPions;
  Int_t NChargedPions;
  Int_t NOtherParticles;

  Int_t NAboveThresholdProtons[kNThreshMax]; //[NThresh]
  Int_t NAboveThresholdGammas[kNThreshMax]; //[NThresh]
  Int_t NAboveThresholdPiPlus[kNThreshMax]; //[NThresh]
  Int_t NAboveThresholdPiMinus[kNThreshMax]; //[NThresh]
  Int_t NAboveThresholdChargedPions[kNThreshMax]; //[NThresh]
  Int_t NAboveThresholdTrackable[kNThreshMax]; //[NThresh]
  Int_t NAboveThresholdNeutrons[kNThreshMax]; //[NThresh]
  Int_t NAboveThresholdPiZero[kNThreshMax]; //[NThresh]
  Int_t NAboveThresholdNeutrals[kNThreshMax]; //[NThresh]
  Int_t NAboveThresholdExotic[kNThreshMax]; //[NThresh]

//******************************************************************************
//                       Tangible Target Traits
//******************************************************************************

  Int_t TargetPDG;
  Int_t TargetZ;

//******************************************************************************
//                       Others and Transients
//******************************************************************************

  TLorentzVector CCQE_deltaP_MeV;


  //Transients
  PartStruct Muon; //!
  PartStruct MuonNeutrino; //!
  PartStruct StruckNucleon; //!
  PartStruct HMProton; //!
  PartStruct HMCPion; //!
  PartStruct FirstProton; //!
  PartStruct HMTrackable; //!
  Double_t TargetBE_MeV; //!


//******************************************************************************
//******************************************************************************

private:
  void HandleProton(TLorentzVector &StdHepPTLV,
    Double_t &StdHepP3Mod, UInt_t &StdHepPosition);
  void HandleCPion(TLorentzVector &StdHepPTLV,
    Double_t &StdHepP3Mod, Int_t pdg);
  void HandleHMTrackable(TLorentzVector &StdHepPTLV,
    Double_t &StdHepP3Mod, Int_t PDG);
  Int_t PartIsAboveThresh(TLorentzVector FourMom);
public:
  bool SetNeutConventionReactionCode(int rc){
    NeutConventionReactionCode = rc;
    return (rc==1);
  }

  bool HandleStdHepParticle(UInt_t &StdHepPosition,
                            Int_t &StdHepPdg,
                            Int_t &StdHepStatus,
                            Double_t * &StdHepP4);
  void HandleStruckNucleon(TLorentzVector &StdHepPTLV,
    Double_t &StdHepP3Mod, Int_t pdg);
  void Finalise();
  void Reset();

  ClassDef(TransversityVars,1);
};

#endif
