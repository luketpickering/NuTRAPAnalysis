#ifndef __TRANSVERSITY_VARIABLE_OBJECTS_SEEN__
#define __TRANSVERSITY_VARIABLE_OBJECTS_SEEN__
#include "TObject.h"
#include "TLorentzVector.h"

#include "TransversityUtils.hxx"

struct PartStruct {
  PartStruct(){
    Reset();
  }
  Int_t PDG;
  Double_t Momentum;
  TLorentzVector FourMomentum;
  void Reset(){
    PDG = 0;
    Momentum = 0;
    FourMomentum = TLorentzVector(0,0,0,0);
  }
};

struct MuonProtonTransversity : public TObject {
private:
  //Transients
  PartStruct Muon; //!
  PartStruct MuonNeutrino; //!
  PartStruct StruckNucleon; //!
  PartStruct HMProton; //!
  PartStruct HMCPion; //!
  PartStruct FirstProton; //!
  PartStruct HMTrackable; //!

  bool IsInGev; //!

  Int_t SmearingMag; //!
  bool DoSmear = false; //!
public:
  MuonProtonTransversity(bool InGeV=true, Int_t MomentumThresh_MeV=10,
  Int_t SmearingMag=0, bool ds=false);

  Double_t ReconNuEnergy;
  Double_t ReconTargetMass;

//DeltaPhiT
  Double_t DeltaPhiT_HMProton;
  Double_t DeltaPhiT_FirstProton;

  Double_t DeltaPhiT_HMProton_deg;
  Double_t DeltaPhiT_FirstProton_deg;
  Double_t DeltaPhiT_HMTrackable_deg;
  Double_t DeltaPhiT_HMTrackable_GSmear_deg;
  Double_t DeltaPhiT_HMTrackable_CSmear_deg;
  Double_t GSmear_mu_deg;
  Double_t GSmear_HMTrackable_deg;
  Double_t CSmear_mu_deg;
  Double_t CSmear_HMTrackable_deg;
//DeltaPt
  TVector3 DeltaPT_HMProton_MeV;
  TVector3 DeltaPT_FirstProton_MeV;
//DeltaAlphat
  Double_t DeltaAlphaT_HMProton;
  Double_t DeltaAlphaT_FirstProton;

  Double_t DeltaAlphaT_HMProton_deg;
  Double_t DeltaAlphaT_FirstProton_deg;
//DeltaP_TT
  Double_t DeltaP_TT;
  Int_t DeltaP_TT_PionPDG;
//Generator reaction code
  Int_t NeutConventionReactionCode;

//Final state contents
  Int_t NFinalStateParticles;

  Int_t NProtons;
  Int_t NAboveThresholdProtons;
  Int_t NGamma;
  Int_t NAboveThresholdGamma;
  Int_t NNeutrons;

  Int_t NPiPlus;
  Int_t NAboveThresholdPiPlus;
  Int_t NPiZero;
  Int_t NPiMinus;
  Int_t NAboveThresholdPiMinus;
  Int_t NPions;
  Int_t NChargedPions;
  Int_t NAboveThresholdChargedPions;

  Int_t NAboveThresholdTrackable;

  Int_t NOtherParticles;

//Selected final state protperties
  Int_t MuonPDG;
  Int_t HMProtonPDG;
  Int_t FirstProtonPDG;
  Int_t HMTrackablePDG;

  TVector3 MuonDirection;
  Double_t MuonMomentum_MeV;

  TVector3 StruckNucleonDirection;
  Double_t StruckNucleonMomentum_MeV;
  Int_t StruckNucleonPDG;

  TVector3 HMProtonDirection;
  TVector3 FirstProtonDirection;
  TVector3 HMTrackableDirection;

  Double_t HMProtonMomentum_MeV;
  Double_t FirstProtonMomentum_MeV;
  Double_t HMTrackableMomentum_MeV;

  TVector3 MuonPt_MeV;
  TVector3 HMProtonPt_MeV;
  TVector3 FirstProtonPt_MeV;
  TVector3 HMTrackablePt_MeV;
  TVector3 StruckNucleonPt_MeV;
//Neutrino properties
  Int_t IncNeutrinoPDG;
  TLorentzVector IncNeutrinoMmtm;
//Target Properties
  Int_t TargetPDG;
  Int_t TargetZ;

  Int_t MomentumThresh_MeV;

private:
  void HandleProton(TLorentzVector &StdHepPTLV,
    Double_t &StdHepP3Mod);
  void HandleCPion(TLorentzVector &StdHepPTLV,
    Double_t &StdHepP3Mod, Int_t pdg);
  void HandleHMTrackable(TLorentzVector &StdHepPTLV,
    Double_t &StdHepP3Mod, Int_t PDG);
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

  ClassDef(MuonProtonTransversity,1);
};

struct PionProductionTransversity : public TObject {

private:
  //Transients
  PartStruct MuonNeutrino; //!
  PartStruct Muon; //!
  PartStruct HMPiPlus; //!
  PartStruct FirstPiPlus; //!
  bool IsInGev; //!

public:
  PionProductionTransversity(bool InGeV=true, Int_t MomentumThresh_MeV=10);

  Double_t ReconNuEnergy;
  Double_t ReconTargetMass;

  Double_t DeltaPhiT_FirstPiPlus;
  Double_t DeltaPhiT_HMPiPlus;

  Double_t DeltaPhiT_FirstPiPlus_deg;
  Double_t DeltaPhiT_HMPiPlus_deg;

  TVector3 DeltaPT_HMPiPlus_MeV;
  TVector3 DeltaPT_FirstPiPlus_MeV;
  Double_t DeltaAlphaT_HMPiPlus;
  Double_t DeltaAlphaT_FirstPiPlus;
  Double_t DeltaAlphaT_HMPiPlus_deg;
  Double_t DeltaAlphaT_FirstPiPlus_deg;

  Int_t NeutConventionReactionCode;

  Int_t NFinalStateParticles;
  Int_t NProtons;
  Int_t NAboveThresholdProtons;
  Int_t NGamma;
  Int_t NAboveThresholdGamma;
  Int_t NNeutrons;

  Int_t NPiPlus;
  Int_t NAboveThresholdPiPlus;
  Int_t NPiZero;
  Int_t NPiMinus;
  Int_t NAboveThresholdPiMinus;
  Int_t NPions;
  Int_t NChargedPions;
  Int_t NAboveThresholdChargedPions;

  Int_t NOtherParticles;

  Int_t MuonPDG;
  Int_t HMPiPlusPDG;
  Int_t FirstPiPlusPDG;

  TVector3 MuonDirection;
  TVector3 HMPiPlusDirection;
  TVector3 FirstPiPlusDirection;

  Double_t MuonMomentum_MeV;
  Double_t HMPiPlusMomentum_MeV;
  Double_t FirstPiPlusMomentum_MeV;

  TVector3 MuonPt_MeV;
  TVector3 HMPiPlusPt_MeV;
  TVector3 FirstPiPlusPt_MeV;

  Int_t IncNeutrinoPDG;
  TLorentzVector IncNeutrinoMmtm;

  Int_t TargetPDG;
  Int_t TargetZ;

  Int_t MomentumThresh_MeV;

private:
  void HandlePiPlus(TLorentzVector &StdHepPTLV,
    Double_t &StdHepP3Mod);

public:
  bool SetNeutConventionReactionCode(int rc){
    NeutConventionReactionCode = rc;
    return (rc==11)||(rc==16);
  }

  bool HandleStdHepParticle(UInt_t &StdHepPosition,
                            int &StdHepPdg,
                            int &StdHepStatus,
                            Double_t * &StdHepP4);
  void Finalise();
  void Reset();

  ClassDef(PionProductionTransversity,1);

};

#endif
