#ifndef __TRANSVERSITY_VARIABLE_OBJECTS_SEEN__
#define __TRANSVERSITY_VARIABLE_OBJECTS_SEEN__

#include <exception>

#include "TObject.h"
#include "TLorentzVector.h"

#include "TransversityUtils.hxx"

static const Int_t kNThreshMax = 10;
static const Int_t kMaxFSMomenta = 20;

struct ENotImplemented : public std::exception {
  std::string Msg;
  ENotImplemented(std::string msg=""){
    Msg = msg;
  }
  virtual const char* what() const throw(){
    return Msg.c_str();
  }
};

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

struct TransversityVarsB : public TObject {

  virtual bool SetNeutConventionReactionCode(int) = 0;
  virtual bool HandleStdHepParticle(UInt_t &StdHepPosition,
                              Int_t &StdHepPdg,
                              Int_t &StdHepStatus,
                              Double_t * &StdHepP4) = 0;
  virtual void HandleRescat(Int_t PDG, Int_t RescatCode) = 0;
  virtual void HandleStruckNucleon(TLorentzVector &StdHepPTLV,
      Double_t &StdHepP3Mod, Int_t pdg) = 0;
  virtual void Finalise() = 0;
  virtual void Reset() = 0;

  virtual Int_t GetIncNeutrino_PDG() = 0;
  virtual Int_t GetStruckNucleonPDG() = 0;
  virtual Int_t GetNeutConventionReactionCode() = 0;
};

struct TransversityVars : public TransversityVarsB {
private:

  bool IsInGev; //!

public:

  TransversityVars(bool InGeV=true, Double_t TargetBE_MeV=0xdeadbeef,
    Int_t NThresh=0, Int_t* Threshs_MeV=0,
    TString generatorName="");
  ~TransversityVars();
//******************************************************************************
//                     Pertinent Particle Properties
//******************************************************************************

  Int_t NThresh; //!
  Int_t* Threshs_MeV; //! //[NThresh]

  TString GeneratorName; //!

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

//First Proton
  Int_t FirstProton_PDG;
  TLorentzVector FirstProton_4Mom_MeV;
  Int_t FirstProton_StdHepPosition;

//Highest Momentum Proton
  Int_t HMProton_PDG;
  TLorentzVector HMProton_4Mom_MeV;
  Int_t HMProton_StdHepPosition;

//Highest Momentum Charged Pion
  Int_t HMCPion_PDG;
  TLorentzVector HMCPion_4Mom_MeV;

//Highest Momentum Trackable
  Int_t HMTrackable_PDG;
  TLorentzVector HMTrackable_4Mom_MeV;

//******************************************************************************
//                       'Verse Variable Values
//******************************************************************************

//DeltaPhiT
  Double_t DeltaPhiT_HMProton_deg;
  Double_t DeltaPhiT_FirstProton_deg;
  Double_t DeltaPhiT_HMTrackable_deg;

//DeltaPt
  TVector3 DeltaPT_HMProton_MeV;
  TVector3 DeltaPT_FirstProton_MeV;

//DeltaPTotal
  TLorentzVector DeltaPTotal_HMProton_MeV;
  TLorentzVector DeltaPTotal_FirstProton_MeV;

//DeltaPNuclearEffect
  TLorentzVector DeltaPNuclearEffect_HMProton_MeV;
  TLorentzVector DeltaPNuclearEffect_FirstProton_MeV;

//DeltaPProton
  TVector3 DeltaPProton_MeV;

//DeltaPTargetNucleon
  TVector3 DeltaPTargetNucleon_MeV;

//DeltaAlphat
  Double_t DeltaAlphaT_HMProton_deg;
  Double_t DeltaAlphaT_FirstProton_deg;

//DeltaP_TT
  Int_t DeltaP_TT_PionPDG;
  Double_t DeltaP_TT;

//ProtonPion Combo Platter
  TVector3 HMProtonPion_3Mom_MeV;
  Double_t DeltaPhiT_HMProtonPion_deg;
  TVector3 DeltaPT_HMProtonPion_MeV;
  Double_t DeltaAlphaT_HMProtonPion_deg;
  TLorentzVector DeltaPTotal_HMProtonPion_MeV;


//******************************************************************************
//                       Subsequent Species Sums
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

  Int_t NAboveThresholdProtons[kNThreshMax];
  Int_t NAboveThresholdGammas[kNThreshMax];
  Int_t NAboveThresholdPiPlus[kNThreshMax];
  Int_t NAboveThresholdPiMinus[kNThreshMax];
  Int_t NAboveThresholdChargedPions[kNThreshMax];
  Int_t NAboveThresholdTrackable[kNThreshMax];
  Int_t NAboveThresholdNeutrons[kNThreshMax];
  Int_t NAboveThresholdPiZero[kNThreshMax];
  Int_t NAboveThresholdNeutrals[kNThreshMax];
  Int_t NAboveThresholdExotic[kNThreshMax];

  Int_t NInEKinBinProtons[kNThreshMax];
  Int_t NInEKinBinGammas[kNThreshMax];
  Int_t NInEKinBinPiPlus[kNThreshMax];
  Int_t NInEKinBinPiMinus[kNThreshMax];
  Int_t NInEKinBinChargedPions[kNThreshMax];
  Int_t NInEKinBinTrackable[kNThreshMax];
  Int_t NInEKinBinNeutrons[kNThreshMax];
  Int_t NInEKinBinPiZero[kNThreshMax];
  Int_t NInEKinBinNeutrals[kNThreshMax];
  Int_t NInEKinBinExotic[kNThreshMax];

//******************************************************************************
//                      FS Particle Stuff
//******************************************************************************
  TLorentzVector __OtherFSPiPlus4Momenta_MeV[kMaxFSMomenta + 1]; //!
  TLorentzVector __OtherFSProton4Momenta_MeV[kMaxFSMomenta + 1]; //!

  Int_t NOtherFSPiPlus4Momenta_MeV;
  Int_t NOtherFSProton4Momenta_MeV;

  Double_t* OtherFSPiPlus4Momenta_MeV_X; //[NOtherFSPiPlus4Momenta_MeV]
  Double_t* OtherFSPiPlus4Momenta_MeV_Y; //[NOtherFSPiPlus4Momenta_MeV]
  Double_t* OtherFSPiPlus4Momenta_MeV_Z; //[NOtherFSPiPlus4Momenta_MeV]
  Double_t* OtherFSPiPlus4Momenta_MeV_T; //[NOtherFSPiPlus4Momenta_MeV]
  Double_t* OtherFSProton4Momenta_MeV_X; //[NOtherFSProton4Momenta_MeV]
  Double_t* OtherFSProton4Momenta_MeV_Y; //[NOtherFSProton4Momenta_MeV]
  Double_t* OtherFSProton4Momenta_MeV_Z; //[NOtherFSProton4Momenta_MeV]
  Double_t* OtherFSProton4Momenta_MeV_T; //[NOtherFSProton4Momenta_MeV]


//******************************************************************************
//                       Tangible Target Traits
//******************************************************************************

  Int_t TargetPDG;
  Int_t TargetZ;

//******************************************************************************
//                       Others and Transients
//******************************************************************************

  TLorentzVector CCQE_deltaP_MeV;

  Bool_t ProtonRescat_contains_NoInt;
  Bool_t ProtonRescat_contains_chrgEx;
  Bool_t ProtonRescat_contains_elastic;
  Bool_t ProtonRescat_contains_inelastic;
  Bool_t ProtonRescat_contains_knockout;

  Double_t CCQ2;

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
  void HandleRescat(Int_t PDG, Int_t RescatCode);

  void HandleStruckNucleon(TLorentzVector &StdHepPTLV,
    Double_t &StdHepP3Mod, Int_t pdg);
  void Finalise();
  void Reset();

  Int_t GetIncNeutrino_PDG(){
    return IncNeutrino_PDG;
  }
  Int_t GetStruckNucleonPDG(){
    return StruckNucleonPDG;
  }
  Int_t GetNeutConventionReactionCode(){
    return NeutConventionReactionCode;
  }

  ClassDef(TransversityVars,1);
};

struct TransversityVarsLite : TransversityVarsB {

  bool IsInGev; //!

public:

  TransversityVarsLite(bool InGeV=true);
  ~TransversityVarsLite();

//******************************************************************************
//                     Event Properties
//******************************************************************************

//Generator reaction code
  Int_t NeutConventionReactionCode;

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

//Highest Momentum Proton
  Int_t HMProton_PDG;
  TLorentzVector HMProton_4Mom_MeV;

//Highest Momentum Charged Pion
  Int_t HMCPion_PDG;
  TLorentzVector HMCPion_4Mom_MeV;

//******************************************************************************
//                       'Verse Variable Values
//******************************************************************************

//DeltaPhiT
  Double_t DeltaPhiT_HMProton_deg;

//DeltaPt
  TVector3 DeltaPT_HMProton_MeV;

//DeltaPTotal
  TLorentzVector DeltaPTotal_HMProton_MeV;

//DeltaPNuclearEffect
  TLorentzVector DeltaPNuclearEffect_HMProton_MeV;

//DeltaPProton
  TVector3 DeltaPProton_MeV;

//DeltaPTargetNucleon
  TVector3 DeltaPTargetNucleon_MeV;

//DeltaAlphat
  Double_t DeltaAlphaT_HMProton_deg;

//DeltaP_TT
  Double_t DeltaP_TT;

//ProtonPion Combo Platter
  TVector3 HMProtonPion_3Mom_MeV;
  Double_t DeltaPhiT_HMProtonPion_deg;
  TVector3 DeltaPT_HMProtonPion_MeV;
  Double_t DeltaAlphaT_HMProtonPion_deg;
  TLorentzVector DeltaPTotal_HMProtonPion_MeV;


//******************************************************************************
//                       Subsequent Species Sums
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

  Double_t CCQ2;

  //Transients
  PartStruct Muon; //!
  PartStruct MuonNeutrino; //!
  PartStruct StruckNucleon; //!
  PartStruct HMProton; //!
  PartStruct HMCPion; //!
  PartStruct FirstProton; //!
  Double_t TargetBE_MeV; //!


//******************************************************************************
//******************************************************************************

private:
  void HandleProton(TLorentzVector &StdHepPTLV,
    Double_t &StdHepP3Mod, UInt_t &StdHepPosition);
  void HandleCPion(TLorentzVector &StdHepPTLV,
    Double_t &StdHepP3Mod, Int_t pdg);
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
  void HandleRescat(Int_t PDG, Int_t RescatCode) {}
  void Finalise();
  void Reset();

  Int_t GetIncNeutrino_PDG(){
    return IncNeutrino_PDG;
  }
  Int_t GetStruckNucleonPDG(){
    return StruckNucleonPDG;
  }
  Int_t GetNeutConventionReactionCode(){
    return NeutConventionReactionCode;
  }

  ClassDef(TransversityVarsLite,1);
};

#endif
