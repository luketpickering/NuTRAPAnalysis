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
protected:
  bool IsInGev; //!
public:

  TransversityVarsB(bool InGeV=true);
  virtual ~TransversityVarsB();

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

//Highest Momentum Trackable
  Int_t HMTrackable_PDG;
  TLorentzVector HMTrackable_4Mom_MeV;

//******************************************************************************
//                       'Verse Variable Values
//******************************************************************************

//deltaphiT
  Double_t deltaphiT_HMProton_deg;

//Deltap
  TLorentzVector Deltap_HMProton_MeV;

//deltap
  TLorentzVector deltap_HMProton_MeV;

//deltapT
  TVector3 deltapT_HMProton_MeV;

//deltaalphaT
  Double_t deltaalphaT_HMProton_deg;

//deltap_TT
  Double_t deltap_TT;

//ProtonPion Combo Platter
  TVector3 HMProtonPion_3Mom_MeV;
  Double_t deltaphiT_HMProtonPion_deg;
  TVector3 deltapT_HMProtonPion_MeV;
  Double_t deltaalphaT_HMProtonPion_deg;
  TLorentzVector Deltap_HMProtonPion_MeV;

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

//******************************************************************************
//                       Tangible Target Traits
//******************************************************************************

  Int_t TargetPDG;
  Int_t TargetZ;

//******************************************************************************
//                       Others and Transients
//******************************************************************************

  Double_t CCQ2;

  //Transients
  PartStruct Muon; //!
  PartStruct MuonNeutrino; //!
  PartStruct StruckNucleon; //!
  PartStruct HMProton; //!
  PartStruct HMCPion; //!
  PartStruct HMTrackable; //!

//******************************************************************************
//******************************************************************************

protected:
  virtual void HandleProton(TLorentzVector &StdHepPTLV,
    Double_t &StdHepP3Mod, UInt_t &StdHepPosition);
  virtual void HandleCPion(TLorentzVector &StdHepPTLV,
    Double_t &StdHepP3Mod, Int_t pdg);
  virtual void HandleHMTrackable(TLorentzVector &StdHepPTLV,
    Double_t &StdHepP3Mod, Int_t PDG);
public:

  virtual bool HandleStdHepParticle(UInt_t &StdHepPosition,
                            Int_t &StdHepPdg,
                            Int_t &StdHepStatus,
                            Double_t * &StdHepP4);
  virtual void HandleStruckNucleon(TLorentzVector &StdHepPTLV,
    Double_t &StdHepP3Mod, Int_t pdg);
  virtual void Finalise();
  virtual void Reset();

  ClassDef(TransversityVarsB,1);

};

struct TransversityVars : public TransversityVarsB {

public:

  TransversityVars(bool InGeV=true, Int_t NThresh=0, Int_t* Threshs_MeV=0,
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

  Double_t ReconNuEnergy;
  Double_t ReconTargetMass;

//******************************************************************************
//                     Pertinent Particle Properties
//******************************************************************************

//Muon
  TVector3 Muon_Pt_MeV;

//First Proton
  Int_t FirstProton_PDG;
  TLorentzVector FirstProton_4Mom_MeV;
  Int_t FirstProton_StdHepPosition;

//Highest Momentum Proton
  Int_t HMProton_StdHepPosition;

//StruckNucleon_3Mom_Recon
  TVector3 StruckNucleon_3Mom_Recon_MeV;

//PreFSINucleon_3Mom_Recon
  TVector3 PreFSINucleon_3Mom_Recon_MeV;

//******************************************************************************
//                       'Verse Variable Values
//******************************************************************************

//deltaphiT
  Double_t deltaphiT_FirstProton_deg;
  Double_t deltaphiT_HMTrackable_deg;

//deltapT
  TVector3 deltapT_FirstProton_MeV;

//Deltap
  TLorentzVector Deltap_FirstProton_MeV;

//deltap
  TLorentzVector deltap_FirstProton_MeV;

//deltaalphaT
  Double_t deltaalphaT_FirstProton_deg;

//deltap_TT
  Int_t deltap_TT_PionPDG;

//******************************************************************************
//                       Subsequent Species Sums
//******************************************************************************

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
//                       Others and Transients
//******************************************************************************

  Bool_t ProtonRescat_contains_NoInt;
  Bool_t ProtonRescat_contains_chrgEx;
  Bool_t ProtonRescat_contains_elastic;
  Bool_t ProtonRescat_contains_inelastic;
  Bool_t ProtonRescat_contains_knockout;

  PartStruct FirstProton; //!

//******************************************************************************
//******************************************************************************

private:
  void HandleProton(TLorentzVector &StdHepPTLV,
    Double_t &StdHepP3Mod, UInt_t &StdHepPosition);
  void HandleCPion(TLorentzVector &StdHepPTLV,
    Double_t &StdHepP3Mod, Int_t pdg);
  Int_t PartIsAboveThresh(TLorentzVector FourMom);
public:

  bool HandleStdHepParticle(UInt_t &StdHepPosition,
                            Int_t &StdHepPdg,
                            Int_t &StdHepStatus,
                            Double_t * &StdHepP4);
  void HandleRescat(Int_t PDG, Int_t RescatCode);

  void Finalise();
  void Reset();

  ClassDef(TransversityVars,1);
};

#endif
