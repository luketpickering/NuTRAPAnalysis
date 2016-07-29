#ifndef __TRANSVERSITY_VARIABLE_OBJECTS_SEEN__
#define __TRANSVERSITY_VARIABLE_OBJECTS_SEEN__

#include <stdexcept>

#include "TObject.h"
#include "TLorentzVector.h"
#include "TTree.h"

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

struct TransversityVarsB {
protected:
  bool IsInGev; //!

  //Transients
  PartStruct Muon; //!
  PartStruct MuonNeutrino; //!
  PartStruct StruckNucleon; //!
  PartStruct HMProton; //!
  PartStruct HMCPion; //!
  PartStruct HMPion; //!
  PartStruct HMTrackable; //!

  TLorentzVector* _IncNeutrino_4Mom_MeV;
  TLorentzVector* _StruckNucleon_4Mom_MeV;
  TLorentzVector* _FourMomentumTransfer;
  TLorentzVector* _Muon_4Mom_MeV;
  TLorentzVector* _HMProton_4Mom_MeV;
  TLorentzVector* _HMCPion_4Mom_MeV;
  TLorentzVector* _HMPion_4Mom_MeV;
  TLorentzVector* _HMTrackable_4Mom_MeV;
  TLorentzVector* _Deltap_HMProton_MeV;
  TLorentzVector* _deltap_HMProton_MeV;
  TVector3* _deltapT_HMProton_MeV;
  TVector3* _HMProtonPion_3Mom_MeV;
  TVector3* _deltapT_HMProtonPion_MeV;
  TLorentzVector* _Deltap_HMProtonPion_MeV;

  TVector3* _FS_PSum;
  TVector3* _ChargedFS_PSum;


public:

  TransversityVarsB(bool InGeV=true);
  virtual ~TransversityVarsB();

//******************************************************************************
//                     Event Properties
//******************************************************************************

//Generator reaction code
  Int_t NeutConventionReactionCode;
  Double_t EvtWght;

//******************************************************************************
//                     Pertinent Particle Properties
//******************************************************************************

//Neutrino
  Int_t IncNeutrino_PDG;
  TLorentzVector IncNeutrino_4Mom_MeV;

//Struck Nucleon
  Int_t StruckNucleonPDG;
  TLorentzVector StruckNucleon_4Mom_MeV;

  TLorentzVector FourMomentumTransfer;

//Muon
  Int_t Muon_PDG;
  TLorentzVector Muon_4Mom_MeV;

//Highest Momentum Proton
  Int_t HMProton_PDG;
  TLorentzVector HMProton_4Mom_MeV;

//Highest Momentum Pion
  Int_t HMPion_PDG;
  TLorentzVector HMPion_4Mom_MeV;

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

  Double_t HadrMass;
  Double_t ChargedHadrMass;

  Double_t FS_ESum;
  TVector3 FS_PSum;
  Double_t ChargedFS_ESum;
  TVector3 ChargedFS_PSum;

//******************************************************************************
//******************************************************************************

protected:
  virtual void HandleProton(TLorentzVector &StdHepPTLV,
    Double_t &StdHepP3Mod, UInt_t &StdHepPosition);
  virtual void HandlePion(TLorentzVector &StdHepPTLV,
    Double_t &StdHepP3Mod, Int_t pdg);
  virtual void HandleCPion(TLorentzVector &StdHepPTLV,
    Double_t &StdHepP3Mod, Int_t pdg);
  virtual void HandleHMTrackable(TLorentzVector &StdHepPTLV,
    Double_t &StdHepP3Mod, Int_t PDG);
public:
  virtual void AddBranches(TTree* tree);
  virtual bool HandleStdHepParticle(UInt_t &StdHepPosition,
                            Int_t &StdHepPdg,
                            Int_t &StdHepStatus,
                            Double_t * &StdHepP4);
  virtual void HandleStruckNucleon(TLorentzVector &StdHepPTLV, Int_t pdg);
  virtual void Finalise();
  virtual void Reset();

};

struct TransversityVars : public TransversityVarsB {

protected:
  PartStruct FirstProton; //!
  TString* _GeneratorName;
  TVector3* _Muon_Pt_MeV;
  TLorentzVector* _FirstProton_4Mom_MeV;
  TVector3* _StruckNucleon_3Mom_Recon_MeV;
  TVector3* _PreFSINucleon_3Mom_Recon_MeV;
  TVector3* _deltapT_FirstProton_MeV;
  TLorentzVector* _Deltap_FirstProton_MeV;
  TLorentzVector* _deltap_FirstProton_MeV;
  TLorentzVector __OtherFSPiPlus4Momenta_MeV[kMaxFSMomenta + 1]; //!
  TLorentzVector __OtherFSProton4Momenta_MeV[kMaxFSMomenta + 1]; //!

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


  Int_t NOtherFSPiPlus4Momenta_MeV;
  Int_t NOtherFSProton4Momenta_MeV;

  Double_t OtherFSPiPlus4Momenta_MeV_X[kMaxFSMomenta]; //[NOtherFSPiPlus4Momenta_MeV]
  Double_t OtherFSPiPlus4Momenta_MeV_Y[kMaxFSMomenta]; //[NOtherFSPiPlus4Momenta_MeV]
  Double_t OtherFSPiPlus4Momenta_MeV_Z[kMaxFSMomenta]; //[NOtherFSPiPlus4Momenta_MeV]
  Double_t OtherFSPiPlus4Momenta_MeV_T[kMaxFSMomenta]; //[NOtherFSPiPlus4Momenta_MeV]
  Double_t OtherFSProton4Momenta_MeV_X[kMaxFSMomenta]; //[NOtherFSProton4Momenta_MeV]
  Double_t OtherFSProton4Momenta_MeV_Y[kMaxFSMomenta]; //[NOtherFSProton4Momenta_MeV]
  Double_t OtherFSProton4Momenta_MeV_Z[kMaxFSMomenta]; //[NOtherFSProton4Momenta_MeV]
  Double_t OtherFSProton4Momenta_MeV_T[kMaxFSMomenta]; //[NOtherFSProton4Momenta_MeV]

//******************************************************************************
//                       Others and Transients
//******************************************************************************

//******************************************************************************
//******************************************************************************

private:
  void HandleProton(TLorentzVector &StdHepPTLV,
    Double_t &StdHepP3Mod, UInt_t &StdHepPosition);
  void HandleCPion(TLorentzVector &StdHepPTLV,
    Double_t &StdHepP3Mod, Int_t pdg);
  Int_t PartIsAboveThresh(TLorentzVector FourMom);
public:
  void AddBranches(TTree* tree);
  bool HandleStdHepParticle(UInt_t &StdHepPosition,
                            Int_t &StdHepPdg,
                            Int_t &StdHepStatus,
                            Double_t * &StdHepP4);

  void Finalise();
  void Reset();
};

TransversityVarsB const * const MakeReadingTransversityVarsB(TTree* tree);
void UnsetBranchAddressesTransversityVarsB(TTree* tree,
  TransversityVarsB const * const tvb);
TransversityVars* MakeReadingTransversityVars(TTree* tree);

#endif
