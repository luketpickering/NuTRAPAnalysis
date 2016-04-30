#include <iostream>

#include "TTree.h"

#include "TransversityVariableObjects.hxx"

using namespace TransversityUtils;

namespace {
constexpr int kStdHepIdxPx = 0;
constexpr int kStdHepIdxPy = 1;
constexpr int kStdHepIdxPz = 2;
constexpr int kStdHepIdxE = 3;

template<size_t N>
void ClearArray(TLorentzVector (&arr)[N]){
  for(size_t i = 0; i < N; ++i){
    arr[i].SetXYZT(0,0,0,0);
  }
}
template<size_t N>
void ClearArray(Double_t (&arr)[N]){
  for(size_t i = 0; i < N; ++i){
    arr[i] = 0.0;
  }
}

void ClearArray(Double_t *arr, size_t N){
  for(size_t i = 0; i < N; ++i){
    arr[i] = 0.0;
  }
}

}

TransversityVarsB::~TransversityVarsB(){}

TransversityVarsB::TransversityVarsB(
  bool InGeV){

  IsInGev = InGeV;

  _IncNeutrino_4Mom_MeV = &IncNeutrino_4Mom_MeV;
  _StruckNucleon_4Mom_MeV = &StruckNucleon_4Mom_MeV;
  _Muon_4Mom_MeV = &Muon_4Mom_MeV;
  _HMProton_4Mom_MeV = &HMProton_4Mom_MeV;
  _HMCPion_4Mom_MeV = &HMCPion_4Mom_MeV;
  _HMPion_4Mom_MeV = &HMPion_4Mom_MeV;
  _HMTrackable_4Mom_MeV = &HMTrackable_4Mom_MeV;
  _Deltap_HMProton_MeV = &Deltap_HMProton_MeV;
  _deltap_HMProton_MeV = &deltap_HMProton_MeV;
  _deltapT_HMProton_MeV = &deltapT_HMProton_MeV;
  _HMProtonPion_3Mom_MeV = &HMProtonPion_3Mom_MeV;
  _deltapT_HMProtonPion_MeV = &deltapT_HMProtonPion_MeV;
  _Deltap_HMProtonPion_MeV = &Deltap_HMProtonPion_MeV;

  Reset();
}

void TransversityVarsB::HandleProton(TLorentzVector &StdHepPTLV,
  Double_t &StdHepP3Mod, UInt_t &StdHepPosition){

  if(!HMProton.PDG){
    HMProton.PDG = 2212;
    HMProton.Momentum = StdHepP3Mod;
    HMProton.FourMomentum = StdHepPTLV;
    HMProton.StdHepPosition = StdHepPosition;
  }
  else if(StdHepP3Mod > HMProton.Momentum){
    HMProton.StdHepPosition = StdHepPosition;
    HMProton.Momentum = StdHepP3Mod;
    HMProton.FourMomentum = StdHepPTLV;
  }
}

void TransversityVarsB::HandleHMTrackable(TLorentzVector &StdHepPTLV,
  Double_t &StdHepP3Mod, Int_t PDG){

  if(!HMTrackable.PDG){
    HMTrackable.PDG = PDG;
    HMTrackable.Momentum = StdHepP3Mod;
    HMTrackable.FourMomentum = StdHepPTLV;
  }
  else if(StdHepP3Mod > HMTrackable.Momentum){
    HMTrackable.PDG = PDG;
    HMTrackable.Momentum = StdHepP3Mod;
    HMTrackable.FourMomentum = StdHepPTLV;
  }
}

void TransversityVarsB::HandlePion(TLorentzVector &StdHepPTLV,
  Double_t &StdHepP3Mod, Int_t pdg){
  if(StdHepP3Mod > HMPion.Momentum){
    HMPion.Momentum = StdHepP3Mod;
    HMPion.FourMomentum = StdHepPTLV;
    HMPion.PDG = pdg;
  }
}


void TransversityVarsB::HandleCPion(TLorentzVector &StdHepPTLV,
  Double_t &StdHepP3Mod, Int_t pdg){

  if(StdHepP3Mod > HMCPion.Momentum){
    HMCPion.Momentum = StdHepP3Mod;
    HMCPion.FourMomentum = StdHepPTLV;
    HMCPion.PDG = pdg;
  }
  HandlePion(StdHepPTLV,StdHepP3Mod,pdg);
}

void TransversityVarsB::HandleStruckNucleon(TLorentzVector &StdHepPTLV,
  Int_t pdg){

  StruckNucleon.Momentum = StdHepPTLV.Vect().Mag();
  StruckNucleon.FourMomentum = StdHepPTLV;
  StruckNucleon.PDG = pdg;
}

bool TransversityVarsB::HandleStdHepParticle(
  UInt_t &StdHepPosition,
  Int_t &StdHepPdg,
  Int_t &StdHepStatus,
  Double_t * &StdHepP4){

  TLorentzVector StdHepPTLV = TLorentzVector(
    StdHepP4[kStdHepIdxPx],
    StdHepP4[kStdHepIdxPy],
    StdHepP4[kStdHepIdxPz],
    StdHepP4[kStdHepIdxE]);
  Double_t StdHepP3Mod = StdHepPTLV.Vect().Mag();

  if(StdHepPosition == 0){ //Incomming Neutrino
    MuonNeutrino.PDG = StdHepPdg;
    MuonNeutrino.FourMomentum = StdHepPTLV;
    MuonNeutrino.Momentum = StdHepP3Mod;
    return true;
  } else if(StdHepPosition == 1){ //Target Nucleus
    TargetPDG = StdHepPdg;
    TargetZ = ((StdHepPdg/10000)%1000);
    return true;
  }

  if(StdHepStatus != 1){return false;}
  if(StdHepPdg >= 1000000000){return false;} //Should catch nuclear PDGs
  if(StdHepPdg >= 2000000000){return false;} //GENIE psuedo particles

  if(StdHepPdg > 100){
    FS_ESum += StdHepPTLV.E();
    FS_PSum += StdHepPTLV.Vect();
  }

  switch(StdHepPdg){
    case -13:
    case 13:{
      Muon.PDG = StdHepPdg;
      Muon.FourMomentum = StdHepPTLV;
      Muon.Momentum = StdHepP3Mod;
      NFinalStateParticles++;
      break;
    }
    //NC
    case -14:
    case 14:{
      return false;
      break;
    }
    case 22:{
      NGammas++;
      NFinalStateParticles++;
      break;
    }
    case 2212:{
      this->HandleProton(StdHepPTLV, StdHepP3Mod, StdHepPosition);
      this->HandleHMTrackable(StdHepPTLV,StdHepP3Mod,StdHepPdg);
      NProtons++;
      NFinalStateParticles++;
      ChargedFS_ESum += StdHepPTLV.E();
      ChargedFS_PSum += StdHepPTLV.Vect();
      break;
    }
    case 2112:{
      NNeutrons++;
      NFinalStateParticles++;
      break;
    }
    case 211:{
      this->HandleCPion(StdHepPTLV,StdHepP3Mod,StdHepPdg);
      this->HandleHMTrackable(StdHepPTLV,StdHepP3Mod,StdHepPdg);
      NPiPlus++;
      NChargedPions++;
      NPions++;
      NFinalStateParticles++;
      ChargedFS_ESum += StdHepPTLV.E();
      ChargedFS_PSum += StdHepPTLV.Vect();
      break;
    }
    case -211:{
      this->HandleCPion(StdHepPTLV,StdHepP3Mod,StdHepPdg);
      this->HandleHMTrackable(StdHepPTLV,StdHepP3Mod,StdHepPdg);
      NPiMinus++;
      NChargedPions++;
      NPions++;
      NFinalStateParticles++;
      ChargedFS_ESum += StdHepPTLV.E();
      ChargedFS_PSum += StdHepPTLV.Vect();
      break;
    }
    case 111:{
      this->HandlePion(StdHepPTLV,StdHepP3Mod,StdHepPdg);
      NPiZero++;
      NPions++;
      NFinalStateParticles++;
      break;
    }
    case 221:
    case 331:
    case 130:
    case 310:
    case 311:
    case 321:
    case -221:
    case -331:
    case -130:
    case -310:
    case -311:
    case -321:
    case 3122: {
      NFinalStateParticles++;
      break;
    }
    default:{
      NOtherParticles++;
      NFinalStateParticles++;
      std::cout << "[INFO]: Found Particle with PDG: " << StdHepPdg
        << std::endl;
      return false;
      break;
    }
  }
  return true;
}

void TransversityVarsB::Finalise(){

  constexpr static Float_t GeVToMeV = 1000.0;
  constexpr static Float_t RadToDeg = 180.0/M_PI;
  constexpr static Float_t DegToRad = M_PI/180.0;

  if(IsInGev){
    Muon.Momentum *= GeVToMeV;
    Muon.FourMomentum *= GeVToMeV;
    MuonNeutrino.Momentum *= GeVToMeV;
    MuonNeutrino.FourMomentum *= GeVToMeV;
    HMProton.Momentum *= GeVToMeV;
    HMProton.FourMomentum *= GeVToMeV;
    HMTrackable.Momentum *= GeVToMeV;
    HMTrackable.FourMomentum *= GeVToMeV;
    HMCPion.Momentum *= GeVToMeV;
    HMCPion.FourMomentum *= GeVToMeV;
    StruckNucleon.Momentum *= GeVToMeV;
    StruckNucleon.FourMomentum *= GeVToMeV;
  }

//******************************************************************************
//                     Pertinent Particle Properties
//******************************************************************************

//Neutrino
  IncNeutrino_PDG = MuonNeutrino.PDG;
  IncNeutrino_4Mom_MeV = MuonNeutrino.FourMomentum;

//Struck Nucleon
  StruckNucleonPDG = StruckNucleon.PDG;
  StruckNucleon_4Mom_MeV = StruckNucleon.FourMomentum;

//Muon
  Muon_PDG = Muon.PDG;
  Muon_4Mom_MeV = Muon.FourMomentum;
  TVector3 Muon_Pt_MeV = TransversityUtils::GetVectorInTPlane(
    Muon.FourMomentum.Vect(),
    MuonNeutrino.FourMomentum.Vect());

//Highest Momentum Proton
  HMProton_PDG = HMProton.PDG;
  HMProton_4Mom_MeV = HMProton.FourMomentum;

//Highest Momentum Pion
  HMPion_PDG = HMPion.PDG;
  HMPion_4Mom_MeV = HMPion.FourMomentum;

//Highest Momentum Charged Pion
  HMCPion_PDG = HMCPion.PDG;
  HMCPion_4Mom_MeV = HMCPion.FourMomentum;

//Highest Momentum Trackable
  HMTrackable_PDG = HMTrackable.PDG;
  HMTrackable_4Mom_MeV = HMTrackable.FourMomentum;

//******************************************************************************
//                       'Verse Variable Values
//******************************************************************************

//Projections
  TVector3 HMProtonPt_MeV = GetVectorInTPlane(HMProton.FourMomentum.Vect(),
    MuonNeutrino.FourMomentum.Vect());


//deltaphiT
  deltaphiT_HMProton_deg = GetDeltaPhiT(
    Muon.FourMomentum.Vect(), HMProton.FourMomentum.Vect(),
    MuonNeutrino.FourMomentum.Vect())*RadToDeg;

//Deltap
  Deltap_HMProton_MeV = (Muon.FourMomentum + HMProton.FourMomentum
    - MuonNeutrino.FourMomentum - StruckNucleon.FourMomentum);

//deltap
  deltap_HMProton_MeV = (Muon.FourMomentum + HMProton.FourMomentum
    - MuonNeutrino.FourMomentum);

//deltapt
  deltapT_HMProton_MeV = GetDeltaPT(Muon_Pt_MeV,
                                    HMProtonPt_MeV,
                                    MuonNeutrino.FourMomentum.Vect());

//deltaalphaT
  deltaalphaT_HMProton_deg = GetDeltaAlphaT(Muon_Pt_MeV,
                                    HMProtonPt_MeV,
                                    MuonNeutrino.FourMomentum.Vect())*RadToDeg;
//deltap_TT
  if((HMPion.Momentum > 1E-3) && (HMProton.Momentum > 1E-3)){

    deltap_TT = TransversityUtils::GetDeltaPTT(
      Muon.FourMomentum.Vect(),
      HMPion.FourMomentum.Vect(),
      HMProton.FourMomentum.Vect(),
      MuonNeutrino.FourMomentum.Vect());

    HMProtonPion_3Mom_MeV = HMProton.FourMomentum.Vect() +
      HMPion.FourMomentum.Vect();

    TVector3 HMProtonPionPt_MeV = GetVectorInTPlane(HMProtonPion_3Mom_MeV,
      MuonNeutrino.FourMomentum.Vect());

    deltaphiT_HMProtonPion_deg = GetDeltaPhiT(
                    Muon.FourMomentum.Vect(),
                    HMProtonPion_3Mom_MeV,
                    MuonNeutrino.FourMomentum.Vect())*RadToDeg;
    deltapT_HMProtonPion_MeV = GetDeltaPT(Muon_Pt_MeV,
                    HMProtonPionPt_MeV,
                    MuonNeutrino.FourMomentum.Vect());
    deltaalphaT_HMProtonPion_deg = GetDeltaAlphaT(Muon_Pt_MeV,
                    HMProtonPionPt_MeV,
                    MuonNeutrino.FourMomentum.Vect())*RadToDeg;

    Deltap_HMProtonPion_MeV = (Muon.FourMomentum + HMProton.FourMomentum
      + HMPion.FourMomentum - MuonNeutrino.FourMomentum
      - StruckNucleon.FourMomentum);
  }

  CCQ2 = (Muon.FourMomentum - MuonNeutrino.FourMomentum).Mag2();

  HadrMass = sqrt(FS_ESum*FS_ESum - FS_PSum.Mag2());
  ChargedHadrMass = sqrt(ChargedFS_ESum*ChargedFS_ESum - ChargedFS_PSum.Mag2());

//******************************************************************************
//******************************************************************************
}

void TransversityVarsB::Reset(){
//******************************************************************************
//                     Event Properties
//******************************************************************************

//Generator reaction code
  NeutConventionReactionCode = 0;

//******************************************************************************
//                     Pertinent Particle Properties
//******************************************************************************

//Neutrino
  IncNeutrino_PDG = 0;
  IncNeutrino_4Mom_MeV = TLorentzVector(0,0,0,0);

//Struck Nucleon
  StruckNucleonPDG = 0;
  StruckNucleon_4Mom_MeV = TLorentzVector(0,0,0,0);

//Muon
  Muon_PDG = 0;
  Muon_4Mom_MeV = TLorentzVector(0,0,0,0);

//Highest Momentum Proton
  HMProton_PDG = 0;
  HMProton_4Mom_MeV = TLorentzVector(0,0,0,0);

//Highest Momentum Pion
  HMPion_PDG = 0;
  HMPion_4Mom_MeV = TLorentzVector(0,0,0,0);

//Highest Momentum Charged Pion
  HMCPion_PDG = 0;
  HMCPion_4Mom_MeV = TLorentzVector(0,0,0,0);

//Highest Momentum Trackable
  HMTrackable_PDG = 0;
  HMTrackable_4Mom_MeV = TLorentzVector(0,0,0,0);

//******************************************************************************
//                       'Verse Variable Values
//******************************************************************************

//deltaphiT
  deltaphiT_HMProton_deg = 0;

//deltapT
  deltapT_HMProton_MeV = TVector3(0,0,0);

//Deltap
  Deltap_HMProton_MeV = TLorentzVector(0,0,0,0);

//deltap
  deltap_HMProton_MeV = TLorentzVector(0,0,0,0);

//deltaalphaT
  deltaalphaT_HMProton_deg = 0;

//deltap_TT
  deltap_TT = 0;

//ProtonPion Combo Platter
  HMProtonPion_3Mom_MeV = TVector3(0,0,0);
  deltaphiT_HMProtonPion_deg = 0;
  deltapT_HMProtonPion_MeV = TVector3(0,0,0);
  deltaalphaT_HMProtonPion_deg = 0;
  Deltap_HMProtonPion_MeV = TLorentzVector(0,0,0,0);

//******************************************************************************
//                       Subsequent Species Sums
//******************************************************************************

  NFinalStateParticles = 0;

  NProtons = 0;
  NGammas = 0;
  NNeutrons = 0;
  NPiPlus = 0;
  NPiZero = 0;
  NPiMinus = 0;
  NPions = 0;
  NChargedPions = 0;
  NOtherParticles = 0;

//******************************************************************************
//                       Tangible Target Traits
//******************************************************************************

  TargetPDG = 0;
  TargetZ = 0;
  HadrMass = 0;
  ChargedHadrMass = 0;

//******************************************************************************
//                       Others and Transients
//******************************************************************************

  CCQ2 = 0;

  //Transients
  Muon.Reset();
  MuonNeutrino.Reset();
  StruckNucleon.Reset();
  HMProton.Reset();
  HMCPion.Reset();
  HMPion.Reset();
  HMTrackable.Reset();

  FS_ESum = 0;
  FS_PSum = TVector3(0,0,0);
  ChargedFS_ESum = 0;
  ChargedFS_PSum = TVector3(0,0,0);
//******************************************************************************
//******************************************************************************

}

void TransversityVarsB::AddBranches(TTree* tree){

//******************************************************************************
//                     Event Properties
//******************************************************************************

//Generator reaction code
  tree->Branch("NeutConventionReactionCode",&NeutConventionReactionCode,
    "NeutConventionReactionCode/I");

//******************************************************************************
//                     Pertinent Particle Properties
//******************************************************************************

//Neutrino
  tree->Branch("IncNeutrino_PDG",&IncNeutrino_PDG, "IncNeutrino_PDG/I");
  tree->Branch("IncNeutrino_4Mom_MeV", &_IncNeutrino_4Mom_MeV);

//Struck Nucleon
  tree->Branch("StruckNucleonPDG",&StruckNucleonPDG, "StruckNucleonPDG/I");
  tree->Branch("StruckNucleon_4Mom_MeV", &_StruckNucleon_4Mom_MeV);

//Muon
  tree->Branch("Muon_PDG",&Muon_PDG, "Muon_PDG/I");
  tree->Branch("Muon_4Mom_MeV", &_Muon_4Mom_MeV);

//Highest Momentum Proton
  tree->Branch("HMProton_PDG",&HMProton_PDG, "HMProton_PDG/I");
  tree->Branch("HMProton_4Mom_MeV", &_HMProton_4Mom_MeV);

//Highest Momentum Charged Pion
  tree->Branch("HMCPion_PDG",&HMCPion_PDG, "HMCPion_PDG/I");
  tree->Branch("HMCPion_4Mom_MeV", &_HMCPion_4Mom_MeV);

//Highest Momentum Pion
  tree->Branch("HMPion_PDG",&HMPion_PDG, "HMPion_PDG/I");
  tree->Branch("HMPion_4Mom_MeV", &_HMPion_4Mom_MeV);

//Highest Momentum Trackable
  tree->Branch("HMTrackable_PDG",&HMTrackable_PDG, "HMTrackable_PDG/I");
  tree->Branch("HMTrackable_4Mom_MeV", &_HMTrackable_4Mom_MeV);

//******************************************************************************
//                       'Verse Variable Values
//******************************************************************************

//deltaphiT
  tree->Branch("deltaphiT_HMProton_deg", &deltaphiT_HMProton_deg,
    "deltaphiT_HMProton_deg/D");

//deltapT
  tree->Branch("deltapT_HMProton_MeV", &_deltapT_HMProton_MeV);

//Deltap
  tree->Branch("Deltap_HMProton_MeV", &_Deltap_HMProton_MeV);

//deltap
  tree->Branch("deltap_HMProton_MeV", &_deltap_HMProton_MeV);

//deltaalphaT
  tree->Branch("deltaalphaT_HMProton_deg", &deltaalphaT_HMProton_deg,
    "deltaalphaT_HMProton_deg/D");

//deltap_TT
  tree->Branch("deltap_TT", &deltap_TT, "deltap_TT/D");

//ProtonPion Combo Platter
  tree->Branch("HMProtonPion_3Mom_MeV", &_HMProtonPion_3Mom_MeV);
  tree->Branch("deltaphiT_HMProtonPion_deg", &deltaphiT_HMProtonPion_deg,
    "deltaphiT_HMProtonPion_deg/D");
  tree->Branch("deltapT_HMProtonPion_MeV", &_deltapT_HMProtonPion_MeV);
  tree->Branch("deltaalphaT_HMProtonPion_deg", &deltaalphaT_HMProtonPion_deg,
    "deltaalphaT_HMProtonPion_deg/D");
  tree->Branch("Deltap_HMProtonPion_MeV", &_Deltap_HMProtonPion_MeV);

//******************************************************************************
//                       Subsequent Species Sums
//******************************************************************************

  tree->Branch("NFinalStateParticles", &NFinalStateParticles,
    "NFinalStateParticles/I");

  tree->Branch("NProtons", &NProtons, "NProtons/I");
  tree->Branch("NGammas", &NGammas, "NGammas/I");
  tree->Branch("NNeutrons", &NNeutrons, "NNeutrons/I");
  tree->Branch("NPiPlus", &NPiPlus, "NPiPlus/I");
  tree->Branch("NPiZero", &NPiZero, "NPiZero/I");
  tree->Branch("NPiMinus", &NPiMinus, "NPiMinus/I");
  tree->Branch("NPions", &NPions, "NPions/I");
  tree->Branch("NChargedPions", &NChargedPions, "NChargedPions/I");
  tree->Branch("NOtherParticles", &NOtherParticles, "NOtherParticles/I");

//******************************************************************************
//                       Tangible Target Traits
//******************************************************************************

  tree->Branch("TargetPDG", &TargetPDG, "TargetPDG/I");
  tree->Branch("TargetZ", &TargetZ, "TargetZ/I");

//******************************************************************************
//                       Others and Transients
//******************************************************************************

  tree->Branch("CCQ2", &CCQ2, "CCQ2/D");
  tree->Branch("HadrMass", &HadrMass, "HadrMass/D");
  tree->Branch("ChargedHadrMass", &ChargedHadrMass, "ChargedHadrMass/D");

//******************************************************************************
//******************************************************************************

}

TransversityVars::~TransversityVars(){
}

TransversityVars::TransversityVars(
  bool InGeV, Int_t NThresh, Int_t* Threshs_MeV,
  TString generatorName) : TransversityVarsB(InGeV){
  GeneratorName = generatorName;

  this->NThresh = NThresh;
  this->Threshs_MeV = Threshs_MeV;

  _GeneratorName = &GeneratorName;
  _Muon_Pt_MeV = &Muon_Pt_MeV;
  _FirstProton_4Mom_MeV = &FirstProton_4Mom_MeV;
  _StruckNucleon_3Mom_Recon_MeV = &StruckNucleon_3Mom_Recon_MeV;
  _PreFSINucleon_3Mom_Recon_MeV = &PreFSINucleon_3Mom_Recon_MeV;
  _deltapT_FirstProton_MeV = &deltapT_FirstProton_MeV;
  _Deltap_FirstProton_MeV = &Deltap_FirstProton_MeV;
  _deltap_FirstProton_MeV = &deltap_FirstProton_MeV;

  Reset();
}

Int_t TransversityVars::PartIsAboveThresh(TLorentzVector FourMom){
  Double_t EKin = FourMom.E() - FourMom.M();
  for(int i = 0; i < NThresh; ++i){
    if(EKin < Threshs_MeV[i]){
      return (i-1);
    }
  }
  return (NThresh-1);
}

template <size_t N>
void IncrementThreshArray(Int_t (&NAboveThreshold) [N], Int_t threshlvl){
  for(size_t i = 0; ((i < N) && (Int_t(i) <= threshlvl)); ++i){
    NAboveThreshold[i]++;
  }
}

void TransversityVars::HandleProton(TLorentzVector &StdHepPTLV,
  Double_t &StdHepP3Mod, UInt_t &StdHepPosition){

  if(!FirstProton.PDG){
    FirstProton.PDG = 2212;
    FirstProton.Momentum = StdHepP3Mod;
    FirstProton.FourMomentum = StdHepPTLV;
    FirstProton.StdHepPosition = StdHepPosition;

    HMProton.PDG = 2212;
    HMProton.Momentum = StdHepP3Mod;
    HMProton.FourMomentum = StdHepPTLV;
    HMProton.StdHepPosition = StdHepPosition;

  }
  else if(StdHepP3Mod > HMProton.Momentum){
    HMProton.StdHepPosition = StdHepPosition;
    HMProton.Momentum = StdHepP3Mod;
    HMProton.FourMomentum = StdHepPTLV;
  }
  if( (StdHepPTLV.E()*(IsInGev?1000.0:1.0)) >
      __OtherFSProton4Momenta_MeV[kMaxFSMomenta].E() ){

    __OtherFSProton4Momenta_MeV[kMaxFSMomenta] =
      IsInGev?(StdHepPTLV*1000.0):(StdHepPTLV);
    std::stable_sort(
      __OtherFSProton4Momenta_MeV,
      __OtherFSProton4Momenta_MeV+kMaxFSMomenta+1,
      [](TLorentzVector const & LHS, TLorentzVector const & RHS) -> bool {
        return (LHS.E() > RHS.E());
      }
    );
  }
}

void TransversityVars::HandleCPion(TLorentzVector &StdHepPTLV,
  Double_t &StdHepP3Mod, Int_t pdg){

  if(StdHepP3Mod > HMCPion.Momentum){
    HMCPion.Momentum = StdHepP3Mod;
    HMCPion.FourMomentum = StdHepPTLV;
    HMCPion.PDG = pdg;
  }
  if((pdg == 211) &&
     ((StdHepPTLV.E()*(IsInGev?1000.0:1.0)) >
           __OtherFSPiPlus4Momenta_MeV[kMaxFSMomenta].E()) ){

    __OtherFSPiPlus4Momenta_MeV[kMaxFSMomenta] =
      IsInGev?(StdHepPTLV*1000.0):(StdHepPTLV);
    std::stable_sort(
      __OtherFSPiPlus4Momenta_MeV,
      __OtherFSPiPlus4Momenta_MeV+kMaxFSMomenta+1,
      [](TLorentzVector const & LHS, TLorentzVector const & RHS) -> bool {
        return (LHS.E() > RHS.E());
      });
  }
}

bool TransversityVars::HandleStdHepParticle(
  UInt_t &StdHepPosition,
  Int_t &StdHepPdg,
  Int_t &StdHepStatus,
  Double_t * &StdHepP4){

  if(!TransversityVarsB::HandleStdHepParticle(StdHepPosition,
                                              StdHepPdg,
                                              StdHepStatus,
                                              StdHepP4) ){
    return false;
  }

  TLorentzVector StdHepPTLV = TLorentzVector(
    StdHepP4[kStdHepIdxPx],
    StdHepP4[kStdHepIdxPy],
    StdHepP4[kStdHepIdxPz],
    StdHepP4[kStdHepIdxE]);
  Double_t StdHepP3Mod = StdHepPTLV.Vect().Mag();

  Int_t threshlvl = PartIsAboveThresh(StdHepPTLV*(IsInGev?1000.0:1.0));

  switch(StdHepPdg){
    case -13:
    case 13:{
      break;
    }
    //NC
    case -14:
    case 14:{
      return false;
      break;
    }
    case 22:{
      IncrementThreshArray(NAboveThresholdGammas,threshlvl);
      if(threshlvl > -1){
        NInEKinBinGammas[threshlvl]++;
      }
      break;
    }
    case 2212:{
      IncrementThreshArray(NAboveThresholdProtons,threshlvl);
      IncrementThreshArray(NAboveThresholdTrackable,threshlvl);
      if(threshlvl > -1){
        NInEKinBinProtons[threshlvl]++;
        NInEKinBinTrackable[threshlvl]++;
      }
      break;
    }
    case 2112:{
      IncrementThreshArray(NAboveThresholdNeutrons,threshlvl);
      IncrementThreshArray(NAboveThresholdNeutrals,threshlvl);
      if(threshlvl > -1){
        NAboveThresholdNeutrons[threshlvl]++;
        NAboveThresholdNeutrals[threshlvl]++;
      }
      break;
    }
    case 211:{
      IncrementThreshArray(NAboveThresholdChargedPions,threshlvl);
      IncrementThreshArray(NAboveThresholdPiPlus,threshlvl);
      IncrementThreshArray(NAboveThresholdTrackable,threshlvl);
      if(threshlvl > -1){
        NInEKinBinChargedPions[threshlvl]++;
        NInEKinBinPiPlus[threshlvl]++;
        NInEKinBinTrackable[threshlvl]++;
      }
      break;
    }
    case -211:{
      HandleHMTrackable(StdHepPTLV,StdHepP3Mod,StdHepPdg);
      HandleCPion(StdHepPTLV,StdHepP3Mod,StdHepPdg);
      IncrementThreshArray(NAboveThresholdChargedPions,threshlvl);
      IncrementThreshArray(NAboveThresholdPiMinus,threshlvl);
      IncrementThreshArray(NAboveThresholdTrackable,threshlvl);
      if(threshlvl > -1){
        NInEKinBinChargedPions[threshlvl]++;
        NInEKinBinPiMinus[threshlvl]++;
        NInEKinBinTrackable[threshlvl]++;
      }
      break;
    }
    case 111:{
      IncrementThreshArray(NAboveThresholdPiZero,threshlvl);
      IncrementThreshArray(NAboveThresholdNeutrals,threshlvl);
      if(threshlvl  > -1){
        NInEKinBinPiZero[threshlvl]++;
        NInEKinBinNeutrals[threshlvl]++;
      }
      break;
    }
    case 221:
    case 331:
    case 130:
    case 310:
    case 311:
    case 321:
    case -221:
    case -331:
    case -130:
    case -310:
    case -311:
    case -321:
    case 3122: {
      IncrementThreshArray(NAboveThresholdExotic,threshlvl);
      if(threshlvl > -1){
        NInEKinBinExotic[threshlvl]++;
      }
      break;
    }
    default:{
      if(false){
        std::cout << "[INFO]: Found Particle with PDG: " << StdHepPdg
          << std::endl;
      }
      return false;
      break;
    }
  }
  return true;
}

void TransversityVars::Finalise(){

  TransversityVarsB::Finalise();

  constexpr static Float_t GeVToMeV = 1000.0;
  constexpr static Float_t RadToDeg = 180.0/M_PI;
  constexpr static Float_t DegToRad = M_PI/180.0;

  if(IsInGev){
    FirstProton.Momentum *= GeVToMeV;
    FirstProton.FourMomentum *= GeVToMeV;
  }

//******************************************************************************
//                     Event Properties
//******************************************************************************

  ReconNuEnergy = GetReconNuEnergy(MuonNeutrino.FourMomentum.Vect(),
                                    Muon.FourMomentum,
                                    HMProton.FourMomentum);
  ReconTargetMass = GetReconTgtMass(ReconNuEnergy,
                                    Muon.FourMomentum,
                                    HMProton.FourMomentum);

//******************************************************************************
//                     Pertinent Particle Properties
//******************************************************************************

//Muon
  Muon_Pt_MeV = TransversityUtils::GetVectorInTPlane(Muon.FourMomentum.Vect(),
    MuonNeutrino.FourMomentum.Vect());

//First Proton
  FirstProton_PDG = FirstProton.PDG;
  FirstProton_4Mom_MeV = FirstProton.FourMomentum;
  FirstProton_StdHepPosition = FirstProton.StdHepPosition;

//Highest Momentum Proton
  HMProton_StdHepPosition = HMProton.StdHepPosition;

//StruckNucleon_3Mom_Recon
  StruckNucleon_3Mom_Recon_MeV = (Muon.FourMomentum + FirstProton.FourMomentum
    - MuonNeutrino.FourMomentum).Vect();
//PreFSINucleon_3Mom_Recon
  PreFSINucleon_3Mom_Recon_MeV = ( MuonNeutrino.FourMomentum
    + StruckNucleon.FourMomentum - Muon.FourMomentum ).Vect();

//******************************************************************************
//                       'Verse Variable Values
//******************************************************************************

//Projections
  TVector3 FirstProtonPt_MeV = GetVectorInTPlane(
    FirstProton.FourMomentum.Vect(), MuonNeutrino.FourMomentum.Vect());

//deltaphiT
  deltaphiT_FirstProton_deg = GetDeltaPhiT(
    Muon.FourMomentum.Vect(), FirstProton.FourMomentum.Vect(),
    MuonNeutrino.FourMomentum.Vect())*RadToDeg;

  deltaphiT_HMTrackable_deg = GetDeltaPhiT(
    Muon.FourMomentum.Vect(), HMTrackable.FourMomentum.Vect(),
    MuonNeutrino.FourMomentum.Vect())*RadToDeg;

//deltapT
  deltapT_FirstProton_MeV = GetDeltaPT(Muon_Pt_MeV,
                                       FirstProtonPt_MeV,
                                       MuonNeutrino.FourMomentum.Vect());

//Deltap
  Deltap_FirstProton_MeV = (Muon.FourMomentum + FirstProton.FourMomentum
    - MuonNeutrino.FourMomentum - StruckNucleon.FourMomentum);

//deltap
  deltap_FirstProton_MeV = (Muon.FourMomentum
    + FirstProton.FourMomentum - MuonNeutrino.FourMomentum);

//deltaalphaT
  deltaalphaT_FirstProton_deg =
    GetDeltaAlphaT(Muon_Pt_MeV,
                   FirstProtonPt_MeV,
                   MuonNeutrino.FourMomentum.Vect())*RadToDeg;

//deltap_TT
  if((HMCPion.Momentum > 1E-6) && (HMProton.Momentum > 1E-6)){
    deltap_TT_PionPDG = HMCPion.PDG;
  }

//******************************************************************************
//                      FS Particle Stuff
//******************************************************************************
  auto GetLastNonZeroIndex = [](TLorentzVector (&tlvs)[kMaxFSMomenta+1]) -> Int_t {
    for(Int_t i = 0; i < (kMaxFSMomenta+1); ++i){
      if(tlvs[i].Vect().Mag() < 1E-6){
        return i;
      }
    }
    return (kMaxFSMomenta+1);
  };

  //Gets the number of non zero momenta in the array ignoring the highest momenta
  // particle which is stored in separate variables.
  NOtherFSPiPlus4Momenta_MeV =
    GetLastNonZeroIndex(__OtherFSPiPlus4Momenta_MeV) - 1;
  NOtherFSProton4Momenta_MeV =
    GetLastNonZeroIndex(__OtherFSProton4Momenta_MeV) - 1;

  NOtherFSPiPlus4Momenta_MeV =
    NOtherFSPiPlus4Momenta_MeV < 0 ? 0 : NOtherFSPiPlus4Momenta_MeV;
  NOtherFSProton4Momenta_MeV =
    NOtherFSProton4Momenta_MeV < 0 ? 0 : NOtherFSProton4Momenta_MeV;

  auto SplitTLV = [](TLorentzVector (&tlvs)[kMaxFSMomenta+1], Double_t* X,
    Double_t* Y, Double_t* Z, Double_t* T, Int_t N = (kMaxFSMomenta+1)){
     //Need to ignore the first proton
     for(Int_t i = 1; i < N; ++i){
      X[i-1] = tlvs[i].X();
      Y[i-1] = tlvs[i].Y();
      Z[i-1] = tlvs[i].Z();
      T[i-1] = tlvs[i].T();
     }
  };

  SplitTLV( __OtherFSPiPlus4Momenta_MeV,
          OtherFSPiPlus4Momenta_MeV_X,
          OtherFSPiPlus4Momenta_MeV_Y,
          OtherFSPiPlus4Momenta_MeV_Z,
          OtherFSPiPlus4Momenta_MeV_T,
          NOtherFSPiPlus4Momenta_MeV);
  SplitTLV( __OtherFSProton4Momenta_MeV,
          OtherFSProton4Momenta_MeV_X,
          OtherFSProton4Momenta_MeV_Y,
          OtherFSProton4Momenta_MeV_Z,
          OtherFSProton4Momenta_MeV_T,
          NOtherFSProton4Momenta_MeV);

//******************************************************************************
//******************************************************************************
}

void TransversityVars::Reset(){

  TransversityVarsB::Reset();

//******************************************************************************
//                     Event Properties
//******************************************************************************

  ReconNuEnergy = 0;
  ReconTargetMass = 0;

//******************************************************************************
//                     Pertinent Particle Properties
//******************************************************************************

//Muon
  Muon_Pt_MeV = TVector3(0,0,0);

//First Proton
  FirstProton_PDG = 0;
  FirstProton_4Mom_MeV = TLorentzVector(0,0,0,0);
  FirstProton_StdHepPosition = -1;

//Highest Momentum Proton
  HMProton_StdHepPosition = -1;

//StruckNucleon_3Mom_Recon
  StruckNucleon_3Mom_Recon_MeV = TVector3(0,0,0);

//PreFSINucleon_3Mom_Recon
  PreFSINucleon_3Mom_Recon_MeV = TVector3(0,0,0);

//******************************************************************************
//                       'Verse Variable Values
//******************************************************************************

//deltaphiT
  deltaphiT_FirstProton_deg = 0;
  deltaphiT_HMTrackable_deg = 0;

//DeltaPt
  deltapT_FirstProton_MeV = TVector3(0,0,0);

//Deltap_
  Deltap_FirstProton_MeV = TLorentzVector(0,0,0,0);

//deltap
  deltap_FirstProton_MeV = TLorentzVector(0,0,0,0);

//DeltaAlphat
  deltaalphaT_FirstProton_deg = 0;

//deltap_TT
  deltap_TT_PionPDG = 0;

//******************************************************************************
//                       Final State Particles
//******************************************************************************

  for(int i = 0; i < kNThreshMax; ++i){
    NAboveThresholdProtons[i] = 0;
    NAboveThresholdGammas[i] = 0;
    NAboveThresholdPiPlus[i] = 0;
    NAboveThresholdPiMinus[i] = 0;
    NAboveThresholdChargedPions[i] = 0;
    NAboveThresholdTrackable[i] = 0;
    NAboveThresholdNeutrons[i] = 0;
    NAboveThresholdPiZero[i] = 0;
    NAboveThresholdNeutrals[i] = 0;
    NAboveThresholdExotic[i] = 0;

    NInEKinBinProtons[i] = 0;
    NInEKinBinGammas[i] = 0;
    NInEKinBinPiPlus[i] = 0;
    NInEKinBinPiMinus[i] = 0;
    NInEKinBinChargedPions[i] = 0;
    NInEKinBinTrackable[i] = 0;
    NInEKinBinNeutrons[i] = 0;
    NInEKinBinPiZero[i] = 0;
    NInEKinBinNeutrals[i] = 0;
    NInEKinBinExotic[i] = 0;
  }

//******************************************************************************
//                      FS Particle Stuff
//******************************************************************************
  ClearArray(__OtherFSPiPlus4Momenta_MeV);
  ClearArray(__OtherFSProton4Momenta_MeV);
  ClearArray(OtherFSPiPlus4Momenta_MeV_X,kMaxFSMomenta);
  ClearArray(OtherFSPiPlus4Momenta_MeV_Y,kMaxFSMomenta);
  ClearArray(OtherFSPiPlus4Momenta_MeV_Z,kMaxFSMomenta);
  ClearArray(OtherFSPiPlus4Momenta_MeV_T,kMaxFSMomenta);
  ClearArray(OtherFSProton4Momenta_MeV_X,kMaxFSMomenta);
  ClearArray(OtherFSProton4Momenta_MeV_Y,kMaxFSMomenta);
  ClearArray(OtherFSProton4Momenta_MeV_Z,kMaxFSMomenta);
  ClearArray(OtherFSProton4Momenta_MeV_T,kMaxFSMomenta);

  NOtherFSPiPlus4Momenta_MeV = 0;
  NOtherFSProton4Momenta_MeV = 0;

//******************************************************************************
//                       Others and Transients
//******************************************************************************

  //Transients
  FirstProton.Reset();

//******************************************************************************
//******************************************************************************
}

void TransversityVars::AddBranches(TTree* tree){

  TransversityVarsB::AddBranches(tree);

//******************************************************************************
//                     Event Properties
//******************************************************************************
  tree->Branch("GeneratorName", &_GeneratorName);

  tree->Branch("ReconNuEnergy", &ReconNuEnergy, "ReconNuEnergy/D");
  tree->Branch("ReconTargetMass", &ReconTargetMass, "ReconTargetMass/D");

//******************************************************************************
//                     Pertinent Particle Properties
//******************************************************************************

//Muon
  tree->Branch("Muon_Pt_MeV", &_Muon_Pt_MeV);

//First Proton
  tree->Branch("FirstProton_PDG", &FirstProton_PDG,
    "FirstProton_PDG/I");
  tree->Branch("FirstProton_4Mom_MeV", &_FirstProton_4Mom_MeV);

  tree->Branch("FirstProton_StdHepPosition", &FirstProton_StdHepPosition,
    "FirstProton_StdHepPosition/I");

//Highest Momentum Proton
  tree->Branch("HMProton_StdHepPosition", &HMProton_StdHepPosition,
    "HMProton_StdHepPosition/I");

//StruckNucleon_3Mom_Recon
  tree->Branch("StruckNucleon_3Mom_Recon_MeV", &_StruckNucleon_3Mom_Recon_MeV);

//PreFSINucleon_3Mom_Recon
  tree->Branch("PreFSINucleon_3Mom_Recon_MeV", &_PreFSINucleon_3Mom_Recon_MeV);

//******************************************************************************
//                       'Verse Variable Values
//******************************************************************************

//deltaphiT
  tree->Branch("deltaphiT_FirstProton_deg", &deltaphiT_FirstProton_deg,
    "deltaphiT_FirstProton_deg/D");
  tree->Branch("deltaphiT_HMTrackable_deg", &deltaphiT_HMTrackable_deg,
    "deltaphiT_HMTrackable_deg/D");

//DeltaPt
  tree->Branch("deltapT_FirstProton_MeV", &_deltapT_FirstProton_MeV);

//Deltap_
  tree->Branch("Deltap_FirstProton_MeV", &_Deltap_FirstProton_MeV);

//deltap
  tree->Branch("deltap_FirstProton_MeV", &_deltap_FirstProton_MeV);

//DeltaAlphat
  tree->Branch("deltaalphaT_FirstProton_deg", &deltaalphaT_FirstProton_deg,
    "deltaalphaT_FirstProton_deg/D");

//deltap_TT
  tree->Branch("deltap_TT_PionPDG", &deltap_TT_PionPDG, "deltap_TT_PionPDG/D");

//******************************************************************************
//                       Final State Particles
//******************************************************************************

  tree->Branch("NThresh", &NThresh, "NThresh/I");
  tree->Branch("Threshs_MeV", Threshs_MeV, "Threshs_MeV[NThresh]/I");

  tree->Branch("NAboveThresholdProtons", NAboveThresholdProtons,
    "NAboveThresholdProtons[NThresh]/I");
  tree->Branch("NAboveThresholdGammas", NAboveThresholdGammas,
    "NAboveThresholdGammas[NThresh]/I");
  tree->Branch("NAboveThresholdPiPlus", NAboveThresholdPiPlus,
    "NAboveThresholdPiPlus[NThresh]/I");
  tree->Branch("NAboveThresholdPiMinus", NAboveThresholdPiMinus,
    "NAboveThresholdPiMinus[NThresh]/I");
  tree->Branch("NAboveThresholdChargedPions", NAboveThresholdChargedPions,
    "NAboveThresholdChargedPions[NThresh]/I");
  tree->Branch("NAboveThresholdTrackable", NAboveThresholdTrackable,
    "NAboveThresholdTrackable[NThresh]/I");
  tree->Branch("NAboveThresholdNeutrons", NAboveThresholdNeutrons,
    "NAboveThresholdNeutrons[NThresh]/I");
  tree->Branch("NAboveThresholdPiZero", NAboveThresholdPiZero,
    "NAboveThresholdPiZero[NThresh]/I");
  tree->Branch("NAboveThresholdNeutrals", NAboveThresholdNeutrals,
    "NAboveThresholdNeutrals[NThresh]/I");
  tree->Branch("NAboveThresholdExotic", NAboveThresholdExotic,
    "NAboveThresholdExotic[NThresh]/I");

  tree->Branch("NInEKinBinProtons", NInEKinBinProtons,
    "NInEKinBinProtons[NThresh]/I");
  tree->Branch("NInEKinBinGammas", NInEKinBinGammas,
    "NInEKinBinGammas[NThresh]/I");
  tree->Branch("NInEKinBinPiPlus", NInEKinBinPiPlus,
    "NInEKinBinPiPlus[NThresh]/I");
  tree->Branch("NInEKinBinPiMinus", NInEKinBinPiMinus,
    "NInEKinBinPiMinus[NThresh]/I");
  tree->Branch("NInEKinBinChargedPions", NInEKinBinChargedPions,
    "NInEKinBinChargedPions[NThresh]/I");
  tree->Branch("NInEKinBinTrackable", NInEKinBinTrackable,
    "NInEKinBinTrackable[NThresh]/I");
  tree->Branch("NInEKinBinNeutrons", NInEKinBinNeutrons,
    "NInEKinBinNeutrons[NThresh]/I");
  tree->Branch("NInEKinBinPiZero", NInEKinBinPiZero,
    "NInEKinBinPiZero[NThresh]/I");
  tree->Branch("NInEKinBinNeutrals", NInEKinBinNeutrals,
    "NInEKinBinNeutrals[NThresh]/I");
  tree->Branch("NInEKinBinExotic", NInEKinBinExotic,
    "NInEKinBinExotic[NThresh]/I");

//******************************************************************************
//                      FS Particle Stuff
//******************************************************************************

  tree->Branch("NOtherFSPiPlus4Momenta_MeV", &NOtherFSPiPlus4Momenta_MeV,
    "NOtherFSPiPlus4Momenta_MeV/I");
  tree->Branch("NOtherFSProton4Momenta_MeV", &NOtherFSProton4Momenta_MeV,
    "NOtherFSProton4Momenta_MeV/I");

  tree->Branch("OtherFSPiPlus4Momenta_MeV_X", OtherFSPiPlus4Momenta_MeV_X,
    "OtherFSPiPlus4Momenta_MeV_X[NOtherFSPiPlus4Momenta_MeV]/D");
  tree->Branch("OtherFSPiPlus4Momenta_MeV_Y", OtherFSPiPlus4Momenta_MeV_Y,
    "OtherFSPiPlus4Momenta_MeV_Y[NOtherFSPiPlus4Momenta_MeV]/D");
  tree->Branch("OtherFSPiPlus4Momenta_MeV_Z", OtherFSPiPlus4Momenta_MeV_Z,
    "OtherFSPiPlus4Momenta_MeV_Z[NOtherFSPiPlus4Momenta_MeV]/D");
  tree->Branch("OtherFSPiPlus4Momenta_MeV_T", OtherFSPiPlus4Momenta_MeV_T,
    "OtherFSPiPlus4Momenta_MeV_T[NOtherFSPiPlus4Momenta_MeV]/D");

  tree->Branch("OtherFSProton4Momenta_MeV_X", OtherFSProton4Momenta_MeV_X,
    "OtherFSProton4Momenta_MeV_X[NOtherFSProton4Momenta_MeV]/D");
  tree->Branch("OtherFSProton4Momenta_MeV_Y", OtherFSProton4Momenta_MeV_Y,
    "OtherFSProton4Momenta_MeV_Y[NOtherFSProton4Momenta_MeV]/D");
  tree->Branch("OtherFSProton4Momenta_MeV_Z", OtherFSProton4Momenta_MeV_Z,
    "OtherFSProton4Momenta_MeV_Z[NOtherFSProton4Momenta_MeV]/D");
  tree->Branch("OtherFSProton4Momenta_MeV_T", OtherFSProton4Momenta_MeV_T,
    "OtherFSProton4Momenta_MeV_T[NOtherFSProton4Momenta_MeV]/D");

//******************************************************************************
//                       Others and Transients
//******************************************************************************
//******************************************************************************
//******************************************************************************
}

struct Proxy {

  Proxy(TransversityVarsB* pf) : proxyFor(pf),
    IncNeutrino_4Mom_MeV(nullptr),
    StruckNucleon_4Mom_MeV(nullptr),
    Muon_4Mom_MeV(nullptr),
    HMProton_4Mom_MeV(nullptr),
    HMCPion_4Mom_MeV(nullptr),
    HMTrackable_4Mom_MeV(nullptr),
    Deltap_HMProton_MeV(nullptr),
    deltap_HMProton_MeV(nullptr),
    deltapT_HMProton_MeV(nullptr),
    HMProtonPion_3Mom_MeV(nullptr),
    deltapT_HMProtonPion_MeV(nullptr),
    Deltap_HMProtonPion_MeV(nullptr){}

  TransversityVarsB* proxyFor;

//******************************************************************************
//                     Pertinent Particle Properties
//******************************************************************************

//Neutrino
  TLorentzVector *IncNeutrino_4Mom_MeV;

//Struck Nucleon
  TLorentzVector *StruckNucleon_4Mom_MeV;

//Muon
  TLorentzVector *Muon_4Mom_MeV;

//Highest Momentum Proton
  TLorentzVector *HMProton_4Mom_MeV;

//Highest Momentum Charged Pion
  TLorentzVector *HMCPion_4Mom_MeV;

//Highest Momentum Trackable
  TLorentzVector *HMTrackable_4Mom_MeV;

//******************************************************************************
//                       'Verse Variable Values
//******************************************************************************

//Deltap
  TLorentzVector *Deltap_HMProton_MeV;

//deltap
  TLorentzVector *deltap_HMProton_MeV;

//deltapT
  TVector3 *deltapT_HMProton_MeV;

//ProtonPion Combo Platter
  TVector3 *HMProtonPion_3Mom_MeV;
  TVector3 *deltapT_HMProtonPion_MeV;
  TLorentzVector *Deltap_HMProtonPion_MeV;

//******************************************************************************
//******************************************************************************
};
constexpr size_t kMaxProxies = 100;
std::vector<Proxy> _prxy;

void UnsetBranchAddressesTransversityVarsB(TTree* tree,
  TransversityVarsB const * const tvb){

  tree->ResetBranchAddresses();

  for(auto pr_it = _prxy.begin(); pr_it < _prxy.end();){
    if(tvb == pr_it->proxyFor){
      pr_it = _prxy.erase(pr_it);
      continue;
    }
    ++pr_it;
  }
  delete tvb;
}

void SetBranchAddressesTransversityVarsB(TTree* tree, TransversityVarsB* tvb){
//******************************************************************************
//                     Event Properties
//******************************************************************************

  if(kMaxProxies == _prxy.size()){
    throw std::overflow_error("Attempted to set up more than 100 TTrees with "
      "separate instances of TransversityVars. Are you attempting to do this"
      "for each entry?");

  }
  _prxy.push_back(Proxy(tvb));

//Generator reaction code
  tree->SetBranchAddress("NeutConventionReactionCode",
    &tvb->NeutConventionReactionCode);

//******************************************************************************
//                     Pertinent Particle Properties
//******************************************************************************

//Neutrino
  tree->SetBranchAddress("IncNeutrino_PDG",&tvb->IncNeutrino_PDG);
  _prxy.back().IncNeutrino_4Mom_MeV = &tvb->IncNeutrino_4Mom_MeV;
  tree->SetBranchAddress("IncNeutrino_4Mom_MeV", &_prxy.back().IncNeutrino_4Mom_MeV);

//Struck Nucleon
  tree->SetBranchAddress("StruckNucleonPDG",&tvb->StruckNucleonPDG);
  _prxy.back().StruckNucleon_4Mom_MeV = &tvb->StruckNucleon_4Mom_MeV;
  tree->SetBranchAddress("StruckNucleon_4Mom_MeV", &_prxy.back().StruckNucleon_4Mom_MeV);

//Muon
  tree->SetBranchAddress("Muon_PDG",&tvb->Muon_PDG);
  _prxy.back().Muon_4Mom_MeV = &tvb->Muon_4Mom_MeV;
  tree->SetBranchAddress("Muon_4Mom_MeV", &_prxy.back().Muon_4Mom_MeV);

//Highest Momentum Proton
  tree->SetBranchAddress("HMProton_PDG",&tvb->HMProton_PDG);
  _prxy.back().HMProton_4Mom_MeV = &tvb->HMProton_4Mom_MeV;
  tree->SetBranchAddress("HMProton_4Mom_MeV", &_prxy.back().HMProton_4Mom_MeV);

//Highest Momentum Charged Pion
  tree->SetBranchAddress("HMCPion_PDG",&tvb->HMCPion_PDG);
  _prxy.back().HMCPion_4Mom_MeV = &tvb->HMCPion_4Mom_MeV;
  tree->SetBranchAddress("HMCPion_4Mom_MeV", &_prxy.back().HMCPion_4Mom_MeV);

//Highest Momentum Trackable
  tree->SetBranchAddress("HMTrackable_PDG",&tvb->HMTrackable_PDG);
  _prxy.back().HMTrackable_4Mom_MeV = &tvb->HMTrackable_4Mom_MeV;
  tree->SetBranchAddress("HMTrackable_4Mom_MeV", &_prxy.back().HMTrackable_4Mom_MeV);

//******************************************************************************
//                       'Verse Variable Values
//******************************************************************************

//deltaphiT
  tree->SetBranchAddress("deltaphiT_HMProton_deg", &tvb->deltaphiT_HMProton_deg);

//deltapT
  _prxy.back().deltapT_HMProton_MeV = &tvb->deltapT_HMProton_MeV;
  tree->SetBranchAddress("deltapT_HMProton_MeV", &_prxy.back().deltapT_HMProton_MeV);

//Deltap
  _prxy.back().Deltap_HMProton_MeV = &tvb->Deltap_HMProton_MeV;
  tree->SetBranchAddress("Deltap_HMProton_MeV", &_prxy.back().Deltap_HMProton_MeV);

//deltap
  _prxy.back().deltap_HMProton_MeV = &tvb->deltap_HMProton_MeV;
  tree->SetBranchAddress("deltap_HMProton_MeV", &_prxy.back().deltap_HMProton_MeV);

//deltaalphaT
  tree->SetBranchAddress("deltaalphaT_HMProton_deg", &tvb->deltaalphaT_HMProton_deg);

//deltap_TT
  tree->SetBranchAddress("deltap_TT", &tvb->deltap_TT);

//ProtonPion Combo Platter
  _prxy.back().HMProtonPion_3Mom_MeV = &tvb->HMProtonPion_3Mom_MeV;
  tree->SetBranchAddress("HMProtonPion_3Mom_MeV", &_prxy.back().HMProtonPion_3Mom_MeV);
  tree->SetBranchAddress("deltaphiT_HMProtonPion_deg", &tvb->deltaphiT_HMProtonPion_deg);
  _prxy.back().deltapT_HMProtonPion_MeV = &tvb->deltapT_HMProtonPion_MeV;
  tree->SetBranchAddress("deltapT_HMProtonPion_MeV", &_prxy.back().deltapT_HMProtonPion_MeV);
  tree->SetBranchAddress("deltaalphaT_HMProtonPion_deg", &tvb->deltaalphaT_HMProtonPion_deg);
  _prxy.back().Deltap_HMProtonPion_MeV = &tvb->Deltap_HMProtonPion_MeV;
  tree->SetBranchAddress("Deltap_HMProtonPion_MeV", &_prxy.back().Deltap_HMProtonPion_MeV);

//******************************************************************************
//                       Subsequent Species Sums
//******************************************************************************

  tree->SetBranchAddress("NFinalStateParticles", &tvb->NFinalStateParticles);

  tree->SetBranchAddress("NProtons", &tvb->NProtons);
  tree->SetBranchAddress("NGammas", &tvb->NGammas);
  tree->SetBranchAddress("NNeutrons", &tvb->NNeutrons);
  tree->SetBranchAddress("NPiPlus", &tvb->NPiPlus);
  tree->SetBranchAddress("NPiZero", &tvb->NPiZero);
  tree->SetBranchAddress("NPiMinus", &tvb->NPiMinus);
  tree->SetBranchAddress("NPions", &tvb->NPions);
  tree->SetBranchAddress("NChargedPions", &tvb->NChargedPions);
  tree->SetBranchAddress("NOtherParticles", &tvb->NOtherParticles);

//******************************************************************************
//                       Tangible Target Traits
//******************************************************************************

  tree->SetBranchAddress("TargetPDG", &tvb->TargetPDG);
  tree->SetBranchAddress("TargetZ", &tvb->TargetZ);

//******************************************************************************
//                       Others and Transients
//******************************************************************************

  tree->SetBranchAddress("CCQ2", &tvb->CCQ2);

//******************************************************************************
//******************************************************************************
}
TransversityVarsB const * const MakeReadingTransversityVarsB(TTree* tree){
  TransversityVarsB* rtn = new TransversityVarsB();
  SetBranchAddressesTransversityVarsB(tree,rtn);
  return rtn;
}
TransversityVars* MakeReadingTransversityVars(TTree* tree){
  TransversityVars* rtn = new TransversityVars();
  SetBranchAddressesTransversityVarsB(tree,rtn);
  //Add other branches here.
  return rtn;
}
