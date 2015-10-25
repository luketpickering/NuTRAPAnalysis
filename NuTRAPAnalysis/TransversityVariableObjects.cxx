#include <iostream>

#include "TRandom.h"
#include "TF1.h"

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

ClassImp(TransversityVars);

TransversityVars::~TransversityVars(){
  delete OtherFSPiPlus4Momenta_MeV_X;
  delete OtherFSPiPlus4Momenta_MeV_Y;
  delete OtherFSPiPlus4Momenta_MeV_Z;
  delete OtherFSPiPlus4Momenta_MeV_T;
  delete OtherFSProton4Momenta_MeV_X;
  delete OtherFSProton4Momenta_MeV_Y;
  delete OtherFSProton4Momenta_MeV_Z;
  delete OtherFSProton4Momenta_MeV_T;
}

TransversityVars::TransversityVars(
  bool InGeV, Double_t TargetBE_MeV, Int_t NThresh, Int_t* Threshs_MeV,
  TString generatorName){
  GeneratorName = generatorName;

  this->NThresh = NThresh;
  this->Threshs_MeV = Threshs_MeV;
  IsInGev = InGeV;
  this->TargetBE_MeV = TargetBE_MeV;

  OtherFSPiPlus4Momenta_MeV_X = new Double_t[kMaxFSMomenta];
  OtherFSPiPlus4Momenta_MeV_Y = new Double_t[kMaxFSMomenta];
  OtherFSPiPlus4Momenta_MeV_Z = new Double_t[kMaxFSMomenta];
  OtherFSPiPlus4Momenta_MeV_T = new Double_t[kMaxFSMomenta];
  OtherFSProton4Momenta_MeV_X = new Double_t[kMaxFSMomenta];
  OtherFSProton4Momenta_MeV_Y = new Double_t[kMaxFSMomenta];
  OtherFSProton4Momenta_MeV_Z = new Double_t[kMaxFSMomenta];
  OtherFSProton4Momenta_MeV_T = new Double_t[kMaxFSMomenta];

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

void TransversityVars::HandleHMTrackable(TLorentzVector &StdHepPTLV,
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

void TransversityVars::HandleStruckNucleon(TLorentzVector &StdHepPTLV,
  Double_t &StdHepP3Mod, Int_t pdg){

  StruckNucleon.Momentum = StdHepP3Mod;
  StruckNucleon.FourMomentum = StdHepPTLV;
  StruckNucleon.PDG = pdg;
}


void TransversityVars::HandleRescat(Int_t PDG, Int_t RescatCode){
  if(PDG == 2212){
    switch(RescatCode){
      case 1:{
        ProtonRescat_contains_NoInt = true;
        break;
      }
      case 2:{
        ProtonRescat_contains_chrgEx = true;
        break;
      }
      case 3:{
        ProtonRescat_contains_elastic = true;
        break;
      }
      case 4:{
        ProtonRescat_contains_inelastic = true;
        break;
      }
      case 5:{
        ProtonRescat_contains_knockout = true;
        break;
      }
    }
  }
}

bool TransversityVars::HandleStdHepParticle(
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

  Int_t threshlvl = PartIsAboveThresh(StdHepPTLV*(IsInGev?1000.0:1.0));

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
      IncrementThreshArray(NAboveThresholdGammas,threshlvl);
      if(threshlvl > -1){
        NInEKinBinGammas[threshlvl]++;
      }
      NFinalStateParticles++;
      break;
    }
    case 2212:{
      HandleProton(StdHepPTLV, StdHepP3Mod, StdHepPosition);
      HandleHMTrackable(StdHepPTLV,StdHepP3Mod,StdHepPdg);
      NProtons++;
      IncrementThreshArray(NAboveThresholdProtons,threshlvl);
      IncrementThreshArray(NAboveThresholdTrackable,threshlvl);
      if(threshlvl > -1){
        NInEKinBinProtons[threshlvl]++;
        NInEKinBinTrackable[threshlvl]++;
      }
      NFinalStateParticles++;
      break;
    }
    case 2112:{
      IncrementThreshArray(NAboveThresholdNeutrons,threshlvl);
      IncrementThreshArray(NAboveThresholdNeutrals,threshlvl);
      NNeutrons++;
      NFinalStateParticles++;
      if(threshlvl > -1){
        NAboveThresholdNeutrons[threshlvl]++;
        NAboveThresholdNeutrals[threshlvl]++;
      }
      break;
    }
    case 211:{
      HandleHMTrackable(StdHepPTLV,StdHepP3Mod,StdHepPdg);
      HandleCPion(StdHepPTLV,StdHepP3Mod,StdHepPdg);
      IncrementThreshArray(NAboveThresholdChargedPions,threshlvl);
      IncrementThreshArray(NAboveThresholdPiPlus,threshlvl);
      IncrementThreshArray(NAboveThresholdTrackable,threshlvl);
      NPiPlus++;
      NChargedPions++;
      NPions++;
      if(threshlvl > -1){
        NInEKinBinChargedPions[threshlvl]++;
        NInEKinBinPiPlus[threshlvl]++;
        NInEKinBinTrackable[threshlvl]++;
      }
      NFinalStateParticles++;
      break;
    }
    case -211:{
      HandleHMTrackable(StdHepPTLV,StdHepP3Mod,StdHepPdg);
      HandleCPion(StdHepPTLV,StdHepP3Mod,StdHepPdg);
      IncrementThreshArray(NAboveThresholdChargedPions,threshlvl);
      IncrementThreshArray(NAboveThresholdPiMinus,threshlvl);
      IncrementThreshArray(NAboveThresholdTrackable,threshlvl);
      NPiMinus++;
      NChargedPions++;
      NPions++;
      if(threshlvl > -1){
        NInEKinBinChargedPions[threshlvl]++;
        NInEKinBinPiMinus[threshlvl]++;
        NInEKinBinTrackable[threshlvl]++;
      }
      NFinalStateParticles++;
      break;
    }
    case 111:{
      IncrementThreshArray(NAboveThresholdPiZero,threshlvl);
      IncrementThreshArray(NAboveThresholdNeutrals,threshlvl);
      NPiZero++;
      NPions++;
      if(threshlvl  > -1){
        NInEKinBinPiZero[threshlvl]++;
        NInEKinBinNeutrals[threshlvl]++;
      }
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
      IncrementThreshArray(NAboveThresholdExotic,threshlvl);
      if(threshlvl > -1){
        NInEKinBinExotic[threshlvl]++;
      }
      NFinalStateParticles++;
      break;
    }
    default:{
      NOtherParticles++;
      NFinalStateParticles++;
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
    FirstProton.Momentum *= GeVToMeV;
    FirstProton.FourMomentum *= GeVToMeV;
    HMTrackable.Momentum *= GeVToMeV;
    HMTrackable.FourMomentum *= GeVToMeV;
    HMCPion.Momentum *= GeVToMeV;
    HMCPion.FourMomentum *= GeVToMeV;
    StruckNucleon.Momentum *= GeVToMeV;
    StruckNucleon.FourMomentum *= GeVToMeV;
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

  if(TargetBE_MeV != 0xdeadbeef){
    ERec_CCQE = GetCCQEERec(MuonNeutrino.FourMomentum.Vect(),
      Muon.FourMomentum, TargetBE_MeV);
    ERec_E_PMu = GetMuPERec(Muon.FourMomentum,
      HMProton.FourMomentum, TargetBE_MeV);
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
  Muon_Pt_MeV = TransversityUtils::GetVectorInTPlane(Muon.FourMomentum.Vect(),
    MuonNeutrino.FourMomentum.Vect());

//First Proton
  FirstProton_PDG = FirstProton.PDG;
  FirstProton_4Mom_MeV = FirstProton.FourMomentum;
  FirstProton_StdHepPosition = FirstProton.StdHepPosition;

//Highest Momentum Proton
  HMProton_PDG = HMProton.PDG;
  HMProton_4Mom_MeV = HMProton.FourMomentum;
  HMProton_StdHepPosition = HMProton.StdHepPosition;

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
  TVector3 FirstProtonPt_MeV = GetVectorInTPlane(
    FirstProton.FourMomentum.Vect(), MuonNeutrino.FourMomentum.Vect());
  TVector3 HMTrackablePt_MeV = GetVectorInTPlane(
    HMTrackable.FourMomentum.Vect(), MuonNeutrino.FourMomentum.Vect());

//DeltaPhiT
  DeltaPhiT_HMProton_deg = GetDeltaPhiT(
    Muon.FourMomentum.Vect(), HMProton.FourMomentum.Vect(),
    MuonNeutrino.FourMomentum.Vect())*RadToDeg;

  DeltaPhiT_FirstProton_deg = GetDeltaPhiT(
    Muon.FourMomentum.Vect(), FirstProton.FourMomentum.Vect(),
    MuonNeutrino.FourMomentum.Vect())*RadToDeg;

  DeltaPhiT_HMTrackable_deg = GetDeltaPhiT(
    Muon.FourMomentum.Vect(), HMTrackable.FourMomentum.Vect(),
    MuonNeutrino.FourMomentum.Vect())*RadToDeg;

//DeltaPt
  DeltaPT_HMProton_MeV = GetDeltaPT(Muon_Pt_MeV,
                                    HMProtonPt_MeV,
                                    MuonNeutrino.FourMomentum.Vect());
  DeltaPT_FirstProton_MeV = GetDeltaPT(Muon_Pt_MeV,
                                       FirstProtonPt_MeV,
                                       MuonNeutrino.FourMomentum.Vect());

  DeltaPTotal_HMProton_MeV = (Muon.FourMomentum + HMProton.FourMomentum
    - MuonNeutrino.FourMomentum - StruckNucleon.FourMomentum);
  DeltaPTotal_FirstProton_MeV = (Muon.FourMomentum + FirstProton.FourMomentum
    - MuonNeutrino.FourMomentum - StruckNucleon.FourMomentum);

  DeltaPNuclearEffect_HMProton_MeV = (Muon.FourMomentum + HMProton.FourMomentum
    - MuonNeutrino.FourMomentum);
  DeltaPNuclearEffect_FirstProton_MeV = (Muon.FourMomentum
    + FirstProton.FourMomentum - MuonNeutrino.FourMomentum);

  DeltaPProton_MeV = (Muon.FourMomentum
    - MuonNeutrino.FourMomentum - StruckNucleon.FourMomentum).Vect();

  DeltaPTargetNucleon_MeV = (Muon.FourMomentum + FirstProton.FourMomentum
    - MuonNeutrino.FourMomentum).Vect();

//DeltaAlphat
  DeltaAlphaT_HMProton_deg = GetDeltaAlphaT(Muon_Pt_MeV,
                                            HMProtonPt_MeV,
                                            MuonNeutrino.FourMomentum.Vect())*RadToDeg;
  DeltaAlphaT_FirstProton_deg =
    GetDeltaAlphaT(Muon_Pt_MeV,
                   FirstProtonPt_MeV,
                   MuonNeutrino.FourMomentum.Vect())*RadToDeg;
//DeltaP_TT
  if((HMCPion.Momentum > 1E-3) && (HMProton.Momentum > 1E-3)){
    DeltaP_TT = TransversityUtils::GetDeltaPTT(
      Muon.FourMomentum.Vect(),
      HMCPion.FourMomentum.Vect(),
      HMProton.FourMomentum.Vect(),
      MuonNeutrino.FourMomentum.Vect());
    DeltaP_TT_PionPDG = HMCPion.PDG;

    HMProtonPion_3Mom_MeV = HMProton.FourMomentum.Vect() +
      HMCPion.FourMomentum.Vect();

    TVector3 HMProtonPionPt_MeV = GetVectorInTPlane(HMProtonPion_3Mom_MeV,
      MuonNeutrino.FourMomentum.Vect());

    DeltaPhiT_HMProtonPion_deg = GetDeltaPhiT(
                    Muon.FourMomentum.Vect(),
                    HMProtonPion_3Mom_MeV,
                    MuonNeutrino.FourMomentum.Vect())*RadToDeg;
    DeltaPT_HMProtonPion_MeV = GetDeltaPT(Muon_Pt_MeV,
                    HMProtonPionPt_MeV,
                    MuonNeutrino.FourMomentum.Vect());
    DeltaAlphaT_HMProtonPion_deg = GetDeltaAlphaT(Muon_Pt_MeV,
                    HMProtonPionPt_MeV,
                    MuonNeutrino.FourMomentum.Vect())*RadToDeg;

    DeltaPTotal_HMProtonPion_MeV = (Muon.FourMomentum + HMProton.FourMomentum
      + HMCPion.FourMomentum - MuonNeutrino.FourMomentum
      - StruckNucleon.FourMomentum);
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
//                       Others and Transients
//******************************************************************************

  CCQE_deltaP_MeV = Muon.FourMomentum + HMProton.FourMomentum -
    MuonNeutrino.FourMomentum - StruckNucleon.FourMomentum;

  CCQ2 = (Muon.FourMomentum - MuonNeutrino.FourMomentum).Mag2();

//******************************************************************************
//******************************************************************************
}

void TransversityVars::Reset(){

//******************************************************************************
//                     Event Properties
//******************************************************************************

//Generator reaction code
  NeutConventionReactionCode = 0;
  ReconNuEnergy = 0;
  ReconTargetMass = 0;

  ERec_CCQE = 0;
  ERec_E_PMu = 0;

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
  Muon_Pt_MeV = TVector3(0,0,0);

//First Proton
  FirstProton_PDG = 0;
  FirstProton_4Mom_MeV = TLorentzVector(0,0,0,0);
  FirstProton_StdHepPosition = -1;

//Highest Momentum Proton
  HMProton_PDG = 0;
  HMProton_4Mom_MeV = TLorentzVector(0,0,0,0);
  HMProton_StdHepPosition = -1;

//Highest Momentum Charged Pion
  HMCPion_PDG = 0;
  HMCPion_4Mom_MeV = TLorentzVector(0,0,0,0);

//Highest Momentum Trackable
  HMTrackable_PDG = 0;
  HMTrackable_4Mom_MeV = TLorentzVector(0,0,0,0);

//******************************************************************************
//                       'Verse Variable Values
//******************************************************************************

//DeltaPhiT
  DeltaPhiT_HMProton_deg = 0;
  DeltaPhiT_FirstProton_deg = 0;
  DeltaPhiT_HMTrackable_deg = 0;

//DeltaPt
  DeltaPT_HMProton_MeV = TVector3(0,0,0);
  DeltaPT_FirstProton_MeV = TVector3(0,0,0);

//DeltaPTotal
  DeltaPTotal_HMProton_MeV = TLorentzVector(0,0,0,0);
  DeltaPTotal_FirstProton_MeV = TLorentzVector(0,0,0,0);

//DeltaPNuclearEffect
  DeltaPNuclearEffect_HMProton_MeV = TLorentzVector(0,0,0,0);
  DeltaPNuclearEffect_FirstProton_MeV = TLorentzVector(0,0,0,0);

//DeltaPProton
  DeltaPProton_MeV = TVector3(0,0,0);
  DeltaPTargetNucleon_MeV = TVector3(0,0,0);

//DeltaAlphat
  DeltaAlphaT_HMProton_deg = 0;
  DeltaAlphaT_FirstProton_deg = 0;
//DeltaP_TT
  DeltaP_TT_PionPDG = 0;
  DeltaP_TT = 0;

//ProtonPion Combo Platter
  HMProtonPion_3Mom_MeV = TVector3(0,0,0);
  DeltaPhiT_HMProton_deg = 0;
  DeltaPT_HMProton_MeV = TVector3(0,0,0);
  DeltaAlphaT_HMProton_deg = 0;
  DeltaPTotal_HMProtonPion_MeV = TLorentzVector(0,0,0,0);

//******************************************************************************
//                       Final State Particles
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
//                       Tangible Target Traits
//******************************************************************************

  TargetPDG = 0;
  TargetZ = 0;

//******************************************************************************
//                       Others and Transients
//******************************************************************************

  CCQE_deltaP_MeV = TLorentzVector(0,0,0,0);

  CCQ2 = 0xdeadbeef;

  ProtonRescat_contains_NoInt = false;
  ProtonRescat_contains_chrgEx = false;
  ProtonRescat_contains_elastic = false;
  ProtonRescat_contains_inelastic = false;
  ProtonRescat_contains_knockout = false;

  //Transients
  Muon.Reset();
  MuonNeutrino.Reset();
  StruckNucleon.Reset();
  HMProton.Reset();
  HMCPion.Reset();
  FirstProton.Reset();
  HMTrackable.Reset();

//******************************************************************************
//******************************************************************************

}

ClassImp(TransversityVarsLite);

TransversityVarsLite::~TransversityVarsLite(){}

TransversityVarsLite::TransversityVarsLite(
  bool InGeV){

  IsInGev = InGeV;

  Reset();
}

void TransversityVarsLite::HandleProton(TLorentzVector &StdHepPTLV,
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
}

void TransversityVarsLite::HandleCPion(TLorentzVector &StdHepPTLV,
  Double_t &StdHepP3Mod, Int_t pdg){

  if(StdHepP3Mod > HMCPion.Momentum){
    HMCPion.Momentum = StdHepP3Mod;
    HMCPion.FourMomentum = StdHepPTLV;
    HMCPion.PDG = pdg;
  }
}

void TransversityVarsLite::HandleStruckNucleon(TLorentzVector &StdHepPTLV,
  Double_t &StdHepP3Mod, Int_t pdg){

  StruckNucleon.Momentum = StdHepP3Mod;
  StruckNucleon.FourMomentum = StdHepPTLV;
  StruckNucleon.PDG = pdg;
}

bool TransversityVarsLite::HandleStdHepParticle(
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
    return true;
  }

  if(StdHepStatus != 1){return false;}
  if(StdHepPdg >= 1000000000){return false;} //Should catch nuclear PDGs
  if(StdHepPdg >= 2000000000){return false;} //GENIE psuedo particles

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
      HandleProton(StdHepPTLV, StdHepP3Mod, StdHepPosition);
      NProtons++;
      NFinalStateParticles++;
      break;
    }
    case 2112:{
      NNeutrons++;
      NFinalStateParticles++;
      break;
    }
    case 211:{
      HandleCPion(StdHepPTLV,StdHepP3Mod,StdHepPdg);
      NPiPlus++;
      NChargedPions++;
      NPions++;
      NFinalStateParticles++;
      break;
    }
    case -211:{
      HandleCPion(StdHepPTLV,StdHepP3Mod,StdHepPdg);
      NPiMinus++;
      NChargedPions++;
      NPions++;
      NFinalStateParticles++;
      break;
    }
    case 111:{
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

void TransversityVarsLite::Finalise(){

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
    FirstProton.Momentum *= GeVToMeV;
    FirstProton.FourMomentum *= GeVToMeV;
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

//Highest Momentum Charged Pion
  HMCPion_PDG = HMCPion.PDG;
  HMCPion_4Mom_MeV = HMCPion.FourMomentum;

//******************************************************************************
//                       'Verse Variable Values
//******************************************************************************

//Projections
  TVector3 HMProtonPt_MeV = GetVectorInTPlane(HMProton.FourMomentum.Vect(),
    MuonNeutrino.FourMomentum.Vect());

//DeltaPhiT
  DeltaPhiT_HMProton_deg = GetDeltaPhiT(
    Muon.FourMomentum.Vect(), HMProton.FourMomentum.Vect(),
    MuonNeutrino.FourMomentum.Vect())*RadToDeg;


//DeltaPt
  DeltaPT_HMProton_MeV = GetDeltaPT(Muon_Pt_MeV,
                                    HMProtonPt_MeV,
                                    MuonNeutrino.FourMomentum.Vect());

  DeltaPTotal_HMProton_MeV = (Muon.FourMomentum + HMProton.FourMomentum
    - MuonNeutrino.FourMomentum - StruckNucleon.FourMomentum);

  DeltaPNuclearEffect_HMProton_MeV = (Muon.FourMomentum + HMProton.FourMomentum
    - MuonNeutrino.FourMomentum);

  DeltaPProton_MeV = (Muon.FourMomentum
    - MuonNeutrino.FourMomentum - StruckNucleon.FourMomentum).Vect();

  DeltaPTargetNucleon_MeV = (Muon.FourMomentum + FirstProton.FourMomentum
    - MuonNeutrino.FourMomentum).Vect();

//DeltaAlphat
  DeltaAlphaT_HMProton_deg = GetDeltaAlphaT(Muon_Pt_MeV,
                                            HMProtonPt_MeV,
                                            MuonNeutrino.FourMomentum.Vect())*RadToDeg;
//DeltaP_TT
  if((HMCPion.Momentum > 1E-3) && (HMProton.Momentum > 1E-3)){

    DeltaP_TT = TransversityUtils::GetDeltaPTT(
      Muon.FourMomentum.Vect(),
      HMCPion.FourMomentum.Vect(),
      HMProton.FourMomentum.Vect(),
      MuonNeutrino.FourMomentum.Vect());

    HMProtonPion_3Mom_MeV = HMProton.FourMomentum.Vect() +
      HMCPion.FourMomentum.Vect();

    TVector3 HMProtonPionPt_MeV = GetVectorInTPlane(HMProtonPion_3Mom_MeV,
      MuonNeutrino.FourMomentum.Vect());

    DeltaPhiT_HMProtonPion_deg = GetDeltaPhiT(
                    Muon.FourMomentum.Vect(),
                    HMProtonPion_3Mom_MeV,
                    MuonNeutrino.FourMomentum.Vect())*RadToDeg;
    DeltaPT_HMProtonPion_MeV = GetDeltaPT(Muon_Pt_MeV,
                    HMProtonPionPt_MeV,
                    MuonNeutrino.FourMomentum.Vect());
    DeltaAlphaT_HMProtonPion_deg = GetDeltaAlphaT(Muon_Pt_MeV,
                    HMProtonPionPt_MeV,
                    MuonNeutrino.FourMomentum.Vect())*RadToDeg;

    DeltaPTotal_HMProtonPion_MeV = (Muon.FourMomentum + HMProton.FourMomentum
      + HMCPion.FourMomentum - MuonNeutrino.FourMomentum
      - StruckNucleon.FourMomentum);
  }

  CCQ2 = (Muon.FourMomentum - MuonNeutrino.FourMomentum).Mag2();

//******************************************************************************
//******************************************************************************
}

void TransversityVarsLite::Reset(){

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

//Highest Momentum Charged Pion
  HMCPion_PDG = 0;
  HMCPion_4Mom_MeV = TLorentzVector(0,0,0,0);


//******************************************************************************
//                       'Verse Variable Values
//******************************************************************************

//DeltaPhiT
  DeltaPhiT_HMProton_deg = 0;

//DeltaPt
  DeltaPT_HMProton_MeV = TVector3(0,0,0);

//DeltaPTotal
  DeltaPTotal_HMProton_MeV = TLorentzVector(0,0,0,0);

//DeltaPNuclearEffect
  DeltaPNuclearEffect_HMProton_MeV = TLorentzVector(0,0,0,0);

//DeltaPProton
  DeltaPProton_MeV = TVector3(0,0,0);

//DeltaPTargetNucleon
  DeltaPTargetNucleon_MeV = TVector3(0,0,0);

//DeltaAlphat
  DeltaAlphaT_HMProton_deg = 0;

  DeltaP_TT = 0;

//ProtonPion Combo Platter
  HMProtonPion_3Mom_MeV = TVector3(0,0,0);
  DeltaPhiT_HMProton_deg = 0;
  DeltaPT_HMProton_MeV = TVector3(0,0,0);
  DeltaAlphaT_HMProton_deg = 0;
  DeltaPTotal_HMProtonPion_MeV = TLorentzVector(0,0,0,0);

//******************************************************************************
//                       Final State Particles
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
//                       Others and Transients
//******************************************************************************

  CCQ2 = 0xdeadbeef;

  //Transients
  Muon.Reset();
  MuonNeutrino.Reset();
  StruckNucleon.Reset();
  HMProton.Reset();
  HMCPion.Reset();
  FirstProton.Reset();

//******************************************************************************
//******************************************************************************

}
