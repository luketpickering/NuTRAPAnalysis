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
}

ClassImp(TransversityVars);

TransversityVars::TransversityVars(
  bool InGeV, Double_t TargetBE_MeV, Int_t NThresh, Int_t* Threshs_MeV,
  Int_t SmearingMag_PhiDeg, Int_t SmearingMag_MuMmtmMeV,
  Int_t SmearingMag_PMmtmMeV, bool ds,
  TString generatorName){
  GeneratorName = generatorName;
  Reset();
  this->NThresh = NThresh;
  this->Threshs_MeV = Threshs_MeV;
  IsInGev = InGeV;
  this->SmearingMag_PhiDeg = SmearingMag_PhiDeg;
  this->SmearingMag_MuMmtmMeV = SmearingMag_MuMmtmMeV;
  this->SmearingMag_PMmtmMeV = SmearingMag_PMmtmMeV;
  this->DoSmear = ds;
  this->TargetBE_MeV = TargetBE_MeV;
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
  for(size_t i = 0; ((i < N) && (Int_t(i) < threshlvl)); ++i){
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
}

void TransversityVars::HandleCPion(TLorentzVector &StdHepPTLV,
  Double_t &StdHepP3Mod, Int_t pdg){

if(StdHepP3Mod > HMCPion.Momentum){
    HMCPion.Momentum = StdHepP3Mod;
    HMCPion.FourMomentum = StdHepPTLV;
    HMCPion.PDG = pdg;
  }
}

void TransversityVars::HandleStruckNucleon(TLorentzVector &StdHepPTLV,
  Double_t &StdHepP3Mod, Int_t pdg){

  StruckNucleon.Momentum = StdHepP3Mod;
  StruckNucleon.FourMomentum = StdHepPTLV;
  StruckNucleon.PDG = pdg;
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
      IncrementThreshArray(NAboveThresholdGammas,
        PartIsAboveThresh(StdHepPTLV*(IsInGev?1000.0:1.0)));
      NFinalStateParticles++;
      break;
    }
    case 2212:{
      HandleProton(StdHepPTLV, StdHepP3Mod, StdHepPosition);
      HandleHMTrackable(StdHepPTLV,StdHepP3Mod,StdHepPdg);
      NProtons++;
      Int_t threshlvl = PartIsAboveThresh(StdHepPTLV*(IsInGev?1000.0:1.0));
      IncrementThreshArray(NAboveThresholdProtons,threshlvl);
      IncrementThreshArray(NAboveThresholdTrackable,threshlvl);
      NFinalStateParticles++;
      break;
    }
    case 2112:{
      Int_t threshlvl = PartIsAboveThresh(StdHepPTLV*(IsInGev?1000.0:1.0));
      IncrementThreshArray(NAboveThresholdNeutrons,threshlvl);
      IncrementThreshArray(NAboveThresholdNeutrals,threshlvl);
      NNeutrons++;
      NFinalStateParticles++;
      break;
    }
    case 211:{
      HandleHMTrackable(StdHepPTLV,StdHepP3Mod,StdHepPdg);
      HandleCPion(StdHepPTLV,StdHepP3Mod,StdHepPdg);
      Int_t threshlvl = PartIsAboveThresh(StdHepPTLV*(IsInGev?1000.0:1.0));
      IncrementThreshArray(NAboveThresholdChargedPions,threshlvl);
      IncrementThreshArray(NAboveThresholdPiPlus,threshlvl);
      IncrementThreshArray(NAboveThresholdTrackable,threshlvl);
      NPiPlus++;
      NChargedPions++;
      NPions++;
      NFinalStateParticles++;
      break;
    }
    case -211:{
      HandleHMTrackable(StdHepPTLV,StdHepP3Mod,StdHepPdg);
      HandleCPion(StdHepPTLV,StdHepP3Mod,StdHepPdg);
      Int_t threshlvl = PartIsAboveThresh(StdHepPTLV*(IsInGev?1000.0:1.0));
      IncrementThreshArray(NAboveThresholdChargedPions,threshlvl);
      IncrementThreshArray(NAboveThresholdPiMinus,threshlvl);
      IncrementThreshArray(NAboveThresholdTrackable,threshlvl);
      NPiMinus++;
      NChargedPions++;
      NPions++;
      NFinalStateParticles++;
      break;
    }
    case 111:{
      Int_t threshlvl = PartIsAboveThresh(StdHepPTLV*(IsInGev?1000.0:1.0));
      IncrementThreshArray(NAboveThresholdPiZero,threshlvl);
      IncrementThreshArray(NAboveThresholdNeutrals,threshlvl);
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
      IncrementThreshArray(NAboveThresholdNeutrals,
        PartIsAboveThresh(StdHepPTLV*(IsInGev?1000.0:1.0)));
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

Double_t GetCauchySmear(TVector3 const & pt, Int_t GaussWidthDeg){
  static bool first = true;
  static TF1 *fa1;
  if(first){
    fa1 = new TF1("CauchyDist","(1/([0])) * ([1]/(x**2 + [1]))",0,1000);
  }
  Float_t gamma = 250.0/pt.Mag();
  Float_t pig = M_PI*gamma;
  fa1->SetParameter(0,pig);
  fa1->SetParameter(1,gamma*gamma);

  return GaussWidthDeg*(fa1->GetRandom());
}

Double_t PullSmear(TVector3 const & pt, Int_t GaussWidthDeg){
  // [1] R. Fr uhwirth et al.: Data Analysis Techniques for High-Energy-Physics,
  // 2nd edition (eds. M. Regler and R. Fr uhwirth), Cambridge University
  //Press (2000).
  if(!GaussWidthDeg){ return 0; }
  static TRandom gaus;
  constexpr static Float_t c1 = 0.5; //Intrinsic smear
  constexpr static Float_t c2 = 0.5/(500.0*500.0); // as a func of dpt[1/MeV]

  Double_t ooptm = 1.0/pt.Mag();

  Double_t DeltaPt_Pt2 = GaussWidthDeg*sqrt(c1 + c2*ooptm*ooptm);

  //500MeV has x degree sigma

  return gaus.Gaus(0.0,DeltaPt_Pt2);
}

Double_t SimpleGSmearMeasure(Double_t const & Meas, Double_t GWidth){

  static TRandom gaus;

  //500MeV has x degree sigma

  return gaus.Gaus(Meas,GWidth);
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
    - MuonNeutrino.FourMomentum - StruckNucleon.FourMomentum).Vect();
  DeltaPTotal_FirstProton_MeV = (Muon.FourMomentum + FirstProton.FourMomentum
    - MuonNeutrino.FourMomentum - StruckNucleon.FourMomentum).Vect();

  DeltaPNuclearEffect_HMProton_MeV = (Muon.FourMomentum + HMProton.FourMomentum
    - MuonNeutrino.FourMomentum).Vect();
  DeltaPNuclearEffect_FirstProton_MeV = (Muon.FourMomentum
    + FirstProton.FourMomentum - MuonNeutrino.FourMomentum).Vect();

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
      - StruckNucleon.FourMomentum).Vect();
  }

//******************************************************************************
//                       Smeared Sample States
//******************************************************************************

  if(DoSmear && (HMTrackablePt_MeV.Mag() > 1E-6)){

    Muon_phi_GSmear_deg = SimpleGSmearMeasure(0,SmearingMag_PhiDeg);
    HMTrackable_phi_GSmear_deg = SimpleGSmearMeasure(0,SmearingMag_PhiDeg);

    TVector3 MuPtSmeared = Muon_Pt_MeV;
    MuPtSmeared.SetPhi(MuPtSmeared.Phi()+
      (Muon_phi_GSmear_deg*DegToRad));

    TVector3 HMTrackableSmeared = HMTrackablePt_MeV;
    HMTrackableSmeared.SetPhi(HMTrackableSmeared.Phi()+
      (HMTrackable_phi_GSmear_deg*DegToRad));

    DeltaPhiT_HMTrackable_GSmear_deg = GetDeltaPhiT(MuPtSmeared,
        HMTrackableSmeared, MuonNeutrino.FourMomentum.Vect())*RadToDeg;

    Muon_MomMag_MeV_GSmear =
      SimpleGSmearMeasure(Muon.Momentum,SmearingMag_MuMmtmMeV);
    HMProton_MomMag_MeV_GSmear =
      SimpleGSmearMeasure(HMProton.Momentum,SmearingMag_PMmtmMeV);
    HMTrackable_MomMag_MeV_GSmear =
      SimpleGSmearMeasure(HMTrackable.Momentum,SmearingMag_PMmtmMeV);
  }

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
  Muon_MomMag_MeV_GSmear = 0;

//First Proton
  FirstProton_PDG = 0;
  FirstProton_4Mom_MeV = TLorentzVector(0,0,0,0);
  FirstProton_StdHepPosition = -1;

//Highest Momentum Proton
  HMProton_PDG = 0;
  HMProton_4Mom_MeV = TLorentzVector(0,0,0,0);
  HMProton_MomMag_MeV_GSmear = 0;
  HMProton_StdHepPosition = -1;

//Highest Momentum Charged Pion
  HMCPion_PDG = 0;
  HMCPion_4Mom_MeV = TLorentzVector(0,0,0,0);

//Highest Momentum Trackable
  HMTrackable_PDG = 0;
  HMTrackable_4Mom_MeV = TLorentzVector(0,0,0,0);
  HMTrackable_MomMag_MeV_GSmear = 0;

//******************************************************************************
//                       'Verse Variable Values
//******************************************************************************

//DeltaPhiT
  DeltaPhiT_HMProton_deg = 0;
  DeltaPhiT_FirstProton_deg = 0;
  DeltaPhiT_HMTrackable_deg = 0;
  DeltaPhiT_HMTrackable_GSmear_deg = 0;

//DeltaPt
  DeltaPT_HMProton_MeV = TVector3(0,0,0);
  DeltaPT_FirstProton_MeV = TVector3(0,0,0);

//DeltaPTotal
  DeltaPTotal_HMProton_MeV = TVector3(0,0,0);
  DeltaPTotal_FirstProton_MeV = TVector3(0,0,0);

//DeltaPNuclearEffect
  DeltaPNuclearEffect_HMProton_MeV = TVector3(0,0,0);
  DeltaPNuclearEffect_FirstProton_MeV = TVector3(0,0,0);

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
  DeltaPTotal_HMProtonPion_MeV = TVector3(0,0,0);

//******************************************************************************
//                       Smeared Sample States
//******************************************************************************

  Muon_phi_GSmear_deg = 0;
  HMTrackable_phi_GSmear_deg = 0;

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
  }

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
