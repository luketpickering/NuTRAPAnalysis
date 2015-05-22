#include <iostream>

#include "TransversityVariableObjects.hxx"

namespace {
constexpr int kStdHepIdxPx = 0;
constexpr int kStdHepIdxPy = 1;
constexpr int kStdHepIdxPz = 2;
constexpr int kStdHepIdxE = 3;
}

ClassImp(MuonProtonTransversity);
ClassImp(PionProductionTransversity);

MuonProtonTransversity::MuonProtonTransversity(
  TransversityUtils::Generators gen){
  Reset();
  Gen = gen;
  if((Gen==TransversityUtils::kGENIE) || (Gen==TransversityUtils::kNuWro)){
    IsInGev = true;
  } else { IsInGev = false; }
}

void MuonProtonTransversity::HandleProton(TLorentzVector &StdHepPTLV,
  Double_t &StdHepP3Mod){

  if(!FirstProton.PDG){
    FirstProton.PDG = 2212;
    FirstProton.Momentum = StdHepP3Mod;
    FirstProton.FourMomentum = StdHepPTLV;

    HMProton.PDG = 2212;
    HMProton.Momentum = StdHepP3Mod;
    HMProton.FourMomentum = StdHepPTLV;
  }
  else if(StdHepP3Mod > HMProton.Momentum){
    HMProton.Momentum = StdHepP3Mod;
    HMProton.FourMomentum = StdHepPTLV;
  }
}

bool MuonProtonTransversity::HandleStdHepParticle(
  UInt_t &StdHepPosition,
  int &StdHepPdg,
  int &StdHepStatus,
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
    case 13:{
      Muon.PDG = StdHepPdg;
      Muon.FourMomentum = StdHepPTLV;
      Muon.Momentum = StdHepP3Mod;
      NFinalStateParticles++;
      break;
    }
    //NC
    case 14:{
      return false;
      break;
    }
    case 22:{
      NGamma++;
      NFinalStateParticles++;
      break;
    }
    case 2212:{
      HandleProton(StdHepPTLV, StdHepP3Mod);
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
      NPiPlus++;
      NChargedPions++;
      NPions++;
      NFinalStateParticles++;
      break;
    }
    case -211:{
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

Double_t GetReconNuEnergy(const TVector3 &nudir,
  const TLorentzVector &mumom, const TLorentzVector &pmom){

  TVector3 m3mom = mumom.Vect();
  TVector3 p3mom = pmom.Vect();

  Double_t p_tot_p = p3mom.Mag();
  Double_t p_tot_m = m3mom.Mag();
  Double_t nudirMag = nudir.Mag();

  Double_t mnangle =
    TMath::ACos((nudir.Dot(m3mom))/(p_tot_m*nudirMag));

  Double_t EnuRec =
    p_tot_m*TMath::Cos(mnangle)
      + sqrt(p_tot_m*TMath::Cos(mnangle)*p_tot_m*TMath::Cos(mnangle)
        - p_tot_m*p_tot_m + p_tot_p*p_tot_p);
  return EnuRec;
}

Double_t GetReconTgtMass(Double_t const &nuE,
  const TLorentzVector &mumom, const TLorentzVector &pmom){
  Double_t TgtMass = pmom.E() + mumom.E() - nuE;

  return TgtMass;
}

void MuonProtonTransversity::Finalise(){

  constexpr Float_t GeVToMeV = 1000.0;
  constexpr Float_t RadToDeg = 180.0/M_PI;

  MuonPDG = Muon.PDG;
  HMProtonPDG = HMProton.PDG;
  FirstProtonPDG = FirstProton.PDG;

  if(IsInGev){
    Muon.Momentum *= GeVToMeV;
    Muon.FourMomentum *= GeVToMeV;
    MuonNeutrino.Momentum *= GeVToMeV;
    MuonNeutrino.FourMomentum *= GeVToMeV;
    HMProton.Momentum *= GeVToMeV;
    HMProton.FourMomentum *= GeVToMeV;
    FirstProton.Momentum *= GeVToMeV;
    FirstProton.FourMomentum *= GeVToMeV;
  }

  MuonDirection = Muon.FourMomentum.Vect().Unit();

  HMProtonDirection  = HMProton.FourMomentum.Vect().Unit();
  FirstProtonDirection  = FirstProton.FourMomentum.Vect().Unit();

  MuonMomentum_MeV = Muon.Momentum;
  HMProtonMomentum_MeV = HMProton.Momentum;
  FirstProtonMomentum_MeV = FirstProton.Momentum;

  IncNeutrinoPDG = MuonNeutrino.PDG;
  IncNeutrinoMmtm = MuonNeutrino.FourMomentum;

  DeltaPhiT_HMProton = TransversityUtils::GetDeltaPhiT(MuonDirection,
    HMProtonDirection, IncNeutrinoMmtm.Vect());
  DeltaPhiT_FirstProton = TransversityUtils::GetDeltaPhiT(MuonDirection,
    FirstProtonDirection, IncNeutrinoMmtm.Vect());

  DeltaPhiT_HMProton_deg = DeltaPhiT_HMProton * RadToDeg;
  DeltaPhiT_FirstProton_deg = DeltaPhiT_FirstProton * RadToDeg;

  MuonPt_MeV = TransversityUtils::GetVectorInTPlane(Muon.FourMomentum.Vect(),
    IncNeutrinoMmtm.Vect());

  HMProtonPt_MeV =
    TransversityUtils::GetVectorInTPlane(HMProton.FourMomentum.Vect(),
      IncNeutrinoMmtm.Vect());

  FirstProtonPt_MeV = TransversityUtils::GetVectorInTPlane(
    FirstProton.FourMomentum.Vect(), IncNeutrinoMmtm.Vect());

  DeltaPT_HMProton_MeV = TransversityUtils::GetDeltaPT(MuonPt_MeV,
                                                        HMProtonPt_MeV,
                                                        IncNeutrinoMmtm.Vect());
  DeltaPT_FirstProton_MeV = TransversityUtils::GetDeltaPT(MuonPt_MeV,
                                                        FirstProtonPt_MeV,
                                                        IncNeutrinoMmtm.Vect());

  DeltaAlphaT_HMProton = TransversityUtils::GetDeltaAlphaT(MuonPt_MeV,
                                                        HMProtonPt_MeV,
                                                        IncNeutrinoMmtm.Vect());
  DeltaAlphaT_FirstProton = TransversityUtils::GetDeltaAlphaT(MuonPt_MeV,
                                                        FirstProtonPt_MeV,
                                                        IncNeutrinoMmtm.Vect());

  DeltaAlphaT_HMProton_deg = DeltaAlphaT_HMProton*RadToDeg;
  DeltaAlphaT_FirstProton_deg = DeltaAlphaT_FirstProton*RadToDeg;

  ReconNuEnergy = GetReconNuEnergy(MuonNeutrino.FourMomentum.Vect(),
                                    Muon.FourMomentum,
                                    HMProton.FourMomentum);
  ReconTargetMass = GetReconTgtMass(ReconNuEnergy,
                                    Muon.FourMomentum,
                                    HMProton.FourMomentum);

}
void MuonProtonTransversity::Reset(){

  ReconNuEnergy = 0;
  ReconTargetMass = 0;

  DeltaPhiT_HMProton = 0;
  DeltaPhiT_FirstProton = 0;

  DeltaPhiT_HMProton_deg = 0;
  DeltaPhiT_FirstProton_deg = 0;

  DeltaPT_HMProton_MeV = TVector3(0,0,0);
  DeltaPT_FirstProton_MeV = TVector3(0,0,0);

  DeltaAlphaT_HMProton = 0;
  DeltaAlphaT_FirstProton = 0;

  DeltaAlphaT_HMProton_deg = 0;
  DeltaAlphaT_FirstProton_deg = 0;

  NeutConventionReactionCode = 0;

  NFinalStateParticles = 0;
  NProtons = 0;
  NGamma = 0;
  NNeutrons = 0;
  NPiPlus = 0;
  NPiZero = 0;
  NPiMinus = 0;
  NPions = 0;
  NChargedPions = 0;
  NOtherParticles = 0;

  MuonPDG = 0;
  HMProtonPDG = 0;
  FirstProtonPDG = 0;

  MuonDirection = TVector3(0,0,0);
  MuonMomentum_MeV = 0;
  HMProtonDirection = TVector3(0,0,0);
  FirstProtonDirection = TVector3(0,0,0);
  HMProtonMomentum_MeV = 0;
  FirstProtonMomentum_MeV = 0;

  MuonPt_MeV = TVector3(0,0,0);
  HMProtonPt_MeV = TVector3(0,0,0);
  FirstProtonPt_MeV = TVector3(0,0,0);

  IncNeutrinoPDG = 0;
  IncNeutrinoMmtm = TLorentzVector(0,0,0,0);

  TargetPDG = 0;
  TargetZ = 0;

  Muon.Reset();
  MuonNeutrino.Reset();
  HMProton.Reset();
  FirstProton.Reset();
}

PionProductionTransversity::PionProductionTransversity(
  TransversityUtils::Generators gen){
  Reset();
  Gen = gen;
  if((Gen==TransversityUtils::kGENIE) || (Gen==TransversityUtils::kNuWro)){
    IsInGev = true;
  } else { IsInGev = false; }
}

void PionProductionTransversity::HandlePiPlus(TLorentzVector &StdHepPTLV,
  Double_t &StdHepP3Mod){

  if(!FirstPiPlus.PDG){
    FirstPiPlus.PDG = 211;
    FirstPiPlus.Momentum = StdHepP3Mod;
    FirstPiPlus.FourMomentum = StdHepPTLV;

    HMPiPlus.PDG = 211;
    HMPiPlus.Momentum = StdHepP3Mod;
    HMPiPlus.FourMomentum = StdHepPTLV;
  }
  else if(StdHepP3Mod > HMPiPlus.Momentum){
    HMPiPlus.Momentum = StdHepP3Mod;
    HMPiPlus.FourMomentum = StdHepPTLV;
  }
}

bool PionProductionTransversity::HandleStdHepParticle(
  UInt_t &StdHepPosition,
  int &StdHepPdg,
  int &StdHepStatus,
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
    case 13:{
      Muon.PDG = StdHepPdg;
      Muon.FourMomentum = StdHepPTLV;
      Muon.Momentum = StdHepP3Mod;
      NFinalStateParticles++;
      break;
    }
    //NC
    case 14:{
      return false;
      break;
    }
    case 22:{
      NGamma++;
      NFinalStateParticles++;
      break;
    }
    case 2212:{
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
      HandlePiPlus(StdHepPTLV, StdHepP3Mod);
      NPiPlus++;
      NChargedPions++;
      NPions++;
      NFinalStateParticles++;
      break;
    }
    case -211:{
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
void PionProductionTransversity::Finalise(){

  constexpr Float_t GeVToMeV = 1000.0;
  constexpr Float_t RadToDeg = 180.0/M_PI;

  MuonPDG = Muon.PDG;
  HMPiPlusPDG = HMPiPlus.PDG;
  FirstPiPlusPDG = FirstPiPlus.PDG;

  if(IsInGev){
    Muon.Momentum *= GeVToMeV;
    Muon.FourMomentum *= GeVToMeV;
    MuonNeutrino.Momentum *= GeVToMeV;
    MuonNeutrino.FourMomentum *= GeVToMeV;
    HMPiPlus.Momentum *= GeVToMeV;
    HMPiPlus.FourMomentum *= GeVToMeV;
    FirstPiPlus.Momentum *= GeVToMeV;
    FirstPiPlus.FourMomentum *= GeVToMeV;
  }

  MuonDirection = Muon.FourMomentum.Vect().Unit();
  HMPiPlusDirection = HMPiPlus.FourMomentum.Vect().Unit();
  FirstPiPlusDirection = FirstPiPlus.FourMomentum.Vect().Unit();

  MuonMomentum_MeV = Muon.Momentum;
  HMPiPlusMomentum_MeV = HMPiPlus.Momentum;
  FirstPiPlusMomentum_MeV = FirstPiPlus.Momentum;

  IncNeutrinoPDG = MuonNeutrino.PDG;
  IncNeutrinoMmtm = MuonNeutrino.FourMomentum;

  DeltaPhiT_HMPiPlus = TransversityUtils::GetDeltaPhiT(MuonDirection,
    HMPiPlusDirection, IncNeutrinoMmtm.Vect());
  DeltaPhiT_FirstPiPlus = TransversityUtils::GetDeltaPhiT(MuonDirection,
    FirstPiPlusDirection, IncNeutrinoMmtm.Vect());

  DeltaPhiT_HMPiPlus_deg = DeltaPhiT_HMPiPlus * RadToDeg;
  DeltaPhiT_FirstPiPlus_deg = DeltaPhiT_FirstPiPlus * RadToDeg;

  MuonPt_MeV = TransversityUtils::GetVectorInTPlane(Muon.FourMomentum.Vect(),
    MuonNeutrino.FourMomentum.Vect());

  FirstPiPlusPt_MeV =
    TransversityUtils::GetVectorInTPlane(FirstPiPlus.FourMomentum.Vect(),
      MuonNeutrino.FourMomentum.Vect());

  HMPiPlusPt_MeV =
    TransversityUtils::GetVectorInTPlane(HMPiPlus.FourMomentum.Vect(),
      MuonNeutrino.FourMomentum.Vect());

  DeltaPT_HMPiPlus_MeV = TransversityUtils::GetDeltaPT(MuonPt_MeV,
                                                        HMPiPlusPt_MeV,
                                                        IncNeutrinoMmtm.Vect());
  DeltaPT_FirstPiPlus_MeV = TransversityUtils::GetDeltaPT(MuonPt_MeV,
                                                        FirstPiPlusPt_MeV,
                                                        IncNeutrinoMmtm.Vect());

  DeltaAlphaT_HMPiPlus = TransversityUtils::GetDeltaAlphaT(MuonPt_MeV,
                                                        HMPiPlusPt_MeV,
                                                        IncNeutrinoMmtm.Vect());
  DeltaAlphaT_FirstPiPlus = TransversityUtils::GetDeltaAlphaT(MuonPt_MeV,
                                                        FirstPiPlusPt_MeV,
                                                        IncNeutrinoMmtm.Vect());

  DeltaAlphaT_HMPiPlus_deg = DeltaAlphaT_HMPiPlus*RadToDeg;
  DeltaAlphaT_FirstPiPlus_deg = DeltaAlphaT_FirstPiPlus*RadToDeg;

}
void PionProductionTransversity::Reset(){

  DeltaPhiT_HMPiPlus = 0;
  DeltaPhiT_FirstPiPlus = 0;

  DeltaPhiT_HMPiPlus_deg = 0;
  DeltaPhiT_FirstPiPlus_deg = 0;

  DeltaPT_HMPiPlus_MeV = TVector3(0,0,0);
  DeltaPT_FirstPiPlus_MeV = TVector3(0,0,0);

  DeltaAlphaT_HMPiPlus = 0;
  DeltaAlphaT_FirstPiPlus = 0;

  DeltaAlphaT_HMPiPlus_deg = 0;
  DeltaAlphaT_FirstPiPlus_deg = 0;

  NeutConventionReactionCode = 0;

  NFinalStateParticles = 0;
  NProtons = 0;
  NGamma = 0;
  NNeutrons = 0;
  NPiPlus = 0;
  NPiZero = 0;
  NPiMinus = 0;
  NPions = 0;
  NChargedPions = 0;
  NOtherParticles = 0;

  MuonPDG = 0;
  HMPiPlusPDG = 0;
  FirstPiPlusPDG = 0;

  MuonDirection = TVector3(0,0,0);
  HMPiPlusDirection = TVector3(0,0,0);
  FirstPiPlusDirection = TVector3(0,0,0);

  MuonMomentum_MeV = 0;
  HMPiPlusMomentum_MeV = 0;
  FirstPiPlusMomentum_MeV = 0;

  MuonPt_MeV = TVector3(0,0,0);

  FirstPiPlusPt_MeV = TVector3(0,0,0);
  HMPiPlusPt_MeV = TVector3(0,0,0);

  IncNeutrinoPDG = 0;
  IncNeutrinoMmtm = TLorentzVector(0,0,0,0);

  TargetPDG = 0;
  TargetZ = 0;

  Muon.Reset();
  MuonNeutrino.Reset();
  HMPiPlus.Reset();
  FirstPiPlus.Reset();
}
