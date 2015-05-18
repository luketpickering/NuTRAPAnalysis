#include <iostream>

#include "TransversityVariableObjects.hxx"

namespace {
constexpr int kStdHepIdxPx = 0;
constexpr int kStdHepIdxPy = 1;
constexpr int kStdHepIdxPz = 2;
constexpr int kStdHepIdxE = 3;
}

ClassImp(CCQEFSITransversity);
ClassImp(PionProductionTransversity);

CCQEFSITransversity::CCQEFSITransversity(TransversityUtils::Generators gen){
  Reset();
  Gen = gen;
}

void CCQEFSITransversity::HandleProton(TLorentzVector &StdHepPTLV,
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

bool CCQEFSITransversity::HandleStdHepParticle(
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
void CCQEFSITransversity::Finalise(){

  constexpr Float_t GeVToMeV = 1.0/1000.0;
  constexpr Float_t RadToDeg = 180.0/M_PI;

  MuonPDG = Muon.PDG;
  HMProtonPDG = HMProton.PDG;
  FirstProtonPDG = FirstProton.PDG;

  MuonDirection = Muon.FourMomentum.Vect().Unit();

  HMProtonDirection  = HMProton.FourMomentum.Vect().Unit();
  FirstProtonDirection  = FirstProton.FourMomentum.Vect().Unit();

  MuonMomentum_MeV = Muon.Momentum;
  HMProtonMomentum_MeV = HMProton.Momentum;
  FirstProtonMomentum_MeV = FirstProton.Momentum;

  if(Gen == TransversityUtils::kGENIE){
    MuonMomentum_MeV *= GeVToMeV;
    HMProtonMomentum_MeV *= GeVToMeV;
    FirstProtonMomentum_MeV *= GeVToMeV;
  }

  IncNeutrinoPDG = MuonNeutrino.PDG;
  IncNeutrinoMmtm = MuonNeutrino.FourMomentum;

  DeltaPhiT_HMProton = TransversityUtils::GetDeltaPhiT(MuonDirection,
    HMProtonDirection, IncNeutrinoMmtm.Vect());
  DeltaPhiT_FirstProton = TransversityUtils::GetDeltaPhiT(MuonDirection,
    FirstProtonDirection, IncNeutrinoMmtm.Vect());

  DeltaPhiT_HMProton_deg = DeltaPhiT_HMProton * RadToDeg;
  DeltaPhiT_FirstProton_deg = DeltaPhiT_FirstProton * RadToDeg;

  MuonPt = TransversityUtils::GetVectorInTPlane(Muon.FourMomentum.Vect(),
    IncNeutrinoMmtm.Vect());

  HMProtonPt = TransversityUtils::GetVectorInTPlane(HMProton.FourMomentum.Vect(),
    IncNeutrinoMmtm.Vect());

  FirstProtonPt = TransversityUtils::GetVectorInTPlane(
    FirstProton.FourMomentum.Vect(), IncNeutrinoMmtm.Vect());

  if(Gen == TransversityUtils::kGENIE){
    MuonPt *= GeVToMeV;
    HMProtonPt *= GeVToMeV;
    FirstProtonPt *= GeVToMeV;
  }

}
void CCQEFSITransversity::Reset(){

  DeltaPhiT_HMProton = 0;
  DeltaPhiT_FirstProton = 0;

  DeltaPhiT_HMProton_deg = 0;
  DeltaPhiT_FirstProton_deg = 0;

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

  MuonPt = TVector3(0,0,0);
  HMProtonPt = TVector3(0,0,0);
  FirstProtonPt = TVector3(0,0,0);

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

  constexpr Float_t GeVToMeV = 1.0/1000.0;
  constexpr Float_t RadToDeg = 180.0/M_PI;

  MuonPDG = Muon.PDG;
  HMPiPlusPDG = HMPiPlus.PDG;
  FirstPiPlusPDG = FirstPiPlus.PDG;

  MuonDirection = Muon.FourMomentum.Vect().Unit();
  HMPiPlusDirection = HMPiPlus.FourMomentum.Vect().Unit();
  FirstPiPlusDirection = FirstPiPlus.FourMomentum.Vect().Unit();

  MuonMomentum_MeV = Muon.Momentum;
  HMPiPlusMomentum_MeV = HMPiPlus.Momentum;
  FirstPiPlusMomentum_MeV = FirstPiPlus.Momentum;

  if(Gen == TransversityUtils::kGENIE){
    MuonMomentum_MeV *= GeVToMeV;
    HMPiPlusMomentum_MeV *= GeVToMeV;
    FirstPiPlusMomentum_MeV *= GeVToMeV;
  }

  IncNeutrinoPDG = MuonNeutrino.PDG;
  IncNeutrinoMmtm = MuonNeutrino.FourMomentum;

  DeltaPhiT_HMPiPlus = TransversityUtils::GetDeltaPhiT(MuonDirection,
    HMPiPlusDirection, IncNeutrinoMmtm.Vect());
  DeltaPhiT_FirstPiPlus = TransversityUtils::GetDeltaPhiT(MuonDirection,
    FirstPiPlusDirection, IncNeutrinoMmtm.Vect());

  DeltaPhiT_HMPiPlus_deg = DeltaPhiT_HMPiPlus * RadToDeg;
  DeltaPhiT_FirstPiPlus_deg = DeltaPhiT_FirstPiPlus * RadToDeg;

  MuonPt = TransversityUtils::GetVectorInTPlane(Muon.FourMomentum.Vect(),
    MuonNeutrino.FourMomentum.Vect());

  FirstPiPlusPt = TransversityUtils::GetVectorInTPlane(FirstPiPlus.FourMomentum.Vect(),
    MuonNeutrino.FourMomentum.Vect());

  HMPiPlusPt = TransversityUtils::GetVectorInTPlane(HMPiPlus.FourMomentum.Vect(),
    MuonNeutrino.FourMomentum.Vect());

  if(Gen == TransversityUtils::kGENIE){
    MuonPt *= GeVToMeV;
    FirstPiPlusPt *= GeVToMeV;
    HMPiPlusPt *= GeVToMeV;
  }
}
void PionProductionTransversity::Reset(){

  DeltaPhiT_HMPiPlus = 0;
  DeltaPhiT_FirstPiPlus = 0;
  DeltaPhiT_HMPiPlus_deg = 0;
  DeltaPhiT_FirstPiPlus_deg = 0;

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

  MuonPt = TVector3(0,0,0);

  FirstPiPlusPt = TVector3(0,0,0);
  HMPiPlusPt = TVector3(0,0,0);

  IncNeutrinoPDG = 0;
  IncNeutrinoMmtm = TLorentzVector(0,0,0,0);

  TargetPDG = 0;
  TargetZ = 0;

  Muon.Reset();
  MuonNeutrino.Reset();
  HMPiPlus.Reset();
  FirstPiPlus.Reset();

}
