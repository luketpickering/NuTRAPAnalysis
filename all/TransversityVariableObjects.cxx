#include <iostream>

#include "TransversityVariableObjects.hxx"

namespace {
const int kStdHepIdxPx = 0;
const int kStdHepIdxPy = 1;
const int kStdHepIdxPz = 2;
const int kStdHepIdxE  = 3;
}

ClassImp(CCQEFSITransversity);

CCQEFSITransversity::CCQEFSITransversity(){
  Reset();
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

  if(StdHepPosition == 0){ //Neutrino
    MuonNeutrino.PDG = StdHepPdg;
    MuonNeutrino.FourMomentum = StdHepPTLV;
    MuonNeutrino.Momentum = StdHepP3Mod;
    return true;
  } else if(StdHepPosition == 1){ //TargetNucleus
    TargetPDG = StdHepPdg;
    TargetZ = ((StdHepPdg/10000)%1000);
    return true;
  }

  if(StdHepStatus != 1){return false;}
  if(StdHepPdg >= 1000000000){return false;}
  if(StdHepPdg >= 2000000000){return false;}

  switch(StdHepPdg){
    case 13:{
      Muon.PDG = StdHepPdg;
      Muon.FourMomentum = StdHepPTLV;
      Muon.Momentum = StdHepP3Mod;
      break;
    }
    case 14:{
      // std::cerr << "Found another Muon Neutrino!!!" << StdHepPosition
      //   << " " << NeutConventionReactionCode << std::endl;
      return false;
      break;
    }
    case 22:{
      NGamma++;
      break;
    }
    case 2212:{
      HandleProton(StdHepPTLV, StdHepP3Mod);
      NProtons++;
      break;
    }
    case 2112:{
      NNeutrons++;
      break;
    }
    case 211:{
      NPiPlus++;
      NChargedPions++;
      NPions++;
      break;
    }
    case -211:{
      NPiMinus++;
      NChargedPions++;
      NPions++;
      break;
    }
    case 111:{
      NPiZero++;
      NPions++;
      break;
    }
    default:{
      NOtherTraj++;
      return false;
      break;
    }
  }
  return true;

}
void CCQEFSITransversity::Finalise(){

  MuonPDG = Muon.PDG;
  HMProtonPDG = HMProton.PDG;
  FirstProtonPDG = FirstProton.PDG;

  MuonDirection = Muon.FourMomentum.Vect().Unit();

  HMProtonDirection  = HMProton.FourMomentum.Vect().Unit();
  FirstProtonDirection  = FirstProton.FourMomentum.Vect().Unit();

  HMProtonMomentum = HMProton.Momentum;
  FirstProtonMomentum = FirstProton.Momentum;

  IncNeutrinoPDG = MuonNeutrino.PDG;
  IncNeutrinoMmtm = MuonNeutrino.FourMomentum;

  TVector3 PlanarUnitMuon = GetUnitVectorInPlane(
    MuonDirection,
    MuonNeutrino.FourMomentum.Vect());

  TVector3 PlanarUnitHMProton = GetUnitVectorInPlane(
    HMProtonDirection,
    MuonNeutrino.FourMomentum.Vect());

  TVector3 PlanarUnitFirstProton = GetUnitVectorInPlane(
    FirstProtonDirection,
    MuonNeutrino.FourMomentum.Vect());

  DeltaPhiT_HMProton =
    (M_PI - acos(PlanarUnitMuon.Dot(PlanarUnitHMProton)));
  DeltaPhiT_FirstProton =
    (M_PI - acos(PlanarUnitMuon.Dot(PlanarUnitFirstProton)));

  MuonPt = GetVectorInPlane(Muon.FourMomentum.Vect(),
    MuonNeutrino.FourMomentum.Vect());
}
void CCQEFSITransversity::Reset(){

  DeltaPhiT_HMProton = 0;
  DeltaPhiT_FirstProton = 0;
  NeutConventionReactionCode = 0;
  NProtons = 0;
  NGamma = 0;
  NNeutrons = 0;
  NPiPlus = 0;
  NPiZero = 0;
  NPiMinus = 0;
  NPions = 0;
  NChargedPions = 0;
  NOtherTraj = 0;
  MuonPDG = 0;
  HMProtonPDG = 0;
  FirstProtonPDG = 0;
  MuonDirection = TVector3(0,0,0);
  HMProtonDirection = TVector3(0,0,0);
  FirstProtonDirection = TVector3(0,0,0);
  HMProtonMomentum = 0;
  FirstProtonMomentum = 0;
  MuonPt = TVector3(0,0,0);

  IncNeutrinoPDG = 0;
  IncNeutrinoMmtm = TLorentzVector(0,0,0,0);

  TargetPDG = 0;
  TargetZ = 0;
}
