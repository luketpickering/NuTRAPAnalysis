#ifndef __TRANSVERSITY_VARIABLE_OBJECTS_SEEN__
#define __TRANSVERSITY_VARIABLE_OBJECTS_SEEN__
#include "TObject.h"
#include "TLorentzVector.h"

namespace {
  TVector3 GetVectorInPlane(const TVector3& inp,
  const TVector3& planarNormal){
  TVector3 pnUnit = planarNormal.Unit();
  double inpProjectPN = inp.Dot(pnUnit);

  TVector3 InPlanarInput = inp-(inpProjectPN*pnUnit);
  return InPlanarInput;
}

TVector3 GetUnitVectorInPlane(const TVector3& inp,
  const TVector3& planarNormal){
  return GetVectorInPlane(inp,planarNormal).Unit();
}
}

struct PartStruct {
  PartStruct(){
    PDG = 0;
    Momentum = 0;
    FourMomentum = TLorentzVector(0,0,0,0);
  }
  Int_t PDG;
  Double_t Momentum;
  TLorentzVector FourMomentum;
};

struct CCQEFSITransversity : public TObject {
private:
  //Transients
  PartStruct Muon; //!
  PartStruct MuonNeutrino; //!
  PartStruct HMProton; //!
  PartStruct FirstProton; //!

public:
  CCQEFSITransversity();

  Double_t DeltaPhiT_HMProton;
  Double_t DeltaPhiT_FirstProton;

  Int_t NeutConventionReactionCode;

  Int_t NProtons;
  Int_t NGamma;
  Int_t NNeutrons;

  Int_t NPiPlus;
  Int_t NPiZero;
  Int_t NPiMinus;
  Int_t NPions;
  Int_t NChargedPions;

  Int_t NOtherTraj;

  Int_t MuonPDG;
  Int_t HMProtonPDG;
  Int_t FirstProtonPDG;

  TVector3 MuonDirection;
  TVector3 HMProtonDirection;
  TVector3 FirstProtonDirection;
  Double_t HMProtonMomentum;
  Double_t FirstProtonMomentum;

  TVector3 MuonPt;

  Int_t IncNeutrinoPDG;
  TLorentzVector IncNeutrinoMmtm;

  Int_t TargetPDG;
  Int_t TargetZ;

  bool SetNeutConventionReactionCode(int rc){
    NeutConventionReactionCode = rc;
    return (rc==1);
  }

  void HandleProton(TLorentzVector &StdHepPTLV,
    Double_t &StdHepP3Mod);

  bool HandleStdHepParticle(UInt_t &StdHepPosition,
                            int &StdHepPdg,
                            int &StdHepStatus,
                            Double_t * &StdHepP4);
  void Finalise();
  void Reset();

  ClassDef(CCQEFSITransversity,1);
};

// struct PionProductionTransversity : public TObject {


// };

#endif
