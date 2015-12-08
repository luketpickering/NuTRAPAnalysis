#ifndef __GENERATOR_SPECIFICS_SEEN__
#define __GENERATOR_SPECIFICS_SEEN__

#include "TTree.h"

#include "TransversityVariableObjects.hxx"
#include "TransversityUtils.hxx"

enum GeneratorID {kNEUT,kGENIE,kNuWro,kEmuNuWro,kGiBUU,kInvalid};

constexpr int kStdHepIdxPx = 0;
constexpr int kStdHepIdxPy = 1;
constexpr int kStdHepIdxPz = 2;
constexpr int kStdHepIdxE = 3;

class Generator {
protected:
  TString GeneratorName;

  TransversityVarsB* OutObjectInfo;

  Int_t* StdHepN = 0;
  Int_t* StdHepPdg = 0;
  Double_t** StdHepP4 = 0;
  Int_t* StdHepStatus = 0;

  bool LiteOutput;

  virtual void StartEvent() = 0;

  virtual void HandleStdHepParticle(
    UInt_t &StdHepPosition,
    Int_t &StdHepPdg,
    Int_t &StdHepStatus,
    Double_t * &StdHepP4) = 0;

  void HandleStruckNucleon(Double_t *&StdHepP4, Int_t pdg);
public:
  int NeutConventionReactionCode;

  virtual void Finalise();

  GeneratorID generator;
  std::string TreeName;
  Generator();
  virtual ~Generator();

  virtual void Init(TTree* tree) = 0;

  virtual void AddOutputBranches(TTree* tree, bool LiteOutput,
    bool MultiplyByGeVToMeV=true, Int_t NThresh=0, Int_t* Threshs_MeV=0);

  void DoEvent();
};

class NEUT : public Generator {
protected:
  constexpr static int kNStdHepNPmax = 100;

  TObjString* NeutReacCode = 0;
  Int_t NStdHepN;
  Int_t NStdHepPdg[kNStdHepNPmax];
  Int_t NStdHepStatus[kNStdHepNPmax];
  Double_t NStdHepP4[kNStdHepNPmax][4];

  void StartEvent();

  void HandleStdHepParticle(
    UInt_t &StdHepPosition,
    Int_t &StdHepPdg,
    Int_t &StdHepStatus,
    Double_t * &StdHepP4);

public:
  NEUT() : Generator() { generator = kNEUT; TreeName = "nRooTracker";
  GeneratorName = "NEUT"; }

  void Init(TTree* tree);

};

class GENIE : public Generator {
  constexpr static int kGStdHepNPmax = 350;

  Bool_t ProtonRescat_contains_NoInt;
  Bool_t ProtonRescat_contains_chrgEx;
  Bool_t ProtonRescat_contains_elastic;
  Bool_t ProtonRescat_contains_inelastic;
  Bool_t ProtonRescat_contains_knockout;

  Int_t G2NeutEvtCode;
  Int_t GStdHepN;
  Int_t GStdHepPdg[kGStdHepNPmax];
  Int_t GStdHepStatus[kGStdHepNPmax];
  Int_t GStdHepRescat[kGStdHepNPmax];
  Double_t GStdHepP4[kGStdHepNPmax][4];

  void StartEvent();

  void HandleRescat(Int_t PDG, Int_t RescatCode);

  void HandleStdHepParticle(
    UInt_t &StdHepPosition,
    Int_t &StdHepPdg,
    Int_t &StdHepStatus,
    Double_t * &StdHepP4);

  public:
    GENIE() : Generator() { generator = kGENIE; TreeName = "gRooTracker";
  GeneratorName = "GENIE"; }

  void Init(TTree* tree);

  void AddOutputBranches(TTree* tree, bool LiteOutput,
    bool MultiplyByGeVToMeV=true, Int_t NThresh=0, Int_t* Threshs_MeV=0);
};

class NuWro : public Generator {
protected:
  constexpr static int kNuStdHepNPmax = 4000;

  TObjString* NuWroEvtCode = 0;
  Int_t NuStdHepN;
  Int_t NuStdHepPdg[kNuStdHepNPmax];
  Int_t NuStdHepStatus[kNuStdHepNPmax];
  Double_t NuStdHepP4[kNuStdHepNPmax][4];

  virtual void StartEvent();

  virtual void HandleStdHepParticle(
    UInt_t &StdHepPosition,
    Int_t &StdHepPdg,
    Int_t &StdHepStatus,
    Double_t * &StdHepP4);

  public:
    NuWro() : Generator() { generator = kNuWro; TreeName = "nRooTracker";
  GeneratorName = "NuWro"; }

  virtual void Init(TTree* tree);
  virtual void Finalise();
};

class EmuNuWro : public NuWro {
protected:
  Int_t StruckNucleonPDG;
  void HandleStdHepParticle(
    UInt_t &StdHepPosition,
    Int_t &StdHepPdg,
    Int_t &StdHepStatus,
    Double_t * &StdHepP4);

public:
  EmuNuWro() : NuWro() { generator = kEmuNuWro; };
  void Init(TTree* tree);

};

class GiBUU : public Generator {
  constexpr static int kGiStdHepNPmax = 100;

  Int_t Gi2NeutEvtCode;
  Int_t GiBUUReactionCode;
  Double_t GiBUUPerWeight;
  Int_t GiStdHepN;
  Int_t GiStdHepPdg[kGiStdHepNPmax];
  Int_t GiStdHepStatus[kGiStdHepNPmax];
  Double_t GiStdHepP4[kGiStdHepNPmax][4];

  void StartEvent();

  void HandleStdHepParticle(
    UInt_t &StdHepPosition,
    Int_t &StdHepPdg,
    Int_t &StdHepStatus,
    Double_t * &StdHepP4);

  public:
    GiBUU() : Generator() { generator = kGiBUU; TreeName = "giRooTracker";
  GeneratorName = "GiBUU"; }

  void Init(TTree* tree);

  void AddOutputBranches(TTree* tree, bool LiteOutput,
    bool MultiplyByGeVToMeV=true, Int_t NThresh=0, Int_t* Threshs_MeV=0);

};

#endif
