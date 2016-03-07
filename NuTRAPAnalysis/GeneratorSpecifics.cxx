#include <iostream>

#include "PureGenUtils.hxx"

#include "GeneratorSpecifics.hxx"

//******************************************************************************
//                     Generic
//******************************************************************************

Generator::Generator(){ generator = kInvalid; TreeName = ""; }
Generator::~Generator(){
  if(OutObjectInfo!=nullptr){
    delete OutObjectInfo;
  }
  if(StdHepP4!=nullptr){
    delete [] StdHepP4;
  }
}

void Generator::AddOutputBranches(TTree* tree, bool LiteOutput,
  bool MultiplyByGeVToMeV, Int_t NThresh, Int_t* Threshs_MeV){

  OutObjectInfo = (LiteOutput?
    new TransversityVarsB(MultiplyByGeVToMeV):
    new TransversityVars(MultiplyByGeVToMeV,
      NThresh, Threshs_MeV, GeneratorName));
  OutObjectInfo->AddBranches(tree);
  this->LiteOutput = LiteOutput;
}

void Generator::DoEvent(){
  OutObjectInfo->Reset();
  this->StartEvent();
  for(UInt_t partNum = 0; partNum < UInt_t(*StdHepN); ++partNum){
    this->HandleStdHepParticle(partNum, StdHepPdg[partNum],
      StdHepStatus[partNum], StdHepP4[partNum]);
  }
  this->Finalise();
}

void Generator::HandleStruckNucleon(Double_t *&StdHepP4, Int_t pdg){
  TLorentzVector StdHepPTLV = TLorentzVector(
    StdHepP4[kStdHepIdxPx],
    StdHepP4[kStdHepIdxPy],
    StdHepP4[kStdHepIdxPz],
    StdHepP4[kStdHepIdxE]);
  OutObjectInfo->HandleStruckNucleon(StdHepPTLV,pdg);
}

void Generator::Finalise(){
  OutObjectInfo->Finalise();
  //Fixes broken antineutrino codes.
  if( (OutObjectInfo->IncNeutrino_PDG < 0) &&
             (NeutConventionReactionCode > 0) ) {
    //want antinu codes to be < 0
    OutObjectInfo->NeutConventionReactionCode =
      (-1*NeutConventionReactionCode);
  } else {
    OutObjectInfo->NeutConventionReactionCode = NeutConventionReactionCode;
  }
}

//******************************************************************************
//                     NEUT
//******************************************************************************

void NEUT::Init(TChain* tree){
  StdHepN = &NStdHepN;
  StdHepPdg = NStdHepPdg;
  StdHepP4 = PGUtils::NewPPOf2DArray(NStdHepP4);
  StdHepStatus = NStdHepStatus;

  tree->SetBranchAddress("EvtCode", &NeutReacCode);
  tree->SetBranchAddress("StdHepN", &NStdHepN);
  tree->SetBranchAddress("StdHepPdg", NStdHepPdg);
  tree->SetBranchAddress("StdHepP4", NStdHepP4);
  tree->SetBranchAddress("StdHepStatus", NStdHepStatus);
}

void NEUT::StartEvent(){
  if(PGUtils::str2int(NeutConventionReactionCode, NeutReacCode->String().Data())
    != PGUtils::STRINT_SUCCESS){
    std::cout << "[WARN]: " << "Couldn't parse reaction code: " <<
      NeutReacCode->String() << std::endl;
  }
}

void NEUT::HandleStdHepParticle(
    UInt_t &StdHepPosition,
    Int_t &StdHepPdg,
    Int_t &StdHepStatus,
    Double_t * &StdHepP4){

  OutObjectInfo->HandleStdHepParticle(StdHepPosition,StdHepPdg,StdHepStatus,
    StdHepP4);

  if(StdHepStatus==11){
    Generator::HandleStruckNucleon(StdHepP4, StdHepPdg);
  }

}
//******************************************************************************
//                     GENIE
//******************************************************************************
void GENIE::Init(TChain* tree){
  StdHepN = &GStdHepN;
  StdHepPdg = GStdHepPdg;
  StdHepP4 = PGUtils::NewPPOf2DArray(GStdHepP4);
  StdHepStatus = GStdHepStatus;

  tree->SetBranchAddress("G2NeutEvtCode", &G2NeutEvtCode);
  tree->SetBranchAddress("StdHepN", &GStdHepN);
  tree->SetBranchAddress("StdHepPdg", GStdHepPdg);
  tree->SetBranchAddress("StdHepP4", GStdHepP4);
  tree->SetBranchAddress("StdHepStatus", GStdHepStatus);
  tree->SetBranchAddress("StdHepRescat", GStdHepRescat);
}

void GENIE::StartEvent(){
  NeutConventionReactionCode = G2NeutEvtCode;
  ProtonRescat_contains_NoInt = false;
  ProtonRescat_contains_chrgEx = false;
  ProtonRescat_contains_elastic = false;
  ProtonRescat_contains_inelastic = false;
  ProtonRescat_contains_knockout = false;
}

void GENIE::HandleStdHepParticle(
    UInt_t &StdHepPosition,
    Int_t &StdHepPdg,
    Int_t &StdHepStatus,
    Double_t * &StdHepP4){

  OutObjectInfo->HandleStdHepParticle(StdHepPosition,StdHepPdg,StdHepStatus,
    StdHepP4);

  if(StdHepStatus == 14 && !LiteOutput){
    HandleRescat(StdHepPdg, GStdHepRescat[StdHepPosition]);
  }

  if(StdHepStatus==11){
    Generator::HandleStruckNucleon(StdHepP4, StdHepPdg);
  }
}

void GENIE::HandleRescat(Int_t PDG, Int_t RescatCode){
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

void GENIE::AddOutputBranches(TTree* tree, bool LiteOutput,
    bool MultiplyByGeVToMeV, Int_t NThresh, Int_t* Threshs_Mev){
  Generator::AddOutputBranches(tree, LiteOutput, MultiplyByGeVToMeV,
    NThresh,Threshs_Mev);
  if(!LiteOutput){
    tree->Branch("ProtonRescat_contains_NoInt",
      &ProtonRescat_contains_NoInt, "ProtonRescat_contains_NoInt/O");
    tree->Branch("ProtonRescat_contains_chrgEx",
      &ProtonRescat_contains_chrgEx, "ProtonRescat_contains_chrgEx/O");
    tree->Branch("ProtonRescat_contains_elastic",
      &ProtonRescat_contains_elastic, "ProtonRescat_contains_elastic/O");
    tree->Branch("ProtonRescat_contains_inelastic",
      &ProtonRescat_contains_inelastic, "ProtonRescat_contains_inelastic/O");
    tree->Branch("ProtonRescat_contains_knockout",
      &ProtonRescat_contains_knockout, "ProtonRescat_contains_knockout/O");
  }
}
//******************************************************************************
//                     NuWro
//******************************************************************************
void NuWro::Init(TChain* tree){
  StdHepN = &NuStdHepN;
  StdHepPdg = NuStdHepPdg;
  StdHepP4 = PGUtils::NewPPOf2DArray(NuStdHepP4);
  StdHepStatus = NuStdHepStatus;

  tree->SetBranchAddress("EvtCode", &NuWroEvtCode);
  tree->SetBranchAddress("StdHepN", &NuStdHepN);
  tree->SetBranchAddress("StdHepPdg", NuStdHepPdg);
  tree->SetBranchAddress("StdHepP4", NuStdHepP4);
  tree->SetBranchAddress("StdHepStatus", NuStdHepStatus);
  tree->SetBranchAddress("EvtWght", &EvtWght);

  ///Need this so that we can keep up to date with the number of entries in
  ///the input tree.
  InpChain = tree;
}

void NuWro::StartEvent(){
  if(PGUtils::str2int(NeutConventionReactionCode, NuWroEvtCode->String().Data())
    != PGUtils::STRINT_SUCCESS){
    std::cout << "[WARN]: " << "Couldn't parse reaction code: " <<
      NuWroEvtCode->String().Data() << std::endl;
  }
}

void NuWro::HandleStdHepParticle(
    UInt_t &StdHepPosition,
    Int_t &StdHepPdg,
    Int_t &StdHepStatus,
    Double_t * &StdHepP4){

  OutObjectInfo->HandleStdHepParticle(StdHepPosition,StdHepPdg,StdHepStatus,
    StdHepP4);

  if(StdHepPosition == 1){

    Int_t StruckNucleonPDGGuess = 0;
    Double_t StruckNucleonMass = sqrt(
      StdHepP4[kStdHepIdxE]*StdHepP4[kStdHepIdxE]
    - StdHepP4[kStdHepIdxPx]*StdHepP4[kStdHepIdxPx]
    - StdHepP4[kStdHepIdxPy]*StdHepP4[kStdHepIdxPy]
    - StdHepP4[kStdHepIdxPz]*StdHepP4[kStdHepIdxPz]);

    if(StruckNucleonMass < 0.939 && StruckNucleonMass > 0.938){
      StruckNucleonPDGGuess = 2212;
    } else if(StruckNucleonMass < 0.940 && StruckNucleonMass > 0.939){
      StruckNucleonPDGGuess = 2112;
    } else if(StruckNucleonMass > 1E-6 &&
              NeutConventionReactionCode == 11){
      std::cout << "[WARN]: Found struck nucleon with mass: "
        << StruckNucleonMass << ", reaction code: "
        << NeutConventionReactionCode << std::endl;
    }
    Generator::HandleStruckNucleon(StdHepP4, StruckNucleonPDGGuess);
  }
}

void NuWro::AddOutputBranches(TTree* tree, bool LiteOutput,
    bool MultiplyByGeVToMeV, Int_t NThresh, Int_t* Threshs_Mev){
  Generator::AddOutputBranches(tree, LiteOutput, MultiplyByGeVToMeV,
    NThresh,Threshs_Mev);
  if(!LiteOutput){
    tree->Branch("EvtWght", &ScaledEvtWght, "EvtWght/D");
  }
}

void NuWro::Finalise(){
  //Moves non Delta++ codes up one.
  if( (OutObjectInfo->IncNeutrino_PDG > 0) &&
      (NeutConventionReactionCode == 11) &&
      (OutObjectInfo->StruckNucleonPDG != 2212) ){

    NeutConventionReactionCode = 12;
  //Forces antineutrino Delta0 code
  } else if( (OutObjectInfo->IncNeutrino_PDG < 0) &&
             (NeutConventionReactionCode == 11) &&
              (OutObjectInfo->StruckNucleonPDG == 2212) ){
    NeutConventionReactionCode = -13;
  }
  ScaledEvtWght = EvtWght/double(InpChain->GetTree()->GetEntries());
  Generator::Finalise();
}
//******************************************************************************
//                     Emulated NuWro
//******************************************************************************

void EmuNuWro::Init(TChain* tree){
  NuWro::Init(tree);
  if(tree->SetBranchAddress("StruckNucleonPDG", &StruckNucleonPDG) !=
    TTree::kMatch){
    std::cout << "[ERROR]: Emulated NuWro tree did not contain auxilliary "
      "StruckNucleonPDG Branch." << std::endl;
    throw 5;
  }
}

void EmuNuWro::HandleStdHepParticle(
    UInt_t &StdHepPosition,
    Int_t &StdHepPdg,
    Int_t &StdHepStatus,
    Double_t * &StdHepP4){

  OutObjectInfo->HandleStdHepParticle(StdHepPosition,StdHepPdg,StdHepStatus,
    StdHepP4);

  if(StdHepPosition == 1){
    Generator::HandleStruckNucleon(StdHepP4, StruckNucleonPDG);
  }

}
//******************************************************************************
//                     GiBUU
//******************************************************************************
void GiBUU::Init(TChain* tree){
  StdHepN = &GiStdHepN;
  StdHepPdg = GiStdHepPdg;
  StdHepP4 = PGUtils::NewPPOf2DArray(GiStdHepP4);
  StdHepStatus = GiStdHepStatus;

  tree->SetBranchAddress("GiBUUReactionCode",&GiBUUReactionCode);
  tree->SetBranchAddress("GiBUUPerWeight",&GiBUUPerWeight);

  tree->SetBranchAddress("GiBUU2NeutCode", &Gi2NeutEvtCode);
  tree->SetBranchAddress("StdHepN", &GiStdHepN);
  tree->SetBranchAddress("StdHepPdg", GiStdHepPdg);
  tree->SetBranchAddress("StdHepP4", GiStdHepP4);
  tree->SetBranchAddress("StdHepStatus", GiStdHepStatus);
}


void GiBUU::AddOutputBranches(TTree* tree, bool LiteOutput,
    bool MultiplyByGeVToMeV, Int_t NThresh, Int_t* Threshs_Mev){
  Generator::AddOutputBranches(tree, LiteOutput, MultiplyByGeVToMeV,
    NThresh,Threshs_Mev);
  tree->Branch("GiBUUReactionCode",
    &GiBUUReactionCode, "GiBUUReactionCode/I");
  tree->Branch("GiBUUPerWeight",
    &GiBUUPerWeight, "GiBUUPerWeight/D");
}

void GiBUU::StartEvent(){
  NeutConventionReactionCode = Gi2NeutEvtCode;
}

void GiBUU::HandleStdHepParticle(
    UInt_t &StdHepPosition,
    Int_t &StdHepPdg,
    Int_t &StdHepStatus,
    Double_t * &StdHepP4){

  OutObjectInfo->HandleStdHepParticle(StdHepPosition,StdHepPdg,StdHepStatus,
    StdHepP4);

  if(StdHepStatus==11){
    Generator::HandleStruckNucleon(StdHepP4, StdHepPdg);
  }

}
