#include <algorithm>

#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TObjString.h"

#include "PureGenUtils.hxx"
#include "CLITools.hxx"

#include "TransversityUtils.hxx"
#include "TransversityVariableObjects.hxx"

using namespace TransversityUtils;

namespace GeneratorDependent {

constexpr int kGStdHepNPmax = 350;
constexpr int kNStdHepNPmax = 100;
constexpr int kNuStdHepNPmax = 4000;
constexpr int kGiStdHepNPmax = 100;

TObjString* NeutReacCode = 0;
Int_t NStdHepN;
Int_t NStdHepPdg[kNStdHepNPmax];
Int_t NStdHepStatus[kNStdHepNPmax];
Double_t NStdHepP4[kNStdHepNPmax][4];

TObjString* NuWroEvtCode = 0;
Int_t NuStdHepN;
Int_t NuStdHepPdg[kNuStdHepNPmax];
Int_t NuStdHepStatus[kNuStdHepNPmax];
Double_t NuStdHepP4[kNuStdHepNPmax][4];

Int_t G2NeutEvtCode;
Int_t GStdHepN;
Int_t GStdHepPdg[kGStdHepNPmax];
Int_t GStdHepStatus[kGStdHepNPmax];
Double_t GStdHepP4[kGStdHepNPmax][4];

Int_t Gi2NeutEvtCode;
Int_t GiStdHepN;
Int_t GiStdHepPdg[kGiStdHepNPmax];
Int_t GiStdHepStatus[kGiStdHepNPmax];
Double_t GiStdHepP4[kGiStdHepNPmax][4];

} // namespace GeneratorDependent


//There'll be no name collisions in my DOJO, might be binary bloat though...
namespace {

  int Verbosity = 1;
  bool OutputInGev;
  int SmearDegrees = 0;
  int Smear_Mom_Mu_MeV = 0;
  int Smear_Mom_P_MeV = 0;
  Double_t TargetBE_MeV = 0;
  bool DoSmear = false;

  Int_t NThresh = 0;
  Int_t Threshs_MeV[kNThreshMax];

  std::string TreeName="nRooTracker";
  TString generatorName;

  Generators Generator = kInvalid;

  constexpr int kStdHepIdxPx = 0;
  constexpr int kStdHepIdxPy = 1;
  constexpr int kStdHepIdxPz = 2;
  constexpr int kStdHepIdxE = 3;

  Int_t* StdHepN = 0;
  Int_t* StdHepPdg = 0;
  Double_t** StdHepP4 = 0;
  Int_t* StdHepStatus = 0;

  int NeutConventionReactionCode;

} // namespace

bool SetUpGeneratorDependence(std::string GeneratorName){

  if((GeneratorName == "NEUT") || (GeneratorName == "neut")){

    generatorName = "NEUT";
    StdHepN = &GeneratorDependent::NStdHepN;
    StdHepPdg = GeneratorDependent::NStdHepPdg;
    StdHepP4 = PGUtils::NewPPOf2DArray(GeneratorDependent::NStdHepP4);
    StdHepStatus = GeneratorDependent::NStdHepStatus;
    TreeName = "nRooTracker";
    Generator = kNEUT;
    std::cout << "Working on NEUT Tree: " << TreeName
      << std::endl;
  } else if((GeneratorName == "GENIE") || (GeneratorName == "genie")){

    generatorName = "GENIE";
    StdHepN = &GeneratorDependent::GStdHepN;
    StdHepPdg = GeneratorDependent::GStdHepPdg;
    StdHepP4 = PGUtils::NewPPOf2DArray(GeneratorDependent::GStdHepP4);
    StdHepStatus = GeneratorDependent::GStdHepStatus;
    TreeName = "gRooTracker";
    Generator = kGENIE;
    std::cout << "Working on GENIE Tree: " << TreeName
      << std::endl;
  } else if((GeneratorName == "NuWro") || (GeneratorName == "nuwro") ||
      (GeneratorName == "NUWRO")){

    generatorName = "NuWro";
    StdHepN = &GeneratorDependent::NuStdHepN;
    StdHepPdg = GeneratorDependent::NuStdHepPdg;
    StdHepP4 = PGUtils::NewPPOf2DArray(GeneratorDependent::NuStdHepP4);
    StdHepStatus = GeneratorDependent::NuStdHepStatus;
    TreeName = "nRooTracker";
    Generator = kNuWro;
    std::cout << "Working on NuWro Tree: " << TreeName
      << std::endl;
  } else if((GeneratorName == "GiBUU") || (GeneratorName == "gibuu") ||
      (GeneratorName == "GIBUU")){

    generatorName = "GiBUU";
    StdHepN = &GeneratorDependent::GiStdHepN;
    StdHepPdg = GeneratorDependent::GiStdHepPdg;
    StdHepP4 = PGUtils::NewPPOf2DArray(GeneratorDependent::GiStdHepP4);
    StdHepStatus = GeneratorDependent::GiStdHepStatus;
    TreeName = "giRooTracker";
    Generator = kGiBUU;
    std::cout << "Working on GiBUU Tree: " << TreeName
      << std::endl;
  } else {
    std::cout << "Unknown Generator string: " << GeneratorName
      << ". Exiting in error." << std::endl;
    return false;
  }

  return true;
}

int ProcessRootrackerToTransversityVariables(
  const char* InputName,
  const char* OutputName="pure_sim_transversity_variables.root",
  std::string GeneratorName="NEUT",
  long long MaxEntries=-1){

  if(!SetUpGeneratorDependence(GeneratorName)){
    return 1;
  }

  TChain* RooTrackerChain = new TChain(TreeName.c_str());
  RooTrackerChain->SetDirectory(0);

  int nFiles = 0, nEntries = 0;
  if( !(nFiles =RooTrackerChain->Add(InputName)) ||
      !(nEntries =RooTrackerChain->GetEntries())){
    std::cout << "[ERROR] Found no files (" << nFiles
      << ") or entries (" << nEntries  << ")" << std::endl;
    return 2;
  }

  switch(Generator){
    case kNEUT:{
      RooTrackerChain->SetBranchAddress("EvtCode",
        &GeneratorDependent::NeutReacCode);
      RooTrackerChain->SetBranchAddress("StdHepN",
        &GeneratorDependent::NStdHepN);
      RooTrackerChain->SetBranchAddress("StdHepPdg",
        GeneratorDependent::NStdHepPdg);
      RooTrackerChain->SetBranchAddress("StdHepP4",
        GeneratorDependent::NStdHepP4);
      RooTrackerChain->SetBranchAddress("StdHepStatus",
        GeneratorDependent::NStdHepStatus);
      break;
    }
    case kGENIE:{
      RooTrackerChain->SetBranchAddress("G2NeutEvtCode",
        &GeneratorDependent::G2NeutEvtCode);
      RooTrackerChain->SetBranchAddress("StdHepN",
        &GeneratorDependent::GStdHepN);
      RooTrackerChain->SetBranchAddress("StdHepPdg",
        GeneratorDependent::GStdHepPdg);
      RooTrackerChain->SetBranchAddress("StdHepP4",
        GeneratorDependent::GStdHepP4);
      RooTrackerChain->SetBranchAddress("StdHepStatus",
        GeneratorDependent::GStdHepStatus);
      break;
    }
    case kNuWro:{
      RooTrackerChain->SetBranchAddress("EvtCode",
        &GeneratorDependent::NuWroEvtCode);
      RooTrackerChain->SetBranchAddress("StdHepN",
        &GeneratorDependent::NuStdHepN);
      RooTrackerChain->SetBranchAddress("StdHepPdg",
        GeneratorDependent::NuStdHepPdg);
      RooTrackerChain->SetBranchAddress("StdHepP4",
        GeneratorDependent::NuStdHepP4);
      RooTrackerChain->SetBranchAddress("StdHepStatus",
        GeneratorDependent::NuStdHepStatus);
      break;
    }
    case kGiBUU:{
      RooTrackerChain->SetBranchAddress("GiBUU2NeutCode",
        &GeneratorDependent::Gi2NeutEvtCode);
      RooTrackerChain->SetBranchAddress("StdHepN",
        &GeneratorDependent::GiStdHepN);
      RooTrackerChain->SetBranchAddress("StdHepPdg",
        GeneratorDependent::GiStdHepPdg);
      RooTrackerChain->SetBranchAddress("StdHepP4",
        GeneratorDependent::GiStdHepP4);
      RooTrackerChain->SetBranchAddress("StdHepStatus",
        GeneratorDependent::GiStdHepStatus);
      break;
    }
    case kInvalid:
    default:{
      std::cerr << "This really shouldn't happen." << std::endl;
      return 4;
      break;
    }
  }

  TFile* outFile = new TFile(OutputName,"CREATE");
  if(!outFile->IsOpen()){
    std::cout << "Failed to open: " << OutputName << std::endl;
    return 8;
  }
  TTree* outTreePureSim = new TTree("TransversitudenessPureSim","");

  TransversityVars* OutInfoCCQEFSI =
    new TransversityVars(OutputInGev, TargetBE_MeV,
      NThresh, Threshs_MeV, SmearDegrees, Smear_Mom_Mu_MeV, Smear_Mom_P_MeV, DoSmear,
      generatorName);
  outTreePureSim->Branch("TransV", &OutInfoCCQEFSI);

  long long doEntries = (MaxEntries==-1) ?
    RooTrackerChain->GetEntries() :
    (std::min(MaxEntries, RooTrackerChain->GetEntries()));

  for(long long i = 0; i < doEntries; ++i){
    RooTrackerChain->GetEntry(i);

    OutInfoCCQEFSI->Reset();
    if(!(i%10000)){
      std::cout << "On entry: " << i  << "/" << doEntries << std::endl;
    }

    switch(Generator){
      case kNuWro:{
        if(PGUtils::str2int(NeutConventionReactionCode,
                  GeneratorDependent::NuWroEvtCode->String().Data())
          != PGUtils::STRINT_SUCCESS){
          std::cout << "[WARN]: " << "Couldn't parse reaction code: " <<
            GeneratorDependent::NuWroEvtCode->String().Data() << std::endl;
        }
        break;
      }
      case kNEUT:{
        if(PGUtils::str2int(NeutConventionReactionCode,
                  GeneratorDependent::NeutReacCode->String().Data())
          != PGUtils::STRINT_SUCCESS){
          std::cout << "[WARN]: " << "Couldn't parse reaction code: " <<
            GeneratorDependent::NeutReacCode->String() << std::endl;
        }
        break;
      }
      case kGENIE:{
        NeutConventionReactionCode = GeneratorDependent::G2NeutEvtCode;
        break;
      }
      case kGiBUU:{
        NeutConventionReactionCode = GeneratorDependent::Gi2NeutEvtCode;
        break;
      }
      case kInvalid:
      default:{
        std::cerr << "This really shouldn't happen." << std::endl;
        break;
      }
    }

    OutInfoCCQEFSI->SetNeutConventionReactionCode(NeutConventionReactionCode);
    for(UInt_t partNum = 0; partNum < UInt_t(*StdHepN); ++partNum){
      OutInfoCCQEFSI->HandleStdHepParticle(partNum, StdHepPdg[partNum],
        StdHepStatus[partNum], StdHepP4[partNum]);

      //Struck Nucleon kinematics are stored differently for the different
      //generators
      switch(Generator){
        case kNuWro:{
          if(partNum == 1){
            TLorentzVector StdHepPTLV = TLorentzVector(
              StdHepP4[partNum][kStdHepIdxPx],
              StdHepP4[partNum][kStdHepIdxPy],
              StdHepP4[partNum][kStdHepIdxPz],
              StdHepP4[partNum][kStdHepIdxE]);
            Double_t StdHepP3Mod = StdHepPTLV.Vect().Mag();
            OutInfoCCQEFSI->HandleStruckNucleon(StdHepPTLV, StdHepP3Mod, 0);
          }
          break;
        }
        case kNEUT:
        case kGENIE:{
          if(StdHepStatus[partNum]==11){
            TLorentzVector StdHepPTLV = TLorentzVector(
              StdHepP4[partNum][kStdHepIdxPx],
              StdHepP4[partNum][kStdHepIdxPy],
              StdHepP4[partNum][kStdHepIdxPz],
              StdHepP4[partNum][kStdHepIdxE]);
            Double_t StdHepP3Mod = StdHepPTLV.Vect().Mag();
            OutInfoCCQEFSI->HandleStruckNucleon(StdHepPTLV, StdHepP3Mod,
              StdHepPdg[partNum]);
          }
          break;
        }
        case kGiBUU:{
          break;
        }
      }

    }
    OutInfoCCQEFSI->Finalise();

    if(Verbosity>1){
      std::cout << "**"
"******************************************************************************"
        << "\n\t#Ev: " << i << "\n"
        << " (OutInfoCCQEFSI->ProtonMom_HighestMomProton == "
        << OutInfoCCQEFSI->HMProton.Momentum << " [MeV/C]): " << "\n**"
"******************************************************************************"
        << std::endl;
      std::cout << "NParticles: " << (*StdHepN) << " - (Neut Reac Code: "
        << NeutConventionReactionCode << ")" << std::endl;
      std::cout << "Struck Nucleon { Momentum: "
        << OutInfoCCQEFSI->StruckNucleon.Momentum
        << ", PDG: " << OutInfoCCQEFSI->StruckNucleon.PDG << "}" << std::endl;

        std::cout << "------\nTransverse Variables"
        << "\n\t DeltaP:" << OutInfoCCQEFSI->DeltaPTotal_HMProton_MeV.Mag()
        << "\n\tP_0^p: " <<  OutInfoCCQEFSI->DeltaPProton_MeV.Mag()
        << "\n------" << std::endl;

      for(int partNum = 0; partNum < (*StdHepN); ++partNum){
        std::cout << "\t" << partNum << ": " << StdHepPdg[partNum]
          << " (Status==" << StdHepStatus[partNum] << ") "
          << PGUtils::PrintArray(StdHepP4[partNum],4) << std::endl;
      }
    }

    outTreePureSim->Fill();
  }
  outTreePureSim->Write();
  outFile->Write();

  std::cout << "Wrote " << outTreePureSim->GetEntries() << " entries to disk."
    << std::endl;

  outFile->Close();

  delete RooTrackerChain;
  delete outFile;
  delete [] StdHepP4;
  delete OutInfoCCQEFSI;
}

namespace {
std::string InputName;
std::string OutputName="TransverseVars.root";
std::string GeneratorName="NEUT";
long MaxEntries=-1;

void SetOpts(){
  CLIArgs::OptSpec.emplace_back("-h","--help", false,
    [&] (std::string const &opt) -> bool {
      CLIArgs::SayRunLike();
      exit(0);
    });

  CLIArgs::OptSpec.emplace_back("-i", "--input-file", true,
    [&] (std::string const &opt) -> bool {
      std::cout << "\t--Reading from file descriptor : " << opt << std::endl;
      InputName = opt;
      return true;
    }, true,[](){},"<TChain::Add descriptor>");

  CLIArgs::OptSpec.emplace_back("-o", "--output-file", true,
    [&] (std::string const &opt) -> bool {
      std::cout << "\t--Writing to File: " << opt << std::endl;
      OutputName = opt;
      return true;
    }, false,
    [&](){OutputName = "TransverseVars.root";},
    "<File Name>{default=TransverseVars.root}");

  CLIArgs::OptSpec.emplace_back("-n", "--nentries", true,
    [&] (std::string const &opt) -> bool {
      long vbhold;
      if(PGUtils::str2int(vbhold,opt.c_str()) == PGUtils::STRINT_SUCCESS){
        if(vbhold != -1){
          std::cout << "\t--Looking at, at most, " << vbhold << " entries."
            << std::endl;
        }
        MaxEntries = vbhold;
        return true;
      }
      return false;
    }, false,
    [&](){MaxEntries = -1;}, "<-1>: Read all {default=-1}");

  CLIArgs::OptSpec.emplace_back("-v", "--verbosity", true,
    [&] (std::string const &opt) -> bool {
      int vbhold;
      if(PGUtils::str2int(vbhold,opt.c_str()) == PGUtils::STRINT_SUCCESS){
        std::cout << "\t--Verbosity: " << vbhold << std::endl;
        Verbosity = vbhold;
        return true;
      }
      return false;
    }, false,
    [&](){Verbosity = 0;}, "<0-4>{default=0}");
  CLIArgs::OptSpec.emplace_back("-M", "--MeV-mode", false,
    [&] (std::string const &opt) -> bool {
      std::cout << "\t--Outputting in MeV." << std::endl;
      OutputInGev = false;
      return true;
    }, false,
    [&](){OutputInGev = true;}, "Use MeV rather than GeV.{default=false}");

  CLIArgs::OptSpec.emplace_back("-g", "--generator", true,
    [&] (std::string const &opt) -> bool {
      std::cout << "\t--Attempting to read generator: " << opt << std::endl;
      GeneratorName = opt;
      return true;
    }, false,
    [&](){GeneratorName = "NEUT";}, "{default=NEUT}");

  CLIArgs::OptSpec.emplace_back("-m", "--EKin-Threshold", true,
    [&] (std::string const &opt) -> bool {
      if(NThresh==kNThreshMax){
        std::cout << "[ERROR]: Tried to add too many momentum thresholds."
          << std::endl;
        return false;
      }
      int vbhold;
      if(PGUtils::str2int(vbhold,opt.c_str()) == PGUtils::STRINT_SUCCESS){
        std::cout << "\t--Added EKin threshold: " << vbhold << " MeV" << std::endl;
        if(NThresh>0 && (vbhold <= Threshs_MeV[NThresh-1])){
          std::cout << "[ERROR]: Attempting to add EKin threshold at " << vbhold
            << " MeV, but the previous one at " << Threshs_MeV[NThresh-1]
            << " is lower.\n\tThresholds must be added in ascending order."
            << std::endl;
          return false;
        }
        Threshs_MeV[NThresh] = vbhold;
        NThresh++;
        return true;
      }
      return false;
    }, false,
    [&](){}, "<int> Add EKin threshold [MeV] {default=N/A}");

  CLIArgs::OptSpec.emplace_back("-P", "--Phi-Smear", true,
    [&] (std::string const &opt) -> bool {
      int vbhold;
      if(PGUtils::str2int(vbhold,opt.c_str()) == PGUtils::STRINT_SUCCESS){
        std::cout << "\t--Smear track angle by: " << vbhold << " degrees" << std::endl;
        SmearDegrees = vbhold;
        DoSmear = true;
        return true;
      }
      return false;
    }, false,
    [&](){SmearDegrees = 0; DoSmear = false;}, "<int> Gaussian smear width [degrees] {default=0}");

  CLIArgs::OptSpec.emplace_back("-Sp", "--Mom-Proton-Smear", true,
    [&] (std::string const &opt) -> bool {
      int vbhold;
      if(PGUtils::str2int(vbhold,opt.c_str()) == PGUtils::STRINT_SUCCESS){
        std::cout << "\t--Smear measured proton momentum by: " << vbhold << " MeV" << std::endl;
        Smear_Mom_P_MeV = vbhold;
        DoSmear = true;
        return true;
      }
      return false;
    }, false,
    [&](){Smear_Mom_P_MeV = 0;}, "<int> Gaussian smear Proton Momenta by [MeV] {default=0.0}");

  CLIArgs::OptSpec.emplace_back("-Smu", "--Mom-Muon-Smear", true,
    [&] (std::string const &opt) -> bool {
      int vbhold;
      if(PGUtils::str2int(vbhold,opt.c_str()) == PGUtils::STRINT_SUCCESS){
        std::cout << "\t--Smear measured muon momentum by: " << vbhold << " MeV" << std::endl;
        Smear_Mom_Mu_MeV = vbhold;
        DoSmear = true;
        return true;
      }
      return false;
    }, false,
    [&](){Smear_Mom_P_MeV = 0;}, "<int> Gaussian smear Muon Momenta by [MeV] {default=0.0}");

  CLIArgs::OptSpec.emplace_back("-B", "--Binding-Energy", true,
    [&] (std::string const &opt) -> bool {
      float vbhold = std::stof(opt);
      if(vbhold){
        std::cout << "\t--Using  " << vbhold << " MeV binding energy in ERec"
          << std::endl;
        TargetBE_MeV = vbhold;
        return true;
      }
      return false;
    }, false,
    [&](){TargetBE_MeV = 0xdeadbeef;},
    "<float> Binding energy used in nu_erec calculations [MeV] {default=25.0}");
}
}


int main(int argc, char const * argv[]){
  SetOpts();

  CLIArgs::AddArguments(argc,argv);
  if(!CLIArgs::GetOpts()){
    CLIArgs::SayRunLike();
    return 1;
  }

  std::cout << "Running with " << NThresh << " set momentum thresholds."
    << std::endl;

  return ProcessRootrackerToTransversityVariables(InputName.c_str(),
                                                  OutputName.c_str(),
                                                  GeneratorName.c_str(),
                                                  MaxEntries);

}
