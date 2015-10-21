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
Int_t GStdHepRescat[kGStdHepNPmax];
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
  bool LiteOutput = false;

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
  Int_t* StdHepRescat = 0;

  Float_t TargetBE_MeV = 0xdeadbeef;

  int NeutConventionReactionCode;

  std::vector<int> ModeIgnores;

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
    StdHepRescat = GeneratorDependent::GStdHepRescat;
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
      RooTrackerChain->SetBranchAddress("StdHepRescat",
        GeneratorDependent::GStdHepRescat);
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

  TransversityVarsB* OutObjectInfo;
  if(LiteOutput) {
    TransversityVarsLite* OutObjectInfo_obj = new TransversityVarsLite(OutputInGev);
    outTreePureSim->Branch("TransV", &OutObjectInfo_obj);
    OutObjectInfo = OutObjectInfo_obj;

  } else {
    TransversityVars* OutObjectInfo_obj =
      new TransversityVars(OutputInGev, TargetBE_MeV,
        NThresh, Threshs_MeV, generatorName);
    outTreePureSim->Branch("TransV", &OutObjectInfo_obj);
    OutObjectInfo = OutObjectInfo_obj;
  }

  long long doEntries = (MaxEntries==-1) ?
    RooTrackerChain->GetEntries() :
    (std::min(MaxEntries, RooTrackerChain->GetEntries()));

  for(long long i = 0; i < doEntries; ++i){
    RooTrackerChain->GetEntry(i);

    OutObjectInfo->Reset();
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

    for(UInt_t partNum = 0; partNum < UInt_t(*StdHepN); ++partNum){
      OutObjectInfo->HandleStdHepParticle(partNum, StdHepPdg[partNum],
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
            Int_t StruckNucleonPDGGuess = 0;
            Double_t StruckNucleonMass = StdHepPTLV.M();
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
            OutObjectInfo->HandleStruckNucleon(StdHepPTLV, StdHepP3Mod,
              StruckNucleonPDGGuess);
          }
          break;
        }
        case kGENIE:{
          if(LiteOutput){
            if(StdHepStatus[partNum] == 14){
              OutObjectInfo->HandleRescat(StdHepPdg[partNum],
                StdHepRescat[partNum]);
            }
          }
          break;
        }
        case kNEUT:{
          if(StdHepStatus[partNum]==11){
            TLorentzVector StdHepPTLV = TLorentzVector(
              StdHepP4[partNum][kStdHepIdxPx],
              StdHepP4[partNum][kStdHepIdxPy],
              StdHepP4[partNum][kStdHepIdxPz],
              StdHepP4[partNum][kStdHepIdxE]);
            Double_t StdHepP3Mod = StdHepPTLV.Vect().Mag();
            OutObjectInfo->HandleStruckNucleon(StdHepPTLV, StdHepP3Mod,
              StdHepPdg[partNum]);
          }
          break;
        }
        case kGiBUU:{
          break;
        }
      }

    }

    OutObjectInfo->Finalise();

    //Bumps up the code for non-Delta++ resonances for incoming numus.
    if((OutObjectInfo->GetIncNeutrino_PDG() == 14) &&
      (Generator == kNuWro) &&
      (NeutConventionReactionCode == 11) &&
      (OutObjectInfo->GetStruckNucleonPDG() != 2212) ){
      OutObjectInfo->SetNeutConventionReactionCode(12);
    } else {
      OutObjectInfo->SetNeutConventionReactionCode(NeutConventionReactionCode);
    }

    if(ModeIgnores.size()){
      bool found = false;
      for(auto const & mi : ModeIgnores){
        if(mi == OutObjectInfo->GetNeutConventionReactionCode()){
          found = true;
          break;
        }
      }
      if(!found){
        continue;
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
  delete OutObjectInfo;
}

namespace {
std::string InputName;
std::string OutputName="TransverseVars.root";
std::string GeneratorName="NEUT";
long MaxEntries=-1;

void SetOpts(){

  CLIArgs::AddOpt("-i", "--input-file", true,
    [&] (std::string const &opt) -> bool {
      std::cout << "\t--Reading from file descriptor : " << opt << std::endl;
      InputName = opt;
      return true;
    }, true,[](){},"<TChain::Add descriptor>");

  CLIArgs::AddOpt("-o", "--output-file", true,
    [&] (std::string const &opt) -> bool {
      std::cout << "\t--Writing to File: " << opt << std::endl;
      OutputName = opt;
      return true;
    }, false,
    [&](){OutputName = "TransverseVars.root";},
    "<File Name>{default=TransverseVars.root}");

  CLIArgs::AddOpt("-n", "--nentries", true,
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

  CLIArgs::AddOpt("-v", "--verbosity", true,
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
  CLIArgs::AddOpt("-M", "--MeV-mode", false,
    [&] (std::string const &opt) -> bool {
      std::cout << "\t--Reading in MeV." << std::endl;
      OutputInGev = false;
      return true;
    }, false,
    [&](){OutputInGev = true;}, "Assume input is in MeV.{default=false}");

  CLIArgs::AddOpt("-g", "--generator", true,
    [&] (std::string const &opt) -> bool {
      std::cout << "\t--Attempting to read generator: " << opt << std::endl;
      GeneratorName = opt;
      return true;
    }, false,
    [&](){GeneratorName = "NEUT";}, "{default=NEUT}");

  CLIArgs::AddOpt("-m", "--EKin-Threshold", true,
    [&] (std::string const &opt) -> bool {
      if(NThresh==kNThreshMax){
        std::cout << "[ERROR]: Tried to add too many momentum thresholds."
          << std::endl;
        return false;
      }
      int vbhold;
      if(PGUtils::str2int(vbhold,opt.c_str()) == PGUtils::STRINT_SUCCESS){
        std::cout << "\t--Added EKin threshold: " << vbhold << " MeV" << std::endl;
        if(NThresh>0 && (vbhold < Threshs_MeV[NThresh-1])){
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

    CLIArgs::AddOpt("-N", "--NEUT-Modes", true,
    [&] (std::string const &opt) -> bool {
      ModeIgnores =
        PGUtils::StringVToIntV(PGUtils::SplitStringByDelim(opt,","));

      if(ModeIgnores.size()){
        std::cout << "\t--Ignoring interactions except of the modes:  "
          << std::flush;
        for(auto const &mi : ModeIgnores){
          std::cout << mi << ", " << std::flush;
        }
        std::cout << std::endl;
        return true;
      }
      return false;
    }, false,
    [](){},
    "<int,int,...> NEUT modes to save output from.");

  CLIArgs::AddOpt("-B", "--Binding-Energy", true,
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

    CLIArgs::AddOpt("-L", "--Lite-Output", false,
    [&] (std::string const &opt) -> bool {
      LiteOutput = true;
      std::cout << "\t--Outputting Lite format." << std::endl;
      return true;
    }, false,
    [](){LiteOutput = false;},
    "Will output in Lite mode which contains less output variables.");

}
}

int main(int argc, char const * argv[]){
  SetOpts();

  CLIArgs::AddArguments(argc,argv);
  if(!CLIArgs::HandleArgs()){
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
