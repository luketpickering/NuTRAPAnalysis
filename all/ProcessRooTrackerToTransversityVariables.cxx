#include <algorithm>

#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include "TransversityUtils.hxx"

#include "TransversityVariableObjects.hxx"

namespace GeneratorDependent {

const int kGStdHepNPmax = 350;

const int kGStdHepIdxPx = 0;
const int kGStdHepIdxPy = 1;
const int kGStdHepIdxPz = 2;
const int kGStdHepIdxE  = 3;

const int kNStdHepNPmax = 100;

const int kNStdHepIdxPx = 0;
const int kNStdHepIdxPy = 1;
const int kNStdHepIdxPz = 2;
const int kNStdHepIdxE = 3;

TString* NeutReacCode = 0;
Int_t NStdHepN;
Int_t NStdHepPdg[kNStdHepNPmax];
Int_t NStdHepStatus[kNStdHepNPmax];
Double_t NStdHepP4[kNStdHepNPmax][4];

Int_t G2NeutEvtCode;
Int_t GStdHepN;
Int_t GStdHepPdg[kGStdHepNPmax];
Int_t GStdHepStatus[kGStdHepNPmax];
Double_t GStdHepP4[kGStdHepNPmax][4];

} // namespace GeneratorDependent


//There'll be no name collisions in my DOJO, might be binary bloat though...
namespace {
  int VERBOSE = 1;

  std::string TreeName="nRooTracker";

  enum Generators {kNEUT,kGENIE,kInvalid};
  Generators Generator = kInvalid;


//The following variables are the ones that the thinking code will use to get
// at the specific generator output.
  int kStdHepIdxPx = 0;
  int kStdHepIdxPy = 0;
  int kStdHepIdxPz = 0;
  int kStdHepIdxE = 0;

  Int_t* StdHepN = 0;
  Int_t* StdHepPdg = 0;
  Double_t** StdHepP4 = 0;
  Int_t* StdHepStatus = 0;

  int NeutConventionReactionCode;

} // namespace

bool SetUpGeneratorDependence(std::string GeneratorName){

  if(GeneratorName == "NEUT"){
    kStdHepIdxPx = GeneratorDependent::kNStdHepIdxPx;
    kStdHepIdxPy = GeneratorDependent::kNStdHepIdxPy;
    kStdHepIdxPz = GeneratorDependent::kNStdHepIdxPz;
    kStdHepIdxE = GeneratorDependent::kNStdHepIdxE;

    StdHepN = &GeneratorDependent::NStdHepN;
    StdHepPdg = GeneratorDependent::NStdHepPdg;
    StdHepP4 = TransversityUtils::NewPPOf2DArray(GeneratorDependent::NStdHepP4);
    StdHepStatus = GeneratorDependent::NStdHepStatus;
    TreeName = "nRooTracker";
    Generator = kNEUT;
    std::cout << "Working on NEUT Tree: " << TreeName
      << std::endl;
  } else if(GeneratorName == "GENIE"){
    kStdHepIdxPx = GeneratorDependent::kGStdHepIdxPx;
    kStdHepIdxPy = GeneratorDependent::kGStdHepIdxPy;
    kStdHepIdxPz = GeneratorDependent::kGStdHepIdxPz;
    kStdHepIdxE = GeneratorDependent::kGStdHepIdxE;

    StdHepN = &GeneratorDependent::GStdHepN;
    StdHepPdg = GeneratorDependent::GStdHepPdg;
    StdHepP4 = TransversityUtils::NewPPOf2DArray(GeneratorDependent::GStdHepP4);
    StdHepStatus = GeneratorDependent::GStdHepStatus;
    TreeName = "gRooTracker";
    Generator = kGENIE;
    std::cout << "Working on GENIE Tree: " << TreeName
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

  RooTrackerChain->Add(InputName);

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
    case kInvalid:
    default:{
      std::cerr << "This really shouldn't happen." << std::endl;
      return 1;
      break;
    }
  }


  TFile* outFile = new TFile(OutputName,"CREATE");
  if(!outFile->IsOpen()){
    std::cout << "Failed to open: " << OutputName << std::endl;
    return 2;
  }
  TTree* outTreePureSim = new TTree("TransversitudenessPureSim","");

  CCQEFSITransversity* OutInfoCCQEFSI = 0;
  outTreePureSim->Branch("CCQEFSITransversity", &OutInfoCCQEFSI);

  PionProductionTransversity* OutInfoPionProduction = 0;
  outTreePureSim->Branch("PionProductionTransversity", &OutInfoPionProduction);

  long long doEntries = (MaxEntries==-1) ?
    RooTrackerChain->GetEntries() :
    (std::min(MaxEntries, RooTrackerChain->GetEntries()));

  for(long long i = 0; i < doEntries; ++i){
    RooTrackerChain->GetEntry(i);

    OutInfoCCQEFSI->Reset();
    OutInfoPionProduction->Reset();

    if(i && !(i%5000)){
      std::cout << "On entry: " << i << std::endl;
    }

    switch(Generator){
      case kNEUT:{
        if(TransversityUtils::str2int(NeutConventionReactionCode,
                  GeneratorDependent::NeutReacCode->Data())
          != TransversityUtils::STRINT_SUCCESS){
          std::cout << "[WARN]: " << "Couldn't parse reaction code: " <<
            GeneratorDependent::NeutReacCode->Data() << std::endl;
        }
        break;
      }
      case kGENIE:{
        NeutConventionReactionCode = GeneratorDependent::G2NeutEvtCode;
        break;
      }
      case kInvalid:
      default:{
        std::cerr << "This really shouldn't happen." << std::endl;
        break;
      }
    }

    OutInfoCCQEFSI->SetNeutConventionReactionCode(NeutConventionReactionCode);
    OutInfoPionProduction->\
      SetNeutConventionReactionCode(NeutConventionReactionCode);

    for(UInt_t partNum = 0; partNum < UInt_t(*StdHepN); ++partNum){
      OutInfoCCQEFSI->HandleStdHepParticle(partNum, StdHepPdg[partNum],
        StdHepStatus[partNum], StdHepP4[partNum]);
      OutInfoPionProduction->HandleStdHepParticle(partNum, StdHepPdg[partNum],
        StdHepStatus[partNum], StdHepP4[partNum]);
    }
    OutInfoCCQEFSI->Finalise();
    OutInfoPionProduction->Finalise();

    if((VERBOSE>1)){
      std::cout << "(NProtons==" << (OutInfoCCQEFSI->NProtons) << "),  "
        << " (OutInfoCCQEFSI->ProtonMom_HighestMomProton == "
        << OutInfoCCQEFSI->HMProtonMomentum << "): " << std::endl;
      std::cout << "NParticles: " << (*StdHepN) << " - ("
        << NeutConventionReactionCode << ")" << std::endl;
      for(int partNum = 0; partNum < (*StdHepN); ++partNum){
        std::cout << "\t" << partNum << ": " << StdHepPdg[partNum]
          << " (Status==" << StdHepStatus[partNum] << ") "
          << PrintArray(StdHepP4[partNum],4) << std::endl;
      }
    }

    outTreePureSim->Fill();
  }
  outTreePureSim->Write();
  outFile->Write();
  outFile->Close();

  delete [] StdHepP4;
}

void SayUsage(char const *runcmd){
  std::cout << "Run Like:\n\t" << runcmd << " <Input-File-Descriptor> "
  "[<OutputName>=pure_sim_transversity_variables.root] [<GeneratorName>=NEUT]"
  " [<int:MaximumEntriesToProcess>=-1]"
  "\n\tMaximumEntriesToProcess == -1 means process all found in input files."
  "\n\tIf you provide a argument you must provide all preceeding ones."
  << std::endl;
}

int main(int argc, char* argv[]){


  if(argc < 2){
    SayUsage(argv[0]);
    return 1;
  }
  const char* InputName = argv[1];
  const char* OutputName="pure_sim_transversity_variables.root";
  const char* GeneratorName="NEUT";
  long MaxEntries=-1;

  if(argc >= 3){
    OutputName = argv[2];
  }
  if(argc >= 4){
    GeneratorName = argv[3];
  }
  if(argc >= 5){
    str2int(MaxEntries,argv[4]);
  }
  if(argc >= 6){
    str2int(VERBOSE,argv[5]);
  }
  return ProcessRootrackerToTransversityVariables(InputName,
                                                  OutputName,
                                                  GeneratorName,
                                                  MaxEntries);

}
