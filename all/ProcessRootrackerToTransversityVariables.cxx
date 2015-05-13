#include <cstdlib>
#include <cerrno>
#include <climits>

#include <algorithm>

#include <sstream>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include "TransversityVariableObjects.hxx"

TChain* RooTrackerChain;

namespace GeneratorDependent {

const int kGStdHepNPmax = 99999;

const int kGStdHepIdxPx = 0;
const int kGStdHepIdxPy = 1;
const int kGStdHepIdxPz = 2;
const int kGStdHepIdxE  = 3;

const int kNStdHepNPmax = 99999;

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

}

enum Generators {kNEUT,kGENIE,kInvalid};

namespace {
  int VERBOSE = 1;

  std::string TreeName="nRooTracker";

  int kStdHepIdxPx = 0;
  int kStdHepIdxPy = 0;
  int kStdHepIdxPz = 0;
  int kStdHepIdxE = 0;

  int MaxStdHepEntries = 0;
  std::string MaxStdHepEntriesStr = "";
  Int_t* StdHepN = 0;
  Int_t* StdHepPdg = 0;
  Double_t** StdHepP4 = 0;
  Int_t* StdHepStatus = 0;

  int NeutConventionReactionCode;
  Generators Generator = kInvalid;
}

//I guess that the pointer array in a statically
//allocated 2D array shouldn't change this is probably alright?
template<typename T, size_t N, size_t M>
T** NewPPOf2DArray(T (&arr)[N][M]){
  T** DynArrOfP = new T*[N];
  for(size_t i =0; i < N; ++i){
    DynArrOfP[i] = arr[i];
  }
  return DynArrOfP;
}


enum STR2INT_ERROR { STRINT_SUCCESS,
                     STRINT_OVERFLOW,
                     STRINT_UNDERFLOW,
                     STRINT_INCONVERTIBLE };

///Converts a string to a long, checking for errors.
///See STR2INT_ERROR for error codes.
STR2INT_ERROR str2int (long &i, char const *s, int base = 0) {
  char *end;
  long  l;
  errno = 0;
  l = strtol(s, &end, base);
  if ((errno == ERANGE && l == LONG_MAX) || l > LONG_MAX) {
      return STRINT_OVERFLOW;
  }
  if ((errno == ERANGE && l == LONG_MIN) || l < LONG_MIN) {
      return STRINT_UNDERFLOW;
  }
  if (*s == '\0' || *end != '\0') {
      return STRINT_INCONVERTIBLE;
  }
  i = l;
  return STRINT_SUCCESS;
}

///Converts a string to a int, checking for errors.
///See STR2INT_ERROR for error codes.
STR2INT_ERROR str2int (int &i, char const *s, int base = 0) {
  long holder;
  STR2INT_ERROR retC = str2int(holder,s,base);
  if(retC != STRINT_SUCCESS){
    return retC;
  }
  if(holder > INT_MAX) {
    return STRINT_OVERFLOW;
  } else if (holder < INT_MIN){
    return STRINT_UNDERFLOW;
  }
  i = holder;
  return retC;
}

template<typename T>
std::string PrintArray(const T* arr, size_t N){
  std::stringstream ss("");
  ss << "[ ";
  for(size_t i = 0; i < N; ++i){
    ss << arr[i] << (((i+1)==N)?"":", ");
  }
  ss << " ]";
  return ss.str();
}

bool SetUpGeneratorDependence(std::string GeneratorName){
  if(GeneratorName == "NEUT"){
    MaxStdHepEntries = GeneratorDependent::kNStdHepNPmax;
    kStdHepIdxPx = GeneratorDependent::kNStdHepIdxPx;
    kStdHepIdxPy = GeneratorDependent::kNStdHepIdxPy;
    kStdHepIdxPz = GeneratorDependent::kNStdHepIdxPz;
    kStdHepIdxE = GeneratorDependent::kNStdHepIdxE;
    StdHepN = &GeneratorDependent::NStdHepN;
    StdHepPdg = GeneratorDependent::NStdHepPdg;
    StdHepP4 = NewPPOf2DArray(GeneratorDependent::NStdHepP4);
    StdHepStatus = GeneratorDependent::NStdHepStatus;
    TreeName = "nRooTracker";
    Generator = kNEUT;
    std::cout << "Working on NEUT Tree: " << TreeName
      << std::endl;
  } else if(GeneratorName == "GENIE"){
    MaxStdHepEntries = GeneratorDependent::kGStdHepNPmax;
    kStdHepIdxPx = GeneratorDependent::kGStdHepIdxPx;
    kStdHepIdxPy = GeneratorDependent::kGStdHepIdxPy;
    kStdHepIdxPz = GeneratorDependent::kGStdHepIdxPz;
    kStdHepIdxE = GeneratorDependent::kGStdHepIdxE;

    StdHepN = &GeneratorDependent::GStdHepN;
    StdHepPdg = GeneratorDependent::GStdHepPdg;
    StdHepP4 = NewPPOf2DArray(GeneratorDependent::GStdHepP4);
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

void ProcessRootrackerToTransversityVariables(
  const char* InputName,
  const char* OutputName="pure_sim_transversity_variables.root",
  std::string GeneratorName="NEUT",
  long long MaxEntries=-1){

  if(!SetUpGeneratorDependence(GeneratorName)){
    return;
  }

  RooTrackerChain = new TChain(TreeName.c_str());
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
      return;
      break;
    }
  }


  TFile* outFile = new TFile(OutputName,"RECREATE");
  TTree* outTreePureSim = new TTree("TransversitudenessPureSim","");

  CCQEFSITransversity* OutInfoCCQEFSI = 0;
  outTreePureSim->Branch("CCQEFSITransversity",&OutInfoCCQEFSI);

  // PionProductionTransversity* OutInfoPionProduction = 0;
  // outTreePureSim->Branch("PionProductionTransversity",&OutInfoPionProduction);

  long long doEntries = (MaxEntries==-1) ?
    RooTrackerChain->GetEntries() :
    (std::min(MaxEntries, RooTrackerChain->GetEntries()));

  for(long long i = 0; i < doEntries; ++i){

    if(!i && (VERBOSE > 0)){
      std::cout << "Getting first entry" << std::endl;
    }
    RooTrackerChain->GetEntry(i);
    if(!i && (VERBOSE > 0)){
      std::cout << "Got first entry" << std::endl;
    }

    if(i && !(i%10000)){
      std::cout << "On entry: " << i << std::endl;
    }

    switch(Generator){
      case kNEUT:{
        // NeutConventionReactionCode = 1;
        str2int(NeutConventionReactionCode,GeneratorDependent::NeutReacCode->Data());
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

    for(UInt_t partNum = 0; partNum < UInt_t(*StdHepN); ++partNum){
      OutInfoCCQEFSI->HandleStdHepParticle(partNum, StdHepPdg[partNum],
        StdHepStatus[partNum], StdHepP4[partNum]);
      // OutInfoPionProduction->HandleStdHepParticle(partNum, StdHepPdg[partNum],
      //   StdHepStatus[partNum], StdHepP4[partNum]);
    }
    OutInfoCCQEFSI->Finalise();

    if((VERBOSE>1) ||
       ( (NeutConventionReactionCode==1) &&
         ( ((VERBOSE>0) && (OutInfoCCQEFSI->NProtons==0)) ||
           ((VERBOSE>0) && (OutInfoCCQEFSI->HMProtonMomentum < 1E-8))
         )
       )
      ){
      std::cout << "%%%%%%(NProtons==" << (OutInfoCCQEFSI->NProtons) << "),  "
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
  std::cout << "Writing tree..." << std::endl;
  outTreePureSim->Write();
  std::cout << "Writing file..." << std::endl;
  outFile->Write();
  std::cout << "Done." << std::endl;
  outFile->Close();
  std::cout << "Shut." << std::endl;
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
  ProcessRootrackerToTransversityVariables(InputName,
                                            OutputName,
                                            GeneratorName,
                                            MaxEntries);

}
