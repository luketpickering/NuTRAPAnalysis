#include <cmath>
#include <iomanip>

#include "TChain.h"
#include "TH1D.h"

#include "CLITools.hxx"
#include "PureGenUtils.hxx"

#include "EventLoop.hxx"

namespace {
  std::string InputFileDescriptor = "";
  Int_t NEvs = 0;
}

void SetOpts(){

  CLIArgs::AddOpt("-i", "--input", true,
    [&] (std::string const &opt) -> bool {
      InputFileDescriptor = opt;
      return true;
    }, true,
    [](){},
    "File to read input tree from.");

  CLIArgs::AddOpt("-n", "--nevs", true,
  [&] (std::string const &opt) -> bool {
    NEvs = std::stoi(opt);
    return true;
  }, false,
  [](){NEvs=-1;},
  "[Read at most N events. \'-1\' means read all. {Default:-1}]");
}

Long64_t LoopEvents(TChain *TRAPChain,
  std::function<bool(Long64_t LoopNum, TransversityVarsB const &tv)> CallBack,
  Long64_t NMax=NEvs){

  TransversityVarsB const * const tv =
    MakeReadingTransversityVarsB(TRAPChain);
  Long64_t SelEntries = 0;
  Long64_t MaxLoop = (NMax==-1)?TRAPChain->GetEntries():
                                (std::min(NMax,TRAPChain->GetEntries()) );
  for(Long64_t ent = 0; ent < MaxLoop; ++ent){
    TRAPChain->GetEntry(ent);
    SelEntries += CallBack(ent, *tv);
  }
  UnsetBranchAddressesTransversityVarsB(TRAPChain,tv);
  return SelEntries;
}

int RunLoop(int argc,
  char const * argv[],
  std::function<bool(Long64_t LoopNum, TransversityVarsB const &tv)> EvLoopFunc){

  try {
    SetOpts();
  } catch (std::exception const & e){
    std::cerr << "[ERROR]: " << e.what() << std::endl;
    return -1;
  }

  CLIArgs::AddArguments(argc,argv);
  if(!CLIArgs::HandleArgs()){
    CLIArgs::SayRunLike();
    return -1;
  }

  TChain *TRAPChain = new TChain("TransversitudenessPureSim");
  TRAPChain->Add(InputFileDescriptor.c_str());

  Long64_t rtn = LoopEvents(TRAPChain,EvLoopFunc);

  delete TRAPChain;
  return rtn;
}
