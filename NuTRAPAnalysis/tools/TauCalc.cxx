#include <cmath>
#include <iomanip>

#include "TChain.h"
#include "TH1D.h"

#include "CLITools.hxx"
#include "PureGenUtils.hxx"

#include "TransversityVariableObjects.hxx"

std::string InputFileDescriptor = "";
Double_t BindingEnergy;
bool IsLiteMode = false;


template<typename TH>
TH* ToPDF(const TH *hraw){

  const Double_t EPSILON = 1e-12;
  const Int_t x0 = 0;
  const Int_t x1 = hraw->GetNbinsX()+1;
  const Double_t Integral = hraw->Integral(x0, x1);

  TH * hist = (TH*)hraw->Clone((std::string(hraw->GetName())+"pdf").c_str());
  hist->Scale(0);

  for(Int_t ib = x0; ib <= x1; ib++){
    const Double_t BinWidth = hraw->GetBinWidth(ib);
    const Double_t BinContent = hraw->GetBinContent(ib);

    //in case of finit number of bins (i.e. eff not always small),
    //Binomial error is more accurate than Poisson error
    const Double_t Scaled = BinContent/Integral;
    const Double_t PDF = Scaled/BinWidth;

    const Double_t PDFErr = sqrt(Scaled*(1-Scaled)/Integral) / BinWidth;
    hist->SetBinContent(ib, PDF);
    hist->SetBinError(ib, PDFErr);
  }

  hist->SetEntries(Integral);

  return hist;
}

void SetOpts(){

  CLIArgs::AddOpt("-i", "--input", true,
    [&] (std::string const &opt) -> bool {
      InputFileDescriptor = opt;
      return true;
    }, true,
    [](){},
    "File to read input tree from.");

    CLIArgs::AddOpt("-B", "--Binding-Energy", true,
    [&] (std::string const &opt) -> bool {
      BindingEnergy = std::stod(opt);
      return true;
    }, true,
    [](){},
    "Input tree is built from TransversityVarsLite instances.");

    CLIArgs::AddOpt("-L", "--Lite-Mode", false,
    [&] (std::string const &opt) -> bool {
      IsLiteMode = PGUtils::str2bool(opt);
      return true;
    }, false,
    [](){IsLiteMode = false;},
    "Input tree is built from TransversityVarsLite instances.");
}

Long64_t LoopEvents(TChain *TRAPChain,
  std::function<bool(TransversityVarsB const &mytv)> CallBack,
  Long64_t NMax=-1){

  TransversityVarsB const * const mytv =
    MakeReadingTransversityVarsB(TRAPChain);
  Long64_t SelEntries = 0;
  Long64_t MaxLoop = (NMax==-1)?TRAPChain->GetEntries():
                                (std::min(NMax,TRAPChain->GetEntries()) );
  for(Long64_t ent = 0; ent < TRAPChain->GetEntries(); ++ent){
    TRAPChain->GetEntry(ent);
    SelEntries += CallBack(*mytv);
  }
  UnsetBranchAddressesTransversityVarsB(TRAPChain,mytv);
  return SelEntries;
}

int main(int argc, char const * argv[]){

  try {
    SetOpts();
  } catch (std::exception const & e){
    std::cerr << "[ERROR]: " << e.what() << std::endl;
    return 1;
  }

  CLIArgs::AddArguments(argc,argv);
  if(!CLIArgs::HandleArgs()){
    CLIArgs::SayRunLike();
    return 1;
  }

  TChain *TRAPChain = new TChain("TransversitudenessPureSim");
  TRAPChain->Add(InputFileDescriptor.c_str());

  auto QESel = [](TransversityVarsB const &mytv) -> bool {
      return (
        (mytv.NeutConventionReactionCode==1) &&
        (mytv.Muon_PDG==13) &&
        (mytv.HMProton_PDG==2212) );
    };

  auto QESel_NegDE = [&](TransversityVarsB const &mytv) -> bool {
      return (
        QESel(mytv) &&
        (-1.0*(mytv.Deltap_HMProton_MeV.E()+BindingEnergy) < 0.01) );
    };

  auto QESel_DP = [&](TransversityVarsB const &mytv) -> bool {
      return (
        QESel(mytv) &&
        (mytv.Deltap_HMProton_MeV.Vect().Mag() > 10.0) );
    };

  auto QESel_NE = [&](TransversityVarsB const &mytv) -> bool {
      return (
        QESel(mytv) &&
        (mytv.NFinalStateParticles!=2) );
    };

  auto QESel_NoNE = [&](TransversityVarsB const &mytv) -> bool {
      return (
        QESel(mytv) &&
        (mytv.NFinalStateParticles==2) );
    };

  auto RESSel = [](TransversityVarsB const &mytv) -> bool {
      return (
        (mytv.NeutConventionReactionCode==11) &&
        (mytv.Muon_PDG==13) &&
        (mytv.HMProton_PDG==2212) &&
        (mytv.HMCPion_PDG==211) );
    };

  auto RESSel_DP = [&](TransversityVarsB const &mytv) -> bool {
      return (
        RESSel(mytv) &&
        (mytv.Deltap_HMProtonPion_MeV.Vect().Mag() > 10.0) );
    };

  Double_t N_QE = LoopEvents(TRAPChain, [&](TransversityVarsB const &mytv){
    if(QESel(mytv)){ return true;}
    return false;
  });
  Double_t N_QE_DP = LoopEvents(TRAPChain, [&](TransversityVarsB const &mytv){
    if(QESel_DP(mytv)){ return true;}
    return false;
  });

  Double_t N_RES = LoopEvents(TRAPChain, [&](TransversityVarsB const &mytv){
    if(RESSel(mytv)){ return true;}
    return false;
  });
  Double_t N_RES_DP = LoopEvents(TRAPChain, [&](TransversityVarsB const &mytv){
    if(RESSel_DP(mytv)){ return true;}
    return false;
  });

  std::cout << "In File: " << InputFileDescriptor << std::endl;
  std::cout << "Taus:" << std::endl;
  std::cout << "\tQE: " << (N_QE_DP/N_QE) << " = (" << N_QE_DP << "/"
    << N_QE << ")" << std::endl;
  std::cout << "\tRES: " << (N_RES_DP/N_RES) << " = (" << N_RES_DP << "/"
    << N_RES << ")" << std::endl;

  delete TRAPChain;
  return 0;
}
