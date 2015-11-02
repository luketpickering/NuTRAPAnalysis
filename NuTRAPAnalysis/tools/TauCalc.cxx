#include "TChain.h"

#include "CLITools.hxx"

std::string InputFileDescriptor = "";

void SetOpts(){

  CLIArgs::AddOpt("-i", "--input", true,
    [&] (std::string const &opt) -> bool {
      InputFileDescriptor = opt;
      return true;
    }, true,
    [](){},
    "File to read input tree from.");
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

  float N_QE = TRAPChain->GetEntries("(TransV.NeutConventionReactionCode==1)&&(TransV.Muon_PDG==13)&&(TransV.HMProton_PDG==2212)");
  float N_QE_DP = TRAPChain->GetEntries("(TransV.DeltaPTotal_HMProton_MeV.Vect().Mag()>10)&&(TransV.NeutConventionReactionCode==1)&&(TransV.Muon_PDG==13)&&(TransV.HMProton_PDG==2212)");


  float N_RES = TRAPChain->GetEntries("(TransV.NeutConventionReactionCode==11)&&(TransV.Muon_PDG==13)&&(TransV.HMProton_PDG==2212)&&(TransV.HMCPion_PDG==211)");
  float N_RES_DP = TRAPChain->GetEntries("(TransV.DeltaPTotal_HMProtonPion_MeV.Vect().Mag()>10)&&(TransV.NeutConventionReactionCode==11)&&(TransV.Muon_PDG==13)&&(TransV.HMProton_PDG==2212)&&(TransV.HMCPion_PDG==211)");


  std::cout << "In File: " << InputFileDescriptor << std::endl;
  std::cout << "Taus:" << std::endl;
  std::cout << "\tQE: " << (N_QE_DP/N_QE) << " = (" << N_QE_DP << "/" << N_QE << ")" << std::endl;
  std::cout << "\tQE: " << (N_RES_DP/N_RES) << " = (" << N_RES_DP << "/" << N_RES << ")" << std::endl;

  delete TRAPChain;
  return 0;
}
