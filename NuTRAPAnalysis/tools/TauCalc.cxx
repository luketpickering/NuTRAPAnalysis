#include <iostream>

#include "CLITools.hxx"

#include "EventLoop.hxx"

int main(int argc, char const * argv[]){

  Double_t DeltaPCut = 0;

  try {
    CLIArgs::AddOpt("-P", "--DeltaP-Threshold", true,
      [&] (std::string const &opt) -> bool {
        DeltaPCut = stod(opt);
        return true;
      }, false,
      [&](){DeltaPCut = 10.0;},
      "DeltaP threshold to count as having undergone FSI.");
  } catch (std::exception const & e){
    std::cerr << "[ERROR]: " << e.what() << std::endl;
    return 1;
  }

  int argc2; char const ** argv2;
  CLIArgs::AddArgumentsAndGetUnexpecteds(argc,argv,argc2,argv2);

  auto QESel = [&](TransversityVarsB const &mytv) -> bool {
      return (
        (mytv.NeutConventionReactionCode==1) &&
        (mytv.Muon_PDG==13) &&
        (mytv.HMProton_PDG==2212) );
    };

  auto QESel_DP = [&](TransversityVarsB const &mytv) -> bool {
      return (
        QESel(mytv) &&
        (mytv.Deltap_HMProton_MeV.Vect().Mag() > DeltaPCut) );
    };

  auto RESSel = [&](TransversityVarsB const &mytv) -> bool {
      return (
        (mytv.NeutConventionReactionCode==11) &&
        (mytv.Muon_PDG==13) &&
        (mytv.HMProton_PDG==2212) &&
        (mytv.HMCPion_PDG==211) );
    };

  auto RESSel_DP = [&](TransversityVarsB const &mytv) -> bool {
      return (
        RESSel(mytv) &&
        (mytv.Deltap_HMProtonPion_MeV.Vect().Mag() > DeltaPCut) );
    };

  Double_t N_QE = 0;
  Double_t N_QE_DP = 0;

  Double_t N_RES = 0;
  Double_t N_RES_DP = 0;

  Long64_t rtn = RunLoop(argc2,argv2,
    [&](Long64_t ent, TransversityVarsB const &tv){
      N_QE += QESel(tv);
      N_QE_DP += QESel_DP(tv);
      N_RES += RESSel(tv);
      N_RES_DP += RESSel_DP(tv);
      return true;
    });
  if(rtn==-1){return 1;}

  std::cout << "Taus:" << std::endl;
  std::cout << "\tQE: " << (N_QE_DP/N_QE) << " = (" << N_QE_DP << "/"
    << N_QE << ")" << std::endl;
  std::cout << "\tRES: " << (N_RES_DP/N_RES) << " = (" << N_RES_DP << "/"
    << N_RES << ")" << std::endl;

}
