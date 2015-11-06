#include <iostream>
#include <iomanip>

#include "CLITools.hxx"

#include "EventLoop.hxx"

namespace {
template<typename T>
inline std::string AddNegSpacer(T const & num){
  std::stringstream ss;
  ss.precision(3);
  ss.flags(std::ios::scientific);
  ss << std::setw(4) << num;
  return std::string((num >= 0)?" ":"")+ss.str();
}
}

int main(int argc, char const * argv[]){

  Double_t DeltaPCut = 0;
  Double_t DeltaECut = 0;

  try {
    CLIArgs::AddOpt("-P", "--DeltaP-Threshold", true,
      [&] (std::string const &opt) -> bool {
        DeltaPCut = stod(opt);
        return true;
      }, false,
      [&](){DeltaPCut = 10.0;},
      "DeltaP threshold to count as having undergone FSI.");
    CLIArgs::AddOpt("-E", "--DeltaE-Threshold", true,
      [&] (std::string const &opt) -> bool {
        DeltaECut = stod(opt);
        return true;
      }, false,
      [&](){DeltaECut = 10.0;},
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
  auto RESSel = [&](TransversityVarsB const &mytv) -> bool {
      return (
        (mytv.NeutConventionReactionCode==11) &&
        (mytv.Muon_PDG==13) &&
        (mytv.HMProton_PDG==2212) &&
        (mytv.HMCPion_PDG==211) );
    };

  std::cout.precision(3);
  std::cout.flags(std::ios::scientific);
  RunLoop(argc2,argv2,[&](Long64_t ent, TransversityVarsB const &tv){

    bool IsQE = QESel(tv);
    if(!IsQE && !RESSel(tv)){ return false; }
    if(IsQE){
      TLorentzVector elePSum = tv.IncNeutrino_4Mom_MeV
                               + tv.StruckNucleon_4Mom_MeV
                               - tv.Muon_4Mom_MeV
                               - tv.HMProton_4Mom_MeV;
      if(elePSum.Vect().Mag() < DeltaPCut || elePSum.E() < DeltaECut ){
        std::cout << "[" << std::setw(5) << ent << "] " << "QE " << std::flush;
        std::cout << " #FS: " << std::setw(2)
          << tv.NFinalStateParticles << std::flush;
        std::cout << " DP: " << AddNegSpacer(elePSum.Vect().Mag())
          << " DE: " << AddNegSpacer(elePSum.E()) << std::endl;
      }
    } else {
      TLorentzVector elePSum = tv.IncNeutrino_4Mom_MeV
                               + tv.StruckNucleon_4Mom_MeV
                               - tv.Muon_4Mom_MeV
                               - tv.HMProton_4Mom_MeV
                               - tv.HMCPion_4Mom_MeV;
      if(elePSum.Vect().Mag() < DeltaPCut || elePSum.E() < DeltaECut ){
        std::cout << "[" << std::setw(5) << ent << "] " << "RES" << std::flush;
        std::cout << " #FS: " << std::setw(2)
          << tv.NFinalStateParticles << std::flush;
        std::cout << " DP: " << AddNegSpacer(elePSum.Vect().Mag())
          << " DE: " << AddNegSpacer(elePSum.E()) << std::endl;
      }
    }

    return true;
  });

}
