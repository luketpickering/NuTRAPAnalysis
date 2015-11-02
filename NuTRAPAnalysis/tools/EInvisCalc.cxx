#include <cmath>

#include "TChain.h"
#include "TH1D.h"

#include "CLITools.hxx"

std::string InputFileDescriptor = "";


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

  Double_t N_QE = TRAPChain->GetEntries("(TransV.NeutConventionReactionCode==1)&&(TransV.Muon_PDG==13)&&(TransV.HMProton_PDG==2212)");
  Double_t N_QE_DP = TRAPChain->GetEntries("(TransV.DeltaPTotal_HMProton_MeV.Vect().Mag()>10)&&(TransV.NeutConventionReactionCode==1)&&(TransV.Muon_PDG==13)&&(TransV.HMProton_PDG==2212)");


  Double_t N_RES = TRAPChain->GetEntries("(TransV.NeutConventionReactionCode==11)&&(TransV.Muon_PDG==13)&&(TransV.HMProton_PDG==2212)&&(TransV.HMCPion_PDG==211)");
  Double_t N_RES_DP = TRAPChain->GetEntries("(TransV.DeltaPTotal_HMProtonPion_MeV.Vect().Mag()>10)&&(TransV.NeutConventionReactionCode==11)&&(TransV.Muon_PDG==13)&&(TransV.HMProton_PDG==2212)&&(TransV.HMCPion_PDG==211)");


  std::cout << "In File: " << InputFileDescriptor << std::endl;
  std::cout << "Taus:" << std::endl;
  std::cout << "\tQE: " << (N_QE_DP/N_QE) << " = (" << N_QE_DP << "/" << N_QE << ")" << std::endl;
  std::cout << "\tRES: " << (N_RES_DP/N_RES) << " = (" << N_RES_DP << "/" << N_RES << ")" << std::endl;


  TH1D* QE_NE = new TH1D("QE_NE","",100,0,100);
  TH1D* QE = new TH1D("QE","",100,0,100);

  TRAPChain->Draw("(-1.0*(TransV.DeltaPTotal_HMProton_MeV.E()+32)) >> QE_NE","(TransV.NeutConventionReactionCode==1)&&(TransV.Muon_PDG==13)&&(TransV.HMProton_PDG==2212)&&(TransV.NFinalStateParticles!=2)");
  TRAPChain->Draw("(-1.0*(TransV.DeltaPTotal_HMProton_MeV.E()+32)) >> QE","(TransV.NeutConventionReactionCode==1)&&(TransV.Muon_PDG==13)&&(TransV.HMProton_PDG==2212)");

  TH1D* QEEtaDeltaE = ToPDF(QE);

  Double_t integral_QE = 0;
  for(Int_t i = 0; i < QE->GetNbinsX()+1; ++i){
    integral_QE = QE->GetBinWidth(i) * QE->GetBinCenter(i) * QEEtaDeltaE->GetBinContent(i) * (1.0 - (QE_NE->GetBinContent(i)/QE->GetBinContent(i)) );
  }
  integral_QE *= (N_QE_DP/N_QE);

  std::cout << "<E_invis^QE> = " << integral_QE << " MeV" << std::endl;


  TH1D* RES_NE = new TH1D("RES_NE","",100,0,100);
  TH1D* RES = new TH1D("RES","",100,0,100);

  TRAPChain->Draw("(-1.0*(TransV.DeltaPTotal_HMProton_MeV.E()+32)) >> RES_NE","(TransV.NeutConventionReactionCode==1)&&(TransV.Muon_PDG==13)&&(TransV.HMProton_PDG==2212)&&(TransV.NFinalStateParticles!=2)");
  TRAPChain->Draw("(-1.0*(TransV.DeltaPTotal_HMProton_MeV.E()+32)) >> RES","(TransV.NeutConventionReactionCode==1)&&(TransV.Muon_PDG==13)&&(TransV.HMProton_PDG==2212)");

  TH1D* RESEtaDeltaE = ToPDF(RES);

  Double_t integral_RES = 0;
  for(Int_t i = 0; i < RES->GetNbinsX()+1; ++i){
    integral_RES = RES->GetBinWidth(i) * RES->GetBinCenter(i) * RESEtaDeltaE->GetBinContent(i) * (1.0 - (RES_NE->GetBinContent(i)/RES->GetBinContent(i)) );
  }
  integral_RES *= (N_RES_DP/N_RES);

  std::cout << "<E_invis^RES> = " << integral_RES << " MeV" << std::endl;

  delete TRAPChain;
  delete QE_NE;
  delete QE;
  delete QEEtaDeltaE;
  delete RES_NE;
  delete RES;
  delete RESEtaDeltaE;
  return 0;
}
