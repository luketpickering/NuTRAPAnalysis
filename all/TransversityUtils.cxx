#include "TransversityUtils.hxx"

namespace TransversityUtils {

TVector3 GetVectorInTPlane(const TVector3& inp,
  const TVector3& planarNormal){
  TVector3 pnUnit = planarNormal.Unit();
  double inpProjectPN = inp.Dot(pnUnit);

  TVector3 InPlanarInput = inp-(inpProjectPN*pnUnit);
  return InPlanarInput;
}

TVector3 GetUnitVectorInTPlane(const TVector3& inp,
  const TVector3& planarNormal){
  return GetVectorInTPlane(inp,planarNormal).Unit();
}

Double_t GetDeltaPhiT(TVector3 const &V_lepton,
  TVector3 const &V_other, TVector3 const &Normal,
  bool PiMinus){

  TVector3 V_lepton_T = GetUnitVectorInTPlane(
    V_lepton,
    Normal);

  TVector3 V_other_T = GetUnitVectorInTPlane(
    V_other,
    Normal);

  return PiMinus ? acos(V_lepton_T.Dot(V_other_T)) :
    (M_PI - acos(V_lepton_T.Dot(V_other_T))) ;
}

TVector3 GetDeltaPT(TVector3 const &V_lepton,
  TVector3 const &V_other, TVector3 const &Normal){

  TVector3 V_lepton_T = GetVectorInTPlane(
    V_lepton,
    Normal);

  TVector3 V_other_T = GetVectorInTPlane(
    V_other,
    Normal);

  return V_lepton_T + V_other_T;
}

Double_t GetDeltaAlphaT(TVector3 const &V_lepton,
  TVector3 const &V_other, TVector3 const &Normal, bool PiMinus){

  TVector3 DeltaPT = GetDeltaPT(V_lepton,V_other,Normal);

  return GetDeltaPhiT(V_lepton,DeltaPT,Normal,PiMinus);
}

Double_t GetDeltaPTT(TVector3 const &V_lepton,
  TVector3 const &V_other1, TVector3 const &V_other2, TVector3 const &Normal){
  TVector3 ZTT_Unit = (Normal.Cross(V_lepton)).Unit();

  return V_other1.Dot(ZTT_Unit) + V_other2.Dot(ZTT_Unit);
}

Double_t GetReconNuEnergy(const TVector3 &nudir,
  const TLorentzVector &mumom, const TLorentzVector &pmom){

  TVector3 m3mom = mumom.Vect();
  TVector3 p3mom = pmom.Vect();

  Double_t p_tot_p = p3mom.Mag();
  Double_t p_tot_m = m3mom.Mag();
  Double_t nudirMag = nudir.Mag();

  Double_t mnangle =
    TMath::ACos((nudir.Dot(m3mom))/(p_tot_m*nudirMag));

  Double_t EnuRec =
    p_tot_m*TMath::Cos(mnangle)
      + sqrt(p_tot_m*TMath::Cos(mnangle)*p_tot_m*TMath::Cos(mnangle)
        - p_tot_m*p_tot_m + p_tot_p*p_tot_p);

  return EnuRec;
}

Double_t GetCCQEERec(const TVector3 &nudir,
  const TLorentzVector &mumom, Double_t E_b){

  constexpr static Double_t m_p = 938.0; //[MeV]
  constexpr static Double_t m_p2 = m_p*m_p; //[MeV^2]
  constexpr static Double_t m_n = 940.0; //[MeV]
  constexpr static Double_t m_mu = 106.0; //[MeV]
  constexpr static Double_t m_mu2 = m_mu*m_mu; //[MeV^2]

  Double_t m_n_red = m_n - E_b;
  Double_t e_mu = mumom.E();
  Double_t mom_mu = mumom.Vect().Mag();
  Double_t mu_ctheta = nudir.Unit().Dot(mumom.Vect().Unit());

  Double_t EnuRec =
    (m_p2 - (m_n_red*m_n_red) - m_mu2 + 2.0*m_n_red*e_mu )/
    //------------------------------------------------------
          (2.0*(m_n_red - e_mu + mom_mu*mu_ctheta));

  return EnuRec;
}

Double_t GetMuPERec(const TLorentzVector &mumom, const TLorentzVector &pmom,
  Double_t E_b){

//Enu = Emu + Ep - m_n + 25 MeV (where 25 MeV) is the binding energy
  constexpr static Double_t m_n = 940.0; //[MeV]
  Double_t m_n_red = m_n - E_b;

  Double_t EnuRec = mumom.E() + pmom.E() - m_n_red;
  return EnuRec;
}

Double_t GetReconTgtMass(Double_t const &nuE,
  const TLorentzVector &mumom, const TLorentzVector &pmom){
  Double_t TgtMass = pmom.E() + mumom.E() - nuE;

  return TgtMass;
}

}
