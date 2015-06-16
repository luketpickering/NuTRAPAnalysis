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

}
