#include "TransversityUtils.hxx"

namespace TransversityUtils {

///Converts a string to a long, checking for errors.
///See STR2INT_ERROR for error codes.
STR2INT_ERROR str2int (long &i, char const *s, int base) {
  char *end;
  long  l;
  errno;
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
STR2INT_ERROR str2int (int &i, char const *s, int base) {
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

}
