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

Double_t GetDeltaPhiT(TVector3 const &V1,
  TVector3 const &V2, TVector3 const &Normal,
  bool PiMinus){

  TVector3 V1_T = GetUnitVectorInTPlane(
    V1,
    Normal);

  TVector3 V2_T = GetUnitVectorInTPlane(
    V2,
    Normal);

  return PiMinus ? acos(V1_T.Dot(V2_T)) : (M_PI - acos(V1_T.Dot(V2_T))) ;
}

}
