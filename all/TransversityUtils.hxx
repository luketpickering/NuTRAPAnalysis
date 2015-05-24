#ifndef __TRANSVERSITY_UTILS_SEEN__
#define __TRANSVERSITY_UTILS_SEEN__
#include <cstdlib>
#include <cerrno>
#include <climits>
#include <cmath>

#include <sstream>

#include <TVector3.h>

namespace TransversityUtils {

enum Generators {kNEUT,kGENIE,kNuWro,kGiBUU,kInvalid};

TVector3 GetVectorInTPlane(const TVector3& inp,
  const TVector3& planarNormal);

TVector3 GetUnitVectorInTPlane(const TVector3& inp,
  const TVector3& planarNormal);

//Get delta phi between V1 and V2 in the plane defined by Normal.
//PiMinus changes the convention between back-to-back is delta phi t == 0
// and delta phi t == pi
Double_t GetDeltaPhiT(TVector3 const &V_lepton,
  TVector3 const &V_other, TVector3 const &Normal,
  bool PiMinus=false);
TVector3 GetDeltaPT(TVector3 const &V_lepton,
  TVector3 const &V_other, TVector3 const &Normal);
Double_t GetDeltaAlphaT(TVector3 const &V_lepton,
  TVector3 const &V_other, TVector3 const &Normal, bool PiMinus=false);

} // namespace Transversityutils
#endif
