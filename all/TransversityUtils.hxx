#ifndef __TRANSVERSITY_UTILS_SEEN__
#define __TRANSVERSITY_UTILS_SEEN__
#include <cstdlib>
#include <cerrno>
#include <climits>
#include <cmath>

#include <sstream>

#include <TVector3.h>

namespace TransversityUtils {

enum Generators {kNEUT,kGENIE,kInvalid};

//I guess that the pointer array in a statically
//allocated 2D array shouldn't change, this is probably alright.
template<typename T, size_t N, size_t M>
T** NewPPOf2DArray(T (&arr)[N][M]){
  T** DynArrOfP = new T*[N];
  for(size_t i =0; i < N; ++i){
    DynArrOfP[i] = arr[i];
  }
  return DynArrOfP;
}

template<typename T, size_t N>
std::string PrintArray(T const (&arr)[N]){
  std::stringstream ss("");
  ss << "[ ";
  for(size_t i = 0; i < N; ++i){
    ss << arr[i] << (((i+1)==N)?"":", ");
  }
  ss << " ]";
  return ss.str();
}

template<typename T>
std::string PrintArray(const T* arr, size_t N){
  std::stringstream ss("");
  ss << "[ ";
  for(size_t i = 0; i < N; ++i){
    ss << arr[i] << (((i+1)==N)?"":", ");
  }
  ss << " ]";
  return ss.str();
}

enum STR2INT_ERROR { STRINT_SUCCESS,
                     STRINT_OVERFLOW,
                     STRINT_UNDERFLOW,
                     STRINT_INCONVERTIBLE };

///Converts a string to a long, checking for errors.
///See STR2INT_ERROR for error codes.
STR2INT_ERROR str2int (long &i, char const *s, int base = 0);

///Converts a string to a int, checking for errors.
///See STR2INT_ERROR for error codes.
STR2INT_ERROR str2int (int &i, char const *s, int base = 0);

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
