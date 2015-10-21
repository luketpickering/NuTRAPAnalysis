#include <iostream>

#include <assert.h>

#include <TVector3.h>

#include "TransversityUtils.hxx"

using namespace TransversityUtils;

void GetVectorInTPlane_Tests(){
/*******************************GetVectorInTPlane*******************************/
  TVector3 norm1(0,0,1);
  TVector3 norm2(1,1,1);

  TVector3 TestV1(1,1,1);
  TVector3 TestV2(0,0,1);
  TVector3 TestV3(1,0,0);
  TVector3 TestV4(0,1,0);
  TVector3 TestV5(-1,0,0);
  TVector3 TestV6(0,-1,0);

  TVector3 TestV1T = GetVectorInTPlane(TestV1,norm1);
  TVector3 TestV2T = GetVectorInTPlane(TestV2,norm1);
  TVector3 TestV3T = GetVectorInTPlane(TestV3,norm1);
  TVector3 TestV4T = GetVectorInTPlane(TestV4,norm1);
  TVector3 TestV5T = GetVectorInTPlane(TestV5,norm1);
  TVector3 TestV6T = GetVectorInTPlane(TestV6,norm1);

  assert(TestV1T == TVector3(1,1,0));
  std::cout << "[PASSED] GetVectorInTPlane test #1" << std::endl;
  assert(TestV2T == TVector3(0,0,0));
  std::cout << "[PASSED] GetVectorInTPlane test #2" << std::endl;
  assert(TestV3T == TVector3(1,0,0));
  std::cout << "[PASSED] GetVectorInTPlane test #3" << std::endl;
  assert(TestV4T == TVector3(0,1,0));
  std::cout << "[PASSED] GetVectorInTPlane test #4" << std::endl;
  assert(TestV5T == TVector3(-1,0,0));
  std::cout << "[PASSED] GetVectorInTPlane test #5" << std::endl;
  assert(TestV6T == TVector3(0,-1,0));
  std::cout << "[PASSED] GetVectorInTPlane test #6" << std::endl;


  std::cout <<
"*******************************************************************************"
"**\n"
"***********************PASSED GetVectorInTPlane tests**************************"
"**\n"
"*******************************************************************************"
"**"  << std::endl;

/******************************************************************************/
}

void GetDeltaPhiT_Tests(){
/***********************************GetDeltaPhiT*******************************/
  TVector3 norm1_dphi_t(0,0,1);

  TVector3 TestV1_dphi_t(1,0,0);
  TVector3 TestV2_dphi_t(0,1,0);
  TVector3 TestV3_dphi_t(-1,0,0);
  TVector3 TestV4_dphi_t(0,-1,0);
  //TVector3 TestV5_dphi_t(1,1,0);

  Double_t Res1 = GetDeltaPhiT(TestV1_dphi_t,TestV1_dphi_t,norm1_dphi_t);
  assert((Res1 - M_PI) < 1E-15);
  std::cout << "[PASSED] GetDeltaPhiT test #1" << std::endl;
  Double_t Res2 = GetDeltaPhiT(TestV1_dphi_t,TestV2_dphi_t,norm1_dphi_t);
  assert((Res2 - (M_PI/2.0)) < 1E-15);
  std::cout << "[PASSED] GetDeltaPhiT test #1" << std::endl;
  Double_t Res3 = GetDeltaPhiT(TestV1_dphi_t,TestV3_dphi_t,norm1_dphi_t);
  assert((Res3 - 0) < 1E-15);
  std::cout << "[PASSED] GetDeltaPhiT test #1" << std::endl;
  Double_t Res4 = GetDeltaPhiT(TestV1_dphi_t,TestV4_dphi_t,norm1_dphi_t);
  assert((Res4 - (M_PI/2.0)) < 1E-15);
  std::cout << "[PASSED] GetDeltaPhiT test #1" << std::endl;

  std::cout <<
"*******************************************************************************"
"**\n"
"**************************PASSED GetDeltaPhiT tests****************************"
"**\n"
"*******************************************************************************"
"**"  << std::endl;

}

void GetDeltaPT_Tests(){
/***********************************GetDeltaPhiT*******************************/
  TVector3 norm1_dphi_t(0,0,1);

  TVector3 TestV1_dphi_t(1,0,0);
  TVector3 TestV2_dphi_t(0,1,0);
  TVector3 TestV3_dphi_t(-1,0,0);
  TVector3 TestV4_dphi_t(0,-1,0);
  TVector3 TestV5_dphi_t(1,2,1);

  TVector3 Res1 = GetDeltaPT(TestV1_dphi_t,TestV1_dphi_t,norm1_dphi_t);
  assert(Res1 == TVector3(2,0,0));
  std::cout << "[PASSED] GetDeltaPT test #1" << std::endl;

  TVector3 Res2 = GetDeltaPT(TestV1_dphi_t,TestV2_dphi_t,norm1_dphi_t);
  assert(Res2 == TVector3(1,1,0));
  std::cout << "[PASSED] GetDeltaPT test #1" << std::endl;

  TVector3 Res3 = GetDeltaPT(TestV1_dphi_t,TestV3_dphi_t,norm1_dphi_t);
  assert(Res3 == TVector3(0,0,0));
  std::cout << "[PASSED] GetDeltaPT test #1" << std::endl;

  TVector3 Res4 = GetDeltaPT(TestV1_dphi_t,TestV4_dphi_t,norm1_dphi_t);
  assert(Res4 == TVector3(1,-1,0));
  std::cout << "[PASSED] GetDeltaPT test #1" << std::endl;

  TVector3 Res5 = GetDeltaPT(TestV1_dphi_t,TestV5_dphi_t,norm1_dphi_t);
  assert(Res5 == TVector3(2,2,0));
  std::cout << "[PASSED] GetDeltaPT test #1" << std::endl;

  TVector3 Res6 = GetDeltaPT(TestV3_dphi_t,TestV5_dphi_t,norm1_dphi_t);
  assert(Res6 == TVector3(0,2,0));
  std::cout << "[PASSED] GetDeltaPT test #1" << std::endl;

  std::cout <<
"*******************************************************************************"
"**\n"
"****************************PASSED GetDeltaPT tests****************************"
"**\n"
"*******************************************************************************"
"**" << std::endl;

}

int main(){
  std::cout <<
"*******************************************************************************"
"**\n"
"*******************************************************************************"
"**\n"
"**********************************TESTS****************************************"
"**\n"
"*******************************************************************************"
"**\n"
"*******************************************************************************"
"**"  << std::endl;
  GetVectorInTPlane_Tests();
  GetDeltaPhiT_Tests();
  GetDeltaPT_Tests();
  return 0;
}
