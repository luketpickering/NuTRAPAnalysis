#include <cstdlib>
#include <cerrno>
#include <climits>

#include <algorithm>

#include <sstream>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include "objects.hxx"

const int kGStdHepNPmax = 350;
const int kGStdHepIdxPx = 0;
const int kGStdHepIdxPy = 1;
const int kGStdHepIdxPz = 2;
const int kGStdHepIdxE  = 3;

TVector3 GetVectorInPlane(const TVector3& inp,
  const TVector3& planarNormal){
  TVector3 pnUnit = planarNormal.Unit();
  double inpProjectPN = inp.Dot(pnUnit);

  TVector3 InPlanarInput = inp-(inpProjectPN*pnUnit);
  return InPlanarInput;
}

TVector3 GetUnitVectorInPlane(const TVector3& inp,
  const TVector3& planarNormal){
  return GetVectorInPlane(inp,planarNormal).Unit();
}


  enum STR2INT_ERROR { STRINT_SUCCESS,
                       STRINT_OVERFLOW,
                       STRINT_UNDERFLOW,
                       STRINT_INCONVERTIBLE };

  ///Converts a string to a long, checking for errors.
  ///See STR2INT_ERROR for error codes.
  STR2INT_ERROR str2int (long &i, char const *s, int base = 0) {
    char *end;
    long  l;
    errno = 0;
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
  STR2INT_ERROR str2int (int &i, char const *s, int base = 0) {
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

namespace {
  int StdHepN;
  int StdHepPdg[kGStdHepNPmax];
  double StdHepP4[kGStdHepNPmax][4];
  int StdHepStatus[kGStdHepNPmax];
  int G2NeutEvtCode;
  int VERBOSE = 0;
}

template<typename T, size_t N>
std::string PrintArray(const T (&arr)[N]){
  std::stringstream ss("");
  ss << "[ ";
  for(size_t i = 0; i < N; ++i){
    ss << arr[i] << (((i+1)==N)?"":", ");
  }
  ss << " ]";
  return ss.str();
}

void genieConv(const char* inpName,
  const char* outName="gn_PureSim_transversity.root",
  long long max=-1){

  TChain* gRooTrackerChain = new TChain("gRooTracker");
  gRooTrackerChain->Add(inpName);
  gRooTrackerChain->SetBranchAddress("StdHepN", &StdHepN);
  gRooTrackerChain->SetBranchAddress("StdHepPdg", &StdHepPdg);
  gRooTrackerChain->SetBranchAddress("StdHepP4", &StdHepP4);
  gRooTrackerChain->SetBranchAddress("G2NeutEvtCode", &G2NeutEvtCode);
  gRooTrackerChain->SetBranchAddress("StdHepStatus", &StdHepStatus);

  TFile* outFile = new TFile(outName,"RECREATE");
  TTree* outTreePureSim = new TTree("TransversitudenessPureSim","");
  TransversityInfoPureSim* outInfoPureSim = 0;
  outTreePureSim->Branch("TransversityInfo",&outInfoPureSim);

  long long doEntries = (max==-1)?gRooTrackerChain->GetEntries() :
    (std::min(max,gRooTrackerChain->GetEntries()));

  for(long long i = 0; i < doEntries; ++i){
    if(!i && (VERBOSE > 0)){
      std::cout << "Getting first entry" << std::endl;
    }
    gRooTrackerChain->GetEntry(i);
    if(!i && (VERBOSE > 0)){
      std::cout << "Got first entry" << std::endl;
    }

    std::pair<int, std::pair<double,TLorentzVector> > highestMomProton(0,
      std::make_pair(0,TLorentzVector(0,0,0,0)));
    std::pair<int, std::pair<double,TLorentzVector> > highestMomOther(0,
      std::make_pair(0,TLorentzVector(0,0,0,0)));
    std::pair<int, std::pair<double,TLorentzVector> > highestMomTrackable(0,
      std::make_pair(0,TLorentzVector(0,0,0,0)));
    std::pair<int, TLorentzVector> firstProton(0,TLorentzVector(0,0,0,0));
    std::pair<int, TLorentzVector> Muon(0,TLorentzVector(0,0,0,0));

    int NProtons = 0;
    int NOtherTraj = 0;
    int NGamma = 0;
    int NNeutrons = 0;
    int NPiPlus = 0;
    int NPiZero = 0;
    int NPiMinus = 0;
    int NPrimTraj = 0;

    for(int partNum = 2; partNum < StdHepN; ++partNum){
      if(StdHepStatus[partNum] != 1){continue;}
      if(StdHepPdg[partNum]>=2000000000){continue;}
      NPrimTraj++;

      switch(StdHepPdg[partNum]){
        case 13:
          break;
        case 22:
          NGamma++;
          break;
        case 2212:
          NProtons++;
          break;
        case 2112:
          NNeutrons++;
          break;
        case 211:
          NPiPlus++;
          break;
        case -211:
          NPiMinus++;
          break;
        case 111:
          NPiZero++;
          break;
        default:
          std::cout << "Found particle: PDG == "
            << StdHepPdg[partNum] << std::endl;
          NOtherTraj++;
          break;
      }

      TLorentzVector momTLV = TLorentzVector(
        StdHepP4[partNum][kGStdHepIdxPx],
        StdHepP4[partNum][kGStdHepIdxPy],
        StdHepP4[partNum][kGStdHepIdxPz],
        StdHepP4[partNum][kGStdHepIdxE]);
      if(StdHepPdg[partNum] == 13){
        Muon.second = momTLV;
        Muon.first = StdHepPdg[partNum];
      }
      else if(StdHepPdg[partNum] == 2212){

        if(!firstProton.first){
          firstProton.second = momTLV;
          firstProton.first = StdHepPdg[partNum];
          highestMomProton.first = StdHepPdg[partNum];
          highestMomProton.second = std::make_pair(
            momTLV.Vect().Mag2(),momTLV);
          if(VERBOSE>0){
            std::cout << "First Proton("<<i<<"): " << highestMomProton.second.second.Vect().Mag() << std::endl;
          }
        }

        if(highestMomProton.second.first < momTLV.Vect().Mag2()){
          if(VERBOSE>0){
            std::cout << "Upping highestMomProton("<<i<<"): " << highestMomProton.second.second.Vect().Mag() << std::endl;
          }
          highestMomProton.second = std::make_pair(
            momTLV.Vect().Mag2(),momTLV);
        }
      } else {

        if(!highestMomOther.first){
          highestMomOther.first = StdHepPdg[partNum];
          highestMomOther.second = std::make_pair(
            momTLV.Vect().Mag2(),momTLV);
        }

        if(highestMomOther.second.first < momTLV.Vect().Mag2()){
          highestMomOther.first = StdHepPdg[partNum];
          highestMomOther.second = std::make_pair(
            momTLV.Vect().Mag2(),momTLV);
        }
      }
      if( (StdHepPdg[partNum] == 2212) ||
          (StdHepPdg[partNum] == 211) ||
          (StdHepPdg[partNum] == -211) ){

        if(highestMomTrackable.second.first < momTLV.Vect().Mag2()){
          highestMomTrackable.first = StdHepPdg[partNum];
          highestMomTrackable.second = std::make_pair(
            momTLV.Vect().Mag2(),momTLV);
        }
      }
    }

    outInfoPureSim->IsValid = bool(Muon.first) && bool(firstProton.first);
    outInfoPureSim->MuonDirn = Muon.second.Vect().Unit();
    outInfoPureSim->HighestMomProtonDirn = highestMomProton.second.second.Vect().Unit();
    outInfoPureSim->MuonPDG = Muon.first;
    outInfoPureSim->HighestMomProtonPDG = highestMomProton.first;
    outInfoPureSim->HighestMomOtherPDG = highestMomOther.first;
    outInfoPureSim->ProtonMom_HighestMomProton = highestMomProton.second.second.Vect().Mag();
    outInfoPureSim->ProtonMom_FirstProton = firstProton.second.Vect().Mag();

    if((VERBOSE>1) ||
       ( (G2NeutEvtCode==1) &&
         ( ((VERBOSE>0) && (NProtons==0)) ||
           ((VERBOSE>0) && (outInfoPureSim->ProtonMom_HighestMomProton < 1E-8))
         )
       )
      ){
      std::cout << "%%Reason("<<i<<") = (NProtons==0): " << (NProtons)
        << " (outInfoPureSim->ProtonMom_HighestMomProton < 1E-8): "
        << outInfoPureSim->ProtonMom_HighestMomProton << std::endl;
      std::cout << "NParticles: " << StdHepN << " - (" << G2NeutEvtCode << ")" << std::endl;
      for(int partNum = 0; partNum < StdHepN; ++partNum){
        std::cout << "\t" << partNum << ": " << StdHepPdg[partNum]
          << " (Status==" << StdHepStatus[partNum] << ") "
          << PrintArray(StdHepP4[partNum]) << std::endl;
      }
    }

    outInfoPureSim->deltaPhiTTrue = 0;

    outInfoPureSim->VtxLocation = TLorentzVector(0,0,0,0);
    outInfoPureSim->IncNeutrinoPDG = StdHepPdg[0];
    outInfoPureSim->IncNeutrinoMmtm = TLorentzVector(
      StdHepP4[0][kGStdHepIdxPx],
      StdHepP4[0][kGStdHepIdxPy],
      StdHepP4[0][kGStdHepIdxPz],
      StdHepP4[0][kGStdHepIdxE]);
    outInfoPureSim->TargetPDG = StdHepPdg[1];
    outInfoPureSim->TargetZ = ((outInfoPureSim->TargetPDG/10000)%1000);

    outInfoPureSim->NPrimTraj = NPrimTraj;

    outInfoPureSim->NProtons = NProtons;
    outInfoPureSim->NNeutrons = NNeutrons;
    outInfoPureSim->NGamma = NGamma;
    outInfoPureSim->NPiPlus = NPiPlus;
    outInfoPureSim->NPiMinus = NPiMinus;
    outInfoPureSim->NPiZero = NPiZero;
    outInfoPureSim->NOtherTraj = NOtherTraj;

    // outInfoPureSim->ReacCode = EvtCode;
    // str2int(outInfoPureSim->NeutReacCode,EvtCode.Data());
    outInfoPureSim->ReacCode = "";
    outInfoPureSim->NeutReacCode = G2NeutEvtCode;

    outInfoPureSim->NUnsavedTrajectories = 0;
    outInfoPureSim->ModBit = 0;

    TVector3 PlanarUnitMuonTrue = GetUnitVectorInPlane(
      outInfoPureSim->MuonDirn,
      outInfoPureSim->IncNeutrinoMmtm.Vect());

    TVector3 PlanarUnitProtonTrue_high = GetUnitVectorInPlane(
      outInfoPureSim->HighestMomProtonDirn,
      outInfoPureSim->IncNeutrinoMmtm.Vect());

    outInfoPureSim->DeltaPhiT_HighestMomProton =
      (M_PI - acos(PlanarUnitMuonTrue.Dot(PlanarUnitProtonTrue_high)));

    TVector3 PlanarUnitProtonTrue_first = GetUnitVectorInPlane(
      firstProton.second.Vect().Unit(),
      outInfoPureSim->IncNeutrinoMmtm.Vect());

    outInfoPureSim->DeltaPhiT_FirstProton =
      (M_PI - acos(PlanarUnitMuonTrue.Dot(PlanarUnitProtonTrue_first)));

    TVector3 PlanarUnitProtonTrue_trackable = GetUnitVectorInPlane(
      highestMomTrackable.second.second.Vect().Unit(),
      outInfoPureSim->IncNeutrinoMmtm.Vect());

    outInfoPureSim->DeltaPhiT_HighestMomTrackable =
      (M_PI - acos(PlanarUnitMuonTrue.Dot(PlanarUnitProtonTrue_trackable)));

    outInfoPureSim->HighestMomTrackablePDG = highestMomTrackable.first;

    outInfoPureSim->MuonPt = GetVectorInPlane(Muon.second.Vect(),
      outInfoPureSim->IncNeutrinoMmtm.Vect());
    outInfoPureSim->ModMuonPt = outInfoPureSim->MuonPt.Mag();

    outTreePureSim->Fill();

  }
  outTreePureSim->Write();
  outFile->Write();
}
