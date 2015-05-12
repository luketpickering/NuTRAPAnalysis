#include <cstdlib>
#include <cerrno>
#include <climits>

#include <algorithm>

#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include "objects.hxx"

#include "ntrac/NRooTrackerVtx.h"

const int kNStdHepIdxPx = 0;
const int kNStdHepIdxPy = 1;
const int kNStdHepIdxPz = 2;
const int kNStdHepIdxE = 3;

const int kNStdHepIdxX = 0;
const int kNStdHepIdxY = 1;
const int kNStdHepIdxZ = 2;
const int kNStdHepIdxT = 3;

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

void neutConv(const char* inpName,
  const char* outName="nt_PureSim_transversity.root",
  long long max=-1){


  TChain* nRooTrackerChain = new TChain("nRooTracker");
  nRooTrackerChain->Add(inpName);
  NRooTrackerVtx* rootrackervtx = new NRooTrackerVtx();
  nRooTrackerChain->SetBranchAddress("nRooTracker", &rootrackervtx);

  TFile* outFile = new TFile(outName,"RECREATE");
  TTree* outTreePureSim = new TTree("TransversitudenessPureSim","");
  TransversityInfoPureSim* outInfoPureSim = 0;
  outTreePureSim->Branch("TransversityInfo",&outInfoPureSim);

  long long doEntries = (max==-1)?nRooTrackerChain->GetEntries() :
    (std::min(max,nRooTrackerChain->GetEntries()));

  for(long long i = 0; i < doEntries; ++i){
    nRooTrackerChain->GetEntry(i);

    // std::cout << rootrackervtx->StdHepN << " Particles" << std::endl;

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

    for(int partNum = 2; partNum < rootrackervtx->StdHepN; ++partNum){
      if(rootrackervtx->StdHepStatus[partNum] != 1){
        std::cout << "Found a particle in NEUT with non FS status." << std::endl;
        continue;
      }

      switch(rootrackervtx->StdHepPdg[partNum]){
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
            << rootrackervtx->StdHepPdg[partNum] << std::endl;
          NOtherTraj++;
          break;
      }

      TLorentzVector momTLV = TLorentzVector(
        rootrackervtx->StdHepP4[partNum][kNStdHepIdxPx],
        rootrackervtx->StdHepP4[partNum][kNStdHepIdxPy],
        rootrackervtx->StdHepP4[partNum][kNStdHepIdxPz],
        rootrackervtx->StdHepP4[partNum][kNStdHepIdxE]);
      if(rootrackervtx->StdHepPdg[partNum] == 13){
        Muon.second = momTLV;
        Muon.first = rootrackervtx->StdHepPdg[partNum];
      }
      else if(rootrackervtx->StdHepPdg[partNum] == 2212){
        if(!firstProton.first){
          firstProton.second = momTLV;
          firstProton.first = rootrackervtx->StdHepPdg[partNum];
          highestMomProton.first = rootrackervtx->StdHepPdg[partNum];
          highestMomProton.second = std::make_pair(
            momTLV.Vect().Mag2(),momTLV);
        }

        if(highestMomProton.second.first < momTLV.Vect().Mag2()){
          highestMomProton.second = std::make_pair(
            momTLV.Vect().Mag2(),momTLV);
        }
      } else {
        if(!highestMomOther.first){
          highestMomOther.first = rootrackervtx->StdHepPdg[partNum];
          highestMomOther.second = std::make_pair(
            momTLV.Vect().Mag2(),momTLV);
        }

        if(highestMomOther.second.first < momTLV.Vect().Mag2()){
          highestMomOther.first = rootrackervtx->StdHepPdg[partNum];
          highestMomOther.second = std::make_pair(
            momTLV.Vect().Mag2(),momTLV);
        }
      }
      if( (rootrackervtx->StdHepPdg[partNum] == 2212) ||
          (rootrackervtx->StdHepPdg[partNum] == 211) ||
          (rootrackervtx->StdHepPdg[partNum] == -211) ){

        if(highestMomTrackable.second.first < momTLV.Vect().Mag2()){
          highestMomTrackable.first = rootrackervtx->StdHepPdg[partNum];
          highestMomTrackable.second = std::make_pair(
            momTLV.Vect().Mag2(),momTLV);
        }
      }
    }

    outInfoPureSim->IsValid = bool(Muon.first) && bool(firstProton.first);
    outInfoPureSim->MuonDirn = Muon.second.Vect().Unit();
    outInfoPureSim->HighestMomProtonDirn =
      highestMomProton.second.second.Vect().Unit();
    outInfoPureSim->MuonPDG = Muon.first;
    outInfoPureSim->HighestMomProtonPDG = highestMomProton.first;
    outInfoPureSim->HighestMomOtherPDG = highestMomOther.first;
    outInfoPureSim->ProtonMom_HighestMomProton =
      highestMomProton.second.second.Vect().Mag();
    outInfoPureSim->ProtonMom_FirstProton = firstProton.second.Vect().Mag();

    outInfoPureSim->deltaPhiTTrue = 0;

    outInfoPureSim->VtxLocation = TLorentzVector(0,0,0,0);
    outInfoPureSim->IncNeutrinoPDG = rootrackervtx->StdHepPdg[0];
    outInfoPureSim->IncNeutrinoMmtm = TLorentzVector(
      rootrackervtx->StdHepP4[0][kNStdHepIdxPx],
      rootrackervtx->StdHepP4[0][kNStdHepIdxPy],
      rootrackervtx->StdHepP4[0][kNStdHepIdxPz],
      rootrackervtx->StdHepP4[0][kNStdHepIdxE]);
    outInfoPureSim->TargetPDG = rootrackervtx->StdHepPdg[1];
    outInfoPureSim->TargetZ = ((outInfoPureSim->TargetPDG/10000)%1000);

    outInfoPureSim->NPrimTraj = (rootrackervtx->StdHepN - 2);

    outInfoPureSim->NProtons = NProtons;
    outInfoPureSim->NNeutrons = NNeutrons;
    outInfoPureSim->NGamma = NGamma;
    outInfoPureSim->NPiPlus = NPiPlus;
    outInfoPureSim->NPiMinus = NPiMinus;
    outInfoPureSim->NPiZero = NPiZero;
    outInfoPureSim->NOtherTraj = NOtherTraj;

    outInfoPureSim->ReacCode = rootrackervtx->EvtCode;

    str2int(outInfoPureSim->NeutReacCode,rootrackervtx->EvtCode.Data());

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
