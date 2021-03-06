#include <algorithm>

#include <iostream>

#include "TChain.h"
#include "TFile.h"
#include "TObjString.h"
#include "TTree.h"

#include "LUtils/CLITools.hxx"
#include "LUtils/Utils.hxx"

#include "GeneratorSpecifics.hxx"
#include "TransversityUtils.hxx"
#include "TransversityVariableObjects.hxx"

using namespace TransversityUtils;

// There'll be no name collisions in my DOJO, might be binary bloat though...
namespace {

int Verbosity = 1;
bool MultiplyByGeVToMeV;
bool LiteOutput = false;

Int_t NThresh = 0;
Int_t Threshs_MeV[kNThreshMax];

std::vector<int> ModesToSave;

Generator *GenProxy;

bool DoRewrite = false;

double ExtraWeight = 1.0;

}  // namespace

bool SetUpGeneratorDependence(std::string GeneratorName) {
  if ((GeneratorName == "NEUT") || (GeneratorName == "neut")) {
    GenProxy = new NEUT();
    std::cout << "Working on NEUT Tree: " << GenProxy->TreeName << std::endl;
  } else if ((GeneratorName == "GENIE") || (GeneratorName == "genie")) {
    GenProxy = new GENIE();
    std::cout << "Working on GENIE Tree: " << GenProxy->TreeName << std::endl;
  } else if ((GeneratorName == "NuWro") || (GeneratorName == "nuwro") ||
             (GeneratorName == "NUWRO")) {
    GenProxy = new NuWro();
    std::cout << "Working on NuWro Tree: " << GenProxy->TreeName << std::endl;
  } else if ((GeneratorName == "ENuWro") || (GeneratorName == "enuwro") ||
             (GeneratorName == "ENUWRO")) {
    GenProxy = new EmuNuWro();
    std::cout << "Working on Emulated NuWro Tree: " << GenProxy->TreeName
              << std::endl;
  } else if ((GeneratorName == "GiBUU") || (GeneratorName == "gibuu") ||
             (GeneratorName == "GIBUU")) {
    GenProxy = new GiBUU();
    std::cout << "Working on GiBUU Tree: " << GenProxy->TreeName << std::endl;
  } else {
    return false;
  }

  return true;
}

int ProcessRootrackerToTransversityVariables(
    const char *InputName,
    const char *OutputName = "pure_sim_transversity_variables.root",
    std::string GeneratorName = "NEUT", long long MaxEntries = -1) {
  if (!SetUpGeneratorDependence(GeneratorName)) {
    return 1;
  }

  TChain *RooTrackerChain = new TChain(GenProxy->TreeName.c_str());
  RooTrackerChain->SetDirectory(0);

  int nFiles = 0, nEntries = 0;
  if (!(nFiles = RooTrackerChain->Add(InputName)) ||
      !(nEntries = RooTrackerChain->GetEntries())) {
    std::cout << "[ERROR] Found no files (" << nFiles << ") or entries ("
              << nEntries << ")" << std::endl;
    return 2;
  }

  GenProxy->Init(RooTrackerChain);
  GenProxy->ExtraWeight = ExtraWeight;

  TFile *outFile = new TFile(OutputName, DoRewrite ? "RECREATE" : "CREATE");
  if (!outFile->IsOpen()) {
    std::cout << "Failed to open: " << OutputName << std::endl;
    return 8;
  }
  TTree *outTreePureSim = new TTree("TransversitudenessPureSim", "");

  GenProxy->AddOutputBranches(outTreePureSim, LiteOutput, MultiplyByGeVToMeV,
                              NThresh, Threshs_MeV);

  long long doEntries =
      (MaxEntries == -1)
          ? RooTrackerChain->GetEntries()
          : (std::min(MaxEntries, RooTrackerChain->GetEntries()));

  for (long long i = 0; i < doEntries; ++i) {
    RooTrackerChain->GetEntry(i);

    if (!(i % 10000)) {
      std::cout << "On entry: " << i << "/" << doEntries << std::endl;
    }

    GenProxy->DoEvent();
    GenProxy->Finalise();

    // Skips modes that we don't want to deal with.
    if (ModesToSave.size()) {
      bool found = false;
      for (auto const &mi : ModesToSave) {
        if (mi == GenProxy->NeutConventionReactionCode) {
          found = true;
          break;
        }
      }
      if (!found) {  // do not save this event
        continue;
      }
    }
    outTreePureSim->Fill();
  }
  outTreePureSim->Write();
  outFile->Write();

  std::cout << "Wrote " << outTreePureSim->GetEntries() << " entries to disk."
            << std::endl;

  outFile->Close();

  delete GenProxy;
  delete RooTrackerChain;
  delete outFile;
  return 0;
}

namespace {
std::string InputName;
std::string OutputName = "TransverseVars.root";
std::string GeneratorName = "NEUT";
long MaxEntries = -1;

void SetOpts() {
  CLIArgs::AddOpt("-i", "--input-file", true,
                  [&](std::string const &opt) -> bool {
                    std::cout << "\t--Reading from file descriptor : " << opt
                              << std::endl;
                    InputName = opt;
                    return true;
                  },
                  true, []() {}, "<TChain::Add descriptor>");

  CLIArgs::AddOpt("-o", "--output-file", true,
                  [&](std::string const &opt) -> bool {
                    std::cout << "\t--Writing to File: " << opt << std::endl;
                    OutputName = opt;
                    DoRewrite = true;
                    return true;
                  },
                  false, [&]() { OutputName = "TransverseVars.root"; },
                  "<File Name>{default=TransverseVars.root}");

  CLIArgs::AddOpt("-n", "--nentries", true,
                  [&](std::string const &opt) -> bool {

                    int ival = 0;
                    try {
                      ival = Utils::str2i(opt, true);
                    } catch (...) {
                      return false;
                    }
                    std::cout << "\t--Looking at, at most, " << ival
                              << " entries." << std::endl;
                    MaxEntries = ival;
                    return true;

                  },
                  false, [&]() { MaxEntries = -1; },
                  "<int>: Read all {default=-1}");

  CLIArgs::AddOpt("-v", "--verbosity", true,
                  [&](std::string const &opt) -> bool {
                    int ival = 0;
                    try {
                      ival = Utils::str2i(opt, true);
                    } catch (...) {
                      return false;
                    }
                    std::cout << "\t--Verbosity: " << ival << std::endl;
                    Verbosity = ival;
                    return true;

                  },
                  false, [&]() { Verbosity = 0; }, "<0-4>{default=0}");

  CLIArgs::AddOpt("-W", "--extra-weight", true,
                  [&](std::string const &opt) -> bool {
                    double ival = 0;
                    bool IsReciprocal = false;
                    std::string inp = opt;
                    if (inp[0] == 'i') {
                      IsReciprocal = true;
                      inp = Utils::Replace(inp, "i", "");
                    }

                    try {
                      ival = Utils::str2d(inp, true);
                      if (IsReciprocal) {
                        ival = (1.0 / ival);
                      }
                    } catch (...) {
                      return false;
                    }
                    std::cout << "\t--Assigning extra weight: " << ival
                              << std::endl;
                    ExtraWeight = ival;
                    return true;
                  },
                  false, []() {}, "[i]<Extra weight [1.0/]'W'>");

  CLIArgs::AddOpt("-M", "--MeV-mode", false,
                  [&](std::string const &opt) -> bool {
                    std::cout
                        << "\t--Multiplying output momenta and energies by 1E3."
                        << std::endl;
                    MultiplyByGeVToMeV = true;
                    return true;
                  },
                  false, [&]() { MultiplyByGeVToMeV = false; },
                  "[Multiply output momenta and energies "
                  "by 1E3.{default=false}]");

  CLIArgs::AddOpt("-g", "--generator", true,
                  [&](std::string const &opt) -> bool {
                    std::cout << "\t--Attempting to read generator: " << opt
                              << std::endl;
                    GeneratorName = opt;
                    return true;
                  },
                  false, [&]() { GeneratorName = "NEUT"; }, "{default=NEUT}");

  CLIArgs::AddOpt(
      "-m", "--EKin-Threshold", true,
      [&](std::string const &opt) -> bool {
        if (NThresh == kNThreshMax) {
          std::cout << "[ERROR]: Tried to add too many momentum thresholds."
                    << std::endl;
          return false;
        }
        int ival = 0;
        try {
          ival = Utils::str2i(opt, true);
        } catch (...) {
          return false;
        }

        std::cout << "\t--Added EKin threshold: " << ival << " MeV"
                  << std::endl;
        if (NThresh > 0 && (ival < Threshs_MeV[NThresh - 1])) {
          std::cout
              << "[ERROR]: Attempting to add EKin threshold at " << ival
              << " MeV, but the previous one at " << Threshs_MeV[NThresh - 1]
              << " is lower.\n\tThresholds must be added in ascending order."
              << std::endl;
          return false;
        }
        Threshs_MeV[NThresh] = ival;
        NThresh++;
        return true;
      },
      false, [&]() {}, "<int> [Add EKin threshold [MeV] {default=N/A}]");

  CLIArgs::AddOpt(
      "-N", "--NEUT-Modes", true,
      [&](std::string const &opt) -> bool {
        ModesToSave =
            Utils::StringVToIntV(Utils::SplitStringByDelim(opt, ","));

        if (ModesToSave.size()) {
          std::cout << "\t--Ignoring interactions except of the modes:  "
                    << std::flush;
          for (auto const &mi : ModesToSave) {
            std::cout << mi << ", " << std::flush;
          }
          std::cout << std::endl;
          return true;
        }
        return false;
      },
      false, []() {}, "<int,int,...> [NEUT modes to save output from.]");

  CLIArgs::AddOpt(
      "-L", "--Lite-Output", false,
      [&](std::string const &opt) -> bool {
        LiteOutput = true;
        std::cout << "\t--Outputting Lite format." << std::endl;
        return true;
      },
      false, []() { LiteOutput = false; },
      "[Will output in Lite mode which contains less output variables.]");
}
}

int main(int argc, char const *argv[]) {
  try {
    SetOpts();
  } catch (std::exception const &e) {
    std::cerr << "[ERROR]: " << e.what() << std::endl;
    return 1;
  }

  CLIArgs::AddArguments(argc, argv);
  if (!CLIArgs::HandleArgs()) {
    CLIArgs::SayRunLike();
    return 1;
  }

  std::cout << "Running with " << NThresh << " set momentum thresholds."
            << std::endl;

  return ProcessRootrackerToTransversityVariables(
      InputName.c_str(), OutputName.c_str(), GeneratorName.c_str(), MaxEntries);
}
