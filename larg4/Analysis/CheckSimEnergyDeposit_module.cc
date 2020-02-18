// art Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "artg4tk/services/DetectorHolder_service.hh"
#include "lardataobj/Simulation/SimEnergyDeposit.h"

// Root includes.
#include "TH1F.h"
#include "TH2F.h"
#include "TNtuple.h"

// STL includes.
#include <cmath>

// Other includes.
#include "CLHEP/Units/SystemOfUnits.h"


typedef std::tuple<std::vector<double>, std::vector<double>, std::vector<double> > calo_t;

calo_t condense_calo(calo_t oldcalo, double step_limit=0.5){ // default 0.5cm
  std::vector<double> old_dE, old_dx, old_currLen; 
  std::tie(old_dE, old_dx, old_currLen) = oldcalo;

  std::vector<double> new_dE, new_dx, new_currLen;
  double sumEnergy=0, sumLength=0;
  for(size_t i=0; i<old_dE.size(); i++) {
    sumEnergy += old_dE.at(i);    
    sumLength += old_dx.at(i);    
    if (sumLength > step_limit) {
      new_dE.push_back(sumEnergy);
      new_dx.push_back(sumLength);
      new_currLen.push_back(old_currLen.at(i));
      sumEnergy = 0.;
      sumLength = 0.;
    }
  }

  // clear up the remaining steps shorter than the step limit
  if (sumLength>0){
      new_dE.push_back(sumEnergy);
      new_dx.push_back(sumLength);
      double current_length=sumLength;
      if (not new_currLen.empty()) 
        current_length += new_currLen.back();
      new_currLen.push_back(current_length);
      sumEnergy = 0.;
      sumLength = 0.;
  }

  return std::make_tuple(new_dE, new_dx, new_currLen);
}

using namespace std;
namespace larg4 {
    class CheckSimEnergyDeposit;
}

class larg4::CheckSimEnergyDeposit : public art::EDAnalyzer {
public:

  explicit CheckSimEnergyDeposit(fhicl::ParameterSet const& p);

private:
  void beginJob() override;
  void analyze(const art::Event& event) override;

  TH1F* _hnHits{nullptr};         // number of SimEnergyDepositHits
  TH1F* _hEdep{nullptr};          // average energy deposition in SimEnergyDepositHits
  TH1F* _hnumPhotons{nullptr};    // number of Photons per SimEnergyDepositHits
  TH1F* _hLandauPhotons{nullptr}; // Edep/cm  SimEnergyDepositHits
  TH1F* _hLandauEdep{nullptr};    // number of Photons/cm SimEnergyDepositHits
  TH1F* _hSteplength{nullptr};    // Geant 4 step length
  TH2F* _hdEdxRR{nullptr};    // dE/dx vs residual range
  TH2F* _hSteplengthRR{nullptr};    // step length vs residual range
  TNtuple* _ntuple{nullptr};
};

larg4::CheckSimEnergyDeposit::CheckSimEnergyDeposit(fhicl::ParameterSet const& p) :
  art::EDAnalyzer(p)
{}

void larg4::CheckSimEnergyDeposit::beginJob()
{
  art::ServiceHandle<art::TFileService const> tfs;
  _hnHits = tfs->make<TH1F>("hnHits", "Number of SimEnergyDeposits", 300, 0, 0);
  _hEdep = tfs->make<TH1F>("hEdep", "Energy deposition in SimEnergyDeposits", 100,0.,0.02);
  _hnumPhotons = tfs->make<TH1F>("hnumPhotons", "number of photons per  SimEnergyDeposit", 100,0.,500.);
  _hLandauPhotons= tfs->make<TH1F>("hLandauPhotons", "number of photons/cm", 100,0.,2000000.);
  _hLandauEdep= tfs->make<TH1F>("hLandauEdep", "Edep/cm", 100,0.,10.);
  _hSteplength= tfs->make<TH1F>("hSteplength", "geant 4 step length", 200,0.,1.0);
  _ntuple = tfs->make<TNtuple>("ntuple","Demo ntuple",
                               "Event:Edep:em_Edep:nonem_Edep:xpos:ypos:zpos:time");
  _hdEdxRR= tfs->make<TH2F>("hdEdxRR", "dEdx vs. RR", 240,0,120, 300,0,15);
  _hSteplengthRR= tfs->make<TH2F>("hSteplengthRR", "StepLength vs. RR", 240,0,120, 200,0,1.0);
} // end beginJob


void larg4::CheckSimEnergyDeposit::analyze(const art::Event& event)
{
  std::vector<art::Handle<sim::SimEnergyDepositCollection>> allSims;
  event.getManyByType(allSims);

  bool has_reinteraction_proton = false;
  for (auto const& sims : allSims) {
    // double sumPhotons=0.0;
    // double sumE = 0.0;
    // _hnHits->Fill(sims->size());
    // for (auto const& hit : *sims) {
    //   // sum up energy deposit in a 1cm slice of liquid Argon.
    //   if (std::abs(hit.EndZ())<0.5) {
    //     sumPhotons= sumPhotons + hit.NumPhotons();
    //     sumE= sumE +hit.Energy();
    //   }
    //   _hnumPhotons->Fill( hit.NumPhotons());
    //   _hEdep->Fill( hit.Energy());   // energy deposit in MeV
    //   _hSteplength->Fill( hit.StepLength()); // step length in cm
    //   /*
    //     _ntuple->Fill(event.event(),
    //     hit.GetEdep(),
    //     hit.GetEdepEM(),
    //     hit.GetEdepnonEM(),
    //     hit.GetXpos(),
    //     hit.GetYpos(),
    //     hit.GetZpos(),
    //     hit.GetTime());
    //   */
    // }
    // _hLandauPhotons->Fill(sumPhotons);
    // _hLandauEdep->Fill(sumE);

    vector<double> dE, dx, currLen; // dE, dx and current length for a step
    double sumLength(0);
    for (auto const& hit : *sims) {
      int trkId = hit.TrackID();
      int pdg = hit.PdgCode();
      if (pdg!=2212) continue; // only care about proton
      if (trkId>1) {
        has_reinteraction_proton = true;
        break; // ignore the event if there are multiple protons
      }

      double e = hit.Energy();
      double l = hit.StepLength();
      _hSteplength->Fill(l);

      sumLength += l;
      dE.push_back(e);
      dx.push_back(l);
      currLen.push_back(sumLength);
    } // end a simEnergyDeposit

    if (has_reinteraction_proton) continue;

    double totalLength = sumLength;
    auto [new_dE, new_dx, new_currLen] = condense_calo(std::make_tuple(dE,dx,currLen));

    for(size_t i=0; i<new_dE.size(); i++){
      double res_range = totalLength - new_currLen.at(i);
      if(res_range>120) std::cout << res_range << " " << new_dE.at(i)/new_dx.at(i) << std::endl;
      _hdEdxRR->Fill(res_range, new_dE.at(i)/new_dx.at(i));
      _hSteplengthRR->Fill(res_range, new_dx.at(i));
    }
  } // end simEnergyDeposits
} // end analyze

DEFINE_ART_MODULE(larg4::CheckSimEnergyDeposit)
