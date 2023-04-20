#include <iostream>
#include <sstream>
#include <string>

//Functions for Isolation and BTagging
#include "Isolation.h"
#include "Pdg_check.h"
#include "Btag.h"

//ROOT Classes
#include "TFile.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TMath.h"

//HepMC2 Classes
#include "HepMC/IO_GenEvent.h"
#include "HepMC/GenEvent.h"

//FastJet3 Classes
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"

using namespace std;
using namespace HepMC;
using namespace fastjet;

int main(int argc, char* argv[]) {

if(argc != 3) {
  cout<<"Usage: main [input_file] [output_file]"<<endl;
  cout<<"input_file -- HepMC2 format"<<endl;
  cout<<"output_file -- ROOT format"<<endl;
  return 0;
}

IO_GenEvent ascii_in(argv[1], std::ios::in);

Int_t nparticles;

//Basic observables
Int_t nleptons;
Int_t njets;
Int_t nbjets;
Int_t ncjets;
Int_t nnus;

Float_t ptj[999];
Float_t ptl[999];
Float_t ptnu[999];

Float_t enerj[999];
Float_t enerl[999];
Float_t enernu[999];

Bool_t btag[999];
Bool_t ctag[999];

Int_t tagged[999];

//Define the output file and tree
TFile *outputFile;
outputFile = TFile::Open(argv[2], "RECREATE");
TTree *tree = new TTree("OutputTree", "ROOT output file");
tree -> Branch("nleptons", &nleptons, "nleptons/I");
tree -> Branch("njets", &njets, "njets/I");
tree -> Branch("nbjets", &nbjets, "nbjets/I");
tree -> Branch("ncjets", &nbjets, "nbjets/I");
tree -> Branch("nnus", &nnus, "nnus/I");
tree -> Branch("ptj1", &ptj[0], "ptj1/F");
tree -> Branch("ptj2", &ptj[1], "ptj2/F");
tree -> Branch("ptl1", &ptl[0], "ptl1/F");
tree -> Branch("ptl2", &ptl[1], "ptl2/F");
tree -> Branch("ptnu1", &ptnu[0], "ptnu1/F");
tree -> Branch("ptnu2", &ptnu[1], "ptnu2/F");
tree -> Branch("enerj1", &enerj[0], "enerj1/F");
tree -> Branch("enerj2", &enerj[1], "enerj2/F");
tree -> Branch("enerl1", &enerl[0], "enerl1/F");
tree -> Branch("enerl2", &enerl[1], "enerl2/F");
tree -> Branch("enernu1", &enernu[0], "enernu1/F");
tree -> Branch("enernu2", &enernu[1], "enernu2/F");
tree -> Branch("tagged1", &tagged[0], "tagged1/I");
tree -> Branch("tagged2", &tagged[1], "tagged2/I");

//Print the Input and Output file names
cout<<"Input_file: "<<argv[1]<<endl;
cout<<"Output_file: "<<argv[2]<<endl;

int i, j;
int entries;
int nEvents = 0;

//Define the b-tagging rate and missing rate of c and light quarks
double b_rate = 0.70;
double c_rate = 0.10;
double q_rate = 0.01;

//Define the R and p parameters for jet algorithm
double RParameter = 0.1;
int pParameter = 1;

GenEvent* evt = ascii_in.read_next_event();

TLorentzVector pvec[999];

TLorentzVector pj[999];
TLorentzVector pl[999];
TLorentzVector pnu[999];

TLorentzVector plcan;
TLorentzVector temp;

while(evt) {//Start reading events
  nEvents++;
  vector<PseudoJet> particles;

  //Initialization
  nparticles = 0; 
  njets = 0; 
  nbjets = 0; 
  ncjets = 0;
  nleptons = 0;
  nnus = 0;
  plcan.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
  temp.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
  for(i=0; i<999; i++) {
     pvec[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
     pj[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
     pl[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
     pnu[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
     ptj[i] = 0.0;
     ptl[i] = 0.0;
     ptnu[i] = 0.0;
     enerj[i] = 0.0;
     enerl[i] = 0.0;
     enernu[i] = 0.0;
     btag[i] = false;
     ctag[i] = false;
     tagged[i] = 0;
  }

  //Loop all particles
  for(GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); ++p) {
    if((*p)->status()==1) {
      pvec[nparticles].SetPxPyPzE((*p)->momentum().px(), (*p)->momentum().py(), (*p)->momentum().pz(), (*p)->momentum().e());
      nparticles++;
    }
  }
  //End of looping all particles

  //Tag objects originated from b/c hadrons
  std::set<int> bobjects = Find_all_objects(5, evt);
  std::set<int> cobjects = Find_all_objects(4, evt);

  class BC4FJ: public PseudoJet::UserInfoBase {
    public:
      BC4FJ(int barcode): _barcode(barcode) {};
      int _barcode;
  };

  //Find isolated lepton and neutrino
  for(GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); ++p) {
    if((*p)->status()==1) {
      //Find isolated leptons
      if((TMath::Abs((*p)->pdg_id())==11 || TMath::Abs((*p)->pdg_id())==13 || TMath::Abs((*p)->pdg_id())==15) && (*p)->momentum().perp()>20.0) { 
        plcan.SetPxPyPzE((*p)->momentum().px(), (*p)->momentum().py(), (*p)->momentum().pz(), (*p)->momentum().e());
        if(IsIsolated(plcan, pvec, nparticles)) {
          pl[nleptons] = plcan;
          ptl[nleptons] = pl[nleptons].Pt();
          enerl[nleptons] = pl[nleptons].E();
          nleptons++;
        }
        else {
          //If a lepton is not isolated, add it to the particles list for running fastjet
          particles.push_back(fastjet::PseudoJet((*p)->momentum().px(), (*p)->momentum().py(), (*p)->momentum().pz(), (*p)->momentum().e()));  
          //Record the barcode
          BC4FJ* info = new BC4FJ((*p)->barcode());
          particles.back().set_user_info(info);
        }
      }
      //Find neutrino
      else if((TMath::Abs((*p)->pdg_id())==12 || TMath::Abs((*p)->pdg_id())==14 || TMath::Abs((*p)->pdg_id())==16)) { 
        pnu[nnus].SetPxPyPzE((*p)->momentum().px(), (*p)->momentum().py(), (*p)->momentum().pz(), (*p)->momentum().e());
        ptnu[nnus] = pnu[nnus].Pt();
        enernu[nnus] = pnu[nnus].E();
        nnus++;
      }
      else {
        //Remained particles are added to the list to run fastjet
        particles.push_back(fastjet::PseudoJet((*p)->momentum().px(), (*p)->momentum().py(), (*p)->momentum().pz(), (*p)->momentum().e())); 
        //Record the barcode
        BC4FJ* info = new BC4FJ((*p)->barcode());
        particles.back().set_user_info(info);
      }
    }
  }
  //End of find isolated lepton and neutrion

  //Define the jet algorithm
  JetDefinition jet_def(ee_genkt_algorithm, RParameter, pParameter);
  ClusterSequence cs_ee(particles, jet_def);
  vector<PseudoJet> jets_ee = sorted_by_E(cs_ee.inclusive_jets());

  entries = jets_ee.size();

  //Jets Observables
  for(i=0; i<entries; i++) {
      vector<PseudoJet> constituents = jets_ee[i].constituents();
      bool is_b_jet = false;
      bool is_c_jet = false;
      for (j=0; j<constituents.size(); j++) {
        assert(constituents[j].has_user_info());
        int barcode = dynamic_cast<const BC4FJ*>(constituents[j].user_info_ptr())->_barcode;
        if (bobjects.find(barcode) != bobjects.end()) {
           is_b_jet = true;
           is_c_jet = false;
           break;
        } 
        if (cobjects.find(barcode) != cobjects.end()) {
           is_b_jet = false;
           is_c_jet = true;
        }
      }
      if (is_b_jet) {
        btag[i] = true; 
        ctag[i] = false;
        nbjets++;
      }
      if (is_c_jet) { 
        btag[i] = false;
        ctag[i] = true;
        ncjets++; 
      }

      pj[njets].SetPxPyPzE(jets_ee[i].px(), jets_ee[i].py(), jets_ee[i].pz(), jets_ee[i].e());
      ptj[njets] = pj[i].Pt();
      enerj[njets] = pj[i].E();
      if(btag[njets] && Btag(b_rate)) tagged[njets] = 1;
      else if(ctag[njets] && Btag(c_rate)) tagged[njets] = 1;
      else if(Btag(q_rate)) tagged[njets] = 1;
      else tagged[njets] = 0;
      njets++;
  }

  //Leptons Observables
  if(nleptons>1) {
    for(i=0; i<nleptons-1; i++) {
      for(j=0; j<nleptons-i; j++) {
        if(ptl[j] < ptl[j+1]) {
          temp = pl[j];
          pl[j] = pl[j+1];
          pl[j+1] = temp;
        }
      }
    }
    for(i=0; i<nleptons; i++) {
      ptl[i] = pl[i].Pt();
      enerl[i] = pl[i].E();
    }
  }

  //Neutrinos Observables
  if(nnus>1) {
    for(i=0; i<nnus-1; i++) {
      for(j=0; j<nnus-i; j++) {
        if(ptnu[j] < ptnu[j+1]) {
          temp = pnu[j];
          pnu[j] = pnu[j+1];
          pnu[j+1] = temp;
        }
      }
    }
    for(i=0; i<nnus; i++) {
      ptnu[i] = pnu[i].Pt();
      enernu[i] = pnu[i].E();
    }
  }

  tree->Fill();
  delete evt;
  ascii_in >> evt; 
}//End of reading events

cout<<"Number of Events"<<nEvents<<endl;
tree->Write();
outputFile->Close();
return 0;
}//end of main
