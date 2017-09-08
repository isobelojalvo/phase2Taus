// -*- C++ -*-
//
// Package:    RecoTauTag/timingTausMiniAOD
// Class:      timingTausMiniAOD
// 
/**\class timingTausMiniAOD timingTausMiniAOD.cc RecoTauTag/timingTausMiniAOD/plugins/timingTausMiniAOD.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Isabel Ojalvo
//         Created:  Tue, 15 Nov 2016 16:00:32 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "iostream"
#include "TMath.h"
#include "TLorentzVector.h"
#include <iostream>

#include "RecoTauTag/phase2Taus/plugins/PATTauClusterVariables.h"

#include "RecoTauTag/RecoTau/interface/RecoTauVertexAssociator.h"
#include "RecoTauTag/RecoTau/interface/ConeTools.h"
#include "RecoTauTag/RecoTau/interface/RecoTauQualityCuts.h"
#include "RecoTauTag/RecoTau/interface/TauDiscriminationProducerBase.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class timingTausMiniAOD : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  typedef reco::tau::cone::DeltaRPtrFilter<reco::PFCandidatePtr> DRFilter;

   public:
      explicit timingTausMiniAOD(const edm::ParameterSet&);
      ~timingTausMiniAOD();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<pat::TauCollection>     tauToken_;
  edm::EDGetTokenT<std::vector<pat::Jet> >    jetSrc_;
  edm::EDGetTokenT<reco::GenJetCollection> genJetSrc_;
  std::string tauID_;
  edm::EDGetTokenT<std::vector <reco::GenParticle> > prunedGenToken_;

  edm::EDGetTokenT<std::vector<pat::PackedCandidate> > pfCand_token;
  edm::ParameterSet qualityCutsPSet_;
  std::auto_ptr<reco::tau::RecoTauQualityCuts> pileupQcutsPUTrackSelection_;

  TTree* tree;
  TTree* jetTree;
  double tauPt_;
  double tauEta_;
  double jetPt_;
  double jetEta_;
  double tauMass_;
  double genTauPt_;
  double genTauEta_;
  double vtxX_, vtxY_, vtxZ_, vtxT_;
  double taudXY_, taudZ_, taudT_;
  int vtxIndex_;
  double tauCHIso_; 
  double tauCHIso0_;
  double tauCHIso1_;
  double tauCHIso2_;
  double tauCHIso3_;
  double tauCHIso4_;
  double tauCHIso5_;
  double tauCHIso6_;
  double tauCHIso7_;
  double tauCHIso8_;
  double tauCHIso9_;
  double tauCHIso10_;
  double tauCHIso11_;
  double tauCHIso12_;
  double tauCHIso13_;
  double tauCHIso14_;
  double tauCHIso15_;
  double tauCHIso16_;
  double tauCHIso17_;
  double tauCHIso18_;
  double tauCHIso19_;
  double tauCHIso20_;
  double tauCHIsoOther_;

  double tauCHIso0_ntracks_;
  double tauCHIso1_ntracks_;
  double tauCHIso2_ntracks_;
  double tauCHIso3_ntracks_;
  double tauCHIso4_ntracks_;
  double tauCHIso5_ntracks_;
  double tauCHIso6_ntracks_;
  double tauCHIso7_ntracks_;
  double tauCHIso8_ntracks_;
  double tauCHIso9_ntracks_;
  double tauCHIso10_ntracks_;
  double tauCHIso11_ntracks_;
  double tauCHIso12_ntracks_;
  double tauCHIso13_ntracks_;
  double tauCHIso14_ntracks_;
  double tauCHIso15_ntracks_;
  double tauCHIso16_ntracks_;
  double tauCHIso17_ntracks_;
  double tauCHIso18_ntracks_;
  double tauCHIso19_ntracks_;
  double tauCHIso20_ntracks_;
  double tauCHIsoOther_ntracks_;
  double z_2_;		   
  double pt_weighted_dr_signal_;
  double pt_weighted_deta_strip_;
  double pt_weighted_dphi_strip_;
  double pt_weighted_dr_iso_;   
  double n_photons_total_;	   
  double nIsoTracks_;	   
  double nIsoNeutral_;	   
  double nIsoGamma_;		   
  double nSigGamma_;            
  double nhIso_;
  double sigPhoIso_;

  int nvtx_;
  int dmf_;
  int goodReco_;
  int genTauMatch_;
  int jetTauMatch_;
  int genJetMatch_;
  double genJetPt_;
  double genJetEta_;

  reco::Candidate::LorentzVector GetVisibleP4(std::vector<const reco::GenParticle*>& daughters);
  void findDaughters(const reco::GenParticle* mother, std::vector<const reco::GenParticle*>& daughters);
  bool isNeutrino(const reco::Candidate* daughter);
  double getTimedCHIsoSum( pat::PackedCandidate const* leadTauCand, std::vector<pat::PackedCandidate const*> isolationCands, float interval, float &n_tracks);
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;


      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
timingTausMiniAOD::timingTausMiniAOD(const edm::ParameterSet& iConfig):
  vtxToken_      (consumes<reco::VertexCollection>          (iConfig.getParameter<edm::InputTag>("vertices"))),
  tauToken_      (consumes<pat::TauCollection>              (iConfig.getParameter<edm::InputTag>("taus"))),
  jetSrc_        (consumes<std::vector<pat::Jet> >    (iConfig.getParameter<edm::InputTag>("jets"))),
  genJetSrc_     (consumes<reco::GenJetCollection>          (iConfig.getParameter<edm::InputTag>("genJets"))),
  prunedGenToken_(consumes<std::vector<reco::GenParticle> > (iConfig.getParameter<edm::InputTag>("pruned"))),
  pfCand_token   (consumes<std::vector<pat::PackedCandidate> > (iConfig.getParameter<edm::InputTag>("pfCands"))),
  qualityCutsPSet_(iConfig.getParameter<edm::ParameterSet>("qualityCuts")) 
{
   //now do what ever initialization is needed
   usesResource("TFileService");

   tauID_    = iConfig.getParameter<std::string>("tauID");
   edm::Service<TFileService> fs;

   tree = fs->make<TTree>("Ntuple", "Ntuple");  
   tree->Branch("tauPt",       &tauPt_,       "tauPt/D"        );
   tree->Branch("tauEta",      &tauEta_,      "tauEta/D"       );
   tree->Branch("genTauPt",    &genTauPt_,    "genTauPt/D"     );
   tree->Branch("genTauEta",   &genTauEta_,   "genTauEta/D"    );
   tree->Branch("genTauMatch", &genTauMatch_, "genTauMatch/I"  );
   tree->Branch("nvtx",        &nvtx_,        "nvtx/I"         );
   tree->Branch("vtxX",        &vtxX_,        "vtxX/D"         );
   tree->Branch("vtxY",        &vtxY_,        "vtxY/D"         );
   tree->Branch("vtxZ",        &vtxZ_,        "vtxZ/D"         );
   tree->Branch("vtxT",        &vtxT_,        "vtxT/D"         );

   tree->Branch("taudXY",      &taudXY_,      "taudXY/D"       );
   tree->Branch("taudZ",       &taudZ_,       "taudZ/D"        );
   tree->Branch("taudT",       &taudT_,       "taudT/D"        );

   tree->Branch("dmf",           &dmf_,           "dmf/I"       );
   tree->Branch("goodReco",      &goodReco_,      "goodReco/I"  );
   tree->Branch("tauMass",       &tauMass_,       "tauMass/D"   );
   tree->Branch("tauCHIso",      &tauCHIso_ ,     "tauCHIso/D"  );
   tree->Branch("tauCHIso0",     &tauCHIso0_,     "tauCHIso0/D" );
   tree->Branch("tauCHIso1",     &tauCHIso1_,     "tauCHIso1/D" );
   tree->Branch("tauCHIso2",     &tauCHIso2_,     "tauCHIso2/D" );
   tree->Branch("tauCHIso3",     &tauCHIso3_,     "tauCHIso3/D" );
   tree->Branch("tauCHIso4",     &tauCHIso4_,     "tauCHIso4/D" );
   tree->Branch("tauCHIso5",     &tauCHIso5_,     "tauCHIso5/D" );
   tree->Branch("tauCHIso6",     &tauCHIso6_,     "tauCHIso6/D" );
   tree->Branch("tauCHIso7",     &tauCHIso7_,     "tauCHIso7/D" );
   tree->Branch("tauCHIso8",     &tauCHIso8_,     "tauCHIso8/D" );
   tree->Branch("tauCHIso9",     &tauCHIso9_,     "tauCHIso9/D" );
   tree->Branch("tauCHIso10",    &tauCHIso10_,    "tauCHIso10/D" );
   tree->Branch("tauCHIso11",    &tauCHIso11_,    "tauCHIso11/D" );
   tree->Branch("tauCHIso12",    &tauCHIso12_,    "tauCHIso12/D" );
   tree->Branch("tauCHIso13",    &tauCHIso13_,    "tauCHIso13/D" );
   tree->Branch("tauCHIso14",    &tauCHIso14_,    "tauCHIso14/D" );
   tree->Branch("tauCHIso15",    &tauCHIso15_,    "tauCHIso15/D" );
   tree->Branch("tauCHIso16",    &tauCHIso16_,    "tauCHIso16/D" );
   tree->Branch("tauCHIso17",    &tauCHIso17_,    "tauCHIso17/D" );
   tree->Branch("tauCHIso18",    &tauCHIso18_,    "tauCHIso18/D" );
   tree->Branch("tauCHIso19",    &tauCHIso19_,    "tauCHIso19/D" );
   tree->Branch("tauCHIso20",    &tauCHIso20_,    "tauCHIso20/D" );
   tree->Branch("tauCHIsoOther", &tauCHIsoOther_, "tauCHIsoOther/D");

   tree->Branch("tauCHIso0_ntracks",     &tauCHIso0_ntracks_,     "tauCHIso0_ntracks/D" );
   tree->Branch("tauCHIso1_ntracks",     &tauCHIso1_ntracks_,     "tauCHIso1_ntracks/D" );
   tree->Branch("tauCHIso2_ntracks",     &tauCHIso2_ntracks_,     "tauCHIso2_ntracks/D" );
   tree->Branch("tauCHIso3_ntracks",     &tauCHIso3_ntracks_,     "tauCHIso3_ntracks/D" );
   tree->Branch("tauCHIso4_ntracks",     &tauCHIso4_ntracks_,     "tauCHIso4_ntracks/D" );
   tree->Branch("tauCHIso5_ntracks",     &tauCHIso5_ntracks_,     "tauCHIso5_ntracks/D" );
   tree->Branch("tauCHIso6_ntracks",     &tauCHIso6_ntracks_,     "tauCHIso6_ntracks/D" );
   tree->Branch("tauCHIso7_ntracks",     &tauCHIso7_ntracks_,     "tauCHIso7_ntracks/D" );
   tree->Branch("tauCHIso8_ntracks",     &tauCHIso8_ntracks_,     "tauCHIso8_ntracks/D" );
   tree->Branch("tauCHIso9_ntracks",     &tauCHIso9_ntracks_,     "tauCHIso9_ntracks/D" );
   tree->Branch("tauCHIso10_ntracks",     &tauCHIso10_ntracks_,     "tauCHIso10_ntracks/D" );
   tree->Branch("tauCHIso11_ntracks",     &tauCHIso11_ntracks_,     "tauCHIso11_ntracks/D" );
   tree->Branch("tauCHIso12_ntracks",     &tauCHIso12_ntracks_,     "tauCHIso12_ntracks/D" );
   tree->Branch("tauCHIso13_ntracks",     &tauCHIso13_ntracks_,     "tauCHIso13_ntracks/D" );
   tree->Branch("tauCHIso14_ntracks",     &tauCHIso14_ntracks_,     "tauCHIso14_ntracks/D" );
   tree->Branch("tauCHIso15_ntracks",     &tauCHIso15_ntracks_,     "tauCHIso15_ntracks/D" );
   tree->Branch("tauCHIso16_ntracks",     &tauCHIso16_ntracks_,     "tauCHIso16_ntracks/D" );
   tree->Branch("tauCHIso17_ntracks",     &tauCHIso17_ntracks_,     "tauCHIso17_ntracks/D" );
   tree->Branch("tauCHIso18_ntracks",     &tauCHIso18_ntracks_,     "tauCHIso18_ntracks/D" );
   tree->Branch("tauCHIso19_ntracks",     &tauCHIso19_ntracks_,     "tauCHIso19_ntracks/D" );
   tree->Branch("tauCHIso20_ntracks",     &tauCHIso20_ntracks_,     "tauCHIso20_ntracks/D" );

   tree->Branch("tauCHIsoOther_ntracks", &tauCHIsoOther_ntracks_, "tauCHIsoOther_ntracks/D");

   tree->Branch("vtxIndex" ,     &vtxIndex_,      "vtxIndex_/I" );

   tree->Branch("pt_weighted_dr_signal_" ,   &pt_weighted_dr_signal_,  "pt_weighted_dr_signal/D"  );
   tree->Branch("pt_weighted_deta_strip_" ,  &pt_weighted_deta_strip_, "pt_weighted_deta_strip/D" );
   tree->Branch("pt_weighted_dphi_strip_" ,  &pt_weighted_dphi_strip_, "pt_weighted_dphi_strip/D" );
   tree->Branch("pt_weighted_dr_iso_" ,      &pt_weighted_dr_iso_,     "pt_weighted_dr_iso/D"     );
   tree->Branch("n_photons_total_" ,         &n_photons_total_,        "n_photons_total/D"        );
   tree->Branch("nIsoTracks_" ,              &nIsoTracks_,             "nIsoTracks/D"             );
   tree->Branch("nIsoNeutral_" ,             &nIsoNeutral_,            "nIsoNeutral/D"            );
   tree->Branch("nIsoGamma_" ,               &nIsoGamma_,              "nIsoGamma/D"              );
   tree->Branch("nSigGamma_" ,               &nSigGamma_,              "nSigGamma/D"              );
   tree->Branch("nhIso_" ,                   &nhIso_,                  "nhIso/D"                  );
   tree->Branch("sigPhoIso_" ,               &sigPhoIso_,              "sigPhoIso/D"              );


   jetTree = fs->make<TTree>(      "jetNtuple",   "jetNtuple"       );
   jetTree->Branch("tauPt",       &tauPt_,       "tauPt/D"        );
   jetTree->Branch("tauEta",      &tauEta_,      "tauEta/D"       );
   jetTree->Branch("genTauPt",    &genTauPt_,    "genTauPt/D"     );
   jetTree->Branch("genTauEta",   &genTauEta_,   "genTauEta/D"    );
   jetTree->Branch("genTauMatch", &genTauMatch_, "genTauMatch/I"  );
   jetTree->Branch("jetPt",        &jetPt_,       "jetPt_/D"        );
   jetTree->Branch("jetEta",       &jetEta_,      "jetEta_/D"       );
   jetTree->Branch("jetTauMatch",  &jetTauMatch_, "jetTauMatch_/I"  );
   jetTree->Branch("genJetMatch",  &genJetMatch_, "genJetMatch_/I"  );
   jetTree->Branch("genJetPt",     &genJetPt_,    "genJetPt_/D"  );
   jetTree->Branch("genJetEta",    &genJetEta_,   "genJetEta_/D"  );
   jetTree->Branch("nvtx",        &nvtx_,        "nvtx/I"         );
   jetTree->Branch("vtxX",        &vtxX_,        "vtxX/D"         );
   jetTree->Branch("vtxY",        &vtxY_,        "vtxY/D"         );
   jetTree->Branch("vtxZ",        &vtxZ_,        "vtxZ/D"         );
   jetTree->Branch("vtxT",        &vtxT_,        "vtxT/D"         );

   jetTree->Branch("taudXY",      &taudXY_,      "taudXY/D"        );
   jetTree->Branch("taudZ",       &taudZ_,       "taudZ/D"        );
   jetTree->Branch("taudT",       &taudT_,       "taudT/D"        );

   jetTree->Branch("dmf",           &dmf_,           "dmf/I"       );
   jetTree->Branch("goodReco",      &goodReco_,      "goodReco/I"  );
   jetTree->Branch("tauMass",       &tauMass_,       "tauMass/D"   );
   jetTree->Branch("tauCHIso",      &tauCHIso_ ,     "tauCHIso/D"  );
   jetTree->Branch("tauCHIso0",     &tauCHIso0_,     "tauCHIso0/D" );
   jetTree->Branch("tauCHIso1",     &tauCHIso1_,     "tauCHIso1/D" );
   jetTree->Branch("tauCHIso2",     &tauCHIso2_,     "tauCHIso2/D" );
   jetTree->Branch("tauCHIso3",     &tauCHIso3_,     "tauCHIso3/D" );
   jetTree->Branch("tauCHIso4",     &tauCHIso4_,     "tauCHIso4/D" );
   jetTree->Branch("tauCHIso5",     &tauCHIso5_,     "tauCHIso5/D" );
   jetTree->Branch("tauCHIso6",     &tauCHIso6_,     "tauCHIso6/D" );
   jetTree->Branch("tauCHIso7",     &tauCHIso7_,     "tauCHIso7/D" );
   jetTree->Branch("tauCHIso8",     &tauCHIso8_,     "tauCHIso8/D" );
   jetTree->Branch("tauCHIso9",     &tauCHIso9_,     "tauCHIso9/D" );
   jetTree->Branch("tauCHIso10",    &tauCHIso10_,    "tauCHIso10/D" );
   jetTree->Branch("tauCHIso11",    &tauCHIso11_,    "tauCHIso11/D" );
   jetTree->Branch("tauCHIso12",    &tauCHIso12_,    "tauCHIso12/D" );
   jetTree->Branch("tauCHIso13",    &tauCHIso13_,    "tauCHIso13/D" );
   jetTree->Branch("tauCHIso14",    &tauCHIso14_,    "tauCHIso14/D" );
   jetTree->Branch("tauCHIso15",    &tauCHIso15_,    "tauCHIso15/D" );
   jetTree->Branch("tauCHIso16",    &tauCHIso16_,    "tauCHIso16/D" );
   jetTree->Branch("tauCHIso17",    &tauCHIso17_,    "tauCHIso17/D" );
   jetTree->Branch("tauCHIso18",    &tauCHIso18_,    "tauCHIso18/D" );
   jetTree->Branch("tauCHIso19",    &tauCHIso19_,    "tauCHIso19/D" );
   jetTree->Branch("tauCHIso20",    &tauCHIso20_,    "tauCHIso20/D" );
   jetTree->Branch("tauCHIsoOther", &tauCHIsoOther_, "tauCHIsoOther/D");
   jetTree->Branch("vtxIndex" ,     &vtxIndex_,      "vtxIndex_/I" );

   jetTree->Branch("tauCHIso0_ntracks",     &tauCHIso0_ntracks_,     "tauCHIso0_ntracks/D" );
   jetTree->Branch("tauCHIso1_ntracks",     &tauCHIso1_ntracks_,     "tauCHIso1_ntracks/D" );
   jetTree->Branch("tauCHIso2_ntracks",     &tauCHIso2_ntracks_,     "tauCHIso2_ntracks/D" );
   jetTree->Branch("tauCHIso3_ntracks",     &tauCHIso3_ntracks_,     "tauCHIso3_ntracks/D" );
   jetTree->Branch("tauCHIso4_ntracks",     &tauCHIso4_ntracks_,     "tauCHIso4_ntracks/D" );
   jetTree->Branch("tauCHIso5_ntracks",     &tauCHIso5_ntracks_,     "tauCHIso5_ntracks/D" );
   jetTree->Branch("tauCHIso6_ntracks",     &tauCHIso6_ntracks_,     "tauCHIso6_ntracks/D" );
   jetTree->Branch("tauCHIso7_ntracks",     &tauCHIso7_ntracks_,     "tauCHIso7_ntracks/D" );
   jetTree->Branch("tauCHIso8_ntracks",     &tauCHIso8_ntracks_,     "tauCHIso8_ntracks/D" );
   jetTree->Branch("tauCHIso9_ntracks",     &tauCHIso9_ntracks_,     "tauCHIso9_ntracks/D" );
   jetTree->Branch("tauCHIso10_ntracks",     &tauCHIso10_ntracks_,     "tauCHIso10_ntracks/D" );
   jetTree->Branch("tauCHIso11_ntracks",     &tauCHIso11_ntracks_,     "tauCHIso11_ntracks/D" );
   jetTree->Branch("tauCHIso12_ntracks",     &tauCHIso12_ntracks_,     "tauCHIso12_ntracks/D" );
   jetTree->Branch("tauCHIso13_ntracks",     &tauCHIso13_ntracks_,     "tauCHIso13_ntracks/D" );
   jetTree->Branch("tauCHIso14_ntracks",     &tauCHIso14_ntracks_,     "tauCHIso14_ntracks/D" );
   jetTree->Branch("tauCHIso15_ntracks",     &tauCHIso15_ntracks_,     "tauCHIso15_ntracks/D" );
   jetTree->Branch("tauCHIso16_ntracks",     &tauCHIso16_ntracks_,     "tauCHIso16_ntracks/D" );
   jetTree->Branch("tauCHIso17_ntracks",     &tauCHIso17_ntracks_,     "tauCHIso17_ntracks/D" );
   jetTree->Branch("tauCHIso18_ntracks",     &tauCHIso18_ntracks_,     "tauCHIso18_ntracks/D" );
   jetTree->Branch("tauCHIso19_ntracks",     &tauCHIso19_ntracks_,     "tauCHIso19_ntracks/D" );
   jetTree->Branch("tauCHIso20_ntracks",     &tauCHIso20_ntracks_,     "tauCHIso20_ntracks/D" );

   jetTree->Branch("tauCHIsoOther_ntracks", &tauCHIsoOther_ntracks_, "tauCHIsoOther_ntracks/D");

   jetTree->Branch("pt_weighted_dr_signal_" ,   &pt_weighted_dr_signal_,  "pt_weighted_dr_signal/D"  );
   jetTree->Branch("pt_weighted_deta_strip_" ,  &pt_weighted_deta_strip_, "pt_weighted_deta_strip/D" );
   jetTree->Branch("pt_weighted_dphi_strip_" ,  &pt_weighted_dphi_strip_, "pt_weighted_dphi_strip/D" );
   jetTree->Branch("pt_weighted_dr_iso_" ,      &pt_weighted_dr_iso_,     "pt_weighted_dr_iso/D"     );
   jetTree->Branch("n_photons_total_" ,         &n_photons_total_,        "n_photons_total/D"        );
   jetTree->Branch("nIsoTracks_" ,              &nIsoTracks_,             "nIsoTracks/D"             );
   jetTree->Branch("nIsoNeutral_" ,             &nIsoNeutral_,            "nIsoNeutral/D"            );
   jetTree->Branch("nIsoGamma_" ,               &nIsoGamma_,              "nIsoGamma/D"              );
   jetTree->Branch("nSigGamma_" ,               &nSigGamma_,              "nSigGamma/D"              );
   jetTree->Branch("nhIso_" ,                   &nhIso_,                  "nhIso/D"                  );
   jetTree->Branch("sigPhoIso_" ,               &sigPhoIso_,              "sigPhoIso/D"              );

}


timingTausMiniAOD::~timingTausMiniAOD()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
timingTausMiniAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   //initialize tau isolation
   std::auto_ptr<reco::tau::RecoTauQualityCuts> pileupQcutsGeneralQCuts_;

   edm::ParameterSet isolationQCuts = qualityCutsPSet_.getParameterSet("isolationQualityCuts");
   std::vector<reco::PFCandidatePtr> chargedPFCandidatesInEvent_;
   /*
   edm::Handle<vector<pat::PackedCandidate> > pfCandidates;
   iEvent.getByToken(pfCand_token, pfCandidates);

   chargedPFCandidatesInEvent_.clear();
   chargedPFCandidatesInEvent_.reserve(pfCandidates->size());

   size_t numPFCandidates = pfCandidates->size();
   for ( size_t i = 0; i < numPFCandidates; ++i ) {
     reco::PFCandidatePtr pfCandidate(pfCandidates, i);
     if ( pfCandidate->charge() != 0 ) {
       chargedPFCandidatesInEvent_.push_back(pfCandidate);
     }
     }*/

   edm::Handle<reco::VertexCollection> vtxs;
   iEvent.getByToken(vtxToken_, vtxs);
   nvtx_=vtxs->size();
   edm::Handle<pat::TauCollection> taus;
   iEvent.getByToken(tauToken_, taus);
   edm::Handle<std::vector<reco::GenParticle> > genParticles;
   iEvent.getByToken(prunedGenToken_, genParticles);
   std::vector<const reco::GenParticle*> GenTaus;

   std::unique_ptr<pat::TauCollection> newTaus(new pat::TauCollection);

   for(std::vector<reco::GenParticle>::const_iterator genParticle = genParticles->begin(); genParticle != genParticles->end(); genParticle++){
     if(TMath::Abs(genParticle->pdgId()) == 15) GenTaus.push_back(&(*genParticle));
   }

   //get objects
   Handle <std::vector<pat::Jet> > JetObjects;
   if(!iEvent.getByToken(jetSrc_, JetObjects))
     std::cout<<"Error getting Jets"<<std::endl;
   std::vector<pat::Jet> Jets;

   for (unsigned int iJet = 0; iJet < JetObjects->size() ; ++iJet){
     //reco::PFJetRef jetCand(JetObjects, iJet);
     pat::Jet jetCand = JetObjects->at(iJet);
     if(jetCand.pt() < 18 )continue;
     bool isATau=false;
     for(auto genTau : GenTaus){
       std::vector<const reco::GenParticle*> genTauDaughters;
       findDaughters(genTau, genTauDaughters);
       reco::Candidate::LorentzVector genTauVis = GetVisibleP4(genTauDaughters);
       genTauPt_  = (float) genTauVis.pt();
       genTauEta_ = (float) genTauVis.eta();
       if (reco::deltaR(jetCand.eta(),jetCand.phi(),genTauVis.eta(),genTauVis.phi()) < 0.5)
	 isATau=true;
       
     }
     if(!isATau)
       Jets.push_back(jetCand);
   }   

   //Build Iso 
   for(unsigned int  i=0;i!=taus->size();++i){
     pat::Tau tau = taus->at(i);

     float dZ  = -999;
     float dXY = -999;
     float z_2=  -999;
     if(vtxs->size()){z_2 = tau.vz() + (130. / tan(tau.theta()));}
     else std::cout<<"VERTICES NULL"<<std::endl;

     float dm = tau.decayMode();
     float pt_weighted_dr_signal  = tau_pt_weighted_dr_signal(tau,dm); 
     float pt_weighted_deta_strip = tau_pt_weighted_deta_strip(tau,dm);
     float pt_weighted_dphi_strip = tau_pt_weighted_dphi_strip(tau,dm);
     float pt_weighted_dr_iso     = tau_pt_weighted_dr_iso(tau,dm);
     float n_photons_total        = tau_n_photons_total(tau);

     float nIsoTracks  = tau.isolationChargedHadrCands().size();
     float nIsoNeutral = tau.isolationNeutrHadrCands().size();
     float nIsoGamma   = tau.isolationGammaCands().size();
     float nSigGamma   = tau.signalGammaCands().size();



     tau.addUserFloat("zIP",                        z_2);
     tau.addUserFloat("tau_pt_weighted_dr_signal" , pt_weighted_dr_signal);
     tau.addUserFloat("tau_pt_weighted_deta_strip", pt_weighted_deta_strip);
     tau.addUserFloat("tau_pt_weighted_dphi_strip", pt_weighted_dphi_strip);
     tau.addUserFloat("tau_pt_weighted_dr_iso",     pt_weighted_dr_iso);
     tau.addUserFloat("n_photons_total",            n_photons_total);
     tau.addUserFloat("nIsoTracks",                 nIsoTracks);
     tau.addUserFloat("nIsoNeutral",                nIsoNeutral);
     tau.addUserFloat("nIsoGamma",                  nIsoGamma);
     tau.addUserFloat("nSigGamma",                  nSigGamma);              
     
    double chargedPtForDB = 0;
    for ( auto const & isoObject : tau.isolationPFChargedHadrCands() ) {
      chargedPtForDB += isoObject->pt();
    }
    double neutralPtForDB = 0.;
    for ( auto const & isoObject : tau.isolationPFGammaCands() ) {
      neutralPtForDB += isoObject->pt();
    }
    double puPtForDB;
    std::pair<edm::ParameterSet, edm::ParameterSet> puFactorizedIsoQCuts;

    std::vector<reco::PFCandidatePtr> allPU;

    /// make pu filtertered 
    std::vector<reco::PFCandidatePtr> cleanPU =
      pileupQcutsGeneralQCuts_->filterCandRefs(allPU);

    // Only select PU tracks inside the isolation cone.
    std::vector<reco::PFCandidatePtr> isoPU_;
    DRFilter deltaBetaFilter(tau.p4(), 0, 0.5); //0.5 is the cone size for isolation
    for ( auto const & cand : cleanPU ) {
      if ( deltaBetaFilter(cand) )  isoPU_.push_back(cand);
    }

    for ( auto const & isoObject : isoPU_ ) {
      puPtForDB += isoObject->pt();
    }

    neutralPtForDB -= (0.2*puPtForDB);

    double totalPtForDB = chargedPtForDB + neutralPtForDB;


     float oldChIso  = 0;
     float newChIso0 = 0;
     float newChIso1 = 0;
     float newChIso2 = 0;
     float newChIso3 = 0;
     float newChIso4 = 0;
     float newChIso5 = 0;
     float newChIso6 = 0;
     float newChIso7 = 0;
     float newChIso8 = 0;
     float newChIso9 = 0;
     float newChIso10 = 0;
     float newChIso11 = 0;
     float newChIso12 = 0;
     float newChIso13 = 0;
     float newChIso14 = 0;
     float newChIso15 = 0;
     float newChIso16 = 0;
     float newChIso17 = 0;
     float newChIso18 = 0;
     float newChIso19 = 0;
     float newChIso20 = 0;
     float newChIsoOther = 0;

     float newChIso0_ntracks = 0;
     float newChIso1_ntracks = 0;
     float newChIso2_ntracks = 0;
     float newChIso3_ntracks = 0;
     float newChIso4_ntracks = 0;
     float newChIso5_ntracks = 0;
     float newChIso6_ntracks = 0;
     float newChIso7_ntracks = 0;
     float newChIso8_ntracks = 0;
     float newChIso9_ntracks = 0;
     float newChIso10_ntracks = 0;
     float newChIso11_ntracks = 0;
     float newChIso12_ntracks = 0;
     float newChIso13_ntracks = 0;
     float newChIso14_ntracks = 0;
     float newChIso15_ntracks = 0;
     float newChIso16_ntracks = 0;
     float newChIso17_ntracks = 0;
     float newChIso18_ntracks = 0;
     float newChIso19_ntracks = 0;
     float newChIso20_ntracks = 0;
     float newChIsoOther_ntracks = 0;

     float nhIso     = 0.; 
     float sigPhoIso = 0.; 
     std::cout<<"N Neutral Hadrons "<< tau.isolationNeutrHadrCands().size()<<std::endl;
     for(auto cand : tau.isolationChargedHadrCands() ) {oldChIso  += cand->pt();}
     for(auto cand : tau.isolationNeutrHadrCands() )   {nhIso     += cand->pt();}
     for(auto cand : tau.signalGammaCands() )          {sigPhoIso += cand->pt();}

     tau.addUserFloat("oldChIso",   oldChIso);
     tau.addUserFloat("nhIsoPt",    nhIso);
     tau.addUserFloat("sigPhoIsoPt",sigPhoIso);

     pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(tau.leadChargedHadrCand().get());

     //std::cout<<"packedLeadTauCand timing "<< packedLeadTauCand->dtimeAssociatedPV()<<std::endl;
     //create CHIso sum for 0,1,2,3,4,5,6,7,8, and all timing intervals
     //for each charge iso track check 1.) difference from central track, 2.) add to CHIso Sum interval
     //central cand, isocands, interval, value to return
     std::vector<pat::PackedCandidate const*> packedCHIsoTauCands; 
     for(auto cand : tau.isolationChargedHadrCands() ){
       pat::PackedCandidate const* tempIsoCand = dynamic_cast<pat::PackedCandidate const*>(cand.get());
       packedCHIsoTauCands.push_back(tempIsoCand);
     }

     const float dz = std::abs(packedLeadTauCand->dz(vtxs->at(0).position()));
     const float dxy = std::abs(packedLeadTauCand->dxy(vtxs->at(0).position()));
     const float dt = packedLeadTauCand->dtimeAssociatedPV();

     tau.addUserFloat("tauDZ",   dz);
     tau.addUserFloat("tauDXY",  dxy);
     tau.addUserFloat("tauDT",   dt);
     const double res_det = 0.030; // ns
     //const float dt = std::abs(time - thevtx.t());

     newChIso0 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 0*res_det, newChIso0_ntracks );
     newChIso1 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 1*res_det, newChIso1_ntracks );
     newChIso2 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 2*res_det, newChIso2_ntracks );
     newChIso3 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 3*res_det, newChIso3_ntracks ); // this is the nominal cut, 3 sigma
     newChIso4 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 4*res_det, newChIso4_ntracks );
     newChIso5 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 5*res_det, newChIso5_ntracks );
     newChIso6 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 6*res_det, newChIso6_ntracks );
     newChIso7 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 7*res_det, newChIso7_ntracks );
     newChIso8 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 8*res_det, newChIso8_ntracks );
     newChIso9 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 9*res_det, newChIso9_ntracks );
     newChIso10 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 10*res_det, newChIso10_ntracks );
     newChIso11 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 11*res_det, newChIso11_ntracks );
     newChIso12 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 12*res_det, newChIso12_ntracks );
     newChIso13 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 13*res_det, newChIso13_ntracks );
     newChIso14 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 14*res_det, newChIso14_ntracks );
     newChIso15 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 15*res_det, newChIso15_ntracks );
     newChIso16 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 16*res_det, newChIso16_ntracks );
     newChIso17 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 17*res_det, newChIso17_ntracks );
     newChIso18 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 18*res_det, newChIso18_ntracks );
     newChIso19 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 19*res_det, newChIso19_ntracks );
     newChIso20 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 20*res_det, newChIso20_ntracks );
     newChIsoOther = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 100*res_det, newChIsoOther_ntracks);

     

     tau.addUserFloat("newChIso0",newChIso0);
     tau.addUserFloat("newChIso1",newChIso1);
     tau.addUserFloat("newChIso2",newChIso2);
     tau.addUserFloat("newChIso3",newChIso3);
     tau.addUserFloat("newChIso4",newChIso4);
     tau.addUserFloat("newChIso5",newChIso5);
     tau.addUserFloat("newChIso6",newChIso6);
     tau.addUserFloat("newChIso7",newChIso7);
     tau.addUserFloat("newChIso8",newChIso8);
     tau.addUserFloat("newChIso9",newChIso9);
     tau.addUserFloat("newChIso10",newChIso10);
     tau.addUserFloat("newChIso11",newChIso11);
     tau.addUserFloat("newChIso12",newChIso12);
     tau.addUserFloat("newChIso13",newChIso13);
     tau.addUserFloat("newChIso14",newChIso14);
     tau.addUserFloat("newChIso15",newChIso15);
     tau.addUserFloat("newChIso16",newChIso16);
     tau.addUserFloat("newChIso17",newChIso17);
     tau.addUserFloat("newChIso18",newChIso18);
     tau.addUserFloat("newChIso19",newChIso19);
     tau.addUserFloat("newChIso20",newChIso20);

     tau.addUserFloat("newChIsoOther",newChIsoOther);

     tau.addUserFloat("newChIso0_ntracks",newChIso0_ntracks );
     tau.addUserFloat("newChIso1_ntracks",newChIso1_ntracks );
     tau.addUserFloat("newChIso2_ntracks",newChIso2_ntracks );
     tau.addUserFloat("newChIso3_ntracks",newChIso3_ntracks );
     tau.addUserFloat("newChIso4_ntracks",newChIso4_ntracks );
     tau.addUserFloat("newChIso5_ntracks",newChIso5_ntracks );
     tau.addUserFloat("newChIso6_ntracks",newChIso6_ntracks );
     tau.addUserFloat("newChIso7_ntracks",newChIso7_ntracks );
     tau.addUserFloat("newChIso8_ntracks",newChIso8_ntracks );
     tau.addUserFloat("newChIso9_ntracks",newChIso9_ntracks );
     tau.addUserFloat("newChIso10_ntracks",newChIso10_ntracks );
     tau.addUserFloat("newChIso11_ntracks",newChIso11_ntracks );
     tau.addUserFloat("newChIso12_ntracks",newChIso12_ntracks );
     tau.addUserFloat("newChIso13_ntracks",newChIso13_ntracks );
     tau.addUserFloat("newChIso14_ntracks",newChIso14_ntracks );
     tau.addUserFloat("newChIso15_ntracks",newChIso15_ntracks );
     tau.addUserFloat("newChIso16_ntracks",newChIso16_ntracks );
     tau.addUserFloat("newChIso17_ntracks",newChIso17_ntracks );
     tau.addUserFloat("newChIso18_ntracks",newChIso18_ntracks );
     tau.addUserFloat("newChIso19_ntracks",newChIso19_ntracks );
     tau.addUserFloat("newChIso20_ntracks",newChIso20_ntracks );
     tau.addUserFloat("newChIsoOther_ntracks",newChIsoOther_ntracks );

     //integrate all ingredients for DB Iso study
     //also check as a function of DZ
     newTaus->push_back(tau);
   }

   for(auto genTau : GenTaus){
      tauPt_     = -10;
      tauEta_    = -10;
      tauMass_   = -10;
      genTauPt_  = -10;
      genTauEta_ = -10;
      tauCHIso_  = -10;
      tauCHIso0_ = -10;
      tauCHIso1_ = -10;
      tauCHIso2_ = -10;     
      tauCHIso3_ = -10;
      tauCHIso4_ = -10;
      tauCHIso5_ = -10;
      tauCHIso6_ = -10;
      tauCHIso7_ = -10;
      tauCHIso8_ = -10;
      tauCHIso9_  = -10;
      tauCHIso10_  = -10;
      tauCHIso11_  = -10;
      tauCHIso12_  = -10;
      tauCHIso13_  = -10;
      tauCHIso14_  = -10;
      tauCHIso15_  = -10;
      tauCHIso16_  = -10;
      tauCHIso17_  = -10;
      tauCHIso18_  = -10;
      tauCHIso19_  = -10;
      tauCHIso20_  = -10;
      tauCHIsoOther_ = -10;

      tauCHIso0_ntracks_ = -10;
      tauCHIso1_ntracks_ = -10;
      tauCHIso2_ntracks_ = -10;     
      tauCHIso3_ntracks_ = -10;
      tauCHIso4_ntracks_ = -10;
      tauCHIso5_ntracks_ = -10;
      tauCHIso6_ntracks_ = -10;
      tauCHIso7_ntracks_ = -10;
      tauCHIso8_ntracks_ = -10;
      tauCHIso9_ntracks_ = -10;
      tauCHIso10_ntracks_ = -10;
      tauCHIso11_ntracks_ = -10;
      tauCHIso12_ntracks_ = -10;
      tauCHIso13_ntracks_ = -10;
      tauCHIso14_ntracks_ = -10;
      tauCHIso15_ntracks_ = -10;
      tauCHIso16_ntracks_ = -10;
      tauCHIso17_ntracks_ = -10;
      tauCHIso18_ntracks_ = -10;
      tauCHIso19_ntracks_ = -10;
      tauCHIso20_ntracks_ = -10;
      tauCHIsoOther_ntracks_ = -10;

      vtxX_ =  -10;
      vtxY_ =  -10;
      vtxZ_ =  -10;
      vtxT_ =  -10;

      taudXY_ =  -10;
      taudZ_ =  -10;
      taudT_ =  -10;

      vtxIndex_ = -10;
      pt_weighted_dr_signal_ = -10; 
      pt_weighted_deta_strip_ = -10;
      pt_weighted_dphi_strip_ = -10;
      pt_weighted_dr_iso_ = -10;    
      n_photons_total_ = -10;       
      nIsoTracks_ = -10;            
      nIsoNeutral_ = -10;           
      nIsoGamma_ = -10;             
      nSigGamma_ = -10;             
      nhIso_ = -10;                 
      sigPhoIso_ = -10;             
      dmf_=-10;
      goodReco_=0;
      genTauMatch_=0;
      
      std::vector<const reco::GenParticle*> genTauDaughters;
      findDaughters(genTau, genTauDaughters);
      reco::Candidate::LorentzVector genTauVis = GetVisibleP4(genTauDaughters);
      genTauPt_  = (float) genTauVis.pt();
      genTauEta_ = (float) genTauVis.eta();
      
     //std::cout<<" pt "<<genTauPt_<<" eta "<<genTauEta_<<std::endl;
      for(const pat::Tau &tau : *newTaus){
	if (reco::deltaR(tau.eta(),tau.phi(),genTauVis.eta(),genTauVis.phi()) < 0.5 && tau.tauID("decayModeFinding")>0.5){

	  genTauMatch_ = 1;
	  tauPt_  = tau.pt();
	  tauEta_ = tau.eta();
	  dmf_ = tau.decayMode();
	  tauMass_ = tau.mass();

	  //get the matched vertex
	  int vtx_index = -1;
	  // find the 4D vertex this muon is best associated to..
	  float max_weight = 0.f;
	  for( unsigned i = 0; i < vtxs->size(); ++i ) {
	    const auto& vtx = (*vtxs)[i];      
	    const float weight = vtx.trackWeight(tau.leadChargedHadrCand());
	    if( weight > max_weight ) {
	      max_weight = weight;
	      vtx_index = i;
	    }
	  }
	  //now do vtx variable filling
	  vtxIndex_ = vtx_index;

	  const reco::Vertex& vtx = (vtx_index == -1 ? (*vtxs)[0] : (*vtxs)[vtx_index]);
	  vtxX_ = vtx.x();
	  vtxY_ = vtx.y();
	  vtxZ_ = vtx.z();
	  vtxT_ = vtx.t();

	  taudXY_ = tau.userFloat("tauDXY");
	  taudZ_ = tau.userFloat("tauDZ");
	  taudT_ = tau.userFloat("tauDT");

	  pt_weighted_dr_signal_ = tau.userFloat("tau_pt_weighted_dr_signal"); 
	  pt_weighted_deta_strip_ = tau.userFloat("tau_pt_weighted_deta_strip");
	  pt_weighted_dphi_strip_ = tau.userFloat("tau_pt_weighted_dphi_strip");
	  pt_weighted_dr_iso_ = tau.userFloat("tau_pt_weighted_dr_iso");	  
	  n_photons_total_ = tau.userFloat("n_photons_total");       
	  nIsoTracks_ = tau.userFloat("nIsoTracks");	  
	  nIsoNeutral_ = tau.userFloat("nIsoNeutral");           
	  nIsoGamma_ = tau.userFloat("nIsoGamma");	  
	  nSigGamma_ = tau.userFloat("nSigGamma");             
	  nhIso_ = tau.userFloat("nhIsoPt");	  
	  sigPhoIso_ = tau.userFloat("sigPhoIsoPt");             

	  tauCHIso_  = tau.userFloat("oldChIso");
	  tauCHIso0_ = tau.userFloat("newChIso0");
	  tauCHIso1_ = tau.userFloat("newChIso1");
	  tauCHIso2_ = tau.userFloat("newChIso2");
	  tauCHIso3_ = tau.userFloat("newChIso3");
	  tauCHIso4_ = tau.userFloat("newChIso4");
	  tauCHIso5_ = tau.userFloat("newChIso5");
	  tauCHIso6_ = tau.userFloat("newChIso6");
	  tauCHIso7_ = tau.userFloat("newChIso7");
	  tauCHIso8_ = tau.userFloat("newChIso8");
	  tauCHIso9_ = tau.userFloat("newChIso9");
	  tauCHIso10_ = tau.userFloat("newChIso10");
	  tauCHIso11_ = tau.userFloat("newChIso11");
	  tauCHIso12_ = tau.userFloat("newChIso12");
	  tauCHIso13_ = tau.userFloat("newChIso13");
	  tauCHIso14_ = tau.userFloat("newChIso14");
	  tauCHIso15_ = tau.userFloat("newChIso15");
	  tauCHIso16_ = tau.userFloat("newChIso16");
	  tauCHIso17_ = tau.userFloat("newChIso17");
	  tauCHIso18_ = tau.userFloat("newChIso18");
	  tauCHIso19_ = tau.userFloat("newChIso19");
	  tauCHIso20_ = tau.userFloat("newChIso20");
	  tauCHIsoOther_ = tau.userFloat("newChIsoOther");
	  
	  //tauCHIso_ntracks_  = tau.userFloat("oldChIso");
	  tauCHIso0_ntracks_ = tau.userFloat("newChIso0_ntracks");
	  tauCHIso1_ntracks_ = tau.userFloat("newChIso1_ntracks");
	  tauCHIso2_ntracks_ = tau.userFloat("newChIso2_ntracks");
	  tauCHIso3_ntracks_ = tau.userFloat("newChIso3_ntracks");
	  tauCHIso4_ntracks_ = tau.userFloat("newChIso4_ntracks");
	  tauCHIso5_ntracks_ = tau.userFloat("newChIso5_ntracks");
	  tauCHIso6_ntracks_ = tau.userFloat("newChIso6_ntracks");
	  tauCHIso7_ntracks_ = tau.userFloat("newChIso7_ntracks");
	  tauCHIso8_ntracks_ = tau.userFloat("newChIso8_ntracks");
	  tauCHIso9_ntracks_ = tau.userFloat("newChIso9_ntracks");
	  tauCHIso10_ntracks_ = tau.userFloat("newChIso10_ntracks");
	  tauCHIso11_ntracks_ = tau.userFloat("newChIso11_ntracks");
	  tauCHIso12_ntracks_ = tau.userFloat("newChIso12_ntracks");
	  tauCHIso13_ntracks_ = tau.userFloat("newChIso13_ntracks");
	  tauCHIso14_ntracks_ = tau.userFloat("newChIso14_ntracks");
	  tauCHIso15_ntracks_ = tau.userFloat("newChIso15_ntracks");
	  tauCHIso16_ntracks_ = tau.userFloat("newChIso16_ntracks");
	  tauCHIso17_ntracks_ = tau.userFloat("newChIso17_ntracks");
	  tauCHIso18_ntracks_ = tau.userFloat("newChIso18_ntracks");
	  tauCHIso19_ntracks_ = tau.userFloat("newChIso19_ntracks");
	  tauCHIso20_ntracks_ = tau.userFloat("newChIso20_ntracks");
	  tauCHIsoOther_ntracks_ = tau.userFloat("newChIsoOther_ntracks");

	  goodReco_ = (bool) tau.tauID(tauID_) >0.5;
	  //std::cout<<"iso0: "<<tauCHIso0_<<" iso1: "<<tauCHIso1_<<" iso8: "<<tauCHIso8_<<" iso: "<<tauCHIso_<<" isoOther: "<<tauCHIsoOther_<<std::endl;
	  break;
	}
      }
     tree->Fill(); 
   }


   edm::Handle<reco::GenJetCollection> genJets;
   if(!iEvent.getByToken(genJetSrc_, genJets))
     std::cout<<"Error getting genJets"<<std::endl;


   //for(auto genTau : genJets){
   for(auto jet : Jets){

      tauPt_     = -10;
      tauEta_    = -10;
      tauMass_   = -10;
      genTauPt_  = -10;
      genTauEta_ = -10;
      tauCHIso_  = -10;
      tauCHIso0_ = -10;
      tauCHIso1_ = -10;
      tauCHIso2_ = -10;     
      tauCHIso3_ = -10;
      tauCHIso4_ = -10;
      tauCHIso5_ = -10;
      tauCHIso6_ = -10;
      tauCHIso7_ = -10;
      tauCHIso8_ = -10;
      tauCHIso9_  = -10;
      tauCHIso10_  = -10;
      tauCHIso11_  = -10;
      tauCHIso12_  = -10;
      tauCHIso13_  = -10;
      tauCHIso14_  = -10;
      tauCHIso15_  = -10;
      tauCHIso16_  = -10;
      tauCHIso17_  = -10;
      tauCHIso18_  = -10;
      tauCHIso19_  = -10;
      tauCHIso20_  = -10;
      tauCHIsoOther_ = -10;

      tauCHIso0_ntracks_ = -10;
      tauCHIso1_ntracks_ = -10;
      tauCHIso2_ntracks_ = -10;     
      tauCHIso3_ntracks_ = -10;
      tauCHIso4_ntracks_ = -10;
      tauCHIso5_ntracks_ = -10;
      tauCHIso6_ntracks_ = -10;
      tauCHIso7_ntracks_ = -10;
      tauCHIso8_ntracks_ = -10;
      tauCHIso9_ntracks_ = -10;
      tauCHIso10_ntracks_ = -10;
      tauCHIso11_ntracks_ = -10;
      tauCHIso12_ntracks_ = -10;
      tauCHIso13_ntracks_ = -10;
      tauCHIso14_ntracks_ = -10;
      tauCHIso15_ntracks_ = -10;
      tauCHIso16_ntracks_ = -10;
      tauCHIso17_ntracks_ = -10;
      tauCHIso18_ntracks_ = -10;
      tauCHIso19_ntracks_ = -10;
      tauCHIso20_ntracks_ = -10;
      tauCHIsoOther_ntracks_ = -10;

      vtxX_ =  -10;
      vtxY_ =  -10;
      vtxZ_ =  -10;
      vtxT_ =  -10;

      taudXY_ =  -10;
      taudZ_ =  -10;
      taudT_ =  -10;

      vtxIndex_ = -10;
      pt_weighted_dr_signal_ = -10; 
      pt_weighted_deta_strip_ = -10;
      pt_weighted_dphi_strip_ = -10;
      pt_weighted_dr_iso_ = -10;    
      n_photons_total_ = -10;       
      nIsoTracks_ = -10;            
      nIsoNeutral_ = -10;           
      nIsoGamma_ = -10;             
      nSigGamma_ = -10;             
      nhIso_ = -10;                 
      sigPhoIso_ = -10;             
      dmf_=-10;
      goodReco_=0;
      genTauMatch_=0;
      genJetMatch_ = 0;
      genJetPt_ = -10;
      genJetEta_ = -10;
      jetPt_  = (float) jet.pt();
      jetEta_ = (float) jet.eta();
      if(jetPt_<18)
	continue;

      for (unsigned int iGenJet = 0; iGenJet < genJets->size() ; ++iGenJet){
	reco::GenJetRef genJet(genJets, iGenJet);
	if (reco::deltaR(genJet->eta(),genJet->phi(),jet.eta(),jet.phi()) < 0.4){
	  genJetMatch_ = 1;
	  genJetPt_ = genJet->pt();
	  genJetEta_ = genJet->eta();
	}
      }
      
      
     //std::cout<<" pt "<<genTauPt_<<" eta "<<genTauEta_<<std::endl;
      for(const pat::Tau &tau : *newTaus){
	if (reco::deltaR(tau.eta(),tau.phi(),jet.eta(),jet.phi()) < 0.5 && tau.tauID("decayModeFinding")>0.5){

	  tauPt_  = tau.pt();
	  tauEta_ = tau.eta();
	  dmf_ = tau.decayMode();
	  tauMass_ = tau.mass();

	  //get the matched vertex
	  int vtx_index = -1;
	  // find the 4D vertex this muon is best associated to..
	  float max_weight = 0.f;
	  for( unsigned i = 0; i < vtxs->size(); ++i ) {
	    const auto& vtx = (*vtxs)[i];      
	    const float weight = vtx.trackWeight(tau.leadChargedHadrCand());
	    if( weight > max_weight ) {
	      max_weight = weight;
	      vtx_index = i;
	    }
	  }
	  //now do vtx variable filling
	  vtxIndex_ = vtx_index;

	  const reco::Vertex& vtx = (vtx_index == -1 ? (*vtxs)[0] : (*vtxs)[vtx_index]);
	  vtxX_ = vtx.x();
	  vtxY_ = vtx.y();
	  vtxZ_ = vtx.z();
	  vtxT_ = vtx.t();

	  taudXY_ = tau.userFloat("tauDXY");
	  taudZ_ = tau.userFloat("tauDZ");
	  taudT_ = tau.userFloat("tauDT");

	  pt_weighted_dr_signal_ = tau.userFloat("tau_pt_weighted_dr_signal"); 
	  pt_weighted_deta_strip_ = tau.userFloat("tau_pt_weighted_deta_strip");
	  pt_weighted_dphi_strip_ = tau.userFloat("tau_pt_weighted_dphi_strip");
	  pt_weighted_dr_iso_ = tau.userFloat("tau_pt_weighted_dr_iso");	  
	  n_photons_total_ = tau.userFloat("n_photons_total");       
	  nIsoTracks_ = tau.userFloat("nIsoTracks");	  
	  nIsoNeutral_ = tau.userFloat("nIsoNeutral");           
	  nIsoGamma_ = tau.userFloat("nIsoGamma");	  
	  nSigGamma_ = tau.userFloat("nSigGamma");             
	  nhIso_ = tau.userFloat("nhIsoPt");	  
	  sigPhoIso_ = tau.userFloat("sigPhoIsoPt");             

	  tauCHIso_  = tau.userFloat("oldChIso");
	  tauCHIso0_ = tau.userFloat("newChIso0");
	  tauCHIso1_ = tau.userFloat("newChIso1");
	  tauCHIso2_ = tau.userFloat("newChIso2");
	  tauCHIso3_ = tau.userFloat("newChIso3");
	  tauCHIso4_ = tau.userFloat("newChIso4");
	  tauCHIso5_ = tau.userFloat("newChIso5");
	  tauCHIso6_ = tau.userFloat("newChIso6");
	  tauCHIso7_ = tau.userFloat("newChIso7");
	  tauCHIso8_ = tau.userFloat("newChIso8");
	  tauCHIso9_ = tau.userFloat("newChIso9");
	  tauCHIso10_ = tau.userFloat("newChIso10");
	  tauCHIso11_ = tau.userFloat("newChIso11");
	  tauCHIso12_ = tau.userFloat("newChIso12");
	  tauCHIso13_ = tau.userFloat("newChIso13");
	  tauCHIso14_ = tau.userFloat("newChIso14");
	  tauCHIso15_ = tau.userFloat("newChIso15");
	  tauCHIso16_ = tau.userFloat("newChIso16");
	  tauCHIso17_ = tau.userFloat("newChIso17");
	  tauCHIso18_ = tau.userFloat("newChIso18");
	  tauCHIso19_ = tau.userFloat("newChIso19");
	  tauCHIso20_ = tau.userFloat("newChIso20");
	  tauCHIsoOther_ = tau.userFloat("newChIsoOther");

	  //tauCHIso_ntracks_  = tau.userFloat("oldChIso");
	  tauCHIso0_ntracks_ = tau.userFloat("newChIso0_ntracks");
	  tauCHIso1_ntracks_ = tau.userFloat("newChIso1_ntracks");
	  tauCHIso2_ntracks_ = tau.userFloat("newChIso2_ntracks");
	  tauCHIso3_ntracks_ = tau.userFloat("newChIso3_ntracks");
	  tauCHIso4_ntracks_ = tau.userFloat("newChIso4_ntracks");
	  tauCHIso5_ntracks_ = tau.userFloat("newChIso5_ntracks");
	  tauCHIso6_ntracks_ = tau.userFloat("newChIso6_ntracks");
	  tauCHIso7_ntracks_ = tau.userFloat("newChIso7_ntracks");
	  tauCHIso8_ntracks_ = tau.userFloat("newChIso8_ntracks");
	  tauCHIso9_ntracks_ = tau.userFloat("newChIso9_ntracks");
	  tauCHIso10_ntracks_ = tau.userFloat("newChIso10_ntracks");
	  tauCHIso11_ntracks_ = tau.userFloat("newChIso11_ntracks");
	  tauCHIso12_ntracks_ = tau.userFloat("newChIso12_ntracks");
	  tauCHIso13_ntracks_ = tau.userFloat("newChIso13_ntracks");
	  tauCHIso14_ntracks_ = tau.userFloat("newChIso14_ntracks");
	  tauCHIso15_ntracks_ = tau.userFloat("newChIso15_ntracks");
	  tauCHIso16_ntracks_ = tau.userFloat("newChIso16_ntracks");
	  tauCHIso17_ntracks_ = tau.userFloat("newChIso17_ntracks");
	  tauCHIso18_ntracks_ = tau.userFloat("newChIso18_ntracks");
	  tauCHIso19_ntracks_ = tau.userFloat("newChIso19_ntracks");
	  tauCHIso20_ntracks_ = tau.userFloat("newChIso20_ntracks");
	  tauCHIsoOther_ntracks_ = tau.userFloat("newChIsoOther_ntracks");

	  goodReco_ = (bool) tau.tauID(tauID_) >0.5;
	  //std::cout<<"iso0: "<<tauCHIso0_<<" iso1: "<<tauCHIso1_<<" iso8: "<<tauCHIso8_<<" iso: "<<tauCHIso_<<" isoOther: "<<tauCHIsoOther_<<std::endl;
	  break;
	}
      }

      if(tauPt_ > 10)
	jetTree->Fill(); 
   }

}


// ------------ method called once each job just before starting event loop  ------------
void 
timingTausMiniAOD::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
timingTausMiniAOD::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
timingTausMiniAOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
}


//Gets visible 4-momentum of a particle from list of daughters
reco::Candidate::LorentzVector timingTausMiniAOD::GetVisibleP4(std::vector<const reco::GenParticle*>& daughters){
  reco::Candidate::LorentzVector p4_vis(0,0,0,0);
  for(size_t i = 0; i < daughters.size(); ++i){
    if (!isNeutrino(daughters[i]) && daughters[i]->status() == 1){  //list should include intermediate daughters, so check for status = 1
      p4_vis += daughters[i]->p4();
    }
  }
  return p4_vis;
}


//Creates a vector of all (including intermediate) daughters for a given mother particle
void timingTausMiniAOD::findDaughters(const reco::GenParticle* mother, std::vector<const reco::GenParticle*>& daughters)
{
  unsigned numDaughters = mother->numberOfDaughters();
  if (numDaughters == 0) std::cout << " none ";
  for (unsigned iDaughter = 0; iDaughter < numDaughters; ++iDaughter ) {
    const reco::GenParticle* daughter = mother->daughterRef(iDaughter).get();
    if (daughter->status() == 1){  //status = 1 is a final state daughter
      daughters.push_back(daughter); 
    }
    if (daughter->status() == 2){  //status = 2 is an intermediate daughter; will decay further
      daughters.push_back(daughter); 
      findDaughters(daughter, daughters);
    }
  }
}

bool timingTausMiniAOD::isNeutrino(const reco::Candidate* daughter)
{
  return (TMath::Abs(daughter->pdgId()) == 12 || TMath::Abs(daughter->pdgId()) == 14 || TMath::Abs(daughter->pdgId()) == 16 || TMath::Abs(daughter->pdgId()) == 18);
}

double timingTausMiniAOD::getTimedCHIsoSum( pat::PackedCandidate const* leadTauCand, std::vector<pat::PackedCandidate const*> isolationCands, float interval, float & n_tracks)
{
  float sum = 0;
  n_tracks = 0;
  for( auto cand : isolationCands ) {
    //std::cout<<"diff time: "<<cand->time()*10-leadTauCand->time()*10<<std::endl;
    //if(fabs(cand->dtimeAssociatedPV()*10-leadTauCand->dtimeAssociatedPV()*10) < interval){
    if(std::abs(cand->time()-leadTauCand->time()) < interval){
      sum += cand->pt();
      n_tracks +=1;
    }
    /*
    if(fabs(cand->dtimeAssociatedPV()-leadTauCand->dtimeAssociatedPV())/0.0010219 < interval){
      sum += cand->pt();
      n_tracks +=1;
      }*/
  }

  return sum;
}


//define this as a plug-in
DEFINE_FWK_MODULE(timingTausMiniAOD);
