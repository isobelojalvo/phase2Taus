// -*- C++ -*-
//
// Package:    RecoTauTag/phase2TausRECO
// Class:      phase2TausRECO
// 
/**\class phase2TausRECO phase2TausRECO.cc RecoTauTag/phase2TausRECO/plugins/phase2TausRECO.cc

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
#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// user include files
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "DataFormats/Candidate/interface/CompositePtrCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Math/interface/Point3D.h"

#include "DataFormats/Provenance/interface/EventAuxiliary.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTau.h"

#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TMath.h"
#include <unordered_map>

using namespace reco;
using std::vector;
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

struct customPFTau{
  reco::PFTau pfTau;
  float chargedIso;
  float signalTracksMatchTimeSlice;
  int goodThreeProngTauT3;
  int goodThreeProngTauT4;
  std::vector<reco::PFCandidatePtr> signalCharged_;
  std::vector<reco::PFCandidatePtr> isoCharged_;
  std::vector<reco::PFCandidatePtr> isoChargedT1_;
  std::vector<reco::PFCandidatePtr> isoChargedT2_;
};

class phase2TausRECO : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit phase2TausRECO(const edm::ParameterSet&);
      ~phase2TausRECO();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  //edm::EDGetTokenT<std::vector<reco::Track> > tracksTag_;
  edm::InputTag  tracksTag_;
  edm::InputTag  timesTag_;
  edm::InputTag  timeResosTag_;
  edm::EDGetTokenT<reco::VertexCollection>   vtxToken_;
  edm::EDGetTokenT<reco::PFTauCollection>    tauToken_;
  edm::EDGetTokenT<reco::PFJetCollection>    jetSrc_;
  edm::EDGetTokenT<reco::GenJetCollection>   genJetSrc_;
  edm::EDGetTokenT<reco::PFTauDiscriminator> discriminatorSrc_;
  edm::EDGetTokenT<reco::PFTauDiscriminator> pfChargedSrc_;
  edm::EDGetTokenT<reco::PFTauDiscriminator> dmfToken_;
  edm::EDGetTokenT<std::vector <reco::GenParticle> > genToken_;

  std::string tauID_;

  TTree* tree;
  TTree* jetTree;
  double tauPt_;
  double tauEta_;
  double tauPhi_;
  double taudXY_, taudZ_, taudT_;
  double tauMass_;
  double genTauPt_;
  double genTauEta_;
  double genTauPhi_;
  double jetPt_;
  double jetEta_;
  double vtxX_, vtxY_, vtxZ_, vtxT_;
  double tauCHIso_;
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

  int genMuonMatch_;
  int genElectronMatch_;

  double genJetPt_;
  double genJetEta_;
  double genJetPhi_;

  double isAMuon_;
  double isAElectron_;


  int nvtx_;
  int dmf_;
  int goodReco_;
  int genTauMatch_;
  int jetTauMatch_;
  int genJetMatch_;
  int vtxIndex_;
  int good3ProngT3_;
  int good3ProngT4_;
  bool cutByDiscriminator_;

  reco::Candidate::LorentzVector GetVisibleP4(std::vector<const reco::GenParticle*>& daughters);
  void findDaughters(const reco::GenParticle* mother, std::vector<const reco::GenParticle*>& daughters);
  bool isNeutrino(const reco::Candidate* daughter);
  double getTimedCHIsoSum( const PFCandidatePtr&  leadTauCand, vector<PFCandidatePtr>& isolationCands, float interval, float & n_tracks, reco::Vertex pv);
  /*
  bool initializeTrackTiming( edm::Handle<edm::View<reco::Track> > tracks, 
			      edm::Handle<reco::VertexCollection> &vtxs, 
			      edm::Handle<edm::ValueMap<float> > times, 
			      edm::Handle<edm::ValueMap<float> > timeResos,
			      std::unordered_multimap<unsigned,reco::TrackBaseRef> &vertices_to_tracks_z,
			      std::unordered_multimap<unsigned,reco::TrackBaseRef> &vertices_to_tracks_zt1, 
			      std::unordered_multimap<unsigned,reco::TrackBaseRef> &vertices_to_tracks_zt2, 
			      std::unordered_multimap<unsigned,reco::TrackBaseRef> &vertices_to_tracks_zt3, 
			      std::unordered_multimap<unsigned,reco::TrackBaseRef> &vertices_to_tracks_zt4, 
			      std::unordered_multimap<unsigned,reco::TrackBaseRef> &vertices_to_tracks_zt5, 
			      std::unordered_multimap<unsigned,reco::TrackBaseRef> &vertices_to_tracks_zt6
			      );

  void checkThreeProngTau(customPFTau tau,
			  int vtx_index,
			  std::unordered_multimap<unsigned,reco::TrackBaseRef> &vertices_to_tracks_zt3,
			  std::unordered_multimap<unsigned,reco::TrackBaseRef> &vertices_to_tracks_zt4,
			  bool &matchT3,
			  bool &matchT4);

void calculateIsoQuantities(customPFTau tau,
			    int vtx_index,
			    std::unordered_multimap<unsigned,reco::TrackBaseRef> &vertices_to_tracks_z, 
			    std::unordered_multimap<unsigned,reco::TrackBaseRef> &vertices_to_tracks_zt1, 
			    std::unordered_multimap<unsigned,reco::TrackBaseRef> &vertices_to_tracks_zt2, 
			    std::unordered_multimap<unsigned,reco::TrackBaseRef> &vertices_to_tracks_zt3, 
			    std::unordered_multimap<unsigned,reco::TrackBaseRef> &vertices_to_tracks_zt4, 
			    std::unordered_multimap<unsigned,reco::TrackBaseRef> &vertices_to_tracks_zt5, 
			    std::unordered_multimap<unsigned,reco::TrackBaseRef> &vertices_to_tracks_zt6,
			    double &ptSum, 
			    double &ptSumT1, 
			    double &ptSumT2, 
			    double &ptSumT3, 
			    double &ptSumT4, 
			    double &ptSumT5, 
			    double &ptSumT6);
  */
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
phase2TausRECO::phase2TausRECO(const edm::ParameterSet& iConfig):
  //tracksTag_       (consumes<std::vector<reco::Track> >  (iConfig.getParameter<edm::InputTag>("tracks"))),
  tracksTag_       ("generalTracks",""),
  timesTag_        ("trackTimeValueMapProducer","generalTracksConfigurableFlatResolutionModel"),
  timeResosTag_    ("trackTimeValueMapProducer","generalTracksConfigurableFlatResolutionModelResolution"),
  vtxToken_        (consumes<reco::VertexCollection>    (iConfig.getParameter<edm::InputTag>("vertices")      )),
  tauToken_        (consumes<reco::PFTauCollection>     (iConfig.getParameter<edm::InputTag>("taus")          )),
  jetSrc_          (consumes<reco::PFJetCollection>     (iConfig.getParameter<edm::InputTag>("jets")          )),
  genJetSrc_       (consumes<reco::GenJetCollection>    (iConfig.getParameter<edm::InputTag>("genJets")       )),
  discriminatorSrc_(consumes<reco::PFTauDiscriminator>  (iConfig.getParameter<edm::InputTag>("discriminator") )),
  pfChargedSrc_    (consumes<reco::PFTauDiscriminator>  (iConfig.getParameter<edm::InputTag>("hpsPFTauChargedIsoPtSum"))),
  dmfToken_        (consumes<reco::PFTauDiscriminator>  (iConfig.getParameter<edm::InputTag>("dmf") )),
  genToken_  (consumes<std::vector<reco::GenParticle> > (iConfig.getParameter<edm::InputTag>("genParticles")  ))
{
  consumes<edm::View<reco::Track> >(tracksTag_);
  consumes<edm::ValueMap<float> >(timesTag_);
  consumes<edm::ValueMap<float> >(timeResosTag_);
  cutByDiscriminator_ = true;
  cutByDiscriminator_ = iConfig.getUntrackedParameter<bool>("cutByDiscriminator");

  //now do what ever initialization is needed
  usesResource("TFileService");

   //tauID_    = iConfig.getParameter<std::string>("tauID");
   edm::Service<TFileService> fs;

   tree = fs->make<TTree>(      "Ntuple",      "Ntuple"          );
   tree->Branch("tauPt",        &tauPt_,       "tauPt_/D"        );
   tree->Branch("tauEta",       &tauEta_,      "tauEta_/D"       );
   tree->Branch("genTauPt",     &genTauPt_,    "genTauPt_/D"     );
   tree->Branch("genTauEta",    &genTauEta_,   "genTauEta_/D"    );
   tree->Branch("genTauPhi",   &genTauPhi_,  "genTauPhi/D"  );
   tree->Branch("genTauMatch",  &genTauMatch_, "genTauMatch_/I"  );
   tree->Branch("good3ProngT3", &good3ProngT3_,"good3ProngT3__/I");
   tree->Branch("good3ProngT4", &good3ProngT4_,"good3ProngT4__/I");
   tree->Branch("nvtx",         &nvtx_,        "nvtx_/I"         );
   tree->Branch("vtxX",         &vtxX_,        "vtxX_/D"         );
   tree->Branch("vtxY",         &vtxY_,        "vtxY_/D"         );
   tree->Branch("vtxZ",         &vtxZ_,        "vtxZ_/D"         );
   tree->Branch("vtxT",         &vtxT_,        "vtxT_/D"         );
   tree->Branch("vtxIndex",     &vtxIndex_,    "vtxIndex_/I"     );
   tree->Branch("taudXY",      &taudXY_,      "taudXY/D"       );
   tree->Branch("taudZ",       &taudZ_,       "taudZ/D"        );
   tree->Branch("taudT",       &taudT_,       "taudT/D"        );
   tree->Branch("dmf",          &dmf_,         "dmf_/I"          );
   tree->Branch("goodReco",     &goodReco_,    "goodReco_/I"     );
   tree->Branch("tauMass",      &tauMass_,     "tauMass_/D"      );
   tree->Branch("tauCHIso",    &tauCHIso_,   "tauCHIso_/D"    );
   tree->Branch("tauCHIso1",  &tauCHIso1_, "tauCHIso1_/D"  );
   tree->Branch("tauCHIso2",  &tauCHIso2_, "tauCHIso2_/D"  );
   tree->Branch("tauCHIso3",  &tauCHIso3_, "tauCHIso3_/D"  );
   tree->Branch("tauCHIso4",  &tauCHIso4_, "tauCHIso4_/D"  );
   tree->Branch("tauCHIso5",  &tauCHIso5_, "tauCHIso5_/D"  );
   tree->Branch("tauCHIso6",  &tauCHIso6_, "tauCHIso6_/D"  );
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


   jetTree = fs->make<TTree>(      "jetNtuple",   "jetNtuple"       );
   jetTree->Branch("tauPt",        &tauPt_,       "tauPt_/D"        );
   jetTree->Branch("tauEta",       &tauEta_,      "tauEta_/D"       );
   jetTree->Branch("jetPt",        &jetPt_,       "jetPt_/D"        );
   jetTree->Branch("jetEta",       &jetEta_,      "jetEta_/D"       );
   jetTree->Branch("jetTauMatch",  &jetTauMatch_, "jetTauMatch_/I"  );
   jetTree->Branch("genJetMatch",  &genJetMatch_, "genJetMatch_/I"  );
   jetTree->Branch("genJetPt",     &genJetPt_,    "genJetPt_/D"  );
   jetTree->Branch("genJetEta",    &genJetEta_,   "genJetEta_/D"  );
   jetTree->Branch("genJetPhi",    &genJetPhi_,   "genJetPhi_/D"  );
   jetTree->Branch("genMuonMatch", &genMuonMatch_, "genMuonMatch/I"  );
   jetTree->Branch("genElectronMatch", &genElectronMatch_, "genElectronMatch/I"  );
   jetTree->Branch("isAMuon",     &isAMuon_,     "isAMuon_/D"     );
   jetTree->Branch("isAElectron", &isAElectron_, "isAElectron_/D" );
   jetTree->Branch("nvtx",         &nvtx_,        "nvtx_/I"         );
   jetTree->Branch("vtxX",         &vtxX_,        "vtxX_/D"         );
   jetTree->Branch("vtxY",         &vtxY_,        "vtxY_/D"         );
   jetTree->Branch("vtxZ",         &vtxZ_,        "vtxZ_/D"         );
   jetTree->Branch("vtxT",         &vtxT_,        "vtxT_/D"         );
   jetTree->Branch("vtxIndex",     &vtxIndex_,    "vtxIndex_/I"     );
   jetTree->Branch("dmf",          &dmf_,         "dmf_/I"          );
   jetTree->Branch("goodReco",     &goodReco_,    "goodReco_/I"     );
   jetTree->Branch("tauMass",      &tauMass_,     "tauMass_/D"      );
   jetTree->Branch("tauCHIso",      &tauCHIso_ ,     "tauCHIso/D"  );
   // jetTree->Branch("tauCHIso0",     &tauCHIso0_,     "tauCHIso0/D" );
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

   //jetTree->Branch("tauCHIso0_ntracks",     &tauCHIso0_ntracks_,     "tauCHIso0_ntracks/D" );
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

}


phase2TausRECO::~phase2TausRECO()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
phase2TausRECO::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   //get track colletions

   edm::Handle<edm::View<reco::Track> >tracks;
   if(!iEvent.getByLabel(tracksTag_,tracks))
     std::cout<<"Error getting Tracks"<<std::endl;

   edm::Handle<edm::ValueMap<float> > times;
   if(!iEvent.getByLabel(timesTag_,times))
     std::cout<<"Error getting Times"<<std::endl;

   edm::Handle<edm::ValueMap<float> > timeResos;
   if(!iEvent.getByLabel(timeResosTag_,timeResos))
     std::cout<<"Error getting TimeResos"<<std::endl;

   Handle<reco::VertexCollection> vtxs;
   if(!iEvent.getByToken(vtxToken_, vtxs))
     std::cout<<"Error getting vtxs"<<std::endl;
   nvtx_ = vtxs->size();

   //get objects
   Handle<reco::PFJetCollection> JetObjects;
   if(!iEvent.getByToken(jetSrc_, JetObjects))
     std::cout<<"Error getting Jets"<<std::endl;

   edm::Handle<reco::PFTauCollection> taus;
   if(!iEvent.getByToken(tauToken_, taus))
     std::cout<<"Error getting Taus"<<std::endl;

   edm::Handle<reco::GenJetCollection> genJets;
   if(!iEvent.getByToken(genJetSrc_, genJets))
     std::cout<<"Error getting genJets"<<std::endl;

   Handle<reco::PFTauDiscriminator> discriminator;
   if(!iEvent.getByToken(discriminatorSrc_, discriminator))
     std::cout<<"Error getting tau discriminator"<<std::endl;

   Handle<reco::PFTauDiscriminator> chargedDiscriminator;
   if(!iEvent.getByToken(pfChargedSrc_, chargedDiscriminator))
     std::cout<<"Error getting Tau charged Iso"<<std::endl;


   Handle<reco::PFTauDiscriminator> DMF; 
   //iEvent.getByLabel("hpsPFTauDiscriminationByDecayModeFindingOldDMs",DMF);
   if(!iEvent.getByToken(dmfToken_,DMF))
     std::cout<<"Error getting DMF disc"<<std::endl;

   //std::vector<const reco::GenParticle*> GenObjects = getGenParticleCollection(iEvent);
   //if (GenObjects.size()!=0) return;

   edm::Handle<std::vector<reco::GenParticle> > genParticles;
   if(!iEvent.getByToken(genToken_, genParticles))
     std::cout<<"Error getting Gen Particles"<<std::endl;

   std::vector<const reco::GenParticle*> GenTaus;
   std::vector<const reco::GenParticle*> GenElectrons;
   std::vector<const reco::GenParticle*> GenMuons;
   for(std::vector<reco::GenParticle>::const_iterator genParticle = genParticles->begin(); genParticle != genParticles->end(); genParticle++){
     if(TMath::Abs(genParticle->pdgId()) == 15) GenTaus.push_back(&(*genParticle));
     if(TMath::Abs(genParticle->pdgId()) == 11) GenElectrons.push_back(&(*genParticle));
     if(TMath::Abs(genParticle->pdgId()) == 13) GenMuons.push_back(&(*genParticle));
   }

   std::vector<reco::PFJet> Jets;

   for (unsigned int iJet = 0; iJet < JetObjects->size() ; ++iJet){
     reco::PFJetRef jetCand(JetObjects, iJet);
     if(jetCand->pt() < 18 )continue;
     bool isATau = false;
     bool isAMuon = false;
     bool isAElectron = false;
     for(auto genTau : GenTaus){
       std::vector<const reco::GenParticle*> genTauDaughters;
       findDaughters(genTau, genTauDaughters);
       reco::Candidate::LorentzVector genTauVis = GetVisibleP4(genTauDaughters);
       genTauPt_  = (float) genTauVis.pt();
       genTauEta_ = (float) genTauVis.eta();
       if (reco::deltaR(jetCand->eta(),jetCand->phi(),genTauVis.eta(),genTauVis.phi()) < 0.4)
	 isATau=true;

     }
     for(auto genMuon: GenMuons)
       if (reco::deltaR(jetCand->eta(),jetCand->phi(),genMuon->eta(),genMuon->phi()) < 0.4)
	 isAMuon = true;

     for(auto genElectron: GenElectrons)
       if (reco::deltaR(jetCand->eta(),jetCand->phi(),genElectron->eta(),genElectron->phi()) < 0.4)
	 isAElectron = true;

     if(!isATau&&!isAMuon&&!isAElectron)
       Jets.push_back(*jetCand);
   }   

   /*
   std::vector<customPFTau> goodTaus;
   for (unsigned int iTau = 0; iTau<taus->size() ; ++iTau){
     reco::PFTauRef tauCandidate(taus, iTau);
     if((*DMF)[tauCandidate] > -1){
       if(cutByDiscriminator_){
	 if((*discriminator)[tauCandidate] > 0.5){
	   customPFTau tempCPFTau;
	   tempCPFTau.pfTau = *tauCandidate;
	   tempCPFTau.chargedIso = (*chargedDiscriminator)[tauCandidate];
	   ///check three prong tau code
	   //get the matched vertex
	   int vtx_index = -1;
	   // find the 4D vertex this tau is best associated to..
	   float max_weight = 0.f;
	   
	   for( unsigned i = 0; i < vtxs->size(); ++i ) {
	     const auto& vtx = (*vtxs)[i];      
	     const float weight = vtx.trackWeight(tempCPFTau.pfTau.leadPFChargedHadrCand()->trackRef());// check me -> Get track ref for charged hadron candidate
	     if( weight > max_weight ) {
	       max_weight = weight;
	       vtx_index = i;
	     }
	   }
	   
	   //check if 3 prong tau tracks all come from the same time slice
	   bool goodThreeProngT3 = false;
	   bool goodThreeProngT4 = false;       
	   
	   if(tempCPFTau.pfTau.decayMode()==10)
	 checkThreeProngTau(tempCPFTau,vtx_index,vertices_to_tracks_zt3,vertices_to_tracks_zt4,goodThreeProngT3,goodThreeProngT4);
	   //////// finish three prong tau code
	   
	   for( auto const & cand : tauCandidate->isolationPFChargedHadrCands() ){
	     tempCPFTau.isoCharged_.push_back(cand);}
	   
	   goodTaus.push_back(tempCPFTau);
	 }
       }
       else{
	 if((*chargedDiscriminator)[tauCandidate] > -1){
	   customPFTau tempCPFTau;
	   tempCPFTau.pfTau = *tauCandidate;
	   tempCPFTau.chargedIso = (*chargedDiscriminator)[tauCandidate];
	   
	   for( auto const & cand : tauCandidate->isolationPFChargedHadrCands() ) 
	     tempCPFTau.isoCharged_.push_back(cand);


	   ///check three prong tau code
	   //get the matched vertex
	   int vtx_index = -1;
	   // find the 4D vertex this tau is best associated to..
	   float max_weight = 0.f;
	   
	   for( unsigned i = 0; i < vtxs->size(); ++i ) {
	     const auto& vtx = (*vtxs)[i];      
	     const float weight = vtx.trackWeight(tempCPFTau.pfTau.leadPFChargedHadrCand()->trackRef());// check me -> Get track ref for charged hadron candidate
	     if( weight > max_weight ) {
	       max_weight = weight;
	       vtx_index = i;
	     }
	   }
	   
	   //check if 3 prong tau tracks all come from the same time slice
	   bool goodThreeProngT3 = false;
	   bool goodThreeProngT4 = false;       
	   
	   if(tempCPFTau.pfTau.decayMode()==10)
	 checkThreeProngTau(tempCPFTau,vtx_index,vertices_to_tracks_zt3,vertices_to_tracks_zt4,goodThreeProngT3,goodThreeProngT4);
	   //////// finish three prong tau code


	   
	   goodTaus.push_back(tempCPFTau);
	 }
       }
       //std::vector<reco::PFCandidatePtr> isoCharged_;

     }
   }

*/
   for(auto genTau : GenTaus){
      tauPt_=-10;
      tauEta_=-10;
      tauMass_=-10;
      genTauPt_=-10;
      genTauEta_=-10;
     
      //nvtx_=-10;
      dmf_=-10;
      goodReco_=0;
      genTauMatch_=0;
      tauCHIso_=0;
      tauCHIso1_=0;
      tauCHIso2_=0;
      tauCHIso3_=0;
      tauCHIso4_=0;
      tauCHIso5_=0;
      tauCHIso6_=0;
      tauCHIso7_=0;
      tauCHIso8_=0;
      tauCHIso9_=0;
      tauCHIso10_=0;
      tauCHIso11_=0;
      tauCHIso12_=0;
      tauCHIso13_=0;
      tauCHIso14_=0;
      tauCHIso15_=0;
      tauCHIso16_=0;
      tauCHIso17_=0;
      tauCHIso18_=0;
      tauCHIso19_=0;
      tauCHIso20_=0;
      tauCHIsoOther_=0;

      tauCHIso1_ntracks_ = 0;
      tauCHIso2_ntracks_ = 0;
      tauCHIso3_ntracks_ = 0;
      tauCHIso4_ntracks_ = 0;
      tauCHIso5_ntracks_ = 0;
      tauCHIso6_ntracks_ = 0;
      tauCHIso7_ntracks_ = 0;
      tauCHIso8_ntracks_ = 0;
      tauCHIso9_ntracks_ = 0;
      tauCHIso10_ntracks_ =0;
      tauCHIso11_ntracks_ =0;
      tauCHIso12_ntracks_ =0;
      tauCHIso13_ntracks_ =0;
      tauCHIso14_ntracks_ =0;
      tauCHIso15_ntracks_ =0;
      tauCHIso16_ntracks_ =0;
      tauCHIso17_ntracks_ =0;
      tauCHIso18_ntracks_ =0;
      tauCHIso19_ntracks_ =0;
      tauCHIso20_ntracks_ =0;
      tauCHIsoOther_ntracks_ =0;

      std::vector<const reco::GenParticle*> genTauDaughters;
      findDaughters(genTau, genTauDaughters);
      reco::Candidate::LorentzVector genTauVis = GetVisibleP4(genTauDaughters);
      genTauPt_  = (float) genTauVis.pt();
      genTauEta_ = (float) genTauVis.eta();
      genTauPhi_ = (float) genTauVis.phi();
      float minDR = 20;

      //std::cout<<" pt "<<genTauPt_<<" eta "<<genTauEta_<<std::endl;
      for (unsigned int iTau = 0; iTau<taus->size() ; ++iTau){
	reco::PFTauRef tauCandidate(taus, iTau);
	reco::PFTau pfTau= *tauCandidate;
	if((*DMF)[tauCandidate] > -1){

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
	  
	  //for(auto tau : goodTaus){
	  if (reco::deltaR(pfTau.eta(),pfTau.phi(),genTauVis.eta(),genTauVis.phi()) < 0.5){
	    genTauMatch_ = 1;
	    tauPt_  = pfTau.pt();
	    tauEta_ = pfTau.eta();
	    dmf_    = pfTau.decayMode();
	    //tauCHIso_ = (*chargedDiscriminator)[tauCandidate];
	    
	    //get the matched vertex
	    int vtx_index = -1;
	    // find the 4D vertex this muon is best associated to..
	    float max_weight = 0.f;
	  for( unsigned i = 0; i < vtxs->size(); ++i ) {
	    const auto& vtx = (*vtxs)[i];      
	    const float weight = vtx.trackWeight(pfTau.leadPFChargedHadrCand()->trackRef());// check me -> Get track ref for charged hadron candidate
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
	  const PFCandidatePtr& packedLeadTauCand = pfTau.leadPFChargedHadrCand();
	  //const std::vector<PFRecoTauChargedHadron>& pfCHIsoTauCands = pfTau.isolationPFChargedHadrCands();
	  std::vector<PFCandidatePtr> packedCHIsoTauCands;
	  //std::cout<<"pfTau.isolationTauChargedHadronCandidates size: "<<pfTau.isolationPFChargedHadrCands().size()<<std::endl;
	  for(auto isocand :  pfTau.isolationPFChargedHadrCands()){
	    packedCHIsoTauCands.push_back(isocand);
	    tauCHIso_ += isocand->pt();
	  }

	  //const float dz = std::abs(packedLeadTauCand->dz(vtxs->at(0).position()));
	  //const float dxy = std::abs(packedLeadTauCand->dxy(vtxs->at(0).position()));
	  //const float dt = packedLeadTauCand->dtime(0);


	  const double res_det = 0.030; // ns

	  //std::cout<<"packedleadtaucand pt: "<<packedLeadTauCand->pt()<<" ncands: "<<packedCHIsoTauCands.size()<<std::endl;
	  newChIso0  = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 0*res_det, newChIso0_ntracks, vtxs->at(0));
	  newChIso1  = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 1*res_det, newChIso1_ntracks, vtxs->at(0));
	  newChIso2  = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 2*res_det, newChIso2_ntracks, vtxs->at(0));
	  newChIso3  = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 3*res_det, newChIso3_ntracks, vtxs->at(0)); // this is the nominal cut, 3 sigma
	  newChIso4  = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 4*res_det, newChIso4_ntracks, vtxs->at(0));
	  newChIso5  = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 5*res_det, newChIso5_ntracks, vtxs->at(0));
	  newChIso6  = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 6*res_det, newChIso6_ntracks, vtxs->at(0));
	  newChIso7  = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 7*res_det, newChIso7_ntracks, vtxs->at(0));
	  newChIso8  = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 8*res_det, newChIso8_ntracks, vtxs->at(0));
	  newChIso9  = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 9*res_det, newChIso9_ntracks, vtxs->at(0));
	  newChIso10 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 10*res_det, newChIso10_ntracks, vtxs->at(0));
	  newChIso11 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 11*res_det, newChIso11_ntracks, vtxs->at(0));
	  newChIso12 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 12*res_det, newChIso12_ntracks, vtxs->at(0));
	  newChIso13 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 13*res_det, newChIso13_ntracks, vtxs->at(0));
	  newChIso14 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 14*res_det, newChIso14_ntracks, vtxs->at(0));
	  newChIso15 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 15*res_det, newChIso15_ntracks, vtxs->at(0));
	  newChIso16 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 16*res_det, newChIso16_ntracks, vtxs->at(0));
	  newChIso17 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 17*res_det, newChIso17_ntracks, vtxs->at(0));
	  newChIso18 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 18*res_det, newChIso18_ntracks, vtxs->at(0));
	  newChIso19 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 19*res_det, newChIso19_ntracks, vtxs->at(0));
	  newChIso20 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 20*res_det, newChIso20_ntracks, vtxs->at(0));
	  newChIsoOther = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 100*res_det, newChIsoOther_ntracks, vtxs->at(0));

	  tauCHIso1_ = newChIso1;
	  tauCHIso2_ = newChIso2 ; 
	  tauCHIso3_ = newChIso3 ; 
	  tauCHIso4_ = newChIso4 ; 
	  tauCHIso5_ = newChIso5 ; 
	  tauCHIso6_ = newChIso6 ; 
	  tauCHIso7_ = newChIso7 ; 
	  tauCHIso8_ = newChIso8 ; 
	  tauCHIso9_ = newChIso9 ; 
	  tauCHIso10_ = newChIso10; 
	  tauCHIso11_ = newChIso11; 
	  tauCHIso12_ = newChIso12; 
	  tauCHIso13_ = newChIso13; 
	  tauCHIso14_ = newChIso14; 
	  tauCHIso15_ = newChIso15; 
	  tauCHIso16_ = newChIso16; 
	  tauCHIso17_ = newChIso17; 
	  tauCHIso18_ = newChIso18; 
	  tauCHIso19_ = newChIso19; 

	  tauCHIso1_ntracks_ = newChIso1_ntracks;
	  tauCHIso2_ntracks_ = newChIso2_ntracks; 
	  tauCHIso3_ntracks_ = newChIso3_ntracks; 
	  tauCHIso4_ntracks_ = newChIso4_ntracks; 
	  tauCHIso5_ntracks_ = newChIso5_ntracks; 
	  tauCHIso6_ntracks_ = newChIso6_ntracks; 
	  tauCHIso7_ntracks_ = newChIso7_ntracks; 
	  tauCHIso8_ntracks_ = newChIso8_ntracks; 
	  tauCHIso9_ntracks_ = newChIso9_ntracks; 
	  tauCHIso10_ntracks_ = newChIso10_ntracks; 
	  tauCHIso11_ntracks_ = newChIso11_ntracks; 
	  tauCHIso12_ntracks_ = newChIso12_ntracks; 
	  tauCHIso13_ntracks_ = newChIso13_ntracks; 
	  tauCHIso14_ntracks_ = newChIso14_ntracks; 
	  tauCHIso15_ntracks_ = newChIso15_ntracks; 
	  tauCHIso16_ntracks_ = newChIso16_ntracks; 
	  tauCHIso17_ntracks_ = newChIso17_ntracks; 
	  tauCHIso18_ntracks_ = newChIso18_ntracks; 
	  tauCHIso19_ntracks_ = newChIso19_ntracks; 

	  break;
	  }
	}
      }
      tree->Fill(); 
   }
   

   for(auto jet : Jets){
     if(jet.pt()<18)
       continue;
     tauPt_=-10;
     tauEta_=-10;
     tauMass_=-10;
     jetPt_=jet.pt();
     jetEta_=jet.eta();
     genJetMatch_ = 0;     
     //nvtx_=-10;
     dmf_=-10;
     goodReco_=0;
     jetTauMatch_=0;

      tauCHIso_=0;
      tauCHIso1_=0;
      tauCHIso2_=0;
      tauCHIso3_=0;
      tauCHIso4_=0;
      tauCHIso5_=0;
      tauCHIso6_=0;
      tauCHIso7_=0;
      tauCHIso8_=0;
      tauCHIso9_=0;
      tauCHIso10_=0;
      tauCHIso11_=0;
      tauCHIso12_=0;
      tauCHIso13_=0;
      tauCHIso14_=0;
      tauCHIso15_=0;
      tauCHIso16_=0;
      tauCHIso17_=0;
      tauCHIso18_=0;
      tauCHIso19_=0;
      tauCHIso20_=0;
      tauCHIsoOther_=0;

      tauCHIso1_ntracks_ = 0;
      tauCHIso2_ntracks_ = 0;
      tauCHIso3_ntracks_ = 0;
      tauCHIso4_ntracks_ = 0;
      tauCHIso5_ntracks_ = 0;
      tauCHIso6_ntracks_ = 0;
      tauCHIso7_ntracks_ = 0;
      tauCHIso8_ntracks_ = 0;
      tauCHIso9_ntracks_ = 0;
      tauCHIso10_ntracks_ =0;
      tauCHIso11_ntracks_ =0;
      tauCHIso12_ntracks_ =0;
      tauCHIso13_ntracks_ =0;
      tauCHIso14_ntracks_ =0;
      tauCHIso15_ntracks_ =0;
      tauCHIso16_ntracks_ =0;
      tauCHIso17_ntracks_ =0;
      tauCHIso18_ntracks_ =0;
      tauCHIso19_ntracks_ =0;
      tauCHIso20_ntracks_ =0;
      tauCHIsoOther_ntracks_ =0;
     
      int vtx_index = -1;
      const reco::Vertex& vtx = (vtx_index == -1 ? (*vtxs)[0] : (*vtxs)[vtx_index]);
      vtxX_ = vtx.x();
      vtxY_ = vtx.y();
      vtxZ_ = vtx.z();
      vtxT_ = vtx.t();
      
      for (unsigned int iGenJet = 0; iGenJet < genJets->size() ; ++iGenJet){
	reco::GenJetRef genJet(genJets, iGenJet);
	if (reco::deltaR(genJet->eta(),genJet->phi(),jet.eta(),jet.phi()) < 0.4)
	  genJetMatch_ = 1;
      }
      
      for (unsigned int iTau = 0; iTau<taus->size() ; ++iTau){
	reco::PFTauRef tauCandidate(taus, iTau);
	reco::PFTau pfTau= *tauCandidate;
	if((*DMF)[tauCandidate] > -1){
	  
	  if (reco::deltaR(pfTau.eta(),pfTau.phi(),jet.eta(),jet.phi()) < 0.3){
	    jetTauMatch_ = 1;
	    tauPt_  = pfTau.pt();
	    tauEta_ = pfTau.eta();
	    dmf_    = pfTau.decayMode();
	    
	    
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

	    tauPt_  = pfTau.pt();
	    tauEta_ = pfTau.eta();
	    dmf_    = pfTau.decayMode();

	    //get the matched vertex
	    int vtx_index = -1;
	    // find the 4D vertex this muon is best associated to..
	    float max_weight = 0.f;
	  for( unsigned i = 0; i < vtxs->size(); ++i ) {
	    const auto& vtx = (*vtxs)[i];      
	    const float weight = vtx.trackWeight(pfTau.leadPFChargedHadrCand()->trackRef());// check me -> Get track ref for charged hadron candidate
	    if( weight > max_weight ) {
	      max_weight = weight;
	      vtx_index = i;
	    }
	  }
	  //now do vtx variable filling
	  vtxIndex_ = vtx_index;

	  const PFCandidatePtr& packedLeadTauCand = pfTau.leadPFChargedHadrCand();
	  //const std::vector<PFRecoTauChargedHadron>& pfCHIsoTauCands = pfTau.isolationPFChargedHadrCands();
	  std::vector<PFCandidatePtr> packedCHIsoTauCands;
	  //std::cout<<"pfTau.isolationTauChargedHadronCandidates size: "<<pfTau.isolationPFChargedHadrCands().size()<<std::endl;
	  for(auto isocand :  pfTau.isolationPFChargedHadrCands()){
	    packedCHIsoTauCands.push_back(isocand);
	    tauCHIso_ += isocand->pt();
	  }

	  //currently this does not work with reco... not sure why
	  //const float dz = std::abs(packedLeadTauCand->dz(vtxs->at(0).position()));
	  //const float dxy = std::abs(packedLeadTauCand->dxy(vtxs->at(0).position()));
	  //const float dt = packedLeadTauCand->dtime(0);


	  const double res_det = 0.030; // ns

	  //std::cout<<"packedleadtaucand pt: "<<packedLeadTauCand->pt()<<" ncands: "<<packedCHIsoTauCands.size()<<std::endl;
	  newChIso0  = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 0*res_det, newChIso0_ntracks, vtxs->at(0));
	  newChIso1  = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 1*res_det, newChIso1_ntracks, vtxs->at(0));
	  newChIso2  = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 2*res_det, newChIso2_ntracks, vtxs->at(0));
	  newChIso3  = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 3*res_det, newChIso3_ntracks, vtxs->at(0)); // this is the nominal cut, 3 sigma
	  newChIso4  = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 4*res_det, newChIso4_ntracks, vtxs->at(0));
	  newChIso5  = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 5*res_det, newChIso5_ntracks, vtxs->at(0));
	  newChIso6  = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 6*res_det, newChIso6_ntracks, vtxs->at(0));
	  newChIso7  = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 7*res_det, newChIso7_ntracks, vtxs->at(0));
	  newChIso8  = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 8*res_det, newChIso8_ntracks, vtxs->at(0));
	  newChIso9  = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 9*res_det, newChIso9_ntracks, vtxs->at(0));
	  newChIso10 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 10*res_det, newChIso10_ntracks, vtxs->at(0));
	  newChIso11 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 11*res_det, newChIso11_ntracks, vtxs->at(0));
	  newChIso12 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 12*res_det, newChIso12_ntracks, vtxs->at(0));
	  newChIso13 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 13*res_det, newChIso13_ntracks, vtxs->at(0));
	  newChIso14 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 14*res_det, newChIso14_ntracks, vtxs->at(0));
	  newChIso15 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 15*res_det, newChIso15_ntracks, vtxs->at(0));
	  newChIso16 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 16*res_det, newChIso16_ntracks, vtxs->at(0));
	  newChIso17 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 17*res_det, newChIso17_ntracks, vtxs->at(0));
	  newChIso18 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 18*res_det, newChIso18_ntracks, vtxs->at(0));
	  newChIso19 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 19*res_det, newChIso19_ntracks, vtxs->at(0));
	  newChIso20 = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 20*res_det, newChIso20_ntracks, vtxs->at(0));
	  newChIsoOther = getTimedCHIsoSum(packedLeadTauCand, packedCHIsoTauCands , 100*res_det, newChIsoOther_ntracks, vtxs->at(0));

	  tauCHIso1_ = newChIso1;
	  tauCHIso2_ = newChIso2 ; 
	  tauCHIso3_ = newChIso3 ; 
	  tauCHIso4_ = newChIso4 ; 
	  tauCHIso5_ = newChIso5 ; 
	  tauCHIso6_ = newChIso6 ; 
	  tauCHIso7_ = newChIso7 ; 
	  tauCHIso8_ = newChIso8 ; 
	  tauCHIso9_ = newChIso9 ; 
	  tauCHIso10_ = newChIso10; 
	  tauCHIso11_ = newChIso11; 
	  tauCHIso12_ = newChIso12; 
	  tauCHIso13_ = newChIso13; 
	  tauCHIso14_ = newChIso14; 
	  tauCHIso15_ = newChIso15; 
	  tauCHIso16_ = newChIso16; 
	  tauCHIso17_ = newChIso17; 
	  tauCHIso18_ = newChIso18; 
	  tauCHIso19_ = newChIso19; 
	  tauCHIso20_ = newChIso20; 

	  tauCHIso1_ntracks_ = newChIso1_ntracks;
	  tauCHIso2_ntracks_ = newChIso2_ntracks; 
	  tauCHIso3_ntracks_ = newChIso3_ntracks; 
	  tauCHIso4_ntracks_ = newChIso4_ntracks; 
	  tauCHIso5_ntracks_ = newChIso5_ntracks; 
	  tauCHIso6_ntracks_ = newChIso6_ntracks; 
	  tauCHIso7_ntracks_ = newChIso7_ntracks; 
	  tauCHIso8_ntracks_ = newChIso8_ntracks; 
	  tauCHIso9_ntracks_ = newChIso9_ntracks; 
	  tauCHIso10_ntracks_ = newChIso10_ntracks; 
	  tauCHIso11_ntracks_ = newChIso11_ntracks; 
	  tauCHIso12_ntracks_ = newChIso12_ntracks; 
	  tauCHIso13_ntracks_ = newChIso13_ntracks; 
	  tauCHIso14_ntracks_ = newChIso14_ntracks; 
	  tauCHIso15_ntracks_ = newChIso15_ntracks; 
	  tauCHIso16_ntracks_ = newChIso16_ntracks; 
	  tauCHIso17_ntracks_ = newChIso17_ntracks; 
	  tauCHIso18_ntracks_ = newChIso18_ntracks; 
	  tauCHIso19_ntracks_ = newChIso19_ntracks; 
	  tauCHIso20_ntracks_ = newChIso20_ntracks; 

	    break;
	  }
	}
      }

      jetTree->Fill(); 
   }


   
}


// ------------ method called once each job just before starting event loop  ------------
void 
phase2TausRECO::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
phase2TausRECO::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
phase2TausRECO::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
}


//Gets visible 4-momentum of a particle from list of daughters
reco::Candidate::LorentzVector phase2TausRECO::GetVisibleP4(std::vector<const reco::GenParticle*>& daughters){
  reco::Candidate::LorentzVector p4_vis(0,0,0,0);
  for(size_t i = 0; i < daughters.size(); ++i){
    if (!isNeutrino(daughters[i]) && daughters[i]->status() == 1){  //list should include intermediate daughters, so check for status = 1
      p4_vis += daughters[i]->p4();
    }
  }
  return p4_vis;
}


//Creates a vector of all (including intermediate) daughters for a given mother particle
void phase2TausRECO::findDaughters(const reco::GenParticle* mother, std::vector<const reco::GenParticle*>& daughters)
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

bool phase2TausRECO::isNeutrino(const reco::Candidate* daughter)
{
  return (TMath::Abs(daughter->pdgId()) == 12 || TMath::Abs(daughter->pdgId()) == 14 || TMath::Abs(daughter->pdgId()) == 16 || TMath::Abs(daughter->pdgId()) == 18);
}


double phase2TausRECO::getTimedCHIsoSum( const PFCandidatePtr&  leadTauCand, vector<PFCandidatePtr>& isolationCands, float interval, float & n_tracks, reco::Vertex pv)
{
  float sum = 0;
  n_tracks = 0;
  for( auto cand : isolationCands ) {
    if(std::abs(cand->time()-leadTauCand->time()) < interval){
      sum += cand->pt();
      n_tracks +=1;
    }

  }

  return sum;
}


//define this as a plug-in
DEFINE_FWK_MODULE(phase2TausRECO);
