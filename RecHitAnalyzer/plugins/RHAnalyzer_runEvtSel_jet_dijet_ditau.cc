#include "MLAnalyzerRun3/RecHitAnalyzer/interface/RecHitAnalyzer.h"

using std::vector;

const unsigned nJets = 50; //TODO: use cfg level nJets_
TH1F *hNpassed_genJetMatch; 
TH1F *hNpassed_minTwoJets; 


std::vector<int> jetIDs_;                          // jet index or unique ID
std::vector<std::vector<int>> matchedGenIDs_;      // vector of matched gen IDs per jet

// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesEvtSel_jet_dijet_ditau ( TTree* tree, edm::Service<TFileService> &fs ) {

  hNpassed_genJetMatch = fs->make<TH1F>("hNpassed_genJetMatch","Jet matched to gen particle (0: No, 1: Yes)", 2, 0, 2);
  hNpassed_minTwoJets  = fs->make<TH1F>("hNpassed_minTwoJets","Atleast two jets in the event (0: No, 1: Yes)", 2, 0, 2);

} // branchesEvtSel_jet_dijet_tau()

// Run jet selection _____________________________________________________//
bool RecHitAnalyzer::runEvtSel_jet_dijet_ditau( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken( genParticleCollectionT_, genParticles );

  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);

  edm::Handle<pat::TauCollection> taus;
  iEvent.getByToken(tauCollectionT_, taus);

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertexCollectionT_, vertices);


  vJetIdxs.clear();
  jetIDs_.clear();
  matchedGenIDs_.clear();

  unsigned int nMatchedJets = 0;
  unsigned int goodVertices = 0;

  if ( debug ) std::cout << " >>>>>>>>>>>>>>>>>>>> evt:" << std::endl;
 
  if (vertices.isValid())
    if (vertices->size() > 0)
      for (auto v : *vertices)
        if (v.ndof() >= 4 && !v.isFake())
          ++goodVertices;
  if ( debug ) std::cout << " good vertices in the event (PU) = " << goodVertices << std::endl;

  if ( debug ) std::cout << " JETS IN THE EVENT = " << jets->size() << " | Selection requires minpT = " << minJetPt_ << " and maxEta = "<< maxJetEta_ << std::endl;

  // Loop over jets
  for ( unsigned iJ(0); iJ != jets->size(); ++iJ ) {

    unsigned int PdgId        = 0;
    float jetdR               = -99.;
    float taupT               = -99.;
    int tauDaughters          = -1;
    int taupi0                = -1;
    //bool JetIsTau             = false;
    float DeepTau             = -1;
    float LooseDeepTau        = -1;
    float MediumDeepTau       = -1;
    float TightDeepTau        = -1;

    pat::Jet iJet = (*jets)[iJ];
    if ( std::abs(iJet.pt())  < minJetPt_ ) continue;
    if ( std::abs(iJet.eta()) > maxJetEta_ ) continue;
    if (debug ) std::cout << "  >>>>>> Jet [" << iJ << "] ->  Pt: " << iJet.pt() << ", Eta: " << iJet.eta() << ", Phi: " << iJet.phi() << std::endl;
    if (isMC_) {
      //bool passedGenSel = false;
      unsigned int iGenParticle = 0;
      
      std::vector<int> matchedGenIDs;
      
      for ( unsigned int iG = 0; iG < genParticles->size(); iG++ ) {
	
	reco::GenParticleRef iGen( genParticles, iG );
	
	float dR = reco::deltaR( iJet.eta(),iJet.phi(), iGen->eta(),iGen->phi() );
        if ( dR > 0.4 ) continue;
	
        if ( iGen->pt() > 20 && (std::abs(iGen->pdgId()) == 11 || std::abs(iGen->pdgId()) == 13) ) break; //only clean jets (lepton veto) 
        if ( std::abs(iGen->pdgId()) == 12 || std::abs(iGen->pdgId()) == 14 || std::abs(iGen->pdgId()) == 16 ) continue;
	
        if (  isSignal_ && !( std::abs(iGen->pdgId()) == 15 && iGen->status() == 2 ) ) continue;  //only for tau signal
        if ( !isSignal_ && !isW_ && !( iGen->status() == 23 ) ) continue;                         //for QCD background
        if ( !isSignal_ &&  isW_ && !( iGen->status() == 71 ) ) continue;                         //only for W + jet background
	
        if ( debug ) std::cout << "   GEN particle " << iGenParticle << " index [" << iG << "] -> status: " << iGen->status() << ", id: " << iGen->pdgId() << ", nDaught: " << iGen->numberOfDaughters() << " nMoms: " <<iGen->numberOfMothers() << " | pt: "<< iGen->pt() << " eta: " <<iGen->eta() << " phi: " <<iGen->phi() << " | dR: " << dR << std::endl;
	
        if (debug ) std::cout << "  >>>>>> Jet [" << iJ << "] ->  Pt: " << iJet.pt() << ", Eta: " << iJet.eta() << ", Phi: " << iJet.phi() << std::endl;
	
        bool isHadronic = true;
        if ( std::abs(iGen->pdgId()) == 15 ) {
          for (unsigned int iDaughter = 0; iDaughter != iGen->numberOfDaughters(); ++iDaughter ){
	    if ( abs(iGen->daughter(iDaughter)->pdgId()) == 11 || abs(iGen->daughter(iDaughter)->pdgId()) == 13 ) isHadronic = false;
          }
          if (isSignal_ && !isHadronic) continue;
          //JetIsTau = true;
	  
	  matchedGenIDs.push_back(iG);
          tauDaughters = 0;
          taupi0 = 0;
          for (unsigned int iDaughter = 0; iDaughter != iGen->numberOfDaughters(); ++iDaughter ){
            if ( debug ) std::cout << "    Tau daughter [" << iDaughter << "] : "<<  std::abs(iGen->daughter(iDaughter)->pdgId()) << std::endl;
            if ( abs(iGen->daughter(iDaughter)->pdgId()) == 111 ) taupi0++;
            if ( iGen->daughter(iDaughter)->charge() == 0 ) continue;
            tauDaughters++;
          }
          if ( debug ) std::cout << "    Tau prongs = " << tauDaughters << " + Tau pi0 = " << taupi0 << std::endl;
	  
          if (!isSignal_){
            //passedGenSel = false;  //only for background
            break;                 //only for background
          }	  
        } 
	//passedGenSel = true;
        ++iGenParticle;
	
      } // primary gen particles
      
      if (!matchedGenIDs.empty()) {
	jetIDs_.push_back(iJ);  // Index of the jet in the collection
	matchedGenIDs_.push_back(matchedGenIDs);  // All matched gen IDs for this jet
      }
    } // is MC selection
    
  } // reco jets
  
  // After you've filled jetIDs_ and matchedGenIDs_ for the current event
  
  if(matchedGenIDs_.empty()){
    hNpassed_genJetMatch->Fill(0);
    return false;
  }
  hNpassed_genJetMatch->Fill(1);
  
  //apply lepton veto on jets here
  
  if(jetIDs_.size()<2){
    hNpassed_minTwoJets->Fill(0);
    return false;
  }
  hNpassed_minTwoJets->Fill(1);

   for (size_t i = 0; i < jetIDs_.size(); ++i) {
    std::cout << "********************Jet ID: " << jetIDs_[i] << " matched to GenParticles IDs: ";
    vJetIdxs.push_back(jetIDs_[i]);
    for (size_t j = 0; j < matchedGenIDs_[i].size(); ++j) {
      std::cout << matchedGenIDs_[i][j];
      if (j != matchedGenIDs_[i].size() - 1) std::cout << ", ";
    }
    std::cout << std::endl;
  }

 if(vJetIdxs.empty()){
	 return false;
 }
  
  // Check jet multiplicity
  if( debug) std::cout << " Matched jets " << jetIDs_.size() << std::endl;
  if ( debug ) std::cout << " >> Event contains a tau candidate" << std::endl;
  return true;
  
} // runEvtSel_jet_dijet_tau()

// Fill branches and histograms _____________________________________________________//
void RecHitAnalyzer::fillEvtSel_jet_dijet_ditau ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);

  //h_tau_jet_nJet->Fill( vJetIdxs.size() );

  for ( size_t i=0; i < vJetIdxs.size(); ++i ) {
int jetIndex = jetIDs_[i];
std::cout << "HBHE passed jetIndex:" << jetIndex << std::endl;

    // Fill histograms 
  
  }

} // fillEvtSel_jet_dijet_tau()
