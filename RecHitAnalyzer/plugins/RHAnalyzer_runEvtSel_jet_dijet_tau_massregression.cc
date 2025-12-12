#include "MLAnalyzerRun3/RecHitAnalyzer/interface/RecHitAnalyzer.h"

using std::vector;
const unsigned nJets = 50; //TODO: use cfg level nJets_

//gen quatities
int v_mr_NGen_a_; //number of Gen pseudoscalars
vector<float> v_mr_Gen_mass_a_;
vector<float> v_mr_Gen_pt_a_;
int v_mr_NGenTaus_; //number of Gen taus
vector<float> v_mr_Gen_tau_pt_;
vector<float> v_mr_Gen_tau_eta_;
vector<float> v_mr_Gen_tau_phi_;
vector<float> v_mr_Gen_tau1_tau2_dR_;

//jet quantities
int v_mr_NJets_; //number of reco'ed jets
vector<float> v_mr_jet_mass_;
vector<float> v_mr_jet_pt_;
vector<float> v_mr_jet_eta_;
vector<float> v_mr_jet_phi_;
vector<float> v_mr_jet_genTau_dR_;
vector<float> v_mr_jet_tau_dR_;
vector<float> v_mr_jet1_jet2_dR_;
int v_mr_NGenTau_JetMatched_; //number of gen taus matched to one jet

//tau quantities
int v_mr_NTaus_;
vector<float> v_mr_tau_mass_;
vector<float> v_mr_tau_pt_;
vector<float> v_mr_tau_eta_;
vector<float> v_mr_tau_phi_;
vector<float> v_mr_tau_genTau_dR_;
vector<float> v_mr_tau1_tau2_dR_;
int v_mr_NTau_JetMatched_; //number of taus matched to one jet

std::map<int, std::vector<int>> v_mr_jetToGenMap_;  // jet index -> matched gen particle indices
std::vector<int> v_mr_jetIDs_;                          // jet index or unique ID
std::vector<std::vector<int>> v_mr_matchedGenIDs_;      // vector of matched gen IDs per jet
std::vector<std::vector<int>> v_mr_matchedTauIDs_;      // vector of matched tau IDs per jet
std::map<int, std::vector<int>> v_mr_jetToTauMap_; 

//////////////////////////////////////////////////////////////////////////////
// Initialize branches _____________________________________________________//
/////////////////////////////////////////////////////////////////////////////
void RecHitAnalyzer::branchesEvtSel_jet_dijet_tau_massregression ( TTree* tree, edm::Service<TFileService> &fs ) {
  tree->Branch("a_NumGen",          &v_mr_NGen_a_);
  tree->Branch("a_mass",            &v_mr_Gen_mass_a_);
  tree->Branch("a_pt",              &v_mr_Gen_pt_a_);
  tree->Branch("NgenTaus",          &v_mr_NGenTaus_);
  tree->Branch("genTau_pt",         &v_mr_Gen_tau_pt_);
  tree->Branch("genTau_eta",        &v_mr_Gen_tau_eta_);
  tree->Branch("genTau_phi",        &v_mr_Gen_tau_phi_);
  tree->Branch("genTau1_Tau2_dr",   &v_mr_Gen_tau1_tau2_dR_);

  tree->Branch("Njets",             &v_mr_NJets_);
  tree->Branch("jet_mass",          &v_mr_jet_mass_);
  tree->Branch("jet_pt",            &v_mr_jet_pt_);
  tree->Branch("jet_eta",           &v_mr_jet_eta_);
  tree->Branch("jet_phi",           &v_mr_jet_phi_);
  tree->Branch("jet_gen_dR",        &v_mr_jet_genTau_dR_);
  tree->Branch("jet_tau_dR",        &v_mr_jet_tau_dR_);
  tree->Branch("jet1_jet2_dR",      &v_mr_jet1_jet2_dR_);
  tree->Branch("NGenTau_JetMatched",&v_mr_NGenTau_JetMatched_);

  tree->Branch("NTaus",             &v_mr_NTaus_);
  tree->Branch("tau_mass",          &v_mr_tau_mass_);
  tree->Branch("tau_pt",            &v_mr_tau_pt_);
  tree->Branch("tau_eta",           &v_mr_tau_eta_);
  tree->Branch("tau_phi",           &v_mr_tau_phi_);
  tree->Branch("tau_gen_dR",        &v_mr_tau_genTau_dR_);
  tree->Branch("tau1_tau2_dR",      &v_mr_tau1_tau2_dR_);
  tree->Branch("NTau_JetMatched",   &v_mr_NTau_JetMatched_);
} // branchesEvtSel_jet_dijet_tau_massregression()

////////////////////////////////////////////////////////////////////////////
// Run jet selection _____________________________________________________//
////////////////////////////////////////////////////////////////////////////
bool RecHitAnalyzer::runEvtSel_jet_dijet_tau_massregression( const edm::Event& iEvent, const edm::EventSetup& iSetup ){

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken( genParticleCollectionT_, genParticles );

  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);

  edm::Handle<pat::TauCollection> taus;
  iEvent.getByToken(tauCollectionT_, taus);

  vJetIdxs.clear();
  v_mr_jetIDs_.clear();
  v_mr_matchedGenIDs_.clear();
  v_mr_matchedTauIDs_.clear();
  v_mr_jetToGenMap_.clear(); 
  v_mr_jetToTauMap_.clear();

  unsigned int nMatchedJets = 0;
  if ( debug ) std::cout << " JETS IN THE EVENT = " << jets->size() << " | Selection requires minpT = " << minJetPt_ << " and maxEta = "<< maxJetEta_ << std::endl;
  
  // Loop over jets
  for ( unsigned iJ(0); iJ != jets->size(); ++iJ ) {
    pat::Jet iJet = (*jets)[iJ];
    if (iJet.pt() < minJetPt_ ) continue;
    if (std::abs(iJet.eta()) > maxJetEta_ ) continue;
    if (debug ) std::cout << "  >>>>>> Jet [" << iJ << "] ->  Pt: " << iJet.pt() << ", Eta: " << iJet.eta() << ", Phi: " << iJet.phi() << " ,mass: "<< iJet.mass()<<  std::endl;
    
    unsigned int iGenParticle = 0;
    std::vector<int> matchedGenIDs;
    for ( unsigned int iG = 0; iG < genParticles->size(); iG++ ) {
      reco::GenParticleRef iGen( genParticles, iG );
      float dR = reco::deltaR( iJet.eta(),iJet.phi(), iGen->eta(),iGen->phi() );
      if ( dR > 0.4 ) continue;
      if ( !( std::abs(iGen->pdgId()) == 15 && iGen->status() == 2 ) ) continue;  //only for tau signal
      if(iGen->numberOfMothers() != 1) continue;
            
      if ( debug ) std::cout << "   GEN particle " << iGenParticle << " index [" << iG << "] -> status: " << iGen->status() << ", id: " << iGen->pdgId() << ", nDaught: " << iGen->numberOfDaughters() << " nMoms: " <<iGen->numberOfMothers() << " | pt: "<< iGen->pt() << " eta: " <<iGen->eta() << " phi: " <<iGen->phi() << " | dR: " << dR << std::endl;
      if (debug ) std::cout << "  >>>>>> Jet [" << iJ << "] ->  Pt: " << iJet.pt() << ", Eta: " << iJet.eta() << ", Phi: " << iJet.phi() << std::endl;
      
      bool isHadronic = true;
      if ( std::abs(iGen->pdgId()) == 15 ) {
	for (unsigned int iDaughter = 0; iDaughter != iGen->numberOfDaughters(); ++iDaughter ){
	  if ( abs(iGen->daughter(iDaughter)->pdgId()) == 11 || abs(iGen->daughter(iDaughter)->pdgId()) == 13 ) isHadronic = false;
	}
	if (isSignal_ && !isHadronic) continue;
	matchedGenIDs.push_back(iG);
      }
      
      ++iGenParticle;
    } // primary gen particles
    
    if (!matchedGenIDs.empty()) {
      std::cout << " Matched gen ID not empty" << std::endl;
      v_mr_jetIDs_.push_back(iJ);  // Index of the jet in the collection
      v_mr_matchedGenIDs_.push_back(matchedGenIDs);  // All matched gen IDs for this jet
      v_mr_jetToGenMap_[iJ] = matchedGenIDs;
    }
  } // reco jets


  // After you've filled jetIDs_ and matchedGenIDs_ for the current event
  for (int jetIdx : v_mr_jetIDs_) {
    const pat::Jet& jet = (*jets)[jetIdx];
    std::vector<int> matchedTauIDs;
    for (unsigned int iT = 0; iT < taus->size(); ++iT) {
      const pat::Tau& tau = (*taus)[iT];
      if (tau.pt() < 10.0 || std::abs(tau.eta()) > 3.0) continue;
      if (reco::deltaR(jet, tau) > 0.4) continue;
      matchedTauIDs.push_back(iT);
      std::cout << "Jet "<< jetIdx << " with pt:" << jet.pt() << " matched to Tau pt:" << tau.pt() <<  " jet tau dR:" << reco::deltaR(jet, tau) << std::endl;
    }
    if (!matchedTauIDs.empty()) {
      v_mr_matchedTauIDs_.push_back(matchedTauIDs); 
      v_mr_jetToTauMap_[jetIdx] = matchedTauIDs;
    }
  }
  
  
  for (size_t i = 0; i < v_mr_jetIDs_.size(); ++i) {
    if( debug ) std::cout << "********************Jet ID: " << v_mr_jetIDs_[i] << " matched to GenParticles IDs: ";
    vJetIdxs.push_back(v_mr_jetIDs_[i]);
    for (size_t j = 0; j < v_mr_matchedGenIDs_[i].size(); ++j) {
      if( debug ) std::cout << v_mr_matchedGenIDs_[i][j];
      if (j != v_mr_matchedGenIDs_[i].size() - 1) { if( debug ) std::cout << ", ";}
    }
    if(debug)std::cout << std::endl;
    for (size_t jj = 0; jj < v_mr_matchedTauIDs_[i].size(); ++jj) {
      if( debug ) std::cout << "Matched Tau ID:" << v_mr_matchedTauIDs_[i][jj];
      if (jj != v_mr_matchedTauIDs_[i].size() - 1) { if( debug ) std::cout << " , ";}
    }
    if(debug)std::cout << std::endl;
  }
  
  if(vJetIdxs.empty()){
    return false;
  }
  
  // Check jet multiplicity
  if ( debug ) std::cout << " Matched jets " << v_mr_jetIDs_.size() << std::endl;
  if ( debug ) std::cout << " >> Event contains a tau candidate" << std::endl;
  return true;
  
} // runEvtSel_jet_dijet_tau()

// Fill branches and histograms _____________________________________________________//
void RecHitAnalyzer::fillEvtSel_jet_dijet_tau_massregression ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {
  
  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);
  
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken( genParticleCollectionT_, genParticles );
 
  edm::Handle<pat::TauCollection> taus;
  iEvent.getByToken(tauCollectionT_, taus);
  
  v_mr_Gen_mass_a_.clear();
  v_mr_Gen_pt_a_.clear();
  v_mr_Gen_tau_pt_.clear();
  v_mr_Gen_tau_eta_.clear();
  v_mr_Gen_tau_phi_.clear();
  v_mr_Gen_tau1_tau2_dR_.clear();
 
  v_mr_jet_mass_.clear();
  v_mr_jet_pt_.clear();
  v_mr_jet_eta_.clear();
  v_mr_jet_phi_.clear();
  v_mr_jet_genTau_dR_.clear();
  v_mr_jet_tau_dR_.clear();
  v_mr_jet1_jet2_dR_.clear();
 
  v_mr_tau_mass_.clear();
  v_mr_tau_pt_.clear();
  v_mr_tau_eta_.clear();
  v_mr_tau_phi_.clear();
  v_mr_tau_genTau_dR_.clear();
  v_mr_tau1_tau2_dR_.clear();

  // Track unique mothers (a and a-bar) by pointer
  std::set<const reco::Candidate*> uniqueMothers;
  
  // Loop over selected jets (in original order from vJetIdxs)
  for (int jetIdx : vJetIdxs) {
    const pat::Jet& jet = (*jets)[jetIdx];
    
    v_mr_jet_mass_.push_back(jet.mass());
    v_mr_jet_pt_.push_back(jet.pt());
    v_mr_jet_eta_.push_back(jet.eta());
    v_mr_jet_phi_.push_back(jet.phi());
    
    std::cout<< "--------------------------------------------------------------" << std::endl;
    std::cout<< "FIlling out jet info pt and eta:" << jet.pt() << "  eta:" << jet.eta() << std::endl;
    float minGenDR = 999.0f;
    const reco::GenParticle* bestGenTau = nullptr;
    std::vector<const reco::GenParticle*> genTausThisJet;

    // =============================================
    // Gen tau loop: collect truth + deduplicate mothers
    // =============================================
    auto genIt = v_mr_jetToGenMap_.find(jetIdx);
    if (genIt != v_mr_jetToGenMap_.end()) {
      for (int gIdx : genIt->second) {
        const reco::GenParticle& genTau = (*genParticles)[gIdx];
        genTausThisJet.push_back(&genTau);

        // Store gen tau kinematics
        v_mr_Gen_tau_pt_.push_back(genTau.pt());
        v_mr_Gen_tau_eta_.push_back(genTau.eta());
        v_mr_Gen_tau_phi_.push_back(genTau.phi());

	std::cout <<" Gen tau pt:" << genTau.pt() << " eta:" << genTau.eta() << std::endl; 

	// === DEDUPLICATED MOTHER INFO ===
        if (genTau.numberOfMothers() > 0) {
          const reco::Candidate* mother = genTau.mother(0);
          if (uniqueMothers.insert(mother).second) {  // true only the first time
            v_mr_Gen_mass_a_.push_back(mother->mass());
            v_mr_Gen_pt_a_.push_back(mother->pt());
	    std::cout << "Pseudoscalar mother mass :" << mother->mass() << " pt:" << mother->pt() << std::endl;
	    
	    // Collect all status=2 tau daughters
	    std::vector<const reco::Candidate*> status2Taus;
	    unsigned int nDaughters = mother->numberOfDaughters();
	    std::cout << "Mother has " << nDaughters << " daughters:" << std::endl;
	    
	    for (unsigned int d = 0; d < nDaughters; ++d) {
	      const reco::Candidate* dau = mother->daughter(d);
	      int pdgId = dau->pdgId();
	      
	      if (std::abs(pdgId) == 15 && dau->status() == 2) {
		status2Taus.push_back(dau);
		std::cout << "  Daughter " << d << ": tau (pdgId=" << pdgId
			  << "), status=" << dau->status()
			  << ", pt=" << dau->pt()
			  << ", eta=" << dau->eta()
			  << ", phi=" << dau->phi() << std::endl;
	      } //status and pdgId
	    }//ndaughters
	    
	    std::cout << "Found " << status2Taus.size() << " status=2 taus among daughters." << std::endl;
	    
	    // Compute and store pairwise dR between status=2 taus (you expect exactly 2)
	    if (status2Taus.size() >= 2) {
	      // Assuming exactly 2 (as you said there will definitely be 2)
	      float tauPairDR = reco::deltaR(*status2Taus[0], *status2Taus[1]);
	      v_mr_Gen_tau1_tau2_dR_.push_back(tauPairDR);
	      
	      std::cout << "dR between the two status=2 taus: " << tauPairDR << std::endl; 
	    }
	  }
	}
	  
        // Update closest gen tau for this jet
        float dR = reco::deltaR(jet, genTau);
        v_mr_jet_genTau_dR_.push_back(dR);
      } //gen loop
    }//jetGen map
    
    //float tau1_tau2_dR = genTausThisJet.size() >= 2 ? reco::deltaR(*genTausThisJet[0], *genTausThisJet[1]) : -1.0f;
    //v_mr_Gen_tau1_tau2_dR_.push_back(tau1_tau2_dR);
    //std::cout <<"jet and gen dR:" << minGenDR << " tau1-tau2 dR gen level: " << tau1_tau2_dR << std::endl;

    // =============================================
    // Reco tau lop: store ALL matched reco taus
    // =============================================
    float minRecoDR = 999.0f;
    std::vector<const pat::Tau*> recoTausThisJet;

    auto tauIt = v_mr_jetToTauMap_.find(jetIdx);
    if (tauIt != v_mr_jetToTauMap_.end()) {
      for (int tIdx : tauIt->second) {
        const pat::Tau& tau = (*taus)[tIdx];
        recoTausThisJet.push_back(&tau);

        v_mr_tau_mass_.push_back(tau.mass());
        v_mr_tau_pt_.push_back(tau.pt());
        v_mr_tau_eta_.push_back(tau.eta());
        v_mr_tau_phi_.push_back(tau.phi());

        float dRjet = reco::deltaR(jet, tau);
	v_mr_jet_tau_dR_.push_back(dRjet);
	std::cout << " tau pt:" << tau.pt() << " eta:" << tau.eta() << " mass:" << tau.mass() << " tau-jet dR:" << dRjet << std::endl;
	//v_mr_tau_genTau_dR_.push_back(bestGenTau ? reco::deltaR(tau, *bestGenTau) : -1.0f);
      }
    }//jetTauMap
  }  

  // Safety: if no mother found (should never happen in signal), push dummies
  if (v_mr_Gen_mass_a_.empty()) {
    v_mr_Gen_mass_a_.push_back(-1.0f);
    v_mr_Gen_pt_a_.push_back(-1.0f);
  }
}//fill  

