#include "MLAnalyzerRun3/RecHitAnalyzer/interface/RecHitAnalyzer.h"

using std::vector;

const unsigned nJets = 50; //TODO: use cfg level nJets_
TH1F *hNpassed_genJetMatch; 
TH1F *hNpassed_minTwoJets; 
TH1F *hNpassed_eleVeto; 
TH1F *hNpassed_muVeto; 
TH1F *hNpassed_hbheCrop; 

//gen variables
vector<float> v_att_genHiggs_M_;
vector<float> v_att_genPS_M_;
vector<float> v_att_genTau_pT_;
vector<float> v_att_genTau_eta_;
vector<float> v_att_genTau_phi_;
vector<float> v_att_genTau_prongs_;
vector<float> v_att_genTau1Tau2_dR_;

//jet variables
int v_att_tau_njet_;
vector<float> v_att_tau_jet_m0_;
vector<float> v_att_tau_jet_pt_;
vector<float> v_att_tau_jet_eta_;
vector<float> v_att_tau_jet_phi_;
vector<float> v_att_tau_jetIsSignal_;
vector<float> v_att_tau_jetdR_;

std::map<int, std::vector<int>> jetToGenMap_;  // jet index -> matched gen particle indices
std::vector<int> jetIDs_;                          // jet index or unique ID
std::vector<std::vector<int>> matchedGenIDs_;      // vector of matched gen IDs per jet
std::vector<std::vector<int>> matchedEleIDs_;      // vector of matched ele IDs per jet
std::vector<std::vector<int>> matchedMuIDs_;      // vector of matched muon IDs per jet

std::map<int, std::vector<int>> jetToEleMap_; 
std::map<int, std::vector<int>> jetToMuMap_;

int jetMatchedEle;
int jetMatchedMu;

// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesEvtSel_jet_dijet_ditau ( TTree* tree, edm::Service<TFileService> &fs ) {

  hNpassed_genJetMatch = fs->make<TH1F>("hNpassed_genJetMatch","Jet matched to gen particle (0: No, 1: Yes)", 2, 0, 2);
  hNpassed_eleVeto     = fs->make<TH1F>("hNpassed_eleVeto","Jets is not matched to electron in the event (0: No, 1: Yes)", 2, 0, 2);
  hNpassed_muVeto      = fs->make<TH1F>("hNpassed_muVeto","Jet is not matched to muon in the event (0: No, 1: Yes)", 2, 0, 2);
  hNpassed_hbheCrop    = fs->make<TH1F>("hNpassed_hbheCrop","Jets in the events passing HB-HE edge cut (0: No, 1: Yes)", 2, 0, 2);
  hNpassed_minTwoJets  = fs->make<TH1F>("hNpassed_minTwoJets","Atleast two jets in the event (0: No, 1: Yes)", 2, 0, 2);

  //gen variables
  tree->Branch("genTauPt",            &v_att_genTau_pT_);
  tree->Branch("genTauEta",           &v_att_genTau_eta_);
  tree->Branch("genTauPhi",           &v_att_genTau_phi_);
  tree->Branch("genTau1Tau2dR",       &v_att_genTau1Tau2_dR_);

  //jet variables
  tree->Branch("nJets",            &v_att_tau_njet_);
  tree->Branch("jetM",             &v_att_tau_jet_m0_);
  tree->Branch("jetPt",            &v_att_tau_jet_pt_);
  tree->Branch("jetEta",           &v_att_tau_jet_eta_);
  tree->Branch("jetPhi",           &v_att_tau_jet_phi_);
  tree->Branch("jetIsSignal",      &v_att_tau_jetIsSignal_);

} // branchesEvtSel_jet_dijet_tau()

// Run jet selection _____________________________________________________//
bool RecHitAnalyzer::runEvtSel_jet_dijet_ditau( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken( genParticleCollectionT_, genParticles );

  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);

  edm::Handle<pat::ElectronCollection> eles;
  iEvent.getByToken(eleCollectionT_, eles);

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertexCollectionT_, vertices);

  edm::Handle<pat::MuonCollection> mus;
  iEvent.getByToken(muCollectionT_, mus);

  vJetIdxs.clear();
  jetIDs_.clear();
  matchedGenIDs_.clear();
  matchedEleIDs_.clear();
  matchedMuIDs_.clear();
  jetToGenMap_.clear(); 
  jetToEleMap_.clear();
  jetToMuMap_.clear();

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
    if (debug ) std::cout << "  >>>>>> Jet [" << iJ << "] ->  Pt: " << iJet.pt() << ", Eta: " << iJet.eta() << ", Phi: " << iJet.phi() << " ,mass: "<< iJet.mass()<<  std::endl;
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
	jetToGenMap_[iJ] = matchedGenIDs;
      }
    } // is MC selection
    
  } // reco jets
  
  // After you've filled jetIDs_ and matchedGenIDs_ for the current event
  
  if(matchedGenIDs_.empty()){
    hNpassed_genJetMatch->Fill(0);
    return false;
  }
  hNpassed_genJetMatch->Fill(1);
 
  jetMatchedEle = 0; 
  //apply lepton veto on jets here
  for (size_t i = 0; i < jetIDs_.size(); ++i) {
    int jetIdx_i = jetIDs_[i];
    pat::Jet iJet = (*jets)[jetIdx_i];

    //check the jet matching to electron
    std::vector<int> matchedEleIDs;
    for (size_t j = 0; j < eles->size(); ++j){
      pat::Electron iEle = (*eles)[j];
      
      if ( std::abs(iEle.pt())  < 10.0 ) continue;
      if ( std::abs(iEle.eta()) > 3.0 ) continue;
      
      float dR = reco::deltaR( iJet.eta(),iJet.phi(), iEle.eta(),iEle.phi() );
      if ( dR > 0.4 ) continue;
      bool eleID  = iEle.electronID("cutBasedElectronID-RunIIIWinter22-V1-veto");
      if (eleID == false) continue;
      std::cout << "------------------------------------------------------------------------------electron pt:" << iEle.pt() << std::endl;
      jetMatchedEle ++;
      
      matchedEleIDs.push_back(j);
    }//electron loop

    //check the jet matching to muon
    std::vector<int> matchedMuIDs;
    for (size_t mm = 0; mm < mus->size(); ++mm){
      pat::Muon iMu = (*mus)[mm];
      
      if ( std::abs(iMu.pt())  < 4.0 ) continue;
      if ( std::abs(iMu.eta()) > 2.4 ) continue;
      
      float dR = reco::deltaR( iJet.eta(),iJet.phi(), iMu.eta(),iMu.phi() );
      if ( dR > 0.4 ) continue;
      bool muonID  = iMu.isMediumMuon();
      if (muonID == false) continue;
      std::cout << "-----------------------------------------------------------------------------------muon pt:" << iMu.pt() << std::endl;
      jetMatchedMu ++;
      
      matchedMuIDs.push_back(mm);
    }//muon loop
    
    if (!matchedEleIDs.empty()) {
      matchedEleIDs_.push_back(matchedEleIDs);  // All matched gen IDs for this jet
      jetToEleMap_[i] = matchedEleIDs;
    }

    if (!matchedMuIDs.empty()) {
      matchedMuIDs_.push_back(matchedMuIDs);  // All matched gen IDs for this jet
      jetToMuMap_[i] = matchedMuIDs;
    }
  }//jet loop
  
  if(!matchedEleIDs_.empty()){
    hNpassed_eleVeto->Fill(0);
    return false;
  }
  hNpassed_eleVeto->Fill(1);
  
   if(!matchedMuIDs_.empty()){
    hNpassed_muVeto->Fill(0);
    return false;
  }
  hNpassed_muVeto->Fill(1);


  for (size_t i = 0; i < jetIDs_.size(); ++i) {
    if( debug ) std::cout << "********************Jet ID: " << jetIDs_[i] << " matched to GenParticles IDs: ";
    vJetIdxs.push_back(jetIDs_[i]);
    for (size_t j = 0; j < matchedGenIDs_[i].size(); ++j) {
      if( debug ) std::cout << matchedGenIDs_[i][j];
      if (j != matchedGenIDs_[i].size() - 1) { if( debug ) std::cout << ", ";}
    }
    if(debug)std::cout << std::endl;
  }
  
  if(vJetIdxs.empty()){
    return false;
  }
  
  // Check jet multiplicity
  if ( debug ) std::cout << " Matched jets " << jetIDs_.size() << std::endl;
  if ( debug ) std::cout << " >> Event contains a tau candidate" << std::endl;
  return true;
  
} // runEvtSel_jet_dijet_tau()

// Fill branches and histograms _____________________________________________________//
void RecHitAnalyzer::fillEvtSel_jet_dijet_ditau ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken( genParticleCollectionT_, genParticles );

  v_att_tau_njet_ = vJetIdxs.size();
  v_att_genTau_pT_.clear();
  v_att_genTau_eta_.clear();
  v_att_genTau_phi_.clear();
  v_att_genTau1Tau2_dR_.clear();
  v_att_tau_jet_m0_.clear();
  v_att_tau_jet_pt_.clear();
  v_att_tau_jet_eta_.clear();
  v_att_tau_jet_phi_.clear();

  //h_tau_jet_nJet->Fill( vJetIdxs.size() );
    /*for (const auto& pair : jetToGenMap_) {
 
      const int jetIdx = pair.first;
      if (std::find(vJetIdxs.begin(), vJetIdxs.end(), jetIdx) == vJetIdxs.end())continue; 
      const std::vector<int>& matchedGenIdxs = pair.second;
      for (int genIdx : matchedGenIdxs) {
        if(debug) std::cout<< "**************************** jet matched to "<< jetIdx << "  genIdx:" << genIdx << std::endl;
      }
    }*/

   if(vJetIdxs.size()<2){
    hNpassed_minTwoJets->Fill(0);
    return;
  }
  hNpassed_minTwoJets->Fill(1);


  std::vector<std::vector<const reco::GenParticle*>> allTauDaughters;
  
  for (size_t i = 0; i < genParticles->size(); ++i) {
    const reco::GenParticle& gen = genParticles->at(i);
    
    if (std::abs(gen.pdgId()) != 25) continue; 
    
    std::vector<const reco::GenParticle*> tauDaughters;
    
    for (unsigned int d = 0; d < gen.numberOfDaughters(); ++d) {
      const reco::GenParticle* dau = dynamic_cast<const reco::GenParticle*>(gen.daughter(d));
      if (!dau) continue;
      if (std::abs(dau->pdgId()) == 15 && dau->status()==2 ) {
	tauDaughters.push_back(dau);
      }
      else if(std::abs(dau->pdgId()) == 15 &&  dau->status()==23){
	// Now look at daughter's daughters
	int nGrandDau = dau->numberOfDaughters();
	for (int j = 0; j < nGrandDau; ++j) {
	  const reco::GenParticle* grandDau = dynamic_cast<const reco::GenParticle*>(dau->daughter(j));
	  if (std::abs(grandDau->pdgId())==15){
	    if( debug )std::cout<<"grand dau pdgID "<<grandDau->pdgId() << " status:" << grandDau->status() << std::endl;
	    tauDaughters.push_back(grandDau);
  	  }	  
	}
      }
    } //no. of daughters
    if (!tauDaughters.empty()) {
      allTauDaughters.push_back(tauDaughters);
    }
  }//genparticles loop

  if( debug ) std::cout << "Found " << allTauDaughters.size() << " pseudoscalars with tau daughters" << std::endl;

  for (size_t psIdx = 0; psIdx < allTauDaughters.size(); ++psIdx) {
    std::cout << " - Pseudoscalar " << psIdx << " has " << allTauDaughters[psIdx].size() << " tau daughters" << std::endl;
    
    const std::vector<const reco::GenParticle*>& tauDaughters = allTauDaughters[psIdx];   
    if (tauDaughters.size() != 2)continue;
    
    const reco::GenParticle* tau1 = tauDaughters[0];
    const reco::GenParticle* tau2 = tauDaughters[1];
    
    bool tau1Matched = false, tau2Matched = false;
    
    for (const auto& pair : jetToGenMap_) {
      const int jetIdx = pair.first;
      if (std::find(vJetIdxs.begin(), vJetIdxs.end(), jetIdx) == vJetIdxs.end())continue;

      const std::vector<int>& matchedGenIdxs = pair.second;
      
      for (int genIdx : matchedGenIdxs) {
	if(debug) std::cout<< "**************************** genIdx:" << genIdx << std::endl;
	const reco::GenParticle& matchedGen = genParticles->at(genIdx);
	if (&matchedGen == tau1){ tau1Matched = true; if( debug )std::cout <<  " matched genID:" << genIdx << " tau 1 Macthed "<<  std::endl;}
	if (&matchedGen == tau2){ tau2Matched = true; if( debug )std::cout <<  " matched genID:" << genIdx << " tau 2 Macthed "<<  std::endl;}
      }
    }
    
    if (tau1Matched || tau2Matched) {
      float dR = reco::deltaR(tau1->eta(), tau1->phi(), tau2->eta(), tau2->phi());
      if(debug) std::cout << " coming in dR loop---------------------------------------------------" << std::endl;
      v_att_genTau_pT_.push_back(tau1->pt());
      v_att_genTau_pT_.push_back(tau2->pt());
      v_att_genTau_eta_.push_back(tau1->eta());
      v_att_genTau_eta_.push_back(tau2->eta());
      v_att_genTau_phi_.push_back(tau1->phi());
      v_att_genTau_phi_.push_back(tau2->phi());
      v_att_genTau1Tau2_dR_.push_back(dR);
      
      if( debug )std::cout << "[dR from gen pseudoscalar daughters] dR = " << dR << " tau1 pt :" << tau1->pt() << " eta:"<< tau1->eta() << " status:" << tau1->status();
      if ( debug ) std::cout << " tau2 pt:"<< tau2->pt() << " eta:"<< tau2->eta() << " status:" << tau2->status() << std::endl;
    }
    else {std::cout << "none of the gen tau matched to gen-jet value map gen particle" << std::endl;}
    
  } //2 PS loop

  // jet loop ///////
  ///////////////////
  for ( size_t i=0; i < vJetIdxs.size(); ++i ) {
    int jetIdx_i = vJetIdxs[i];
    if(debug)std::cout << "HBHE passed jetIdx_i:" << jetIdx_i << std::endl;
    pat::Jet iJet = (*jets)[jetIdx_i];
    if(debug)std::cout << "-------------------------->Jet pt after HBHE cut " << iJet.pt() << "  eta:" << iJet.eta() << " phi:"<< iJet.phi() << " mass:" << iJet.mass()<< std::endl;
    // Fill histograms
    v_att_tau_jet_m0_.push_back(iJet.mass());
    v_att_tau_jet_pt_.push_back(iJet.pt());
    v_att_tau_jet_eta_.push_back(iJet.eta());
    v_att_tau_jet_phi_.push_back(iJet.phi());
    

  } // jetID loop
  
} // fillEvtSel_jet_dijet_tau()
