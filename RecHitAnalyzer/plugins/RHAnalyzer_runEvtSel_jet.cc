
#include "MLAnalyzerRun3/RecHitAnalyzer/interface/RecHitAnalyzer.h"

using std::vector;

// Run jet event selection ////////////////////////////////
// Only the jet seed finding is explicitly done here.
// The explicit jet selection routines are contained in
// the individual branchesEvtSel_jet_*() and runEvtSel_jet_*()

const int search_window = 7;
const int image_padding = 12;
unsigned int jet_runId_;
unsigned int jet_lumiId_;
unsigned long long jet_eventId_;
vector<float> vJetSeed_iphi_;
vector<float> vJetSeed_ieta_;
vector<int>   vFailedJetIdx_;
int hltAccept_;
int HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1_;
int HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60_;
int HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet75_;
int HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_OneProng_M5to80_;

TH1F * hNpassed_hlt;

//const std::string task_ = ""; // TODO: put switch at cfg level
//const std::string task_ = "dijet_gg_qq"; // TODO: put switch at cfg level
//const std::string task_ = "jet_tau"; //for GluGluHToTauTau
//const std::string task_ = "dijet_ditau";
//const std::string task_ = "dijet_tau_massregression";
//const std::string task_ = "dijet_ele_massregression";
//const std::string task_ = "jet_ele_classification";
//const std::string task_ = "jet_background";


// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesEvtSel_jet ( TTree* tree, edm::Service<TFileService> &fs ) {

  tree->Branch("hltAccept",      &hltAccept_);
  tree->Branch("eventId",        &jet_eventId_);
  tree->Branch("runId",          &jet_runId_);
  tree->Branch("HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1",          &HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1_);
  tree->Branch("HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60",  &HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60_);
  tree->Branch("HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet75",  &HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet75_);
  tree->Branch("HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_OneProng_M5to80", &HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_OneProng_M5to80_);
  tree->Branch("lumiId",         &jet_lumiId_);
  tree->Branch("jetSeed_iphi",   &vJetSeed_iphi_);
  tree->Branch("jetSeed_ieta",   &vJetSeed_ieta_);

  hNpassed_hlt      = fs->make<TH1F>("hNpassed_hlt", "isPassed;isPassed;N", 2, 0., 2);


  // Fill branches in explicit jet selection
  if ( task_ == "dijet_ditau" ) {
    branchesEvtSel_jet_dijet_ditau( tree, fs );
    if(isMC_)branchesEvtSel_jet_dijet_ditau_h2aa4Tau( tree, fs );
  } 
  
} // branchesEvtSel_jet()

// Run event selection ___________________________________________________________________//
bool RecHitAnalyzer::runEvtSel_jet ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {


  edm::Handle<edm::TriggerResults> hltresults;
  iEvent.getByToken(triggerResultsToken_,hltresults);

 //if(iEvent.id().event()<21930)return false;
  if (!hltresults.isValid()) {
   std::cout << "!!! Error in getting TriggerResults product from Event !!!" << std::endl;
   //triger_valid = false;
  }

  int ntrigs = hltresults->size();
  edm::TriggerNames const& triggerNames = iEvent.triggerNames(*hltresults);

  if(debug)std::cout << " N triggers:" << ntrigs << std::endl;

  //int passTrigger = 0;
  /*for (int itrig = 0; itrig != ntrigs; ++itrig)
    {
   std::string trigName = triggerNames.triggerName(itrig);
   // std::cout << ">>>>>>>>>>>>>>>>Available Trigger: " << trigName << std::endl;
   bool accept = hltresults->accept(itrig);
   if (!(accept)) continue;
   std::cout << " Accept >>>>>>>>>>>>>>>>>>>>" << trigName << std::endl;
    }
  */
  int hltAccept = -1;
  //std::string trgName = "HLT_PFHT280_QuadPFJet30_v1";
  std::string trgName1 = "HLT_DoubleMediumDeepTauPFTauHPS*";
  std::string trgName2 = "HLT_PFMET*";

  std::string trgPattern = "(" + trgName1 + "|" + trgName2 + ")";
  std::vector< std::vector<std::string>::const_iterator > trgMatches = edm::regexMatch( triggerNames.triggerNames(), trgPattern );
  //std::vector< std::vector<std::string>::const_iterator > trgMatches = edm::regexMatch( triggerNames.triggerNames(), trgName );
  if(debug)std::cout << " N matches: " << trgMatches.size() << std::endl;

  if ( !trgMatches.empty() ) {
    //std::vector<std::string>  HLTPathsByName_;
    //    //std::vector<unsigned int> HLTPathsByIndex_;
    hltAccept = 0;
    for ( auto const& iT : trgMatches ) {
	 if(debug)   std::cout << triggerNames.triggerIndex(*iT)<<"]:"<< *iT << std::endl;
      //HLTPathsByName_.push_back( *iT );
      //HLTPathsByIndex_.push_back( triggerNames.triggerIndex(*iT) );
      if ( hltresults->accept(triggerNames.triggerIndex(*iT)) ){
	    hltAccept = 1;
            if(debug)std::cout << " name["<<triggerNames.triggerIndex(*iT)<<"]:"<< *iT << " -> " << hltresults->accept(triggerNames.triggerIndex(*iT)) << std::endl;
      	  break;
	}
    }
  }
  if(debug)std::cout << "*************** hltAccept:" << hltAccept << std::endl;
  hltAccept_ = hltAccept;
  // Ensure trigger acceptance
  if(hltAccept==1) {
    hNpassed_hlt->Fill(1);
  }
  else hNpassed_hlt->Fill(0);
  if ( hltAccept_ != 1 ) return false;
  
  // Each jet selection must fill vJetIdxs with good jet indices
  
  // Run explicit jet selection
  bool hasPassed;
  if ( task_ == "dijet_ditau" ) {
    hasPassed = runEvtSel_jet_dijet_ditau( iEvent, iSetup );
    if ( debug && hasPassed ) std::cout << "!!!!!!   JET SELECTION HAS PASSED! " << std::endl;
  } 

  if ( !hasPassed ) return false;
  if (task_ == "dijet_ditau" && isMC_) runEvtSel_jet_dijet_ditau_h2aa4Tau( iEvent, iSetup );

  std::sort(vJetIdxs.begin(), vJetIdxs.end());
  if ( debug ) {
    for ( int thisJetIdx : vJetIdxs ) {
      if(debug)std::cout << " index order:" << thisJetIdx << std::endl;
    }
  }

  // edm::ESHandle<CaloGeometry> caloGeomH_;
  // iSetup.get<CaloGeometryRecord>().get( caloGeomH_ );
  // const CaloGeometry* caloGeom = caloGeomH_.product();
  auto const& caloGeom = iSetup.getData(caloGeomToken_);



  edm::Handle<HBHERecHitCollection> HBHERecHitsH_;
  iEvent.getByToken( HBHERecHitCollectionT_, HBHERecHitsH_ );

  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);
  if ( debug ) std::cout << " >> PATJetCol.size: " << jets->size() << std::endl;

  float seedE;
  int iphi_, ieta_, ietaAbs_;

  int nJet = 0;
  vJetSeed_iphi_.clear();
  vJetSeed_ieta_.clear();
  vFailedJetIdx_.clear();
  passedJetIdxs.clear();

  // Loop over jets
  for ( int thisJetIdx : vJetIdxs ) {

    pat::Jet iJet = (*jets)[thisJetIdx];
    
    if ( debug ) std::cout << " >> jet[" << thisJetIdx << "]Pt:" << iJet.pt()  << " Eta:" << iJet.eta()  << " Phi:" << iJet.phi()
			   << " jetE:" << iJet.energy() << " jetM:" << iJet.mass() << std::endl;

    // Get closest HBHE tower to jet position
    // This will not always be the most energetic deposit
    HcalDetId hId( spr::findDetIdHCAL(&caloGeom, iJet.eta(), iJet.phi(), false ) );
    if ( hId.subdet() != HcalBarrel && hId.subdet() != HcalEndcap ){
      vFailedJetIdx_.push_back(thisJetIdx);
      if(debug)std::cout << "Fail getting HBHE tower to jet position" << std::endl;
      continue;
    }
    HBHERecHitCollection::const_iterator iRHit( HBHERecHitsH_->find(hId) );
    seedE = ( iRHit == HBHERecHitsH_->end() ) ? 0. : iRHit->energy() ;
    HcalDetId seedId = hId;
    if ( debug ) std::cout << " >> hId.ieta:" << hId.ieta() << " hId.iphi:" << hId.iphi() << " E:" << seedE << std::endl;

    // Look for the most energetic HBHE tower deposit within a search window
    for ( int ieta = 0; ieta < search_window; ieta++ ) {

      ieta_ = hId.ieta() - (search_window/2)+ieta;

      if ( std::abs(ieta_) > HBHE_IETA_MAX_HE-1 ) continue;
      if ( std::abs(ieta_) < HBHE_IETA_MIN_HB ) continue;

      HcalSubdetector subdet_ = std::abs(ieta_) > HBHE_IETA_MAX_HB ? HcalEndcap : HcalBarrel;

      for ( int iphi = 0; iphi < search_window; iphi++ ) {

        iphi_ = hId.iphi() - (search_window/2)+iphi;

        // iphi should wrap around
        if ( iphi_ > HBHE_IPHI_MAX ) {
          iphi_ = iphi_-HBHE_IPHI_MAX;
        } else if ( iphi_ < HBHE_IPHI_MIN ) {
          iphi_ = HBHE_IPHI_MAX-abs(iphi_);
        }

        // Skip non-existent and lower energy towers
        HcalDetId hId_( subdet_, ieta_, iphi_, 1 );
        HBHERecHitCollection::const_iterator iRHit( HBHERecHitsH_->find(hId_) );
        if ( iRHit == HBHERecHitsH_->end() ) continue;
        if ( iRHit->energy() <= seedE ) continue;
        if ( debug ) std::cout << " !! hId.ieta:" << hId_.ieta() << " hId.iphi:" << hId_.iphi() << " E:" << iRHit->energy() << std::endl;

        seedE = iRHit->energy();
        seedId = hId_;

      } // iphi
    } // ieta

    // NOTE: HBHE iphi = 1 does not correspond to EB iphi = 1!
    // => Need to shift by 2 HBHE towers: HBHE::iphi: [1,...,71,72]->[3,4,...,71,72,1,2]
    iphi_  = seedId.iphi() + 2; // shift
    iphi_  = iphi_ > HBHE_IPHI_MAX ? iphi_-HBHE_IPHI_MAX : iphi_; // wrap-around
    iphi_  = iphi_ - 1; // make histogram-friendly
    ietaAbs_  = seedId.ietaAbs() == HBHE_IETA_MAX_HE ? HBHE_IETA_MAX_HE-1 : seedId.ietaAbs(); // last HBHE ieta embedded
    ieta_  = seedId.zside() > 0 ? ietaAbs_-1 : -ietaAbs_;
    ieta_  = ieta_+HBHE_IETA_MAX_HE-1;

    // If the seed is too close to the edge of HE, discard event
    // Required to keep the seed at the image center
    if ( HBHE_IETA_MAX_HE-1 - ietaAbs_ < image_padding ) {
      if ( debug ) std::cout << " Fail HE edge cut " << std::endl;
      vFailedJetIdx_.push_back(thisJetIdx);
      continue;
    }

    // Save position of most energetic HBHE tower
    // in EB-aligned coordinates
    if ( debug ) std::cout << " !! ieta_:" << ieta_ << " iphi_:" << iphi_ << " ietaAbs_:" << ietaAbs_ << " E:" << seedE << std::endl;
    vJetSeed_iphi_.push_back( iphi_ );
    vJetSeed_ieta_.push_back( ieta_ );
    nJet++;

  } // good jets

  // Remove jets that failed the Seed cuts
  for(int failedJetIdx : vFailedJetIdx_){
    vJetIdxs.erase(std::remove(vJetIdxs.begin(),vJetIdxs.end(),failedJetIdx),vJetIdxs.end());
    if(debug)std::cout << "Failed jets ID:" << failedJetIdx << std::endl;
  }
  if ( vJetIdxs.size() == 0){
    if ( debug ) std::cout << " No passing jets...  " << std::endl;
    if(debug) std::cout << " >> analyze failed: no passing jets" << std::endl;
    return false;
  }

  for (int passJetIdx : vJetIdxs){
    passedJetIdxs.push_back(passJetIdx);
    if(debug)std::cout << "passed jet index is :" << passJetIdx << std::endl;
 }

  if ((nJets_ > 0) && nJet == nJets_) std::cout << " >> analyze failed: " << nJets_ << " passing jets" << std::endl;
  if ( (nJets_ > 0) && nJet != nJets_ ) return false;
  if ( debug ) std::cout << " >> analyze: passed" << std::endl;

  jet_eventId_ = iEvent.id().event();
  jet_runId_ = iEvent.id().run();
  jet_lumiId_ = iEvent.id().luminosityBlock();

  if ( task_ == "dijet_ditau" ) {
    fillEvtSel_jet_dijet_ditau( iEvent, iSetup );
    if(isMC_)fillEvtSel_jet_dijet_ditau_h2aa4Tau( iEvent, iSetup );
  } 

  return true;

} // runEvtSel_jet()
