

import FWCore.ParameterSet.Config as cms

from RecoMET.METProducers.METSignificanceParams_cfi import METSignificanceParams

fevt = cms.EDAnalyzer('RecHitAnalyzer'
    , task                           = cms.string("tau_massregression")
    # , task                           = cms.string("dijet_ditau")
    #, task                           = cms.string("tau_classification")
    #, task                           = cms.string("jet_ele_classification")
    #, task                           = cms.string("qcd")
    #, task                           = cms.string("boostedTop")
    , isDebug                        = cms.bool(True)
    , isMC                           = cms.bool(True)
    , isSignal                       = cms.bool(False)
    , isW                            = cms.bool(False)
    , isBoostedTop                   = cms.bool(False)

    #, tracks = cms.untracked.InputTag('ctfWithMaterialTracks')
    , EBRecHitCollection             = cms.InputTag('ecalRecHit:EcalRecHitsEB') #RAW
    , reducedEBRecHitAODCollection   = cms.InputTag('reducedEcalRecHitsEB') #AOD
    , reducedEBRecHitminiAODCollection = cms.InputTag('reducedEgamma:reducedEBRecHits')#miniAOD
    , EERecHitCollection             = cms.InputTag('ecalRecHit:EcalRecHitsEE')#RAW
    , reducedEERecHitAODCollection   = cms.InputTag('reducedEcalRecHitsEE')#AOD
    , reducedEERecHitminiAODCollection = cms.InputTag('reducedEgamma:reducedEERecHits')#miniAOD
    
    , HBHERecHitCollection    = cms.InputTag('hbhereco')#RAW
    , reducedHBHERecHitAODCollection    = cms.InputTag('reducedHcalRecHits:hbhereco')#AOD
    , reducedHBHERecHitminiAODCollection= cms.InputTag('reducedEgamma:reducedHBHEHits')#miniAOD reduced Egamma HCAL rechits
    #, reducedHBHERecHitminiAODCollection= cms.InputTag('slimmedHcalRecHits:reducedHcalRecHits')#miniAOD
    , genParticleCollection          = cms.InputTag('prunedGenParticles')
    , ak4PFJetCollection             = cms.InputTag('slimmedJetsPuppi')
    , genJetCollection               = cms.InputTag('slimmedGenJets')
    , trackRecHitCollection          = cms.InputTag('generalTracks')
    , trackCollection                = cms.InputTag("generalTracks")
    , vertexCollection               = cms.InputTag("offlineSlimmedPrimaryVertices")
    , secVertexCollection            = cms.InputTag("slimmedSecondaryVertices")
    , siPixelRecHitCollection        = cms.InputTag("siPixelRecHits")
    , siStripMatchedRecHitCollection = cms.InputTag("siStripMatchedRecHits", "matchedRecHit")
    , siStripRphiRecHits             = cms.InputTag("siStripMatchedRecHits", "rphiRecHit")
    , siStripStereoRecHits           = cms.InputTag("siStripMatchedRecHits", "stereoRecHit")
    #, pfCollection                   = cms.InputTag("particleFlow")
    #, srcPFCandidates                = cms.InputTag("particleFlow")
    #, srcPfJets                      = cms.InputTag("ak4PFJets")
    , metCollection                  = cms.InputTag("slimmedMETs")
    , tauCollection                  = cms.InputTag("slimmedTaus")
    , eleCollection                  = cms.InputTag("slimmedElectrons")
    , muonCollection                 = cms.InputTag("slimmedMuons")

    , triggerResultsTag              = cms.InputTag("TriggerResults", "", "HLT")
    , triggerSummaryTag              = cms.InputTag("hltTriggerSummaryAOD","","HLT")
    #, ipTagInfoCollection            = cms.InputTag("pfImpactParameterTagInfos")
    , mode                           = cms.string("JetLevel")
    #, rhoLabel                       = cms.InputTag('fixedGridRhoAll')
    #, srcJetSF                       = cms.string('AK4PFchs')
    #, srcJetResPt                    = cms.string('AK4PFchs_pt')
    #, srcJetResPhi                   = cms.string('AK4PFchs_phi')
    #, srcLeptons                     = cms.VInputTag("slimmedElectrons","slimmedMuons","slimmedPhotons"),
    #, srcLeptons                     = cms.VInputTag("gedGsfElectrons","muons","gedPhotons")
    , transTrackBuilder              = cms.ESInputTag("", "TransientTrackBuilder")

    # Jet level cfg
    , nJets     = cms.int32(-1)
    , minJetPt  = cms.double(20.)
    , maxJetEta = cms.double(3.0)
    , z0PVCut   = cms.double(0.1)

    # MET parameter
    , parameters = METSignificanceParams

    #granularity multiplier wrt ECAL maps for tracker and tracking RH images
    , granularityMultiPhi = cms.int32(5)
    , granularityMultiEta = cms.int32(5)
    )
