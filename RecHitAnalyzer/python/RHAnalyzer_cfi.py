

import FWCore.ParameterSet.Config as cms

from RecoMET.METProducers.METSignificanceParams_cfi import METSignificanceParams

fevt = cms.EDAnalyzer('RecHitAnalyzer'
    #, task                           = cms.string("dijet_tau_massregression")
     , task                           = cms.string("dijet_ditau")
    #, task                           = cms.string("tau_classification")
    #, task                           = cms.string("jet_ele_classification")
    #, task                           = cms.string("qcd")
    #, task                           = cms.string("boostedTop")
    , isDebug                        = cms.bool(True)
    , isMC                           = cms.bool(True)
    , isSignal                       = cms.bool(True)
    , isW                            = cms.bool(False)
    , isBoostedTop                   = cms.bool(False)

    #, tracks = cms.untracked.InputTag('ctfWithMaterialTracks')
    #, EBRecHitCollection             = cms.InputTag('ecalRecHit:EcalRecHitsEB')
    , reducedEBRecHitCollection      = cms.InputTag('ecalRecHit:EcalRecHitsEB')
    #, reducedEBRecHitCollection      = cms.InputTag('reducedEcalRecHitsEB')
    #, EERecHitCollection             = cms.InputTag('ecalRecHit:EcalRecHitsEE')
    , reducedEERecHitCollection      = cms.InputTag('ecalRecHit:EcalRecHitsEE')
    #, reducedEERecHitCollection      = cms.InputTag('reducedEcalRecHitsEE')
    #, EBDigiCollection               = cms.InputTag('simEcalDigis:ebDigis')
    #, selectedEBDigiCollection       = cms.InputTag('selectDigi:selectedEcalEBDigiCollection')
    , reducedHBHERecHitCollection    = cms.InputTag('hbhereco')
    # , reducedHBHERecHitCollection    = cms.InputTag('reducedHcalRecHits:hbhereco')
    , genParticleCollection          = cms.InputTag('prunedGenParticles')
    , ak4PFJetCollection             = cms.InputTag('slimmedJets')
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
    #, eleCollection                  = cms.InputTag("gedGsfElectrons")
    #, muonCollection                 = cms.InputTag("muons")
    , tauCollection                  = cms.InputTag("slimmedTaus")

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
    , maxJetEta = cms.double(2.4)
    , z0PVCut   = cms.double(0.1)

    # MET parameter
    , parameters = METSignificanceParams

    #granularity multiplier wrt ECAL maps for tracker and tracking RH images
    , granularityMultiPhi = cms.int32(5)
    , granularityMultiEta = cms.int32(5)
    )
