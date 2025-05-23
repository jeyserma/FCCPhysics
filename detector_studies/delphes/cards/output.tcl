
module EDM4HepOutput EDM4HepOutput {   
    add ReconstructedParticleCollections EFlowTrack EFlowPhoton EFlowNeutralHadron
    add GenParticleCollections           Particle
    add JetCollections                   Jet
    add MuonCollections                  Muon 
    add ElectronCollections              Electron
    add PhotonCollections                Photon
    set RecoParticleCollectionName       ReconstructedParticles
    set MCRecoAssociationCollectionName  MCRecoAssociations
} 