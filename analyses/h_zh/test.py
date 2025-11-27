
import ROOT
ROOT.TH1.SetDefaultSumw2(ROOT.kTRUE)



processList = {
    'wz3p6_ee_qqH_Hinv_ecm240':         {'fraction':1},
    'wz3p6_ee_qqH_Hinv_ecm365':         {'fraction':1},

    'wz3p6_ee_ccH_Hinv_ecm240':         {'fraction':1},
    'wz3p6_ee_ccH_Hinv_ecm365':         {'fraction':1},


    'wzp6_ee_numunumuH_Huu_ecm240':         {'fraction':1},
    'wzp6_ee_numunumuH_Huu_ecm365':         {'fraction':1},
}




inputDir = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"
procDict = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/samplesDict.json"

# additional/custom C++ functions
includePaths = ["../../functions/functions.h", "../../functions/functions_gen.h"]


# output directory
outputDir   = f"output/test/histmaker//"

nCPUS       = 8

# scale the histograms with the cross-section and integrated luminosity
doScale = True
intLumi = 1

# define histograms
bins_m = (500, 0, 500) # 100 MeV bins
bins_maa = (100, 120, 130) # 100 MeV bins
bins_p = (500, 0, 500) # 100 MeV bins
bins_m_zoom = (100, 120, 130) # 100 MeV

bins_theta = (500, 0, 5)
bins_phi = (400, -4, 4)

bins_count = (100, 0, 100)
bins_pdgid = (60, -30, 30)
bins_charge = (10, -5, 5)

bins_resolution = (10000, 0.95, 1.05)
bins_resolution_1 = (20000, 0, 2)

jet_energy = (5000, 0, 500) # 100 MeV bins
dijet_m = (5000, 0, 500) # 100 MeV bins
visMass = (5000, 0, 500) # 100 MeV bins
missEnergy  = (5000, 0, 500) # 100 MeV bins

dijet_m_final = (500, 50, 100) # 100 MeV bins

bins_cos = (100, -1, 1)
bins_cos_abs = (100, 0, 1)
bins_iso = (1000, 0, 10)
bins_aco = (1000,0,5)
bins_cosThetaMiss = (10000, 0, 1)

bins_m_fine = (500, 110, 130) # 100 MeV bins


def build_graph(df, dataset):

    hists, cols = [], []

    df = df.Define("weight", "1.0")
    weightsum = df.Sum("weight")




    # define collections
    df = df.Alias("Particle0", "Particle#0.index")
    df = df.Alias("Particle1", "Particle#1.index")
    df = df.Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
    df = df.Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")

    # all photons
    df = df.Alias("Photon0", "Photon#0.index")
    df = df.Define("photons_all", "FCCAnalyses::ReconstructedParticle::get(Photon0, ReconstructedParticles)")
    df = df.Define("photons_all_p", "FCCAnalyses::ReconstructedParticle::get_p(photons_all)")
    df = df.Define("photons_all_theta", "FCCAnalyses::ReconstructedParticle::get_theta(photons_all)")
    df = df.Define("photons_all_no", "FCCAnalyses::ReconstructedParticle::get_n(photons_all)")

    hists.append(df.Histo1D(("photons_all_p", "", *bins_p), "photons_all_p"))
    hists.append(df.Histo1D(("photons_all_no", "", *bins_count), "photons_all_no"))

    df = df.Define("photons", "FCCAnalyses::sel_range(5, 170, false)(photons_all, photons_all_p)")
    df = df.Define("photons_p", "FCCAnalyses::ReconstructedParticle::get_p(photons)")
    df = df.Define("photons_theta", "FCCAnalyses::ReconstructedParticle::get_theta(photons)")
    df = df.Define("photons_n", "FCCAnalyses::ReconstructedParticle::get_n(photons)")

    hists.append(df.Histo1D(("photons_n", "", *bins_count), "photons_n"))

    # define PF candidates collection by removing the muons
    #df = df.Define("rps_no_photons", "FCCAnalyses::ReconstructedParticle::remove(ReconstructedParticles, photons)")
    df = df.Alias("rps_no_photons", "ReconstructedParticles")
    df = df.Define("RP_px", "FCCAnalyses::ReconstructedParticle::get_px(rps_no_photons)")
    df = df.Define("RP_py", "FCCAnalyses::ReconstructedParticle::get_py(rps_no_photons)")
    df = df.Define("RP_pz","FCCAnalyses::ReconstructedParticle::get_pz(rps_no_photons)")
    df = df.Define("RP_e", "FCCAnalyses::ReconstructedParticle::get_e(rps_no_photons)")
    df = df.Define("RP_m", "FCCAnalyses::ReconstructedParticle::get_mass(rps_no_photons)")
    df = df.Define("RP_q", "FCCAnalyses::ReconstructedParticle::get_charge(rps_no_photons)")
    df = df.Define("pseudo_jets", "FCCAnalyses::JetClusteringUtils::set_pseudoJets(RP_px, RP_py, RP_pz, RP_e)")

    df = df.Define("clustered_jets", "JetClustering::clustering_ee_kt(2, 2, 1, 0)(pseudo_jets)")
    df = df.Define("jets", "FCCAnalyses::JetClusteringUtils::get_pseudoJets(clustered_jets)")
    df = df.Define("jetconstituents", "FCCAnalyses::JetClusteringUtils::get_constituents(clustered_jets)")
    df = df.Define("jets_e", "FCCAnalyses::JetClusteringUtils::get_e(jets)")
    df = df.Define("jets_px", "FCCAnalyses::JetClusteringUtils::get_px(jets)")
    df = df.Define("jets_py", "FCCAnalyses::JetClusteringUtils::get_py(jets)")
    df = df.Define("jets_pz", "FCCAnalyses::JetClusteringUtils::get_pz(jets)")
    df = df.Define("jets_m", "FCCAnalyses::JetClusteringUtils::get_m(jets)")

    df = df.Define("jet1", "ROOT::Math::PxPyPzEVector(jets_px[0], jets_py[0], jets_pz[0], jets_e[0])")
    df = df.Define("jet2", "ROOT::Math::PxPyPzEVector(jets_px[1], jets_py[1], jets_pz[1], jets_e[1])")
    df = df.Define("jet1_p", "jet1.P()")
    df = df.Define("jet2_p", "jet2.P()")
    df = df.Define("jet1_phi", "jet1.Phi()")
    df = df.Define("jet2_phi", "jet2.Phi()")
    df = df.Define("jet1_eta", "jet1.Eta()")
    df = df.Define("jet2_eta", "jet2.Eta()")
    df = df.Define("jet1_theta", "jet1.Theta()")
    df = df.Define("jet2_theta", "jet2.Theta()")
    df = df.Define("dijet", "jet1+jet2")
    df = df.Define("dijet_tlv", "TLorentzVector ret; ret.SetPxPyPzE(dijet.Px(), dijet.Py(), dijet.Pz(), dijet.E()); return ret;")
    df = df.Define("dijet_m", "dijet.M()")
    df = df.Define("dijet_p", "dijet.P()")
    df = df.Define("dijet_theta", "dijet.Theta()")
    df = df.Define("costheta1", "abs(cos(dijet.Theta()))")

    df = df.Define("dphi", "ROOT::Math::VectorUtil::DeltaPhi(jet1, jet2);")
    df = df.Define("deta", "jet1_eta-jet2_eta")
    df = df.Define("dr", "std::sqrt(deta * deta + dphi * dphi)")

    hists.append(df.Histo1D(("dphi", "", *(100, -4, 4)), "dphi"))
    hists.append(df.Histo1D(("deta", "", *(100, -4, 4)), "deta"))
    hists.append(df.Histo1D(("dr", "", *(100,0,10)), "dr"))

    hists.append(df.Histo1D(("dijet_m", "", *bins_m), "dijet_m"))
    hists.append(df.Histo1D(("dijet_p", "", *bins_p), "dijet_p"))

    hists.append(df.Histo1D(("zqq_jet1_p_nOne", "", *bins_p), "jet1_p"))
    hists.append(df.Histo1D(("zqq_jet2_p_nOne", "", *bins_p), "jet2_p"))


    return hists, weightsum




