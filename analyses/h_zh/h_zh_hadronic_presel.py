
import ROOT
ROOT.TH1.SetDefaultSumw2(ROOT.kTRUE)

ecm = 240 # 240 365
flavor = "mumu" # ee mumu
analysisType = "mass"

# list of all processes
fraction = 0.1
processList = {

    'wzp6_ee_qqH_Hbb_ecm240':           {'fraction':fraction},
    'wzp6_ee_qqH_Hcc_ecm240':           {'fraction':fraction},
    'wzp6_ee_qqH_Hss_ecm240':           {'fraction':fraction},
    'wzp6_ee_qqH_Hgg_ecm240':           {'fraction':fraction},
    'wzp6_ee_qqH_Haa_ecm240':           {'fraction':fraction},
    'wzp6_ee_qqH_HZa_ecm240':           {'fraction':fraction},
    'wzp6_ee_qqH_HWW_ecm240':           {'fraction':fraction},
    'wzp6_ee_qqH_HZZ_ecm240':           {'fraction':fraction},
    'wzp6_ee_qqH_Hmumu_ecm240':          {'fraction':fraction},
    'wzp6_ee_qqH_Htautau_ecm240':        {'fraction':fraction},
    'wz3p6_ee_qqH_Hinv_ecm240':        {'fraction':fraction},

    'wzp6_ee_ssH_Hbb_ecm240':           {'fraction':fraction},
    'wzp6_ee_ssH_Hcc_ecm240':           {'fraction':fraction},
    'wzp6_ee_ssH_Hss_ecm240':           {'fraction':fraction},
    'wzp6_ee_ssH_Hgg_ecm240':           {'fraction':fraction},
    'wzp6_ee_ssH_Haa_ecm240':           {'fraction':fraction},
    'wzp6_ee_ssH_HZa_ecm240':           {'fraction':fraction},
    'wzp6_ee_ssH_HWW_ecm240':           {'fraction':fraction},
    'wzp6_ee_ssH_HZZ_ecm240':           {'fraction':fraction},
    'wzp6_ee_ssH_Hmumu_ecm240':          {'fraction':fraction},
    'wzp6_ee_ssH_Htautau_ecm240':        {'fraction':fraction},
    'wz3p6_ee_ssH_Hinv_ecm240':        {'fraction':fraction},

    'wzp6_ee_ccH_Hbb_ecm240':           {'fraction':fraction},
    'wzp6_ee_ccH_Hcc_ecm240':           {'fraction':fraction},
    'wzp6_ee_ccH_Hss_ecm240':           {'fraction':fraction},
    'wzp6_ee_ccH_Hgg_ecm240':           {'fraction':fraction},
    'wzp6_ee_ccH_Haa_ecm240':           {'fraction':fraction},
    'wzp6_ee_ccH_HZa_ecm240':           {'fraction':fraction},
    'wzp6_ee_ccH_HWW_ecm240':           {'fraction':fraction},
    'wzp6_ee_ccH_HZZ_ecm240':           {'fraction':fraction},
    'wzp6_ee_ccH_Hmumu_ecm240':          {'fraction':fraction},
    'wzp6_ee_ccH_Htautau_ecm240':        {'fraction':fraction},
    'wz3p6_ee_ccH_Hinv_ecm240':        {'fraction':fraction},


    'wzp6_ee_bbH_Hbb_ecm240':           {'fraction':fraction},
    'wzp6_ee_bbH_Hcc_ecm240':           {'fraction':fraction},
    'wzp6_ee_bbH_Hss_ecm240':           {'fraction':fraction},
    'wzp6_ee_bbH_Hgg_ecm240':           {'fraction':fraction},
    'wzp6_ee_bbH_Haa_ecm240':           {'fraction':fraction},
    'wzp6_ee_bbH_HZa_ecm240':           {'fraction':fraction},
    'wzp6_ee_bbH_HWW_ecm240':           {'fraction':fraction},
    'wzp6_ee_bbH_HZZ_ecm240':           {'fraction':fraction},
    'wzp6_ee_bbH_Hmumu_ecm240':          {'fraction':fraction},
    'wzp6_ee_bbH_Htautau_ecm240':        {'fraction':fraction},
    'wz3p6_ee_bbH_Hinv_ecm240':        {'fraction':fraction},


    'p8_ee_WW_ecm240':            {'fraction':fraction},
    'p8_ee_ZZ_ecm240':            {'fraction':fraction},
    'p8_ee_Zqq_ecm240':           {'fraction':fraction},
}




inputDir = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"
procDict = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/samplesDict.json"

# additional/custom C++ functions
includePaths = ["../../functions/functions.h", "../../functions/functions_gen.h", "utils.h", "utils_hadronic.h"]


# output directory
outputDir   = "output/h_zh_hadronic/treemaker/"


# optional: ncpus, default is 4, -1 uses all cores available
nCPUS       = 32

# scale the histograms with the cross-section and integrated luminosity
doScale = True
intLumi = 10.8e6

# define histograms

# define histograms
bins_p_mu = (2000, 0, 200) # 100 MeV bins
bins_m_ll = (2000, 0, 200) # 100 MeV bins
bins_p_ll = (200, 0, 200) # 1 GeV bins
bins_recoil = (20000, 0, 200) # 10 MeV bins 
bins_recoil_fine = (20000, 120, 140) # 1 MeV bins 
bins_cosThetaMiss = (10000, 0, 1)

bins_theta = (500, 0, 5)
bins_phi = (500, -5, 5)
bins_aco = (1000, -10, 10)

bins_count = (50, 0, 50)
bins_pdgid = (60, -30, 30)
bins_charge = (10, -5, 5)
bins_iso = (500, 0, 5)
bins_dR = (1000, 0, 10)

bins_cat = (10, 0, 10)
bins_resolution = (10000, 0.95, 1.05)

bins_m_fine = (100, 120, 130) # 100 MeV bins
bins_cos_abs = (100, 0, 1)

def exclusive_clustering(df, njets):
    if njets==0: # inclusive
        df = df.Define("clustered_jets_N0", "JetClustering::clustering_kt(0.6, 0, 5, 1, 0)(pseudo_jets)")
    else:
        df = df.Define(f"clustered_jets_N{njets}", f"JetClustering::clustering_ee_kt(2, {njets}, 0, 10)(pseudo_jets)")
        #df = df.Define(f"clustered_jets_N{njets}", f"JetClustering::clustering_ee_genkt(0.5, 2, {njets}, 0, 10, 0)(pseudo_jets)")
    df = df.Define(f"jets_N{njets}", f"FCCAnalyses::JetClusteringUtils::get_pseudoJets(clustered_jets_N{njets})")
    df = df.Define(f"njets_N{njets}", f"jets_N{njets}.size()")
    df = df.Define(f"jetconstituents_N{njets}", f"FCCAnalyses::JetClusteringUtils::get_constituents(clustered_jets_N{njets})")
    df = df.Define(f"jets_e_N{njets}", f"FCCAnalyses::JetClusteringUtils::get_e(jets_N{njets})")
    df = df.Define(f"jets_px_N{njets}", f"FCCAnalyses::JetClusteringUtils::get_px(jets_N{njets})")
    df = df.Define(f"jets_py_N{njets}", f"FCCAnalyses::JetClusteringUtils::get_py(jets_N{njets})")
    df = df.Define(f"jets_pz_N{njets}", f"FCCAnalyses::JetClusteringUtils::get_pz(jets_N{njets})")
    df = df.Define(f"jets_m_N{njets}", f"FCCAnalyses::JetClusteringUtils::get_m(jets_N{njets})")
    df = df.Define(f"jets_rp_N{njets}", f"FCCAnalyses::jets2rp(jets_px_N{njets}, jets_py_N{njets}, jets_pz_N{njets}, jets_e_N{njets}, jets_m_N{njets})")
    df = df.Define(f"jets_rp_cand_N{njets}", f"FCCAnalyses::select_jets(jets_rp_N{njets}, jetconstituents_N{njets}, {njets}, ReconstructedParticles)") # reduces potentially the jet multiplicity
    df = df.Define(f"njets_cand_N{njets}", f"jets_rp_cand_N{njets}.size()")
    return df

def define_variables(df, njets):
    #if njets == 0:
    #df = df.Filter(f"if(njets_cand_N{njets} < 2) std::cout << \"FILTER={njets} \" << njets_cand_N{njets}  <<  std::endl; return njets_cand_N{njets} >= 2") #### tricky need to fix it

    df = df.Define(f"zbuilder_result_N{njets}", f"FCCAnalyses::resonanceBuilder_mass_recoil_hadronic(91.2, 125, 0.0, ecm)(jets_rp_cand_N{njets})")
    df = df.Define(f"zqq_N{njets}", f"ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{{zbuilder_result_N{njets}[0]}}") # the Z
    df = df.Define(f"zqq_jets_N{njets}", f"ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{{zbuilder_result_N{njets}[1],zbuilder_result_N{njets}[2]}}") # the leptons
    df = df.Define(f"zqq_m_N{njets}", f"FCCAnalyses::ReconstructedParticle::get_mass(zqq_N{njets})[0]")
    df = df.Define(f"zqq_p_N{njets}", f"FCCAnalyses::ReconstructedParticle::get_p(zqq_N{njets})[0]")
    df = df.Define(f"zqq_recoil_N{njets}", f"FCCAnalyses::ReconstructedParticle::recoilBuilder(ecm)(zqq_N{njets})")
    df = df.Define(f"zqq_recoil_m_N{njets}", f"FCCAnalyses::ReconstructedParticle::get_mass(zqq_recoil_N{njets})[0]")


    return df


def veto_selection_leptonic(df, flavor):

    if flavor == "muon":
        df = df.Define(f"leps_{flavor}", "FCCAnalyses::ReconstructedParticle::sel_p(20)(muons_all)")
    else:
        df = df.Define(f"leps_{flavor}", "FCCAnalyses::ReconstructedParticle::sel_p(20)(electrons_all)")

    df = df.Define(f"leps_{flavor}_q", f"FCCAnalyses::ReconstructedParticle::get_charge(leps_{flavor})")
    df = df.Define(f"leps_{flavor}_no", f"FCCAnalyses::ReconstructedParticle::get_n(leps_{flavor})")
    df = df.Define(f"leps_{flavor}_iso", f"FCCAnalyses::coneIsolation(0.01, 0.5)(leps_{flavor}, ReconstructedParticles)")
    df = df.Define(f"leps_{flavor}_sel_iso", f"FCCAnalyses::sel_iso(0.25)(leps_{flavor}, leps_{flavor}_iso)") # 0.25




    df = df.Define(f"zbuilder_result_{flavor}", f"FCCAnalyses::resonanceBuilder_mass_recoil(91.2, 125, 0.4, ecm, false)(leps_{flavor}, MCRecoAssociations0, MCRecoAssociations1, ReconstructedParticles, Particle, Particle0, Particle1)")
    df = df.Define(f"zll_{flavor}", f"ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{{zbuilder_result_{flavor}[0]}}") # the Z
    df = df.Define(f"zll_leps_{flavor}", f"ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{{zbuilder_result_{flavor}[1],zbuilder_result_{flavor}[2]}}") # the leptons
    df = df.Define(f"zll_m_{flavor}", f"FCCAnalyses::ReconstructedParticle::get_mass(zll_{flavor})[0]")
    df = df.Define(f"zll_p_{flavor}", f"FCCAnalyses::ReconstructedParticle::get_p(zll_{flavor})[0]")
    df = df.Define(f"zll_recoil_{flavor}", f"FCCAnalyses::ReconstructedParticle::recoilBuilder(ecm)(zll_{flavor})")
    df = df.Define(f"zll_recoil_m_{flavor}", f"FCCAnalyses::ReconstructedParticle::get_mass(zll_recoil_{flavor})[0]")


    #df = df.Define(f"missingEnergy_{flavor}", f"FCCAnalyses::missingEnergy(ecm, ReconstructedParticles)")
    #df = df.Define(f"cosTheta_miss_{flavor}", f"FCCAnalyses::get_cosTheta_miss(missingEnergy_{flavor})")


    sel_leps = f"leps_{flavor}_no >= 2 && leps_{flavor}_sel_iso.size() > 0 && abs(Sum(leps_{flavor}_q)) < leps_{flavor}_q.size()"
    sel_mll = f"zll_m_{flavor} > 86 && zll_m_{flavor} < 96"
    sel_pll = f"zll_p_{flavor} > 20 && zll_p_{flavor} < 70"
    sel_recoil = f"zll_recoil_m_{flavor} < 140 && zll_recoil_m_{flavor} > 120"
    
    #df = df.Filter(f"!({sel_leps} && {sel_mll} && {sel_pll} && {sel_recoil})")
    df = df.Filter(f"!({sel_leps} && {sel_mll} && {sel_pll})")
    #df = df.Filter(f"!(({sel_leps}) && {sel_mll})")
    return df
    #df = df.Filter("cosTheta_miss < 0.98")






class RDFanalysis():

    # encapsulate analysis logic, definitions and filters in the dataframe
    def analysers(df):

        hists = []
        sigProcs = ["wzp6_ee_mumuH_ecm240", "wzp6_ee_eeH_ecm240"]

        df = df.Define("ecm", "240" if ecm == 240 else "365")
        df = df.Define("weight", "1.0")
        weightsum = df.Sum("weight")

        df = df.Define("cut0", "0")
        df = df.Define("cut1", "1")
        df = df.Define("cut2", "2")
        df = df.Define("cut3", "3")
        df = df.Define("cut4", "4")
        df = df.Define("cut5", "5")
        df = df.Define("cut6", "6")
        df = df.Define("cut7", "7")
        df = df.Define("cut8", "8")
        df = df.Define("cut9", "9")
        df = df.Define("cut10", "10")
        df = df.Define("cut11", "11")
        df = df.Define("cut12", "12")
        df = df.Define("cut13", "13")
        df = df.Define("cut14", "14")

        hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut0")) ## all events

        df = df.Alias("Particle0", "Particle#0.index")
        df = df.Alias("Particle1", "Particle#1.index")
        df = df.Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
        df = df.Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")

        df = df.Alias("Muon0", "Muon#0.index")
        df = df.Alias("Electron0", "Electron#0.index")
        df = df.Alias("Photon0", "Photon#0.index")
        df = df.Define("muons_all", "FCCAnalyses::ReconstructedParticle::get(Muon0, ReconstructedParticles)")
        df = df.Define("electrons_all", "FCCAnalyses::ReconstructedParticle::get(Electron0, ReconstructedParticles)")
        df = df.Define("photons_all", "FCCAnalyses::ReconstructedParticle::get(Photon0, ReconstructedParticles)")

        ## orthogonality w.r.t. leptonic channels
        df = veto_selection_leptonic(df, "muon")
        df = veto_selection_leptonic(df, "electron")
        hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut1"))







        #df = df.Filter("ReconstructedParticles.size() > 30")
        
        #### remove isolated photons from clustering
        df = df.Define("photons_all_p", "FCCAnalyses::ReconstructedParticle::get_p(photons_all)")
        df = df.Define("photons", "FCCAnalyses::sel_range(40, 95, false)(photons_all, photons_all_p)")
        df = df.Define("rps_no_photons", "FCCAnalyses::ReconstructedParticle::remove(ReconstructedParticles, photons)")


        #### remove isolated muons from clustering
        # ensure orthogonality with leptonic channel
        df = df.Define("muons_all_p", "FCCAnalyses::ReconstructedParticle::get_p(muons_all)")
        df = df.Define("muons", "FCCAnalyses::sel_range(40, 95, false)(muons_all, muons_all_p)")

        #### remove isolated electrons from clustering
        # ensure orthogonality with leptonic channel
        df = df.Define("electrons_all_p", "FCCAnalyses::ReconstructedParticle::get_p(electrons_all)")
        df = df.Define("electrons", "FCCAnalyses::sel_range(40, 95, false)(electrons_all, electrons_all_p)")


        ### remove isolated electrons from clustering
        # ensure orthogonality with leptonic channel
        df = df.Define("rps_no_photons_muons", "FCCAnalyses::ReconstructedParticle::remove(rps_no_photons, muons)")
        df = df.Define("rps_no_photons_muons_electrons", "FCCAnalyses::ReconstructedParticle::remove(rps_no_photons_muons, electrons)")


        # define PF candidates collection by removing the muons
        df = df.Alias("rps_sel", "rps_no_photons_muons_electrons")
        df = df.Define("RP_px", "FCCAnalyses::ReconstructedParticle::get_px(rps_sel)")
        df = df.Define("RP_py", "FCCAnalyses::ReconstructedParticle::get_py(rps_sel)")
        df = df.Define("RP_pz","FCCAnalyses::ReconstructedParticle::get_pz(rps_sel)")
        df = df.Define("RP_e", "FCCAnalyses::ReconstructedParticle::get_e(rps_sel)")
        df = df.Define("RP_m", "FCCAnalyses::ReconstructedParticle::get_mass(rps_sel)")
        df = df.Define("RP_q", "FCCAnalyses::ReconstructedParticle::get_charge(rps_sel)")
        df = df.Define("pseudo_jets", "FCCAnalyses::JetClusteringUtils::set_pseudoJets(RP_px, RP_py, RP_pz, RP_e)")


        # perform possible clusterings
        
        df = exclusive_clustering(df, 6)
        df = exclusive_clustering(df, 4)
        df = exclusive_clustering(df, 2)
        df = exclusive_clustering(df, 0)

        df = define_variables(df, 6)
        df = define_variables(df, 4)
        df = define_variables(df, 2)
        df = define_variables(df, 0)


        df = df.Define("zqq", "std::vector<Vec_rp> r = {zqq_N0, zqq_N2, zqq_N4, zqq_N6}; return r;")
        df = df.Define("zqq_jets", "std::vector<Vec_rp> r = {jets_rp_cand_N0, jets_rp_cand_N2, jets_rp_cand_N4, jets_rp_cand_N6}; return r;")
        df = df.Define("zqq_m", "Vec_f r = {zqq_m_N0, zqq_m_N2, zqq_m_N4, zqq_m_N6}; return r;")
        df = df.Define("zqq_p", "Vec_f r = {zqq_p_N0, zqq_p_N2, zqq_p_N4, zqq_p_N6}; return r;")
        df = df.Define("zqq_recoil_m", "Vec_f r = {zqq_recoil_m_N0, zqq_recoil_m_N2, zqq_recoil_m_N4, zqq_recoil_m_N6}; return r;")
        df = df.Define("njets", "Vec_i r = {(int)njets_cand_N0, (int)njets_cand_N2, (int)njets_cand_N4, (int)njets_cand_N6}; return r;")
        df = df.Define("njets_target", "Vec_i r = {0, 2, 4, 6}; return r;")
        df = df.Define("best_clustering_idx", "FCCAnalyses::best_clustering_idx(zqq_m, zqq_p, zqq_recoil_m, njets, njets_target)")

        hists.append(df.Histo1D(("best_clustering_idx_nosel", "", *(15, -5, 10)), "best_clustering_idx"))

        ## njets for inclusive clustering
        hists.append(df.Histo1D(("njets_inclusive", "", *bins_count), "njets_cand_N0"))
        df_incl = df.Filter("best_clustering_idx == 0")
        hists.append(df_incl.Histo1D(("njets_inclusive_sel", "", *bins_count), "njets_cand_N0"))




        df = df.Filter("best_clustering_idx >= 0")
        hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut2")) ## after clustering

        df = df.Define("zqq_best", "zqq[best_clustering_idx]")
        df = df.Define("zqq_jets_best", "zqq_jets[best_clustering_idx]")
        df = df.Define("zqq_m_best", "zqq_m[best_clustering_idx]")
        df = df.Define("zqq_recoil_m_best", "zqq_recoil_m[best_clustering_idx]")
        df = df.Define("zqq_p_best", "zqq_p[best_clustering_idx]")
        df = df.Define("zqq_jets_p_best", "FCCAnalyses::ReconstructedParticle::get_p(zqq_jets_best)")
        df = df.Define("zqq_jets_theta_best", "FCCAnalyses::ReconstructedParticle::get_theta(zqq_jets_best)")
        df = df.Define("zqq_jets_costheta_best", "Vec_f ret; for(auto & theta: zqq_jets_theta_best) ret.push_back(std::abs(cos(theta))); return ret;")

        df = df.Define("z_theta", "FCCAnalyses::ReconstructedParticle::get_theta(zqq_best)")
        df = df.Define("z_costheta", "std::abs(cos(z_theta[0]))")
        hists.append(df.Histo1D(("z_costheta_nosel", "", *bins_cos_abs), "z_costheta"))
        

        hists.append(df.Histo1D(("zqq_m_best_nosel", "", *(200, 0, 200)), "zqq_m_best"))
        hists.append(df.Histo1D(("zqq_p_best_nosel", "", *(200, 0, 200)), "zqq_p_best"))
        hists.append(df.Histo1D(("zqq_recoil_m_best_nosel", "", *(200, 0, 200)), "zqq_recoil_m_best"))
        hists.append(df.Histo2D(("zqq_recoil_m_mqq_nosel", "", *((25, 100, 150)+(35, 50, 120))), "zqq_recoil_m_best", "zqq_m_best"))

        ## jet kinematics
        df = df.Define("leading_idx", "(zqq_jets_best[0].energy > zqq_jets_best[1].energy) ? 0 : 1")
        df = df.Define("subleading_idx", "(zqq_jets_best[0].energy > zqq_jets_best[1].energy) ? 1 : 0")
        df = df.Define("leading_jet_p", "zqq_jets_p_best[leading_idx]")
        df = df.Define("subleading_jet_p", "zqq_jets_p_best[subleading_idx]")
        df = df.Define("leading_jet_costheta", "zqq_jets_costheta_best[leading_idx]")
        df = df.Define("subleading_jet_costheta", "zqq_jets_costheta_best[subleading_idx]")
        
        hists.append(df.Histo1D(("leading_jet_p", "", *(200, 0, 200)), "leading_jet_p"))
        hists.append(df.Histo1D(("subleading_jet_p", "", *(200, 0, 200)), "subleading_jet_p"))
        hists.append(df.Histo1D(("leading_jet_costheta", "", *bins_cos_abs), "leading_jet_costheta"))
        hists.append(df.Histo1D(("subleading_jet_costheta", "", *bins_cos_abs), "subleading_jet_costheta"))

        ## attempt to reconstruct WW with 4 jets
        df = df.Define("ww_pairs", "FCCAnalyses::pair_W(jets_rp_cand_N4)")
        df = df.Define("w1", "ww_pairs[0]")
        df = df.Define("w2", "ww_pairs[1]")
        df = df.Define("dr_mww", "std::sqrt((w1-78)*(w1-78) + (w2-78)*(w2-78))")
        # dphi between W candidates
        df = df.Define("dphi_mww", "FCCAnalyses::pair_W_dphi(jets_rp_cand_N4)")
        
        df = df.Define("ww_pairs_p", "FCCAnalyses::pair_W_p(jets_rp_cand_N4)")
        df = df.Define("w1_p", "ww_pairs_p[0]")
        df = df.Define("w2_p", "ww_pairs_p[1]")
        
        hists.append(df.Histo1D(("dphi_mww", "", *(400, 0, 4)), "dphi_mww"))
        hists.append(df.Histo1D(("dr_mww", "", *(1000, 0, 100)), "dr_mww"))
        hists.append(df.Histo1D(("w1_nosel", "", *(200, 0, 200)), "w1"))
        hists.append(df.Histo1D(("w2_nosel", "", *(200, 0, 200)), "w2"))

        hists.append(df.Histo1D(("w1_p_nosel", "", *(200, 0, 200)), "w1_p"))
        hists.append(df.Histo1D(("w2_p_nosel", "", *(200, 0, 200)), "w2_p"))

        ##############################################################
        ## START SELECTION
        ##############################################################



        ## CUT ON m(qq)
        
        hists.append(df.Histo1D(("zqq_m_best_nOne", "", *(200, 0, 200)), "zqq_m_best"))
        ##df = df.Filter("zqq_m_best > 60 && zqq_m_best < 110") # tighter at 80, but aa/mumu problem
        #df = df.Filter("zqq_m_best > 60 && zqq_m_best < 110") ## loose
        #df = df.Filter("zqq_m_best > 40 && zqq_m_best < 120") ## loose
        df = df.Filter("zqq_m_best > 40 && zqq_m_best < 140") ## loose
        #df = df.Filter("zqq_m_best > 80 && zqq_m_best < 110") ## aggressive 80-60 seems not too much of a difference --> because we fit it!
        hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut3"))


        ## CUT ON p(qq)
        hists.append(df.Histo1D(("zqq_p_best_nOne", "", *(200, 0, 200)), "zqq_p_best"))
        ##df = df.Filter("zqq_p_best < 60 && zqq_p_best > 30")
        df = df.Filter("zqq_p_best < 60") # see significance
        hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut4"))


        ## cut on cos(qq)
        hists.append(df.Histo1D(("z_costheta_nOne", "", *bins_cos_abs), "z_costheta"))
        df = df.Filter("z_costheta < 0.95")
        hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut5"))

        ## acolinearity
        df = df.Define("acolinearity", "FCCAnalyses::acolinearity(zqq_jets_best)")
        hists.append(df.Histo1D(("acolinearity_nOne", "", *bins_aco), "acolinearity"))
        df = df.Filter("acolinearity > 0.35")
        hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut6"))


        ## acoplanarity
        df = df.Define("acoplanarity", "FCCAnalyses::acoplanarity(zqq_jets_best)")
        hists.append(df.Histo1D(("acoplanarity_nOne", "", *bins_aco), "acoplanarity"))
        df = df.Filter("acoplanarity < 5")
        hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut7"))




        hists.append(df.Histo2D(("w1_w2_m_nOne", "", *((150, 0, 150)+(150, 0, 150))), "w1", "w2"))

        hists.append(df.Histo1D(("dphi_mww_nOne", "", *(4000, 0, 4)), "dphi_mww"))

        #df =df.Filter("(w1 < 75 || w1 > 80) && (w2 < 75 || w2 > 80)")
        df =df.Filter("((w1-78)*(w1-78) + (w2-78)*(w2-78)) > 5*5")
        hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut8"))


        hists.append(df.Histo1D(("w1_p_nOne", "", *(200, 0, 200)), "w1_p"))
        hists.append(df.Histo1D(("w2_p_nOne", "", *(200, 0, 200)), "w2_p"))


        ####
        ## CUT 7: cos theta(miss)
        ####
        df = df.Define("missingEnergy_rp", "FCCAnalyses::missingEnergy(ecm, ReconstructedParticles)")
        df = df.Define("missingEnergy", "missingEnergy_rp[0].energy")
        df = df.Define("cosThetaMiss", "FCCAnalyses::get_cosTheta_miss(missingEnergy_rp)")
        hists.append(df.Histo1D(("cosThetaMiss_nOne", "", *bins_cosThetaMiss), "cosThetaMiss"))
        
        df = df.Filter("cosThetaMiss < .99")
        hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut9"))

        ## final histograms
        hists.append(df.Histo1D(("zqq_recoil_m", "", *(200, 0, 200)), "zqq_recoil_m_best"))
        hists.append(df.Histo1D(("zqq_m", "", *(200, 0, 200)), "zqq_m_best"))
        #hists.append(df.Histo2D(("zqq_recoil_m_mqq", "", *((50, 100, 150)+(60, 60, 120))), "zqq_recoil_m_best", "zqq_m_best"))
        hists.append(df.Histo3D(("zqq_recoil_m_mqq_pqq", "", *((50, 100, 150)+(100, 40, 140)+(40, 20, 60))), "zqq_recoil_m_best", "zqq_m_best", "zqq_p_best"))
        hists.append(df.Histo2D(("zqq_recoil_m_mqq", "", *((50, 100, 150)+(100, 40, 140))), "zqq_recoil_m_best", "zqq_m_best"))
        hists.append(df.Histo2D(("zqq_recoil_m_pqq", "", *((50, 100, 150)+(40, 20, 60))), "zqq_recoil_m_best", "zqq_p_best"))





        return df

    # define output branches to be saved
    def output():
        branchList = ["zqq_recoil_m_best", "zqq_p_best", "zqq_m_best", "leading_jet_costheta", "subleading_jet_costheta", "leading_jet_p", "subleading_jet_p", "acolinearity", "acoplanarity", "W1_p", "W2_p", "W1_m", "W2_m", "z_costheta"]
        #if doInference:
        #    branchList.append("mva_score")
        return branchList