#include <random>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <algorithm>

namespace FCCAnalyses {

float frac_mz_240 = 1.0;
float frac_pz_240 = 1.0;
float frac_rec_240 = 1.0;
float pz_240 = 52;

float frac_mz_365 = 5.0;
float frac_pz_365 = 1.0;
float frac_rec_365 = 1.0;
float pz_365 = 143;


Vec_tlv pair_WW_N4(Vec_rp in) {
    // assume 4 input jets
    Vec_tlv ret;

    TLorentzVector j1, j2, j3, j4, W1, W2;
    j1.SetXYZM(in[0].momentum.x, in[0].momentum.y, in[0].momentum.z, in[0].mass);
    j2.SetXYZM(in[1].momentum.x, in[1].momentum.y, in[1].momentum.z, in[1].mass);
    j3.SetXYZM(in[2].momentum.x, in[2].momentum.y, in[2].momentum.z, in[2].mass);
    j4.SetXYZM(in[3].momentum.x, in[3].momentum.y, in[3].momentum.z, in[3].mass);

    float chi2_1 = std::pow((j1+j2).M()-80.0, 2) + std::pow((j3+j4).M()-80.0, 2);
    float chi2_2 = std::pow((j1+j3).M()-80.0, 2) + std::pow((j2+j4).M()-80.0, 2);
    float chi2_3 = std::pow((j1+j4).M()-80.0, 2) + std::pow((j2+j3).M()-80.0, 2);

    if(chi2_1<chi2_2 && chi2_1<chi2_3) {
        W1 = j1+j2;
        W2 = j3+j4;
    }
    else if(chi2_2<chi2_1 && chi2_2<chi2_3) {
        W1 = j1+j3;
        W2 = j2+j4;
    }
    else if(chi2_3<chi2_1 && chi2_3<chi2_2) {
        W1 = j1+j4;
        W2 = j2+j3;
    }

    ret.push_back(W1);
    ret.push_back(W2);
    return ret;
}


Vec_f pair_W_p(Vec_rp in) {
    // assume 4 input jets
    Vec_f ret;

    TLorentzVector j1, j2, j3, j4;
    j1.SetXYZM(in[0].momentum.x, in[0].momentum.y, in[0].momentum.z, in[0].mass);
    j2.SetXYZM(in[1].momentum.x, in[1].momentum.y, in[1].momentum.z, in[1].mass);
    j3.SetXYZM(in[2].momentum.x, in[2].momentum.y, in[2].momentum.z, in[2].mass);
    j4.SetXYZM(in[3].momentum.x, in[3].momentum.y, in[3].momentum.z, in[3].mass);

    float chi2_1 = std::pow((j1+j2).M()-80.0, 2) + std::pow((j3+j4).M()-80.0, 2);
    float chi2_2 = std::pow((j1+j3).M()-80.0, 2) + std::pow((j2+j4).M()-80.0, 2);
    float chi2_3 = std::pow((j1+j4).M()-80.0, 2) + std::pow((j2+j3).M()-80.0, 2);
    
    if(chi2_1<chi2_2 && chi2_1<chi2_3) {
        float w1 = (j1+j2).P();
        float w2 = (j3+j4).P();
        ret.push_back(w1);
        ret.push_back(w2);
    }
    else if(chi2_2<chi2_1 && chi2_2<chi2_3) {
        float w1 = (j1+j3).P();
        float w2 = (j2+j4).P();
        ret.push_back(w1);
        ret.push_back(w2);
    }
    else if(chi2_3<chi2_1 && chi2_3<chi2_2) {
        float w1 = (j1+j4).P();
        float w2 = (j2+j3).P();
        ret.push_back(w1);
        ret.push_back(w2);
    }
    return ret;
}

float pair_W_dphi(Vec_rp in) {
    TLorentzVector j1, j2, j3, j4;
    j1.SetXYZM(in[0].momentum.x, in[0].momentum.y, in[0].momentum.z, in[0].mass);
    j2.SetXYZM(in[1].momentum.x, in[1].momentum.y, in[1].momentum.z, in[1].mass);
    j3.SetXYZM(in[2].momentum.x, in[2].momentum.y, in[2].momentum.z, in[2].mass);
    j4.SetXYZM(in[3].momentum.x, in[3].momentum.y, in[3].momentum.z, in[3].mass);

    float chi2_1 = std::pow((j1+j2).M()-80.0, 2) + std::pow((j3+j4).M()-80.0, 2);
    float chi2_2 = std::pow((j1+j3).M()-80.0, 2) + std::pow((j2+j4).M()-80.0, 2);
    float chi2_3 = std::pow((j1+j4).M()-80.0, 2) + std::pow((j2+j3).M()-80.0, 2);
    
    float ret = -999;
    if(chi2_1<chi2_2 && chi2_1<chi2_3) {
        ret = (j1+j2).DeltaPhi((j3+j4));
    }
    else if(chi2_2<chi2_1 && chi2_2<chi2_3) {
        ret = (j1+j3).DeltaPhi((j2+j4));
    }
    else if(chi2_3<chi2_1 && chi2_3<chi2_2) {
        float w1 = (j1+j4).M();
        float w2 = (j2+j3).M();
        ret = (j1+j4).DeltaPhi((j2+j3));
    }

    return std::abs(ret);
}

Vec_f pair_W(Vec_rp in) {
    // assume 4 input jets
    Vec_f ret;

    TLorentzVector j1, j2, j3, j4;
    j1.SetXYZM(in[0].momentum.x, in[0].momentum.y, in[0].momentum.z, in[0].mass);
    j2.SetXYZM(in[1].momentum.x, in[1].momentum.y, in[1].momentum.z, in[1].mass);
    j3.SetXYZM(in[2].momentum.x, in[2].momentum.y, in[2].momentum.z, in[2].mass);
    j4.SetXYZM(in[3].momentum.x, in[3].momentum.y, in[3].momentum.z, in[3].mass);

    float chi2_1 = std::pow((j1+j2).M()-80.0, 2) + std::pow((j3+j4).M()-80.0, 2);
    float chi2_2 = std::pow((j1+j3).M()-80.0, 2) + std::pow((j2+j4).M()-80.0, 2);
    float chi2_3 = std::pow((j1+j4).M()-80.0, 2) + std::pow((j2+j3).M()-80.0, 2);
    
    if(chi2_1<chi2_2 && chi2_1<chi2_3) {
        float w1 = (j1+j2).M();
        float w2 = (j3+j4).M();
        ret.push_back(w1);
        ret.push_back(w2);
    }
    else if(chi2_2<chi2_1 && chi2_2<chi2_3) {
        float w1 = (j1+j3).M();
        float w2 = (j2+j4).M();
        ret.push_back(w1);
        ret.push_back(w2);
    }
    else if(chi2_3<chi2_1 && chi2_3<chi2_2) {
        float w1 = (j1+j4).M();
        float w2 = (j2+j3).M();
        ret.push_back(w1);
        ret.push_back(w2);
    }

    /*
    float m12 = std::abs((j1+j2).M()-80.0);
    float m13 = std::abs((j1+j3).M()-80.0);
    float m14 = std::abs((j1+j4).M()-80.0);

    float m23 = std::abs((j2+j3).M()-80.0);
    float m24 = std::abs((j2+j4).M()-80.0);

    float m34 = std::abs((j3+j4).M()-80.0);
    

    if(m12<m13 && m12<m14 && m12<m23 && m12<m24 && m12<m34) {
        float w1 = (j1+j2).M();
        float w2 = (j3+j4).M();
        ret.push_back(w1);
        ret.push_back(w2);
    }
    else if(m13<m12 && m13<m14 && m13<m23 && m13<m24 && m13<m34) {
        float w1 = (j1+j3).M();
        float w2 = (j2+j4).M();
        ret.push_back(w1);
        ret.push_back(w2);
    }
    else if(m14<m13 && m14<m12 && m14<m23 && m14<m24 && m14<m34) {
        float w1 = (j1+j4).M();
        float w2 = (j2+j3).M();
        ret.push_back(w1);
        ret.push_back(w2);
    }
    else if(m23<m13 && m23<m14 && m23<m12 && m23<m24 && m23<m34) {
        float w1 = (j2+j3).M();
        float w2 = (j1+j4).M();
        ret.push_back(w1);
        ret.push_back(w2);
    }
    else if(m24<m13 && m24<m14 && m24<m23 && m24<m12 && m24<m34) {
        float w1 = (j2+j4).M();
        float w2 = (j1+j3).M();
        ret.push_back(w1);
        ret.push_back(w2);
    }
    else if(m34<m13 && m34<m14 && m34<m23 && m34<m24 && m34<m12) {
        float w1 = (j3+j4).M();
        float w2 = (j1+j2).M();
        ret.push_back(w1);
        ret.push_back(w2);
    }
    */
    return ret;
}

Vec_rp jets2rp(Vec_f px, Vec_f py, Vec_f pz, Vec_f e, Vec_f m) {
    // checked std::sqrt(px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i] + m[i]*m[i]) equals the energy
    Vec_rp ret;
    for(int i = 0; i < px.size(); i++) {
        edm4hep::ReconstructedParticleData p;
        p.momentum.x = px[i];
        p.momentum.y = py[i];
        p.momentum.z = pz[i];
        p.mass = m[i];
        p.energy = e[i];
        p.charge = 0;
        ret.push_back(p);
    }
    return ret;
}

Vec_rp jets2rp_smeared(
    const Vec_f &px,
    const Vec_f &py,
    const Vec_f &pz,
    const Vec_f &e,
    const Vec_f &m,
    double jes_factor,
    double jer_base,
    double jer_sf,
    std::uint64_t event_id,
    std::uint64_t seed = 12345
) {
    /*
      Jet energy scale correction + jet energy resolution smearing.

      jes_factor: multiplicative scale correction
                  1.00 = nominal
                  1.01 = JES +1%
                  0.99 = JES -1%

      jer_base:   existing relative jet energy resolution
                  0.04 = 4%

      jer_sf:     target resolution / existing resolution
                  1.00 = no additional smearing
                  1.10 = 10% worse resolution
                  1.20 = 20% worse resolution

      event_id:   unique event identifier, e.g. rdfentry_

      seed:       fixed random-seed offset

      Each jet receives an independent Gaussian fluctuation.
      The complete four-momentum is scaled by the same factor.
    */

    if (jes_factor <= 0.0 || jer_base < 0.0 || jer_sf < 1.0) {
        throw std::invalid_argument(
            "Invalid JES or JER parameters"
        );
    }

    if (py.size() != px.size() ||
        pz.size() != px.size() ||
        e.size()  != px.size() ||
        m.size()  != px.size()) {
        throw std::invalid_argument(
            "Jet component vectors have different sizes"
        );
    }

    Vec_rp ret;
    ret.reserve(px.size());

    // Additional relative resolution.
    // Variances are added in quadrature.
    const double sigma_add =
        jer_base * std::sqrt(jer_sf * jer_sf - 1.0);

    for (std::size_t i = 0; i < px.size(); ++i) {

        edm4hep::ReconstructedParticleData p{};

        // Empty/padded jets should remain empty.
        if (e[i] <= 0.0f) {
            ret.push_back(p);
            continue;
        }

        double z = 0.0;

        if (sigma_add > 0.0) {

            // Deterministic seed for this event and jet.
            // No shared RNG: safe for RDataFrame multithreading.
            std::uint64_t key =
                event_id ^
                (seed + 0x9e3779b97f4a7c15ULL *
                        (static_cast<std::uint64_t>(i) + 1));

            // SplitMix64 finalizer to decorrelate nearby seeds.
            key ^= key >> 30;
            key *= 0xbf58476d1ce4e5b9ULL;
            key ^= key >> 27;
            key *= 0x94d049bb133111ebULL;
            key ^= key >> 31;

            std::mt19937_64 rng(key);
            std::normal_distribution<double> gauss(0.0, 1.0);

            z = gauss(rng);
        }

        // Scale correction followed by additional smearing.
        // The lower bound protects against an unphysical
        // negative energy in an extreme Gaussian fluctuation.
        const double smear_factor =
            std::max(0.01, 1.0 + z * sigma_add);

        const double factor =
            jes_factor * smear_factor;

        // Apply the same factor to the entire four-vector.
        p.momentum.x = px[i] * factor;
        p.momentum.y = py[i] * factor;
        p.momentum.z = pz[i] * factor;

        p.energy = e[i] * factor;
        p.mass   = m[i] * factor;

        p.charge = 0;

        ret.push_back(p);
    }

    return ret;
}

Vec_rp select_jets(Vec_rp in, std::vector<std::vector<int>> constituents, int njets_sel, Vec_rp reco) {
    // njets_sel = the current njets clustering algo
    Vec_rp ret;
    for(int i = 0; i < in.size(); i++) {
        float p = std::sqrt(in[i].momentum.x*in[i].momentum.x + in[i].momentum.y*in[i].momentum.y + in[i].momentum.z*in[i].momentum.z);
        //if(p < 0) continue; // at least 10 GeV
        //if(njets_sel == 0 and p < 10) continue;
        if(p < 5) continue; // at least 5 GeV momemtum
        //if(constituents[i].size() < 5) continue; // has impact on m(qq) for H(aa) H(mumu)??
        // might need to reject isolated muons/electrons/photons
        /*
        bool hasTarget = false;
        if(constituents[i].size() < 5) {
            
            for(int k=0; k<constituents[i].size();  k++) {
                int pdgid = std::abs(reco[k].type);
                if((pdgid == 13 || pdgid == 11 || pdgid == 22)) { // 
                    hasTarget = true;
                    cout << "REMOVE FROM JET " << njets_sel << " "  << pdgid << " " << p << endl;
                    break;
                }
            }
        }
        */
        //if(!hasTarget) ret.push_back(in[i]);
        ret.push_back(in[i]);
    }
    return ret;
}

Vec_rp expand_jets(Vec_rp in, int njets) {
    if(in.size() < njets) in.resize(njets); // add empty jets
    return in;
}

int best_clustering_idx(Vec_f mz, Vec_f pz, Vec_f mrec, Vec_i njets, Vec_i njets_target, int ecm=240) {

    float frac_mz = 1.0;
    float frac_pz = 1.0;
    float frac_rec = 1.0;
    float pz_ = 52;
    if(ecm == 240) {
        frac_mz = frac_mz_240;
        frac_pz = frac_pz_240;
        frac_rec = frac_rec_240;
        pz_ = pz_240;
    }
    if(ecm == 365) {
        frac_mz = frac_mz_365;
        frac_pz = frac_pz_365;
        frac_rec = frac_rec_365;
        pz_ = pz_365;
    }

    float mz_ = 91.2;
    float mh_ = 125.0;

    Vec_f chi2;
    for(int i = 0; i < mz.size(); i++) {

        float c = frac_mz*std::pow(mz[i] - mz_, 2) + frac_rec*std::pow(mrec[i] - mh_, 2) + frac_pz*std::pow(pz[i] - pz_, 2);
        if(i==0 && njets[0] < 2) c = 9e99;
        else if(i > 0 && njets[i] != njets_target[i]) c = 9e99;
        chi2.push_back(c);
    }

    // extra constraints: the number of good candidate jets must be at least 2 to form a good z-candidate
    // sometimes exclusive clustering gives wrong njets
    int min_dx = std::distance(std::begin(chi2), std::min_element(std::begin(chi2), std::end(chi2)));
    if(min_dx == 0 and njets[0] >= 2) return 0; // requirement for inclusive clustering
    else if(njets[min_dx] == njets_target[min_dx] && min_dx != 0) return min_dx; // requirement for exlusive clustering
    else {
        //std::cout << "CANNOT FIND GOOD CLUSTERING" <<  std::endl;
        return -1; // could not cluster
    }

    /*
    float mz_ = 91.2;
    float mh_ = 125.0;

    Vec_f chi2;
    for(int i = 0; i < mz.size(); i++) {

        float c = std::pow(mz[i] - mz_, 2) + std::pow(mrec[i] - mh_, 2);
        chi2.push_back(c);
    }

    // extra constraints: the number of good candidate jets must be at least 2 to form a good z-candidate
    // sometimes exclusive clustering gives wrong njets
    int min_dx = std::distance(std::begin(chi2), std::min_element(std::begin(chi2), std::end(chi2)));
    if(min_dx == 0 and njets[0] >= 2) return 0; // requirement for inclusive clustering
    else if(njets[min_dx] == njets_target[min_dx] && min_dx != 0) return min_dx; // requirement for exlusive clustering
    else {
        std::cout << "CANNOT FIND GOOD CLUSTERING" <<  std::endl;
        return -1; // could not cluster
    }
    */
}


// build the Z resonance based on the available leptons. Returns the best lepton pair compatible with the Z mass and recoil at 125 GeV
// technically, it returns a ReconstructedParticleData object with index 0 the di-lepton system, index and 2 the leptons of the pair
struct resonanceBuilder_mass_recoil_hadronic {
    float m_resonance_mass;
    float m_recoil_mass;
    float chi2_recoil_frac;
    float ecm;
    resonanceBuilder_mass_recoil_hadronic(float arg_resonance_mass, float arg_recoil_mass, float arg_chi2_recoil_frac, float arg_ecm);
    Vec_rp operator()(Vec_rp legs);
};

resonanceBuilder_mass_recoil_hadronic::resonanceBuilder_mass_recoil_hadronic(float arg_resonance_mass, float arg_recoil_mass, float arg_chi2_recoil_frac, float arg_ecm) {m_resonance_mass = arg_resonance_mass, m_recoil_mass = arg_recoil_mass, chi2_recoil_frac = arg_chi2_recoil_frac, ecm = arg_ecm;}

Vec_rp resonanceBuilder_mass_recoil_hadronic::resonanceBuilder_mass_recoil_hadronic::operator()(Vec_rp legs) {
    float frac_mz = 1.0;
    float frac_pz = 1.0;
    float frac_rec = 1.0;
    float pz_ = 52;
    if(ecm == 240) {
        frac_mz = frac_mz_240;
        frac_pz = frac_pz_240;
        frac_rec = frac_rec_240;
        pz_ = pz_240;
    }
    if(ecm == 365) {
        frac_mz = frac_mz_365;
        frac_pz = frac_pz_365;
        frac_rec = frac_rec_365;
        pz_ = pz_365;
    }


    Vec_rp result;
    result.reserve(3);
    std::vector<std::vector<int>> pairs; // for each permutation, add the indices of the muons
    int n = legs.size();
    //cout << "******************* " << n << endl;
    if(n > 1) {
        ROOT::VecOps::RVec<bool> v(n);
        std::fill(v.end() - 2, v.end(), true); // helper variable for permutations
        do {
            std::vector<int> pair;
            rp reso;
            reso.charge = 0;
            TLorentzVector reso_lv;
            Vec_f jet_momenta;
            for(int i = 0; i < n; ++i) {
                if(v[i]) {
                    pair.push_back(i);
                    TLorentzVector leg_lv;

                    leg_lv.SetXYZM(legs[i].momentum.x, legs[i].momentum.y, legs[i].momentum.z, legs[i].mass);
                    reso_lv += leg_lv;
                    jet_momenta.push_back(leg_lv.P());
                }
            }
            
            // condition on jet momenta: leading between 20-90 GeV, subleading less than 60 GeV
            //if(!((jet_momenta[0] > jet_momenta[1]) && (jet_momenta[0] < 90 && jet_momenta[0] > 20))) continue;
            //if(!((jet_momenta[0] > jet_momenta[1]) && (jet_momenta[1] < 60))) continue;

            reso.momentum.x = reso_lv.Px();
            reso.momentum.y = reso_lv.Py();
            reso.momentum.z = reso_lv.Pz();
            reso.mass = reso_lv.M();
            
            
            result.emplace_back(reso);
            pairs.push_back(pair);

        } while(std::next_permutation(v.begin(), v.end()));
    }
    else {
        //std::cout << "ERROR: resonanceBuilder_mass_recoil_hadronic, at least two jets required. Return dummy." << std::endl;
        //exit(1);
        return result;
    }

    if(result.size() > 1) {
        //cout << "*******************" << endl;
        Vec_rp bestReso;
        int idx_min = -1;
        float d_min = 9e9;
        for (int i = 0; i < result.size(); ++i) {
            
            // calculate recoil
            auto recoil_p4 = TLorentzVector(0, 0, 0, ecm);
            TLorentzVector tv1;
            tv1.SetXYZM(result.at(i).momentum.x, result.at(i).momentum.y, result.at(i).momentum.z, result.at(i).mass);
            recoil_p4 -= tv1;

            auto recoil_fcc = edm4hep::ReconstructedParticleData();
            recoil_fcc.momentum.x = recoil_p4.Px();
            recoil_fcc.momentum.y = recoil_p4.Py();
            recoil_fcc.momentum.z = recoil_p4.Pz();
            recoil_fcc.mass = recoil_p4.M();

            TLorentzVector tg;
            tg.SetXYZM(result.at(i).momentum.x, result.at(i).momentum.y, result.at(i).momentum.z, result.at(i).mass);

            float momentum = tg.P();
            float mass = frac_mz*std::pow(result.at(i).mass - m_resonance_mass, 2); // mass
            float rec = frac_rec*std::pow(recoil_fcc.mass - m_recoil_mass, 2); // recoil
            float p = frac_pz*std::pow(momentum - pz_, 2); // recoil
            //float d = (1.0-chi2_recoil_frac)*mass + chi2_recoil_frac*rec;
            float d = mass + rec + p;
            if(d < d_min) {
                d_min = d;
                idx_min = i;
            }
            //cout << ecm << " " << pz_ << " masss=" << result.at(i).mass << " recoil_mass=" << recoil_fcc.mass << " momentum=" << momentum << endl;
            //if(result.at(i).mass < 127 && result.at(i).mass > 123) {
            //    cout << "!!!!! PREVIOUS" << endl;
            //}
        }
        if(idx_min > -1) { 
            bestReso.push_back(result.at(idx_min));
            //cout << "SELECTED " << result.at(idx_min).mass << endl;
            auto & l1 = legs[pairs[idx_min][0]];
            auto & l2 = legs[pairs[idx_min][1]];
            bestReso.emplace_back(l1);
            bestReso.emplace_back(l2);
        }
        else {
            std::cout << "ERROR: resonanceBuilder_mass_recoil, no mininum found." << std::endl;
            exit(1);
        }
        return bestReso;
    }
    else {
        auto & l1 = legs[0];
        auto & l2 = legs[1];
        result.emplace_back(l1);
        result.emplace_back(l2);
        return result;
    }
}


Vec_rp energyReconstructFourJet(float ecm, Vec_f px, Vec_f py, Vec_f pz, Vec_f e) {

    // assume 4 jets as input
    float p0 = std::sqrt(px[0]*px[0] + py[0]*py[0] + pz[0]*pz[0]);
    float p1 = std::sqrt(px[1]*px[1] + py[1]*py[1] + pz[1]*pz[1]);
    float p2 = std::sqrt(px[2]*px[2] + py[2]*py[2] + pz[2]*pz[2]);
    float p3 = std::sqrt(px[3]*px[3] + py[3]*py[3] + pz[3]*pz[3]);

    TMatrixD mtrx(4, 4);
    mtrx(0, 0) = 1;
    mtrx(0, 1) = 1;
    mtrx(0, 2) = 1;
    mtrx(0, 3) = 1;

    mtrx(1, 0) = px[0]/e[0];
    mtrx(1, 1) = px[1]/e[1];
    mtrx(1, 2) = px[2]/e[2];
    mtrx(1, 3) = px[3]/e[3];
    
    mtrx(2, 0) = py[0]/e[0];
    mtrx(2, 1) = py[1]/e[1];
    mtrx(2, 2) = py[2]/e[2];
    mtrx(2, 3) = py[3]/e[3];

    mtrx(3, 0) = pz[0]/e[0];
    mtrx(3, 1) = pz[1]/e[1];
    mtrx(3, 2) = pz[2]/e[2];
    mtrx(3, 3) = pz[3]/e[3];

    TMatrixD inv = mtrx.Invert();

    TVectorD vec(4);
    vec(0) = ecm;
    vec(1) = 0;
    vec(2) = 0;
    vec(3) = 0;

    TVectorD res = inv*vec;

    bool isValid = false;
    
    if(res[0]<0 or res[1]<0 or res[2]<0 or res[3]<0 or res[0]>240 or res[1]>240 or res[2]>240 or res[3]>240) {
        isValid = false;
    }
    if(!isValid) {
        /*
        cout << "***************" << endl;
        cout << px[0] << " " << py[0] << " " << pz[0] << " " << p0 << " " << e[0] << " " << res[0] << endl;
        cout << px[1] << " " << py[1] << " " << pz[1] << " " << p1 << " " << e[1] << " " << res[1] << endl;
        cout << px[2] << " " << py[2] << " " << pz[2] << " " << p2 << " " << e[2] << " " << res[2] << endl;
        cout << px[3] << " " << py[3] << " " << pz[3] << " " << p3 << " " << e[3] << " " << res[3] << endl;
        */
    }

    //Vec_tlv ret;
    Vec_rp ret;
    float chi2 = 0;
    for(int i=0; i<4; i++) {
        
        auto rp = edm4hep::ReconstructedParticleData();
        
        
        TLorentzVector tlv;
        if(isValid) {
            tlv.SetPxPyPzE(px[i]*res[i]/e[i], py[i]*res[i]/e[i], pz[i]*res[i]/e[i], res[i]);
            //rp.momentum.x = px[i]*res[i]/e[i];
            //rp.momentum.y = py[i]*res[i]/e[i];
            //rp.momentum.z = pz[i]*res[i]/e[i];
            //rp.energy = res[i];
        }
        else {
            tlv.SetPxPyPzE(px[i], py[i], pz[i], e[i]);
            //rp.momentum.x = px[i];
            //rp.momentum.y = py[i];
            //rp.momentum.z = pz[i];
            //rp.energy = e[i];
        }
        //float m = std::sqrt(rp.energy*rp.energy - rp.momentum.x*rp.momentum.x - rp.momentum.y*rp.momentum.y - rp.momentum.z*rp.momentum.z);
        //rp.mass = m;
        
        rp.momentum.x = tlv.Px();
        rp.momentum.y = tlv.Py();
        rp.momentum.z = tlv.Pz();
        rp.mass = tlv.M();
        rp.energy = tlv.E();
        
        ret.push_back(rp);

        /*
        if(res[i] > 0) {
            float uncert = 0.5*std::sqrt(e[i]) + 0.05*e[i];
            float delta = (e[i]-res[i])/uncert;
            chi2 += delta*delta;
        }
        else {
            chi2 += 1000.;
        }
        */
    }
    
    // add chi2 as dummy to the list of Lorentz vectors
    //TLorentzVector chi2_;
    //chi2_.SetPxPyPzE(0, 0, 0, chi2);
    //ret.push_back(chi2_);
    
    
    return ret;
}



Vec_rp energyReconstructFourJet_rp(float ecm, Vec_rp in) {

    // assume 4 jets as input
    if(in[0].energy <= 0 || in[1].energy <= 0 || in[2].energy <= 0 || in[3].energy <= 0) {
        return in;
    }

    TMatrixD mtrx(4, 4);
    mtrx(0, 0) = 1;
    mtrx(0, 1) = 1;
    mtrx(0, 2) = 1;
    mtrx(0, 3) = 1;

    mtrx(1, 0) = in[0].momentum.x/in[0].energy;
    mtrx(1, 1) = in[1].momentum.x/in[1].energy;
    mtrx(1, 2) = in[2].momentum.x/in[2].energy;
    mtrx(1, 3) = in[3].momentum.x/in[3].energy;
    
    mtrx(2, 0) = in[0].momentum.y/in[0].energy;
    mtrx(2, 1) = in[1].momentum.y/in[1].energy;
    mtrx(2, 2) = in[2].momentum.y/in[2].energy;
    mtrx(2, 3) = in[3].momentum.y/in[3].energy;

    mtrx(3, 0) = in[0].momentum.z/in[0].energy;
    mtrx(3, 1) = in[1].momentum.z/in[1].energy;
    mtrx(3, 2) = in[2].momentum.z/in[2].energy;
    mtrx(3, 3) = in[3].momentum.z/in[3].energy;

    TMatrixD inv = mtrx.Invert();

    TVectorD vec(4);
    vec(0) = ecm;
    vec(1) = 0;
    vec(2) = 0;
    vec(3) = 0;

    TVectorD res = inv*vec;

    bool isValid = true;
    
    if(!(std::isfinite(res[0]) && std::isfinite(res[1]) && std::isfinite(res[2]) && std::isfinite(res[3]))) isValid = false;
    if(res[0] <= 0 or res[1] <= 0 or res[2] <= 0 or res[3] <= 0 or res[0] > ecm or res[1] > ecm or res[2] > ecm or res[3] > ecm) {
        isValid = false;
    }
    if(!isValid) {

       // cout << "***************" << endl;
       // cout << px[0] << " " << py[0] << " " << pz[0] << " " << p0 << " " << e[0] << " " << res[0] << endl;
       // cout << px[1] << " " << py[1] << " " << pz[1] << " " << p1 << " " << e[1] << " " << res[1] << endl;
       // cout << px[2] << " " << py[2] << " " << pz[2] << " " << p2 << " " << e[2] << " " << res[2] << endl;
       // cout << px[3] << " " << py[3] << " " << pz[3] << " " << p3 << " " << e[3] << " " << res[3] << endl;
 
    }

    Vec_rp ret;
    for(int i=0; i<4; i++) {
        TLorentzVector tlv;
        if(isValid) {
            //cout << in[i].energy << " " << res[i] << endl;
            tlv.SetPxPyPzE(in[i].momentum.x*res[i]/in[i].energy, in[i].momentum.y*res[i]/in[i].energy, in[i].momentum.z*res[i]/in[i].energy, res[i]);
            //tlv.SetPxPyPzE(in[i].momentum.x, in[i].momentum.y, in[i].momentum.z, in[i].energy);
        }
        else {
            tlv.SetPxPyPzE(in[i].momentum.x, in[i].momentum.y, in[i].momentum.z, in[i].energy);
            //tlv.SetXYZM(in[i].momentum.x, in[i].momentum.y, in[i].momentum.z, in[i].mass);
        }

        auto rp = edm4hep::ReconstructedParticleData();
        rp.momentum.x = tlv.Px();
        rp.momentum.y = tlv.Py();
        rp.momentum.z = tlv.Pz();
        rp.mass = tlv.M();
        rp.energy = tlv.E();
        ret.push_back(rp);
    }

    return ret;
}











// Function that runs the fit for the thrust axis determination
struct thrustFit_mc {
    public:
        thrustFit_mc(const ROOT::VecOps::RVec<float> & arg_px,
        const ROOT::VecOps::RVec<float> & arg_py,
        const ROOT::VecOps::RVec<float> & arg_pz);
    float operator()(const double *par);

    private:
        ROOT::VecOps::RVec<float> _px; // vector of px
        ROOT::VecOps::RVec<float> _py; // vector of py
        ROOT::VecOps::RVec<float> _pz; // vector of pz
};

thrustFit_mc::thrustFit_mc(const ROOT::VecOps::RVec<float> & arg_px, const ROOT::VecOps::RVec<float> & arg_py, const ROOT::VecOps::RVec<float> & arg_pz) {
    _px=arg_px;
    _py=arg_py;
    _pz=arg_pz;
}

float thrustFit_mc::operator()(const double *pars){
    double num = 0.;
    double den = 0.;
    double mag = sqrt(pars[0]*pars[0] + pars[1]*pars[1] + pars[2]*pars[2]);

    for(unsigned int i =0; i<_px.size(); i++) {
        num += std::abs(_px[i]*(pars[0]/mag) + _py[i]*(pars[1]/mag) + _pz[i]*(pars[2]/mag));
        den += sqrt(_px[i]*_px[i] + _py[i]*_py[i] + _pz[i]*_pz[i]);
    }
    if(den>0.){
        double val = num / den;
        return -val;
    }
    return 0.;
};

// Finds the thrust axis based on a list of px, py, pz
// MC based, i.e. 
struct minimize_thrust_mc {
    minimize_thrust_mc(float ipx=1.0, float ipy=1.0, float ipz=1.0, std::string arg_minname="Minuit2", std::string arg_algoname="Migrad", int arg_maxcalls=10000, float arg_tolerance=0.001);
    ROOT::VecOps::RVec<float> operator()(const ROOT::VecOps::RVec<float> & px, const ROOT::VecOps::RVec<float> & py, const ROOT::VecOps::RVec<float> & pz);

    char const *_minname; // Minimizer to use, Minuit2 default
    char const *_algoname; // Optimisation algorithm, Migrad default
    int _maxcalls; // Maximum call to minimization function, default=100000
    float _tolerance; //Tolerance for minimization, default=0.001
    ROOT::Math::Minimizer *_min; //internal ROOT minimizer
    double _step[3]={0.001,0.001,0.001};
    double _variable[3];
};



minimize_thrust_mc::minimize_thrust_mc(float ipx=1.0, float ipy=1.0, float ipz=1.0, std::string arg_minname, std::string arg_algoname, int arg_maxcalls, float arg_tolerance) {
    _minname=arg_minname.c_str();
    _algoname=arg_algoname.c_str();
    _maxcalls=arg_maxcalls;
    _tolerance=arg_tolerance;

    _min = ROOT::Math::Factory::CreateMinimizer(_minname, _algoname);
    _min->SetMaxFunctionCalls(_maxcalls); // for Minuit/Minuit2
    _min->SetMaxIterations(10000);  // for GSL
    _min->SetTolerance(_tolerance);
    _min->SetPrintLevel(0);

    //std::cout << "INIT: f(" << ipx << "," << ipy << "," << ipz << "): " << std::endl;
    _variable[0] = ipx;
    _variable[1] = ipy;
    _variable[2] = ipz;
}

ROOT::VecOps::RVec<float> minimize_thrust_mc::operator()(const ROOT::VecOps::RVec<float> & px, const ROOT::VecOps::RVec<float> & py, const ROOT::VecOps::RVec<float> & pz) {
    _min->SetVariable(0,"x",_variable[0], _step[0]);
    _min->SetVariable(1,"y",_variable[1], _step[1]);
    _min->SetVariable(2,"z",_variable[2], _step[2]);
    // create functon wrapper for minmizer
    // a IMultiGenFunction type
    ROOT::Math::Functor f(thrustFit_mc(px,py,pz),3);
    _min->SetFunction(f);

    //min->SetValidError(true);
    //min->ProvidesError();
    _min->Minimize();
    //std::cout << "is valid error before hesse " << min->IsValidError() <<std::endl;
    //min->Hesse();
    //std::cout << "is valid error after hesse  " << min->IsValidError() <<std::endl;
    //std::cout << "Ncalls  " << _min->NCalls() << "  Niter " << _min->NIterations() <<std::endl;
    //_min->PrintResults();
    const double *xs = _min->X();
    const double *xs_err = _min->Errors();

    //std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "," << xs[2] << "): " << _min->MinValue()  << std::endl;

    ROOT::VecOps::RVec<float> result;
    result.push_back(-1.*_min->MinValue());
    result.push_back(xs[0]);
    result.push_back(xs_err[0]);
    result.push_back(xs[1]);
    result.push_back(xs_err[1]);
    result.push_back(xs[2]);
    result.push_back(xs_err[2]);

    return result;
}







}