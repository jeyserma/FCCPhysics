
#include <math.h>
#include "edm4hep/SimTrackerHitData.h"
#include "edm4hep/MCParticleData.h"

// definitions here: https://github.com/key4hep/EDM4hep/blob/main/edm4hep.yaml

// generic definitions
using Vec_b = ROOT::VecOps::RVec<bool>;
using Vec_d = ROOT::VecOps::RVec<double>;
using Vec_f = ROOT::VecOps::RVec<float>;
using Vec_i = ROOT::VecOps::RVec<int>;
using Vec_ui = ROOT::VecOps::RVec<unsigned int>;


// detector-specific collections
using Vec_SimTrackerHitData = ROOT::VecOps::RVec<edm4hep::SimTrackerHitData>;
using Vec_MCParticleData = ROOT::VecOps::RVec<edm4hep::MCParticleData>;

Vec_MCParticleData getMCParticle(Vec_SimTrackerHitData in, Vec_MCParticleData mc, Vec_i idx) {
    Vec_MCParticleData ret;
    for(int i=0; i<in.size(); i++) {
        ret.push_back(mc[idx[i]]);
    }
    return ret;
}



Vec_f getSimHitPosition_x(Vec_SimTrackerHitData in) {
    Vec_f ret;
    for(auto & hit: in) {
        ret.push_back(hit.position[0]);
    }
    return ret;
}

Vec_f getSimHitPosition_y(Vec_SimTrackerHitData in) {
    Vec_f ret;
    for(auto & hit: in) {
        ret.push_back(hit.position[1]);
    }
    return ret;
}

Vec_f getSimHitPosition_z(Vec_SimTrackerHitData in) {
    Vec_f ret;
    for(auto & hit: in) {
        ret.push_back(hit.position[2]);
    }
    return ret;
}

Vec_f getSimHitPosition_r(Vec_SimTrackerHitData in) {
    Vec_f ret;
    for(auto & hit: in) {
        ret.push_back(std::sqrt(hit.position[0]*hit.position[0] + hit.position[1]*hit.position[1]));
    }
    return ret;
}

Vec_f getSimHitPosition_theta(Vec_SimTrackerHitData in, bool deg=false) {
    Vec_f ret;
    for(auto & hit: in) {
        float x = hit.position[0];
        float y = hit.position[1];
        float z = hit.position[2];
        float theta = std::acos(z/std::sqrt(x*x + y*y + z*z));
        if(deg) theta = theta * 180 / M_PI;
        ret.push_back(theta);
    }
    return ret;
}

Vec_f getSimHitPosition_phi(Vec_SimTrackerHitData in, bool deg=false) {
    Vec_f ret;
    for(auto & hit: in) {
        float x = hit.position[0];
        float y = hit.position[1];
        float phi = std::atan2(y, x);
        if(deg) phi = phi * 180 / M_PI;
        ret.push_back(phi);
    }
    return ret;
}

Vec_i getSimHitLayer(Vec_f hits_r, Vec_f radii) {
    Vec_i ret;
    for(auto & hit_r: hits_r) {
        int layer = -1;
        for(int i=0; i<radii.size(); i++) {
            if(hit_r > 13.0000 && hit_r < 14.2850) layer = 0;
            //if(abs(hit_r-radii[i]) < 4) layer = i;
        }
        ret.push_back(layer);
    }
    return ret;
}

Vec_b isProducedBySecondary(Vec_SimTrackerHitData in) {
    Vec_b ret;
    for(auto & hit: in) {
        // see https://github.com/key4hep/EDM4hep/blob/v00-10-05/edm4hep.yaml#L249
        ret.push_back(hit.quality & (1 << 30));
    }
    return ret;
}

Vec_f getEnergyDeposition(Vec_SimTrackerHitData in) {
    Vec_f ret;
    for(auto & hit: in) {
        ret.push_back(hit.eDep);
    }
    return ret;
}

Vec_i getCellID(Vec_SimTrackerHitData in) {
    Vec_i ret;
    for(auto & hit: in) {
        ret.push_back(hit.cellID);
    }
    return ret;
}


Vec_f get_theta(Vec_MCParticleData in, bool deg=false) {
    Vec_f ret;
    for(auto & p: in) {
        TLorentzVector tlv;
        tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
        float theta = tlv.Theta();
        if(deg) theta = theta * 180 / M_PI;
        ret.push_back(theta);
    }
    return ret;
}


Vec_f get_p(Vec_MCParticleData in) {
    Vec_f ret;
    for(auto & p: in) {
        TLorentzVector tlv;
        tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
        ret.push_back(tlv.P());
    }
    return ret;
}

Vec_f get_pt(Vec_MCParticleData in) {
    Vec_f ret;
    for(auto & p: in) {
        TLorentzVector tlv;
        tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
        ret.push_back(tlv.Pt());
    }
    return ret;
}

Vec_i get_pdgid(Vec_MCParticleData in, bool abs=false) {
    Vec_i ret;
    for(auto & p: in) {
        if(abs) ret.push_back(std::abs(p.PDG));
        else ret.push_back(p.PDG);
    }
    return ret;
}

Vec_i get_generatorStatus(Vec_MCParticleData in) {
    Vec_i ret;
    for(auto & p: in) {
        ret.push_back(p.generatorStatus);
    }
    return ret;
}

Vec_f tf_theta(Vec_f in) {
    Vec_f ret;
    for(auto & t: in) {
        if(t > M_PI/2.) ret.push_back(M_PI - t);
        else ret.push_back(t);
    }
    return ret;
}


