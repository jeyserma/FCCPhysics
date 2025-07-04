
#include <math.h>
#include "edm4hep/SimTrackerHitData.h"

// definitions here: https://github.com/key4hep/EDM4hep/blob/main/edm4hep.yaml

// generic definitions
using Vec_b = ROOT::VecOps::RVec<bool>;
using Vec_d = ROOT::VecOps::RVec<double>;
using Vec_f = ROOT::VecOps::RVec<float>;
using Vec_i = ROOT::VecOps::RVec<int>;
using Vec_ui = ROOT::VecOps::RVec<unsigned int>;


// detector-specific collections
using Vec_SimTrackerHitData = ROOT::VecOps::RVec<edm4hep::SimTrackerHitData>;



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
        ret.push_back(hit.EDep);
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

