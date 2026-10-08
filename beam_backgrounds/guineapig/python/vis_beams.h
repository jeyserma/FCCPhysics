
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
using Vec_tlv = ROOT::VecOps::RVec<TLorentzVector>;
using Vec_mc = ROOT::VecOps::RVec<edm4hep::MCParticleData>;



Vec_tlv makeLorentzVectorsMomentum(Vec_mc in) {
    Vec_tlv result;
    for(auto & p: in) {
        TLorentzVector tlv;
        tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
        result.push_back(tlv);
    }
    return result;
}

Vec_tlv makeLorentzVectorsVertex(Vec_mc in) {
    Vec_tlv result;
    for(auto & p: in) {
        TLorentzVector tlv;
        tlv.SetXYZT(p.vertex.x, p.vertex.y, p.vertex.z, p.time);
        result.push_back(tlv);
    }
    return result;
}

Vec_i getSlice(Vec_mc in) {
    Vec_i result;
    for(auto & p: in) {
        result.push_back(p.time);
    }
    return result;
}

Vec_f getPos(Vec_mc in, int ax) {
    Vec_f result;
    for(auto & p: in) {
        float k;
        if(ax == 0) k = p.vertex.x;
        else if(ax == 1) k = p.vertex.y;
        else k = p.vertex.z;
        result.push_back(k);
        //if(p.time == 1 and ax==0) std::cout << p.vertex.x << " " << p.vertex.y << " " << p.vertex.z << std::endl;
    }
    return result;
}

Vec_f getMom(Vec_mc in, int ax) {
    Vec_f result;
    for(auto & p: in) {
        float k;
        if(ax == 0) k = p.momentum.x;
        else if(ax == 1) k = p.momentum.y;
        else k = p.momentum.z;
        result.push_back(k);
        //if(p.time == 1 and ax==0) std::cout << p.vertex.x << " " << p.vertex.y << " " << p.vertex.z << std::endl;
    }
    return result;
}

bool print(Vec_f in) {
    for(auto & p: in) {
        std::cout << p <<  std::endl;
    }
    return true;
}


TH2D* SliceXY(const TH3* h3, int zb, const char* name, const char* title="") {
    auto x = h3->GetXaxis();
    auto y = h3->GetYaxis();
    auto h2 = new TH2D(name, title, x->GetNbins(), x->GetXmin(), x->GetXmax(),  y->GetNbins(), y->GetXmin(), y->GetXmax());
    h2->SetDirectory(nullptr);
    for(int ix=1; ix<=x->GetNbins(); ++ix) {
        for(int iy=1; iy<=y->GetNbins(); ++iy) {
            h2->SetBinContent(ix, iy, h3->GetBinContent(ix, iy, zb));
        }
    }
    return h2;
}