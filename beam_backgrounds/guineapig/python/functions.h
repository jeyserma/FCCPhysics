
#include <math.h>
#include "edm4hep/MCParticleData.h"
#include "ROOT/RVec.hxx"

#include <vector>
#include <cmath>
#include <algorithm>

// definitions here: https://github.com/key4hep/EDM4hep/blob/main/edm4hep.yaml

// generic definitions
using Vec_b = ROOT::VecOps::RVec<bool>;
using Vec_d = ROOT::VecOps::RVec<double>;
using Vec_f = ROOT::VecOps::RVec<float>;
using Vec_i = ROOT::VecOps::RVec<int>;
using Vec_ui = ROOT::VecOps::RVec<unsigned int>;
using Vec_tlv = ROOT::VecOps::RVec<TLorentzVector>;
using Vec_mc = ROOT::VecOps::RVec<edm4hep::MCParticleData>;

using P4 = ROOT::Math::PxPyPzMVector;
using Vec_p4 = ROOT::VecOps::RVec<P4>;

Vec_i select_pairs(const Vec_mc& in, int ptype) {

    Vec_i result;
    for (const auto& p : in) {
        if(p.generatorStatus == ptype) result.push_back(1);
        else result.push_back(0);
    }
    return result;
}

Vec_f get_t(const Vec_mc& in) {

    Vec_f result;
    for (const auto& p : in) {
        result.push_back(p.time);
    }
    return result;
}

Vec_p4 makeP4Vector(const Vec_mc& in, double xing = 0.0) {
    const double beta_boost = std::sin(xing); // xing in radians
    ROOT::Math::BoostX boost(beta_boost);

    Vec_p4 result;
    for (const auto& p : in) {
        P4 p4(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
        result.push_back(boost(p4));
    }
    return result;
}

Vec_i get_q(Vec_mc in) {
    Vec_i ret;
    for(auto & p: in) {
        ret.push_back(p.charge);
    }
    return ret;
}

Vec_f get_theta(Vec_p4 in, bool deg=false, bool tf=false) {
    Vec_f ret;
    for(auto & p: in) {
        float theta = p.Theta();
        if(deg) theta = theta * 180 / M_PI;
        if(!deg && tf && theta > M_PI/2.) theta = M_PI - theta; // transform to 0..pi/2 (assume symmetric)
        ret.push_back(theta);
    }
    return ret;
}

Vec_f get_costheta(Vec_p4 in, bool tf=false) {
    Vec_f ret;
    for(auto & p: in) {
        float theta = p.Theta();
        if(tf && theta > M_PI/2.) theta = M_PI - theta; // transform to 0..pi/2 (assume symmetric)
        ret.push_back(std::cos(theta));
    }
    return ret;
}


Vec_f get_phi(Vec_p4 in, bool deg=false) {
    Vec_f ret;
    for(auto & p: in) {
        float phi = p.Phi();
        if(deg) phi = phi * 180 / M_PI;
        ret.push_back(phi);
    }
    return ret;
}



Vec_f get_E(Vec_p4 in) {
    Vec_f ret;
    for(auto & p: in) {
        ret.push_back(p.E());
    }
    return ret;
}

Vec_f get_p(Vec_p4 in) {
    Vec_f ret;
    for(auto & p: in) {
        ret.push_back(p.P());
    }
    return ret;
}

Vec_f get_px(Vec_p4 in) {
    Vec_f ret;
    for(auto & p: in) {
        ret.push_back(p.Px());
    }
    return ret;
}

Vec_f get_py(Vec_p4 in) {
    Vec_f ret;
    for(auto & p: in) {
        ret.push_back(p.Py());
    }
    return ret;
}

Vec_f get_pz(Vec_p4 in) {
    Vec_f ret;
    for(auto & p: in) {
        ret.push_back(p.Pz());
    }
    return ret;
}

Vec_f get_pt(Vec_p4 in) {
    Vec_f ret;
    for(auto & p: in) {
        ret.push_back(p.Pt());
    }
    return ret;
}


Vec_f get_x(const Vec_mc& in) {
    Vec_f result;
    result.reserve(in.size());
    for (const auto& p : in) {
        result.push_back(p.vertex.x);
    }
    return result;
}

Vec_f get_y(const Vec_mc& in) {
    Vec_f result;
    result.reserve(in.size());
    for (const auto& p : in) {
        result.push_back(p.vertex.y);
    }
    return result;
}

Vec_f get_z(const Vec_mc& in) {
    Vec_f result;
    result.reserve(in.size());
    for (const auto& p : in) {
        result.push_back(p.vertex.z);
    }
    return result;
}













double wrap_0_2pi(double a) {
    const double twopi = 2.0 * M_PI;
    a = std::fmod(a, twopi);
    if (a < 0.0) a += twopi;
    return a;
}

bool hitsCylinderOne(
    double x0, double y0, double z0,
    double px, double py, double pz,
    double q, double Bz,
    double R,
    double zmax = -1.0
) {
    const double pT = std::hypot(px, py);

    // Straight-line case: B = 0 or pT = 0
    if (pT == 0.0 || Bz == 0.0) {
        const double a = px*px + py*py;
        if (a == 0.0) return false;

        const double b = 2.0 * (x0*px + y0*py);
        const double c = x0*x0 + y0*y0 - R*R;

        const double disc = b*b - 4.0*a*c;
        if (disc < 0.0) return false;

        const double sqrtDisc = std::sqrt(disc);
        const double roots[2] = {
            (-b - sqrtDisc) / (2.0*a),
            (-b + sqrtDisc) / (2.0*a)
        };

        for (double t : roots) {
            if (t >= 0.0) {
                const double z = z0 + t*pz;
                if (zmax < 0.0 || std::abs(z) <= zmax) return true;
            }
        }
        return false;
    }

    // pT in MeV, B in T, rho in mm
    const double rho = pT / (0.3 * std::abs(q) * std::abs(Bz));

    // Rotation sign
    const double s = (q * Bz > 0.0) ? 1.0 : -1.0;

    // Circle center
    const double xc = x0 + s * rho * py / pT;
    const double yc = y0 - s * rho * px / pT;

    const double d = std::hypot(xc, yc);

    // No intersection between track circle and detector cylinder
    if (d > rho + R) return false;
    if (d < std::abs(rho - R)) return false;
    if (d == 0.0) return false; // degenerate concentric case

    // Circle-circle intersection: detector circle centered at (0,0), radius R
    // track circle centered at (xc,yc), radius rho
    const double a = (R*R - rho*rho + d*d) / (2.0*d);
    double h2 = R*R - a*a;
    if (h2 < 0.0) h2 = 0.0;
    const double h = std::sqrt(h2);

    const double x2 = a * xc / d;
    const double y2 = a * yc / d;

    const double rx = -yc * h / d;
    const double ry =  xc * h / d;

    const double xints[2] = {x2 + rx, x2 - rx};
    const double yints[2] = {y2 + ry, y2 - ry};

    const double phi0 = std::atan2(y0 - yc, x0 - xc);

    // Important: for qB > 0, the motion is clockwise in xy,
    // so the circle angle decreases.
    const double direction = (q * Bz > 0.0) ? -1.0 : +1.0;

    for (int i = 0; i < 2; ++i) {
        const double phii = std::atan2(yints[i] - yc, xints[i] - xc);

        double dphi;
        if (direction > 0.0) {
            dphi = wrap_0_2pi(phii - phi0);  // CCW
        } else {
            dphi = wrap_0_2pi(phi0 - phii);  // CW
        }

        const double zhit = z0 + (pz / pT) * rho * dphi;

        if (zmax < 0.0 || std::abs(zhit) <= zmax) {
            return true;
        }
    }

    return false;
}



Vec_b hitsCylinder(
    Vec_f x,
    Vec_f y,
    Vec_f z,
    Vec_f px,
    Vec_f py,
    Vec_f pz,
    Vec_i q,
    double Bz,
    double R,
    double zmax = -1.0
) {
    const std::size_t n = x.size();

    Vec_b result;
    result.reserve(n);

    for (std::size_t i = 0; i < n; ++i) {
        result.push_back(
            hitsCylinderOne(
                x[i], y[i], z[i],
                px[i], py[i], pz[i],
                q[i], Bz, R, zmax
            )
        );
    }

    return result;
}