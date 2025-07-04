
#include <math.h>

struct PixelDigi {
    int row;
    int col;
    float adc; // ADC counts or energy
};

struct BarrelGeometry {
    float R; // layer radius, mm
    float z; // z extend, mm (single side)
    float pX = 25; // X pitch in um
    float pY = 25; // Y pitch in um
    float c; // circumference of layer
    BarrelGeometry() {}
    BarrelGeometry(float R_, float z_, float pX_, float pY_) : R(R_), z(z_), pX(pX_), pY(pY_), c(2*M_PI*R_) {}
};

using Vec_PixelDigi = ROOT::VecOps::RVec<PixelDigi>;





Vec_PixelDigi DigitizerSimHitBarrel(const BarrelGeometry& geo, Vec_f hits_x, Vec_f hits_y, Vec_f hits_z, Vec_f hits_e, float thrs) {

    Vec_PixelDigi digis; // vector holding all the digis

    for(int i=0; i<hits_e.size(); i++) {
        float depositEnergy = hits_e[i] * 1000; // MeV
        if(depositEnergy < thrs) {
            continue;
        }

        // local coordinates
        float phi = std::atan2(hits_y[i], hits_x[i]);
        phi = (phi > 0 ? phi : (2*M_PI + phi)); // between 0 and 2PI
        float x_pos = geo.R*phi*1000.; // um
        float y_pos = hits_z[i]*1000.; // um
        float z_pos = 0; // always 0 in local coordinates

        int row = static_cast<int>(x_pos / geo.pX); // in rphi direction (unrolled)
        int col = static_cast<int>(y_pos / geo.pY); // in z direction

        digis.push_back({row, col, depositEnergy});
    }

    return digis;
}


Vec_f getOccupancyBarrel(const BarrelGeometry& geo, Vec_PixelDigi digis, int nR, int nC) {
    // nR and nC are the window sizes in rows() and cols()

    int n_row_pixels_per_window = geo.c * 1000 / nR / geo.pX; // number of pixel sides per row
    int n_col_pixels_per_window = geo.z * 1000 * 2.0 / nC / geo.pY; // number of pixel sides per col
    int n_col_pixels_half_z = geo.z * 1000 / geo.pY; // number of pixels in half z side

    Vec_i hitmap = Vec_i(nR*nC, 0);
    for(const auto& b : digis) {
        int r = b.row/n_row_pixels_per_window;
        int c = (b.col + n_col_pixels_half_z)/n_col_pixels_per_window; // avoid negative cols, add n_col_pixels_half_z, so c>=0
        hitmap[r*nC+c] += 1;
    }

    int maxCount = *std::max_element(hitmap.begin(), hitmap.end());
    float avgCount = std::reduce(hitmap.begin(), hitmap.end())*1.0 / hitmap.size();

    Vec_f ret;
    ret.push_back(maxCount);
    ret.push_back(avgCount);
    return ret;
}

