#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include "QEGenerator_3N.hh"
#include "constants.hh"
#include "TVector3.h"

using namespace std;

eNCrossSection * myCS;
QEGenerator_3N * myGen;
TRandom3 * myRand;

/*
Generates theta12-theta23 heatmap based only on wavefunction.
Loops over k1,k2,k3 values s.t k1+k2+k3=0 (vector sum).
Calculate theta12, theta23 using cosine law and adds them to a histogram,
with weight according to wf.
*/

// Parameters
const double k_min = 0.35/GeVfm; // p_min 350MeV/c
const double ktot_max = 5/GeVfm; // arbitrary 
int k_bins = 50; // Default value
int theta_bins = 50; // Number of bins for theta

int main(int argc, char **argv) {
    if (argc > 2) {
        k_bins = std::stoi(argv[1]);
        theta_bins = std::stoi(argv[2]);
    }

    myRand = new TRandom3(0);
    ffModel ffMod = kelly;
    csMethod csMeth = cc1;
    myCS = new eNCrossSection(csMeth, ffMod);
    myGen = new QEGenerator_3N(6.0, myCS, 1, myRand);

    // Create 2D bin array initialized to zero
    vector<vector<double>> heatmap(theta_bins, vector<double>(theta_bins, 0.0));

    double dk = (ktot_max - k_min) / k_bins;
    double dtheta = M_PI / theta_bins;
    int counter = 0;
    for (int i = 0; i < k_bins; ++i) {
        double k1 = k_min + i * dk;
        for (int j = 0; j < i; ++j) {
            double k2 = k_min + j * dk;
            for (int m = 0; m < j; ++m) {
                double k3 = k_min + m * dk;
                counter++;

                if (k1 > k2+k3) {continue;}
                // Jacobian and angles
                double c12 = (k3*k3 - k1*k1 - k2*k2) / (2.0 * k1 * k2); // cos( theta12 )
                double c23 = (k1*k1 - k2*k2 - k3*k3) / (2.0 * k2 * k3); // cos( theta23 )
                double c13 = (k2*k2 - k1*k1 - k3*k3) / (2.0 * k1 * k3); // cos( theta13 )
                // guard against round-off
                if (c12 < -1.0) c12 = -1.0; 
                if (c12 >  1.0) c12 =  1.0;
                if (c23 < -1.0) c23 = -1.0;
                if (c23 >  1.0) c23 =  1.0;
                double theta12 = acos(c12);
                double theta23 = acos(c23);
                double theta13 = acos(c13);

                

                if (theta12 > M_PI && theta23 > M_PI) {continue;}

                // loop over all combinations of (k1,k2,k3) and add weights
                double k[3] = {k1, k2, k3};
                for (int a = 0; a < 3; ++a) {
                    for (int b = 0; b < 3; ++b) {
                        if (a != b) {
                            int c = 3 - a - b;
                            // Uniform (k1,k2,k3) grid: the density per dk1 dk2 dk3 is
                            // the per-(p_tot,f1,f2) density divided by p_tot^2.
                            double p_tot = (k[a] + k[b] + k[c]) * GeVfm;
                            double weight = myGen->get_rho_ptot_f1f2f3(k[a], k[b], k[c]) / (p_tot * p_tot);
                            if (weight > 0) {
                                int bin12 = std::min(int(theta12 / dtheta), theta_bins - 1);
                                int bin23 = std::min(int(theta23 / dtheta), theta_bins - 1);
                                heatmap[bin12][bin23] += weight;
                            }
                        }
                    }
                }
            }
        }
    }

    cout << counter << endl;

    // Write binned result
    std::ofstream out("theta_heatmap_binned.txt");
    out << "theta12_bin_center theta23_bin_center weight\n";
    for (int i = 0; i < theta_bins; ++i) {
        double theta12_center = (i + 0.5) * dtheta;
        for (int j = 0; j < theta_bins; ++j) {
            double theta23_center = (j + 0.5) * dtheta;
            out << theta12_center << " " << theta23_center << " " << heatmap[i][j] << "\n";
        }
    }
    out.close();

    std::cout << "Binned heatmap written to theta_heatmap_binned.txt" << std::endl;
    return 0;
}