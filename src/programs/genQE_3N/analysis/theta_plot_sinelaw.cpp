#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "QEGenerator_3N.hh"
#include "constants.hh"
#include "TVector3.h"

using namespace std;

eNCrossSection * myCS;
QEGenerator_3N * myGen;
TRandom3 * myRand;

// Parameters
const double k_min = 0.35/GeVfm; // p_min 350MeV/c
const double ktot_max = 5/GeVfm; // arbitrary 
int k_bins = 50; // Default value
int theta_bins = 50; // Number of bins for theta

// Angles between nucleons 1-2 and 2-3 when nucleon 1 has |k| = ka, nucleon 2
// has kb and nucleon 3 has kc (ka + kb + kc = 0 as vectors), from the law of
// cosines. Same labelling as get_rho_ptot_f1f2f3(ka, kb, kc).
static void labelled_angles(double ka, double kb, double kc, double &theta12, double &theta23) {
    double c12 = (kc*kc - ka*ka - kb*kb) / (2.0 * ka * kb);
    double c23 = (ka*ka - kb*kb - kc*kc) / (2.0 * kb * kc);
    theta12 = acos(std::max(-1.0, std::min(1.0, c12)));
    theta23 = acos(std::max(-1.0, std::min(1.0, c23)));
}

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

    // Sample k1, theta12, and theta23 directly
    for (int i = 0; i < k_bins; ++i) {
        double k1 = k_min + i * dk;
        
        for (int j = 0; j < theta_bins; ++j) {
            double theta12 = j * dtheta;
            
            for (int m = 0; m < theta_bins; ++m) {
                double theta23 = m * dtheta;
                counter++;

                // Use sine law in the triangle formed by k1, k2, k3
                // The triangle has sides k1, k2, k3 and opposite angles theta23, theta13, theta12
                // From momentum conservation: k1 + k2 + k3 = 0 (vector sum)
                // The third angle theta13 = π - theta12 - theta23
                double phi12 = M_PI - theta12;
                double phi23 = M_PI - theta23;
                double phi13 = M_PI - phi12 - phi23;
                
                bool found_solution = false;
                double k2, k3;
                
                // Check if angles are valid (must sum to π and be positive)
                if (phi13 > 0 && phi12 > 0 && phi23 > 0) {
                    // Apply sine law:
                    // k1/sin(theta23) = k2/sin(theta13) = k3/sin(theta12)

                    double sin12 = sin(phi12);
                    double sin23 = sin(phi23);
                    double sin13 = sin(phi13);

                    // Avoid division by zero
                    if (sin23 > 1e-10 && sin13 > 1e-10 && sin12 > 1e-10) {
                        k2 = k1 * sin13 / sin23;
                        k3 = k1 * sin12 / sin23;
                        
                        // Check if k2 and k3 are within valid ranges
                        if (k2 >= k_min && k2 <= ktot_max && k3 >= k_min && k3 <= ktot_max) {
                            found_solution = true;
                        }
                    }
                }
                
                if (found_solution) {
                    // Sort momenta in descending order and recalculate angles
                    if (k1 < k2) swap(k1, k2);
                    if (k2 < k3) swap(k2, k3);
                    if (k1 < k2) swap(k1, k2);

                    // Recalculate angles after sorting using cosine law
                    theta12 = acos((k3*k3 - k1*k1 - k2*k2) / (2.0 * k1 * k2));
                    theta23 = acos((k1*k1 - k2*k2 - k3*k3) / (2.0 * k2 * k3));

                    double J = (k2 * k3); 

                    // Loop over all combinations of (k1,k2,k3) and add weights with Jacobian
                    double k[3] = {k1, k2, k3};
                    for (int a = 0; a < 3; ++a) {
                        for (int b = 0; b < 3; ++b) {
                            if (a != b) {
                                int c = 3 - a - b;
                                double weight = 1; //myGen->get_rho_ptot_f1f2f3(k[a], k[b], k[c]);
                                if (weight > 0) {
                                    // Bin the angles of this labelling and accumulate:
                                    // the sorted theta12 is always > 90 deg, and '='
                                    // kept only the last permutation.
                                    double perm12, perm23;
                                    labelled_angles(k[a], k[b], k[c], perm12, perm23);
                                    int bin12 = std::min(int(perm12 / dtheta), theta_bins - 1);
                                    int bin23 = std::min(int(perm23 / dtheta), theta_bins - 1);
                                    heatmap[bin12][bin23] += weight * J * dk * dtheta * dtheta;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    cout << counter << endl;

    // Write binned result
    std::ofstream out("theta_heatmap_binned_2.txt");
    out << "theta12_bin_center theta23_bin_center weight\n";
    for (int i = 0; i < theta_bins; ++i) {
        double theta12_center = (i + 0.5) * dtheta;
        for (int j = 0; j < theta_bins; ++j) {
            double theta23_center = (j + 0.5) * dtheta;
            out << theta12_center << " " << theta23_center << " " << heatmap[i][j] << "\n";
        }
    }
    out.close();

    std::cout << "Binned heatmap written to theta_heatmap_binned_2.txt" << std::endl;
    return 0;
}
