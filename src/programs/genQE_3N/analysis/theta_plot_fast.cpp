#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <thread>
#include <mutex>
#include <atomic>
#include "QEGenerator_3N.hh"
#include "constants.hh"
#include "TVector3.h"

using namespace std;

eNCrossSection * myCS;
QEGenerator_3N * myGen;
TRandom3 * myRand;

/*
Works the same as theta_plot, but uses threads to run faster.
*/

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

// Thread-safe data structures
mutex heatmap_mutex;
atomic<int> counter(0);

void compute_chunk(int start_i, int end_i, vector<vector<double>>& global_heatmap) {
    vector<vector<double>> local_heatmap(theta_bins, vector<double>(theta_bins, 0.0));
    double dk = (ktot_max - k_min) / k_bins;
    double dtheta = M_PI / theta_bins;
    int local_counter = 0;

    // Create thread-local random number generator
    TRandom3 local_rand(0);
    
    for (int i = start_i; i < end_i; ++i) {
        double k1 = k_min + i * dk;
        for (int j = 0; j < i; ++j) {
            double k2 = k_min + j * dk;
            for (int m = 0; m < j; ++m) {
                double k3 = k_min + m * dk;
                local_counter++;

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
                                // Bin the angles of this labelling (k[a], k[b], k[c]),
                                // not the sorted (k1, k2, k3): theta12 of the sorted
                                // triangle is always > 90 deg.
                                double perm12, perm23;
                                labelled_angles(k[a], k[b], k[c], perm12, perm23);
                                int bin12 = std::min(int(perm12 / dtheta), theta_bins - 1);
                                int bin23 = std::min(int(perm23 / dtheta), theta_bins - 1);
                                local_heatmap[bin12][bin23] += weight;
                            }
                        }
                    }
                }
            }
        }
    }

    // Merge local results into global heatmap
    lock_guard<mutex> lock(heatmap_mutex);
    for (int i = 0; i < theta_bins; ++i) {
        for (int j = 0; j < theta_bins; ++j) {
            global_heatmap[i][j] += local_heatmap[i][j];
        }
    }
    counter += local_counter;
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

    // Determine number of threads (use hardware concurrency)
    unsigned int num_threads = std::thread::hardware_concurrency();
    if (num_threads == 0) num_threads = 4; // fallback
    
    cout << "Using " << num_threads << " threads for parallel computation" << endl;

    // Divide work among threads
    vector<thread> threads;
    int chunk_size = k_bins / num_threads;
    
    for (unsigned int t = 0; t < num_threads; ++t) {
        int start_i = t * chunk_size;
        int end_i = (t == num_threads - 1) ? k_bins : (t + 1) * chunk_size;
        
        threads.emplace_back(compute_chunk, start_i, end_i, ref(heatmap));
    }

    // Wait for all threads to complete
    for (auto& t : threads) {
        t.join();
    }

    cout << counter << endl;

    // Write binned result
    std::ofstream out("theta_heatmap_binned.txt");
    out << "theta12_bin_center theta23_bin_center weight\n";
    double dtheta = M_PI / theta_bins;
    for (int i = 0; i < theta_bins; ++i) {
        double theta12_center = (i + 0.5) * dtheta;
        for (int j = 0; j < theta_bins; ++j) {
            double theta23_center = (j + 0.5) * dtheta;
            out << theta12_center << " " << theta23_center << " " << heatmap[i][j] << "\n";
        }
    }
    out.close();

    std::cout << "Binned heatmap written to theta_heatmap_binned.txt" << std::endl;
    
    // Cleanup
    delete myGen;
    delete myCS;
    delete myRand;
    
    return 0;
}