#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <map>
#include <sstream>
#include <numeric>
#include "QEGenerator_3N.hh"
#include "constants.hh"
#include "TVector3.h"
#include "TLorentzVector.h"

using namespace std;

eNCrossSection * myCS;
QEGenerator_3N * myGen;
TRandom3 * myRand;

int main(int argc, char **argv) {
    // Defaults
    const double Ebeam = 6.0; // GeV
    int theta_bins = 50;
    long long n_events = 1000000; // number of accepted events to accumulate
    bool use_CM = false;

    // Independent variable: can be "Q2" (default) or "xB" (Bjorken x)
    std::string indep_var = "Q2";

    // Default boundaries for the independent variable (Q2 by default)
    std::vector<double> bounds = {1.0, 4.0, 8.0, 10.0};

    // Command-line argument formats supported (flexible):
    // 1) theta_scan [n_events] [bounds] [theta_bins] [use_CM]
    //    - bounds: comma-separated numbers (assumed Q2)
    // 2) theta_scan [n_events] [var:bounds] [theta_bins] [use_CM]
    //    - var:bounds e.g. "xB:0.1,0.2,0.5" or "Q2:1,4,8,10"
    // 3) theta_scan [n_events] [var] [bounds] [theta_bins] [use_CM]
    //    - var: "xB" and bounds in the next argument

    auto parse_bounds = [&](const std::string &s)->std::vector<double> {
        std::vector<double> parsed;
        std::stringstream ss(s);
        std::string token;
        while (getline(ss, token, ',')) {
            try { parsed.push_back(stod(token)); } catch(...) { }
        }
        return parsed;
    };

    if (argc > 1) {
        n_events = atoll(argv[1]);
    }

    if (argc > 2) {
        std::string arg2 = argv[2];
        // case: var:bounds
        auto pos = arg2.find(':');
        if (pos != std::string::npos) {
            std::string var = arg2.substr(0, pos);
            std::string bstr = arg2.substr(pos+1);
            if (!var.empty()) indep_var = var;
            auto p = parse_bounds(bstr);
            if (p.size() >= 2) bounds = p;
        } else {
            // if arg2 contains any letter, treat as var name and expect bounds in argv[3]
            bool has_alpha = false;
            for (char c : arg2) if (isalpha((unsigned char)c)) { has_alpha = true; break; }
            if (has_alpha) {
                indep_var = arg2;
                if (argc > 3) {
                    auto p = parse_bounds(argv[3]);
                    if (p.size() >= 2) bounds = p;
                }
            } else {
                // otherwise arg2 is bounds (legacy behavior: assume Q2)
                auto p = parse_bounds(arg2);
                if (p.size() >= 2) bounds = p;
            }
        }
    }

    // theta_bins and use_CM may be after bounds depending on format
    // Find a numeric argument that is an integer for theta_bins and use_CM flag
    // For simplicity: assume theta_bins is next positional integer after bounds if provided
    // and use_CM is the following argument
    // Determine which argv index we should read for theta_bins and use_CM
    int arg_idx_for_theta_bins = 3; // default
    // If argv[2] was var with bounds in argv[3], theta_bins is argv[4]
    if (argc > 2) {
        std::string arg2 = argv[2];
        if (arg2.find(':') != std::string::npos) {
            arg_idx_for_theta_bins = 3;
        } else {
            bool has_alpha = false;
            for (char c : arg2) if (isalpha((unsigned char)c)) { has_alpha = true; break; }
            if (has_alpha) arg_idx_for_theta_bins = 4;
            else arg_idx_for_theta_bins = 3;
        }
    }
    if (argc > arg_idx_for_theta_bins) {
        theta_bins = stoi(argv[arg_idx_for_theta_bins]);
    }
    if (argc > arg_idx_for_theta_bins + 1) {
        use_CM = (stoi(argv[arg_idx_for_theta_bins + 1]) != 0);
    }

    // Validate bounds: sort and unique
    sort(bounds.begin(), bounds.end());
    bounds.erase(unique(bounds.begin(), bounds.end()), bounds.end());
    if (bounds.size() < 2) {
        cerr << "Need at least two boundary values.\n";
        return 1;
    }

    int num_regions = static_cast<int>(bounds.size()) - 1;

    cout << "theta_scan: n_events=" << n_events << ", theta_bins=" << theta_bins
        << ", var=" << indep_var << " regions=" << num_regions << "\n";
    cout << indep_var << " boundaries:";
    for (double b : bounds) cout << " " << b;
    cout << "\n";

    // Initialize generator and cross section
    myRand = new TRandom3(0);
    ffModel ffMod = kelly;
    csMethod csMeth = cc1;
    myCS = new eNCrossSection(csMeth, ffMod);
    myGen = new QEGenerator_3N(Ebeam, myCS, 1, myRand);

    const double dThetaDeg = 180.0 / theta_bins;

    // histograms: one theta_bins x theta_bins heatmap per Q2 region
    std::vector<std::vector<std::vector<double>>> hist_regions;
    hist_regions.resize(num_regions, std::vector<std::vector<double>>(theta_bins, std::vector<double>(theta_bins, 0.0)));

    // Event loop: accumulate weights until we have n_events accepted in EACH region
    std::vector<long long> success_counts(num_regions, 0);
    long long event_count = 0;

    // Print initial per-region status
    cout << "Starting generation: target n_events per region=" << n_events << "\n";
    cout << "Initial counts: ";
    for (int r = 0; r < num_regions; ++r) {
        cout << "R" << (r+1) << ":" << success_counts[r];
        if (r < num_regions-1) cout << ", ";
    }
    cout << "\n";

    // continue generating until every region has reached n_events
    while (!std::all_of(success_counts.begin(), success_counts.end(), [&](long long c){ return c >= n_events; })) {
        event_count++;

        double weight;
        int N1_type, N2_type, N3_type;
        TLorentzVector v_k_target, v_Lead_target, v_2_target, v_3_target, v_Am3_target;
        myGen->generate_event(weight, N1_type, N2_type, N3_type, v_k_target, v_Lead_target, v_2_target, v_3_target, v_Am3_target, use_CM);

        if (weight <= 0.0) continue;

        // calculate q
        TLorentzVector vbeam_target(0.0, 0.0, Ebeam, Ebeam);
        TLorentzVector q = vbeam_target - v_k_target;

        // incoming momenta (reconstructed struck and spectators)
        TVector3 p1 = v_Lead_target.Vect() - q.Vect();
        TVector3 p2 = v_2_target.Vect();
        TVector3 p3 = v_3_target.Vect();

        // Basic triangle check (same as other code) to remove unphysical events
        if (p1.Mag() > p2.Mag() + p3.Mag() || p2.Mag() > p1.Mag() + p3.Mag() || p3.Mag() > p1.Mag() + p2.Mag()) {
            continue;
        }

        // sort by magnitude
        std::vector<std::pair<TVector3, int>> particles = {{p1,1}, {p2,2}, {p3,3}};
        std::sort(particles.begin(), particles.end(), [](const auto &a, const auto &b){ return a.first.Mag() > b.first.Mag(); });
        TVector3 p1_sorted = particles[0].first;
        TVector3 p2_sorted = particles[1].first;
        TVector3 p3_sorted = particles[2].first;

        double theta12 = p1_sorted.Angle(p2_sorted) * 180.0 / M_PI;
        double theta23 = p2_sorted.Angle(p3_sorted) * 180.0 / M_PI;

        double Q2 = -q.Mag2();

        // compute independent variable value based on indep_var
        double indep_val = 0.0;
        if (indep_var == "Q2" || indep_var == "q2" || indep_var == "Q^2") {
            indep_val = Q2;
        } else if (indep_var == "xB" || indep_var == "xb" || indep_var == "Xb") {
            // Bjorken x = Q2 / (2 * mN * nu)
            double nu = q.E();
            // Protect against zero/negative nu
            if (nu <= 0.0) continue;
            indep_val = Q2 / (2.0 * mN * nu);
        } else {
            // unknown variable: default to Q2
            indep_val = Q2;
        }

        // find region for indep_val
        int region = -1;
        for (int r = 0; r < num_regions; ++r) {
            if (indep_val >= bounds[r] && indep_val < bounds[r+1]) { region = r; break; }
        }
        // If outside all regions, skip
        if (region == -1) continue;

        // If this region already has enough accepted events, skip storing it
        if (success_counts[region] >= n_events) continue;

        // bin theta12/theta23
        if (theta12 >= 0.0 && theta12 <= 180.0 && theta23 >= 0.0 && theta23 <= 180.0) {
            int it1 = static_cast<int>(theta12 / dThetaDeg);
            int it2 = static_cast<int>(theta23 / dThetaDeg);
            if (it1 >= theta_bins) it1 = theta_bins - 1;
            if (it2 >= theta_bins) it2 = theta_bins - 1;
            hist_regions[region][it2][it1] += weight;
        }

        success_counts[region]++;
        // announce when a region first reaches the target
        if (success_counts[region] == n_events) {
            cout << "Region " << (region+1) << " reached target n_events=" << n_events << "\n";
        }

        // optional: small progress print every 10k generated showing per-region counts
        if ((event_count % 10000) == 0) {
            long long total_accepted = std::accumulate(success_counts.begin(), success_counts.end(), 0LL);
            cout << "events generated=" << event_count << ", accepted(total)=" << total_accepted << " (";
            for (int rr = 0; rr < num_regions; ++rr) {
                cout << "R" << (rr+1) << ":" << success_counts[rr] << "/" << n_events;
                if (rr < num_regions-1) cout << ", ";
            }
            cout << ")\r" << flush;
        }
    }

    // Write out one file per region for the chosen independent variable
    for (int r = 0; r < num_regions; ++r) {
        double v_min = bounds[r];
        double v_max = bounds[r+1];
        std::ostringstream fname;
        // sanitize variable name for filename
        std::string varname = indep_var;
        for (auto &c : varname) if (!isalnum((unsigned char)c)) c = '_';
        fname << "hist_theta12_theta23_" << varname << "_region" << (r+1) << ".txt";
        std::ofstream hout(fname.str());
        hout << std::fixed << std::setprecision(25);
        // Header line that the plotting script expects
        hout << "# " << varname << " region " << (r+1) << ": " << v_min << " ≤ " << varname << " < " << v_max << "\n";
        hout << "# theta_bins " << theta_bins << " range_deg [0,180]\n";
        hout << "# Columns: theta12_center theta23_center weight\n";
        // include how many accepted events were accumulated for this region
        hout << "# accepted_events " << success_counts[r] << " target_n_events " << n_events << "\n";
        for (int it2 = 0; it2 < theta_bins; ++it2) {
            double th23_center = (it2 + 0.5) * dThetaDeg;
            for (int it1 = 0; it1 < theta_bins; ++it1) {
                double th12_center = (it1 + 0.5) * dThetaDeg;
                double w = hist_regions[r][it2][it1];
                hout << th12_center << " " << th23_center << " " << w << "\n";
            }
        }
        cout << "Wrote " << fname.str() << " (" << varname << " in [" << v_min << "," << v_max << "))\n";
    }

    cout << "theta_scan: finished. events generated=" << event_count << ", accepted per region: ";
    for (int r = 0; r < num_regions; ++r) {
        cout << success_counts[r];
        if (r < num_regions-1) cout << ", ";
    }
    cout << "\n";
    return 0;
}
