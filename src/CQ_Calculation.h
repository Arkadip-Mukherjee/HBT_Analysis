#ifndef CQ_CALCULATION_H
#define CQ_CALCULATION_H

#include <vector>
#include <atomic>
#include <string>

struct Particle {
    double e, px, py, pz, x, y, z, t;
};

struct Config {
    int TotalEvents;
    int MixedEvents;
    int MaxParticles;
    int Q_num_bins;
    int kT_num_bins;
    int mT_num_bins;
    double Q_min;
    double Q_max;
    double kT_min;
    double kT_max;
    double mT_min;
    double mT_max;
    double rapidity_limit;
    std::string input_filename;
    std::vector<int> pid_list;
};

int get_PID_from_urqmd_MCID(int mcid, int iso);

double rapidity_calculation(double e, double pz);

std::vector<std::vector<Particle>> readParticleData(const Config& config);

void calculateInvCorrelation(const std::vector<std::vector<Particle>>& particles, const Config& config);

void calculate3DCorrelation(const std::vector<std::vector<Particle>>& particles, const Config& config);

#endif // CQ_CALCULATION_H
