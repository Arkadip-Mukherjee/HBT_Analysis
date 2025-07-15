#include "CQ_Calculation.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <omp.h>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <atomic>
#include <random>

using namespace std;

const double hbar_c = 0.197326980;

// Function to get PID from UrQMD MCID
int get_PID_from_urqmd_MCID(int mcid, int iso) {
	if (mcid == 101 && iso == 2) return 211;   // pi^{+}
	if (mcid == 101 && iso == 0) return 111;   // pi^{0}
	if (mcid == 101 && iso == -2) return -211; // pi^{-}
	if (mcid == 106 && iso == 1) return 321;   // k^{+}
	if (mcid == 106 && iso == -1) return 311;  // k^{0}
	if (mcid == -106 && iso == 1) return -311; // anti k^{0}
	if (mcid == -106 && iso == -1) return -321; // k^{-}
	if (mcid == 109 && iso == 0) return 333;   // phi(1020)
	if (mcid == 102 && iso == 0) return 221;   // eta
	if (mcid == 100 && iso == 0) return 22;    // photon
	if (mcid == 1 && iso == 1) return 2212;    // p
	if (mcid == 1 && iso == -1) return 2112;   // n
	if (mcid == -1 && iso == -1) return -2212; // anti p
	if (mcid == -1 && iso == 1) return -2112;  // anti n
	if (mcid == 40 && iso == 2) return 3222;   // Sigma^{+}
	if (mcid == -40 && iso == -2) return -3222; // anti Sigma^{-}
	if (mcid == 40 && iso == 0) return 3212;   // Sigma^{0}
	if (mcid == -40 && iso == 0) return -3212; // anti Sigma^{0}
	if (mcid == 40 && iso == -2) return 3112;  // Sigma^{-}
	if (mcid == -40 && iso == 2) return -3112; // anti Sigma^{+}
	if (mcid == 49 && iso == -1) return 3312;  // Xi^{-}
	if (mcid == 49 && iso == 1) return 3322;   // Xi^{0}
	if (mcid == -49 && iso == -1) return -3322; // anti Xi^{0}
	if (mcid == -49 && iso == 1) return -3312; // anti Xi^{+}
	if (mcid == 27 && iso == 0) return 3122;   // Lambda
	if (mcid == -27 && iso == 0) return -3122; // anti Lambda
	if (mcid == 55 && iso == 0) return 3334;   // Omega
	if (mcid == -55 && iso == 0) return -3334; // anti Omega
	return 0;
}

// Function to Calculate Rapidity
double rapidity_calculation(double e, double pz){
	double rap = 0.5 * std::log(abs((e + pz)/(e - pz)));
	return std::abs(rap);
}

// Function to read the Particle Data Set from input_file
std::vector<std::vector<Particle>> readParticleData(const Config& config) {
	std::ifstream file(config.input_filename, std::ios::binary | std::ios::in);
	if (!file) {
	throw std::runtime_error("File not found: " + config.input_filename);
	}

	std::vector<std::vector<Particle>> particles(config.TotalEvents);
	for (int event = 0; event < config.TotalEvents; event++) {
		int nch, count = 0;
		file.read(reinterpret_cast<char*>(&nch), sizeof(int));
		if (file.eof()) break;

		particles[event].reserve(nch);
		for (int i = 0; i < nch; i++) {
			float particle_array[11];
			file.read(reinterpret_cast<char*>(particle_array), sizeof(float) * 11);

			int pid_store = get_PID_from_urqmd_MCID(static_cast<int>(particle_array[1]), static_cast<int>(particle_array[2]));
			if (std::find(config.pid_list.begin(), config.pid_list.end(), pid_store) != config.pid_list.end()) {
				count++;
				Particle p;
				p.pid = pid_store;
				p.mass = particle_array[0]
				p.t = particle_array[3];
				p.x = particle_array[4];
				p.y = particle_array[5];
				p.z = particle_array[6];	
				p.e = particle_array[7];
				p.px = particle_array[8];
				p.py = particle_array[9];
				p.pz = particle_array[10];
				
				if (rapidity_calculation(p.e, p.pz) > config.rapidity_limit) continue;
				particles[event].push_back(p);
			}
		} 
		if (event%100 == 0) std::cout<< "Particle count for event " << event << " is:   " << count << std::endl ;
	}
	return particles;
}

// Function to compute the 1-D Correlation function 
void calculateInvCorrelation(const std::vector<std::vector<Particle>>& particles, const Config& config) {

	const double Q_bin_width = (config.Q_max - config.Q_min) / config.Q_num_bins;
	const double kT_bin_width = (config.kT_max - config.kT_min) / config.kT_num_bins;
	
	std::atomic<long long> N_same(0), N_mixed(0);
	
	std::vector<std::vector<double>> numerator(config.kT_num_bins, std::vector<double>(config.Q_num_bins, 0.0));
	std::vector<std::vector<double>> denominator(config.kT_num_bins, std::vector<double>(config.Q_num_bins, 0.0));
	std::vector<std::vector<double>> Qinv_Store(config.kT_num_bins, std::vector<double>(config.Q_num_bins, 0.0));
	std::vector<std::vector<double>> Num_Pairs(config.kT_num_bins, std::vector<double>(config.Q_num_bins, 0.0));

	// Declare thread-local versions
	std::vector<std::vector<std::vector<double>>> local_numerators(omp_get_max_threads(), numerator);
	std::vector<std::vector<std::vector<double>>> local_denominators(omp_get_max_threads(), denominator);
	std::vector<std::vector<std::vector<double>>> local_Qinv_Store(omp_get_max_threads(), Qinv_Store);
	std::vector<std::vector<std::vector<double>>> local_Num_Pairs(omp_get_max_threads(), Num_Pairs);
	
	
	// =============================================  FOR SAME EVENT  ============================================== //


	#pragma omp parallel for
	for (int event = 0; event < config.TotalEvents; event++) {

		int tid = omp_get_thread_num();
		auto& local_num = local_numerators[tid];
		auto& local_Qinvstore = local_Qinv_Store[tid];
		auto& local_NPairs = local_Num_Pairs[tid];

		const auto& ev_particles = particles[event];
		for (size_t i = 0; i < ev_particles.size(); i++) {
			const auto& p1 = ev_particles[i];
			for (size_t j = i+1; j < ev_particles.size(); j++) {
				const auto& p2 = ev_particles[j];

				if (p1.pid != p2.pid) continue; //To ensure pairing between identical particles only

				double Kx = 0.5*(p1.px + p2.px);
				double Ky = 0.5*(p1.py + p2.py);
				double Kz = 0.5*(p1.pz + p2.pz);
				double Ke = 0.5*(p1.e + p2.e);

				double K_perp_sq = Kx * Kx + Ky * Ky;
				double K_perp = std::sqrt(K_perp_sq);
				
				double Mt_sq = Ke * Ke - Kz * Kz;
				double Mt = std::sqrt(Mt_sq);

				if (K_perp < config.kT_min || K_perp > config.kT_max) continue;
				int Kperp_index = static_cast<int>((K_perp - config.kT_min) / kT_bin_width);
				if (Kperp_index >= config.kT_num_bins) continue;
				
				double qx = p1.px - p2.px;
				double qy = p1.py - p2.py;
				double qz = p1.pz - p2.pz;
				double qe = p1.e  - p2.e;
				
				double qinv = std::sqrt(-(qe*qe - qx*qx - qy*qy - qz*qz));
				if (qinv < config.Q_min || qinv > config.Q_max) continue;
				int qinv_index = static_cast<int>((qinv - config.Q_min) / Q_bin_width);
				if (qinv_index >= config.Q_num_bins) continue;

				double cosqx = std::cos((qe * (p1.t - p2.t) - qx * (p1.x - p2.x) 
				- qy * (p1.y - p2.y) - qz * (p1.z - p2.z)) / hbar_c);
				
				local_NPairs[Kperp_index][qinv_index] += 1.0;
				local_num[Kperp_index][qinv_index] += cosqx;
				local_Qinvstore[Kperp_index][qinv_index] += qinv;
				
				N_same++;
			}
		}
		if (event%100 == 0) std::cout << "Calculations for Same Event Number:  " << event << "  done." << std::endl ;
	}
	
	
	// =============================================  FOR MIXED EVENTS  ============================================= //
		
	#pragma omp parallel
	{
	    thread_local std::mt19937 gen(std::random_device{}());
	    thread_local std::uniform_real_distribution<double> dist(0.0, 2*M_PI);

	#pragma omp for
	for (int event = 0; event < config.TotalEvents; event++) {

		int tid = omp_get_thread_num();
		auto& local_den = local_denominators[tid];
		const auto& ev_particles = particles[event];
		
	for (int mixedev = 0; mixedev < config.MixedEvents; mixedev++){
	
		const double random_angle_phi = dist(gen);
		const double cosphi = std::cos(random_angle_phi);
		const double sinphi = std::sin(random_angle_phi);
	
		for (size_t i = 0; i < ev_particles.size(); i++) {
			const auto& p1 = ev_particles[i];
			
			double new_px1 = p1.px * cosphi - p1.py * sinphi;
			double new_py1 = p1.px * sinphi + p1.py * cosphi;
			
			for (size_t j = 0; j < ev_particles.size(); j++) {
				const auto& p2 = ev_particles[j];

				if (p1.pid != p2.pid) continue; //To ensure pairing between identical particles only
				
				double Kx = 0.5*(new_px1 + p2.px);
				double Ky = 0.5*(new_py1 + p2.py);
				double Kz = 0.5*(p1.pz + p2.pz);
				double Ke = 0.5*(p1.e + p2.e);

				double K_perp_sq = Kx * Kx + Ky * Ky;
				double K_perp = std::sqrt(K_perp_sq);
				
				double Mt_sq = Ke * Ke - Kz * Kz;
				double Mt = std::sqrt(Mt_sq);

				if (K_perp < config.kT_min || K_perp > config.kT_max) continue;
				int Kperp_index = static_cast<int>((K_perp - config.kT_min) / kT_bin_width);
				if (Kperp_index >= config.kT_num_bins) continue;
				
				double qx = new_px1 - p2.px;
				double qy = new_py1 - p2.py;
				double qz = p1.pz - p2.pz;
				double qe = p1.e  - p2.e;

				double qinv = std::sqrt(-(qe*qe - qx*qx - qy*qy - qz*qz));
				if (qinv < config.Q_min || qinv > config.Q_max) continue;
				int qinv_index = static_cast<int>((qinv - config.Q_min) / Q_bin_width);
				if (qinv_index >= config.Q_num_bins) continue;
				
				local_den[Kperp_index][qinv_index] += 1.0;	
	
				N_mixed++;		
			}
		}
		if (event%100 == 0) std::cout << "Calculations for Mixed Event Number:  " << event << "  done." << std::endl ;
	}
	}
	
	// Merge thread-local results
	for (int tid = 0; tid < omp_get_max_threads(); tid++) {
		for (int kt = 0; kt < config.kT_num_bins; kt++) {
			for (int ii = 0; ii < config.Q_num_bins; ii++) {
				numerator[kt][ii] += local_numerators[tid][kt][ii];
				denominator[kt][ii] += local_denominators[tid][kt][ii];
				Qinv_Store[kt][ii] += local_Qinv_Store[tid][kt][ii];
				Num_Pairs[kt][ii] += local_Num_Pairs[tid][kt][ii];
			}
		}
	}

	// Write output files
	for (int kTbin = 0; kTbin < config.kT_num_bins; kTbin++) { 
		std::string filename = "Cqinv_kT_" + std::to_string(kTbin) + ".data"; // CHANGE THE FILENAME AS PER CONVENIENCE
		std::ofstream file(filename);
		for (int ii = 0; ii < config.Q_num_bins; ++ii) {
			if (denominator[kTbin][ii] == 0) continue;
			if (numerator[kTbin][ii] == 0) continue;
			if (Num_Pairs[kTbin][ii] == 0) continue;

			double correlfunc = 1 + ((numerator[kTbin][ii] / N_same) / (denominator[kTbin][ii] / N_mixed));
			double Qinv_val = Qinv_Store[kTbin][ii] / Num_Pairs[kTbin][ii] ;

			file << std::scientific << Qinv_val << "   " << correlfunc << "\n";
		}
		file.close();
	}
}


// Function to compute the 3-D Correlation function 
void calculate3DCorrelation(const std::vector<std::vector<Particle>>& particles, const Config& config) {

	const double Q_bin_width = (config.Q_max - config.Q_min) / config.Q_num_bins;
	const double kT_bin_width = (config.kT_max - config.kT_min) / config.kT_num_bins;
	
	std::atomic<long long> N_same(0), N_mixed(0);
	
	std::vector<std::vector<std::vector<std::vector<double>>>> numerator(config.kT_num_bins,
        std::vector<std::vector<std::vector<double>>>(config.Q_num_bins,
        std::vector<std::vector<double>>(config.Q_num_bins,
        std::vector<double>(config.Q_num_bins, 0.0))));

    	std::vector<std::vector<std::vector<std::vector<double>>>> denominator(config.kT_num_bins,
        std::vector<std::vector<std::vector<double>>>(config.Q_num_bins,
        std::vector<std::vector<double>>(config.Q_num_bins,
        std::vector<double>(config.Q_num_bins, 0.0))));
        
        std::vector<std::vector<std::vector<std::vector<double>>>> Num_Pairs(config.kT_num_bins,
        std::vector<std::vector<std::vector<double>>>(config.Q_num_bins,
        std::vector<std::vector<double>>(config.Q_num_bins,
        std::vector<double>(config.Q_num_bins, 0.0))));
        
        std::vector<std::vector<std::vector<std::vector<double>>>> Q_out(config.kT_num_bins,
        std::vector<std::vector<std::vector<double>>>(config.Q_num_bins,
        std::vector<std::vector<double>>(config.Q_num_bins,
        std::vector<double>(config.Q_num_bins, 0.0))));
        
        std::vector<std::vector<std::vector<std::vector<double>>>> Q_side(config.kT_num_bins,
        std::vector<std::vector<std::vector<double>>>(config.Q_num_bins,
        std::vector<std::vector<double>>(config.Q_num_bins,
        std::vector<double>(config.Q_num_bins, 0.0))));
        
        std::vector<std::vector<std::vector<std::vector<double>>>> Q_long(config.kT_num_bins,
        std::vector<std::vector<std::vector<double>>>(config.Q_num_bins,
        std::vector<std::vector<double>>(config.Q_num_bins,
        std::vector<double>(config.Q_num_bins, 0.0))));
        

	// Declare thread-local versions
	std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> local_numerators(omp_get_max_threads(),
	numerator);
	std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> local_denominators(omp_get_max_threads(),
	denominator);
	std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> local_Num_Pairs(omp_get_max_threads(), Num_Pairs);
	std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> local_Q_out(omp_get_max_threads(), Q_out);
	std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> local_Q_side(omp_get_max_threads(), Q_side);
	std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>> local_Q_long(omp_get_max_threads(), Q_long);

	
	// =============================================  FOR SAME EVENT  ============================================== //


	#pragma omp parallel for
	for (int event = 0; event < config.TotalEvents; event++) {

		int tid = omp_get_thread_num();
		auto& local_num = local_numerators[tid];
		auto& local_NPairs = local_Num_Pairs[tid];
		auto& local_Qo = local_Q_out[tid];
		auto& local_Qs = local_Q_side[tid];
		auto& local_Ql = local_Q_long[tid];

		const auto& ev_particles = particles[event];
		for (size_t i = 0; i < ev_particles.size(); i++) {
			const auto& p1 = ev_particles[i];
			for (size_t j = i+1; j < ev_particles.size(); j++) {
				const auto& p2 = ev_particles[j];

				if (p1.pid != p2.pid) continue; //To ensure pairing between identical particles only
				
				double Kx = 0.5*(p1.px + p2.px);
				double Ky = 0.5*(p1.py + p2.py);
				double Kz = 0.5*(p1.pz + p2.pz);
				double Ke = 0.5*(p1.e + p2.e);

				double K_perp_sq = Kx * Kx + Ky * Ky;
				double K_perp = std::sqrt(K_perp_sq);
				
				double Mt_sq = Ke * Ke - Kz * Kz;
				double Mt = std::sqrt(Mt_sq);

				if (K_perp < config.kT_min || K_perp > config.kT_max) continue;
				int Kperp_index = static_cast<int>((K_perp - config.kT_min) / kT_bin_width);
				if (Kperp_index >= config.kT_num_bins) continue;
				
				double phi = std::atan2(Ky, Kx); // Angle of K_perp in lab frame
				double cos_phi = std::cos(-phi);
				double sin_phi = std::sin(-phi);

				// Rotating transverse momentum of p1
				double p1_out = p1.px * cos_phi - p1.py * sin_phi;
				double p1_side = p1.px * sin_phi + p1.py * cos_phi;

				// Rotating transverse momentum of p2
				double p2_out = p2.px * cos_phi - p2.py * sin_phi;
				double p2_side = p2.px * sin_phi + p2.py * cos_phi;

				// Rotating position vectors of p1
				double out1 = p1.x * cos_phi - p1.y * sin_phi;
				double side1 = p1.x * sin_phi + p1.y * cos_phi;

				// Rotating position vectors of p2
				double out2 = p2.x * cos_phi - p2.y * sin_phi;
				double side2 = p2.x * sin_phi + p2.y * cos_phi;

				// Transformation to LCMS frame
				double beta = Kz / Ke;
				double gamma = Ke / Mt;
				
				double pz1LCMS = gamma * (p1.pz - beta * p1.e);
				double pz2LCMS = gamma * (p2.pz - beta * p2.e);
				double e1LCMS = gamma * (p1.e - beta * p1.pz);
				double e2LCMS = gamma * (p2.e - beta * p2.pz);

				double z1LCMS = gamma * (p1.z - beta * p1.t);
				double z2LCMS = gamma * (p2.z - beta * p2.t);
				double t1LCMS = gamma * (p1.t - beta * p1.z);
				double t2LCMS = gamma * (p2.t - beta * p2.z);
				
				double qo = p1_out - p2_out;
				double qs = p1_side - p2_side;
				double ql = pz1LCMS - pz2LCMS;
				double qe = e1LCMS - e2LCMS;

				int qout_index = static_cast<int>((qo - config.Q_min) / Q_bin_width);
				int qside_index = static_cast<int>((qs - config.Q_min) / Q_bin_width);
				int qlong_index = static_cast<int>((ql - config.Q_min) / Q_bin_width);

				if (qout_index  < 0 || qout_index  >= config.Q_num_bins) continue;
				if (qside_index < 0 || qside_index >= config.Q_num_bins) continue;
				if (qlong_index < 0 || qlong_index >= config.Q_num_bins) continue;
				
				double dout = out1 - out2;
				double dside = side1 - side2;
				double dlong = z1LCMS - z2LCMS;
				double dtime = t1LCMS - t2LCMS;

				double cosqx = std::cos((qe*dtime - qo*dout - qs*dside - ql*dlong) / hbar_c);
				
				local_NPairs[Kperp_index][qout_index][qside_index][qlong_index] += 1.0;
				local_num[Kperp_index][qout_index][qside_index][qlong_index] += cosqx;
				local_Qo[Kperp_index][qout_index][qside_index][qlong_index] += qo;
				local_Qs[Kperp_index][qout_index][qside_index][qlong_index] += qs;
				local_Ql[Kperp_index][qout_index][qside_index][qlong_index] += ql;
				
				N_same++;	
			}
		}
		if (event%100 == 0) std::cout << "Calculations for Same Event Number:  " << event << "  done." << std::endl ;
	}
	
	
	// =============================================  FOR MIXED EVENTS  ============================================= //
		
	#pragma omp parallel
	{
	    thread_local std::mt19937 gen(std::random_device{}());
	    thread_local std::uniform_real_distribution<double> dist(0.0, 2*M_PI);

	#pragma omp for
	for (int event = 0; event < config.TotalEvents; event++) {

		int tid = omp_get_thread_num();
		auto& local_den = local_denominators[tid];
		const auto& ev_particles = particles[event];
		
	for (int mixedev = 0; mixedev < config.MixedEvents; mixedev++){

		const double random_angle_phi = dist(gen);
		const double cosphi = std::cos(random_angle_phi);
		const double sinphi = std::sin(random_angle_phi);

		for (size_t i = 0; i < ev_particles.size(); i++) {
			const auto& p1 = ev_particles[i];
			
			double new_px1 = p1.px * cosphi - p1.py * sinphi;
			double new_py1 = p1.px * sinphi + p1.py * cosphi;
			
			for (size_t j = 0; j < ev_particles.size(); j++) {
				const auto& p2 = ev_particles[j];

				if (p1.pid != p2.pid) continue; //To ensure pairing between identical particles only

				double Kx = 0.5*(new_px1 + p2.px);
				double Ky = 0.5*(new_py1 + p2.py);
				double Kz = 0.5*(p1.pz + p2.pz);
				double Ke = 0.5*(p1.e + p2.e);

				double K_perp_sq = Kx * Kx + Ky * Ky;
				double K_perp = std::sqrt(K_perp_sq);
				
				double Mt_sq = Ke * Ke - Kz * Kz;
				double Mt = std::sqrt(Mt_sq);

				if (K_perp < config.kT_min || K_perp > config.kT_max) continue;
				int Kperp_index = static_cast<int>((K_perp - config.kT_min) / kT_bin_width);
				if (Kperp_index >= config.kT_num_bins) continue;
				
				double phi = std::atan2(Ky, Kx); // Angle of K_perp in lab frame
				double cos_phi = std::cos(-phi);
				double sin_phi = std::sin(-phi);

				// Rotating transverse momentum of p1
				double p1_out = p1.px * cos_phi - p1.py * sin_phi;
				double p1_side = p1.px * sin_phi + p1.py * cos_phi;

				// Rotating transverse momentum of p2
				double p2_out = p2.px * cos_phi - p2.py * sin_phi;
				double p2_side = p2.px * sin_phi + p2.py * cos_phi;

				// Transformation to LCMS frame
				double beta = Kz / Ke;
				double gamma = Ke / Mt;
				
				double pz1LCMS = gamma * (p1.pz - beta * p1.e);
				double pz2LCMS = gamma * (p2.pz - beta * p2.e);
				double e1LCMS = gamma * (p1.e - beta * p1.pz);
				double e2LCMS = gamma * (p2.e - beta * p2.pz);

				double qo = p1_out - p2_out;
				double qs = p1_side - p2_side;
				double ql = pz1LCMS - pz2LCMS;
				double qe = e1LCMS - e2LCMS;

				int qout_index  = static_cast<int>((qo - config.Q_min) / Q_bin_width);
				int qside_index = static_cast<int>((qs - config.Q_min) / Q_bin_width);
				int qlong_index = static_cast<int>((ql - config.Q_min) / Q_bin_width);

				if (qout_index  < 0 || qout_index  >= config.Q_num_bins) continue;
				if (qside_index < 0 || qside_index >= config.Q_num_bins) continue;
				if (qlong_index < 0 || qlong_index >= config.Q_num_bins) continue;
				
				local_den[Kperp_index][qout_index][qside_index][qlong_index] += 1.0;	
				
				N_mixed++;		
			}
		}
	}
	if (event%100 == 0) std::cout << "Calculations for Mixed Event Number:  " << event << "  done." << std::endl ;
	}
	}
	
	// Merge thread-local results
	for (int tid = 0; tid < omp_get_max_threads(); tid++) {
		for (int kt = 0; kt < config.kT_num_bins; kt++) {
			for (int ii = 0; ii < config.Q_num_bins; ii++) {
				for (int jj = 0; jj < config.Q_num_bins; jj++) {
					for (int kk = 0; kk < config.Q_num_bins; kk++) {
						numerator[kt][ii][jj][kk] += local_numerators[tid][kt][ii][jj][kk];
						denominator[kt][ii][jj][kk] += local_denominators[tid][kt][ii][jj][kk];
						Num_Pairs[kt][ii][jj][kk] += local_Num_Pairs[tid][kt][ii][jj][kk];
						Q_out[kt][ii][jj][kk] += local_Q_out[tid][kt][ii][jj][kk];
						Q_side[kt][ii][jj][kk] += local_Q_side[tid][kt][ii][jj][kk];
						Q_long[kt][ii][jj][kk] += local_Q_long[tid][kt][ii][jj][kk];
					}
				}
			}
		}
	}

	// Write output files
	for (int kTbin = 0; kTbin < config.kT_num_bins; kTbin++) { 
		std::string filename = "Cq3D_kT_" + std::to_string(kTbin) + ".data"; // CHANGE THE FILENAME AS PER CONVENIENCE
		std::ofstream file(filename);
		for (int ii = 0; ii < config.Q_num_bins; ii++) {
			for (int jj = 0; jj < config.Q_num_bins; jj++) {
				for (int kk = 0; kk < config.Q_num_bins; kk++) {
					if (denominator[kTbin][ii][jj][kk] == 0) continue;
					if (numerator[kTbin][ii][jj][kk] == 0) continue;
					if (Num_Pairs[kTbin][ii][jj][kk] == 0) continue;
					
					double correlfunc = 1 + ((numerator[kTbin][ii][jj][kk] / N_same) / (denominator[kTbin][ii][jj][kk] / N_mixed));
								
					double Qout_val = Q_out[kTbin][ii][jj][kk] / Num_Pairs[kTbin][ii][jj][kk] ;
					double Qside_val = Q_side[kTbin][ii][jj][kk] / Num_Pairs[kTbin][ii][jj][kk] ;
					double Qlong_val = Q_long[kTbin][ii][jj][kk] / Num_Pairs[kTbin][ii][jj][kk] ;

					file << std::scientific << Qout_val  << "   " << Qside_val  << "   " 
								<< Qlong_val << "   " << correlfunc << " \n";
				}
			}
		}
		file.close();
	}
}
