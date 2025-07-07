#include "CQ_Calculation.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

Config readConfig(const std::string& configFile) {
	Config config;
	std::ifstream file(configFile);
	std::string line;

	while (std::getline(file, line)) {
		if (line.empty() || line[0] == '#') continue;

		std::istringstream iss(line);
		std::string key;
		if (std::getline(iss, key, '=')) {
			std::string value;
			if (std::getline(iss, value)) {
				if (key == "TotalEvents") config.TotalEvents = std::stoi(value);
				else if (key == "MixedEvents") config.MixedEvents = std::stoi(value);
				else if (key == "MaxParticles") config.MaxParticles = std::stoi(value);
				else if (key == "Q_num_bins") config.Q_num_bins = std::stoi(value);
				else if (key == "kT_num_bins") config.kT_num_bins = std::stoi(value);
				else if (key == "mT_num_bins") config.mT_num_bins = std::stoi(value);
				else if (key == "Q_min") config.Q_min = std::stod(value);
				else if (key == "Q_max") config.Q_max = std::stod(value);
				else if (key == "kT_min") config.kT_min = std::stod(value);
				else if (key == "kT_max") config.kT_max = std::stod(value);
				else if (key == "mT_min") config.mT_min = std::stod(value);
				else if (key == "mT_max") config.mT_max = std::stod(value);
				else if (key == "rapidity_limit") config.rapidity_limit = std::stod(value);
				else if (key == "input_filename") config.input_filename = value;
				else if (key == "pid_list") {
					std::istringstream pid_stream(value);
					int pid;
					while (pid_stream >> pid) {
						config.pid_list.push_back(pid);
						if (pid_stream.peek() == ',') pid_stream.ignore();
					}
				}
			}
		}
	}
	return config;
}

int main() {
	Config config = readConfig("/HBT_Analysis/config/config.ini"); // Check the correct filepath
	auto particles = readParticleData(config);
	cout << "\nWhich correlation analysis would you like to run?\n"
	     << "  1) 1D  [C(Qinv)]\n"
	     << "  2) 3D  [C(Qout, Qside, Qlong)]\n"
	     << "Enter choice [1 or 2]: ";
	int choice = 0;
	if (!(cin >> choice)) {
		cerr << "Invalid input. Exiting.\n";
		return 1;
	}

	switch (choice) {
		case 1:
			cout << "Running 1D correlation [C(Qinv)]…\n";
			calculateInvCorrelation(particles, config);
			std::cout << "\nCompilation complete. Output stored in Cqinv_kT_*.data\n" << std::endl;
		break;

		case 2:
			cout << "Running 3D correlation [C(Qout, Qside, Qlong)]…\n";
			calculate3DCorrelation(particles, config);
			std::cout << "\nCompilation complete. Output stored in Cq3D_kT_*.data\n" << std::endl;
		break;

		default:
			cerr << "Choice must be 1 or 2. Exiting.\n";
		return 1;
	}
	return 0;
}
