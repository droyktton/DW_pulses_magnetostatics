#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <tuple>


class ElasticLine {
public:
    double dt;
    double t;
    int L;
    std::vector<double> monomer;
    std::vector<double> fuerza;
    std::vector<double> disorder;

    double Temp = 0.0; // temperature
    double f0 = 0.6; // pulse amplitude
    double k0 = 0.0001; // magnetostatic constant
    double DisAmp = 1.0; // disorder amplitude
    double tau = 100.0; // pulse duration
    double tau0 = 10*tau; // ramp duration
    double f1 = 0.1; // ramp rate

    ElasticLine(int L_in) : L(L_in), dt(0.1), t(0.0) {
        monomer.resize(L, 0.0);
        fuerza.resize(L, 0.0);
        disorder.resize(L * L, 0.0);

        std::srand(0); // for reproducibility
        for (int i = 0; i < L; ++i) {
            for (int j = 0; j < L; ++j) {
                int n = i * L + j;
                disorder[n] = ((double)rand() / RAND_MAX * 2.0 - 1.0) * DisAmp;
            }
        }
    }

    // ramp and pulsated field
    double field() {
        int tmax = 3*tau0;
        if (t < tmax) {
            return f1 * t/ tmax; // ramp
        } else {
            return f0 * ((std::cos(t * 2*M_PI/tau) > 0) ? 1 : 0); // pulses
        }
    }

    void fuerzas() {
        for (int i = 0; i < L; ++i) {
            int ip1 = (i + 1) % L;
            int im1 = (i - 1 + L) % L;

            fuerza[i] = monomer[ip1] + monomer[im1] - 2 * monomer[i]
                        + field() - k0 * monomer[i];

            double u = monomer[i];
            u = u - std::floor(u / L) * L;
            int j = static_cast<int>(u) % L;

            int n = i * L + j;
            int np1 = i * L + (j + 1) % L;

            double interp_force = disorder[n] + (u - j) * (disorder[np1] - disorder[n]);
            fuerza[i] += interp_force;
        }
    }

    void update(int nrun) {
        for (int n = 0; n < nrun; ++n) {
            fuerzas();
            for (int i = 0; i < L; ++i) {
                monomer[i] += fuerza[i] * dt +
                              std::sqrt(Temp * dt) * (2.0 * rand() / RAND_MAX - 1.0);
            }
            t += dt;
        }
    }

    void reset() {
        std::fill(monomer.begin(), monomer.end(), 0.0);
        for (int i = 0; i < L; ++i) {
            for (int j = 0; j < L; ++j) {
                int n = i * L + j;
                disorder[n] = ((double)rand() / RAND_MAX * 2.0 - 1.0) * DisAmp;
            }
        }
        t = 0.0;
    }

    void print_config(std::ofstream &outputFile){
        for(int i = 0; i < L; ++i) {
            outputFile << monomer[i] << "\n";
        }
        outputFile << "\n" << std::endl;
    }

    void print_monitor(std::ofstream &monitorFile) {
      std::tuple<double, double> tup = cm_displacement();
      double cm = std::get<0>(tup);
      double var = std::get<1>(tup);
      monitorFile << t << " " << cm << " " << var << " " << f0 << "\n";
    }

    std::tuple<double, double> cm_displacement()
    {
        double cm = 0.0;
        for (int i = 0; i < L; ++i) cm += monomer[i];
        cm /= L;
        double var = 0.0;
        for (int i = 0; i < L; ++i) var += (monomer[i] - cm) * (monomer[i] - cm);
        var /= L;
        return std::tuple(cm, var);
    }
};

int main(int argc, char **argv) {

    int L = atoi(argv[1]);
    ElasticLine cuerda(L);

    std::ofstream outputFile("cuerda.dat");
    std::ofstream monitorFile("monitor.dat");
    std::ofstream monitorFile2("monitor2.dat");
    std::ofstream monitorstrobFile("monitorstrob.dat");

    // Simulation parameters
    cuerda.f0 = atof(argv[2]); // 0.4
    cuerda.k0 = atof(argv[3]); // 0.0001
    cuerda.DisAmp = atof(argv[4]); // 1.0
    cuerda.tau = atof(argv[5]); // 100.0
    cuerda.Temp = 0.0; // 0.0
    int nrun = 1;

    // Run simulation
    for (int step = 0; step < 200000; ++step) {

        double f0=cuerda.field();
        cuerda.update(nrun);
        double f1=cuerda.field();

	bool pulse = (f0<0.000001 && f1>0);

        //std::cout << f0 << " " << f1 << std::endl;
        if(pulse && step>0){
            //std::cout << "Step " << step << ", field = " << cuerda.field() << std::endl;
            cuerda.print_config(outputFile);
            cuerda.print_monitor(monitorstrobFile);
        }

        if(step%100==0)
        cuerda.print_monitor(monitorFile);

        if(pulse==1)
        cuerda.print_monitor(monitorFile2);
    }

    // Reset simulation
    // cuerda.reset();
}
