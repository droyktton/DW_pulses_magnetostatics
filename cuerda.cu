#include <iostream>
#include <fstream>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/transform.h>
#include <thrust/transform_reduce.h>
#include <thrust/fill.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>
#include <curand_kernel.h>
#include <cmath>
#include <tuple>

__global__ void kernel_fuerzas(
    double* monomer,
    double* fuerza,
    const double* disorder,
    int L,
    double field,
    double k0)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= L) return;

    int ip1 = (i + 1) % L;
    int im1 = (i - 1 + L) % L;

    double f = monomer[ip1] + monomer[im1] - 2 * monomer[i]
               + field - k0 * monomer[i];

    double u = monomer[i];
    u = u - floor(u / L) * L;
    int j = static_cast<int>(u) % L;

    int n = i * L + j;
    int np1 = i * L + (j + 1) % L;

    double interp_force = disorder[n] + (u - j) * (disorder[np1] - disorder[n]);
    fuerza[i] = f + interp_force;
}

__global__ void kernel_update(
    double* monomer,
    const double* fuerza,
    double dt,
    double Temp,
    int L,
    curandState* states)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= L) return;

    double eta = curand_uniform_double(&states[i]) * 2.0 - 1.0;
    monomer[i] += fuerza[i] * dt + sqrt(Temp * dt) * eta;
}

__global__ void init_curand(curandState* states, unsigned long seed, int L)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < L)
        curand_init(seed, i, 0, &states[i]);
}

class ElasticLine {
public:
    double dt = 0.1;
    double t = 0.0;
    int L;

    thrust::device_vector<double> monomer;
    thrust::device_vector<double> fuerza;
    thrust::device_vector<double> disorder;

    double Temp = 0.0;
    double f0 = 0.6;
    double k0 = 0.0001;
    double DisAmp = 1.0;
    double tau = 100.0;
    double tau0;
    double f1 = 0.1;
    //bool pulses_on = false;
    double tmax = 3 * tau0;

    curandState* devStates;

    ElasticLine(int L_in) : L(L_in), tau0(10 * tau),
        monomer(L, 0.0),
        fuerza(L, 0.0),
        disorder(L * L)
    {
        // Initialize disorder on host
        thrust::host_vector<double> h_disorder(L * L);
        srand(0);
        for (int i = 0; i < L * L; ++i)
            h_disorder[i] = ((double)rand() / RAND_MAX * 2.0 - 1.0) * DisAmp;
        disorder = h_disorder;

        cudaMalloc(&devStates, L * sizeof(curandState));
        init_curand<<<(L+255)/256, 256>>>(devStates, 1234, L);

        //pulses_on = false;
        tmax = 3 * tau0;

    }

    ~ElasticLine() {
        cudaFree(devStates);
    }

    double field() {
        // tmax = 3 * tau0;
        // pulses_on = (t > tmax);
        if (t < tmax)
            return f1 * t / tmax;
        else{
            //pulses_on = true;
            return f0 * ((cos(t * 2 * M_PI / tau) > 0) ? 1 : 0);
        }
    }

    void fuerzas() {
        kernel_fuerzas<<<(L+255)/256, 256>>>(
            thrust::raw_pointer_cast(monomer.data()),
            thrust::raw_pointer_cast(fuerza.data()),
            thrust::raw_pointer_cast(disorder.data()),
            L, field(), k0
        );
    }

    void update(int nrun) {
        for (int n = 0; n < nrun; ++n) {
            fuerzas();
            kernel_update<<<(L+255)/256, 256>>>(
                thrust::raw_pointer_cast(monomer.data()),
                thrust::raw_pointer_cast(fuerza.data()),
                dt, Temp, L, devStates
            );
            t += dt;
        }
    }

    void print_config(std::ofstream &outputFile) {
        thrust::host_vector<double> h_monomer = monomer;
        for (int i = 0; i < L; ++i)
            outputFile << h_monomer[i] << "\n";
        outputFile << "\n" << std::endl;
    }

    std::tuple<double, double> cm_displacement() {
        thrust::host_vector<double> h_monomer = monomer;
        double cm = thrust::reduce(h_monomer.begin(), h_monomer.end()) / L;

        double var = thrust::transform_reduce(
            h_monomer.begin(), h_monomer.end(),
            [=] __host__ __device__ (double u) {
                return (u - cm) * (u - cm);
            },
            0.0, thrust::plus<double>()
        ) / L;

        return std::tuple(cm, var);
    }

    void print_monitor(std::ofstream &monitorFile) {
        auto [cm, var] = cm_displacement();
        monitorFile << t << " " << cm << " " << var << " " << f0 << "\n";
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
    cuerda.Temp = atof(argv[6]); // 0.0
    int nrun = 1;

    std::cout << "L = " << L << std::endl;
    std::cout << "f0 = " << cuerda.f0 << std::endl;
    std::cout << "k0 = " << cuerda.k0 << std::endl;
    std::cout << "DisAmp = " << cuerda.DisAmp << std::endl;
    std::cout << "tau = " << cuerda.tau << std::endl;
    std::cout << "Temp = " << cuerda.Temp << std::endl;
    std::cout << "tmax = " << cuerda.tmax << std::endl;
    std::cout << std::endl;

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

        //if(step%100==0)
        //cuerda.print_monitor(monitorFile);

        if(pulse==true)
        cuerda.print_monitor(monitorFile2);
    }

    // Reset simulation
    // cuerda.reset();
}



