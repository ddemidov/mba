#include <iostream>
#include <fstream>
#include <string>
#include <random>

#include <mpi.h>

#include <mba.hpp>

#include "profiler.hpp"

class chunk {
    public:
        chunk(int size, int rank, int n) {
            int chunk_size = (n + size - 1) / size;
            chunk_start = rank * chunk_size;
            chunk_end   = std::min(n, chunk_start + chunk_size);
        }

        int begin() const { return chunk_start; }
        int end()   const { return chunk_end;   }
        int size()  const { return chunk_end - chunk_start; }

    private:
        int chunk_start;
        int chunk_end;

};

int main(int argc, char *argv[]) {
    const size_t n = argc < 2 ? 1024 * 1024 : std::stoi(argv[1]);

    profiler<> prof;

    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    chunk chunk(size, rank, n);

    prof.tic("generate data");
    std::default_random_engine rng(0);
    std::uniform_real_distribution<double> rnd(0.0, 1.0);

    std::vector< std::array<double,2> > p(n);
    std::generate(p.begin(), p.end(), [&]() { return mba::make_array<double>(rnd(rng), rnd(rng)); });

    std::vector<double> v(n);
    std::transform(p.begin(), p.end(), v.begin(), [](const std::array<double,2> &c) {
            double x = c[0] - 0.5;
            double y = c[1] - 0.5;
            return x * x + y * y;
            });
    prof.toc("generate data");

    std::array<double, 2> xmin = {-0.01, -0.01};
    std::array<double, 2> xmax = { 1.01,  1.01};

    prof.tic("create cloud");
    mba::cloud<2> c(
            MPI_COMM_WORLD, xmin, xmax,
            p.begin() + chunk.begin(), p.begin() + chunk.end(),
            v.begin() + chunk.begin(),
            mba::default_grid(xmin, xmax)
            );
    prof.toc("create cloud");

    if (rank == 0) {
        prof.tic("interpolate");
        std::ofstream dat("c.dat");
        dat << std::scientific;
        for(double y = 0; y <= 1; y += 0.01) {
            for(double x = 0; x <= 1; x += 0.01)
                dat << c(x, y) << " ";
            dat << std::endl;
        }
        prof.toc("interpolate");

        std::cout << prof << std::endl;
    }

    MPI_Finalize();
}
