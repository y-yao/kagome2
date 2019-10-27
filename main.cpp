#include <fstream>
#include <iostream>
#include <iterator>
#include <chrono>
#include "matrix.h"
#include <Spectra/SymEigsSolver.h>
#include <hps/src/hps.h>

using namespace Spectra;

int main(int argc, char *argv[]) {
  if (argc != 5) {
    std::cout << "Usage: [executable_name] [datafile_name] [dim_of_H] [nev] [ncv]\n";
    std::exit(EXIT_FAILURE);
  }

  FILE* pFile;
  pFile = fopen("output", "w");

  fprintf(pFile, "executable: %s\tdata_file: %s\tdim_of_H: %s\tnev: %s\tncv: %s\n", argv[0], argv[1], argv[2], argv[3], argv[4]);
  fflush(stdout);

  // Load Matrix
  size_t H_dim = std::stoull(argv[2]);
  HamiltonianOp H(H_dim);
  H.load_matrix(argv[1]);

  // Diagonalization
  int nev = std::atoi(argv[3]); // Number of eigenvalues requested; require 1 <= nev <= N-1
  int ncv = std::atoi(argv[4]); // larger ncv <-> faster convergence <-> more costly iters;
               // require nev < ncv <= N; recommend ncv >= 2*nev

  HamiltonianOp op(H);
  SymEigsSolver<double, SMALLEST_ALGE, HamiltonianOp> eigs(&op, nev, ncv);
  
  auto start = std::chrono::high_resolution_clock::now();
  eigs.init();
  int nconv = eigs.compute(100, 1e-8);
  auto end = std::chrono::high_resolution_clock::now();
  std::cout << "\tDiagonalization takes " << std::chrono::duration<double>(end-start).count() / 3600 << "hrs." << std::endl;

  H.clear_matrix();

  Eigen::IOFormat Fmt(8);
  if (eigs.info() == SUCCESSFUL) {
    Eigen::VectorXd evals = eigs.eigenvalues();
    Eigen::MatrixXd evecs = eigs.eigenvectors();
    fprintf(pFile, "\nConverged eigenvalues:\n");
    for (int i = nev - 1; i >= 0; i--) { // ground state before excited states
        fprintf(pFile, "%11.7f  ", evals(i));
    }
    fprintf(pFile, "\n");
    for (int i = nev - 1; i >= 0; i--) {
      std::vector<long int> indices;
      std::vector<double> values;
      for (long int j = 0; j < H_dim; j++) {
        if (abs(evecs(j, i)) > 1e-4) {
          indices.push_back(j);
          values.push_back(evecs(j, i));
        }
      }
      int result_file_index = nev - 1 - i; // 0: ground state
      std::ofstream indices_file("indices"+std::to_string(result_file_index)+".dat", std::ofstream::binary);
      hps::to_stream(indices, indices_file);
      std::ofstream values_file("values"+std::to_string(result_file_index)+".dat", std::ofstream::binary);
      hps::to_stream(values, values_file);
    }
  }
  fclose(pFile);

  return 0;
}
