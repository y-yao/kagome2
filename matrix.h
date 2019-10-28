#include <iostream>
#include <fstream>
#include <iterator>
#include <chrono>
#include <Spectra/SymEigsSolver.h>
//#include "sparse_vec.h"

using namespace Spectra;

class SparseVec {
public:
  std::vector<size_t> inds;
  double diag_entry;
};

class HamiltonianOp {
public:
  HamiltonianOp(size_t H_dim_) {
    H_dim = H_dim_;
    H.resize(H_dim);
  };

  void load_matrix_old(std::string Hamiltonian_file_name) {
    std::ifstream file(Hamiltonian_file_name); // input file name

    auto t_start = std::chrono::high_resolution_clock::now();

    unsigned n_nonzero;
    size_t row;
  
    std::string line;
    while (std::getline(file, line)) {
      std::istringstream iss(line);
      std::vector<std::string> results(std::istream_iterator<std::string>{iss},
                                       std::istream_iterator<std::string>());
      n_nonzero = stoul(results[0]);
      row = stoull(results[1]);
      H[row].inds.push_back(stoull(results[1]));
      H[row].diag_entry = stod(results[1 + n_nonzero]);
      for (unsigned i = 2; i <= n_nonzero; i++) {
        size_t ind = stoull(results[i]);
        if (ind < row) continue; // only store upper triangle
        H[row].inds.push_back(ind);
      }
      H[row].inds.shrink_to_fit();
    }
    auto t_end = std::chrono::high_resolution_clock::now();
    std::cout << "\tReading Hamiltonian takes " << std::chrono::duration<double>(t_end-t_start).count() << "s." << std::endl;
  
    if (row+1 != H_dim) {
      std::cout<<"\nH_dim given is wrong!"<<std::endl;
      std::exit(0);
    }
  }

  void load_matrix(std::string Hamiltonian_file_name) {
    std::ifstream file(Hamiltonian_file_name); // input file name

    auto t_start = std::chrono::high_resolution_clock::now();

    unsigned n_nonzero;
    size_t row;
  
    std::string line;
    while (std::getline(file, line)) {
      std::istringstream iss(line);
      std::vector<std::string> results(std::istream_iterator<std::string>{iss},
                                       std::istream_iterator<std::string>());
      n_nonzero = stoul(results[0]);
      row = stoull(results[1]);
      H[row].inds.resize(n_nonzero);
      H[row].diag_entry = stod(results[1 + n_nonzero]);
      H[row].inds[0] = row;
      for (unsigned i = 2; i <= n_nonzero; i++) {
        H[row].inds[i - 1] = stoull(results[i]);
      }
    }
    auto t_end = std::chrono::high_resolution_clock::now();
    std::cout << "\tReading Hamiltonian takes " << std::chrono::duration<double>(t_end-t_start).count() << "s." << std::endl;
  
    if (row+1 != H_dim) {
      std::cout<<"\nH_dim given is wrong!"<<std::endl;
      std::exit(0);
    }
  }

  size_t rows() {return H_dim;}

  size_t cols() {return H_dim;}

//    void perform_op0(const double* x_in, double* y_out) {
//      for (size_t i = 0; i < rows(); i++) {
//      y_out[i] = 0.;
//      }
//      for (size_t i = 0; i < rows(); i++) {
//        const auto& row_vals = H[i].vals;
//        const auto& row_inds = H[i].inds;
//      for (size_t j_id = 0; j_id < row_vals.size(); j_id++) {
//        size_t j = row_inds[j_id];    
//          y_out[i] += x_in[j] * row_vals[j_id];
//      }
//      }
//    }

  void perform_op(const double* x_in, double* y_out) {
    auto t_start = std::chrono::high_resolution_clock::now();
#pragma omp parallel for schedule(static, 1)  
    for (size_t i = 0; i < rows(); i++) {
      y_out[i] = 0.;
    }
#pragma omp parallel for schedule(static, 1)
    for (size_t i = 0; i < rows(); i++) {
      const auto& row_inds = H[i].inds;
      double diff_i = H[i].diag_entry * x_in[i];
      for (size_t j_id = 1; j_id < row_inds.size(); j_id++) {
        size_t j = row_inds[j_id];
        diff_i += 0.5 * x_in[j];
        if (i != j) {
          double diff_j = 0.5 * x_in[i];
#pragma omp atomic          
          y_out[j] += diff_j;
        }
      }
#pragma omp atomic          
      y_out[i] += diff_i;
    }
    auto t_end = std::chrono::high_resolution_clock::now();
    std::cout << "\tMat_vec multiplication takes " << std::chrono::duration<double>(t_end-t_start).count() << "s." << std::endl;
  }

  void clear_matrix() {
    H.clear();
    H.shrink_to_fit();
  }

private:
  size_t H_dim;

  std::vector<SparseVec> H;

};
