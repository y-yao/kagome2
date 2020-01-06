#include <fstream>
#include <hps/src/hps.h>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <chrono>

unsigned long long site[100], size, *bits;
int sites, pair[100][2], edges, edge[100][2], defects, defect[100][3];
double alpha = -.333333333333;

int mz2total(unsigned long long b) {
  int mz2, s;

  mz2 = -sites;
  for (s = 0; s < sites; ++s)
    if ((b & site[s]) == 0)
      mz2 += 2;

  return mz2;
}

unsigned long long state(unsigned long long b) {
  unsigned long long s1, s2, s;

  s1 = 0;
  s2 = size - 1;
  s = (s1 + s2) / 2;

  while (bits[s] != b) {
    if (s == s1)
      return s2;

    if (bits[s] > b) {
      s2 = s;
      s = (s1 + s2) / 2;
    } else {
      s1 = s;
      s = (s1 + s2) / 2;
    }
  }

  return s;
}

void setbits() {
  int s;
  unsigned long long b;

  for (s = 0; s < sites; ++s)
    site[s] = 1ULL << s;

  size = 0;
  for (b = 0; b < 1ULL << sites; ++b)
    if (mz2total(b) == 0)
      ++size;

  bits = (unsigned long long *)malloc(size * sizeof(unsigned long long));

  size = 0;
  for (b = 0; b < 1ULL << sites; ++b)
    if (mz2total(b) == 0)
      bits[size++] = b;
}

void dimerize(double *psi, unsigned long long b, int sgn, int d) {
  if (d) {
    dimerize(psi, b, sgn, d - 1);
    dimerize(psi, (b ^ site[pair[d - 1][0]]) ^ site[pair[d - 1][1]], -sgn,
             d - 1);
  } else
    psi[state(b)] = (double)sgn;
}

void hpair(int i, int j, double *psi, double *hpsi) {
  unsigned long long s, b, si, sj;

  si = site[i];
  sj = site[j];

  for (s = 0; s < size; ++s) {
    b = bits[s];

    if (((b & si) == 0) == ((b & sj) == 0)) {
      hpsi[s] += .25 * psi[s];
    } else {
      hpsi[s] += -.25 * psi[s];
      hpsi[s] += .5 * psi[state((b ^ si) ^ sj)];
    }
  }
}

void h(double *psi, double *hpsi) {
  unsigned long long s;
  int e;

  for (s = 0; s < size; ++s)
    hpsi[s] = 0.;

  for (e = 0; e < edges; ++e)
    hpair(edge[e][0], edge[e][1], psi, hpsi);
}

void hdefect(int i, int j, int k, double *psi, double *hpsi) {
  unsigned long long s;

  for (s = 0; s < size; ++s)
    hpsi[s] = 0.;

  hpair(i, j, psi, hpsi);
  hpair(j, k, psi, hpsi);
  hpair(k, i, psi, hpsi);

  for (s = 0; s < size; ++s)
    hpsi[s] = (1. + .75 * alpha) * psi[s] + alpha * hpsi[s];
}

double dot(double *psi1, double *psi2) {
  double d;
  unsigned long long s;

  d = 0.;
  for (s = 0; s < size; ++s)
    d += psi1[s] * psi2[s];

  return d;
}

void normalize(double *psi) {
  double norm;
  unsigned long long s;

  norm = sqrt(dot(psi, psi));

  for (s = 0; s < size; ++s)
    psi[s] /= norm;
}

void projout(double *psi, double *psiout) {
  double d;
  unsigned long long s;

  d = dot(psi, psiout);

  for (s = 0; s < size; ++s)
    psi[s] -= d * psiout[s];
}

int main(int argc, char *argv[]) {
  char *edgefile, *dimerfile, *eigenfile;
  FILE *fp;
  int e, dimers, states, sites2, i, j, d, t, t1, t2, t3, num;
  unsigned long long b, s, size2;
  double **basis, *tmp, *evec;

  if (argc == 4) {
    edgefile = argv[1];
    dimerfile = argv[2];
    num = std::atoi(argv[3]);
  } else {
    fprintf(stderr,
            "expected three arguments: edgefile, dimerfile, n_ev\n");
    return 1;
  }

  fp = fopen(edgefile, "r");

  fscanf(fp, "%d%d", &edges, &sites);

  for (e = 0; e < edges; ++e)
    fscanf(fp, "%d%d", &edge[e][0], &edge[e][1]);

  fclose(fp);

  setbits();

  fp = fopen(dimerfile, "r");

  fscanf(fp, "%d%d%d", &states, &sites2, &defects);

  if (sites2 != sites) {
    fprintf(stderr, "edgefile and dimerfile are inconsistent\n");
    return 1;
  }

  auto start = std::chrono::high_resolution_clock::now();

  basis = (double **)malloc(states * sizeof(double *));
  for (i = 0; i < states; ++i)
    basis[i] = (double *)malloc(size * sizeof(double));

  tmp = (double *)malloc(size * sizeof(double));

  dimers = sites / 2;

  for (i = 0; i < states; ++i) {
    for (s = 0; s < size; ++s)
      basis[i][s] = 0.;

    b = 0;
    for (d = 0; d < dimers; ++d) {
      fscanf(fp, "%d%d", &pair[d][0], &pair[d][1]);
      b = b ^ site[pair[d][0]];
    }

    dimerize(basis[i], b, 1, dimers);

    for (t = 0; t < defects; ++t) {
      fscanf(fp, "%d%d%d", &t1, &t2, &t3);

      hdefect(t1, t2, t3, basis[i], tmp);

      for (s = 0; s < size; ++s)
        basis[i][s] = tmp[s];
    }
  }

  fclose(fp);

  for (i = 0; i < states; ++i) {
    for (j = 0; j < i; ++j)
      projout(basis[i], basis[j]);

    normalize(basis[i]);
  }

  for (i = 0; i < states; ++i) {
    h(basis[i], tmp);

    for (j = 0; j < states; ++j)
      printf("%15.10lf", dot(basis[j], tmp));

    printf("\n");
  }

  printf("\n");

  auto end = std::chrono::high_resolution_clock::now();
  std::cout << "\tCalculating low-energy Hamiltonian takes " << std::chrono::duration<double>(end-start).count() / 3600 << "hrs." << std::endl;

  std::cout << "\nsize = " << size << "\n";
  
  start = std::chrono::high_resolution_clock::now();
  
  evec = (double *)malloc(size * sizeof(double));

  std::ifstream eigenvalues_file("eigenvalues.dat", std::ofstream::binary);
  auto eval = hps::from_stream<std::vector<double>>(eigenvalues_file);

  for (i = 0; i < num; i++) {
    std::ifstream indices_file("indices" + std::to_string(i) + ".dat",
                               std::ofstream::binary);
    auto indices = hps::from_stream<std::vector<long int>>(indices_file);
    std::ifstream values_file("values" + std::to_string(i) + ".dat",
                              std::ofstream::binary);
    auto values = hps::from_stream<std::vector<double>>(values_file);
    for (size_t s = 0; s < size; s++) evec[s] = 0.;
    for (size_t s = 0; s < indices.size(); s++) {
      evec[indices[s]] = values[s];
    }
    normalize(evec);
    for (j = 0; j < states; ++j)
      projout(evec, basis[j]);
    printf("%15.10lf%15.10lf\n", eval[i], 1. - dot(evec, evec));
  }
  
  end = std::chrono::high_resolution_clock::now();
  std::cout << "\tProjection takes " << std::chrono::duration<double>(end-start).count() / 3600 << "hrs." << std::endl;

  return 0;
}
