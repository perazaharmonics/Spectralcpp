#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include "EigenSpectral.h"

using std::vector;

// Helper to compare vectors
template <typename T>
bool VectorsClose(const vector<T>& a, const vector<T>& b, T tol = 1e-6) {
  if (a.size() != b.size()) return false;
  for (size_t i = 0; i < a.size(); ++i)
    if (std::abs(a[i] - b[i]) > tol) return false;
  return true;
}

// Helper to compare matrices
template <typename T>
bool MatricesClose(const vector<vector<T>>& A, const vector<vector<T>>& B, T tol = 1e-6) {
  if (A.size() != B.size()) return false;
  for (size_t i = 0; i < A.size(); ++i)
    if (!VectorsClose(A[i], B[i], tol)) return false;
  return true;
}

int main() {
  // --- Test Setup: Two synthetic signal sources with sinusoidal patterns ---
  const int N = 8; // Number of samples
  const int signalSources=2; // Number of signal sources
  std::vector<std::vector<double>> signal(signalSources, std::vector<double>(N));
  signal = {
    {1, 0.707, 0, -0.707, -1, -0.707, 0, 0.707}, // sin wave
    {0, 0.707, 1, 0.707, 0, -0.707, -1, -0.707}  // cos wave
  };

  EigenSpectral<double> spectral(signal);

  // --- Test 1: Covariance Matrix ---
  auto cov = spectral.GetCovarianceMatrix();
  std::cout << "[TEST] Covariance Matrix:\n";
  spectral.PrintMatrix(cov);

  assert(cov.size() == signalSources && cov[0].size() == signalSources);

  // --- Test 2: Eigen Decomposition ---
  auto [eigenvalues, eigenvectors] = spectral.GetEigenDecomposition();
  std::cout << "[TEST] Eigenvalues:\n";
  spectral.PrintVector(eigenvalues);
  std::cout << "[TEST] Eigenvectors:\n";
  spectral.PrintMatrix(eigenvectors);

  // --- Test 3: Determinant and Inverse ---
  double det = spectral.ComputeDeterminant(cov);
  auto inverse = spectral.InverseMatrix(cov);
  std::cout << "[TEST] Determinant: " << det << std::endl;
  std::cout << "[TEST] Inverse Matrix:\n";
  spectral.PrintMatrix(inverse);

  // Sanity: Multiply inverse * original = Identity?
  auto identity = spectral.MultiplyMatrices(cov, inverse);
  for (size_t i = 0; i < identity.size(); ++i)
    for (size_t j = 0; j < identity[i].size(); ++j)
      assert(std::abs(identity[i][j] - (i == j ? 1.0 : 0.0)) < 1e-3);

  // --- Test 4: PCA ---
  auto pca = spectral.PCA(1);
  std::cout << "[TEST] PCA with 1 component:\n";
  spectral.PrintMatrix(pca);
  assert(pca.size() == 1);
  assert(pca[0].size() == signal[0].size());

  // --- Test 5: Phase Unwrapping ---
  auto unwrapped = spectral.UnwrapPhase();
  std::cout << "[TEST] Unwrapped Phase:\n";
  spectral.PrintVector(unwrapped);
  assert(unwrapped.size() == signal[0].size());

  // --- Test 6: MUSIC Frequency Estimation ---
  auto music = spectral.MUSICPhaseEstimation(1, 0.1);
  std::cout << "[TEST] MUSIC Power Spectrum:\n";
  spectral.PrintVector(music);
  assert(music.size() == static_cast<size_t>(spectral.GetMaxIterations()));

  // --- Test 7: ESPRIT Frequency Estimation ---
  auto esprit = spectral.ESPRITPhaseEstimation(1);
  std::cout << "[TEST] ESPRIT Estimated Phases:\n";
  spectral.PrintVector(esprit);
  if (esprit.empty()) 
    std::cerr << "[WARNING] ESPRIT returned empty result.\n";
  else 
    spectral.PrintVector(esprit);
  std::cout << "[SUCCESS] All tests passed.\n";
  return 0;
}
