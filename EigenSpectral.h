/**
 * Filename:
 *   EigenSignalProcessing.h
 * 
 * Description:
 *   This file contains spectral methods for processing signals using Eigenvectors.
 *  Eigenvectors and Eigenvalues are used in Signal Processing for their ability to
 *  model the underlying signal structure, smooth out noisy phase transitions, and
 *  detect phase discontinuities. This class contains applications of: Eigenvectors
 *  for:
 * 1. Eigenvector-Based Smoothing Phase Unwrapping. By projecting the phase 
 *  information onto the dominant eigenvectors (i.e. those associated with the largest
 *  eigenvalues), the phase information can be smoothed and phase discontinuities can
 *  be detected. Thus reducing noise and smoothing out phase transitions, by 
 *  extracting a more continuous phase function before unwrapping.
 * 
 * 2. Principal Component Analysis (PCA) for Phase Regularization. Eigenvalue
 * decomposition of the phase gradient matrix can be used to determine the primary
 * directions of phase changes. The dominant eigenvectors (corresponding to the
 * largest eigenvalues) represent the most significant phase transitions, which
 * can be used to reconstruct a continuous phase estimate. By filtering out lower
 * eigenvalue components (i.e. those associated with noise), noise-induced
 * phase discontinuities can be reduced before applying unwrapping
 * algorithms like the Goldstein branch-cut method, or minimum Lp-norm unwrapping.
 * 
 * 3. Multiple Signal Classification (MUSIC) for Phase Unwrapping. The MUSIC
 * algorithm can be used to estimate the number and locations of phase discontinuities
 * by projecting phase vectors onto the noise subspace (spanned by the eigenvectors
 * that correspond to the smallest eigenvalues), sharp phase transitions and harmonic
 * components can be identified. This method is particularly useful in multi-component
 * signals with close spectral content. In array signal processing, MUSIC is used to
 * estimate the direction of arrival of signals in the presence of noise, which
 * inherently includes phase information of arriving signals.
 * 
 * 4. Estimation of Signal Parameters via Rotational Invariance Techniques (ESPRIT)
 * for Phase Unwrapping. The ESPRIT algorithm provides a high-resolution approach to
 * estimating phase shifts and frequency components by exploiting the rotational
 * invariance within the signal subspace. The signal is decomposed into overlapping segments
 * and the phase difference between successive segments is computed using eigenvectors of
 * the signal covariance matrix. ESPRIT is an effictient algorithm because it does not
 * require spectral search, making it ideal for real-time DSP or resource constrained phase tracking
 * applications.
 * 
 * This class provides efficient C++ implementations of these techniques using fundamental
 * linear algebra routines.
 */

#ifndef EIGENSPECTRAL_H
#define EIGENSPECTRAL_H

#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>
#include <numeric>

template <typename T>
class EigenSpectral
{
  public:
    // Constructors
    EigenSpectral (const std::vector<std::vector<T>> &signal);
    virtual ~EigenSpectral () = default;

    // Eigen Decomposition using Power Iteration
    void ComputeEigenDecomposition (void);

    // Music algorithm for phase estimation.
    std::vector<T> MUSICPhaseEstimation (int signalSources, double frequencySpacing);

    // ESPRIT algorithm for phase estimation.
    std::vector<T> ESPRITPhaseEstimation (int signalSources);

    // PCA for noise reduction and phase regularization.
    std::vector<std::vector<T>> PCA(int retainedComponents);

    // Graph-based Phase Unwrapping using eigenvectors.
    std::vector<T> UnwrapPhase (const std::vector<T>& phaseData) const;

    // Getter functions
    const inline std::vector<std::vector<T>>& GetSignalMatrix (void) const {return signalMatrix;}
    const inline std::vector<std::vector<T>>& GetCovarianceMatrix (void) const {return covarianceMatrix;}
    const inline std::vector<T>& GetEigenvalues (void) const {return eigenvalues;}
    const inline std::vector<std::vector<T>>& GetEigenvectors (void) const {return eigenvectors;}
    const inline std::vector<T>& GetPhaseUnwrapped (void) const {return phaseUnwrapped;}
    const inline std::pair<std::vector<T>, std::vector<std::vector<T>>> & GetEigenDecomposition(void) const 
    {
      return eigenDecomposition;
    }
    const inline int GetMaxIterations(void) const {return maxIterations;}
    const inline double GetTolerance(void) const {return tolerance;}
    
      
    // Print functions
    void PrintMatrix (const std::vector<std::vector<T>> &matrix);
    void PrintVector (const std::vector<T> &vector);
    void PrintEigenvectors (void);
    void PrintEigenvalues (void);
    void PrintEigenDecomposition (void);
    void PrintCovarianceMatrix (void);
    void PrintPhaseUnwrapped (void);

  private:
    std::vector<std::vector<T>>                           signalMatrix;
    std::vector<std::vector<T>>                           eigenvectors;
    std::vector<T>                                        eigenvalues;
    std::vector<std::vector<T>>                           orthonormalEigenvectors;
    std::vector<std::vector<T>>                           covarianceMatrix;
    std::pair<std::vector<T>, std::vector<std::vector<T>>> eigenDecomposition;
    std::vector<T>                                        normalizedVector;
    std::vector<T>                                        phaseEstimate;
    std::vector<T>                                        phaseUnwrapped;
    std::vector<T>                                        noiseSubspace;
    std::vector<T>                                        signalSubspace;
    int                                                   signalLength{0};
    int                                                   signalDimension{0};
    int                                                   signalSources{0};
    int                                                   retainedComponents{0};
    int                                                   maxIterations{1000};
    T                                                     frequencySpacing=static_cast<T>(0);
    double                                                tolerance{1e-6};
  
      // Calculate the covariance matrix of the input signal.
    void ComputeCovarianceMatrix (void);

    // Find the largest off diagonal element of the matrix.
    void FindLargestOffDiagonalElement (const std::vector<std::vector<T>> &matrix, int &p, int &q);

    // Power Iteration for Eigen Decomposition (approximate eigenvalues)
    std::pair<std::vector<T>,std::vector<std::vector<T>>> PowerIteration (int maxIterations=1000,double tolerance=1e-6);
  public:
    // Matrix operation helper functions
    std::vector<std::vector<T>> MultiplyMatrices (const std::vector<std::vector<T>> &A, const std::vector<std::vector<T>> &B);
    std::vector<std::vector<T>> ScaleMatrix (const std::vector<std::vector<T>> &A, T scalar);
    std::vector<T>              MultiplyMatrixVector (const std::vector<std::vector<T>> &A, const std::vector<T> &B);
    std::vector<std::vector<T>> InverseMatrix (const std::vector<std::vector<T>> &A);
    std::vector<std::vector<T>> TransposeMatrix (const std::vector<std::vector<T>> &A);
    std::vector<std::vector<T>> AddMatrices (const std::vector<std::vector<T>> &A, const std::vector<std::vector<T>> &B);
    std::vector<std::vector<T>> SubtractMatrices (const std::vector<std::vector<T>> &A, const std::vector<std::vector<T>> &B);
    T                           VectorNorm (const std::vector<T> &v);
    T                           MatrixNorm(const std::vector<std::vector<T>> &A);
    T                           VectorDotProduct (const std::vector<T> &v1, const std::vector<T> &v2);
    T                           ComputeDeterminant (const std::vector<std::vector<T>> &A);
    std::vector<std::vector<T>> ComputeAdjoint (const std::vector<std::vector<T>> &A);
    std::vector<T>              NormalizeVector (const std::vector<T> &v);
    std::vector<std::vector<T>> NormalizeMatrix(const std::vector<std::vector<T>> &A); 
    void                        Orthonormalize(std::vector<std::vector<T>> &vectors);
    std::pair<std::vector<T>, std::vector<std::vector<T>>> SolveGeneralizedEigenvalueProblem (const std::vector<std::vector<T>> &A, const std::vector<std::vector<T>> &B);

};


template <typename T>
// Constructor
EigenSpectral<T>::EigenSpectral(
    const std::vector<std::vector<T>> &s)// The input signal.
 : signalMatrix(s)                      // Initialize the signal matrix.
{                                       // ------- Constructor -----------
  signalLength=signalMatrix[0].size();  // Set the signal length.
  signalDimension=signalMatrix.size();  // Set the signal dimension.
  signalSources=signalDimension;        // Set the signal sources.
  frequencySpacing=1;                   // Set the frequency spacing.
  retainedComponents=signalSources;     // Set the retained components.
  eigenvectors.resize(signalSources, std::vector<T>(signalSources, 0.0));
  ComputeCovarianceMatrix();            // Compute the covariance matrix.
  ComputeEigenDecomposition();                 // Compute the eigen decomposition.
}                                       // ------- Constructor -----------

// Compute the covariance matrix of the input signal or pairwise covariance between
// the columns of the signal matrix or signal sources.
// Where Covariance(x)=E[(x-E[x])(x-E[x])^T]
// Where E[x] is the expected value of x.
// Where x is the signal matrix.
// Where X^T is the transpose of X.
// Where E[x] is the mean of x or Expected value of x.
// Where E[x]=1/N*sum(x_i) where i=1 to N.
// Where N is the number of samples.
// Where x_i is the ith sample.
template <typename T>
void EigenSpectral<T>::ComputeCovarianceMatrix (void)
{                                       // ----- ComputeCovarianceMatrix -----
  if (signalMatrix.empty() || signalMatrix[0].empty()) {
    std::cerr << "[ERROR] Empty signal matrix!\n";
    return;
  }  
  int rows = signalMatrix.size();       // Number of rows in the signal matrix.
  int cols = signalMatrix[0].size();    // Number of columns in the signal matrix.
  auto s=signalMatrix;                  // Our signal matrix.
  // Compute the mean of each row (signal)
  std::vector<T> means(rows,0.0);       // Where to store the mean of each signal.
  for (int i=0;i<rows;i++)              // For each signal...
  {                                     // Compute the mean.
    T sum{0.0};                         // Where to store the sum of the elements.
    for (int k=0;k<cols;k++)            // For the amount of samples in the signal.
    {                                   // Begin to collect the sum of the samples.
      sum+=s[i][k];                     // Sum the samples (cols) per signal (rows)
    }                                   // Done aggregating the samples.
    means[i]=sum/cols;                  // Calculate the mean.
  }                                     // Done calculating the mean.
  covarianceMatrix.resize(rows, std::vector<T>(rows, 0.0)); // signals*signals
  // --------------------------------- //
  // Compute the covariance matrix. We must subtract each row's mean before
  // taking the dos product. We also loop twice over the rows because each entry
  // in the convariance matrix corresponds to a pair of signals (row i vs. row j)
  // while the innermost loops sums over the samples (columns).
  // ----------------------------------- //
  for (int i=0;i<rows;i++)               // For the amount of signal sources (rows)
  {                                      // Compute the covariance matrix.
    for (int j=0;j<rows;j++)             // For the amount of signal sources.
    {                                    // Compute pairwise conv between col(i) and col(j)
      T sum{0.0};                        // Initialize the sum.
      for (int k=0;k<cols;k++)           // For the amount of samples (cols) in the matrix.
      {                                  // Compute sum of products el elements from col(i) and col(j)
        T centered_i=s[i][k]-mean;       // Substract the mean from this signal's samples
        T centered_j=s[j][k]-mean;       // Substract the mean from this signal's samples
        sum+=centered_i*centered_j;      // Summate the products.
      }                                  // End of the loop over the rows.
      // Sample covariance uses (cols -1)
      covarianceMatrix[i][j]=sum/(cols-1); // Compute the covariance matrix.
    }                                    // End of the loop over the columns.
  }                                      // End of the loop over the columns.
}                                        // ----- ComputeCovarianceMatrix -----



// --------------------------------------------------------------------------------
// This is actually much harder than I originally thought. I need to
// make a separate class for the Generalized Eigenvalue Problem. This class would
// and should take tha pattern of a dispatcher, where different algorithms are applied
// depending on the conditions defined by the type of input matrix to ensure numerical
// stability and flexibility to handle all cases. Thus yielding two
// new classes. One for Linear Algebraic methods, and another just for the Generalized
// Eigenvalue problem. Finally, the third class. EigenSpectral. Which also needs the 
// ---------------------------------------------------------------------------------
// Compute the eigen decomposition of the covariance matrix using the power iteration
// method, followed by sorting the eigenpairs by descending eigenvalue magnitude.
// The power iteration method is an iterative approach for computing the dominant
// eigenvectors and eigenvalues of a matrix. It starts with a random vector and repeatedly
// multiplies it by the matrix, normalizing after each multiplication. This causes the
// vector to converge toward the direction of the dominant eigenvector.
// The eigenvalue corresponding to this eigenvector is estimated using the Rayleigh
// quotient: λ = (v^T * A * v), where A is the matrix and v is the normalized eigenvector.
// After computing each eigenpair, the matrix is deflated by subtracting λ * (v * v^T)
// to remove the contribution of the found eigenvector and allow convergence toward
// the next dominant one.
// Once all eigenpairs are found, they are sorted in descending order by eigenvalue
// magnitude, ensuring that the first eigenvectors represent the most significant signal
// components. This ordering is essential for subspace methods such as PCA, MUSIC,
// and ESPRIT, where separation of signal and noise relies on eigenvalue dominance.
// Compute the eigen decomposition of the covariance matrix using the power iteration
// method. This method estimates eigenvectors and eigenvalues one at a time by iteratively
// multiplying a matrix by a vector until it converges in direction. The associated eigenvalue
// is found using the Rayleigh quotient. After extracting the dominant eigenvector and eigenvalue,
// their contribution is subtracted (deflated) from the covariance matrix so the next dominant
// pair can be extracted. This process is repeated until all eigenpairs are computed. The method
// assumes the matrix is symmetric and positive semi-definite.


/*template <typename T>
void EigenSpectral<T>::ComputeEigenDecomposition(void)
{                                       // ------- ComputeEigenDecomposition -----------
  auto [eVals, eVecs] = PowerIteration(maxIterations, tolerance); // Use PowerIteration
  eigenDecomposition = std::make_pair(eVals, eVecs); // Store the eigendecomposition
  eigenvalues = eVals;                  // Save the eigenvalues
  eigenvectors = eVecs;                 // Save the eigenvectors
  orthonormalEigenvectors = eVecs;      // Assume orthonormal for now
}                                       // ------- ComputeEigenDecomposition -----------

// Perform power iteration to compute all eigenpairs of the covariance matrix.
// This method estimates one eigenpair at a time, subtracting each contribution
// from the matrix before finding the next (deflation).
template <typename T>
std::pair<std::vector<T>, std::vector<std::vector<T>>>
EigenSpectral<T>::PowerIteration(int maxIterations, double tolerance)
{                                       // ------- PowerIteration -----------
  int size = covarianceMatrix.size();   // Matrix is square: size x size
  std::vector<T> eigenvalues(size, 0.0); // Store eigenvalues
  std::vector<std::vector<T>> eigenvectors(size, std::vector<T>(size, 0.0)); // Store eigenvectors
  std::vector<std::vector<T>> A = covarianceMatrix; // Copy of matrix to deflate

  for (int i = 0; i < size; ++i)        // For each eigenpair to compute
  {
    std::vector<T> b(size, 0.0);        // Random starting vector
    for (int j = 0; j < size; ++j)
      b[j] = static_cast<T>(rand()) / RAND_MAX;

    b = NormalizeVector(b);             // Normalize initial guess
    std::vector<T> prev_b(size, 0.0);   // Store previous guess

    for (int iter = 0; iter < maxIterations; ++iter)
    {
      prev_b = b;                       // Save previous vector
      b = MultiplyMatrixVector(A, b);  // Multiply by matrix
      b = NormalizeVector(b);          // Normalize result

      T diff = 0.0;                     // Check convergence
      for (int k = 0; k < size; ++k)
        diff += std::pow(b[k] - prev_b[k], 2);

      if (std::sqrt(diff) < tolerance) break; // Converged
    }

    std::vector<T> Ab = MultiplyMatrixVector(A, b); // A * b
    T lambda = VectorDotProduct(b, Ab);  // Rayleigh quotient

    eigenvalues[i] = lambda;            // Store eigenvalue
    eigenvectors[i] = b;                // Store eigenvector

    // Deflate the matrix A = A - λ * (v * v^T)
    std::vector<std::vector<T>> outer(size, std::vector<T>(size));
    for (int r = 0; r < size; ++r)
      for (int c = 0; c < size; ++c)
        outer[r][c] = lambda * b[r] * b[c];

    A = SubtractMatrices(A, outer);     // Update matrix for next iteration
  }

  return std::make_pair(eigenvalues, eigenvectors); // Return all eigenpairs
}                                       // ------- PowerIteration -----------
*/

// MUSIC algorithm for phase estimation. Used to estimate the frequencies
// of multiple signals in the presence of noise. Using the eigen-decomposition
// of the covariance matrix of the signal. By separating the signal subspace
// spanned by the eigenvector corresponding to the largest eigenvalues, from
// the noise subspace spanned by the eigenvectors corresponding to the smallest
// eigenvalues, the MUSIC algorithm can be used to estimate the number of signal
// sources and their frequencies.
// The algorithm works by projecting the signal onto the noise subspace and
// computing the power spectrum of the projection. The peaks in the power spectrum
// correspond to the frequencies of the signal sources.
// The eigenvalues corresponding to the noise subspace are typically smaller
// than those corresponding to the signal subspace. By taking the reciprocal of
// the eigenvalues, the algorithm emphasizes the noise subspace and suppresses
// the signal subspace.
// The sum of reciprocals of the noise subspace eigenvalues is used to estimate
// the power spectrum of the signal. This sum is inversely related to the power
// of the signal at a given frequency.
// The power spectrum at each frequency is given by the reciprocal of the sum of 
// the reciprocals of the noise subspace eigenvalues. Thus, the peaks in the 
// power spectrum correspond to the frequencies of where the signal is strong, as
// the contribution from the noise space is minimized.
template <typename T>
std::vector<T> EigenSpectral<T>::MUSICPhaseEstimation(
  int signalSources,                   // Number of signal sources
  double frequencySpacing)             // Frequency spacing
{                                      // ------- MUSICPhaseEstimation -----------
  auto [eigenvalues,eigenvectors]=GetEigenDecomposition(); // Get the eigen decomposition
  const int dim=signalMatrix.size();   // Get the dimension of the signal
  const int N=GetMaxIterations();      // Number of frequency steps in the spectrum
  const int noiseDim=dim - signalSources; // Dimensionality of noise subspace

  std::vector<T> pwrSpectrum(N, 0.0);  // Initialize the power spectrum

  // Sanity check: return flat spectrum if invalid setup
  if (signalSources >= dim || eigenvectors.empty()) 
  {
    std::cerr << "[ERROR] MUSIC: signalSources too large or empty eigenvectors.\n";
    return pwrSpectrum;                // Return all-1s power spectrum
  }

  // Build noise subspace matrix E_n from eigenvectors[signalSources..]
  std::vector<std::vector<T>> En(noiseDim); // Initialize noise subspace
  for (int j = 0; j < noiseDim; ++j)
  {
    En[j] = eigenvectors[signalSources + j]; // Copy the noise eigenvectors
  }

  // Evaluate MUSIC spectrum over N frequency bins
  for (int i = 0; i < N; ++i) 
  {
    T f = i * frequencySpacing;        // Current test frequency
    std::vector<std::complex<T>> a(dim); // Complex steering vector

    // Build the complex steering vector: a(f) = [1, e^{-j2πf}, ..., e^{-j2πf(M-1)}]^T
    for (int k = 0; k < dim; ++k) 
    {
      T phase = -2.0 * M_PI * f * k;   // Compute phase angle
      a[k] = std::complex<T>(std::cos(phase), std::sin(phase)); // e^{-j2πfk}
    }

    // Project steering vector onto noise subspace
    std::complex<T> denom(0.0, 0.0);   // Accumulator for projection magnitude

    for (int j = 0; j < noiseDim; ++j) // For each noise basis vector
    {
      std::complex<T> dot(0.0, 0.0);   // Inner product accumulator

      for (int k = 0; k < dim; ++k)    // Dot product a^H * En[j]
      {
        dot += std::conj(a[k]) * En[j][k]; // Hermitian projection
      }

      denom += dot * std::conj(dot);   // Accumulate |a^H * En[j]|^2
    }

    // Compute MUSIC spectrum value at this frequency
    if (std::abs(denom) > 1e-12)       // Avoid division by zero
      pwrSpectrum[i] = 1.0 / std::abs(denom); // Power is inverse of projection magnitude
    else
      pwrSpectrum[i] = 0.0;            // Otherwise, set to zero
  }                                     // Done looping over frequency bins

  return pwrSpectrum;                   // Return MUSIC power spectrum
}                                       // ------- MUSICPhaseEstimation -----------
template <typename T>
// ESPRIT algorithm for phase estimation. Used to estimate the freq of signal in the
// presence of noise. Useful in applications such as Direction of Arrival (DOA) estimation
// and spectral analysis. The algorithm exploits the rotational invariance property of
// the signal subspace to estimate the phase shifts and frequency components of the signal.
// The signal subspace is spanned by eigenvectors corresponding to the largest eigenvalues
// of the signal's covariance matrix. These eigenvectors represent the directions in which
// the signal energy is concentrated.
// By exploitiong the fact that the signal subspace exhibits a rotational invariance property,
// the algorithm can estimate signal parameters. Rotational invariance means that the signal
// subspace is invariant under certain linear transformations, such as rotations or reflections.
// The ESPRIT algorithm assumes that the signal has a shift-invariant structure, which allows
// decomposition of the signal into two overlapping subspaces. These subspaces are then used
// to form a matrix equation that can be solved to estimate the signal parameters.
std::vector<T> EigenSpectral<T>::ESPRITPhaseEstimation(
  int signalSources)                    // Number of signal sources
{                                       // ------- ESPIRITPhaseEstimation -----------
  auto [eigenvalues,eigenvectors]=GetEigenDecomposition(); // Compute the eigenvalues and eigenvectors.
  // Extract the signal subspace (eigenvectors corresponding to the largest eigenvalues).
  std::vector<std::vector<T>> signalSubspace(signalSources,std::vector<T>(eigenvectors[0].size()));
  for (int i=0;i<signalSources;i++)     // For the # of signal sources.
  {
    signalSubspace[i]=eigenvectors[i];  // Extract the signal subspace.
  }                                     // Done extracting the signal subspace.
  if (signalSources < 1 || signalSubspace[0].size() == 0) {
    std::cerr << "[ERROR] ESPRIT: signalSources too small or empty eigenvectors.\n";
    return {};
  }
  // Only attempt to build Y and Z if signalSources > 1
  if (signalSources == 1) 
  {
    std::cerr << "[INFO] ESPRIT skipped subspace shift (only one source).\n";
    return eigenvectors[0]; // or some trivial estimate
  }
  // Form the submatrices by shifting the signal subspace.
  size_t L = signalSubspace[0].size();
  std::vector<std::vector<T>> Y(signalSources - 1, std::vector<T>(L));
  std::vector<std::vector<T>> Z(signalSources - 1, std::vector<T>(L));
  for (int i=0;i<signalSources-1;++i)   // For the # of signal sources - 1.
  {
    
    Y[i]=signalSubspace[i];             // Form the submatrix Y.
    Z[i]=signalSubspace[i+1];           // Form the submatrix Z.
  }                                     // Done forming the submatrices.
  /// Solve the matrix equation to estimate the phase shifts and frequency components.
  // This involves solving the generalized eigenvalue problem Y^T*Y*alpha=Z^T*Z*alpha*lambda
  // where alpha is the phase shift vector and lambda is the frequency vector.
  auto [estimatedPhases,estimatedFrequencies]=SolveGeneralizedEigenvalueProblem(Y,Z);
  return estimatedPhases;
}

// Solve the Generalized Eigenvalue Problem using the Jacobi method.
// The Jacobi method is an iterative algorithm that diagonalizes a symmetric matrix
// by applying a series of rotations to the matrix. The rotations are chosen to
// eliminate off-diagonal elements, resulting in a diagonal matrix.
// The algorithm works by iteratively applying rotations to the matrix until the
// off-diagonal elements are sufficiently small. At each iteration, the algorithm
// selects the largest off-diagonal element and rotates the corresponding row and
// column to eliminate it. The process is repeated until the off-diagonal elements
// are below a specified tolerance.
template <typename T>
std::pair<std::vector<T>, std::vector<std::vector<T>>> EigenSpectral<T>::SolveGeneralizedEigenvalueProblem(
    const std::vector<std::vector<T>> &A, // Matrix A
    const std::vector<std::vector<T>> &B) // Matrix B
{                                       // ------- SolveGeneralizedEigenvalueProblem -----------
  if (A.size() != B.size() || A[0].size() != B[0].size()) 
  {                                     // If the sizes of the matrices do not match.
    std::cerr << "[ERROR] Mismatched matrix sizes in SolveGeneralizedEigenvalueProblem.\n";
    return {};                          // Return an empty pair.
  }                                     // Done checking matrix sizes.
    int n = A.size();                   // Size of the matrix A.
    std::vector<T> eigenvalues(n, 0.0); // Initialize the eigenvalues.
    std::vector<std::vector<T>> eigenvectors(n, std::vector<T>(n, 0.0)); // Initialize the eigenvectors.
    // Initialize the matrix V as the identity matrix.
    std::vector<std::vector<T>> V(n, std::vector<T>(n, 0.0));
    for (int i = 0; i < n; i++) 
    {                                   // For the size of the matrix.
        V[i][i] = 1.0;                  // Set the diagonal elements to 1.
    }                                   // Done initializing the matrix V.
    // Compute the maximum number of iterations.
    int maxIters = GetMaxIterations();  // Get the maximum number of iterations.
    // Compute the matrix C = B^-1 * A.
    std::vector<std::vector<T>> C = MultiplyMatrices(InverseMatrix(B), A); // Compute the matrix C.
    // Compute the eigenvalues and eigenvectors using the Jacobi method.
    for (int i = 0; i < maxIters; i++) { // For the maximum number of iterations.
        // Find the largest off-diagonal element of the matrix C.
        int p{0};                       // Initialize the indices.
        int q{1};                        // Initialize the indices.
        FindLargestOffDiagonalElement(C, p, q); // Find the largest off-diagonal element.
        // Compute the rotation angle.
        T theta = 0.5 * std::atan2(2.0 * C[p][q], C[p][p] - C[q][q]); // Compute the rotation angle.
        // Compute the rotation matrix.
        std::vector<std::vector<T>> R(n, std::vector<T>(n, 0.0)); // Initialize the rotation matrix.
        for (int j = 0; j < n; j++) {   // For the size of the matrix.
            R[j][j] = 1.0;              // Set the diagonal elements to 1.
        }                               // Done initializing the rotation matrix.
        R[p][p] = R[q][q] = std::cos(theta);
        R[p][q] = -std::sin(theta);
        R[q][p] = std::sin(theta);
        // Apply the rotation to the matrix C.
        C = MultiplyMatrices(TransposeMatrix(R), MultiplyMatrices(C, R));
        // Apply the rotation to the eigenvector matrix V.
        V = MultiplyMatrices(V, R);
    }
    // Extract the eigenvalues from the diagonal of the matrix C.
    for (int i = 0; i < n; i++) {
        eigenvalues[i] = C[i][i];
    }
    return std::make_pair(eigenvalues, V); // Return the eigenvalues and eigenvectors.
}

// Find the largest off-diagonal element of the matrix.
template <typename T>
void EigenSpectral<T>::FindLargestOffDiagonalElement(
    const std::vector<std::vector<T>> &matrix, // Input matrix
    int &p,                             // Row index of the largest off-diagonal element
    int &q)                             // Column index of the largest off-diagonal element
{
    int n = matrix.size();              // Size of the matrix
    T maxElement = 0;                   // Initialize the maximum element
    // Iterate over all off-diagonal elements
    for (int i = 0; i < n; ++i)         // For the size of the matrix (# of rows).
    {                                   // Loop over the matrix's elements.
      for (int j = 0; j < n; ++j)       // For the # of columns in the matrix.
      {                                 // Loop over the columns.
        // If the element is off-diagonal..
        if (i != j && std::abs(matrix[i][j]) > maxElement) // and larger than the current maximum,  
        {                               // Then
          maxElement = std::abs(matrix[i][j]);// Update the maximum element.
          p = i;                        // Update the row index to where we found the maximum.
          q = j;                        // Update the column index to where we found the maximum.
        }                               // Done updating the maximum element.
      }                                 // Done looping over the columns.
  }                                     // Done looping over the rows.
}                                       // Done finding the largest off-diagonal element.

// PCA for noise reduction and phase regularization.
// Principal Component Analysis (PCA) is a technique used to reduce the dimensionality
// of data by finding the principal components that capture the most variance in the data.
// The principal components are the eigenvectors of the covariance matrix of the data.
// By projecting the data onto the principal components, the data can be represented in
// a lower-dimensional space while preserving the most important information.
// In the context of phase estimation, PCA can be used to reduce the noise in the phase
// estimates by filtering out the components associated with noise. By retaining only the
// components corresponding to the largest eigenvalues of the phase gradient matrix, the
// most significant phase transitions can be preserved, while the noise-induced phase
// discontinuities can be reduced.
template <typename T>
std::vector<std::vector<T>> EigenSpectral<T>::PCA(
  int retainedComponents)               // Number of components to retain.
{                                       // ------- PCA ------------------
  auto [eigenvalues,eigenvectors]=GetEigenDecomposition(); // Compute the eigenvalues and eigenvectors.
  std::vector<std::vector<T>> principalComponents; // Initialize the principal components.
  for (int i = 0; i < retainedComponents; ++i) // For the # of retained components.
  {                                     // Loop over the retained components.
    std::vector<T> projected(signalMatrix[0].size(), 0.0);// Initialize the projected vector.
    for (size_t t = 0; t < signalMatrix[0].size(); ++t) 
    {                                   // Loop over the signal matrix columns.
      for (size_t j = 0; j < signalMatrix.size(); ++j) 
      {                                 // Get the projected vector.
        projected[t] += orthonormalEigenvectors[i][j] * signalMatrix[j][t];
      }                                 // Done computing the projected vector.
    }                                   // Done looping over the signal matrix columns.
    principalComponents.push_back(projected);// Store the projected vector.
  }                                     // Done looping over the retained components.
  return principalComponents;           // Return the principal components.
}  

/**
 * @brief Unwraps a sequence of phase angles to produce a continuous phase progression.
 *
 * This function removes 2π ambiguities from a wrapped phase vector by detecting 
 * discontinuities between consecutive phase samples and compensating for them. 
 * It ensures that the difference between successive phase values in the output 
 * does not exceed π radians in magnitude (i.e., no jumps larger than 180 degrees), 
 * by adding or subtracting the appropriate multiple of 2π.
 *
 * The algorithm is based on Itoh's 1D phase unwrapping method [oai_citation_attribution:4‡pmc.ncbi.nlm.nih.gov](https://pmc.ncbi.nlm.nih.gov/articles/PMC3156002/#:~:text=Itoh%E2%80%99s%20method%20%28Itoh%2C%201982%29,Then%20it%20wraps):
 * starting from the first phase value (assumed to be within [-π, π)), it iteratively
 * adds or subtracts 2π whenever a jump > π (or < -π) is encountered, effectively 
 * "unwrapping" the phase. This is extended to handle large jumps that may span 
 * multiple 2π intervals, which can occur with under-sampled signals or abrupt 
 * phase shifts in eigenvector data.
 *
 * **Assumptions:** The input vector `phaseData` contains phase angles (in radians) 
 * that are *wrapped* to a principal range, typically [-π, π) or [0, 2π). The function 
 * does not internally wrap the input, so inputs outside these ranges should be 
 * normalized beforehand. The output will be a continuous phase sequence where any 
 * 2π discontinuities are removed. Note that true phase changes of ±π (180°) 
 * remain, as they are inherently ambiguous (the algorithm chooses not to alter 
 * differences of exactly π to avoid arbitrary shifts).
 *
 * **Edge Cases:** 
 * - If `phaseData` is empty, the result is an empty vector.
 * - If `phaseData` has a single element, it is returned as-is (no unwrapping needed).
 * - Large phase discontinuities (> 2π) between samples are handled by applying multiple 
 *   2π corrections iteratively until the jump is within [-π, π]. This means even if the 
 *   phase suddenly jumps several revolutions, the algorithm will unwrap it to the nearest 
 *   equivalent angle, preserving continuity.
 * - Noisy data: Small oscillations below the π threshold are left untouched. If noise causes 
 *   a jump slightly over the threshold, a 2π correction will be applied (which may or may not 
 *   reflect the true continuous phase – users should pre-filter noise if needed). The algorithm 
 *   uses a threshold of π radians for jumps (with a tiny epsilon margin) to decide unwrapping, 
 *   aligning with the standard MATLAB/NumPy approach [oai_citation_attribution:5‡ncl.ucar.edu](https://www.ncl.ucar.edu/Document/Functions/Contributed/unwrap_phase.shtml#:~:text=This%20function%20mimics%20the%20behavior,jump%20tolerance%20of%20pi%20radians).
 *
 * @param phaseData  Vector of wrapped phase values (in radians) to be unwrapped.
 * @return std::vector<T>  New vector of the same length with unwrapped phase values.
 */
template <typename T>
std::vector<T> EigenSpectral<T>::UnwrapPhase(const std::vector<T>& phaseData) const
{
    std::vector<T> unwrapped;
    unwrapped.reserve(phaseData.size());  // reserve output size for efficiency

    if (phaseData.empty()) {
        return unwrapped;  // return empty result for empty input
    }

    // Constants for 2π and π (cast to T to match vector's value type).
    const T TWO_PI = static_cast<T>(2 * M_PI);
    const T PI     = static_cast<T>(M_PI);
    const T TOL    = static_cast<T>(1e-9);  // small tolerance to avoid floating-point edge issues

    // Start unwrapping from the first value.
    unwrapped.push_back(phaseData[0]);  
    T offset = 0;  
    T prevPhase = phaseData[0];  // last input phase value seen (for computing difference)

    // Iterate through the phase data, unwrapping each subsequent value.
    for (size_t i = 1; i < phaseData.size(); ++i) {
        T currPhase = phaseData[i];
        // Calculate raw phase difference between current and previous sample
        T diff = currPhase - prevPhase;

        // Check for discontinuities greater than π radians.
        if (diff > PI + TOL) {
            // **Positive jump exceeds +π**: we likely wrapped downward (e.g., went from just below +π to just above -π).
            // Compute how many 2π to subtract to reduce the jump. At least one 2π will be subtracted.
            // We use a while-loop to handle cases of multiple wraps (difference much larger than π).
            while (diff > PI + TOL) {
                // Subtract 2π from the current offset to compensate for a wrap-around.
                offset -= TWO_PI;
                diff   -= TWO_PI;
            }
            // After this loop, diff should be <= π (the jump has been reduced by the necessary 2π multiples).
        }
        else if (diff < -PI - TOL) {
            // **Negative jump exceeds -π**: likely wrapped upward (e.g., from just above -π to just below +π).
            while (diff < -PI - TOL) {
                // Add 2π to the offset to compensate for the negative wrap-around.
                offset += TWO_PI;
                diff   += TWO_PI;
            }
            // Now diff should be >= -π.
        }

        // Compute the unwrapped phase for the current sample by applying the accumulated offset.
        // The offset is zero most of the time and changes only when a wrap is detected, ensuring continuity.
        T unwrappedPhase = currPhase + offset;
        unwrapped.push_back(unwrappedPhase);

        // Update the reference for previous phase (leave `prevPhase` as the original wrapped phase of the last sample).
        prevPhase = currPhase;
    }

    return unwrapped;
}

// Matrix operation helper functions
template <typename T>
std::vector<std::vector<T>> EigenSpectral<T>::MultiplyMatrices(
  const std::vector<std::vector<T>> &A, // Matrix A
  const std::vector<std::vector<T>> &B) // Matrix B
{                                       // ------- MultiplyMatrices -----------
  size_t rows=A.size();                 // Number of rows in matrix A.
  size_t cols=B[0].size();              // Number of columns in matrix B.
  std::vector<std::vector<T>> C(rows,std::vector<T>(cols,0.0)); // Initialize the result matrix.
  // Perform the MxN * NxP matrix multiplication.
  for (size_t i=0;i<rows;i++)           // For the number of rows in matrix A.
  {                                     // Loop over the rows of matrix A.
    for (size_t j=0;j<cols;j++)         // For the number of columns in matrix B.
    {                                   // Loop over the columns of matrix B.
      for (size_t k=0;k<B.size();k++)   // For the number of rows in matrix B.
      {                                 // Loop over the columns of matrix A.
        C[i][j]+=A[i][k]*B[k][j];       // Compute the matrix product.
      }                                 // Done computing the matrix product.
    }                                   // Done looping over the columns of matrix B.
  }                                     // Done looping over the rows of matrix A.
  return C;                             // Return the result matrix.
}                                       // ------- MultiplyMatrices -----------

template <typename T>
// Cm = Amxn*Bmxn
std::vector<T> EigenSpectral<T>::MultiplyMatrixVector(
  const std::vector<std::vector<T>> &A, // Matrix A
  const std::vector<T> &v)              // Vector v
{                                       // ------- MultiplyMatrixVector -----------
  // The size of the resultant vector will be the number of rows in matrix A.
  std::vector<T> C(A.size(),0.0);       // Initialize the result vector
  for (size_t i=0;i<A.size();i++)       // For the # of rows in matrix A.
  {                                     // Loop over the elements.
    for (size_t j=0;j<v.size();j++)     // For the # of rows in vector v
    {                                   // Loop over the elements.
      C[i]+=A[i][j]*v[j];               // Compute the matrix-vector product.
    }                                   // Done traversing vector.
  }                                     // Done traversing matrix.
  return C;                             // Return resulting matrix
}                                       // ------- MultiplyMatrixVector -----------

template <typename T>
T EigenSpectral<T>::VectorNorm(
  const std::vector<T> &v)              // Vector v
{                                       // ------- VectorNorm -----------
  T sum{0.0};                           // Initialize the sum.
  for (T val: v) sum+=val*val;          // Compute the sum of squares.
  return std::sqrt(sum);                // Return the square root of the sum.
}                                       // ------- VectorNorm -----------

template <typename T>
std::vector<std::vector<T>> EigenSpectral<T>::InverseMatrix(
  const std::vector<std::vector<T>> &A) // Input matrix A
{                                       // ------- InverseMatrix -----------
  T det = ComputeDeterminant(A);        // Compute the determinant of A.
  if (std::abs(det) < 1e-12)            // Check for singular matrix (non-invertible).
  {                                     // and return a zero matrix - its singular, and cannot be inverted.
    return std::vector<std::vector<T>>(A.size(), std::vector<T>(A[0].size(), 0.0));
  }                                      // Done checking for singular matrix.
  // Compute the adjoint (transpose of cofactor matrix).
  std::vector<std::vector<T>> adjoint = ComputeAdjoint(A);
  // Allocate space for the inverse matrix.
  std::vector<std::vector<T>> inverse(
    A.size(), std::vector<T>(A[0].size(), 0.0));
  // Divide each element of the adjoint by the determinant to get the inverse.
  for (size_t i = 0; i < A.size(); ++i) // For the number of rows
  {                                     // ... and ..
    for (size_t j = 0; j < A[0].size(); ++j)// For the number of cols
    {                                   // 
      inverse[i][j] = adjoint[i][j] / det;// Compute inverse element.
    }                                   // Done computing the inverse element.
  }                                     // Done traversing the matrix.
  return inverse;                       // Return the inverse matrix.
}                                       // ------- InverseMatrix -----------


// The determinant of a matrix is a scalar value that can be computed from the elements
// of the matrix. It is a measure of the "volume" of the matrix, representing the
// product of the eigenvalues of the matrix. The determinant is used in various
// mathematical operations, such as solving systems of linear equations, computing
// the inverse of a matrix, and finding the eigenvalues of a matrix.
// The determinant of a matrix can be computed using various methods, such as cofactor
// expansion, LU decomposition, or the product of the eigenvalues.
// In the context of eigenvalue decomposition, the determinant of a matrix is used to
// compute the characteristic polynomial, which is then used to find the eigenvalues
// of the matrix.
template <typename T>
T EigenSpectral<T>::ComputeDeterminant(
  const std::vector<std::vector<T>> &A) // Matrix A
{                                       // ------- ComputeDeterminant -----------
  size_t n=A.size();                    // # of rows (assuming square matrix) is matrix size
  std::vector<std::vector<T>> U=A;      // Initialize Upper triangular matrix copy.
  T det{1.0};                           // Initialize the determinant.
  for (size_t i=0;i<n;++i)              // For # of rows in input matrix.
  {                                     // Loop over pivot rows.
    if (std::abs(U[i][i])<1e-12)        // Check for 0 pivot (singular matrix)
    {  
      return static_cast<T>(0);         // matrix is singular so det(A)=0.
    }
    for (size_t k=i+1;k<n;++k)          // For row i +1 to row n.
    {                                   // Eliminate entries below the pivot.
      T factor=U[k][i]/U[i][i];         // Compute the elimination factor.
      for (size_t j=i;j<n;++j)          // Update row k from col i to end.
        U[k][j]-=factor*U[i][j];        // Substract scaled pivot from target row.
    }                                   // Done with the LU decomposition.
    det*=U[i][i];                       // Multiply diagonal entries for determinant.
  }                                     // Done computing the determinant.
  return det;                           // Return the determinant.
}                                       // ------- ComputeDeterminant -----------
// Compute the dot product of two vectors.
// The dot product is defined as: sum(a_i * b_i)
// It measures the similarity (projection) between two vectors.
template <typename T>
T EigenSpectral<T>::VectorDotProduct(
  const std::vector<T> &v1,             // Vector 1
  const std::vector<T> &v2)             // Vector 2
{                                       // ------- VectorDotProduct -----------
  T sum{0.0};                           // Initialize the sum
  for (size_t i = 0; i < v1.size(); ++i)// For each element of the vectors
  {
    sum += v1[i] * v2[i];               // Accumulate the product
  }                                     // Done traversing both vectors
  return sum;                           // Return the dot product
}                                       // ------- VectorDotProduct -----------
// The Adjunt matrix is the transpose of the cofactor matrix.
// The cofactor matrix is the matrix of determinants of the minors of the original matrix.
// The minor of a matrix is the determinant of the matrix obtained by removing the ith row
// and jth column of the original matrix.
// The cofactor of an element is the determinant of the minor of the element multiplied by
// (-1)^(i+j), where i and j are the row and column indices of the element.
template <typename T>
std::vector<std::vector<T>> EigenSpectral<T>::ComputeAdjoint(
  const std::vector<std::vector<T>> &A) // Matrix A
{                                       // ------- ComputeAdjoint -----------
  // Compute the cofactor matrix.
  std::vector<std::vector<T>> cofactorMatrix(A.size(),std::vector<T>(A[0].size(),0.0)); // Initialize the cofactor matrix.
  for (size_t i=0;i<A.size();i++)       // For the # of rows in matrix A.
  {                                     // Loop over the rows.
    for (size_t j=0;j<A[0].size();j++)  // For the # of columns in matrix A.
    {                                   // Loop over the columns.
      // Compute the cofactor of the element.
      std::vector<std::vector<T>> cofactor(A.size()-1,std::vector<T>(A[0].size()-1,0.0)); // Initialize the cofactor matrix.
      for (size_t k=0;k<A.size();k++)   // For the # of rows in the cofactor matrix.
      {                                 // Loop over the rows.
        for (size_t l=0;l<A[0].size();l++) // For the # of columns in the cofactor matrix.
        {                               // Loop over the columns.
          if (k<i && l<j) cofactor[k][l]=A[k][l]; // Compute the cofactor.
          else if (k<i && l>j) cofactor[k][l-1]=A[k][l]; // Compute the cofactor.
          else if (k>i && l<j) cofactor[k-1][l]=A[k][l]; // Compute the cofactor.
          else if (k>i && l>j) cofactor[k-1][l-1]=A[k][l]; // Compute the cofactor.
        }                               // Done computing the cofactor.
      }                                 // Done looping over the rows.
      // Compute the determinant of the cofactor matrix.
      cofactorMatrix[i][j]=std::pow(-1,i+j)*ComputeDeterminant(cofactor); // Compute the determinant.
    }                                   // Done looping over the columns.
  }                                     // Done looping over the rows.
  // Compute the adjoint matrix by taking the transpose of the cofactor matrix.
  std::vector<std::vector<T>> adjoint=TransposeMatrix(cofactorMatrix); // Compute the adjoint matrix
  return adjoint;                       // Return the adjoint matrix.
}                                       // ------- ComputeAdjoint -----------

template <typename T>
std::vector<std::vector<T>> EigenSpectral<T>::TransposeMatrix(
  const std::vector<std::vector<T>> &A) // Matrix A
{                                       // ------- TransposeMatrix -----------
  // Initialize the transpose matrix.
  std::vector<std::vector<T>> transpose(A[0].size(),std::vector<T>(A.size(),0.0));
  for (size_t i=0;i<A.size();i++)       // For the # of rows in matrix A.
  {                                     // Loop over the rows.
    for (size_t j=0;j<A[0].size();j++)  // For the # of columns in matrix A.
    {                                   // Loop over the columns.
      transpose[j][i]=A[i][j];          // Compute the transpose.
    }                                   // Done computing the transpose.
  }                                     // Done looping over the rows.
  return transpose;                     // Return the transpose matrix.
}                                       // ------- TransposeMatrix -----------

template <typename T>
std::vector<std::vector<T>> EigenSpectral<T>::AddMatrices(
  const std::vector<std::vector<T>> &A, // Matrix A
  const std::vector<std::vector<T>> &B) // Matrix B
{                                       // ------- AddMatrices -----------
  // Initialize the result matrix.
  std::vector<std::vector<T>> C(A.size(),std::vector<T>(A[0].size(),0.0));
  for (size_t i=0;i<A.size();i++)       // For the # of rows in matrix A.
  {                                     // Loop over the rows.
    for (size_t j=0;j<A[0].size();j++)  // For the # of columns in matrix A.
    {                                   // Loop over the columns.
      C[i][j]=A[i][j]+B[i][j];          // Compute the sum.
    }                                   // Done computing the sum.
  }                                     // Done looping over the rows.
  return C;                             // Return the result matrix.
}                                       // ------- AddMatrices -----------

template <typename T>
std::vector<std::vector<T>> EigenSpectral<T>::SubtractMatrices(
  const std::vector<std::vector<T>> &A, // Matrix A
  const std::vector<std::vector<T>> &B) // Matrix B
{                                       // ------- SubtractMatrices -----------
  // Initialize the result matrix.
  std::vector<std::vector<T>> C(A.size(),std::vector<T>(A[0].size(),0.0));
  for (size_t i=0;i<A.size();i++)       // For the # of rows in matrix A.
  {                                     // Loop over the rows.
    for (size_t j=0;j<A[0].size();j++)  // For the # of columns in matrix A.
    {                                   // Loop over the columns.
      C[i][j]=A[i][j]-B[i][j];          // Compute the difference.
    }                                   // Done computing the difference.
  }                                     // Done looping over the rows.
  return C;                             // Return the result matrix.
}                                       // ------- SubtractMatrices -----------

template <typename T>
std::vector<std::vector<T>> EigenSpectral<T>::ScaleMatrix(
  const std::vector<std::vector<T>> &A, // Matrix A
  T scalar)                            // Scalar value
{                                       // ------- ScaleMatrix -----------
  // Initialize the result matrix.
  std::vector<std::vector<T>> C(A.size(),std::vector<T>(A[0].size(),0.0));
  for (size_t i=0;i<A.size();i++)       // For the # of rows in matrix A.
  {                                     // Loop over the rows.
    for (size_t j=0;j<A[0].size();j++)  // For the # of columns in matrix A.
    {                                   // Loop over the columns.
      C[i][j]=A[i][j]*scalar;           // Compute the scaled matrix.
    }                                   // Done computing the scaled matrix.
  }                                     // Done looping over the rows.
  return C;                             // Return the result matrix.
}                                       // ------- ScaleMatrix -----------

template <typename T>
std::vector<std::vector<T>> EigenSpectral<T>::NormalizeMatrix(
  const std::vector<std::vector<T>> &A) // Matrix A
{                                       // ------- NormalizeMatrix -----------
  // Compute the norm of the matrix.
  T norm=MatrixNorm(A);                 // Compute the norm.
  // Normalize the matrix.
  std::vector<std::vector<T>> C(A.size(),std::vector<T>(A[0].size(),0.0));
  for (size_t i=0;i<A.size();i++)       // For the # of rows in matrix A.
  {                                     // Loop over the rows.
    for (size_t j=0;j<A[0].size();j++)  // For the # of columns in matrix A.
    {                                   // Loop over the columns.
      C[i][j]=A[i][j]/norm;             // Compute the normalized matrix.
    }                                   // Done computing the normalized matrix.
  }                                     // Done looping over the rows.
  return C;                             // Return the normalized matrix.
}                                       // ------- NormalizeMatrix -----------

// Perform Gram-Schmidt orthogonalization on a set of vectors.
template <typename T>
void EigenSpectral<T>::Orthonormalize(std::vector<std::vector<T>> &vectors) {
  for (size_t i = 0; i < vectors.size(); ++i) // For # of vectors
  {                                     // Loop over the vectors.
    for (size_t j = 0; j < i; ++j)      // For columns less than i
    { // Calculate the projection of vector i onto vector j and subtract it.
      T dot = std::inner_product(vectors[i].begin(), vectors[i].end(),
                                 vectors[j].begin(), 0.0);
      for (size_t k = 0; k < vectors[i].size(); ++k) 
      {
        vectors[i][k] -= dot * vectors[j][k];
      }                        // Done subtracting the projection.
    }                          // Done looping over the columns.
    vectors[i] = NormalizeVector(vectors[i]);// Normalize the vector.
  }
}

template <typename T>
T EigenSpectral<T>::MatrixNorm(
  const std::vector<std::vector<T>> &A) // Matrix A
{                                       // ------- MatrixNorm -----------
  T sum{0.0};                           // Initialize the sum.
  for (size_t i=0;i<A.size();i++)       // For the # of rows in matrix A.
  {                                     // Loop over the rows.
    for (size_t j=0;j<A[0].size();j++)  // For the # of columns in matrix A.
    {                                   // Loop over the columns.
      sum+=A[i][j]*A[i][j];             // Compute the sum of squares.
    }                                   // Done computing the sum of squares.
  }                                     // Done looping over the rows.
  return std::sqrt(sum);                // Return the square root of the sum.
}                                       // ------- MatrixNorm -----------

template <typename T>
std::vector<T> EigenSpectral<T>::NormalizeVector(
  const std::vector<T> &v)              // Vector v
{
  T norm=VectorNorm(v);                 // Compute the norm of the vector.
  std::vector<T> normalizedVect(v.size()); // Initialize the normalized vector.
  for (size_t i=0;i<v.size();i++)       // For the length of the vector
  {
    normalizedVect[i]=v[i]/norm;        // Normalize the vector.
  }                                     // Done normalizing vector.
  return normalizedVect;                // Return the normalized vector.
}                                       // ------- NormalizeVector -----------

// Variety of Print functions
template <typename T>
void EigenSpectral<T>::PrintMatrix(
  const std::vector<std::vector<T>> &A) // Matrix A
{                                       // ------- PrintMatrix -----------
  for (size_t i=0;i<A.size();i++)       // For the # of rows in matrix A.
  {                                     // Loop over the rows.
    for (size_t j=0;j<A[0].size();j++)  // For the # of columns in matrix A.
    {                                   // Loop over the columns.
      std::cout<<A[i][j]<<" ";           // Print the element.
    }                                   // Done printing the element.
    std::cout<<std::endl;                // Print a new line.
  }                                     // Done looping over the rows.
}                                       // ------- PrintMatrix -----------

template <typename T>
void EigenSpectral<T>::PrintVector(
  const std::vector<T> &v)              // Vector v
{                                       // ------- PrintVector -----------
  for (size_t i=0;i<v.size();i++)       // For the length of the vector.
  {                                     // Loop over the vector.
    std::cout<<v[i]<<" ";                // Print the element.
  }                                     // Done printing the element.
  std::cout<<std::endl;                  // Print a new line.
}                                       // ------- PrintVector -----------

template <typename T>
void EigenSpectral<T>::PrintEigenvalues(void)
{                                       // ------- PrintEigenvalues -----------
  std::cout<<"Eigenvalues: ";           // Print the eigenvalues.
  for (size_t i=0;i<eigenvalues.size();i++) // For the size of the eigenvalues.
  {                                     // Loop over the eigenvalues.
    std::cout<<eigenvalues[i]<<" ";     // Print the eigenvalue.
  }                                     // Done printing the eigenvalues.
  std::cout<<std::endl;                 // Print a new line.
}                                       // ------- PrintEigenvalues -----------

template <typename T>
void EigenSpectral<T>::PrintEigenvectors(void)
{                                       // ------- PrintEigenvectors -----------
  std::cout<<"Eigenvectors: "<<std::endl; // Print the eigenvectors.
  for (size_t i=0;i<eigenvectors.size();i++) // For the size of the eigenvectors.
  {                                     // Loop over the eigenvectors.
    for (size_t j=0;j<eigenvectors[0].size();j++) // For the size of the eigenvectors.
    {                                   // Loop over the eigenvectors.
      std::cout<<eigenvectors[i][j]<<" "; // Print the eigenvector.
    }                                   // Done printing the eigenvector.
    std::cout<<std::endl;                // Print a new line.
  }                                     // Done looping over the eigenvectors.
}                                       // ------- PrintEigenvectors -----------

template <typename T>
void EigenSpectral<T>::PrintCovarianceMatrix(void)
{                                       // ------- PrintCovarianceMatrix -----------
  std::cout<<"Covariance Matrix: "<<std::endl; // Print the covariance matrix.
  for (size_t i=0;i<covarianceMatrix.size();i++) // For the size of the covariance matrix.
  {                                     // Loop over the covariance matrix.
    for (size_t j=0;j<covarianceMatrix[0].size();j++) // For the size of the covariance matrix.
    {                                   // Loop over the covariance matrix.
      std::cout<<covarianceMatrix[i][j]<<" "; // Print the covariance matrix.
    }                                   // Done printing the covariance matrix.
    std::cout<<std::endl;                // Print a new line.
  }                                     // Done looping over the covariance matrix.
}                                       // ------- PrintCovarianceMatrix -----------

template <typename T>
void EigenSpectral<T>::PrintPhaseUnwrapped(void)
{                                       // ------- PrintPhaseUnwrapped -----------
  std::cout<<"Unwrapped Phase: "<<std::endl; // Print the unwrapped phase.
  const int n=GetPhaseUnwrapped().size(); // Get the size of the unwrapped phase.
  auto v=GetPhaseUnwrapped();           // Get the unwrapped phase.
  for (size_t i=0;i<n;i++) // For the size of the unwrapped phase.
  {                                     // Loop over the unwrapped phase.
    std::cout<<v[i]<<" ";  // Print the unwrapped phase.
  }                                     // Done printing the unwrapped phase.
  std::cout<<std::endl;                 // Print a new line.
}                                       // ------- PrintPhaseUnwrapped -----------

template <typename T>
void EigenSpectral<T>::PrintEigenDecomposition(void)
{
  std::cout<<"Eigenvalues: ";           // Print the eigenvalues.
  for (size_t i=0;i<eigenvalues.size();i++) // For the size of the eigenvalues.
  {                                     // Loop over the eigenvalues.
    std::cout<<eigenvalues[i]<<" ";     // Print the eigenvalue.
  }                                     // Done printing the eigenvalues.
  std::cout<<std::endl;                 // Print a new line.
  std::cout<<"Eigenvectors: "<<std::endl; // Print the eigenvectors.
  for (size_t i=0;i<eigenvectors.size();i++) // For the size of the eigenvectors.
  {                                     // Loop over the eigenvectors.
    for (size_t j=0;j<eigenvectors[0].size();j++) // For the size of the eigenvectors.
    {                                   // Loop over the eigenvectors.
      std::cout<<eigenvectors[i][j]<<" "; // Print the eigenvector.
    }                                   // Done printing the eigenvector.
    std::cout<<std::endl;                // Print a new line.
  }                                     // Done looping over the eigenvectors.
}


#endif // EIGENSPECTRAL_H
