/**
 * Description: A object that represents the Fourier Matrix, and its operations.
 * This is to be used in the context of spectral methods, and shall contains the
 * matrix, its properties and operations. More importantly so in the context of
 * Digital Signal Processing (DSP) and Fourier analysis.
 *
 * Author:
 *  J.EP J. Enrique Peraza
 */
#ifndef FOURIER_MATRIX_H
#define FOURIER_MATRIX_H

#include "Spectral.h"       // For Spectral operations.
#include "Matrices.h"       // Matrices<T> for storing full matrices
#include "MatrixSolver.h"   // CirculantEigen via MatrixSolver<T>
#include "Vectors.h"        // For vector operations.

#include <complex>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm> // for std::copy
#include <numeric>   // for std::iota

namespace dsp
{

template <typename T = double>
class FourierMatrix : public Matrices<std::complex<T>>
{
public:
    using Complex      = std::complex<T>;
    using BaseMatrix   = Matrices<Complex>;
    using VectorC      = std::vector<Complex>;

    //--------------------------------------------------------------------------------
    // Factory: Generate the N×N DFT matrix (orthonormal if requested).
    //--------------------------------------------------------------------------------
    static FourierMatrix<T> Generate(   // Generate the N×N Fourier matrix F.
      std::size_t N,                    // The size of the Fourier matrix (N×N).
      bool orthonormal = true)          // If false, scale by 1/sqrt(N) for orthonormal DFT.
    {                                   // ----------- Generate ------------- //
      // ------------------------------ //
      // Generate the N×N Fourier matrix F.
      // Scale = 1/vN for an orthonormal DFT, or =1 for the unnormalized version.
      // ------------------------------ //
      const T scale = orthonormal?T(1)/std::sqrt(static_cast<T>(N)):T(1);
      FourierMatrix<T> F(N,N);          // Create an N×N Fourier matrix F.
      // ------------------------------ //
      // Build F[m,n] = scale * exp(-j·2p·m·n / N).
      // ------------------------------ //
      for (size_t m=0;m<N;m++)          // For each row m in the DFT matrix
        for (size_t n=0;n<N;n++)        // For each column n in the DFT matrix
        {                               // Build the DFT matrix entry.
          T angle=T(-2)*M_PI*static_cast<T>(m)*static_cast<T>(n)/static_cast<T>(N);
          F(m,n)=scale*Complex{ std::cos(angle), std::sin(angle) };
        }                               // Done building the DFT matrix.
      return F;                         // Return the generated Fourier matrix.
    }                                   // ----------- Generate ------------- //

    // Inherit all constructors (size-ctor, initializer-list, copy, move, etc.)
    using BaseMatrix::BaseMatrix;
    //--------------------------------------------------------------------------------
    // 2) Forward/Inverse (FFT-based) routines
    //    Instead of explicitly doing y = F * x, you just run an FFT on x.
    //    - Forward:  X[k] = SUM_{n=0..N-1} x[n]·exp(-j2pnk/N)
    //    - Inverse: x[n] = (1/N) SUM_{k=0..N-1} X[k]·exp(+j2pnk/N)
    //    (If you built F via Generate(N, orthonormal=true), you won't need a 1/N factor here.
    //     But in general, SpectralOps<T> can do both unnormalized or normalized versions.)
    //--------------------------------------------------------------------------------
    /// FFT-style forward (i.e. multiply by DFT matrix).  
    /// Expects x.size() == N.  Returns a length-N vector X.
    std::vector<Complex> Forward(       // Perform the FFT on the input vector x.
      const std::vector<Complex>& x) const// The input signal vector x.
    {                                   // ------------ Forward ------------- //
      const size_t N=this->Rows();      // Get the # of rows in Fourier matrix.
      if (x.size()!=N)                  // Is the input vector x the correct size?
        throw std::invalid_argument{"FourierMatrix::Forward: size(x) != N"};
      SpectralOps<T> spec(static_cast<int>(N));// Create a SpectralOps object for FFT operations.
      return spec.FFT(x);               // Perform the FFT on the input vector x.
    }                                   // ------------ Forward ------------- //
    // -------------------------------- //
    // FFT-style inverse (i.e. multiply by F^H).  
    // Expects X.size() == N.  Returns a length-N vector x.
    // -------------------------------- //
    std::vector<Complex> Inverse(       // Perform IFFT on the input vector X.
      const std::vector<Complex>& X) const// The input spectral vector X.
    {                                   // ------------ Inverse ------------- //
      const size_t N=this->Rows();      // Get the # of rows in Fourier matrix.
      if (X.size()!=N)                  // Is the input vector X the correct size?
          throw std::invalid_argument{"FourierMatrix::Inverse: size(X) != N"};
      SpectralOps<T> spec(static_cast<int>(N));// Our SpectralOps object for IFFT operations.
      return spec.IFFT(X);              // Perform the IFFT on the input vector X.
    }                                   // ------------ Inverse ------------- //

    //--------------------------------------------------------------------------------
    // 3) Circulant diagonalization: if C is an N×N circulant, its eigenvalues = FFT(first column),
    //    and eigenvectors = the (columns of the) N×N DFT matrix.  We can simply call
    //      MatrixSolver<T>::CirculantEigen(C)
    //    to get {lambda, U}, but we demonstrate how to wrap that here.
    //--------------------------------------------------------------------------------
    /// Diagonalize a circulant matrix C by calling MatrixSolver.  Returns (eigvals, eigvecs).
    std::pair<std::vector<T>, BaseMatrix>
    DiagonalizeCirculant(               // Diagonalize a circulant matrix C.
      const Matrices<Complex>& C,       // The circulant matrix C to diagonalize.
      T tol=T(1e-9)) const              // Tolerance for checking circulant property.
    {                                   // ------- DiagonalizeCirculant ----- //
        size_t N = C.Rows();            // Get the number of rows in C.
        if (N != C.Cols())              // Is C square?
            throw std::invalid_argument{"FourierMatrix::DiagonalizeCirculant: must be square"};
        if (!C.IsCircular(tol))         // Is C circulant within the tolerance?
            throw std::invalid_argument{"FourierMatrix::DiagonalizeCirculant: not circulant"};
        MatrixSolver<T> solver;         // Our matrix solver object.
        return solver.CirculantEigen(C);// Call the circulant eigenvalue solver.
    }                                   // ------- DiagonalizeCirculant ----- //

    //--------------------------------------------------------------------------------
    // 4) Fast Toeplitz × vector multiply via circulant embedding & FFT.
    //
    //    If you have an N×N Toeplitz matrix T whose first column is h[0..N-1],
    //    then for any x[0..N-1], the product y = T * x can be computed by:
    //
    //      * Pick M = next power-of-two = 2N-1
    //      * form h_padded[M] = [h[0], h[1], ..., h[N-1], 0, ..., 0]  
    //      * form x_padded[M] = [x[0], x[1], ..., x[N-1], 0, ..., 0]
    //      * H = FFT( h_padded )  (length-M DFT)  
    //      * X = FFT( x_padded )  
    //      * Y[k] = H[k]*X[k]  (pointwise multiply)  
    //      * y_full = IFFT(Y)  
    //      * Then y[0..N-1] = first N entries of y_full.
    //
    //    Complexity:  O(M log M) instead of O(N^2).
    //--------------------------------------------------------------------------------
    /// Return \(\lfloor T * x \rfloor\) in O(M log M), given first column `h` of T.
    static std::vector<Complex>
    FastToeplitzMultiply(               // Operate on x using h.
      const std::vector<T>& h,          // The first col of Toeplitz (LTI) matrix h (a filter, maybe)
      const std::vector<T>& x)          // The input signal vector x.
    {                                   // ------- FastToeplitzMultiply ----- //
      size_t N=h.size();                // Get the number of filter coeffs?
      if (x.size()!=N)                  // Is the input vector x the same size as h?
        throw std::invalid_argument{"FastToeplitzMultiply: h.size()!=x.size()"};
        // ----------------------------- //
        // Compute M = next power-of-two = 2N-1
        // ----------------------------- //
        int M=NextPowerOfTwo(static_cast<int>(2*N-1));
        // ---------------------------- //
        // Zero-pad h[0..N-1] into length-M complex array h_pad;
        //    same for x_pad:
        // ---------------------------- //
        std::vector<Complex> h_pad(M,Complex{0,0}),x_pad(M, Complex{0,0});
        for (size_t i=0;i<N;i++)        // For each element in h and x...
        {                               // Zero-pad them into length-M arrays.
            h_pad[i] = Complex{h[i], 0};// Copy h[i] into h_pad[i] as a complex number.
            x_pad[i] = Complex{x[i], 0};// Copy x[i] into x_pad[i] as a complex number.
        }                               // Done zero-padding h and x.
        // ---------------------------- //
        // Now perform the FFT on the padded vectors:
        // ---------------------------- //
        SpectralOps<T> spec(M);         // Our Spectral Objects for FFT operations.
        auto H = spec.FFT(h_pad);       // Perform FFT on the padded filter h.
        auto X = spec.FFT(x_pad);       // Perform FFT on the padded input x.
        // ---------------------------- //
        // Obtain the output spectrum Y by pointwise multiplying H and X:
        // ---------------------------- //
        std::vector<Complex> Y(M);      // Prepare the output spectrum Y.
        for (int k=0;k<M;k++)           // For each freq bin k in the spectrum...
          Y[k]=H[k]*X[k];               // Pointwise multiply the FFTs H and X.
        // ---------------------------- //
        // Now perform the IFFT on the output spectrum Y:
        // This will give us the full convolution result in y_full.
        // The first N entries of y_full are the Toeplitz product.
        // ---------------------------- //
        auto y_full = spec.IFFT(Y);     // Perform the IFFT on the output spectrum Y.
        std::vector<Complex> y(N);      // Prepare the output vector y of size N.
        for (size_t n=0;n<N;n++)        // For each entry n in the output vector y...
          y[n] = y_full[n];             // Copy the first N entries from y_full to y.
        return y;                       // Return the first N entries of the full convolution result.
    }                                   // ------- FastToeplitzMultiply ----- //
    //-------------------------------------------------------------------------------
    // (B) Overload: Slice out the first column from a full Matrices<T> and then
    //     call the core FFT routine.  But first *verify* that A.IsToeplitz(tol)==true.
    //
    //     If A is not square or not Toeplitz within `tol`, throw an exception.
    //-------------------------------------------------------------------------------
    static std::vector<Complex>         //
    FastToeplitzMultiply(               // Operate on x using h.
      const std::vector<T>& h,          // The first col of Toeplitz (LTI) matrix h (a filter, maybe)
      const std::vector<T>& x)          // The input signal vector x.
    {                                   // ------- FastToeplitzMultiply ----- //
      size_t N=A.Rows();                // Get the number of rows in the matrix A.
      if (A.Cols()!=N)                  // Is A square?
        throw std::invalid_argument{"FastToeplitzMultiply(A,x): A must be square" };
      if (!A.IsToeplitz(tol))           // Is A toeplitz (LTI) within the tolerance?
        throw std::invalid_argument{"FastToeplitzMultiply(A,x): Matrix is not Toeplitz"};
      std::vector<T> h(N);              // Prepare a vector h to hold the first column of A.
      for (size_t i=0;i<N;i++)          // For each row i in the matrix A...
        h[i]=A(i,0);                    // Copy the first col of A into h[i].
      return FastToeplitzMultiply(h, x);// Call the core FastToeplitzMultiply with h and x.
    }                                   // ------- FastToeplitzMultiply ----- //
private:
    // Helper: next power-of-two = n
    static int NextPowerOfTwo(int n) 
    {                                   // -------- NextPowerOfTwo --------- //
        int p=1;                        // 1* 2^k = n, where k is the smallest integer such that 2^k >= n.
        while (p<n) p<<=1;              // Keep doubling p until it is >= n.
        return p;                       // Return the next power of two that is >= n.
    }                                   // -------- NextPowerOfTwo --------- //
};

} // namespace dsp

#endif // FOURIER_MATRIX_H
