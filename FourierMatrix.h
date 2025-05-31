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

#include "Spectral.h"       // your FFT/IFFT routines (SpectralOps<T>)
#include "Matrices.h"       // Matrices<T> for storing full matrices
#include "MatrixSolver.h"   // CirculantEigen via MatrixSolver<T>
#include "Vectors.h"        // Vectors<T> if you need a vector type

#include <complex>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <algorithm> // for std::copy
#include <numeric>   // for std::iota

namespace dsp {

template <typename T = double>
class FourierMatrix : public Matrices<std::complex<T>>
{
public:
    using Complex      = std::complex<T>;
    using BaseMatrix   = Matrices<Complex>;
    using VectorC      = std::vector<Complex>;

    //--------------------------------------------------------------------------------
    // 1) Factory: Generate the N×N DFT matrix (orthonormal if requested).
    //--------------------------------------------------------------------------------
    static FourierMatrix<T> Generate(std::size_t N, bool orthonormal = true)
    {
        // Scale = 1/√N for an orthonormal DFT, or =1 for the “unnormalized” version.
        const T scale = orthonormal
                         ? T(1)/std::sqrt(static_cast<T>(N))
                         : T(1);

        FourierMatrix<T> F(N, N);
        // We want F[m,n] = scale * exp(−j·2π·m·n / N).  (Some conventions use 2π, some π;
        // here we'll assume the “standard” full‐length DFT: exp(−j·2π·m·n/N).)
        for (std::size_t m = 0; m < N; ++m) {
            for (std::size_t n = 0; n < N; ++n) {
                T angle = T(-2) * M_PI * static_cast<T>(m) * static_cast<T>(n)
                          / static_cast<T>(N);
                F(m,n) = scale * Complex{ std::cos(angle), std::sin(angle) };
            }
        }
        return F;
    }

    // Inherit all constructors (size‐ctor, initializer‐list, copy, move, etc.)
    using BaseMatrix::BaseMatrix;

    //--------------------------------------------------------------------------------
    // 2) Forward/Inverse (FFT‐based) routines
    //    Instead of explicitly doing y = F * x, you just run an FFT on x.
    //    - Forward:  X[k] = ∑_{n=0..N−1} x[n]·exp(−j2πnk/N)
    //    - Inverse: x[n] = (1/N) ∑_{k=0..N−1} X[k]·exp(+j2πnk/N)
    //    (If you built F via Generate(N, orthonormal=true), you won’t need a 1/N factor here.
    //     But in general, SpectralOps<T> can do both unnormalized or normalized versions.)
    //--------------------------------------------------------------------------------

    /// FFT‐style “forward” (i.e. multiply by DFT matrix).  
    /// Expects x.size() == N.  Returns a length‐N vector X.
    std::vector<Complex> Forward(const std::vector<Complex>& x) const
    {
        const std::size_t N = this->Rows();
        if (x.size() != N) {
            throw std::invalid_argument{"FourierMatrix::Forward: size(x) != N"};
        }
        // SpectralOps<T> is assumed to have a method FFT(vector<Complex>&) → vector<Complex>
        SpectralOps<T> spec(static_cast<int>(N));
        return spec.FFT(x);
    }

    /// FFT‐style “inverse” (i.e. multiply by F^H).  
    /// Expects X.size() == N.  Returns a length‐N vector x.
    std::vector<Complex> Inverse(const std::vector<Complex>& X) const
    {
        const std::size_t N = this->Rows();
        if (X.size() != N) {
            throw std::invalid_argument{"FourierMatrix::Inverse: size(X) != N"};
        }
        SpectralOps<T> spec(static_cast<int>(N));
        return spec.IFFT(X);
    }

    //--------------------------------------------------------------------------------
    // 3) Circulant diagonalization: if C is an N×N circulant, its eigenvalues = FFT(first column),
    //    and eigenvectors = the (columns of the) N×N DFT matrix.  We can simply call
    //      MatrixSolver<T>::CirculantEigen(C)
    //    to get {lambda, U}, but we demonstrate how to wrap that here.
    //--------------------------------------------------------------------------------

    /// Diagonalize a circulant matrix C by calling MatrixSolver.  Returns (eigvals, eigvecs).
    std::pair<std::vector<T>, BaseMatrix>
    DiagonalizeCirculant(const Matrices<Complex>& C, T tol = T(1e-9)) const
    {
        std::size_t N = C.Rows();
        if (N != C.Cols())
            throw std::invalid_argument{"FourierMatrix::DiagonalizeCirculant: must be square"};
        if (!C.IsCircular(tol))
            throw std::invalid_argument{"FourierMatrix::DiagonalizeCirculant: not circulant"};

        // MatrixSolver<T> expects real‐valued T.  Here our matrix is Complex<T>,
        // but CirculantEigen in MatrixSolver only calls FFT under the hood and returns
        // complex eigenvalues/vectors.  We just forward the call:
        MatrixSolver<T> solver;
        return solver.CirculantEigen(C);
    }

    //--------------------------------------------------------------------------------
    // 4) Fast Toeplitz × vector multiply via circulant embedding & FFT.
    //
    //    If you have an N×N Toeplitz matrix T whose first column is h[0..N−1],
    //    then for any x[0..N−1], the product y = T * x can be computed by:
    //
    //      * Pick M = next power‐of‐two ≥ 2N−1
    //      * form h_padded[M] = [h[0], h[1], …, h[N−1], 0, …, 0]  
    //      * form x_padded[M] = [x[0], x[1], …, x[N−1], 0, …, 0]
    //      * H = FFT( h_padded )  (length‐M DFT)  
    //      * X = FFT( x_padded )  
    //      * Y[k] = H[k]*X[k]  (pointwise multiply)  
    //      * y_full = IFFT(Y)  
    //      * Then y[0..N−1] = first N entries of y_full.
    //
    //    Complexity:  O(M log M) instead of O(N^2).
    //--------------------------------------------------------------------------------

    /// Return \(\lfloor T * x \rfloor\) in O(M log M), given first column `h` of T.
    static std::vector<Complex>
    FastToeplitzMultiply(const std::vector<T>& h, const std::vector<T>& x)
    {
        std::size_t N = h.size();
        if (x.size() != N) {
            throw std::invalid_argument{"FastToeplitzMultiply: h.size()!=x.size()"};
        }

        // 1) Compute M = next power‐of‐two ≥ 2N−1
        int M = NextPowerOfTwo(static_cast<int>(2*N - 1));
        
        // 2) Zero‐pad h[0..N−1] into length‐M complex array h_pad;
        //    same for x_pad:
        std::vector<Complex> h_pad(M, Complex{0,0}), x_pad(M, Complex{0,0});
        for (std::size_t i = 0; i < N; ++i) {
            h_pad[i] = Complex{h[i], 0};
            x_pad[i] = Complex{x[i], 0};
        }

        // 3) FFT(h_pad), FFT(x_pad)
        SpectralOps<T> spec(M);
        auto H = spec.FFT(h_pad);
        auto X = spec.FFT(x_pad);
        
        // 4) Y[k] = H[k] * X[k]
        std::vector<Complex> Y(M);
        for (int k = 0; k < M; ++k) {
            Y[k] = H[k] * X[k];
        }

        // 5) y_full = IFFT(Y)
        auto y_full = spec.IFFT(Y);

        // 6) The first N entries of y_full are the Toeplitz product
        std::vector<Complex> y(N);
        for (std::size_t n = 0; n < N; ++n) {
            y[n] = y_full[n];
        }
        return y;
    }
  //-------------------------------------------------------------------------------
    // (B) Overload: “Slice out” the first column from a full Matrices<T> and then
    //     call the core FFT routine.  But first *verify* that A.IsToeplitz(tol)==true.
    //
    //     If A is not square or not Toeplitz within `tol`, we throw an exception.
    //-------------------------------------------------------------------------------
    static std::vector<Complex>
    FastToeplitzMultiply(const Matrices<T>& A,
                         const std::vector<T>& x,
                         T tol = T(1e-9))
    {
        // 1) Check squareness
        std::size_t N = A.Rows();
        if (A.Cols() != N) {
            throw std::invalid_argument{
                "FastToeplitzMultiply(A,x): A must be square"
            };
        }

        // 2) Check Toeplitz‐ness
        if (!A.IsToeplitz(tol)) {
            throw std::invalid_argument{
                "FastToeplitzMultiply(A,x): Matrix is not Toeplitz"
            };
        }

        // 3) Extract first column into a vector h[0..N−1]
        std::vector<T> h(N);
        for (std::size_t i = 0; i < N; ++i) {
            h[i] = A(i, 0);
        }

        // 4) Forward to core FFT‐embedding routine
        return FastToeplitzMultiply(h, x);
    }
private:
    // Helper: next power‐of‐two ≥ n
    static int NextPowerOfTwo(int n) {
        int p = 1;
        while (p < n) p <<= 1;
        return p;
    }
};

} // namespace dsp

#endif // FOURIER_MATRIX_H