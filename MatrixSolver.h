/**
 * A generalized matrix solver used for DSP applications.
 * 
 * Author:
 * J.EP J. Enrique Peraza
 */
#ifndef MATRIX_SOLVER_H
#define MATRIX_SOLVER_H
#include "Matrices.h"
#include "Vectors.h"
#include "Spectral.h"
#include <tuple>
#include <complex>
#include <cmath>
#include <algorithm>
#include <stdexcept>
using std::size_t;

template <typename T>
class MatrixSolver{
public:
  // Default constructor:
  MatrixSolver(void)=default;
  MatrixSolver(int iters=1000,T defaultTol=T(1e-6))
    :maxiters(iters),tolerance(defaultTol){} // Constructor with max iterations and tolerance.
  ~MatrixSolver(void)=default;

  // Our dispatcher. SolveEigen(). If B==nullptr, solves A*x=lambda*x. Otherwise,
  // solves A*x=lambda*B*x.
  std::pair<std::vector<T>,Matrices<T>> SolveEigen(
    const Matrices<T>& A,               // The mtrix to solve for eigen-values/vectors
    const Matrices<T>* B=nullptr)       // The matrix to solve for generalized eigen-values/vectors)                    
  {                                     // ----------- SolveEigen ------------ //
    // 1. The generalized case:         //
    if (B)                              // Is B nullptr?
      return GeneralizedEigen(A,*B,maxiters,tolerance); // Yes, solve the generalized eigenvalue problem.
    // -------------------------------- //
    // 2. If A is circulant, its eigenvalues are simply the DFT of the
    // first row, and eigenvectors are (normalized) Fourier basis vectors.
    // So we can do an FFT to get the eigenvalues and eigenvectors.
    // -------------------------------- //
    if (A.IsCircular(tol))              // Is A circulant?
      return CirculantEigen(A);         // Yes, use the circulant eigenvalue solver.
    // -------------------------------- //
    // 3. If A is Toeplitz, we can use the Levinson-Durbin algorithm to solve
    // the Toeplitz system and get the eigenvalues and eigenvectors.
    // -------------------------------- //
    if (A.IsToeplitz(tol))              // Is A Toeplitz?
      return ToeplitzEigen(A,maxiters,tolerance);// Yes, use the Toeplitz eigenvalue solver.
    // -------------------------------- //
    // 4. If A is positive definite will do a Tridiagonalization and then
    //   and use the QL algorithm to compute the eigenvalues and eigenvectors.
    // -------------------------------- //
    if (A.IsPositiveDefinite(tolerance))      // Is A positive definite?
      return TridiagonalEigen(A,maxIter,tol); // Yes, use the tridiagonal eigenvalue solver.
    // -------------------------------- //
    // 5. If A is postive semi-definite we will do a "rank-revealing" Cholesky + QR
    // on the rank deficient part of A, and then use the QL algorithm to compute
    // the eigenvalues and eigenvectors.
    // ------------------------------- //
    if (A.IsSemiPositiveDefinite(tolerance))  // Is A positive semi-definite?
      return SemiPositiveDefiniteEigen(A,maxIter,tol); // Yes, use the semi-positive definite eigenvalue solver.
    // -------------------------------- //
    // 6. Real symmetric / Hermitian?   //
    // If A is real symmetric or Hermitian, we will use the Jacobi method for eigenvalue computation.
    // -------------------------------- //
    if (A.IsSymmetric(tolerance))             // Is A Real symmetric or Hermitian?
      return JacobiEigen(A,maxiters,tolerance);// Yes, use the Jacobi method for eigenvalue computation.
    // -------------------------------- //
    // Fallback to the QR iteration algorithm:
    // -------------------------------- //
    return QREigen(A,maxiters,tolerance);      // No, use the QR iteration algorithm for eigenvalue computation.
  }                                     // ----------- SolveEigen ------------ //

  // CirculantEigen: Uses the FFT to compute the eigenvalues and eigenvectors
  std::pair<std::vector<T>, Matrices<std::complex<T>>>
  CirculantEigen(const Matrices<T>& A)  // Use FFT to compute eigenvalues/vectors.
  {                                     // ---------- CirculantEigen --------- //
    const size_t N = A.Rows();          // The number of samples (rows or cols) in matrix A.
    if (A.Cols() != N) 
      throw std::invalid_argument{"CirculantEigen: Matrix must be square!"};
    if (!A.IsCircular())
      throw std::invalid_argument{"CirculantEigen: Matrix must be circulant!"};
    // -------------------------------- //
    // Grab the first row, build complex vector.
    //   (real part is the first row, imaginary part is zero)
    // This is the input to the FFT.
    // -------------------------------- //
    std::vector<std::complex<T>> firstRow(N);
    for (size_t j = 0; j < N; j++)
      firstRow[j] = std::complex<T>(A(0, j),T{0});
    // -------------------------------- //
    // Do a length N FFT on the first row to get the Fourier basis vectors,
    // and for the real egeinvalues vector.
    // -------------------------------- //  
    auto [eigvalsC, eigvecsC]=SpectralOps<T>::FFTStrideEig(firstRow);
    std::vector<T> realVals(N);
    for (size_t k = 0; k < N; k++)
      realVals[k] = static_cast<T>(std::real(eigvalsC[k]));
    // -------------------------------- //
    // Pack eigenvector-matrix (complex) into Matrices< complex<T> >
    // The eigenvectors are the Fourier basis vectors, which are complex.
    // Each column of the eigenvector matrix corresponds to a Fourier basis vector.
    // --------------------------------- //
    Matrices<std::complex<T>> Vc(N, N);  // Our complex matrix of eigenvectors.
    for (size_t i=0;i<N;i++)             // For each eigenvector...
      for (size_tj=0;j<N;j++)            // For each element in the eigenvector...
        Vc(i,j)=eigvecsC[i][j];          // Assign the complex eigenvector element to the matrix.
    return { realVals, Vc };             // Return the real eigenvalues and complex eigenvectors.
  }                                      // ---------- CirculantEigen --------- //
  // ToeplitzEigen: Uses the Tridiagonalization algorithm to compute the eigenvalues and eigenvectors
  std::pair<std::vector<T>, Matrices<T>>
  ToeplitzEigen(const Matrices<T>& A,   // The input Toeplitz matrix A.
    int maxiters=1000,                  // Max number of iterations for convergence.
    T tol = T(1e-6))                    // Tolerance for convergence.
  {                                     // ---------- ToeplitzEigen --------- //
    const size_t N=A.Rows();            // The number of samples (rows or cols) in matrix A.
    if (A.Cols()!=N)                    // Is the matrix square?
      throw std::invalid_argument{"ToeplitzEigen: Matrix must be square!"};
    if (!A.IsToeplitz(tol))             // Is the matrix Toeplitz?
      throw std::invalid_argument{"ToeplitzEigen: Matrix is not Toeplitz!"};
    if (N==1)                           // Is the matrix 1x1 (trivial case)?
      return { std::vector<T>{ A(0,0) }, Matrices<T>(1,1) };
    return TridiagonalEigen(A, maxiters, tol);// Tridiagonalize, get eigval/vecs.
  }                                     // ---------- ToeplitzEigen --------- //
  // TridiagonalEigen: Uses the Householder reduction to tridiagonal form and then
  // applies the QL algorithm to compute the eigenvalues and eigenvectors.  
  std::pair<std::vector<T>, Matrices<T>>
  TridiagonalEigen(const Matrices<T>& A,// The input matrix A.
    int maxiters= 1000,                 // Max number of iterations for convergence.
    T tol=T(1e-6))                      // Tolerance for convergence.
  {                                     // ---------- TridiagonalEigen --------- //
    const size_t N=A.Rows();            // The number of samples (rows or cols) in matrix A.
    if (A.Cols()!=N)                    // Is the matrix square?
      throw std::invalid_argument{"TridiagonalEigen: Matrix must be square!"};
    if (!A.IsSymmetric(tol))            // Is the matrix symmetric or Hermitian?
      throw std::invalid_argument{"TridiagonalEigen: Matrix is not symmetric/Hermitian!"};
    // --------------------------------- //
    // Initialize the diagonal and off-diagonal vectors.
    // d[i] = T(i,i), e[i] = T(i, i+1) (off-diagonal)
    // where T is the tridiagonal matrix.
    // -------------------------------- //
    std::vector<T> d(N, T{}),e(N,T{});  // Diagonal and off-diagonal vectors.
    Matrices<T> V(N, N);                // Initialize the eigenvector matrix V.
    V.fill(T{});                        // Fill V with zeros.
    for (size_t i=0;i<N;i++)            // For each element... 
      V(i,i)=T{1};                      // Set the diagonal of V to 1 (identity matrix).
    Matrices<T> C = A;                  // Copy the input matrix A to C for manipulation.
    // -------------------------------- //
    // Begin performing Householder reduction (for k=0..N-2)
    // -------------------------------- //
    for (size_t k=0;k <N-1;k++)         // For each column k in the matrix...
    {                                   // Start Householder reduction.
      size_t len=N-k;                   // Length of the current column vector.
      std::vector<T> x(len);            // Vector to hold the current column vector.
      for (size_t i=k;i<N;i++)          // Fill the vector x with the current column of C.
        x[i-k]=C(i,k);                  // Copy the k-th column of C into x. 
      // ------------------------------ //
      // Compute the Householder vector v and the Householder transformation.
      // We want to zero out the sub-diagonal elements of the k-th column.
      // The Householder vector v is computed as:
      // v = x - alpha * e1, where e1 is the first basis vector.
      // alpha = -sign(x[0]) * ||x||2, where ||x||2 is the 2-norm of x.
      // We will use the Householder transformation H_k = I - 2 * v * v^H
      // to zero out the sub-diagonal elements of the k-th column.
      // ------------------------------ //
      long double sumnorm = 0.0L;       // Sum of squares for the 2-norm.
      for (auto const& xv : x)          // Compute the 2-norm of the vector x.
        sumnorm += std::norm(xv);       // Sum of squares of the elements in x.
      long double normx = std::sqrt(sumnorm);// Compute the 2-norm of x.
      if (normx<tol)                    // Is the 2-norm of x less than tol?
      {                                 // Yes, so skip the column
        // Already zero ? T(i,k)=0 for i>k
        e[k] = T{};                     // Set the off-diagonal element to zero.
        d[k] = C(k, k);                 // Set the diagonal element to the k-th diagonal of C.
        continue;                       // Skip to the next column.
      }                                 // Compute alpha = -sign(x[0]) * ||x||2
      T alpha;                          // Initialize alpha.
      if (x[0]==T{})                    // Is the first element of x zero?
        alpha = static_cast<T>(-normx); // Yes, set alpha to -||x||2.
      else                              // Else compute alpha based on the sign of x[0].
        alpha = (x[0]<T{})?static_cast<T>(normx):static_cast<T>(-normx);
      // ------------------------------ //
      // Get the Householder vector v.
      // v = x - alpha e1
      // ------------------------------ //
      std::vector<T> v=x;               // Copy x to v.
      v[0]-=alpha;                      // Subtract alpha from the first element of v.
      // ------------------------------ //
      // normalize v so ||v||2 = 1
      // ------------------------------ //
      long double vnorm2 = 0.0L;        // Initialize the 2-norm of v.
      for (auto const& vi : v)          // Compute the 2-norm of v.
        vnorm2 += std::norm(vi);        // Sum of squares of the elements in v.
      if (vnorm2 < tol)                 // Is the 2-norm of v less than tol?
      {                                 // Yes, so skip the column.
        e[k] = T{};                     // Set the off-diagonal element to zero.
        d[k] = C(k, k);                 // Set the diagonal element to the k-th diagonal of C.
        continue;                       // Skip to the next column.
      }                                 //
      // ------------------------------ //
      // Normalize v to have unit 2-norm.
      // ------------------------------ //
      long double vnorm=std::sqrt(vnorm2);// Compute the 2-norm of v.
      for (auto& vi : v)                // For each element in v...
        vi = static_cast<T>(vi/static_cast<T>(vnorm));// Normalize the element by dividing by the 2-norm.
      // ------------------------------ //
      // Apply H_k = I - 2 v v^H to C on rows/cols k..N-1
      // This will zero out the sub-diagonal elements of the k-th column.
      // ------------------------------ //
      for (size_t i=k;i<N;i++)          // For each row i in matrix C starting from k...
      {                                 // ...and
        for (size_t j=k;j<N;j++)        // For each col j in matrix C starting from k...
        {                               // Apply Householder transformation.
          std::complex<T> dot=static_cast<std::complex<T>>(v[i-k]) 
            *std::conj(static_cast<std::complex<T>>(v[j-k]));
          C(i,j)-=T(2)*static_cast<T>(dot);// Subtract 2 * v * v^H from C(i,j).
        }                               // Done with this column.
      }                                 // 
      // ------------------------------ //
      // Update the eigenvector matrix V with the Householder transformation.
      // V = V * H_k, where H_k is the Householder transformation.
      // Accumulate V := V * H_k (Householder update).
      // ------------------------------ //
      for (size_t i=0;i<N;i++)          // For each row i in matrix V...
      {                                 // ...and
        std::complex<T> tau=T{};        // Initializing first tau to zero...
        for (size_t j=k;j<N;j++)        // For each col j in mat V starting from K.
        {                               // Compute dot product of V(i,j) and v[j-k].
          tau+=static_cast<std::complex<T>>(V(i,j)) 
            *std::conj(static_cast<std::complex<T>>(v[j-k]));
        }                               // Done computing dot product.
        // Now tau is the dot product of the i-th row of V and the Householder vector v.
        tau*=T(2);                      // Scale tau by 2.
        for (size_t j=k;j<N;j++)        // For each col j in mat V starting from K.
          V(i,j)-=tau*std::conj(v[j-k]);// Subtract tau * v[j-k] from V(i,j).
      }                                 // Done updating the eigenvector matrix V.
      // ------------------------------ //
      // Now C is partially tridiagonal in rows/cols = k+1
      // ------------------------------ //
      e[k] = C(k,k+1);                  // Set off-diagonal element e[k] to C(k,k+1).
      d[k] = C(k,k);                    // Set diagonal element d[k] to C(k,k).
    }                                   // Done with Householder reduction.
    // Last diagonal                    //
    d[N-1] = C(N-1, N-1);               // Set last diagonal elem d[N-1] to C(N-1,N-1).
    e[N-1] = T{};                       // Set last off-diagonal elem e[N-1] to zero.
    // -------------------------------- //
    // QL iteration on tridiagonal (d[], e[]) to rfine eigenvalues and eigenvectors.
    //  in d[], accumulate rotations into U
    // -------------------------------- //
    Matrices<T> U = V;                  // Initialize the eigenvector matrix U with V.
    for (size_t i=0;i<N;i++)            // For each element in the diagonal vector...
      if (std::abs(e[i]) < tol)         // Is the off-diagonal elem less than tol?
        e[i] = T{};                     // Set the off-diagonal elem to zero.
    for (size_t l=0;l<N;l++)            // For each element l in d...
    {                                   // Begin QL iteration.
      int iter = 0;                     // Initialize iteration counter.
      while (true)                      // While not done QL iteration...
      {                                 // Find the first sub-diagonal element e[l].
        // Find m >= l such that e[m] > 0
        size_t m = l;                   // Initialize m to l.
        while (m<N-1&&std::abs(e[m])>tol)// M is within bounds and e[m] is not zero.
          ++m;                          // Increment m until we find a zero or end of the vector.
        if (m==l)                       // Have we found an eigenvalue?
          break;                        // d[l] is already an eigenvalue, so break.
        if (i++ter > maxiters)
          throw std::runtime_error{"TridiagonalEigen: QL iteration did not converge."};
        // ---------------------------- //
        // If we reach here, we need to perform a QL step.
        // We will use Givens rotations to zero out the sub-diagonal elements.
        // We will form a shift mu to improve convergence.
        // Form shift mu
        // ---------------------------- //
        long double dl= static_cast<long double>(d[l]);// Get the diagonal element d[l] as a long double.
        long double dl1=static_cast<long double>(d[l+1]);// Get the next diagonal element d[l+1] as a long double.
        long double el=static_cast<long double>(e[l]);// Get the off-diagonal element e[l] as a long double.
        // ---------------------------- //
        // Compute the shift mu using the formula:
        // mu = (d[l] + d[l+1]) / 2 - sign(d[l+1] - d[l]) * sqrt((d[l] - d[l+1])^2 + 4*e[l]^2) / 2
        // This shift helps to improve convergence of the QL iteration.
        // The sign function is used to determine the direction of the shift.
        // The shift is computed as the average of the two diagonal elements minus a correction term.
        // The correction term is based on the difference between the two diagonal elements
        // and the off-diagonal element.
        // The correction term is computed as the square root of the sum of the squares of the difference
        // between the two diagonal elements and four times the square of the off-diagonal element.
        // ---------------------------- //
        long double diff=dl1-dl;        // Difference between the two diagonal elements.
        long double mu=(dl + dl1)/2.0L  // Compute the average of the two diagonal elements.
          - static_cast<long double>(std::copysign(1.0L, diff))
          * std::sqrt((dl-dl1)*(dl-dl1)+4.0L*el*el)/2.0L;
        long double p=dl-mu;            // Compute the shift p = d[l] - mu.
        long double c=1.0L, s=0.0L;     // Initialize Givens rotation parameters.
        for (size_t i=l;i<N-1;i++)      // For each element i in the diagonal vector....
        {                               // Compute Given rotation to zero out e[i].
          long double r   = std::sqrt(p*p + static_cast<long double>(e[i])*e[i]);
          if (r<tol)                    // Is the norm r less than tol?
          {                             // Yes so considen element i as zero.
            e[i] = T{};                 // Set the off-diagonal element e[i] to zero.
            break;                      // Break out of the loop.
          }                             // No, so continue with the Givens rotation.
          // -------------------------- //
          // Compute the Givens rotation parameters c and s.
          // c = p / r, s = e[i] / r
          // -------------------------- //
          long double c_=p/r;           // Compute the cosine of the Givens rotation.
          long double s_=e[i]/r;        // Compute the sine of the Givens rotation.
          if (i>l)                      // Not first element?
            e[i-1]=static_cast<T>(r);   // Set the previous off-diagonal element e[i-1] to r.
          // Update the diagonal element d[i] using the Givens rotation.
          long double d_next = static_cast<long double>(d[i+1]);
          long double temp   = c_*d_next- s_*static_cast<long double>(e[i+1]);
          e[i + 1]           = static_cast<T>(s_*d_next+c_*static_cast<long double>(e[i+ ]));
          d[i]               = static_cast<T>(mu+r);
          p                  = temp;
          // Update U's columns i, i+1 by Givens [c_, s_]
          for (size_t row=0;row<N;row++)// For each row in the matrix U...
          {                             // Apply the Givens rotation to the i-th and (i+1)-th columns of U.
            // U(row, i) = c_ * U(row, i) - s_ * U(row, i + 1)
            // U(row, i + 1) = s_ * U(row, i) + c_ * U(row, i + 1)
            T Uik=U(row,i);
            T Uik1=U(row,i+1);
            U(row, i)=static_cast<T>(c_*static_cast<long double>(Uik)-s_*static_cast<long double>(Uik1));
            U(row, i+1)=static_cast<T>(s_*static_cast<long double>(Uik)+c_*static_cast<long double>(Uik1));
          }    // Done applying Givens rotation to U(row,i) and U(row,i+1).
        }                               // Done with Givens rotation for this column.
        d[l] = static_cast<T>(mu+p);    // Update the diagonal element d[l] with the shift and the norm.
        if (std::abs(e[l])<tol)         // Is the off-diagonal element e[l] less than tol?
          e[l]=T{};                     // Yes, set it to zero.
      }                                 // Done with QL iteration.
    }
    // -------------------------------- //
    // Now we have the eigenvalues in d[] and the eigenvectors in U.
    // The eigenvalues are the diagonal elements of the tridiagonal matrix.
    // The eigenvectors are the columns of the matrix U.
    // We will sort the eigenvalues and eigenvectors in ascending order.
    // Sort (optional)
    // -------------------------------- //
    for (size_t i=0;i<N-1;i++)          // For each element i in the diagonal vector...
    {                                   //
      size_t minIdx = i;                // Initialize the minimum index to i.
      for (size_t j = i + 1; j < N; j++)// For each element j in the diagonal vector...
        if (d[j]<d[minIdx])             // If d[j] is less than d[minIdx], update minIdx.
          minIdx=j;                     // Update the minimum index to j.
      if (minIdx!=i)                    // If the minimum index is not i, swap d[i] and d[minIdx].
      {                                 // Swap the eigenvalue with the minimum eigenvalue.
        std::swap(d[i], d[minIdx]);     // Swap the eigenvalues.
        for (size_t row = 0; row < N; row++)// For each row in the matrix U...
          std::swap(U(row, i), U(row, minIdx));// Swap the corresponding eigenvectors.
      }                                 // Done swapping eigenvalues and eigenvectors.
    }                                   // Done sorting eigenvalues and eigenvectors.
    return { d, U };                    // Return the sorted eigenvalues and eigenvectors.
  }                                     // ---------- TridiagonalEigen --------- //

  // --------------------------------------------------------------------------------
  // (J)  SemiPositiveDefiniteEigen:  rank-revealing pivoted Cholesky + QL
  //       Steps:
  //        1) Attempt a “pivoted” Cholesky on A.  If A is PSD, no negative pivots will occur.
  //        2) Count how many pivots (diagonal entries of L) are > tol; that number is r = rank(A).
  //        3) Extract the leading r×r principal block of A (call it A1).  A1 is strictly PD.
  //        4) Run TridiagonalEigen(A1) to get its eigenpairs (?1..?_r, U1).  
  //        5) The remaining (N-r) eigenvalues are all zero.  Build an orthonormal basis 
  //           for Null(A) via Gram-Schmidt on the last (N-r) standard basis vectors, or 
  //           simply complete U = [U1, some orthonormal null-basis].  
  //       6) Return:  eigenvalue-vector = [?1..?_r, 0..0] (length N), eigenvector-matrix U (N×N).
  // --------------------------------------------------------------------------------
  std::pair<std::vector<T>, Matrices<T>>
  SemiPositiveDefiniteEigen(
    const Matrices<T>& A,               // The input matrix A.
    int maxIter = 1000,                 // Max number of iterations for convergence.
    T   tol     = T(1e-6))              // Tolerance for convergence.
  {
    const size_t N = A.Rows();
    if (A.Cols() != N)
      throw std::invalid_argument{"SemiPositiveDefiniteEigen: Matrix must be square!"};
    if (!A.IsSymmetric(tol))
      throw std::invalid_argument{"SemiPositiveDefiniteEigen: Matrix must be symmetric/Hermitian!"};
    // 1) Do a pivoted Cholesky directly on A, but keep track of the diagonal pivots.
    Matrices<T> C=A;                    // working copy
    Matrices<T> L(N, N);  L.fill(T{});  // L will hold the Cholesky factor.
    // pivoted: we will store a vector `diag` of length N to hold diagonals of A,
    // then each time pick the largest remaining pivot to proceed.  However, a simpler
    // approach (assuming A is PSD) is just do a normal Cholesky and watch for any zero/negative pivot.
    // We do a no-pivot Cholesky but watch diag < tol ? rank deficiency.
    size_t r{0};                        // count of strictly positive pivots
    for (size_t k=0;k<N;k++)            // For each col k in the matrix...
    {                                   // Start Cholesky decomposition.
      long double sum{0.0L};            // Initialize sum to zero.
      for (size_t m=0;m<k;m++)          // For each row m in the matrix...
        sum += std::norm(L(k, m));      // Compute sum of squares of the elements in the k-th row of L.
      long double diagk=static_cast<long double>(C(k,k))-sum;// Compute the diagonal element diagk = C(k,k) - sum.
      if (diagk<=tol)                   // Is it effectively zero (rank deficient)?
        break;                          // Yes so break out of the loop.
      // ------------------------------ //
      // If we reach here, we have a strictly positive pivot.
      // We will store it in L(k,k) and update the rest of the k-th column.
      // ------------------------------ //
      L(k, k)=static_cast<T>(std::sqrt(diagk));// Store the square root of diagk in L(k,k).
      ++r;                              // we have found the k-th pivot
      for (size_t i=k+1;i<N;i++)        // Starting from k+1 for each row in Matrix...
      {                                 // Update the k-th colum of L.
        // Compute the factor for the i-th row in the k-th column of L.
        // We will use the formula:
        // L(i,k) = (C(i,k) - sum(L(i,m) * conj(L(k,m)))) / L(k,k)
        // where m = 0..k-1, and C(i,k) is the k-th column of C.
        // We will compute the sum of the products of the elements in the i-th row of L
        // and the conjugate of the elements in the k-th row of L.
        // This will give us the factor to update the k-th column of L.
        long double sum2{0.0L};         // Initialize sum2 to zero.
        for (size_t m=0;m<k;m++)        // For each row m in the matrix...
          sum2+=static_cast<long double>(L(i, m))*static_cast<long double>(std::conj(L(k, m)));
        long double factor = (static_cast<long double>(C(i,k))-sum2) 
                             / static_cast<long double>(L(k,k));
        L(i,k)=static_cast<T>(factor);  // Store the factor in L(i,k).
      }                                 // Done updating the k-th column of L.
    }                                   // Done with Cholesky decomposition.
    // -------------------------------- //
    // If r == 0, that means A is (numerically) zero ? all eigenvalues = 0, eigenvectors = I
    // -------------------------------- //
    if (r==0)                           // Is r zero?
      // Yes, so A is numerically zero.
      // All eigenvalues are zero, and the eigenvectors are the identity matrix.
      // So we return a vector of zeros and an identity matrix.
    {                                   //
      std::vector<T> eigvals(N, T{0});  // Eigenvalues are all zero.
      Matrices<T> U(N, N); U.fill(T{}); // Eigenvectors are the identity matrix.
      for (size_t i=0;i<N;i++)          // Fill the identity matrix U.
        U(i,i)=T{1};                    // Set the diagonal elements to 1.
      return { eigvals,U };             // Return the eigenvalues and eigenvectors.
    }                                   // Done checking for r == 0.
    // -------------------------------- //
    // Now we know the top-left r×r block of A is strictly PD.  So we can
    //   extract it and run a strictly-PD eigen solver on it.
    // -------------------------------- //
    Matrices<T> A1(r,r);                // Create a new matrix A1 of size r×r.
    for (size_t i=0;i<r;i++)            // For each row i in the matrix...
      for (size_t j=0;j<r;j++)          // For each column j in the matrix... 
        A1(i, j) = A(i, j);             // Copy the top-left r×r block of A into A1.
    // -------------------------------- //
    // Run a strictly-PD eigen solver on A1
    // -------------------------------- //
    auto [eigvals1, U1]=TridiagonalEigen(A1,maxIter,tol);
    // -------------------------------- //
    // eigvals1.size() == r;  U1 is r×r (orthonormal)
    // So build the final eigenvalue vector of length N:
    // -------------------------------- //
    std::vector<T> eigvals(N,T{0});     // Initialize the eigenvalues vector with zeros.
    for (size_t i=0;i<r;i++)            // For each eigenvalue in eigvals1...
      eigvals[i] = eigvals1[i];         // Copy the eigenvalues from eigvals1 to eigvals.
    Matrices<T> U(N,N);                 // Initialize the eigenvector matrix U of size N×N.
    U.fill(T{});                        // Fill U with zeros.
    for (size_t i=0;i <r;i++)           // For each row i in the matrix... 
      for (size_t j=0;j<r;j++)          // For each column j in the matrix...
        U(i,j) =U1(i,j);                // Copy the eigenvectors from U1 to U.
    for (size_t i=r;i<N;i++)            // For each row i in the matrix U starting from r...
      U(i,i)=T{1};                      // Choose coordinate axes for the null space of A.
    return { eigvals,U};                // Return the eigenvalues and eigenvectors.
  }                                     // -- SemiPositiveDefiniteEigen ----- //
    // LU Decomposition with partial pivoting. Returns (L,U,pivots) such that:
    // P*A=L*U, where P is a permutation matrix.
    std::tuple<Matrices<T>,Matrices<T>,std::vector<size_t>>
    LUDecomposition(const Matrices<T>& A)
    {                                   // ----------- LUDecomposition ------- //
      size_t N=A.Rows();                // Number of rows in the input matrix A.
      Matrices<T> U=A;                  // U is a copy of A, which hold upper triangular matrix.
      Matrices<T> L(N,N);L.fill(T{});   // L is an NxN matrix initialized to zero. (lower triangular matrix).
      std::vector<size_t> P(N);         // P is a vector of size N to hold the pivot indices.
      std::iota(P.begin(),P.end(),0);   // Initialize P with indices from 0 to N-1.
      for (size_t k=0;k<N;k++)          // For each row in the matrix...
      {                                 //
      // ------------------------------ //
      // Find the pivot element in the current column.
      // ------------------------------ //
        size_t pivot=k;                 // Assume the pivot is the current row element.
        T maxv = std::abs(U(k,k));      // Initialize maxv to the absolute value of the pivot element.
        for (size_t i=k+1;i<N;i++)      // For each row below the current row...
        {                               // Test for pivot element.
          T val=std::abs(U(i,k));       // Get the absolute value of the current element.
          if (val>maxv)                 // Is the current element greater than maxv?
          {                             // Yes, update the pivot.
            maxv=val;                   // Update the maximum value.
            pivot=i;                    // Save the index of the pivot element.
          }                             // Done checking for pivot element.
        }                               // Done checking all rows for pivot element.
        if (maxv<tol)                   // Is the pivot element less than the tolerance?
          throw std::runtime_error{"LUDecomposition: Matrix is singular!"}; // Yes, the matrix is singular.
        if (pivot!=k)                   // Is the pivot element the current row?
        {                               // No, we need to swap rows.
          auto &Ud=U.data();            // Get the data of U.
          std::swap(Ud[k],Ud[pivot]);   // Swap the pivot row with the current row in U.
          std::swap(P[k],P[pivot]);     // Swap the pivot indices in P.
          if (k>0)                      // Are we in the first row?
           std::swap_ranges(L.data()[k].begin(),
          L.data()[k].begin()+k,        // No so we will...
          L.data()[pivot].begin());     // Swap the elements in L up to the current row.
        }                               // Done swapping rows if needed.
        // ---------------------------- //
        // Now we will eliminate the entries below the pivot.
        // ---------------------------- //
        for (size_t i=k+1;i<N;i++)      // For each row below the current row...
        {                               // LU elimination step.
          T factor=U(i,k)/U(k,k);       // Compute the factor to eliminate the entry.
          L(i,k)=factor;                // Store the factor in L.
          for (size_t j=k;j<N;j++)      // For each column in the current row...
            U(i,j)-=factor*U(k,j);      // Eliminate the entry by subtracting the factor times the pivot row.
        }                               // Done eliminating entries below the pivot.
      }                                 // Done processing all rows.
      return {L,U,P};                   // Return the matrices L, U and the pivot indices P.
    }                                   // ----------- LUDecomposition ------- //
    // Householder QR Decomposition: A=Q*R, where Q is an orthogonal matrix and 
    // R is upper triangular.
    std::pair<Matrices<T>,Matrices<T>>
    QRDecomposition(const Matrices<T>& A,// The matrix to decompose.
    T tol=T(1e-6))                      // The tolerance for convergence.
    {                                   // -------- QRDecomposition ---------- //
      size_t M=A.Rows(),N=A.Cols();     // Get the number of rows and columns in A.
      Matrices<T> Q(M,M);Q.fill(T{});   // Initialize Q as an MxM matrix filled with zeros.
      for (size_t i=0;i<M;i++)          // For each row in A...
        Q(i,i)=T{1};                    // Set the diagonal elements of Q to 1.  
      Matrices<T> R=A;                  // Copy A to R, which will hold the upper (right) triangular matrix.
      for (size_t k=0;k<std::min(M,N);k++)// For the min of either rows, or columns...
      {                                 //  we will...
        std::vector<T> x(M-k);          // Build the Householder vector x from the k-th column of R.
        for (size_t i=k;i<M;i++)        // For each row below the k-th row...
          x[i-k]=R(i,k);                // Copy the k-th column of R to x.
        // ---------------------------- //
        // alpha = -sign(x0)*||x||
        // ---------------------------- //
        T normx=T{};                    // Initialize norm of x to zero.
        for (auto &xi: x) normx+=std::norm(xi); // Compute the norm of x.
        T alpha=x[0]>=T{}?-normx:normx; // Compute the sign of x[0] and negate the norm.
        // ---------------------------- //
        // v = x-alhpa*e1
        // ---------------------------- //  
        std::vector<T>v=x;              // Initialize v as a copy of x.
        v[0]-=alpha;                    // Subtract alpha from the first element of v.
        // ---------------------------- //
        // Compute v of the Householder reflector H = I - 2*(v*v^T)/(v^T*v)
        // ---------------------------- //
        T vnorm=T{};                    // Compute the norm of v.
        for (auto &vi:v) vnorm+=std::norm(vi); // Compute the norm of v.
        if (vnorm<tol)                  // Is the norm of v less than the tolerance?
          continue;                     // Yes, skip this iteration.
        vnorm=std::sqrt(vnorm);         // Compute the square root of the norm.
        for (auto &vi:v)                // For each element in v...
          vi/=vnorm;                    // Normalize v by dividing each element by the norm.
        // ---------------------------- //
        // Apply the Householder reflector to R.
        // R := H*R
        // ---------------------------- //
        for (size_t j=k;j<N;j++)        // For each column in R...
        {                               //   we will..
          T tau{};                      // Initialize tau to zero.
          for (size_t i=k;i<M;i++)      // For each row in R...
            tau+=std::conj(v[i-k])*R(i,j); // Compute the inner product of v and the j-th column of R.
          tau*=T(2);                    // Scale tau by 2.
          for (size_t i=k;i<M;i++)      // For each row in R...
            R(i,j)-=v[i-k]*tau;         // Apply the Householder reflector to the j-th column of R.
        }                               // Done applying the Householder reflector to R.
        // ---------------------------- //
        // Apply the Householder reflector to Q.
        // Q := Q*H
        // ---------------------------- //
        for (size_t i=0;i<M;i++)        // For each row in Q...
        {                               //   we will...
          T tau{};                      // Initialize tau to zero...
          for (size_t j=k;j<M;j++)      // For each column in Q...
            tau+=Q(i,j)*v[j-k];         // Compute the inner product of the i-th row of Q and v.
          tau*=T(2);                    // Scale tau by 2.
          for (size_t j=k;j<M;j++)      // For each column in Q...
            Q(i,j)-=tau*std::conj(v[j-k]);// Apply the Householder reflector to the i-th row of Q.
        }                               // Done applying the Householder reflector to Q.
      }                                 // Done processing all rows.
      return {Q,R};                     // Return the matrices Q and R.
    }                                   // -------- QRDecomposition ---------- //

    // CholeskyDecomposition: A=LL^H, where L is a lower triangular matrix.
    // for Symmetric Positive Definite matrices (A = L*L^H)
    Matrices<T> CholeskyDecomposition(const Matrices<T>& A)
    {                                   // -------- CholeskyDecomposition ------ //
      size_t N=A.Rows();                // Get the number of rows in A.
      if (A.Cols()!=N)                  // Is A square?
        throw std::invalid_argument{"CholeskyDecomposition: Matrix must be square!"}; // No, so throw an error.
      Matrices<T> C=A;                  // Copy A to C, which will hold the lower triangular matrix.
      Matrices<T> L(N,N);L.fill(T{});   // Initialize L as an NxN matrix filled with zeros.
      for (size_t i=0;i<N;i++)          // For each row in A...
      {                                 //
        for (size_t j=0;j<=i;j++)       // For each column in the current row...
        {                               // Compute the elements of L.
          std::complex<long double> sum=0.0L;// Initialize the sum to zero.
          for (size_t k=0;k<j;k++)       // For each column in the current row...
            sum+=static_cast<long double>L(i,k)*std::conj(static_cast<std::complex<long double>>(L(i,k)));
          if (i==j)                     // Are we on the diagonal?
          {                             // Yes, compute the diagonal element.
            long double diag=static_cast<long double>(C(i,i))-sum.real();
            if (diag<=0.0L)             // Is the diagonal element non-positive?
              throw std::runtime_error{"CholeskyDecomposition: Matrix is not positive definite!"}; // No, the matrix is not positive definite.
            L(i,i)=static_cast<T>(std::sqrt(diag));// Set the diagonal element of L to the square root of the diagonal element of A.
          }                             // Done computing the diagonal element.
          else                          // Else we are in off-diagonal element.                          
            L(i,j)=(A(i,j)-sum)/L(j,j); // Compute the off-diagonal element of L.
        }                               // Done computing the elements of L.
      }                                 // Done processing all rows.
      return L;                         // Return the lower triangular matrix L.
    }                                   // -------- CholeskyDecomposition ------ //
    // JacobiEigen: Computes the eigenvalues and eigenvectors of a real symmetric matrix
    std::pair<std::vector<T>,Matrices<T>>
    JacobiEigen(
      const Matrices<T>& A,             // The matrix to get eigen's from
      int maxIter=1000,                 // The maximum number of iterations
      T tol=T(1e-6))                    // The tolerance for convergence.
    {                                   // ----------- JacobiEigen ------------ //
      size_t N=A.Rows();                // Get the number of rows in A.             
      Matrices<T> D=A;                  // Initialize D as a copy of A, which will hold the diagonal matrix.
      Matrices<T> V(N,N);V.fill(T{});   // Initialize V as an NxN matrix filled with zeros.
      for (size_t i=0;i<N;i++)          // For each row in V...
        V(i,i)=T{1};                    // Set the diagonal elements of V to 1.
      for (int iter=0;iter<maxIter;iter++)// For the max number of iterations
      {                                 //   we will...
        // ---------------------------- //
        // Find the largest off-diagonal element in D.
        // ---------------------------- //
        size_t p=1,q=1;                 // Initialize p and q to 1.
        T maxv=std::abs(D(0,1));        // Initialize maxv to the absolute value of the first off-diagonal element.
        for (size_t i=0;i<N;i++)        // For each row in D...
          for (size_t j=i+1;j<N;j++)    // For each element in the current row.
          {                             // Test for the largest off-diagonal element.
            T val=std::abs(D(i,j));     // Get the absolute value of current element.
            if (val>maxv)               // Current element larger than maxv?
            {                           // Yes, update the indices and maxv.
              maxv=val;                 // Update the maximum value.
              p=i;                      // Update the row index.
              q=j;                      // Update the column index.
            }                           // Done checking for the largest off-diagonal element.
          }                             // Done checking all elems in current row.
          if (maxv<tol)                 // Is the largest off-diagonal element less than the tolerance?
            break;                      // Yes, we have converged. So break.
          T theta=T{0.5}*std::atan2(T{2}*D(p,q),D(q,q)-D(p,p)); // Compute the rotation angle.
          T c=std::cos(theta);          // Compute the cosine of the rotation angle.
          T s=std::sin(theta);          // Compute the sine of the rotation angle.
          // -------------------------- //
          // Now apply the rotation to D.
          // -------------------------- //
          for (size_t i=0;i<N;i++)      // For each row in D...
          {                             //  we will...
            auto Dip=D(i,p);            // Get the current element in D.
            auto Diq=D(i,q);            // Get the current element in D.
            D(i,p)=c*Dip-s*Diq;         // Update the element in D.
            D(p,i)=D(i,p);              // Update the symmetric element in D.
            D(i,q)=s*Dip+c*Diq;         // Update the element in D.
            D(q,i)=D(i,q);              // Update the symmetric element in D.
          }                             // Done updating the elements in D.
          // -------------------------- //
          // Update the diagonals.
          // -------------------------- //
          T dpp=c*c*D(p,p)+s*s*D(q,q)-T(2)*s*c*D(p,q); // Compute the new diagonal element.
          T dqq=s*s*D(p,p)+c*c*D(q,q)+T(2)*s*c*D(p,q); // Compute the new diagonal element.
          D(p,p)=dpp;                   // Update the diagonal element in D.
          D(q,q)=dqq;                   // Update the diagonal element in D.
          D(p,q)=D(q,p)=T{0};           // Set the off-diagonal elements to zero.
          // -------------------------- //
          // Accumulate the rotation in V.
          // -------------------------- //
          for (size_t i=0;i<N;i++)      // For each row in V...
          {                             //  we will...
            auto Vip=V(i,p);            // Get the current element in V.
            auto Viq=V(i,q);            // Get the current element in V.
            V(i,p)=c*Vip-s*Viq;         // Update the element in V.
            V(i,q)=s*Vip+c*Viq;         // Update the element in V.
          }                             // Done updating the elements in V.
      }                                 // Done processing all iterations.
      std::vector<T> vals(N);           // Initialize the eigenvalues vector.
      for (size_t i=0;i<N;i++)          // For each row in D...
        vals[i]=D(i,i);                 // Set the eigenvalue to the diagonal element in D.
      return {vals,V};                  // Return the eigenvalues and eigenvectors.
    }                                   // ----------- JacobiEigen ------------ //
    // QREigen: Computes the eigenvalues and eigenvectors of a matrix using the QR iteration method.
    std::pair<std::vector<T>,Matrices<T>>
    QREigen(
      const Matrices<T>& A,             // To matrix to get eigen's from.
      int maxIter=1000,                 // The maximum number of iterations.
      T tol=T(1e-6))                    // The tolerance for convergence.
      {                                 // ---------- QREigen ---------------- //
        size_t N=A.Rows();              // Get the number of rows in A.
        if (N!=A.Cols())                // Is A a square matrix?
          throw std::invalid_argument{"QREigen: Matrix must be square!"}; // No, throw an exception.
        Matrices<T> X=A;                // Initialize X as a copy of A, which will hold the QR decomposition.
        Matrices<T> Qtot(N,N);Qtot.fill(T{}); // Initialize Qtot as an NxN matrix filled with zeros.
        for (size_t i=0;i<N;i++)        // For each row in Qtot...
          Qtot(i,i)=T{1};               // Set the diagonal elements of Qtot to 1.
        for (int iter=0;iter<maxIter;iter++) // For the max number of iterations...
        {                               //   we will...
          // -------------------------- //
          // Compute the QR decomposition of X.
          // -------------------------- //
          auto [Q,R]=QRDecomposition(X,tol);// Compute QR Decomposition of X.
          X=R*Q;                        // Update X with the product of R and Q.
          Qtot=Qtot*Q;                  // Update Qtot with the product of Qtot and Q.
          // Is the sub-diagonal elements of X less than the tolerance?
          T off{T{}};                   // Initialize off-diagonal sum to zero.
          for (size_t i=1;i<N;i++)      // For each row in X starting from the second row..
            off+=std::abs(X(i,i-1));    // Sum abs value of the sub-diagonal elements.
          if (off<tol)                  // Is the off-diagonal sum less than the tolerance?
            break;                      // Yes, we have converged. So break.
        }                               // Done processing all iterations.
        std::vector<T> vals(N);         // Initialize the eigenvalues vector.
        for (size_t i=0;i<N;i++)        // For each row in X...
           vals[i]=X(i,i);              // Set the eigenvalue to the diagonal element in X.
        return {vals,Qtot};             // Return the eigenvalues and the accumulated orthogonal matrix Qtot.
      }                                 // ---------- QREigen ---------------- //
    // GeneralizedEigen: Computes the generalized eigenvalues and eigenvectors of a matrix.
    std::pair<std::vector<T>,Matrices<T>>
    GeneralizedEigen(
        const Matrices<T>& A,           // The matrix to solve for eigen-values/vectors
        const Matrices<T>& B,           // The matrix to solve for generalized eigen-values/vectors
        int maxIter=1000,               // The algorithm's maximum number of iterations
        T tol = T(1e-6))                // The tolerance for convergence
    {                                   // ----------- GeneralizedEigen ---------- //
      // Attempt Symmetric Positive Definite Cholesky Decomposition if possible.
      if (A.IsSymmetric(tol) && B.IsSymmetric(tol))// Are they both symmetric?
      {                                 // Yes.
        try                             // Try Cholesky on B.
        {                               // Attempt Cholesky Decomposition.
          auto L=CholeskyDecomposition(B);// Decompose B into L*L^H.
          // Form: C = L^{-1}*A*L^{-H}  //
          auto LU=LUDecomposition(L);   // Decompose L into L*U.
          auto Linv=invertFromLU(LU);   // Invert L using LU decomposition.
          Matrices<T> C=Linv*A*Linv.conjugateTranspose(); // Compute C = L^{-1}*A*L^{-H}.
          // Now we can use JacobiEigen to compute the eigenvalues and eigenvectors of C.
          auto [w,Y]=JacobiEigen(C,maxIter,tol); // Compute the eigenvalues and eigenvectors of C.
          // Now we need to transform the eigenvectors back to the original space.
          Matrices<T> X=Linv.conjugateTranspose()*Y; // Transform the eigenvectors back to the original space.
          return {w,X};                 // Return the eigenvalues and the transformed eigenvectors.
        }                               // Done trying Cholesky Decomposition.
        catch (const std::runtime_error& e) // If Cholesky fails, we fallback to QR.
        {                               // Catch the runtime error.
          std::cerr<<"GeneralizedEigen: "<<e.what()<<std::endl; // Print the error message.
        }                               // Done catching the error.
      }                                 // Done checking for symmetric positive definite matrices.
      // Fallback to the QR iteration algorithm for generalized eigenvalue problem.
      auto [L,U,P]=LUDecomposition(B);  // Decompose B into L*U.
      auto Binv=invertFromLU({L,U,P});  // Invert B using LU decomposition.
      return SolveEigen(Binv*A,nullptr,maxIter,tol); // Solve the generalized eigenvalue problem.
    }                                   // ----------- GeneralizedEigen ---------- //
private:
  // Invert a matrix given its LU decomposition.
  Matrices<T> invertFromLU(
    const std::tuple<Matrices<T>,Matrices<T>,std::vector<size_t>>& lures)
  {                                     // ----------- invertFromLU ------------ //
    auto& [L,U,P]=lures;                // Unpack the LU decomposition.
    size_t N=L.Rows();                  // Get the number of rows in L.
    Matrices<T> Pmat(N,N);Pmat.fill(T{}); // Initialize the permutation matrix P as an NxN matrix filled with zeros.
    if (U.Cols()!=N || L.Cols()!=N || L.Rows()!=N || U.Rows()!=N) // Check if L and U are square matrices.
      throw std::invalid_argument{"invertFromLU: L and U must be square matrices!"}; // If not, throw an error.
    if (P.size()!=N)                   // Check if P has the correct size.
      throw std::invalid_argument{"invertFromLU: P must have size N!"}; // If not, throw an error.
    Matrices<T> Pmat(N,N);
    Pmat.fill(T{});
    // Build the P matrix.              // P is a permutation matrix.
    for (size_t i=0;i<N;i++)            // For each row in P...
      Pmat(i,P[i])=T{1};                // Set the permutation matrix elements.
    // -------------------------------- //
    // Solve for each column of I.
    // -------------------------------- //
    Matrices inv(N,N);                  // Initialize the inverse matrix as an NxN matrix filled with zeros.
    for (size_t i=0;i<N;i++)            // For each column of the identity matrix...
    {
      std::vector<T> b(N,T{});          // Initialize the right-hand side vector b as an Nx1 vector filled with zeros.
      b[P[i]]=T{1};                     // Set the j-th element of b to 1.
      // ------------------------------ //
      // Solve L*y=b using forward substitution.
      // ------------------------------ //
      std::vector<T> y(N);              // Initialize the intermediate vector y as an Nx1 vector.
      for (size_t j=0;j<N;j++)          // For each row in L...
      {                                 // we will...
        T s=T{};                        // Initialize the sum s to zero.
        for (size_t k=0;k<j;k++)        // For each column in the current row...
          s+=L(j,k)*y[k];               // Compute the sum of products of the elements in L and y.
        y[j]=(b[j]-s)/L(j,j);           // Compute the i-th element of y.
      }                                 // Done solving for y.
      // ------------------------------ //
      // Solve U*x=y using backward substitution.
      // ------------------------------ //
      std::vector<T> x(N);              // Initialize the solution vector x as an Nx1 vector.
      for (int i=int(N)-1;i>=0;i--)     // For each row in U from the last row to the first...
      {                                 //
        T s=T{};                        // Initialize the sum s to zero.
        for (size_t k=i+1;k<N;k++)      // For each column in the current row...
          s+=U(i,k)*x[k];               // Compute the sum of products of the elements in U and x.
        x[i]=(y[i]-s)/U(i,i);           // Compute the i-th element of x.
      }                                 // Done solving for x.
      // ------------------------------ //
      // Insert the column into the inverse matrix.
      // ------------------------------ //
      for (size_t i=0;i<N;i++)          // For each row in the inverse matrix...
        inv(i,j)=x[i];                  // Set the i-th element of the j-th column of the inverse matrix to x[i].
    }                                   // Done processing all columns.
    return inv;                         // Return the inverse matrix.
  }                                     // ----------- invertFromLU ------------ //

private:
  int maxiters{1000}; // Maximum number of iterations for the algorithms.
  T tolerance{T(1e-6)}; // Tolerance for convergence.

};
#endif 
