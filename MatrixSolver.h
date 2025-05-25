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
#include <tuple>
#include <complex>
#include <cmath>
#include <algorithm>
#include <stdexcept>

template <typename T>
class MatrixSolver{
public:
  // Default constructor:
  MatrixSolver(void)=default;
  MatrixSolver(int iters=1000,T defaultTol=T(1e-6))
    :maxIter(iters),tol(defaultTol){} // Constructor with max iterations and tolerance.
  ~MatrixSolver(void)=default;

  // Our dispatcher. SolveEigen(). If B==nullptr, solves A*x=lambda*x. Otherwise,
  // solves A*x=lambda*B*x.
  std::pair<std::vector<T>,Matrices<T>> SolveEigen(
    const Matrices<T>& A,               // The mtrix to solve for eigen-values/vectors
    const Matrices<T>* B=nullptr,       // The matrix to solve for generalized eigen-values/vectors
    int maxIter=1000,                   // The algorithm's maximum number of iterations
    T tol = T(1e-6))                    // The tolerance for convergence
  {                                     // ----------- SolveEigen ------------ //
    // 1. The generalized case:         //
    if (B)                              // Is B nullptr?
      return GeneralizedEigen(A,*B,maxIter,tol); // Yes, solve the generalized eigenvalue problem.
    // 2. Real symmetric / Hermitian?   //
    if (A.IsSymmetric(tol))             // Is A Real symmetric or Hermitian?
      return JacobiEigen(A,maxIter,tol);// Yes, use the Jacobi method for eigenvalue computation.
    // Fallback to the QR iteration algorithm:
    return QREigen(A,maxIter,tol);      // No, use the QR iteration algorithm for eigenvalue computation.
  }                                     // ----------- SolveEigen ------------ //

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
      Matrices<T> L(N,N);L.fill(T{});   // Initialize L as an NxN matrix filled with zeros.
      for (size_t i=0;i<N;i++)          // For each row in A...
      {                                 //
        for (size_t j=0;j<=i;j++)       // For each column in the current row...
        {                               // Compute the elements of L.
          std::complex<T> sum{};        // Initialize the sum to zero.
          for (size_t k=0;k<j;k++)       // For each column in the current row...
            sum+=L(i,k)*std::conj(L(j,k));// Compute the sum of products of the elements in L.
          if (i==j)                     // Are we on the diagonal?
          {                             // Yes, compute the diagonal element.
            T diag=A(i,i)-sum.real();   // Compute the diagonal element of L.
            if (diag<=T{})              // Is the diagonal element non-positive?
              throw std::runtime_error{"CholeskyDecomposition: Matrix is not positive definite!"}; // No, the matrix is not positive definite.
            L(i,i)=std::sqrt(diag);     // Set the diagonal element of L to the square root of the diagonal element of A.
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