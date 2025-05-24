/**
 * A class that contains methods for vector operations. 
 * 
 * Author: 
 *   J.EP J. Enrique Peraza
 */
#ifndef VECTORS_H
#define VECTORS_H
#include <vector>
#include <initializer_list>
#include <stdexcept>
#include <type_traits>
#include <complex>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <utility>

template <typename T>
class Vectors{
public:
  // Constructors and destructors:
  Vectors()=default;                    // Empty vector constructor.
  explicit Vectors(size_t n):           // Zero initialized length n 
    vect(n,T{}){}                       // vector constructor.
  Vectors(std::initializer_list<T> init):// Initializer list constructor allows
    vect(init){}                        // to initialize the vector with a list of values.
  Vectors(const Vectors&)=default;      // Copy constructor to initialize the vector with another vector.
  Vectors(Vectors&&)=default            // Move constructor to initialize the vector with another vector.
  virtual ~Vectors()=default;
  
  
  // Element accessors:
  T &operator[](size_t i) noexcept { return vect[i]; }
  const T &operator[](size_t i) const  noexcept { return vect[i]; }
  T& at(size_t i)
  {
    if (i>=vect.size())
      throw std::out_of_range{"Vectors::at(): Index out of range!"};
    return vect[i];                   // Return the element at index i.
  }
  const T& at(size_t i) const
  {
    if (i>=vect.size())
      throw std::out_of_range{"Vectors::at(): Index out of range!"};
    return vect[i];                   // Return the element at index i.
  }
  // Front and back accessors:
  T &front() noexcept { return vect.front(); } // Returns a reference to the first element.
  const T &front() const noexcept { return vect.front(); } // Returns a const reference to the first element.
  T &back() noexcept { return vect.back(); } // Returns a reference to the last element.
  const T &back() const noexcept { return vect.back(); } // Returns a const reference to the last element.
  // Push and pop methods:
  void push_back(const T &value) { vect.push_back(value); } // Adds an element to the end of the vector.
  void push_back(T &&value) { vect.push_back(std::move(value)); } // Adds an rvalue element to the end of the vector.
  void pop_back() { vect.pop_back(); }  // Removes the last element from the vector.
  void clear() noexcept { vect.clear(); }// Clears the vector, removing all elements.
  void resize(size_t n, const T &value = T{}) // Resizes the vector to n elements, initializing new elements with value.
  {                                     // ----------- resize ------------ //
    vect.resize(n, value);              // Resize the vector to n elements.
  }                                     // ----------- resize ------------ //
  void reserve(size_t n) { vect.reserve(n); } // Reserves space for n elements, but does not change the size of the vector.
  // Raw data accessor:
  const std::vector<T>& data() const noexcept { return vect; } // Returns a const reference to the underlying vector data.
  std::vector<T>& data() noexcept { return vect; } // Returns a reference to the underlying vector data.
  // Capacity and size:
  size_t size() const noexcept { return vect.size(); } // Returns the size of the vector.
  bool empty() const noexcept { return vect.empty(); } // Returns true if the vector is empty.

  // Iterators:
  auto begin() noexcept { return vect.begin(); } // Returns an iterator to the beginning of the vector.
  auto end() noexcept { return vect.end(); }     // Returns an iterator to the end of the vector.
  auto begin() const noexcept { return vect.begin(); } // Returns a const iterator to the beginning of the vector.
  auto end() const noexcept { return vect.end(); }     // Returns a const iterator to the end of the vector.

  Vectors &operator=(const Vectors &v)  // Copy assignment operator to copy a vector.
  {                                     // ----------- operator= ------------ //
    if (this!=&v)                       // Are we not assigning to ourselves?
      vect=v.vect;                      // Yes, copy the vector.
    return *this;                       // Return the vector.
  }                                     // ----------- operator= ------------ //
  Vectors &operator=(Vectors &&v) noexcept // Move assignment operator to move a vector.
  {                                     // ----------- operator= ------------ //
    if (this!=&v)                       // Are we not assigning to ourselves?
      vect=std::move(v.vect);           // Yes, move the vector.
    return *this;                       // Return the vector.
  }                                     // ----------- operator= ------------ //
  Vectors &operator+=(const Vectors &v) // Addition assignment operator to add a vector.
  {                                     // ----------- operator+= ------------ //
    if (vect.size()!=v.size())          // Are the vectors the same length?
      throw std::invalid_argument{"Vectors:: -= size mismatch: Vectors must have the same length for addition!"};
    for (size_t i=0;i<vect.size();i++)  // For each element in the vector...
      vect[i]+=v[i];                    // Add the elements of the vectors.
    return *this;                       // Return the sum of the vectors.
  }                                     // ----------- operator+= ------------ //
  Vectors &operator-=(const Vectors &v) // Subtraction assignment operator to subtract a vector.
  {                                     // ----------- operator-= ------------ //
    if (vect.size()!=v.size())          // Are the vectors the same length?
      throw std::invalid_argument{"Vectors:: -= size mismatch: Vectors must have the same length for subtraction!"};
    for (size_t i=0;i<vect.size();i++)  // For each element in the vector...
      vect[i]-=v[i];                    // Subtract the elements of the vectors.
    return *this;                       // Return the difference of the vectors.
  }                                     // ----------- operator-= ------------ //
  Vectors &operator*=(const T &scalar)  // Multiplication assignment operator to scale a vector.
  {                                     // ----------- operator*= ------------ //
    if (vect.empty())                   // Is the vector emtpy?
    {                                   // Yes, return the vector.
      std::cerr<<"[ERROR] Vector is emtpy!\n";
      return vect;                      // Return the vector if it is empty.
    }                                   // Done checking if the vector is empty.
    for (size_t i=0;i<vect.size();i++)  // For each element in the vector...
    {                                   // Loop over the elements of the vector.
      vect[i]*=scalar;                  // Scale each element by the scalar value.
    }                                   // Done scaling the vector.
    return *this;                       // Return the scaled vector.
  }                                     // ----------- operator*= ------------ //
  Vectors &operator/=(const T &scalar)      // Division assignment operator to scale a vector.
  {                                     // ----------- operator/= ------------ //
    if (scalar==T{})
      throw std::invalid_argument{"Vectors:: /= division by zero: Cannot divide by zero!"};
    for (size_t i=0;i<vect.size();i++)  // For each element in the vector...
      vect[i]/=scalar;                  // Divide each element by the scalar value.
    return *this;                       // Return the scaled vector.
  }                                     // ----------- operator/= ------------ //
  // ---------------------------------- //
  // Getter for the vector.
  // ---------------------------------- //
  const std::vector<T> &GetVector() const 
  {                                     // ----------- GetVector ------------ //
    return *this->vect;                 // Return the vector.
  }                                     // ----------- GetVector ------------ //
  // ---------------------------------- //
  // Setter for the vector.
  // ---------------------------------- //
  void SetVector(const std::vector<T> &v)    
  {                                     // ----------- SetVector ------------ //
    *this->vect=v;                      // Set the vector to the input vector.
  }                                     // ----------- SetVector ------------ //

  // Non-member friends:
  friend Vectors operator+(Vectors v1, const Vectors& v2) { return v1+=v2; }
  friend Vectors operator-(Vectors v1, const Vectors& v2) { return v1-=v2; }
  friend Vectors operator*(Vectors v, const T& s) { return v*=s; } // Scale vector by scalar.
  friend Vectors operator*(const T&s, Vectors v) { return v*=s; } // Scale vector by scalar.
  friend Vectors operator/(Vectors v, const T& s) { return v/=s; } // Scale vector by scalar.
  // ---------------------------------- //
  // Compute the dot product of two vectors. The vector dot product tell us
  // The similarity between two vectors. The angle between the two vectors
  // can be computed from the dot product using the formula:
  // cos(theta) = (v1 . v2) / (||v1|| * ||v2||)
  // where v1 and v2 are the vectors, ||v1|| and ||v2|| are the norms of the vectors,
  // and theta is the angle between the vectors.
  // ---------------------------------- //
  T VectorDotProduct(const Vectors& v2) const{
    if (v1.size()!=v2.size())           // Are they the same length?
      throw std::invalid_argument{"Vectors:: dot product size mismatch: Vectors must have the same length for dot product!"};
    T sum = T{};                        // Initialize the sum to zero.
    for (size_t i=0;i<size();i++)       // For each element in the vector...
      sum+=std::conj(vect[i])*v2.vect[i];// Compute the dot product.
    return sum;                         // Return the dot product.
  }                                     // --------- VectorDotProduct ------- //
  
  // ---------------------------------- //
  // Compute the norm of a vector. The vector norm is a measure of the length
  // of the vector. It is defined as the square root of the sum of the squares
  // of the elements of the vector. The norm is used to normalize vectors and
  // to compute the angle between vectors.
  // The norm is also known as the Euclidean norm or L2 norm. It is the equivalent
  // of the magnitute of a vector in physics. Thus it is also the power calculation.
  // ---------------------------------- //

  T VectorNorm() const                  // Returns norm of a vector.
  {                                     // ----------- VectorNorm ------------ //
    return std::sqrt(VectorDotProduct(*this));
  }                                     // ----------- VectorNorm ------------ //   
  
  // ---------------------------------- //
  // Transform a vector into a unit vector by dividing each element by the vector
  // norm. This is useful for normalizing vectors in spectral analysis.
  // The function returns a new vector that is the normalized version of the input vector.
  // ---------------------------------- //
  Vectors NormalizeVector()  const      // Retuns a normalized vector.  
  {                                     // ----------- NormalizeVector --------- //
    auto n=VectorNorm();                // Compute the norm of the vector.
    if (n==T{})                         // Is the norm zero?
      throw std::runtime_error{"Vectors:: NormalizeVector(): Cannot normalize a zero vector!"};
    Vectors v2 = *this;                 // Create a copy of the vector to normalize.
    v2/=n;                              // Scale the vector by the norm to normalize it.
    return v2;                          // Return the normalized vector.
  }                                     // ----------- NormalizeVector --------- //
  // ---------------------------------- //
  // Compute the angle between two vectors using the dot product and norms.
  // The angle is computed using the formula:
  // theta = acos((v1 . v2) / (||v1|| * ||v2||))
  // where v1 and v2 are the vectors, ||v1|| and ||v2|| are the norms of the vectors,
  // and theta is the angle between the vectors.
  // The angle is returned in radians.
  // ---------------------------------- //
  T VectorAngle(const Vectors &v2) const//S
  {                                     // ----------- VectorAngle ------------ //
    if (size()!=v2.size())              // Are the vectors the same length?
      throw std::invalid_argument{"Vectors:: angle size mismatch: Vectors must have the same length for angle computation!"};
    T n1=VectorNorm(),n2=v2.VectorNorm();// Compute the norms of the vectors.
    if (n1==T{} || n2==T{})             // Are the norms zero?
      throw std::runtime_error{"Vectors:: VectorAngle(): Cannot compute angle with zero vector!"};
   auto inner=VectorDotProduct(v2)/(n1*n2);// Compute the dot product and normalize it.
   if constexpr (std::is_floating_point_v<T>)// Is T floating point?
   {                                    // Yes so clamp the value to [-1, 1] 
     inner=std::max<T>(T{-1}, std::min<T>(T{1}, inner)); //  to avoid NaN.
     return std::acos(inner);           // Return the angle in radians.
   }                                    // Done checking if floating point.
   constexpr bool iscomplex{            // Check if T is complex.
     std::is_same_v<T, std::complex<float>> ||
     std::is_same_v<T, std::complex<double>> ||
     std::is_same_v<T, std::complex<long double>>
   };                                   //
   else if (iscomplex)                  // Is it complex?
   {                                    // Yes
    return std::arg(inner);             // Return the phase of the complex number.
   }                                    //
   else                                 // Else we don't know how to hadnle it.
   {                                    // So throw an error if it is not
    static_assert(                      //  any of these types:
      std::is_floating_point_v<T> ||
      std::is_same_v<T,std::complex<float>>   ||
      std::is_same_v<T,std::complex<double>>  ||
      std::is_same_v<T,std::complex<long double>>,
      "Vectors::VectorAngle<T> requires T be float, double, or std::complex thereof"
    );                                  //
   }                                    // Done throwing error if not valid type.   
  }                                     // ----------- VectorAngle ------------ //
  // ---------------------------------- //
  // Compute the cross product of two vectors. The cross product is only defined
  // for 3-dimensional vectors. It is a vector that is orthogonal to both input vectors.
  // The cross product is used in physics to compute the torque and angular momentum.
  // It is used in spectral analysis to compute the orthogonal basis of a set of vectors.
  // The cross product is defined as:
  // v1 x v2 = (v1[1]*v2[2] - v1[2]*v2[1], v1[2]*v2[0] - v1[0]*v2[2], v1[0]*v2[1] - v1[1]*v2[0])
  // where v1 and v2 are the vectors.
  // The function returns a new vector that is the cross product of the input vectors.
  // ---------------------------------- //
  Vectors VectorCrossProduct(const Vectors& v2) const
  {                                     // ----------- VectorCrossProduct ------ //
    if (size()!=3||v2.size()!=3)        // Are the vectors 3D?
      throw std::invalid_argument{"Vectors:: cross product size mismatch: Vectors must be 3D for cross product!"};
    return Vectors{                     // Create a new vector for the cross product.
      vect[1]*v2.vect[2] - vect[2]*v2.vect[1], // i component
      vect[2]*v2.vect[0] - vect[0]*v2.vect[2], // j component
      vect[0]*v2.vect[1] - vect[1]*v2.vect[0]  // k component
    };                                  // Return the cross product vector.
  }                                     // ----------- VectorCrossProduct ------ //

private:
  std::vector<T> vect;                             // The vector of type T.
};
#endif // VECTORS_H
