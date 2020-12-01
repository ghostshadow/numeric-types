/// \file

#ifndef NUMERICTYPES_VEC_HPP
#define NUMERICTYPES_VEC_HPP

#include <type_traits>
#include <stdexcept>
#include <cmath>

namespace NumericTypes {
/* ******** Vector ******** */

// Named Accessors for up to 4 elements
/**
 * \brief Named vector element access
 */
enum class VecElement {
	X=0 /**< \brief X-element */,
	Y=1 /**< \brief Y-element */,
	Z=2 /**< \brief Z-element */,
	W=3 /**< \brief W-element (homogeneous coordinates) */
};

/**
 * \brief General N-dimensional Vector
 * \tparam T Elements Type
 * \tparam N Dimensions
 * \pre The elements type has to support common arithmetic operations
 */
template<class T=double, size_t N=3>
class Vec {
	static_assert(std::is_arithmetic<T>::value,
			"vector elements must support arithmetic operations");
public:
	// Constructors
	/**
	 * \brief Default Constructor (produces 0 vector)
	 */
	Vec()=default;

	/**
	 * \brief Construct vector from array
	 * \tparam T2 Array elements type
	 * \param _v The Array
	 * \pre 1. The array elements type has to be convertible to the vector elements type
	 * \pre 2. The array has to have as many elements as the vector does dimensions
	 */
	template<class T2=T>
	explicit Vec(const T2 (&_v)[N]) noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		for(size_t i=0; i<N; ++i) {
			v[i]=_v[i];
		}
	}

	/**
	 * \brief Construct vector from initializer list
	 * \details The elements of the initializer list are copied to the vector elements.
	 * All none filled elements are set to 0, any excess elements are discarded.
	 * \tparam T2 Initializer list elements type
	 * \param il The initializer list
	 * \pre The initializer list elements type has to be convertible to the vector
	 * elements type
	 */
	template<class T2=T>
	Vec(std::initializer_list<T2> il) noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		for(int i=0; i<(il.size()>=N?N:il.size()); ++i) {
			v[i]=*(il.begin()+i);
		}
	}

	/**
	 * \brief Construct vector from another Vec with a different elements type
	 * \tparam T2 Elements type of the other Vec
	 * \param other The other Vec
	 * \pre 1. The elements type of the other Vec has to be convertible to the vector
	 * elements type
	 * \pre 2. The other Vec has to have as many dimensions as the vector
	 */
	template<class T2=T>
	explicit Vec(const Vec<T2, N> &other) noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		for(size_t i=0; i<N; ++i) {
			v[i]=other[i];
		}
	}

	/**
	 * \brief Vector copy constructor
	 * \param other Other instance to be copied
	 */
	Vec(const Vec<T, N> &other) noexcept {
		for(size_t i=0; i<N; ++i) {
			v[i]=other[i];
		}
	}

	// Index Accessors
	/**
	 * \brief Index access for modifyable lvalue
	 * \param i Index (has to be less than vector dimensions, starting at 0)
	 * \return Modifyable lvalue reference of element at index
	 * \throws std::out_of_range if index is greater or equal the vector dimension
	 */
	T &operator[](const size_t i) {
		if(i>=N) {
			throw std::out_of_range("index must be between 0 and size-1");
		}
		return v[i];
	}

	/**
	 * \brief Index access for copy of element
	 * \param i Index (has to be less tan vector dimensions, starting at 0)
	 * \return Copy of element at index
	 * \throws std::out_of_range if index is greater or equal the vector dimension
	 */
	T operator[](const size_t i) const {
		if(i>=N) {
			throw std::out_of_range("index must be between 0 and size-1");
		}
		return v[i];
	}

	/**
	 * \brief Named element access for modifyable lvalue
	 * \details Translates a named index to a numerical index and calls the
	 * numerical index access.
	 * \param e Named index (see VecElement)
	 * \return Modifyable lvalue reference of element at index
	 * \throws std::out_of_range if index is greater or equal the vector dimension
	 */
	T &operator[](const VecElement e) {
		switch(e) {
		case VecElement::X:
			return (*this)[0];
		case VecElement::Y:
			return (*this)[1];
		case VecElement::Z:
			return (*this)[2];
		case VecElement::W:
			return (*this)[3];
		}
	}

	/**
	 * \brief Named element access for copy of element
	 * \details Translates a named index to a numerical index and calls the
	 * numerical index access.
	 * \param e Named index (see VecElement)
	 * \return Copy of element at index
	 * \throws std::out_of_range if index is greater or equal the vector dimension
	 */
	T operator[](const VecElement e) const {
		switch(e) {
		case VecElement::X:
			return (*this)[0];
		case VecElement::Y:
			return (*this)[1];
		case VecElement::Z:
			return (*this)[2];
		case VecElement::W:
			return (*this)[3];
		}
	}

	// Modifying operators
	/**
	 * \brief Asignment from array
	 * \details The array elements are copied to the corresponding vector elements
	 * \tparam T2 Array elements type
	 * \param other The Array
	 * \return lvalue reference of self
	 * \pre 1. The array elements type has to be convertible to the vector elements type
	 * \pre 2. The array has to have as many elements as the vector does dimensions
	 */
	template<class T2=T>
	Vec<T, N> &operator=(const T2 (&other)[N]) noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		for(size_t i=0; i<N; ++i) {
			v[i]=other[i];
		}
		return *this;
	}

	/**
	 * \brief Asignment from Vec with other element type
	 * \tparam T2 Elements type of other Vec
	 * \param other The other Vec
	 * \return lvalue reference of self
	 * \pre 1. The elements type of the other Vec has to be convertible to the vector
	 * elements type
	 * \pre 2. The other Vec has to have as many dimensions as the vector
	 */
	template<class T2=T>
	Vec<T, N> &operator=(const Vec<T2, N> &other) noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		for(size_t i=0; i<N; ++i) {
			v[i]=other[i];
		}
		return *this;
	}

	/**
	 * \brief Asignment copy
	 * \param other Other instance of Vec
	 * \return lvalue reference of self
	 */
	Vec<T, N> &operator=(const Vec<T, N> &other) noexcept {
		for(size_t i=0; i<N; ++i) {
			v[i]=other[i];
		}
		return *this;
	}

	/**
	 * \brief Short asign sum
	 * \details Asigns each element the sum of itself and the corresponding element of
	 * the rhs Vec.
	 * \tparam T2 Elements type of the rhs Vec
	 * \param rhs The rhs Vec
	 * \return lvalue reference of self
	 * \pre 1. The elements type of the other Vec has to be convertible to the vector
	 * elements type
	 * \pre 2. The other Vec has to have as many dimensions as the vector
	 */
	template<class T2=T>
	Vec<T, N> &operator+=(const Vec<T2, N> &rhs) noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		for(size_t i=0; i<N; ++i) {
			v[i]+=rhs[i];
		}
		return *this;
	}

	/**
	 * \brief Short asign difference
	 * \details Asigns each element the difference of itself and the corresponding element of
	 * the rhs Vec.
	 * \tparam T2 Elements type of the rhs Vec
	 * \param rhs The rhs Vec
	 * \return lvalue reference of self
	 * \pre 1. The elements type of the other Vec has to be convertible to the vector
	 * elements type
	 * \pre 2. The other Vec has to have as many dimensions as the vector
	 */
	template<class T2=T>
	Vec<T, N> &operator-=(const Vec<T2, N> &rhs) noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		for(size_t i=0; i<N; ++i) {
			v[i]-=rhs[i];
		}
		return *this;
	}

	/**
	 * \brief Short asign multiplicative scale
	 * \details Asign each element the product of itself and the scaling value
	 * \tparam T2 Type of the scaling value
	 * \param rhs The scaling value
	 * \return lvalue reference of self
	 * \pre The type of the scaling value has to be a scalar, which is convertible
	 * to the vector elements type
	 */
	template<class T2=T, typename std::enable_if<std::is_scalar<T2>::value,
		bool>::type=true>
	Vec<T, N> &operator*=(const T2 &rhs) noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		for(size_t i=0; i<N; ++i) {
			v[i]*=rhs;
		}
		return *this;
	}

	/**
	 * \brief Short asign dividing scale
	 * \details Asign each element the fraction of itself divided by the scaling value
	 * \tparam T2 Type of the scaling value
	 * \param rhs The scaling value
	 * \return lvalue reference of self
	 * \pre The type of the scaling value has to be a scalar, which is convertible
	 * to the vector elements type
	 */
	template<class T2=T, typename std::enable_if<std::is_scalar<T2>::value,
		bool>::type=true>
	Vec<T, N> &operator/=(const T2 &rhs) {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		for(size_t i=0; i<N; ++i) {
			v[i]/=rhs;
		}
		return *this;
	}

	// Non-Modifying operators
	/**
	 * \brief Equality test
	 * \tparam T2 Elements type of the rhs Vec
	 * \param rhs The rhs Vec
	 * \return true if all elements of the vector are equal to the corresponding elements
	 * of the rhs Vec
	 * \pre 1. The elements type of the other Vec has to be convertible to the vector
	 * elements type
	 * \pre 2. The other Vec has to have as many dimensions as the vector
	 */
	template<class T2=T>
	bool operator==(const Vec<T2, N> &rhs) const noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		for(int i=0; i<N; ++i) {
			if(v[i]!=rhs[i]) {
				return false;
			}
		}
		return true;
	}

	/**
	 * \brief Unary Negaive
	 * \return A new Vec instance with each element being the negative of the original
	 * element
	 */
	Vec<T, N> operator-() const noexcept {
		Vec<T, N> neg{};
		for(size_t i=0; i<N; ++i) {
			neg[i]=-v[i];
		}
		return neg;
	}

	/**
	 * \brief Summation
	 * \tparam T2 Elements type of rhs Vec
	 * \param rhs The rhs Vec
	 * \return A new Vec instance with each element being the sum of the corresponding
	 * elements of the vector and the rhs Vec.
	 * \pre 1. The elements type of the other Vec has to be convertible to the vector
	 * elements type
	 * \pre 2. The other Vec has to have as many dimensions as the vector
	 */
	template<class T2=T>
	Vec<T, N> operator+(const Vec<T2, N> &rhs) const noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		Vec<T, N> sum{};
		for(size_t i=0; i<N; ++i) {
			sum[i]=v[i]+rhs[i];
		}
		return sum;
	}

	/**
	 * \brief Difference
	 * \tparam T2 Elements type of rhs Vec
	 * \param rhs The rhs Vec
	 * \return A new Vec instance with each element being the difference of the
	 * corresponding elements of the vector and the rhs Vec.
	 * \pre 1. The elements type of the other Vec has to be convertible to the vector
	 * elements type
	 * \pre 2. The other Vec has to have as many dimensions as the vector
	 */
	template<class T2=T>
	Vec<T, N> operator-(const Vec<T2, N> &rhs) const noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		Vec<T, N> diff{};
		for(size_t i=0; i<N; ++i) {
			diff[i]=v[i]-rhs[i];
		}
		return diff;
	}

	/**
	 * \brief Multiplicative scaled
	 * \tparam T2 The type of the scaling value
	 * \param rhs the scaling value
	 * \return A copy of the vector multiplicative scaled by the scaling value
	 * \pre The type of the scaling value has to be a scalar, which is convertible
	 * to the vector elements type
	 */
	template<class T2=T, typename std::enable_if<std::is_scalar<T2>::value,
		bool>::type=true>
	Vec<T, N> operator*(const T2 &rhs) const noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		Vec<T, N> scaled{};
		for(size_t i=0; i<N; ++i) {
			scaled[i]=v[i]*rhs;
		}
		return scaled;
	}

	/**
	 * \brief Dot product
	 * \tparam T2 Elements type of the rhs Vec
	 * \param rhs The rhs Vec
	 * \return The sum of the product of each corresponding element in the vector and
	 * the rhs Vec
	 * \pre 1. The elements type of the other Vec has to be convertible to the vector
	 * elements type
	 * \pre 2. The other Vec has to have as many dimensions as the vector
	 */
	template<class T2=T>
	T operator*(const Vec<T2, N> &rhs) const noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		T dot{0.};
		for(size_t i=0; i<N; ++i) {
			dot+=v[i]*rhs[i];
		}
		return dot;
	}

	/**
	 * \brief Dividing scaled
	 * \tparam T2 The type of the scaling value
	 * \param rhs the scaling value
	 * \return A copy of the vector scaled by the inverse of the scaling value
	 * \pre The type of the scaling value has to be a scalar, which is convertible
	 * to the vector elements type
	 */
	template<class T2=T, typename std::enable_if<std::is_scalar<T2>::value,
		bool>::type=true>
	Vec<T, N> operator/(const T2 &rhs) const {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		Vec<T, N> scaled{};
		for(size_t i=0; i<N; ++i) {
			scaled[i]=v[i]/rhs;
		}
		return scaled;
	}

	// Named Operations
	/**
	 * \brief Squared magnitude
	 * \return Squared magnitude (i.e. the sum of each element squared)
	 */
	T magnitude2() const noexcept {
		T mag2{0.};
		for(size_t i=0; i<N; ++i) {
			mag2+=v[i]*v[i];
		}
		return mag2;
	}

	/**
	 * \brief Vector magnitude (AKA length of the Vector)
	 * \return Vector magnitude
	 */
	T magnitude() const {
		return std::sqrt(magnitude2());
	}

	/**
	 * \brief Euclidean norm of the vector (AKA magnitude of the vector)
	 * \return Euclidean norm of the vector
	 */
	T norm() const { return magnitude(); }

	/**
	 * \brief Normalized vector
	 * \return A new Vec instance being the unit vector with the same direction
	 * as the vector
	 */
	Vec<T, N> normalized() const {
		return Vec<T, N>(*this)/norm();
	};

	/**
	 * \brief Normalize the vector
	 * \details Scale the vector to a magnitude of 1
	 * \return lvalue reference of self
	 */
	Vec<T, N> &normalize() {
		(*this)/=norm();
		return *this;
	}

	/**
	 * \brief Dot product
	 * \tparam T2 Elements type of the rhs Vec
	 * \param rhs The rhs Vec
	 * \return The sum of the product of each corresponding element in the vector and
	 * the rhs Vec
	 * \pre 1. The elements type of the other Vec has to be convertible to the vector
	 * elements type
	 * \pre 2. The other Vec has to have as many dimensions as the vector
	 */
	template<class T2=T>
	T dot(const Vec<T2, N> &rhs) const noexcept {
		return ((*this) * rhs);
	}

	/**
	 * \brief Cross product
	 * \tparam T2 Elements type of the rhs Vec
	 * \param rhs The rhs Vec
	 * \return A new Vec instance which is orthogonal to the vector and the rhs Vec,
	 * according to the right hand rule, with the magnitude equal to the area of the
	 * parallelogram spanned by the vector and the rhs Vec.
	 * \pre 1. The elements type of the other Vec has to be convertible to the vector
	 * elements type
	 * \pre 2. The other Vec has to have as many dimensions as the vector
	 */
	template<class T2=T>
	Vec<T, N> cross(const Vec<T2, N> &rhs) const noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		static_assert(N==3, "cross product is only defined for 3-dimensions");
		Vec<T, N> cross{};
		cross[0]=v[1]*rhs[2]-rhs[1]*v[2];
		cross[1]=v[2]*rhs[0]-rhs[2]*v[0];
		cross[2]=v[0]*rhs[1]-rhs[0]*v[1];
		return cross;
	}

	// Type Info
	/** \brief Size of the vector (AKA its dimensions) */
	const size_t size=N;
	/** \brief Elements type */
	using element_type=T;

private:
	/** \brief instance elements storage */
	T v[N]{0.};
};


/**
 * \brief Multiplicative scaled
 * \tparam T Elements Type
 * \tparam N Dimensions
 * \tparam T2 The type of the scaling value
 * \param lhs The scaling value
 * \param rhs The vector
 * \return A copy of the vector multiplicative scaled by the scaling value
 * \pre The type of the scaling value has to be a scalar, which is convertible
 * to the vector elements type
 */
template<class T=double,
	class T2=T, size_t N=3, typename std::enable_if<std::is_scalar<T2>::value,
	bool>::type=true>
inline Vec<T, N> operator*(const T2 &lhs, const Vec<T, N> &rhs) noexcept {
	static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
	return rhs*lhs;
}

/** \brief Shorthand for a 3D-Vector of type double */
using Vec3D=Vec<double, 3>;

}

#endif //NUMERICTYPES_VEC_HPP
