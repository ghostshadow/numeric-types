/// \file

#ifndef NUMERICTYPES_QUAT_HPP
#define NUMERICTYPES_QUAT_HPP

#include <type_traits>
#include <stdexcept>
#include <cmath>
#include "Vec.hpp"

namespace NumericTypes {
/* ********** Quaternion ********* */

/**
 * \brief Named quaternion element access
 */
enum class QuatElement {
	X=0 /**< \brief First vector component (AKA i) */,
	Y=1 /**< \brief Second vector component (AKA j) */,
	Z=2 /**< \brief Third vector component (AKA k) */,
	S=3 /**< \brief Scalar component */
};

/**
 * \brief General Quaternion
 * \tparam T Elements type
 * \pre The elements type as to support common arithmetic operations
 */
template<class T=double>
class Quat {
	static_assert(std::is_arithmetic<T>::value,
			"quaternion elements must support arithmetic operations");
public:
	// Constructors
	/** \brief Default constructor (all elements 0) */
	Quat()=default;

	/**
	 * \brief Construct Quat from array
	 * \details Array elements 0-2 are used as the vector components and array
	 * element 3 is used as the scalar component.
	 * \tparam T2 Array elements type
	 * \param q The Array
	 * \pre 1. The array elements type has to be convertible to the quaternion elements
	 * type
	 * \pre 2. The array has to have 4 elements
	 */
	template<class T2=T>
	explicit Quat(const T2 (&q)[4]) noexcept :S(q[3]), V{q[0], q[1], q[2]} {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
	}

	/**
	 * \brief Construct Quat from initializer list
	 * \details The first 3 elements are used for the vector components and the
	 * 4th element is used for the scalar component. All none filled elements
	 * are set to 0, any excess elements are discarded.
	 * \tparam T2 Initializer list elements type
	 * \param il The initializer list
	 * \pre 1. The initializer list elements type has to be convertible to the
	 * quaternion elements type
	 */
	template<class T2=T>
	Quat(std::initializer_list<T2> il) noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		if(il.size()>=1) {
			V[0]=*(il.begin());
		}
		if(il.size()>=2) {
			V[1]=*(il.begin()+1);
		}
		if(il.size()>=3) {
			V[2]=*(il.begin()+2);
		}
		if(il.size()>=4) {
			S=*(il.begin()+3);
		}
	}

	/**
	 * \brief Construct Quat from scalar and vector
	 * \tparam T2 Vector elements type
	 * \tparam T3 Scalar type
	 * \param S The scalar
	 * \param V The vector
	 * \pre 1. The scalar type and the vector elements type have to be convertible
	 * to the quaternion elements type
	 * \pre 2. The vector has to have 3 elements
	 */
	template<class T2=T, class T3=T, typename std::enable_if<std::is_scalar<T3>::value,
		bool>::type=true>
	Quat(const T3 &S, const Vec<T2, 3> &V) noexcept :S(S), V(V) {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		static_assert(std::is_convertible<T3, T>::value, "can not convert element types");
	}

	/**
	 * \brief Construct Quat from scalar and vector
	 * \tparam T2 Vector elements type
	 * \tparam T3 Scalar type
	 * \param S The scalar
	 * \param V The vector
	 * \pre 1. The scalar type and the vector elements type have to be convertible
	 * to the quaternion elements type
	 * \pre 2. The vector has to have 3 elements
	 */
	template<class T2=T, class T3=T, typename std::enable_if<std::is_scalar<T3>::value,
		bool>::type=true>
	Quat(const Vec<T2, 3> &V, const T3 &S) noexcept :Quat(S, V) {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		static_assert(std::is_convertible<T3, T>::value, "can not convert element types");
	}

	/**
	 * \brief Construst pure scalar Quat
	 * \tparam T2 Scalar type
	 * \param S The scalar
	 * \pre The scalar type has to be convertible to the Quat elements type
	 */
	template<class T2=T, typename std::enable_if<std::is_scalar<T2>::value,
		bool>::type=true>
	explicit Quat(const T2 &S) noexcept :S(S) {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
	}

	/**
	 * \brief Constuct pure imaginary Quat
	 * \tparam T2 Vec elements type
	 * \param V The Vec
	 * \pre The Vec elements type has to be convertible to the Quat elements type
	 */
	template<class T2=T>
	explicit Quat(const Vec <T2> &V) noexcept :V(V) {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
	}

	/**
	 * \brief Construct Quat from other Quat with different elements type
	 * \tparam T2 Elements type of the other Quat
	 * \param other The other Quat
	 * \pre The elements type of the other Quat has to be convertible to the
	 * quaternion elements type
	 */
	template<class T2=T>
	explicit Quat(const Quat<T2> &other) noexcept :S(other.scalar()), V(other.vector()) {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
	}

	/**
	 * \brief Construct Quat copy
	 * \param other The other Quat instance
	 */
	Quat(const Quat<T> &other) noexcept :S(other.scalar()), V(other.vector()) {}

	// Prepared Instances
	/**
	 * \brief Create a identiy quaternion
	 * \details Creates a quaternion with the scalar component being 1 and
	 * the vector component being a 0-vector
	 * \return A identity Quat
	 */
	static Quat<T> Identity() {
		return Quat<T>(T{1.}, Vec<T, 3>{0., 0., 0.});
	}

	// Accessors
	/**
	 * \brief Index access for modifyable lvalue
	 * \param i Index (0-3; 0-2 are the vector components, 3 is the scalar component)
	 * \return Modifyable lvalue reference of element at index
	 * \throws std::out_of_range for indices greater or equal to 4
	 */
	T &operator[](const size_t i) {
		if(i>=4) {
			throw std::out_of_range("index out of bound");
		}
		if(i==3) {
			return S;
		}
		return V[i];
	}

	/**
	 * \brief Index access for copy of element
	 * \param i Index (0-3; 0-2 are the vector components, 3 is the scalar component)
	 * \return Copy of element at index
	 * \throws std::out_of_range for indices greater or equal to 4
	 */
	T operator[](const size_t i) const {
		if(i>=4) {
			throw std::out_of_range("index out of bound");
		}
		if(i==3) {
			return S;
		}
		return V[i];
	}

	/**
	 * \brief Named index access for modifyable lvalue
	 * \param e Element name (see QuatElement)
	 * \return Modifyable lvalue reference of element at index
	 */
	T &operator[](const QuatElement e) noexcept {
		switch(e) {
		case QuatElement::X:
			return V[0];
		case QuatElement::Y:
			return V[1];
		case QuatElement::Z:
			return V[2];
		case QuatElement::S:
			return S;
		}
	}

	/**
	 * \brief Named index access for copy of element
	 * \param e Element name (see QuatElement)
	 * \return Copy of element at index
	 */
	T operator[](const QuatElement e) const noexcept {
		switch(e) {
		case QuatElement::X:
			return V[0];
		case QuatElement::Y:
			return V[1];
		case QuatElement::Z:
			return V[2];
		case QuatElement::S:
			return S;
		}
	}

	/**
	 * \brief Get the scalar component as modifyable lvalue
	 * \return Modifyable lvalue reference of the scalar component
	 */
	T &scalar() noexcept {
		return S;
	}

	/**
	 * \brief Get a copy of the scalar component
	 * \return Copy of the scalar component
	 */
	T scalar() const noexcept {
		return S;
	}

	/**
	 * \brief Get the vector components as modifyable lvalue Vec
	 * \return Modifyable lvalue reference of the vector component Vec
	 */
	Vec<T, 3> &vector() noexcept {
		return V;
	}

	/**
	 * \brief Get a copy of the vector components as Vec
	 * \return A Vec containing a copy of the vector components
	 */
	Vec<T, 3> vector() const noexcept {
		return V;
	}

	// Modifying Operators
	/**
	 * \brief Asignment from Quat with other elements type
	 * \tparam T2 Elements type of the other Quat
	 * \param other The other Quat
	 * \return lvalue reference of self
	 * \pre The elements type of the other Quat has to be convertible to the
	 * quaternion elements type
	 */
	template<class T2=T>
	Quat<T> &operator=(const Quat<T2> &other) noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		S=other.scalar();
		V=other.vector();
		return *this;
	}

	/**
	 * \brief Asignment from other Quat
	 * \param other The other Quat
	 * \return lvalue reference of self
	 */
	Quat<T> &operator=(const Quat<T> &other) noexcept {
		S=other.scalar();
		V=other.vector();
		return *this;
	}

	/**
	 * \brief Short asign sum
	 * \details Asign each element the sum of itself and the corresponding element
	 * of the rhs Quat
	 * \tparam T2 Elements type of the rhs Quat
	 * \param rhs The rhs Quat
	 * \return lvalue reference of self
	 * \pre The elements type of the other Quat has to be convertible to the
	 * quaternion elements type
	 */
	template<class T2=T>
	Quat<T> &operator+=(const Quat<T2> &rhs) noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		S+=rhs.scalar();
		V+=rhs.vector();
		return *this;
	}

	/**
	 * \brief Short asign difference
	 * \details Asign each element the difference of itself and the corresponding
	 * element of the rhs Quat
	 * \tparam T2 Elements type of the rhs Quat
	 * \param rhs The rhs Quat
	 * \return lvalue reference of self
	 * \pre The elements type of the other Quat has to be convertible to the
	 * quaternion elements type
	 */
	template<class T2=T>
	Quat<T> &operator-=(const Quat<T2> &rhs) noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		S-=rhs.scalar();
		V-=rhs.vector();
		return *this;
	}

	/**
	 * \brief Short asign multiplicative scale
	 * \details Asign each element the product of itself and the scaling value
	 * \tparam T2 Scaling value type
	 * \param rhs The scaling value
	 * \return lvalue reference of self
	 * \pre The scaling value type has to be a scalar, which is convertible to the quaternion
	 * elements type
	 */
	template<class T2=T, typename std::enable_if<std::is_scalar<T2>::value,
		bool>::type=true>
	Quat<T> &operator*=(const T2 &rhs) noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		S*=rhs;
		V*=rhs;
		return *this;
	}

	/**
	 * \brief Short asign division scale
	 * \details Asign each element the product of itself and inverse of the scaling value
	 * \tparam T2 Scaling value type
	 * \param rhs The scaling value
	 * \return lvalue reference of self
	 * \pre The scaling value type has to be a scalar, which is convertible to the quaternion
	 * elements type
	 */
	template<class T2=T, typename std::enable_if<std::is_scalar<T2>::value,
		bool>::type=true>
	Quat<T> &operator/=(const T2 &rhs) {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		S/=rhs;
		V/=rhs;
		return *this;
	}

	// Non-Modifying Operators
	/**
	 * \brief Equality test
	 * \tparam T2 Elements type of the rhs Quat
	 * \param rhs The rhs Quat
	 * \return true if all elements of the quaternion are equal the corresponding
	 * elements of rhs Quat
	 * \pre The elements type of the other Quat has to be convertible to the
	 * quaternion elements type
	 */
	template<class T2=T>
	bool operator==(const Quat<T2> &rhs) const noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		return S==rhs.scalar() && V==rhs.vector();
	}

	/**
	 * \brief Unary Negative
	 * \details Negates only the scalar part of the quaternion leaving the
	 * vector part uneffected. (This is basically the reversion of the rotation
	 * direction in the geometric interpretation of the quaternion)
	 * \return A new Quat with the scalar part being the negative of the
	 * scalar part of ther quaternion and the vector part being copied
	 */
	Quat<T> operator-() const noexcept {
		return Quat<T>(-S, V);
	}

	/**
	 * \brief Summation
	 * \tparam T2 Elements type of the rhs Quat
	 * \param rhs The rhs Quat
	 * \return A new Quat with each element being the sum of the corresponding
	 * elements from the quaternion and the rhs Quat
	 * \pre The elements type of the other Quat has to be convertible to the
	 * quaternion elements type
	 */
	template<class T2=T>
	Quat<T> operator+(const Quat<T2> &rhs) const noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		Quat<T> sum{};
		sum.scalar()=S+rhs.scalar();
		sum.vector()=V+rhs.vector();
		return sum;
	}

	/**
	 * \brief Difference
	 * \tparam T2 Elements type of the rhs Quat
	 * \param rhs The rhs Quat
	 * \return A new Quat with each element being the difference of the corresponding
	 * elements from the quaternion and the rhs Quat
	 * \pre The elements type of the other Quat has to be convertible to the
	 * quaternion elements type
	 */
	template<class T2=T>
	Quat<T> operator-(const Quat<T2> &rhs) const noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		Quat<T> sum{};
		sum.scalar()=S-rhs.scalar();
		sum.vector()=V-rhs.vector();
		return sum;
	}

	/**
	 * \brief Multiplicative scaled
	 * \tparam T2 The scaling value type
	 * \param rhs The scaling value
	 * \return A new Quat with each element being the product of the quaternion
	 * and the scaling value
	 * \pre The scaling value type has to be a scalar, which is convertible to the quaternion
	 * elements type
	 */
	template<class T2=T, typename std::enable_if<std::is_scalar<T2>::value,
		bool>::type=true>
	Quat<T> operator*(const T2 &rhs) const noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		Quat<T> scaled{};
		scaled.scalar()=S*rhs;
		scaled.vector()=V*rhs;
		return scaled;
	}

	/**
	 * \brief Product
	 * \details Quaternion product formula:
	 * \f{eqnarray*}{
	 * q_{n}\left(S,\vec{V}\right) &=& q_{l}\left(S,\vec{V}\right)\cdot{}
	 *   q_{r}\left(S,\vec{V}\right)\\
	 * S_{n} &=& S_{l}\cdot{} S_{r}-\vec{V}_{l}\cdot{}\vec{V}_{r}\\
	 * \vec{V}_{n} &=& \vec{V}_{r}\cdot{} S_{l}+\vec{V}_{l}\cdot{} S_{r}+
	 *   \vec{V}_{l}\times{}\vec{V}_{r}
	 * \f}
	 * \tparam T2 Elements type of the rhs Quat
	 * \param rhs The rhs Quat
	 * \return A new Quat with the product of the quaternion and the rhs Quat
	 * according to the quaternion product formula
	 * \pre The elements type of the other Quat has to be convertible to the
	 * quaternion elements type
	 */
	template<class T2=T>
	Quat<T> operator*(const Quat<T2> &rhs) const noexcept {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		return Quat<T>{S*rhs.scalar()-V*rhs.vector(),
			rhs.vector()*S+V*rhs.scalar()+V.cross(rhs.vector())};
	}

	/**
	 * \brief Dividing scaled
	 * \tparam T2 The scaling value type
	 * \param rhs The scaling value
	 * \return A new Quat with each element being the product of the quaternion
	 * and inverse of the scaling value
	 * \pre The scaling value type has to be a scalar, which is convertible to the quaternion
	 * elements type
	 */
	template<class T2=T, typename std::enable_if<std::is_scalar<T2>::value,
		bool>::type=true>
	Quat<T> operator/(const T2 &rhs) const {
		static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
		Quat<T> scaled{};
		scaled.scalar()=S/rhs;
		scaled.vector()=V/rhs;
		return scaled;
	}

	// Named Operations
	/**
	 * \brief Conjungate of the quaternion
	 * \return A new Quat with the scalar part being a copy and the vector part
	 * being the negaive of the vector part of the quaternion
	 */
	Quat<T> conjungate() const noexcept {
		return Quat<T>(S, -V);
	}

	/**
	 * \brief Euclidean norm of the quaternion (AKA magnitude)
	 * \return Euclidean norm of the quaternion (i.e. the square root of the
	 * sum of the squares of all its elements)
	 */
	T norm() const {
		return std::sqrt(S*S+V.magnitude2());
	}

	/**
	 * \brief Squared euclidean norm of the quaternion (AKA magnitude squared)
	 * \return Squared euclidean norm of the quaternion (i.e. the sum of the
	 * squares of all its elements)
	 */
	T norm2() const {
		return S*S+V.magnitude2();
	}


	/**
	 * \brief The versor of the quaternion
	 * \return A new Quat being the quaternion scaled to a magnitude of 1
	 */
	Quat<T> versor() const {
		return Quat<T>(*this)/norm();
	}

	/**
	 * \brief Normalized quaternion (AKA the versor)
	 * \return The versor of the quaternion (see versor())
	 */
	Quat<T> normalized() const { return versor(); }

	/**
	 * \brief Normalize the quaternion
	 * \details Scales the quaternion to a magnitude of 1
	 * \return lvalue reference of self
	 */
	Quat<T> &normalize() {
		(*this)/=norm();
		return *this;
	}

	/**
	 * \brief Reciprocal of the quaternion
	 * \return A new Quat being the reciprocal of the quaternion (i.e. the
	 * conjungate divided by the square of the magnitude)
	 */
	Quat<T> reciprocal() const {
		return conjungate()/(S*S+V.magnitude2());
	}

	/**
	 * \brief Dot (inner) product of two quaternions
	 * \param rhs Right hand side quaternion
	 * \return Scalar result of dot product
	 */
	template<class T2=T>
	T dot(const Quat<T2> &rhs) const noexcept {
		static_assert(std::is_convertible<T2, T>::value,
				"element type is not convertible");
		return 0.5 * ((*this) * rhs.conjungate() + this->conjungate() * rhs).scalar();
	}

	/**
	 * \brief Cross product of two quaternions
	 * \param rhs Right hand side quaternion
	 * \return Quaternion result of dot product
	 */
	template<class T2=T>
	Quat<T> cross(Quat<T2> &rhs) const noexcept {
		static_assert(std::is_convertible<T2, T>::value,
				"element type is not convertible");
		return Quat<T>(0.5 * ((*this) * rhs - rhs * (*this)));
	}


	// Type Info
	/** \brief Elements Type */
	using element_type=T;

private:
	/** \brief instance scalar part storage */
	T S{0.};
	/** \brief instance vector part storage */
	Vec<T, 3> V{};
};

/**
 * \brief Multiplicative scaled
 * \tparam T Element type
 * \tparam T2 The scaling value type
 * \param lhs The scaling value
 * \param rhs The quaternion
 * \return A new Quat with each element being the product of the quaternion
 * and the scaling value
 * \pre The scaling value type has to be convertible to the quaternion
 * elements type
 */
template<class T=double, class T2=T, typename std::enable_if<std::is_scalar<T2>::value,
	bool>::type=true>
inline Quat<T> operator*(const T2 &lhs, const Quat<T> &rhs) noexcept {
	static_assert(std::is_convertible<T2, T>::value, "can not convert element types");
	return rhs*lhs;
}

/** \brief Shorthand for a Quat of type double */
using QuatD=Quat<double>;

}

#endif //NUMERICTYPES_QUAT_HPP
