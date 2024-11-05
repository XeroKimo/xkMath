module;

#include <concepts>
#include <cmath>
#include <algorithm>
#include <numbers>

export module xk.Math.Quaternion;
import xk.Math.Matrix;
import xk.Math.Angles;

namespace xk::Math
{
	export template<std::floating_point T>
	struct Quaternion
	{
	public:
		using value_type = T;
		using data_type = Vector<T, 4>;
		using axis_type = Vector<T, 3>;

	private:
		data_type data;

	public:
		constexpr Quaternion() :
			data{ 0, 0, 0, 1 }
		{

		}
		constexpr Quaternion(value_type x, value_type y, value_type z, value_type w) :
			data{ x, y, z, w }
		{

		}
		constexpr Quaternion(Degree<value_type> x, Degree<value_type> y, Degree<value_type> z) :
			Quaternion(Radian{ x }, Radian{ y }, Radian{ z })
		{
		}

		//Referenced equations https://www.euclideanspace.com/maths/geometry/rotations/conversions/eulerToQuaternion/index.htm
		constexpr Quaternion(Radian<value_type> x, Radian<value_type> y, Radian<value_type> z)
		{
			T cx = std::cos(x._value / 2);
			T cy = std::cos(y._value / 2);
			T cz = std::cos(z._value / 2);
			T sx = std::sin(x._value / 2);
			T sy = std::sin(y._value / 2);
			T sz = std::sin(z._value / 2);

			this->W() = (cx * cy * cz) - (sx * sy * sz);
			this->X() = (cx * sy * sz) + (sx * cy * cz);
			this->Y() = (cx * sy * cz) + (sx * cy * sz);
			this->Z() = (cx * cy * sz) - (sx * sy * cz);
		}

		constexpr Quaternion(axis_type axis, Degree<value_type> angle) :
			Quaternion(axis, Radian{ angle })
		{

		}
		constexpr Quaternion(axis_type axis, Radian<value_type> angle)
		{
			angle._value /= 2;
			value_type sinAngle = sin(angle._value);
			W() = std::cos(angle._value);
			X() = axis.X() * sinAngle;
			Y() = axis.Y() * sinAngle;
			Z() = axis.Z() * sinAngle;
		}

	public:
		constexpr Quaternion& operator+=(const Quaternion& rh)
		{
			std::transform(data.begin(), data.end(), rh.data.begin(), data.begin(), std::plus{});
			return *this;
		}
		constexpr Quaternion& operator-=(const Quaternion& rh)
		{
			std::transform(data.begin(), data.end(), rh.data.begin(), data.begin(), std::minus{});
			return *this;
		}
		constexpr Quaternion& operator*=(const Quaternion& rh)
		{
			data = rh.ToMatrix() * data;

			return *this;
		}
		constexpr Quaternion operator-() const
		{
			Quaternion copy;
			copy.data = data;
			copy.w() *= -1;
			return copy;
		}
		//constexpr Quaternion& operator/=(const Quaternion& rh)
		//{
		//    //std::transform(data.begin(), data.end(), rh.data.begin(), data.begin(), std::minus{});
		//    return *this;
		//}

		friend constexpr Quaternion operator+(Quaternion lh, const Quaternion& rh)
		{
			return (lh += rh);
		}
		friend constexpr Quaternion operator-(Quaternion lh, const Quaternion& rh)
		{
			return (lh -= rh);
		}
		friend constexpr Quaternion operator*(Quaternion lh, const Quaternion& rh)
		{
			return (lh *= rh);
		}
		//friend constexpr Quaternion operator/(Quaternion lh, const Quaternion& rh)
		//{
		//    return (lh /= rh);
		//}

		friend constexpr bool operator==(const Quaternion& lh, const Quaternion& rh) = default;
		friend constexpr bool operator!=(const Quaternion& lh, const Quaternion& rh) = default;

		constexpr value_type& operator[](size_t index) { return data[index]; }
		constexpr const value_type& operator[](size_t index) const { return data[index]; }

	public:
		constexpr value_type& X() { return data[0]; }
		constexpr value_type& Y() { return data[1]; }
		constexpr value_type& Z() { return data[2]; }
		constexpr value_type& W() { return data[3]; }

		constexpr const value_type& X() const { return data[0]; }
		constexpr const value_type& Y() const { return data[1]; }
		constexpr const value_type& Z() const { return data[2]; }
		constexpr const value_type& W() const { return data[3]; }

		constexpr auto begin() { return data.begin(); }
		constexpr auto rbegin() { return data.rbegin(); }
		constexpr auto end() { return data.end(); }
		constexpr auto rend() { return data.rend(); }

		constexpr auto begin() const { return data.begin(); }
		constexpr auto rbegin() const { return data.rbegin(); }
		constexpr auto end() const { return data.end(); }
		constexpr auto rend() const { return data.rend(); }

	public:
		void Normalize()
		{
			Normalize(data);
		}

		Vector<Degree<value_type>, 3> ToEulerDegree() const
		{
			auto Radian = ToEulerRadian();

			return
			{
				 Degree{ Radian.X() },
				 Degree{ Radian.Y() },
				 Degree{ Radian.Z() }
			};
		}

		//Referenced equations https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
		Vector<Radian<value_type>, 3> ToEulerRadian() const
		{
			value_type yValue = 2 * (W() * Y() - Z() * X());
			return
			{
				Radian{ std::atan2(2 * (W() * X() + Y() * Z()), 1 - 2 * (X() * X() + Y() * Y()))                       },
				Radian{ (std::abs(yValue) >= 1) ? std::copysign(std::numbers::pi_v<T> / 2, -yValue) : std::asin(yValue)  },
				Radian{ std::atan2(2 * (W() * Z() + Y() * X()),  1 - 2 * (Z() * Z() + Y() * Y()))                        },
			};
		}

		//Referenced equations https://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/code/index.htm
		//The matrix form of the quaternion, if you're looking to make an actually rotation matrix, use ToRotationMatrix()
		Matrix<value_type, 4, 4> ToMatrix() const
		{
			return
			{
				 W(),  Z(), -Y(),  X(),
				-Z(),  W(),  X(),  Y(),
				 Y(), -X(),  W(),  Z(),
				-X(), -Y(), -Z(),  W()
			};
		}


		//Referenced equations https://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToMatrix/index.htm
		//Creates a Rotation Matrix from the quaternion, if you're looking to view the matrix form of the quaternion, use ToMatrix()
		Matrix<value_type, 4, 4> ToRotationMatrix() const
		{
			Matrix<float, 4, 4> lh =
			{
				 W(), -Z(),  Y(), -X(),
				 Z(),  W(), -X(), -Y(),
				-Y(),  X(),  W(), -Z(),
				 X(),  Y(),  Z(), W()
			};
			Matrix<float, 4, 4> rh =
			{
				 W(), -Z(),  Y(),  X(),
				 Z(),  W(), -X(),  Y(),
				-Y(),  X(),  W(),  Z(),
				-X(), -Y(), -Z(),  W()
			};

			return lh * rh;
		}
	};
}