module;

#include <numbers>

export module xk.Math.Angles;
export import <compare>;

namespace xk::Math
{
	export template<class Ty>
	struct Radian;

	export template<class Ty>
		struct Degree
	{
		using value_type = Ty;
		Ty _value{};

		Degree() = default;
		Degree(value_type value) : _value{ value } {}
		explicit Degree(Radian<value_type> angle);

		auto operator<=>(const Degree&) const noexcept = default;

		Degree& operator=(const Degree& value) noexcept
		{
			_value = value._value;
			return *this;
		}

		Degree& operator+=(const Degree& rh) noexcept
		{
			_value += rh._value;
			return *this;
		}

		Degree& operator-=(const Degree& rh) noexcept
		{
			_value -= rh._value;
			return *this;
		}

		friend Degree operator+(Degree lh, const Degree& rh) noexcept
		{
			return lh += rh;
		}

		friend Degree operator-(Degree lh, const Degree& rh) noexcept
		{
			return lh -= rh;
		}
	};

	export template<class Ty>
		struct Radian
	{
		using value_type = Ty;
		Ty _value{};

		Radian() = default;
		Radian(value_type value) : _value{ value } {}
		explicit Radian(Degree<Ty> angle) :
			_value{ angle._value * static_cast<Ty>(std::numbers::pi) / static_cast<Ty>(180.0) }
		{
		}

		auto operator<=>(const Radian&) const noexcept = default;

		Radian& operator=(const Radian& value) noexcept
		{
			_value = value._value;
			return *this;
		}

		Radian& operator+=(const Radian& rh) noexcept
		{
			_value += rh._value;
			return *this;
		}

		Radian& operator-=(const Radian& rh) noexcept
		{
			_value -= rh._value;
			return *this;
		}

		friend Radian operator+(Radian lh, const Radian& rh) noexcept
		{
			return lh += rh;
		}

		friend Radian operator-(Radian lh, const Radian& rh) noexcept
		{
			return lh -= rh;
		}
	};

	template<class Ty>
	Degree<Ty>::Degree(Radian<Ty> angle) :
		_value{ angle._value / static_cast<Ty>(std::numbers::pi) * static_cast<Ty>(180.0) }
	{
	}
}