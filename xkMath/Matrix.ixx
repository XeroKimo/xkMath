module;

#include <array>
#include <algorithm>
#include <numeric>
#include <concepts>
#include <cstdint>
#include <cmath>
#include <functional>
#include <cassert>

#pragma warning(disable:4244)
export module xk.Math.Matrix;
import xk.Math.Angles;

//Testing resyncing across submodules

namespace xk::Math
{
	export template<class Ty, size_t M, size_t N>
		struct Matrix;

	export template<class Ty, size_t M, size_t N>
		constexpr Matrix<Ty, N, M> Transpose(const Matrix<Ty, M, N>& mat) noexcept;

	template<class Ty, bool IsConst>
	struct BackingMatrixRef
	{
		using matrix = std::remove_const_t<Ty>;
		using matrix_pointer = std::remove_const_t<Ty>*;
		using matrix_reference = std::remove_const_t<Ty>&;
		using reference = matrix::reference;

		matrix_pointer m_matrix;
	};

	template<class Ty>
	struct BackingMatrixRef<Ty, true>
	{
		using matrix = const std::remove_const_t<Ty>;
		using matrix_pointer = const std::remove_const_t<Ty>*;
		using matrix_reference = const std::remove_const_t<Ty>&;
		using reference = matrix::const_reference;

		matrix_pointer m_matrix;
	};

	export template<class Ty, size_t M, size_t N, bool IsConst>
		class ColumnRef : private BackingMatrixRef<Matrix<Ty, M, N>, IsConst>
	{
		using base_type = BackingMatrixRef<Matrix<Ty, M, N>, IsConst>;
		using matrix = base_type::matrix;
		using matrix_pointer = base_type::matrix_pointer;
		using matrix_reference = base_type::matrix_reference;
		using reference = base_type::reference;

		using base_type::m_matrix;
		size_t m_selectedColumn;

	public:
		constexpr ColumnRef(matrix_reference matrix, size_t column) :
			base_type(&matrix),
			m_selectedColumn(column)
		{

		}

	public:
		constexpr reference operator[](size_t index) const;
	};

	export template<class Ty, size_t M, size_t N, bool IsConst>
		class RowRef : private BackingMatrixRef<Matrix<Ty, M, N>, IsConst>
	{
		using base_type = BackingMatrixRef<Matrix<Ty, M, N>, IsConst>;
		using matrix = base_type::matrix;
		using matrix_pointer = base_type::matrix_pointer;
		using matrix_reference = base_type::matrix_reference;

		using reference = base_type::reference;
		using base_type::m_matrix;
		size_t m_selectedRow;

	public:
		constexpr RowRef(matrix_reference matrix, size_t row) :
			base_type(&matrix),
			m_selectedRow(row)
		{

		}

	public:
		constexpr reference operator[](size_t index) const;
	};

	template<class Ty, size_t M, size_t N>
	ColumnRef(Matrix<Ty, M, N>&, size_t) -> ColumnRef<Ty, M, N, false>;

	template<class Ty, size_t M, size_t N>
	ColumnRef(const Matrix<Ty, M, N>&, size_t) -> ColumnRef<Ty, M, N, true>;

	template<class Ty, size_t M, size_t N>
	RowRef(Matrix<Ty, M, N>&, size_t) -> RowRef<Ty, M, N, false>;

	template<class Ty, size_t M, size_t N>
	RowRef(const Matrix<Ty, M, N>&, size_t) -> RowRef<Ty, M, N, true>;

	template<class Ty1, class Ty2, size_t M, size_t N, size_t M2, bool IsConst1, bool IsConst2>
	constexpr auto operator*(RowRef<Ty1, M, N, IsConst1> lh, ColumnRef<Ty2, N, M2, IsConst2> rh);

	export enum class MatrixMemoryLayout
	{
		Row_Major,
		Column_Major
	};

	template<class Ty, size_t M, size_t N>
	struct Matrix
	{
		static constexpr size_t row_count = M;
		static constexpr size_t column_count = N;
		static constexpr size_t element_count = row_count * column_count;
		static constexpr bool is_square_matrix = row_count == column_count;
		using value_type = Ty;
		using reference = Ty&;
		using const_reference = const Ty&;

		std::array<value_type, element_count> _values{};
		static constexpr MatrixMemoryLayout nativeLayout = MatrixMemoryLayout::Row_Major;
		static constexpr MatrixMemoryLayout memoryLayout = MatrixMemoryLayout::Column_Major;
	public:
		constexpr Matrix() = default;

		template<std::convertible_to<Ty> Ty2>
		constexpr Matrix(const Matrix<Ty2, M, N>& other)
		{
			std::copy(other._values.begin(), other._values.end(), _values.begin());
		}

		template<std::convertible_to<Ty>... Ty2>
			requires (sizeof...(Ty2) <= element_count)
		constexpr Matrix(Ty2... values) :
			_values({ static_cast<Ty>(values)... })
		{
			if constexpr(memoryLayout != nativeLayout)
			{
				auto copy = _values;
				for(size_t i = 0; i < row_count; i++)
				{
					for(size_t j = 0; j < column_count; j++)
					{
						_values[GetIndex(i, j)] = copy[NativeMemoryIndex(i, j)];
					}
				}
			}
		}

		static constexpr Matrix Identity() requires (is_square_matrix)
		{
			Matrix mat;
			for(size_t i = 0; i < column_count; i++)
			{
				mat.At(i, i) = 1;
			}
			return mat;
		}

	public:
		constexpr reference At(size_t row, size_t column) { return _values[GetIndex(row, column)]; }
		constexpr const_reference At(size_t row, size_t column) const { return _values[GetIndex(row, column)]; }

	public:
		template<class Ty2>
		constexpr Matrix& operator+=(const Matrix<Ty2, row_count, column_count>& rh)
		{
			for(size_t i = 0; i < element_count; i++)
			{
				_values[i] += rh._values[i];
			}

			return *this;
		}

		template<class Ty2>
		constexpr Matrix& operator-=(const Matrix<Ty2, row_count, column_count>& rh)
		{
			for(size_t i = 0; i < element_count; i++)
			{
				_values[i] -= rh._values[i];
			}

			return *this;
		}

		constexpr Matrix operator-() const
		{
			Matrix temp = *this;
			temp *= -1;
			return temp;
		}

		template<class Ty2, size_t M2>
		friend constexpr Matrix<decltype(std::declval<Ty>()* std::declval<Ty2>()), M, M2> operator*(const Matrix<Ty, M, N>& lh, const Matrix<Ty2, N, M2>& rh)
		{
			Matrix<decltype(std::declval<Ty>()* std::declval<Ty2>()), M, M2> result;
			for(size_t row = 0; row < result.row_count; row++)
			{
				for(size_t column = 0; column < result.column_count; column++)
				{
					result.At(row, column) = RowRef(lh, row) * ColumnRef(rh, column);
				}
			}
			return result;
		}

		template<class Ty2>
		friend constexpr Matrix operator+(Matrix<Ty2, row_count, column_count> lh, const Matrix<Ty2, row_count, column_count>& rh)
		{
			return lh += rh;
		}

		template<class Ty2>
		friend constexpr Matrix operator-(Matrix<Ty2, row_count, column_count> lh, const Matrix<Ty2, row_count, column_count>& rh)
		{
			return lh -= rh;
		}

		template<class Ty2, size_t M2, size_t N2>
		constexpr Matrix& operator*=(const Matrix<Ty2, M2, N2>& rh) requires (M2 == N) && (N2 == N)
		{
			*this = (*this) * rh;
			return *this;
		}

		template<class Ty2>
			requires std::is_arithmetic_v<Ty2>
		constexpr Matrix& operator*=(Ty2 scalar)
		{
			for(reference value : _values)
			{
				value *= scalar;
			}
			return *this;
		}

		template<class Ty2>
			requires std::is_arithmetic_v<Ty2>
		constexpr Matrix& operator/=(Ty2 scalar)
		{
			for(reference value : _values)
			{
				value /= scalar;
			}
			return *this;
		}

		template<class Ty2>
			requires std::is_arithmetic_v<Ty2>
		friend constexpr Matrix operator*(Matrix lh, Ty2 scalar)
		{
			return lh *= scalar;
		}

		template<class Ty2>
			requires std::is_arithmetic_v<Ty2>
		friend constexpr Matrix operator*(Ty2 scalar, Matrix lh)
		{
			return lh *= scalar;
		}

		template<class Ty2>
			requires std::is_arithmetic_v<Ty2>
		friend constexpr Matrix operator/(Matrix lh, Ty2 scalar)
		{
			return lh /= scalar;
		}

		template<class Ty2>
		friend constexpr bool operator==(const Matrix& lh, const Matrix<Ty2, M, N>& rh)
		{
			return std::equal(lh._values.begin(), lh._values.end(), rh._values.begin());
		}

		operator Ty() const requires (element_count == 1) { return _values[0]; }

		constexpr size_t GetIndex(size_t row, size_t column) const noexcept
		{
			assert(row < row_count);
			assert(column < column_count);
			return (memoryLayout == MatrixMemoryLayout::Column_Major) ?
				ColumnMajorIndex(row, column) :
				RowMajorIndex(row, column);
		}

	private:
		constexpr size_t NativeMemoryIndex(size_t row, size_t column) const noexcept
		{
			return (nativeLayout == MatrixMemoryLayout::Column_Major) ?
				ColumnMajorIndex(row, column) :
				RowMajorIndex(row, column);
		}
		constexpr size_t ColumnMajorIndex(size_t row, size_t column) const noexcept { return column * row_count + row; }
		constexpr size_t RowMajorIndex(size_t row, size_t column) const noexcept { return row * column_count + column; }
	};

	export template<class Ty, size_t ElementCount>
		using SquareMatrix = Matrix<Ty, ElementCount, ElementCount>;

	export template<class Ty, size_t ElementCount>
		struct Vector;

	export template<class Ty, size_t ElementCount>
		constexpr Matrix<Ty, ElementCount, 1> ToColumnMatrix(Vector<Ty, ElementCount> vector)
	{
		return[vector]<size_t... Is>(std::index_sequence<Is...>)
		{
			return Matrix<Ty, ElementCount, 1>{ vector[Is]... };
		}(std::make_index_sequence<ElementCount>{});
	}

	export template<class Ty, size_t ElementCount>
		constexpr Matrix<Ty, 1, ElementCount> ToRowMatrix(Vector<Ty, ElementCount> vector)
	{
		return[vector]<size_t... Is>(std::index_sequence<Is...>)
		{
			return Matrix<Ty, 1, ElementCount>{ vector[Is]... };
		}(std::make_index_sequence<ElementCount>{});
	}

	export template<class Ty>
		struct Uniform
	{
		Ty value;
	};

	export template<class Ty, size_t ElementCount>
		struct Vector
	{
		static constexpr size_t element_count = ElementCount;
		using value_type = Ty;
		using reference = Ty&;
		using const_reference = const Ty&;

		std::array<value_type, ElementCount> _values{};

	public:
		constexpr Vector() = default;

		template<std::convertible_to<Ty> Ty2>
		constexpr Vector(Uniform<Ty2> value)
		{
			_values.fill(value.value);
		}

		constexpr Vector(Matrix<Ty, ElementCount, 1> values)
		{
			std::copy(values._values.begin(), values._values.end(), _values.begin());
		}

		constexpr Vector(Matrix<Ty, 1, ElementCount> values)
		{
			std::copy(values._values.begin(), values._values.end(), _values.begin());
		}

		template<std::convertible_to<Ty>... OtherTy2>
			requires (sizeof...(OtherTy2) <= element_count)
		constexpr Vector(OtherTy2... values) :
			_values{ static_cast<Ty>(values)... }
		{

		}

		template<std::convertible_to<Ty> OtherTy, size_t OtherElementCount, std::convertible_to<Ty>... OtherTy2>
			requires (sizeof...(OtherTy2) + OtherElementCount <= element_count)
		constexpr Vector(Vector<OtherTy, OtherElementCount> v, OtherTy2... values)
		{
			std::copy(v._values.begin(), v._values.end(), _values.begin());
			if constexpr(sizeof...(OtherTy2) > 0)
			{
				std::array<value_type, sizeof...(OtherTy2)> vs{ values... };
				std::copy(vs.begin(), vs.end(), _values.begin() + OtherElementCount);
			}
		}

		template<class Ty2>
		constexpr Vector& operator+=(const Vector<Ty2, element_count>& rh)
		{
			std::transform(_values.begin(), _values.end(), rh._values.begin(), _values.begin(), std::plus{});
			return *this;
		}

		template<class Ty2>
		constexpr Vector& operator-=(const Vector<Ty2, element_count>& rh)
		{
			std::transform(_values.begin(), _values.end(), rh._values.begin(), _values.begin(), std::minus{});
			return *this;
		}

		constexpr Vector operator-() const
		{
			Vector temp = *this;
			temp *= -1;
			return temp;
		}

		template<class Ty2>
		friend constexpr Vector operator+(Vector<Ty2, element_count> lh, const Vector<Ty2, element_count>& rh)
		{
			return lh += rh;
		}

		template<class Ty2>
		friend constexpr Vector operator-(Vector<Ty2, element_count> lh, const Vector<Ty2, element_count>& rh)
		{
			return lh -= rh;
		}

		template<class Ty2>
			requires std::is_arithmetic_v<Ty2>
		constexpr Vector& operator*=(Ty2 scalar)
		{
			std::transform(_values.begin(), _values.end(), _values.begin(), [scalar](value_type value) { return value *= scalar; });
			return *this;
		}

		template<class Ty2>
			requires std::is_arithmetic_v<Ty2>
		constexpr Vector& operator/=(Ty2 scalar)
		{
			std::transform(_values.begin(), _values.end(), _values.begin(), [scalar](value_type value) { return value /= scalar; });
			return *this;
		}

		template<class Ty2>
			requires std::is_arithmetic_v<Ty2>
		friend constexpr Vector operator*(Vector lh, Ty2 scalar)
		{
			return lh *= scalar;
		}

		template<class Ty2>
			requires std::is_arithmetic_v<Ty2>
		friend constexpr Vector operator*(Ty2 scalar, Vector lh)
		{
			return lh *= scalar;
		}

		template<class Ty2, size_t N>
			requires std::is_arithmetic_v<Ty2>
		friend constexpr Matrix<std::common_type_t<Ty, Ty2>, 1, N> operator*(Vector lh, Matrix<Ty2, element_count, N> rh)
		{
			return ToRowMatrix(lh) * rh;
		}

		template<class Ty2, size_t M>
			requires std::is_arithmetic_v<Ty2>
		friend constexpr Matrix<std::common_type_t<Ty, Ty2>, M, 1> operator*(Matrix<Ty2, M, element_count> lh, Vector rh)
		{
			return lh * ToColumnMatrix(rh);
		}

		template<class Ty2>
			requires std::is_arithmetic_v<Ty2>
		friend constexpr Vector operator/(Vector lh, Ty2 scalar)
		{
			return lh /= scalar;
		}

		template<class Ty2>
		friend constexpr bool operator==(const Vector& lh, const Vector<Ty2, element_count>& rh)
		{
			return std::equal(lh._values.begin(), lh._values.end(), rh._values.begin());
		}

		template<std::convertible_to<Ty> Ty2>
		Vector& operator=(Uniform<Ty2> value)
		{
			_values.fill(value.value);
			return *this;
		}

		auto begin() { return _values.begin(); }
		auto end() { return _values.end(); }

		auto begin() const { return _values.begin(); }
		auto end() const { return _values.end(); }

		operator Ty() const requires (element_count == 1) { return _values[0]; }

		constexpr reference operator[](size_t index) { return _values[index]; }
		constexpr const_reference operator[](size_t index) const { return _values[index]; }

		constexpr reference X() requires (element_count >= 1) { return _values[0]; }
		constexpr const_reference X() const requires (element_count >= 1) { return _values[0]; }

		constexpr reference Y() requires (element_count >= 2) { return _values[1]; }
		constexpr const_reference Y() const requires (element_count >= 2) { return _values[1]; }

		constexpr reference Z() requires (element_count >= 3) { return _values[2]; }
		constexpr const_reference Z() const requires (element_count >= 3) { return _values[2]; }

		constexpr reference W() requires (element_count >= 4) { return _values[3]; }
		constexpr const_reference W() const requires (element_count >= 4) { return _values[3]; }

		template<size_t... Index>
			requires ((Index < element_count) && ...) && (sizeof...(Index) < element_count)
		constexpr Vector<Ty, sizeof...(Index)> Swizzle()
		{
			return { _values[Index]... };
		}
	};

	template<class... Ty>
	Vector(Ty...) -> Vector<std::tuple_element_t<0, std::tuple<Ty...>>, sizeof...(Ty)>;

	template<class Ty1, class Ty2, size_t M, size_t N, size_t M2, bool IsConst1, bool IsConst2>
	constexpr auto operator*(RowRef<Ty1, M, N, IsConst1> lh, ColumnRef<Ty2, N, M2, IsConst2> rh)
	{
		using value_type = decltype(std::declval<Ty1>()* std::declval<Ty2>());
		value_type value = 0;
		for(size_t i = 0; i < N; i++)
		{
			value += lh[i] * rh[i];
		}
		return value;
	}

	export template<class Ty, size_t M, size_t N, bool IsConst>
		constexpr ColumnRef<Ty, M, N, IsConst>::reference ColumnRef<Ty, M, N, IsConst>::operator[](size_t index) const
	{
		return (m_matrix->At(index, m_selectedColumn));
	}

	export template<class Ty, size_t M, size_t N, bool IsConst>
		constexpr RowRef<Ty, M, N, IsConst>::reference RowRef<Ty, M, N, IsConst>::operator[](size_t index) const
	{
		return (m_matrix->At(m_selectedRow, index));
	}

	template<class Ty, size_t M, size_t N>
	constexpr Matrix<Ty, N, M> Transpose(const Matrix<Ty, M, N>& mat) noexcept
	{
		Matrix<Ty, N, M> result;
		for(size_t row = 0; row < mat.row_count; row++)
		{
			for(size_t column = 0; column < mat.column_count; column++)
			{
				result.At(column, row) = mat.At(row, column);
			}
		}
		return result;
	}

	export template<class Ty, size_t ElementCount>
		constexpr Ty Dot(const Vector<Ty, ElementCount>& lh, const Vector<Ty, ElementCount>& rh)
	{
		return std::transform_reduce(lh.begin(), lh.end(), rh.begin(), Ty{});
	}

	export template<class Ty, size_t ElementCount>
		constexpr Ty MagnitudeSquared(const Vector<Ty, ElementCount>& v)
	{
		return std::accumulate(v.begin(), v.end(), static_cast<Ty>(0), [](Ty total, Ty val)
		{
			return total + val * val;
		});
	}

	export template<class Ty, size_t ElementCount>
		constexpr Ty Magnitude(const Vector<Ty, ElementCount>& v)
	{
		return std::sqrt(MagnitudeSquared(v));
	}

	export template<class Ty, size_t ElementCount>
		constexpr Vector<Ty, ElementCount> Normalize(Vector<Ty, ElementCount> v)
	{
		std::transform(v.begin(), v.end(), v.begin(), [dividor = Magnitude(v)](auto& e)
		{
			return e /= dividor;
		});
		return v;
	}

	export template<class Ty>
		constexpr Matrix<Ty, 4, 4> TransformMatrix(const Vector<Ty, 3>& vector)
	{
		return
		{
			1, 0, 0, vector.X(),
			0, 1, 0, vector.Y(),
			0, 0, 1, vector.Z(),
			0, 0, 0, 1
		};
	}

	export template<class Ty>
		constexpr Matrix<Ty, 4, 4> ScaleMatrix(const Vector<Ty, 3>& vector)
	{
		return
		{
			vector.X(), 0, 0, 0,
			0, vector.Y(), 0, 0,
			0, 0, vector.Z(), 0,
			0, 0, 0, 1
		};
	}

	export template<class Ty, class Ty2, std::size_t ElementCount>
	constexpr Vector<Ty, ElementCount> HadamardProduct(Vector<Ty, ElementCount> lh, const Vector<Ty2, ElementCount>& rh)
	{
		std::transform(lh.begin(), lh.end(), rh.begin(), lh.begin(), std::multiplies{});
		return lh;
	}

	export template<class Ty, class Ty2, std::size_t ElementCount>
	constexpr Vector<Ty, ElementCount> HadamardDivision(Vector<Ty, ElementCount> lh, const Vector<Ty2, ElementCount>& rh)
	{
		std::transform(lh.begin(), lh.end(), rh.begin(), lh.begin(), std::divides{});
		return lh;
	}

	export template<class Ty, class Ty2, std::size_t ElementCount>
	constexpr Vector<Ty, ElementCount> HadamardSafeDivision(Vector<Ty, ElementCount> lh, const Vector<Ty2, ElementCount>& rh)
	{
		std::transform(lh.begin(), lh.end(), rh.begin(), lh.begin(), [](auto lh, auto rh) { return rh == 0 ? lh : lh / rh; });
		return lh;
	}

	export constexpr Matrix<float, 4, 4> OrthographicProjectionLH(Vector<float, 2> resolution, float zNear, float zFar)
	{
		return
		{
			2 / resolution.X(), 0,      0,                              0,
			0,      2 / resolution.Y(), 0,                              0,
			0,      0,     1 / (zFar - zNear),          -zNear / (zFar - zNear),
			0,      0,      0 , 1
		};
	}

	export constexpr Matrix<float, 4, 4> OrthographicProjectionAspectRatioLH(Vector<float, 2> aspectRatio, float viewSize, float zNear, float zFar)
	{
		return OrthographicProjectionLH({ aspectRatio.X() / aspectRatio.Y() * viewSize, viewSize }, zNear, zFar);
	}

	export constexpr Matrix<float, 4, 4> OrthographicProjectionRH(Vector<float, 2> resolution, float zNear, float zFar)
	{
		return
		{
			2 / resolution.X(), 0,      0,                              0,
			0,      2 / resolution.Y(), 0,                              0,
			0,      0,     1 / (zNear - zFar),          zNear / (zNear - zFar),
			0,      0,      0 , 1
		};
	}

	export constexpr Matrix<float, 4, 4> OrthographicProjectionAspectRatioRH(Vector<float, 2> aspectRatio, float viewSize, float zNear, float zFar)
	{
		return OrthographicProjectionRH({ aspectRatio.X() / aspectRatio.Y() * viewSize, viewSize }, zNear, zFar);
	}

	export Matrix<float, 4, 4> PerspectiveProjectionLH(Radian<float> fovY, float aspectRatio, float zNear, float zFar)
	{
		float yScale = 1 / tanf(fovY._value / 2);
		float xScale = yScale / aspectRatio;

		return
		{
			xScale, 0,      0,                              0,
			0,      yScale, 0,                              0,
			0,      0,      zFar / (zFar - zNear),          -zNear * zFar / (zFar - zNear),
			0,      0,      1 , 0
		};
	}

	export Matrix<float, 4, 4> PerspectiveProjectionRH(Radian<float> fovY, float aspectRatio, float zNear, float zFar)
	{
		float yScale = 1 / tanf(fovY._value / 2);
		float xScale = yScale / aspectRatio;

		return
		{
			xScale, 0,      0,                              0,
			0,      yScale, 0,                              0,
			0,      0,      zFar / (zNear - zFar),          zNear * zFar / (zNear - zFar),
			0,      0,      -1 , 0
		};
	}

	export Matrix<float, 4, 4>  PerspectiveProjectionLH(Degree<float> fovY, float aspectRatio, float zNear, float zFar)
	{
		return PerspectiveProjectionLH(Radian{ fovY }, aspectRatio, zNear, zFar);
	}

	export Matrix<float, 4, 4>  PerspectiveProjectionRH(Degree<float> fovY, float aspectRatio, float zNear, float zFar)
	{
		return PerspectiveProjectionRH(Radian{ fovY }, aspectRatio, zNear, zFar);
	}

	export template<class Ty>
	constexpr Vector<Ty, 3> Cross(Vector<Ty, 3> lh, Vector<Ty, 3> rh)
	{
		return
		{
			lh.Y() * rh.Z() - lh.Z() * rh.Y(),
			lh.Z() * rh.X() - lh.X() * rh.Z(),
			lh.X() * rh.Y() - lh.Y() * rh.X()
		};
	}

	export SquareMatrix<float, 4> LookAt(Vector<float, 3> position, Vector<float, 3> target, Vector<float, 3> upVector)
	{
		auto direction = Normalize(target - position);
		auto right = Normalize(Cross(upVector, direction));
		auto up = Normalize(Cross(direction, right));

		return SquareMatrix<float, 4>
		{
			right.X(), right                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 .Y(), right.Z(), 0,
			up.X(), up.Y(), up.Z(), 0,
			direction.X(), direction.Y(), direction.Z(), 0,
			0, 0, 0, 1
		} * TransformMatrix(-position);
	}

	export namespace Aliases
	{
		using Vector2 = Vector<float, 2>;
		using Vector3 = Vector<float, 3>;
		using Vector4 = Vector<float, 4>;

		using iVector2 = Vector<int, 2>;
		using iVector3 = Vector<int, 3>;
		using iVector4 = Vector<int, 4>;

		using uVector2 = Vector<unsigned int, 2>;
		using uVector3 = Vector<unsigned int, 3>;
		using uVector4 = Vector<unsigned int, 4>;

		using dVector2 = Vector<double, 2>;
		using dVector3 = Vector<double, 3>;
		using dVector4 = Vector<double, 4>;

		using i8Vector2 = Vector<std::int8_t, 2>;
		using i8Vector3 = Vector<std::int8_t, 3>;
		using i8Vector4 = Vector<std::int8_t, 4>;

		using u8Vector2 = Vector<std::uint8_t, 2>;
		using u8Vector3 = Vector<std::uint8_t, 3>;
		using u8Vector4 = Vector<std::uint8_t, 4>;

		using i16Vector2 = Vector<std::int16_t, 2>;
		using i16Vector3 = Vector<std::int16_t, 3>;
		using i16Vector4 = Vector<std::int16_t, 4>;

		using u16Vector2 = Vector<std::uint16_t, 2>;
		using u16Vector3 = Vector<std::uint16_t, 3>;
		using u16Vector4 = Vector<std::uint16_t, 4>;

		using Matrix2x2 = Matrix<float, 2, 2>;
		using Matrix3x3 = Matrix<float, 3, 3>;
		using Matrix3x4 = Matrix<float, 3, 4>;
		using Matrix4x4 = Matrix<float, 4, 4>;
	}
}

#pragma warning(default:4244)