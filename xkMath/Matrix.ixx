module;

#include <array>
#include <algorithm>
#include <numeric>
#include <concepts>

export module xkMath.Matrix;

namespace xk
{
	export template<class Ty, size_t M, size_t N>
	struct Matrix;

	export template<class Ty, size_t M, size_t N, template<class Ty, size_t M, size_t N> class MatrixTy>
	class ColumnView
	{
		using value = Ty;
		using reference = Ty&;
		using const_reference = Ty&;

		using matrix = MatrixTy<Ty, M, N>;
		using matrix_pointer = matrix*;
		using matrix_reference = matrix&;

		matrix_pointer m_matrix;
		size_t m_selectedColumn;

	public:
		constexpr ColumnView(matrix_reference matrix, size_t column) :
			m_matrix(&matrix),
			m_selectedColumn(column)
		{

		}

	public:
		constexpr auto operator[](size_t index) const;
	};

	export template<class Ty, size_t M, size_t N, template<class Ty, size_t M, size_t N> class MatrixTy>
	class RowView
	{
		using value = Ty;
		using reference = Ty&;
		using const_reference = Ty&;

		using matrix = MatrixTy<Ty, M, N>;
		using matrix_pointer = matrix*;
		using matrix_reference = matrix&;

		matrix_pointer m_matrix;
		size_t m_selectedRow;

	public:
		constexpr RowView(matrix_reference matrix, size_t row) :
			m_matrix(&matrix),
			m_selectedRow(row)
		{

		}

	public:
		constexpr auto operator[](size_t index) const;
	};

	template<class Ty, class Ty2, size_t M, size_t N, size_t M2, template<class Ty, size_t M, size_t N> class MatrixTy>
	constexpr Ty operator*(RowView<Ty, M, N, MatrixTy> lh, ColumnView<Ty2, N, M2, MatrixTy> rh);

	template<class Ty, size_t M, size_t N>
	struct Matrix
	{
		static constexpr size_t row_count = M;
		static constexpr size_t column_count = N;
		static constexpr size_t element_count = row_count * column_count;
		static constexpr bool is_square_matrix = row_count == column_count;
		static constexpr bool is_vector = M > 0 && N == 1;
		using value = Ty;
		using reference = Ty&;
		using const_reference = Ty&;

		std::array<value, element_count> _values{};

	public:
		constexpr Matrix() = default;
		static constexpr Matrix Identity() requires (is_square_matrix)
		{
			Matrix mat;
			for (size_t i = 0; i < column_count; i++)
			{
				mat.At(i, i) = 1;
			}
			return mat;
		}

	public:
		constexpr reference At(size_t row, size_t column) { return _values[ColumnMajorIndex(row, column)]; }
		constexpr const_reference At(size_t row, size_t column) const { return _values[ColumnMajorIndex(row, column)]; }

		constexpr reference operator[](size_t index) { return _values[index]; }
		constexpr reference operator[](size_t index) const { return _values[index]; }

		constexpr reference X() requires (element_count >= 1) && (is_vector) { return _values[0]; }
		constexpr const_reference X() const requires (element_count >= 1) && (is_vector) { return _values[0]; }

		constexpr reference Y() requires (element_count >= 2) && (is_vector) { return _values[1]; }
		constexpr const_reference Y() const requires (element_count >= 2) && (is_vector) { return _values[1]; }

		constexpr reference Z() requires (element_count >= 3) && (is_vector) { return _values[2]; }
		constexpr const_reference Z() const requires (element_count >= 3) && (is_vector) { return _values[2]; }

		constexpr reference W() requires (element_count >= 4) && (is_vector) { return _values[3]; }
		constexpr const_reference W() const requires (element_count >= 4) && (is_vector) { return _values[3]; }

	public:
		template<class Ty2>
		constexpr Matrix& operator+=(const Matrix<Ty2, row_count, column_count>& rh)
		{
			for (size_t i = 0; i < element_count; i++)
			{
				_values[i] += rh._values[i];
			}

			return *this;
		}

		template<class Ty2>
		constexpr Matrix& operator-=(const Matrix<Ty2, row_count, column_count>& rh)
		{
			for (size_t i = 0; i < element_count; i++)
			{
				_values[i] -= rh._values[i];
			}

			return *this;
		}

		template<class Ty2, size_t M2>
		friend constexpr Matrix<Ty, M, M2> operator*(const Matrix<Ty, M, N>& lh, const Matrix<Ty2, N, M2>& rh)
		{
			Matrix<M, M2> result;
			for (size_t row = 0; row < result.row_count; row++)
			{
				for (size_t column = 0; column < result.column_count; column++)
				{
					result.At(row, column) = RowView(lh, row) * ColumnView(rh, column);
				}
			}
			return result;
		}

		template<class Ty2>
		Matrix<Ty, row_count, column_count>& operator*=(const Matrix<Ty2, row_count, column_count>& rh) requires(column_count == row_count)
		{
			this = (*this) * rh;
			return *this;
		}

		template<class Ty2>
		Matrix<Ty, row_count, column_count>& operator*=(Ty2 scalar)
		{
			for (reference value : _values)
			{
				value *= scalar;
			}
			return *this;
		}

		template<class Ty2>
		Matrix& operator/=(Ty2 scalar)
		{
			for (reference value : _values)
			{
				value /= scalar;
			}
			return *this;
		}

		template<class Ty2>
		friend Matrix& operator*(Matrix lh, Ty2 scalar)
		{
			return lh *= scalar;
		}

		template<class Ty2>
		friend Matrix& operator/(Matrix lh, Ty2 scalar)
		{
			return lh /= scalar;
		}

	private:
		constexpr size_t ColumnMajorIndex(size_t row, size_t column) { return column * row_count + row; }
	};

	export template<class Ty, size_t ElementCount>
	using Vector = Matrix<Ty, ElementCount, 1>;
	
	template<class Ty, class Ty2, size_t M, size_t N, size_t M2, template<class Ty, size_t M, size_t N> class MatrixTy>
	constexpr Ty operator*(RowView<Ty, M, N, MatrixTy> lh, ColumnView<Ty2, N, M2, MatrixTy> rh)
	{
		Ty value = 0;
		for (size_t i = 0; i < N; i++)
		{
			value += lh[i] * rh[i];
		}
		return value;
	}

	template<class Ty, size_t M, size_t N, template<class Ty, size_t M, size_t N> class MatrixTy>
	constexpr auto ColumnView<Ty, M, N, MatrixTy>::operator[](size_t index) const
	{
		return (m_matrix->At(index, m_selectedColumn));
	}

	template<class Ty, size_t M, size_t N, template<class Ty, size_t M, size_t N> class MatrixTy>
	constexpr auto RowView<Ty, M, N, MatrixTy>::operator[](size_t index) const
	{
		return (m_matrix->At(m_selectedRow, index));
	}

	template<class Ty, size_t M, size_t N>
	Matrix<Ty, N, M> Transpose(const Matrix<Ty, M, N>& mat)
	{
		Matrix<Ty, N, M> result;
		for (size_t row = 0; row < mat.row_count; row++)
		{
			for (size_t column = 0; column < mat.column_count; column++)
			{
				result.At(column, row) = mat.At(row, column);
			}
		}
		return result;
	}

	template<class Ty, size_t ElementCount>
	Ty Dot(const Vector<Ty, ElementCount>& lh, const Vector<Ty, ElementCount>& rh)
	{
		Ty value = 0;
		for (size_t i = 0; i < ElementCount; i++)
		{
			value += lh[i] * rh[i];
		}
		return value;
	}
}
