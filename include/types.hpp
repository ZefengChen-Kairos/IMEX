#ifndef TYPES_HPP
#define TYPES_HPP

// Tensor / Grid / Array containers and the Eigen matrix-free operator glue.

#include "config.hpp"

struct Tensor
{
	double* data;
	int nxCells;
	int nyCells;
	int nRKS;
	
	// Access to data as l−value
	
	Tensor(int RKS,int xcell, int ycell){
		nxCells = xcell;
		nyCells = ycell;
		nRKS = RKS;
		data = new double[nRKS * (nxCells + 2 * ghost_cells) * (nyCells + 2 * ghost_cells) * NUM_VARS];
	}
	
	double operator()(int r,int i, int j, int k)const
	{
		return data[r * (NUM_VARS * ((nyCells + 2 * ghost_cells) * (nxCells + 2 * ghost_cells))) + NUM_VARS * (i + j * (nxCells + 2 * ghost_cells)) + k];
	}
	double& operator()(int r,int i, int j, int k){
		// return data[NUM_VARSS * (i + j * (nxCells + 2 * BCells)) + k];
		return data[r * (NUM_VARS * ((nyCells + 2 * ghost_cells) * (nxCells + 2 * ghost_cells))) + NUM_VARS * (i + j * (nxCells + 2 * ghost_cells)) + k];
	}
	~Tensor() {
		delete[] data;
	}
};

struct Grid
{
	double* data;
	int nxCells;
	int nyCells;
	
	// Access to data as l−value
	
	Grid(int xcell, int ycell){
		nxCells = xcell;
		nyCells = ycell;
		data = new double[(nxCells + 2 * ghost_cells) * (nyCells + 2 * ghost_cells) * NUM_VARS];
	}
	
	double operator()(int i, int j, int k)const
	{
		// return data[NUM_VARSS * (i + j * (nxCells + 2 * BCells)) + k];
		return data[NUM_VARS * (i + j * (nxCells + 2 * ghost_cells)) + k];
	}
	double& operator()(int i, int j, int k){
		// return data[NUM_VARSS * (i + j * (nxCells + 2 * BCells)) + k];
		return data[NUM_VARS * (i + j * (nxCells + 2 * ghost_cells)) + k];
	}
	~Grid() {
		delete[] data;
	}
};

struct Array
{
	double* data;
	int nxCells;
	int nyCells;
	
	// Access to data as l−value
	
	Array(int xcell, int ycell){
		nxCells = xcell;
		nyCells = ycell;
		data = new double[(nxCells + 2 * ghost_cells) * (nyCells + 2 * ghost_cells)];
	}
	
	double operator()(int i, int j)const
	{
		// return data[NUM_VARSS * (i + j * (nxCells + 2 * BCells)) + k];
		return data[i + j * (nxCells + 2 * ghost_cells)];
	}
	double& operator()(int i, int j){
		// return data[NUM_VARSS * (i + j * (nxCells + 2 * BCells)) + k];
		return data[i + j * (nxCells + 2 * ghost_cells)];
	}
	~Array() {
		delete[] data;
	}
};

inline Array *gFluxBy_p0 = nullptr;
inline Array *gFluxBx_0p = nullptr;
inline bool   g_ct_selfcheck = true;


class MatrixReplacement;
using Eigen::SparseMatrix;

namespace Eigen {
	namespace internal {
		// MatrixReplacement looks-like a SparseMatrix, so let's inherit its traits:
		template <>
		struct traits<MatrixReplacement> : public Eigen::internal::traits<Eigen::SparseMatrix<double> > {};
	}  // namespace internal
}  // namespace Eigen

// Example of a matrix-free wrapper from a user type to Eigen's compatible type
// For the sake of simplicity, this example simply wrap a Eigen::SparseMatrix.
class MatrixReplacement : public Eigen::EigenBase<MatrixReplacement> {
public:
	// Required typedefs, constants, and method:
	typedef double Scalar;
	typedef double RealScalar;
	typedef int StorageIndex;
	enum { ColsAtCompileTime = Eigen::Dynamic, MaxColsAtCompileTime = Eigen::Dynamic, IsRowMajor = false };
	
	Index rows() const { return mp_mat->rows(); }
	Index cols() const { return mp_mat->cols(); }
	
	template <typename Rhs>
	Eigen::Product<MatrixReplacement, Rhs, Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const {
		return Eigen::Product<MatrixReplacement, Rhs, Eigen::AliasFreeProduct>(*this, x.derived());
	}
	
	// Custom API:
	MatrixReplacement() : mp_mat(0) {}
	
	void attachMyMatrix(const SparseMatrix<double>& mat) { mp_mat = &mat; }
	const SparseMatrix<double> my_matrix() const { return *mp_mat; }
	
private:
	const SparseMatrix<double>* mp_mat;
};

// Implementation of MatrixReplacement * Eigen::DenseVector though a specialization of internal::generic_product_impl:
namespace Eigen {
	namespace internal {
		
		template <typename Rhs>
		struct generic_product_impl<MatrixReplacement, Rhs, SparseShape, DenseShape,
		GemvProduct>  // GEMV stands for matrix-vector
		: generic_product_impl_base<MatrixReplacement, Rhs, generic_product_impl<MatrixReplacement, Rhs> > {
			typedef typename Product<MatrixReplacement, Rhs>::Scalar Scalar;
			
			template <typename Dest>
			static void scaleAndAddTo(Dest& dst, const MatrixReplacement& lhs, const Rhs& rhs, const Scalar& alpha) {
				// This method should implement "dst += alpha * lhs * rhs" inplace,
				// however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
				eigen_assert(alpha == Scalar(1) && "scaling is not implemented");
				EIGEN_ONLY_USED_FOR_DEBUG(alpha);
				
				// Here we could simply call dst.noalias() += lhs.my_matrix() * rhs,
				// but let's do something fancier (and less efficient):
				for (Index i = 0; i < lhs.cols(); ++i) dst += rhs(i) * lhs.my_matrix().col(i);
			}
		};
		
	}  // namespace internal
}  // namespace Eigen

#endif // TYPES_HPP
