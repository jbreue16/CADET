// =============================================================================
//  CADET
//  
//  Copyright © 2008-2022: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

/**
 * @file 
 * Implements the High Resolution Koren method
 */

#ifndef LIBCADET_HIGHRESKOREN_HPP_
#define LIBCADET_HIGHRESKOREN_HPP_

#include "AutoDiff.hpp"
#include "MathUtil.hpp"
#include "Memory.hpp"
#include "common/CompilerSpecific.hpp"
#include "cadet/Exceptions.hpp"

#include <algorithm>

namespace cadet
{

/**
 * @brief Implements the High Resolution Koren scheme for convection
 * @details This scheme is based on upwind stencils and provides WENO methods 1-1, 2-3, and 3-5.
 *          In general, WENO achieves order \f$ r \f$ using a stencil with \f$ (2r-1) \f$ points
 *          that is subdivided into \f$ r \f$ smaller stencils having \f$ r \f$ points each.
 *          WENO combines all substencils with an estimate of the smoothness of the solution (also obtained from the
 *          substencils) in order to achieve a non-oscillatory high order reconstruction of the face values given
 *          volume averages (cell boundary fluxes in finite volume schemes).
 *          For details see \cite Liu1994 and \cite Jiang1996.
 */
class HighResolutionKoren
{
public:

	/**
	 * @brief Creates the HighResolutionKoren scheme
	 */
	HighResolutionKoren() { }

	/**
	 * @brief Returns the maximum order \f$ r \f$ of the implemented schemes
	 * @return Maximum HR Koren order \f$ r \f$
	 */
	CADET_CONSTEXPR static inline unsigned int maxOrder() CADET_NOEXCEPT { return 2; }

	/**
	 * @brief Returns the maximum stencil size for the implemented schemes
	 * @return Maximum stencil size
	 */
	CADET_CONSTEXPR static inline unsigned int maxStencilSize() CADET_NOEXCEPT { return 2 * maxOrder() - 1; }

	/**
	 * @brief Reconstructs a cell face value from volume averages
	 * @param [in] cellIdx Index of the current cell
	 * @param [in] numCells Number of cells
	 * @param [in] w Stencil that contains the \f$ 2r-1 \f$ volume averages from which the cell face values are reconstructed centered at the 
	 *               current cell (i.e., index 0 is the current cell, -2 the next to previous cell, 2 the next but one cell)
	 * @param [out] result Reconstructed cell face value
	 * @param [out] Dvm Gradient of the reconstructed cell face value (array has to be of size \f$ 2r-1\f$ where \f$ r \f$ is the WENO order)
	 * @tparam StateType Type of the state variables
	 * @tparam StencilType Type of the stencil (can be a dedicated class with overloaded operator[] or a simple pointer)
	 * @return Order of the scheme that was used in the computation
	 */
	template <typename StateType, typename StencilType>
	int reconstruct(unsigned int cellIdx, unsigned int numCells, const StencilType& w, StateType& result, double* const Dvm)
	{
		return reconstruct<StateType, StencilType, true>(cellIdx, numCells, w, result, Dvm);
	}

	/**
	 * @brief Reconstructs a cell face value from volume averages
	 * @param [in] cellIdx Index of the current cell
	 * @param [in] numCells Number of cells
	 * @param [in] w Stencil that contains the \f$ 2r-1 \f$ volume averages from which the cell face values are reconstructed centered at the 
	 *               current cell (i.e., index 0 is the current cell, -2 the next to previous cell, 2 the next but one cell)
	 * @param [out] result Reconstructed cell face value
	 * @tparam StateType Type of the state variables
	 * @tparam StencilType Type of the stencil (can be a dedicated class with overloaded operator[] or a simple pointer)
	 * @return Order of the scheme that was used in the computation
	 */
	template <typename StateType, typename StencilType>
	int reconstruct(unsigned int cellIdx, unsigned int numCells, const StencilType& w, StateType& result)
	{
		return reconstruct<StateType, StencilType, false>(cellIdx, numCells, w, result, nullptr);
	}

	/**
	 * @brief Sets the order
	 * @param [in] order Order of the method
	 */
	inline void order(int order)
	{
		cadet_assert(order <= static_cast<int>(maxOrder()));
		cadet_assert(order > 0);
		_order = order;
	}

	/**
	 * @brief Returns the WENO order
	 * @return Order of the WENO method
	 */
	inline int order() const CADET_NOEXCEPT { return _order; }

	/**
	 * @brief Returns the number of upper diagonals required in the Jacobian
	 * @return Number of required Jacobian upper diagonals
	 */
	inline unsigned int upperBandwidth() const CADET_NOEXCEPT { return _order - 1; }

	/**
	 * @brief Returns the number of lower diagonals required in the Jacobian
	 * @return Number of required Jacobian lower diagonals
	 */
	inline unsigned int lowerBandwidth() const CADET_NOEXCEPT { return _order - 1; }

	/**
	 * @brief Returns the size of the stencil (i.e., the number of required elements)
	 * @return Size of the stencil
	 */
	inline unsigned int stencilSize() const CADET_NOEXCEPT { return 2 * _order - 1; }

private:

	/**
	 * @brief Reconstructs a cell face value from volume averages
	 * @param [in] cellIdx Index of the current cell
	 * @param [in] numCells Number of cells
	 * @param [in] w Stencil that contains the \f$ 2r-1 \f$ volume averages from which the cell face values are reconstructed centered at the 
	 *               current cell (i.e., index 0 is the current cell, -2 the next to previous cell, 2 the next but one cell)
	 * @param [out] result Reconstructed cell face value
	 * @param [out] Dvm Gradient of the reconstructed cell face value (array has to be of size \f$ 2r-1\f$ where \f$ r \f$ is the WENO order)
	 * @tparam StateType Type of the state variables
	 * @tparam StencilType Type of the stencil (can be a dedicated class with overloaded operator[] or a simple pointer)
	 * @tparam wantJac Determines if the gradient is computed (@c true) or not (@c false)
	 * @return Order of the scheme that was used in the computation
	 */
	template <typename StateType, typename StencilType, bool wantJac>
	int reconstruct(unsigned int cellIdx, unsigned int numCells, const StencilType& w, StateType& result, double* const Dvm)
	{
#if defined(ACTIVE_SETFAD) || defined(ACTIVE_SFAD)
		using cadet::sqr;
		using sfad::sqr;
#endif

		// Local order of the scheme that is actually used (may be changed by treatment of boundaries)
		int order = _order;

		// Lower order such that maximum order is used at all points
		// This very statement selects the max. order for the current column cell
		// order = min(maxOrderleft, maxOrderright)
		order = std::min(std::min(static_cast<int>(cellIdx) + 1, _order), std::min(static_cast<int>(numCells - cellIdx), _order));

		// Total stencil size
		const int sl = 2 * order - 1;

		// Simple upwind scheme
		if (order == 1)
		{
			result = w[0];
			if (wantJac)
				*Dvm = 1.0;
			return order;
		}

		return order;
	}

	int _order; //!< Selected order
};

} // namespace cadet

#endif  // LIBCADET_HIGHRESKOREN_HPP_
