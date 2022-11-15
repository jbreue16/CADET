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

#include "TimeIntegrator.hpp"
#include "cadet/Logging.hpp"
#include "Logging.hpp"

#include "model/parts/RadialConvectionDispersionKernel.hpp"
#include "linalg/BandMatrix.hpp"

#include <algorithm>
#include <iostream>
#include <vector>

// Uncomment the next line to enable logging output of CADET in unit tests
#define CADETTEST_ENABLE_LOG

#ifdef CADETTEST_ENABLE_LOG
	#include "cadet/Logging.hpp"
	#include <iostream>

	class LogReceiver : public cadet::ILogReceiver
	{
	public:
		LogReceiver() { }

		virtual void message(const char* file, const char* func, const unsigned int line, cadet::LogLevel lvl, const char* lvlStr, const char* message)
		{
			std::cout << '[' << lvlStr << ": " << func << "::" << line << "] " << message << std::flush;
		}
	};
#endif

class RadialFlowModel : public cadet::test::IDiffEqModel
{
public:
	RadialFlowModel(int nComp, int nCol) : _nComp(nComp), _nCol(nCol)
	{
		const int nDof = (_nCol + 1) * _nComp;
		_jacDisc.resize(nDof, 2 * nComp, 2 * nComp);
		_jac.resize(nDof, 2 * nComp, 2 * nComp);

		_params.offsetToBulk = _nComp;
		_params.nCol = _nCol;
		_params.nComp = _nComp;
		_params.offsetToInlet = 0;
		_params.strideCell = _nComp;
	}
	virtual ~RadialFlowModel() CADET_NOEXCEPT { }

	virtual int numDofs() const CADET_NOEXCEPT { return _nComp * (_nCol + 1); }

	virtual void notifyDiscontinuousSectionTransition(double t, int secIdx, double* vecStateY, double* vecStateYdot)
	{
		if (t == 0.0)
		{
			// Consistent init
			std::fill_n(vecStateY, numDofs(), 0.0);
			std::fill_n(vecStateYdot, numDofs(), 0.0);
		}
	}

	virtual int residual(double time, int secIdx, double const* vecStateY, double const* vecStateYdot, double* res)
	{
		return residualWithJacobian(time, secIdx, vecStateY, vecStateYdot, res);
	}

	virtual int residualWithJacobian(double time, int secIdx, double const* vecStateY, double const* vecStateYdot, double* res)
	{
		_jac.setAll(0.0);

		// Inlet block: val - c_i = 0
		for (int i = 0; i < _nComp; ++i)
		{
			_jac.centered(i, 0) = -1.0;
			res[i] = inlet(time, secIdx, i) - vecStateY[i];
		}

		// Inlet-Bulk coupling Jacobian
		const int idxInletCell = (_params.u >= 0.0) ? 0 : _nCol - 1;
		const double factor = _params.u / (_params.cellCenters[idxInletCell] * _params.cellSizes[idxInletCell]);
		for (int i = 0; i < _nComp; ++i)
			_jac.centered(i + _nComp, -_nComp) = -factor;

		return cadet::model::parts::convdisp::residualKernelRadial<double, double, double, cadet::linalg::BandedRowIterator, true>(
			cadet::SimulationTime{time, secIdx},
			vecStateY, vecStateYdot, res, _jac.row(), _params
		);
	}

	virtual double residualNorm(double time, int secIdx, double const* vecStateY, double const* vecStateYdot)
	{
		return 0.0;
	}

	virtual int linearSolve(double t, double alpha, double tol, double* rhs, double const* weight,
		double const* vecStateY, double const* vecStateYdot)
	{
		_jacDisc.copyOver(_jac);

		// Add time derivative
		cadet::linalg::FactorizableBandMatrix::RowIterator jac = _jacDisc.row(_nComp);
		for (unsigned int i = 0; i < _nCol; ++i)
		{
			for (unsigned int j = 0; j < _nComp; ++j, ++jac)
			{
				// Add time derivative to main diagonal
				jac[0] += alpha;
			}
		}

		if (!_jacDisc.factorize())
			return 1;
		if (!_jacDisc.solve(rhs))
			return 1;
		return 0;
	}

	virtual void applyInitialCondition(double* vecStateY, double* vecStateYdot) const
	{
		// Consistent init
		std::fill_n(vecStateY, numDofs(), 0.0);
		std::fill_n(vecStateYdot, numDofs(), 0.0);
	}

	virtual void saveSolution(double t, double const* vecStateY, double const* vecStateYdot)
	{
		_solTimes.push_back(t);

		const int nDof = numDofs();
		_solution.resize(_solution.size() + nDof);
		std::copy(vecStateY, vecStateY + nDof, _solution.data() + _solPos);
		_solPos += nDof;
	}

protected:
	int _nComp;
	int _nCol;

	cadet::model::parts::convdisp::RadialFlowParameters<double> _params;
	cadet::linalg::BandMatrix _jac;
	cadet::linalg::FactorizableBandMatrix _jacDisc;

	std::vector<double> _solTimes;
	std::vector<double> _solution;
	int _solPos;

	double inlet(double t, int secIdx, int comp) const CADET_NOEXCEPT
	{
		return 0.0;
	}
};


int main(int argc, char* argv[])
{
#ifdef CADETTEST_ENABLE_LOG
	// Set LogLevel in CADET library
	const cadet::LogLevel logLevel = cadet::LogLevel::Trace;
	LogReceiver lr;
	cadet::setLogReceiver(&lr);
	cadet::setLogLevel(logLevel);
#endif

	const double tEnd = 100.0;
	std::vector<double> secTimes = {0.0, tEnd};
	std::vector<double> solTimes(101, 0.0);

	for (int i = 0; i <= tEnd; ++i)
		solTimes[i] = i;

	RadialFlowModel model;

	cadet::test::TimeIntegrator sim;	
	sim.configureTimeIntegrator(1e-6, 1e-8, 1e-4, 100000, 0.0);
	sim.setSectionTimes(secTimes);
	sim.setSolutionTimes(solTimes);
	sim.initializeModel(model);

	sim.integrate();

	return 0;
}
