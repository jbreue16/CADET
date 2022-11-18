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
#include "Memory.hpp"
#include "AutoDiff.hpp"

#include "io/hdf5/HDF5Writer.hpp"

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
	RadialFlowModel(int nComp, int nCol) : _nComp(nComp), _nCol(nCol), _stencilMemory(sizeof(cadet::active) * 5)
	{
		const int nPureDof = _nCol * _nComp;
		_jacDisc.resize(nPureDof, 2 * nComp, 2 * nComp);
		_jac.resize(nPureDof, 2 * nComp, 2 * nComp);

		_radDispersion = std::vector<cadet::active>(_nComp, 1e-7);

		const double colLen = 0.1;

		equidistantCells(0.1, 0.5, _nCol);

		_params.u = 1.0 * fromVolumetricFlowRate(8e-2, colLen);
		_params.d_rad = _radDispersion.data();
		_params.cellBounds = _cellBounds.data();
		_params.cellCenters = _cellCenters.data();
		_params.cellSizes = _cellCenters.data();
		_params.stencilMemory = &_stencilMemory;
		_params.offsetToBulk = _nComp;
		_params.nCol = _nCol;
		_params.nComp = _nComp;
		_params.offsetToInlet = 0;
		_params.strideCell = _nComp;
	}
	virtual ~RadialFlowModel() CADET_NOEXCEPT { }

	int numPureDofs() const CADET_NOEXCEPT { return _nComp * _nCol; }
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

		// Inlet block: c_i - val = 0
		for (int i = 0; i < _nComp; ++i)
			res[i] = vecStateY[i] - inlet(time, secIdx, i);

		return cadet::model::parts::convdisp::residualKernelRadial<double, double, double, cadet::linalg::BandMatrix::RowIterator, true>(
			cadet::SimulationTime{time, static_cast<unsigned int>(secIdx)},
			vecStateY, vecStateYdot, res, _jac.row(0), _params
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
		cadet::linalg::FactorizableBandMatrix::RowIterator jac = _jacDisc.row(0);
		for (int i = 0; i < _nCol; ++i)
		{
			for (int j = 0; j < _nComp; ++j, ++jac)
			{
				// Add time derivative to main diagonal
				jac[0] += alpha;
			}
		}

		if (!_jacDisc.factorize())
			return 1;

		// Inlet-Bulk coupling Jacobian
		// A * inlet + J * x = rhs
		// J * x = rhs - A * inlet
		const int idxInletCell = (_params.u >= 0.0) ? 0 : _nCol - 1;
		const double factor = static_cast<double>(_params.u) / (static_cast<double>(_params.cellCenters[idxInletCell]) * static_cast<double>(_params.cellSizes[idxInletCell]));
		double* const rhsBulkInlet = rhs + _nComp * (idxInletCell + 1);
		for (int i = 0; i < _nComp; ++i)
		{
			rhsBulkInlet[i] -= -factor * rhs[i];
		}

		if (!_jacDisc.solve(rhs + _nComp))
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
		const int nTimeSaved = _solTimes.size();
		_solTimes.push_back(t);

		const int nDof = numPureDofs();
		_solution.resize(_solution.size() + nDof);
		std::copy(vecStateY + _nComp, vecStateY + nDof + _nComp, _solution.data() + nTimeSaved * nDof);

		_solutionInlet.resize(_solutionInlet.size() + _nComp);
		std::copy(vecStateY, vecStateY + _nComp, _solutionInlet.data() + nTimeSaved * _nComp);

		_solutionOutlet.resize(_solutionOutlet.size() + _nComp);
		if (_params.u > 0.0)
		{
			// Flow from inner to outer
			std::copy(vecStateY + nDof, vecStateY + nDof + _nComp, _solutionOutlet.data() + nTimeSaved * _nComp);
		}
		else
		{
			// Flow from outer to inner
			std::copy(vecStateY + _nComp, vecStateY + 2 * _nComp, _solutionOutlet.data() + nTimeSaved * _nComp);
		}
	}

	const std::vector<double>& solutionTimes() const CADET_NOEXCEPT { return _solTimes; }
	const std::vector<double>& solution() const CADET_NOEXCEPT { return _solution; }
	const std::vector<double>& solutionInlet() const CADET_NOEXCEPT { return _solutionInlet; }
	const std::vector<double>& solutionOutlet() const CADET_NOEXCEPT { return _solutionOutlet; }
	int numComp() const CADET_NOEXCEPT { return _nComp; }
	int numCol() const CADET_NOEXCEPT { return _nCol; }

	std::vector<double> coordinates() const CADET_NOEXCEPT
	{
		std::vector<double> coords(_cellCenters.size(), 0.0);
		for (int i = 0; i < _cellCenters.size(); ++i)
			coords[i] = static_cast<double>(_cellCenters[i]);
		
		return coords;
	}

protected:
	int _nComp;
	int _nCol;

	cadet::model::parts::convdisp::RadialFlowParameters<double> _params;
	cadet::linalg::BandMatrix _jac;
	cadet::linalg::FactorizableBandMatrix _jacDisc;

	std::vector<cadet::active> _radDispersion;
	std::vector<cadet::active> _cellCenters;
	std::vector<cadet::active> _cellSizes;
	std::vector<cadet::active> _cellBounds;
	cadet::ArrayPool _stencilMemory;

	std::vector<double> _solTimes;
	std::vector<double> _solution;
	std::vector<double> _solutionInlet;
	std::vector<double> _solutionOutlet;

	double inlet(double t, int secIdx, int comp) const CADET_NOEXCEPT
	{
		if (t <= 10.0)
			return 0.1 * t;
		if (t <= 50.0)
			return 1.0;
		if (t <= 60.0)
			return 1.0 - (t-50.0) * 0.1;
		return 0.0;
//		return 1.0;
	}

	void equidistantCells(double inner, double outer, int nCol)
	{
		const double dr = (outer - inner) / nCol;
		std::vector<cadet::active> centers(nCol, 0.0);
		_cellSizes = std::vector<cadet::active>(nCol, dr);
		std::vector<cadet::active> bounds(nCol + 1, 0.0);

		for (int i = 0; i < nCol; ++i)
		{
			centers[i] = (i + 0.5) * dr;
			bounds[i] = i * dr;
		}
		bounds[nCol] = outer;

		_cellCenters = std::move(centers);
		_cellBounds = std::move(bounds);
	}

	double fromVolumetricFlowRate(double volRate, double len)
	{
		const double pi = 3.14159265358979323846;
		return volRate / (pi * 2.0 * len);
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

	const double tEnd = 160.0;
	std::vector<double> secTimes = {0.0, tEnd};
	std::vector<double> solTimes(161, 0.0);

	for (int i = 0; i <= tEnd; ++i)
		solTimes[i] = i;

	RadialFlowModel model(1, 200);

	cadet::test::TimeIntegrator sim;	
	sim.configureTimeIntegrator(1e-6, 1e-8, 1e-4, 100000, 0.0);
	sim.setSectionTimes(secTimes);
	sim.setSolutionTimes(solTimes);
	sim.initializeModel(model);

	sim.integrate();

	cadet::io::HDF5Writer writer;
	writer.openFile("radial.h5", "co");
	writer.vector("SOLUTION_TIMES", model.solutionTimes());
	const std::vector<std::size_t> dims = {model.solutionTimes().size(), static_cast<std::size_t>(model.numCol()), static_cast<std::size_t>(model.numComp())};
	writer.template tensor<double>("SOLUTION", 3, dims.data(), model.solution());
	writer.template matrix<double>("SOLUTION_INLET", model.solutionTimes().size(), model.numComp(), model.solutionInlet());
	writer.template matrix<double>("SOLUTION_OUTLET", model.solutionTimes().size(), model.numComp(), model.solutionOutlet());
	writer.template vector<double>("COORDS", model.coordinates());
	writer.closeFile();
	return 0;
}
