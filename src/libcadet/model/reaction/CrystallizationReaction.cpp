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

#include "model/reaction/ReactionModelBase.hpp"
#include "model/ModelUtils.hpp"
#include "cadet/Exceptions.hpp"
#include "model/Parameters.hpp"
#include "LocalVector.hpp"
#include "SimulationTypes.hpp"
#include "linalg/ActiveDenseMatrix.hpp"
#include "Memory.hpp"

#include <functional>
#include <algorithm>
#include <unordered_map>
#include <string>
#include <iterator>
#include <vector>

namespace cadet
{

	namespace model
	{

		/**
		 * @brief Defines the crystallization reaction model
		 */
		class CrystallizationReaction : public IDynamicReactionModel
		{
		public:

			CrystallizationReaction() : _nComp(0), _nBins(0), _bins(0), _binCenters(0), _binSizes(0) { }
			virtual ~CrystallizationReaction() CADET_NOEXCEPT { }

			static const char* identifier() { return "CRYSTALLIZATION"; }
			virtual const char* name() const CADET_NOEXCEPT { return identifier(); }

			virtual bool requiresConfiguration() const CADET_NOEXCEPT { return true; }
			virtual bool usesParamProviderInDiscretizationConfig() const CADET_NOEXCEPT { return false; }

			virtual bool configure(IParameterProvider& paramProvider, UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx)
			{
				readScalarParameterOrArray(_bins, paramProvider, "CRY_BINS", 1);

				if (_bins.size() != _nBins + 1)
					throw InvalidParameterException("Expected CRY_BINS to have " + std::to_string(_nBins + 1) + " elements (got " + std::to_string(_bins.size()) + ")");

				registerParam1DArray(_parameters, _bins, [=](bool multi, unsigned int idx) { return makeParamId(hashString("CRY_BINS"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, idx, SectionIndep); });

				_binCenters.resize(_nBins);
				_binSizes.resize(_nBins);
				_binCenterDists.resize(_nBins);
				updateBinCoords();


				_nucleiMassDensity = paramProvider.getDouble("CRY_NUCLEI_MASS_DENSITY");
				_parameters[makeParamId(hashString("CRY_NUCLEI_MASS_DENSITY"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_nucleiMassDensity;

				_volShapeFactor = paramProvider.getDouble("CRY_VOL_SHAPE_FACTOR");
				_parameters[makeParamId(hashString("CRY_VOL_SHAPE_FACTOR"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_volShapeFactor;

				_primaryNucleationRate = paramProvider.getDouble("CRY_PRIMARY_NUCLEATION_RATE");
				_parameters[makeParamId(hashString("CRY_PRIMARY_NUCLEATION_RATE"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_primaryNucleationRate;

				_secondaryNucleationRate = paramProvider.getDouble("CRY_SECONDARY_NUCLEATION_RATE");
				_parameters[makeParamId(hashString("CRY_SECONDARY_NUCLEATION_RATE"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_secondaryNucleationRate;

				_growthRateConstant = paramProvider.getDouble("CRY_GROWTH_RATE_CONSTANT");
				_parameters[makeParamId(hashString("CRY_GROWTH_RATE_CONSTANT"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_growthRateConstant;

				_growthConstant = paramProvider.getDouble("CRY_GROWTH_CONSTANT");
				_parameters[makeParamId(hashString("CRY_GROWTH_CONSTANT"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_growthConstant;

				_growthDispersionRate = paramProvider.getDouble("CRY_GROWTH_DISPERSION_RATE");
				_parameters[makeParamId(hashString("CRY_GROWTH_DISPERSION_RATE"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_growthDispersionRate;

				_growthSchemeOrder = paramProvider.getInt("CRY_GROWTH_SCHEME_ORDER");

				_a = paramProvider.getDouble("CRY_A");
				_parameters[makeParamId(hashString("CRY_A"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_a;

				_b = paramProvider.getDouble("CRY_B");
				_parameters[makeParamId(hashString("CRY_B"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_b;

				_g = paramProvider.getDouble("CRY_G");
				_parameters[makeParamId(hashString("CRY_G"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_g;

				_p = paramProvider.getDouble("CRY_P");
				_parameters[makeParamId(hashString("CRY_P"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_p;

				_k = paramProvider.getDouble("CRY_K");
				_parameters[makeParamId(hashString("CRY_K"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_k;

				_u = paramProvider.getDouble("CRY_U");
				_parameters[makeParamId(hashString("CRY_U"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_u;

				return true;
			}

			virtual bool configureModelDiscretization(IParameterProvider& paramProvider, unsigned int nComp, unsigned int const* nBound, unsigned int const* boundOffset)
			{
				_nComp = nComp;

				// Comp 0 is substrate, last comp is equilibrium
				_nBins = _nComp - 2;

				if (_nBins < 1)
					throw InvalidParameterException("Expected at least 3 components (got " + std::to_string(_nComp) + ")");

				return true;
			}

			std::unordered_map<ParameterId, double> getAllParameterValues() const
			{
				std::unordered_map<ParameterId, double> data;
				std::transform(_parameters.begin(), _parameters.end(), std::inserter(data, data.end()),
					[](const std::pair<const ParameterId, active*>& p) { return std::make_pair(p.first, static_cast<double>(*p.second)); });
				return data;
			}

			bool hasParameter(const ParameterId& pId) const
			{
				return _parameters.find(pId) != _parameters.end();
			}

			bool setParameter(const ParameterId& pId, int value) { return false; }
			bool setParameter(const ParameterId& pId, bool value) { return false; }

			bool setParameter(const ParameterId& pId, double value)
			{
				auto paramHandle = _parameters.find(pId);
				if (paramHandle != _parameters.end())
				{
					paramHandle->second->setValue(value);

					// TODO: This does not handle parameter sensitivities wrt. to bin size
					if (pId.name == hashString("CRY_BINS"))
						updateBinCoords();

					return true;
				}

				return false;
			}

			active* getParameter(const ParameterId& pId)
			{
				auto paramHandle = _parameters.find(pId);
				if (paramHandle != _parameters.end())
				{
					return paramHandle->second;
				}

				return nullptr;
			}

			virtual void setExternalFunctions(IExternalFunction** extFuns, unsigned int size) { }
			virtual bool dependsOnTime() const CADET_NOEXCEPT { return false; }
			virtual bool requiresWorkspace() const CADET_NOEXCEPT { return false; }
			virtual unsigned int workspaceSize(unsigned int nComp, unsigned int totalNumBoundStates, unsigned int const* nBoundStates) const CADET_NOEXCEPT
			{
				return 0;
			}

			virtual unsigned int numReactionsLiquid() const CADET_NOEXCEPT { return 1; }
			virtual unsigned int numReactionsCombined() const CADET_NOEXCEPT { return 1; }

			CADET_DYNAMICREACTIONMODEL_BOILERPLATE

		protected:

			std::unordered_map<ParameterId, active*> _parameters; //!< Map used to translate ParameterIds to actual variables
			int _nComp; //!< Number of components
			int _nBins; //!< Number of crystal size bins

			std::vector<active> _bins;
			std::vector<active> _binCenters;
			std::vector<active> _binSizes;
			std::vector<active> _binCenterDists;
			active _nucleiMassDensity; //!< rho
			active _volShapeFactor; //!< k_v
			active _primaryNucleationRate; //!< k_p
			active _secondaryNucleationRate; //!< k_b
			active _growthRateConstant; //!< k_g
			active _growthConstant; //!< gamma
			active _growthDispersionRate; //!< D_g
			active _growthSchemeOrder; // either 1 or 2
			active _a; //!< System constant
			active _b; //!< System constant
			active _g; //!< System constant
			active _p; //!< System constant
			active _k; //!< System constant
			active _u; //!< System constant

			void updateBinCoords() CADET_NOEXCEPT
			{
				for (int i = 0; i < _nBins; ++i)
				{
					if (cadet_likely(i + 1 < _nBins))
						_binCenterDists[i] = 0.5 * (_bins[i + 2] - _bins[i]);

					_binCenters[i] = 0.5 * (_bins[i] + _bins[i + 1]);
					_binSizes[i] = _bins[i + 1] - _bins[i];
				}
			}

			template <typename StateType, typename ResidualType, typename ParamType, typename FactorType>
			int residualLiquidImpl(double t, unsigned int secIdx, const ColumnPosition& colPos,
				StateType const* y, ResidualType* res, const FactorType& factor, LinearBufferAllocator workSpace) const
			{
				typedef typename DoubleActivePromoter<StateType, ParamType>::type StateParam;
				// c_0 is component 0
				// c_eq is last component
				// x_c = bins[0] (i.e. left boundary of first bin is critical nuclei size x_c)

				// Pointer to crystal bins
				StateType const* const yCrystal = y + 1;
				ResidualType* const resCrystal = res + 1;

				// s = (c_0 - c_eq) / c_eq = c_0 / c_eq - 1
				const StateType s = y[0] / y[_nComp - 1] - 1.0;
				const ParamType massDensityShapeFactor = static_cast<ParamType>(_nucleiMassDensity) * static_cast<ParamType>(_volShapeFactor);

				// Numerical approximation of integrals via midpoint rule as volume averages are
				// second order accurate approximation of value at midpoint / center of mass of volume

				const StateParam k_g_times_s_g = static_cast<ParamType>(_growthRateConstant) * pow(s, static_cast<ParamType>(_g));

				StateParam M = 0.0;
				StateParam substrateConversion = 0.0;
				for (int i = 0; i < _nBins; ++i)
				{
					M += yCrystal[i] * static_cast<ParamType>(_binCenters[i]) * static_cast<ParamType>(_binCenters[i]) * static_cast<ParamType>(_binCenters[i]) * static_cast<ParamType>(_binSizes[i]);
					substrateConversion += yCrystal[i] * static_cast<ParamType>(_binCenters[i]) * static_cast<ParamType>(_binCenters[i]) * (static_cast<ParamType>(_a) + static_cast<ParamType>(_growthConstant) * pow(static_cast<ParamType>(_binCenters[i]), static_cast<ParamType>(_p))) * static_cast<ParamType>(_binSizes[i]);
				}
				M *= massDensityShapeFactor;
				substrateConversion *= 3.0 * k_g_times_s_g;

				// B_0 = primary + secondary nucleation rate
				const StateParam B_0 = static_cast<ParamType>(_primaryNucleationRate) * pow(s, static_cast<ParamType>(_u)) + static_cast<ParamType>(_secondaryNucleationRate) * pow(s, static_cast<ParamType>(_b)) * pow(M, static_cast<ParamType>(_k));
				const ParamType x_c_3 = static_cast<ParamType>(_binCenters[0]) * static_cast<ParamType>(_binCenters[0]) * static_cast<ParamType>(_binCenters[0]);

				// adjust c_feed for each time step
				// res[0] -= factor * massDensityShapeFactor * (B_0 * x_c_3 + substrateConversion);

				StateParam v_g = 0.0;

				// upwind scheme
				if (static_cast<ParamType>(_growthSchemeOrder) == 1)
				{
					for (int i = 0; i < _nBins; ++i)
					{
						// Flux through left face
						if (cadet_likely((i > 0) && (i + 1 < _nBins)))
						{
							// flux through the left face
							resCrystal[i] += factor / static_cast<ParamType>(_binSizes[i]) * (v_g * yCrystal[i - 1]) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i] - yCrystal[i - 1]) / static_cast<ParamType>(_binCenterDists[i - 1]);
							// flux through the right face
							v_g = k_g_times_s_g * (static_cast<ParamType>(_a) + static_cast<ParamType>(_growthConstant) * pow(static_cast<ParamType>(_bins[i + 1]), static_cast<ParamType>(_p)));
							resCrystal[i] -= factor / static_cast<ParamType>(_binSizes[i]) * (v_g * yCrystal[i]) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i + 1] - yCrystal[i]) / static_cast<ParamType>(_binCenterDists[i]);
						}
						else if (i == 0)
						{
							// Left boundary condition
							resCrystal[i] += factor / static_cast<ParamType>(_binSizes[i]) * B_0;
							// first order approximation
							v_g = k_g_times_s_g * (static_cast<ParamType>(_a) + static_cast<ParamType>(_growthConstant) * pow(static_cast<ParamType>(_bins[i + 1]), static_cast<ParamType>(_p)));
							resCrystal[i] -= factor / static_cast<ParamType>(_binSizes[i]) * v_g * (0.5 * yCrystal[i] + 0.5 * yCrystal[i + 1]) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i + 1] - yCrystal[i]) / static_cast<ParamType>(_binCenterDists[i]);
						}
						else
						{
							// first order approximation
							resCrystal[i] += factor / static_cast<ParamType>(_binSizes[i]) * v_g * (0.5 * yCrystal[i - 1] + 0.5 * yCrystal[i]) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i] - yCrystal[i - 1]) / static_cast<ParamType>(_binCenterDists[i - 1]);
							// no flux
						}
					}
				}
				// high-resolutuion scheme
				else if (static_cast<ParamType>(_growthSchemeOrder) == 2)
				{
					StateParam r_x_i = 0.0;
					StateParam phi = 0.0;
					StateParam R_i = 0.0;
					for (int i = 0; i < _nBins; ++i)
					{
						if (cadet_likely((i > 1) && (i + 1 < _nBins)))
						{
							// Flux through left face, modified van Leer flux limiter
							r_x_i = (yCrystal[i - 1] - yCrystal[i - 2] + 1e-10) * (static_cast<ParamType>(_binSizes[i]) + static_cast<ParamType>(_binSizes[i - 1])) / (yCrystal[i] - yCrystal[i - 1] + 1e-10) / (static_cast<ParamType>(_binSizes[i - 1]) + static_cast<ParamType>(_binSizes[i - 2]));
							R_i = (static_cast<ParamType>(_binSizes[i]) + static_cast<ParamType>(_binSizes[i - 1])) / static_cast<ParamType>(_binSizes[i - 1]);
							if (cadet_likely(r_x_i > 0))
							{
								phi = R_i * r_x_i / (R_i - 1.0 + r_x_i);
							}
							else
							{
								phi = 0.0;
							}
							// Use v_g from previous iteration, this results in the same flux => conservative; apply first order approximation to the diffusion, 2nd order accurate F_{i-1}
							resCrystal[i] += factor / static_cast<ParamType>(_binSizes[i]) * v_g * (yCrystal[i - 1] + phi / R_i * (yCrystal[i] - yCrystal[i - 1])) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i] - yCrystal[i - 1]) / static_cast<ParamType>(_binCenterDists[i - 1]);
							// Flux through left face, modified van Leer flux limiter, update r_x_i, R_i, and growth rate
							r_x_i = (yCrystal[i] - yCrystal[i - 1] + 1e-10) * (static_cast<ParamType>(_binSizes[i + 1]) + static_cast<ParamType>(_binSizes[i])) / (yCrystal[i + 1] - yCrystal[i] + 1e-10) / (static_cast<ParamType>(_binSizes[i]) + static_cast<ParamType>(_binSizes[i - 1]));
							R_i = (static_cast<ParamType>(_binSizes[i + 1]) + static_cast<ParamType>(_binSizes[i])) / static_cast<ParamType>(_binSizes[i]);
							if (cadet_likely(r_x_i > 0))
							{
								phi = R_i * r_x_i / (R_i - 1.0 + r_x_i);
							}
							else
							{
								phi = 0.0;
							}
							v_g = k_g_times_s_g * (static_cast<ParamType>(_a) + static_cast<ParamType>(_growthConstant) * pow(static_cast<ParamType>(_bins[i + 1]), static_cast<ParamType>(_p)));
							resCrystal[i] -= factor / static_cast<ParamType>(_binSizes[i]) * v_g * (yCrystal[i] + phi / R_i * (yCrystal[i + 1] - yCrystal[i])) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i + 1] - yCrystal[i]) / static_cast<ParamType>(_binCenterDists[i]);
						}
						else if (i == 1)
						{
							// apply first order approximation to F_1 and the diffusion, ideally HR scheme should be applied to F_{1+1/2}, but this saves some code and logical if statement
							resCrystal[i] += factor / static_cast<ParamType>(_binSizes[i]) * v_g * (0.5 * yCrystal[i - 1] + 0.5 * yCrystal[i]) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i] - yCrystal[i - 1]) / static_cast<ParamType>(_binCenterDists[i - 1]);
							v_g = k_g_times_s_g * (static_cast<ParamType>(_a) + static_cast<ParamType>(_growthConstant) * pow(static_cast<ParamType>(_bins[i + 1]), static_cast<ParamType>(_p)));
							resCrystal[i] -= factor / static_cast<ParamType>(_binSizes[i]) * v_g * (0.5 * yCrystal[i] + 0.5 * yCrystal[i + 1]) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i + 1] - yCrystal[i]) / static_cast<ParamType>(_binCenterDists[i]);
						}
						else if (i + 1 == _nBins)
						{
							// first order approximation
							resCrystal[i] += factor / static_cast<ParamType>(_binSizes[i]) * v_g * (0.5 * yCrystal[i - 1] + 0.5 * yCrystal[i]) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i] - yCrystal[i - 1]) / static_cast<ParamType>(_binCenterDists[i - 1]);
							// no flux
						}
						else
						{
							// left boundary condition
							resCrystal[i] += factor / static_cast<ParamType>(_binSizes[i]) * B_0;
							// first order approximation to F_1 and the diffusion
							v_g = k_g_times_s_g * (static_cast<ParamType>(_a) + static_cast<ParamType>(_growthConstant) * pow(static_cast<ParamType>(_bins[i + 1]), static_cast<ParamType>(_p)));
							resCrystal[i] -= factor / static_cast<ParamType>(_binSizes[i]) * v_g * (0.5 * yCrystal[i] + 0.5 * yCrystal[i + 1]) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i + 1] - yCrystal[i]) / static_cast<ParamType>(_binCenterDists[i]);
						}
					}
				}
				// weno 3 on random grids
				else if (static_cast<ParamType>(_growthSchemeOrder) == 3)
				{
					// declare the vectors to store the coefficients
					std::vector<active> q_1_right_coeff(_nBins - 1);
					std::vector<active> q_0_right_coeff(_nBins - 1);
					std::vector<active> C_right_coeff(_nBins - 1);
					std::vector<active> IS_0_coeff(_nBins - 1);
					std::vector<active> IS_1_coeff(_nBins - 1);
					ParamType delta_sum = 0.0;
					ParamType delta_right_sum = 0.0;
					ParamType delta_left_sum = 0.0;
					// calculate the coefficients and store them, the first entry is 0
					for (int i = 1; i + 1 < _nBins; ++i)
					{
						delta_sum = static_cast<ParamType>(_binSizes[i]) + static_cast<ParamType>(_binSizes[i - 1]) + static_cast<ParamType>(_binSizes[i + 1]);
						delta_left_sum = static_cast<ParamType>(_binSizes[i]) + static_cast<ParamType>(_binSizes[i - 1]);
						delta_right_sum = static_cast<ParamType>(_binSizes[i]) + static_cast<ParamType>(_binSizes[i + 1]);
						q_0_right_coeff[i] = static_cast<ParamType>(_binSizes[i + 1]) / delta_right_sum;
						q_1_right_coeff[i] = 1.0 + static_cast<ParamType>(_binSizes[i]) / delta_left_sum;
						C_right_coeff[i] = delta_left_sum / delta_sum;
						IS_0_coeff[i] = (2.0 * static_cast<ParamType>(_binSizes[i]) / delta_right_sum) * (2.0 * static_cast<ParamType>(_binSizes[i]) / delta_right_sum);
						IS_1_coeff[i] = (2.0 * static_cast<ParamType>(_binSizes[i]) / delta_left_sum) * (2.0 * static_cast<ParamType>(_binSizes[i]) / delta_left_sum);
					}
					// ode input
					StateParam IS_0 = 0.0;
					StateParam IS_1 = 0.0;
					StateParam alpha_0 = 0.0;
					StateParam alpha_1 = 0.0;
					StateParam W_0 = 0.0;
					StateParam W_1 = 0.0;
					StateParam q_0 = 0.0;
					StateParam q_1 = 0.0;
					for (int i = 0; i < _nBins; ++i)
					{
						if (cadet_likely((i > 1) && (i + 1 < _nBins)))
						{
							// flux through left face
							IS_0 = static_cast<ParamType>(IS_0_coeff[i - 1]) * (yCrystal[i] - yCrystal[i - 1]) * (yCrystal[i] - yCrystal[i - 1]);
							IS_1 = static_cast<ParamType>(IS_1_coeff[i - 1]) * (yCrystal[i - 1] - yCrystal[i - 2]) * (yCrystal[i - 1] - yCrystal[i - 2]);
							alpha_0 = static_cast<ParamType>(C_right_coeff[i - 1]) / (1e-10 + IS_0) / (1e-10 + IS_0);
							alpha_1 = (1.0 - static_cast<ParamType>(C_right_coeff[i - 1])) / (1e-10 + IS_1) / (1e-10 + IS_1);
							W_0 = alpha_0 / (alpha_0 + alpha_1);
							W_1 = 1.0 - W_0;
							q_0 = static_cast<ParamType>(q_0_right_coeff[i]) * yCrystal[i - 1] + (1.0 - static_cast<ParamType>(q_0_right_coeff[i])) * yCrystal[i];
							q_1 = static_cast<ParamType>(q_1_right_coeff[i]) * yCrystal[i - 1] + (1.0 - static_cast<ParamType>(q_1_right_coeff[i])) * yCrystal[i - 2];
							resCrystal[i] += factor / static_cast<ParamType>(_binSizes[i]) * v_g * (W_0 * q_0 + W_1 * q_1) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i] - yCrystal[i - 1]) / static_cast<ParamType>(_binCenterDists[i - 1]);
							// flux through right face, update IS, alpha, W, q and growth rate
							IS_0 = static_cast<ParamType>(IS_0_coeff[i]) * (yCrystal[i + 1] - yCrystal[i]) * (yCrystal[i + 1] - yCrystal[i]);
							IS_1 = static_cast<ParamType>(IS_1_coeff[i]) * (yCrystal[i] - yCrystal[i - 1]) * (yCrystal[i] - yCrystal[i - 1]);
							alpha_0 = static_cast<ParamType>(C_right_coeff[i]) / (1e-10 + IS_0) / (1e-10 + IS_0);
							alpha_1 = (1.0 - static_cast<ParamType>(C_right_coeff[i])) / (1e-10 + IS_1) / (1e-10 + IS_1);
							W_0 = alpha_0 / (alpha_0 + alpha_1);
							W_1 = 1.0 - W_0;
							q_0 = static_cast<ParamType>(q_0_right_coeff[i]) * yCrystal[i] + (1.0 - static_cast<ParamType>(q_0_right_coeff[i])) * yCrystal[i + 1];
							q_1 = static_cast<ParamType>(q_1_right_coeff[i]) * yCrystal[i] + (1.0 - static_cast<ParamType>(q_1_right_coeff[i])) * yCrystal[i - 1];
							v_g = k_g_times_s_g * (static_cast<ParamType>(_a) + static_cast<ParamType>(_growthConstant) * pow(static_cast<ParamType>(_bins[i + 1]), static_cast<ParamType>(_p)));
							resCrystal[i] -= factor / static_cast<ParamType>(_binSizes[i]) * v_g * (W_0 * q_0 + W_1 * q_1) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i + 1] - yCrystal[i]) / static_cast<ParamType>(_binCenterDists[i]);
						}
						// boundary condition
						else if (i + 1 == _nBins)
						{
							// first order approximation
							resCrystal[i] += factor / static_cast<ParamType>(_binSizes[i]) * v_g * (0.5 * yCrystal[i - 1] + 0.5 * yCrystal[i]) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i] - yCrystal[i - 1]) / static_cast<ParamType>(_binCenterDists[i - 1]);
							// no flux, regularity boundary condition
						}
						else if (i == 1)
						{
							// first order approximation, ideally weno scheme should be applied to F_{1+1/2}, but this saves some code and logical if statement
							resCrystal[i] += factor / static_cast<ParamType>(_binSizes[i]) * v_g * (0.5 * yCrystal[i - 1] + 0.5 * yCrystal[i]) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i] - yCrystal[i - 1]) / static_cast<ParamType>(_binCenterDists[i - 1]);
							v_g = k_g_times_s_g * (static_cast<ParamType>(_a) + static_cast<ParamType>(_growthConstant) * pow(static_cast<ParamType>(_bins[i + 1]), static_cast<ParamType>(_p)));
							resCrystal[i] -= factor / static_cast<ParamType>(_binSizes[i]) * v_g * (0.5 * yCrystal[i] + 0.5 * yCrystal[i + 1]) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i + 1] - yCrystal[i]) / static_cast<ParamType>(_binCenterDists[i]);
						}
						else
						{
							// nucleation boundary condition
							resCrystal[i] += factor / static_cast<ParamType>(_binSizes[i]) * B_0;
							// first order approximation
							v_g = k_g_times_s_g * (static_cast<ParamType>(_a) + static_cast<ParamType>(_growthConstant) * pow(static_cast<ParamType>(_bins[i + 1]), static_cast<ParamType>(_p)));
							resCrystal[i] -= factor / static_cast<ParamType>(_binSizes[i]) * v_g * (0.5 * yCrystal[i] + 0.5 * yCrystal[i + 1]) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i + 1] - yCrystal[i]) / static_cast<ParamType>(_binCenterDists[i]);
						}
					}
				}
				// weno 5 on random grids
				else if (static_cast<ParamType>(_growthSchemeOrder) == 4)
				{
					// declare the vectors to store the coefficients
					std::vector<active> q_2_coeff_1(_nBins - 2);
					std::vector<active> q_2_coeff_2(_nBins - 2);
					std::vector<active> q_1_coeff_1(_nBins - 2);
					std::vector<active> q_1_coeff_2(_nBins - 2);
					std::vector<active> q_0_coeff_1(_nBins - 2);
					std::vector<active> q_0_coeff_2(_nBins - 2);
					std::vector<active> C_0(_nBins - 2);
					std::vector<active> C_1(_nBins - 2);
					std::vector<active> C_2(_nBins - 2);
					std::vector<active> IS_0_coeff_1(_nBins - 2);
					std::vector<active> IS_0_coeff_2(_nBins - 2);
					std::vector<active> IS_0_coeff_3(_nBins - 2);
					std::vector<active> IS_1_coeff_1(_nBins - 2);
					std::vector<active> IS_1_coeff_2(_nBins - 2);
					std::vector<active> IS_1_coeff_3(_nBins - 2);
					std::vector<active> IS_2_coeff_1(_nBins - 2);
					std::vector<active> IS_2_coeff_2(_nBins - 2);
					std::vector<active> IS_2_coeff_3(_nBins - 2);
					ParamType delta_sum = 0.0;
					ParamType delta_sum_I0 = 0.0;
					ParamType delta_sum_I1 = 0.0;
					ParamType delta_sum_I2 = 0.0;
					ParamType IS_0_pre = 0.0;
					ParamType IS_1_pre = 0.0;
					ParamType IS_2_pre = 0.0;
					// calculate the coefficients and store them, the first and second entry is 0, all positive
					for (int i = 2; i + 2 < _nBins; ++i)
					{
						delta_sum = static_cast<ParamType>(_binSizes[i - 2]) + static_cast<ParamType>(_binSizes[i - 1]) + static_cast<ParamType>(_binSizes[i]) + static_cast<ParamType>(_binSizes[i + 1]) + static_cast<ParamType>(_binSizes[i + 2]);
						delta_sum_I0 = static_cast<ParamType>(_binSizes[i - 2]) + static_cast<ParamType>(_binSizes[i - 1]) + static_cast<ParamType>(_binSizes[i]);
						delta_sum_I1 = static_cast<ParamType>(_binSizes[i - 1]) + static_cast<ParamType>(_binSizes[i]) + static_cast<ParamType>(_binSizes[i + 1]);
						delta_sum_I2 = static_cast<ParamType>(_binSizes[i]) + static_cast<ParamType>(_binSizes[i + 1]) + static_cast<ParamType>(_binSizes[i + 2]);
						q_0_coeff_1[i] = static_cast<ParamType>(_binSizes[i + 1]) * (delta_sum_I2 - static_cast<ParamType>(_binSizes[i])) / (delta_sum_I2 - static_cast<ParamType>(_binSizes[i + 2])) / delta_sum_I2;
						q_0_coeff_2[i] = static_cast<ParamType>(_binSizes[i + 1]) * static_cast<ParamType>(_binSizes[i]) / delta_sum_I2 / (delta_sum_I2 - static_cast<ParamType>(_binSizes[i]));
						q_1_coeff_1[i] = static_cast<ParamType>(_binSizes[i]) * (delta_sum_I1 - static_cast<ParamType>(_binSizes[i + 1])) / delta_sum_I1 / (delta_sum_I1 - static_cast<ParamType>(_binSizes[i - 1]));
						q_1_coeff_2[i] = static_cast<ParamType>(_binSizes[i]) * static_cast<ParamType>(_binSizes[i + 1]) / delta_sum_I1 / (delta_sum_I1 - static_cast<ParamType>(_binSizes[i + 1]));
						q_2_coeff_1[i] = static_cast<ParamType>(_binSizes[i]) * (delta_sum_I0 - static_cast<ParamType>(_binSizes[i - 2])) / delta_sum_I0 / (delta_sum_I0 - static_cast<ParamType>(_binSizes[i]));
						q_2_coeff_2[i] = 1.0 + static_cast<ParamType>(_binSizes[i]) / (static_cast<ParamType>(_binSizes[i - 1]) + static_cast<ParamType>(_binSizes[i])) + static_cast<ParamType>(_binSizes[i]) / delta_sum_I0;
						C_0[i] = delta_sum_I0 * (delta_sum_I0 - static_cast<ParamType>(_binSizes[i - 2])) / delta_sum / (delta_sum - static_cast<ParamType>(_binSizes[i - 2]));
						C_1[i] = delta_sum_I0 / delta_sum * (delta_sum_I2 - static_cast<ParamType>(_binSizes[i])) / (static_cast<ParamType>(_binSizes[i - 1]) + delta_sum_I2) * (1.0 + (delta_sum - static_cast<ParamType>(_binSizes[i - 2])) / (delta_sum - static_cast<ParamType>(_binSizes[i + 2])));
						C_2[i] = static_cast<ParamType>(_binSizes[i + 1]) * (static_cast<ParamType>(_binSizes[i + 1]) + static_cast<ParamType>(_binSizes[i + 2])) / delta_sum / (delta_sum - static_cast<ParamType>(_binSizes[i + 2]));
						IS_0_pre = 4.0 * (static_cast<ParamType>(_binSizes[i]) / delta_sum_I2) * (static_cast<ParamType>(_binSizes[i]) / delta_sum_I2);
						IS_0_coeff_1[i] = IS_0_pre * (10.0 * static_cast<ParamType>(_binSizes[i]) * static_cast<ParamType>(_binSizes[i]) + static_cast<ParamType>(_binSizes[i + 1]) * (static_cast<ParamType>(_binSizes[i]) + static_cast<ParamType>(_binSizes[i + 1]))) / (static_cast<ParamType>(_binSizes[i + 1]) + static_cast<ParamType>(_binSizes[i + 2])) / (static_cast<ParamType>(_binSizes[i + 1]) + static_cast<ParamType>(_binSizes[i + 2]));
						IS_0_coeff_2[i] = IS_0_pre * (20.0 * static_cast<ParamType>(_binSizes[i]) * static_cast<ParamType>(_binSizes[i]) + 2.0 * static_cast<ParamType>(_binSizes[i + 1]) * (static_cast<ParamType>(_binSizes[i]) + static_cast<ParamType>(_binSizes[i + 1])) + (2.0 * static_cast<ParamType>(_binSizes[i + 1]) + static_cast<ParamType>(_binSizes[i])) * delta_sum_I2) / (static_cast<ParamType>(_binSizes[i + 1]) + static_cast<ParamType>(_binSizes[i + 2])) / (static_cast<ParamType>(_binSizes[i + 1]) + static_cast<ParamType>(_binSizes[i]));
						IS_0_coeff_3[i] = IS_0_pre * (10.0 * static_cast<ParamType>(_binSizes[i]) * static_cast<ParamType>(_binSizes[i]) + (2.0 * delta_sum_I2 - static_cast<ParamType>(_binSizes[i + 2])) * (delta_sum_I2 + static_cast<ParamType>(_binSizes[i + 1]))) / (static_cast<ParamType>(_binSizes[i]) + static_cast<ParamType>(_binSizes[i + 1])) / (static_cast<ParamType>(_binSizes[i]) + static_cast<ParamType>(_binSizes[i + 1]));
						IS_1_pre = 4.0 * (static_cast<ParamType>(_binSizes[i]) / delta_sum_I1) * (static_cast<ParamType>(_binSizes[i]) / delta_sum_I1);
						IS_1_coeff_1[i] = IS_1_pre * (10.0 * static_cast<ParamType>(_binSizes[i]) * static_cast<ParamType>(_binSizes[i]) + static_cast<ParamType>(_binSizes[i + 1]) * (static_cast<ParamType>(_binSizes[i]) + static_cast<ParamType>(_binSizes[i + 1]))) / (static_cast<ParamType>(_binSizes[i - 1]) + static_cast<ParamType>(_binSizes[i])) / (static_cast<ParamType>(_binSizes[i - 1]) + static_cast<ParamType>(_binSizes[i]));
						IS_1_coeff_2[i] = IS_1_pre * (20.0 * static_cast<ParamType>(_binSizes[i]) * static_cast<ParamType>(_binSizes[i]) - static_cast<ParamType>(_binSizes[i + 1]) * static_cast<ParamType>(_binSizes[i - 1]) - (static_cast<ParamType>(_binSizes[i]) + static_cast<ParamType>(_binSizes[i + 1])) * (static_cast<ParamType>(_binSizes[i]) + static_cast<ParamType>(_binSizes[i - 1]))) / (static_cast<ParamType>(_binSizes[i]) + static_cast<ParamType>(_binSizes[i + 1])) / (static_cast<ParamType>(_binSizes[i]) + static_cast<ParamType>(_binSizes[i - 1]));
						IS_1_coeff_3[i] = IS_1_pre * (10.0 * static_cast<ParamType>(_binSizes[i]) * static_cast<ParamType>(_binSizes[i]) + static_cast<ParamType>(_binSizes[i - 1]) * (static_cast<ParamType>(_binSizes[i - 1]) + static_cast<ParamType>(_binSizes[i]))) / (static_cast<ParamType>(_binSizes[i]) + static_cast<ParamType>(_binSizes[i + 1])) / (static_cast<ParamType>(_binSizes[i]) + static_cast<ParamType>(_binSizes[i + 1]));
						IS_2_pre = 4.0 * (static_cast<ParamType>(_binSizes[i]) / delta_sum_I0) * (static_cast<ParamType>(_binSizes[i]) / delta_sum_I0);
						IS_2_coeff_1[i] = IS_2_pre * (10.0 * static_cast<ParamType>(_binSizes[i]) * static_cast<ParamType>(_binSizes[i]) + static_cast<ParamType>(_binSizes[i - 1]) * (static_cast<ParamType>(_binSizes[i - 1]) + static_cast<ParamType>(_binSizes[i]))) / (static_cast<ParamType>(_binSizes[i - 2]) + static_cast<ParamType>(_binSizes[i - 1])) / (static_cast<ParamType>(_binSizes[i - 2]) + static_cast<ParamType>(_binSizes[i - 1]));
						IS_2_coeff_2[i] = IS_2_pre * (20.0 * static_cast<ParamType>(_binSizes[i]) * static_cast<ParamType>(_binSizes[i]) + 2.0 * static_cast<ParamType>(_binSizes[i - 1]) * (static_cast<ParamType>(_binSizes[i - 1]) + static_cast<ParamType>(_binSizes[i])) + delta_sum_I0 * (2.0 * static_cast<ParamType>(_binSizes[i - 1]) + static_cast<ParamType>(_binSizes[i]))) / (static_cast<ParamType>(_binSizes[i - 1]) + static_cast<ParamType>(_binSizes[i])) / (static_cast<ParamType>(_binSizes[i - 2]) + static_cast<ParamType>(_binSizes[i - 1]));
						IS_2_coeff_3[i] = IS_2_pre * (10.0 * static_cast<ParamType>(_binSizes[i]) * static_cast<ParamType>(_binSizes[i]) + (2.0 * delta_sum_I0 - static_cast<ParamType>(_binSizes[i - 2])) * (delta_sum_I0 + static_cast<ParamType>(_binSizes[i - 1]))) / (static_cast<ParamType>(_binSizes[i - 1]) + static_cast<ParamType>(_binSizes[i])) / (static_cast<ParamType>(_binSizes[i - 1]) + static_cast<ParamType>(_binSizes[i]));
					}
					// define local variables
					StateParam IS_0 = 0.0;
					StateParam IS_1 = 0.0;
					StateParam IS_2 = 0.0;
					StateParam alpha_0 = 0.0;
					StateParam alpha_1 = 0.0;
					StateParam alpha_2 = 0.0;
					StateParam W_0 = 0.0;
					StateParam W_1 = 0.0;
					StateParam W_2 = 0.0;
					StateParam q_0 = 0.0;
					StateParam q_1 = 0.0;
					StateParam q_2 = 0.0;
					for (int i = 0; i < _nBins; ++i)
					{
						if (cadet_likely((i > 2) && (i + 2 < _nBins)))
						{
							// flux through the left face 
							IS_0 = static_cast<ParamType>(IS_0_coeff_1[i - 1]) * (yCrystal[i + 1] - yCrystal[i]) * (yCrystal[i + 1] - yCrystal[i]) + static_cast<ParamType>(IS_0_coeff_2[i - 1]) * (yCrystal[i + 1] - yCrystal[i]) * (yCrystal[i - 1] - yCrystal[i]) + static_cast<ParamType>(IS_0_coeff_3[i - 1]) * (yCrystal[i - 1] - yCrystal[i]) * (yCrystal[i - 1] - yCrystal[i]);
							IS_1 = static_cast<ParamType>(IS_1_coeff_1[i - 1]) * (yCrystal[i - 2] - yCrystal[i - 1]) * (yCrystal[i - 2] - yCrystal[i - 1]) + static_cast<ParamType>(IS_1_coeff_2[i - 1]) * (yCrystal[i] - yCrystal[i - 1]) * (yCrystal[i - 2] - yCrystal[i - 1]) + static_cast<ParamType>(IS_1_coeff_3[i - 1]) * (yCrystal[i] - yCrystal[i - 1]) * (yCrystal[i] - yCrystal[i - 1]);
							IS_2 = static_cast<ParamType>(IS_2_coeff_1[i - 1]) * (yCrystal[i - 3] - yCrystal[i - 2]) * (yCrystal[i - 3] - yCrystal[i - 2]) + static_cast<ParamType>(IS_2_coeff_2[i - 1]) * (yCrystal[i - 1] - yCrystal[i - 2]) * (yCrystal[i - 3] - yCrystal[i - 2]) + static_cast<ParamType>(IS_2_coeff_3[i - 1]) * (yCrystal[i - 1] - yCrystal[i - 2]) * (yCrystal[i - 1] - yCrystal[i - 2]);
							alpha_0 = static_cast<ParamType>(C_0[i - 1]) / (1e-10 + IS_0) / (1e-10 + IS_0);
							alpha_1 = static_cast<ParamType>(C_1[i - 1]) / (1e-10 + IS_1) / (1e-10 + IS_1);
							alpha_2 = static_cast<ParamType>(C_2[i - 1]) / (1e-10 + IS_2) / (1e-10 + IS_2);
							W_0 = alpha_0 / (alpha_0 + alpha_1 + alpha_2);
							W_1 = alpha_1 / (alpha_0 + alpha_1 + alpha_2);
							W_2 = 1.0 - W_0 - W_1;
							q_0 = yCrystal[i] + static_cast<ParamType>(q_0_coeff_1[i - 1]) * (yCrystal[i - 1] - yCrystal[i]) + static_cast<ParamType>(q_0_coeff_2[i - 1]) * (yCrystal[i] - yCrystal[i + 1]);
							q_1 = yCrystal[i - 1] + static_cast<ParamType>(q_1_coeff_1[i - 1]) * (yCrystal[i] - yCrystal[i - 1]) + static_cast<ParamType>(q_1_coeff_2[i - 1]) * (yCrystal[i - 1] - yCrystal[i - 2]);
							q_2 = yCrystal[i - 2] + static_cast<ParamType>(q_2_coeff_1[i - 1]) * (yCrystal[i - 3] - yCrystal[i - 2]) + static_cast<ParamType>(q_2_coeff_2[i - 1]) * (yCrystal[i - 1] - yCrystal[i - 2]);
							resCrystal[i] += factor / static_cast<ParamType>(_binSizes[i]) * v_g * (W_0 * q_0 + W_1 * q_1 + W_2 * q_2) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i] - yCrystal[i - 1]) / static_cast<ParamType>(_binCenterDists[i - 1]);
							// flux through the right face
							IS_0 = static_cast<ParamType>(IS_0_coeff_1[i]) * (yCrystal[i + 2] - yCrystal[i + 1]) * (yCrystal[i + 2] - yCrystal[i + 1]) + static_cast<ParamType>(IS_0_coeff_2[i]) * (yCrystal[i + 2] - yCrystal[i + 1]) * (yCrystal[i] - yCrystal[i + 1]) + static_cast<ParamType>(IS_0_coeff_3[i]) * (yCrystal[i] - yCrystal[i + 1]) * (yCrystal[i] - yCrystal[i + 1]);
							IS_1 = static_cast<ParamType>(IS_1_coeff_1[i]) * (yCrystal[i - 1] - yCrystal[i]) * (yCrystal[i - 1] - yCrystal[i]) + static_cast<ParamType>(IS_1_coeff_2[i]) * (yCrystal[i + 1] - yCrystal[i]) * (yCrystal[i - 1] - yCrystal[i]) + static_cast<ParamType>(IS_1_coeff_3[i]) * (yCrystal[i + 1] - yCrystal[i]) * (yCrystal[i + 1] - yCrystal[i]);
							IS_2 = static_cast<ParamType>(IS_2_coeff_1[i]) * (yCrystal[i - 2] - yCrystal[i - 1]) * (yCrystal[i - 2] - yCrystal[i - 1]) + static_cast<ParamType>(IS_2_coeff_2[i]) * (yCrystal[i] - yCrystal[i - 1]) * (yCrystal[i - 2] - yCrystal[i - 1]) + static_cast<ParamType>(IS_2_coeff_3[i]) * (yCrystal[i] - yCrystal[i - 1]) * (yCrystal[i] - yCrystal[i - 1]);
							alpha_0 = static_cast<ParamType>(C_0[i]) / (1e-10 + IS_0) / (1e-10 + IS_0);
							alpha_1 = static_cast<ParamType>(C_1[i]) / (1e-10 + IS_1) / (1e-10 + IS_1);
							alpha_2 = static_cast<ParamType>(C_2[i]) / (1e-10 + IS_2) / (1e-10 + IS_2);
							W_0 = alpha_0 / (alpha_0 + alpha_1 + alpha_2);
							W_1 = alpha_1 / (alpha_0 + alpha_1 + alpha_2);
							W_2 = 1.0 - W_0 - W_1;
							q_0 = yCrystal[i + 1] + static_cast<ParamType>(q_0_coeff_1[i]) * (yCrystal[i] - yCrystal[i + 1]) + static_cast<ParamType>(q_0_coeff_2[i]) * (yCrystal[i + 1] - yCrystal[i + 2]);
							q_1 = yCrystal[i] + static_cast<ParamType>(q_1_coeff_1[i]) * (yCrystal[i + 1] - yCrystal[i]) + static_cast<ParamType>(q_1_coeff_2[i]) * (yCrystal[i] - yCrystal[i - 1]);
							q_2 = yCrystal[i - 1] + static_cast<ParamType>(q_2_coeff_1[i]) * (yCrystal[i - 2] - yCrystal[i - 1]) + static_cast<ParamType>(q_2_coeff_2[i]) * (yCrystal[i] - yCrystal[i - 1]);
							v_g = k_g_times_s_g * (static_cast<ParamType>(_a) + static_cast<ParamType>(_growthConstant) * pow(static_cast<ParamType>(_bins[i + 1]), static_cast<ParamType>(_p)));
							resCrystal[i] -= factor / static_cast<ParamType>(_binSizes[i]) * v_g * (W_0 * q_0 + W_1 * q_1 + W_2 * q_2) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i + 1] - yCrystal[i]) / static_cast<ParamType>(_binCenterDists[i]);
						}
						else if ((i == 1) || (i + 2 == _nBins) || (i == 2))
						{
							// first order approximation
							resCrystal[i] += factor / static_cast<ParamType>(_binSizes[i]) * v_g * (0.5 * yCrystal[i - 1] + 0.5 * yCrystal[i]) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i] - yCrystal[i - 1]) / static_cast<ParamType>(_binCenterDists[i - 1]);
							v_g = k_g_times_s_g * (static_cast<ParamType>(_a) + static_cast<ParamType>(_growthConstant) * pow(static_cast<ParamType>(_bins[i + 1]), static_cast<ParamType>(_p)));
							resCrystal[i] -= factor / static_cast<ParamType>(_binSizes[i]) * v_g * (0.5 * yCrystal[i] + 0.5 * yCrystal[i + 1]) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i + 1] - yCrystal[i]) / static_cast<ParamType>(_binCenterDists[i]);
						}
						else if (i == 0)
						{
							// BC
							resCrystal[i] += factor / static_cast<ParamType>(_binSizes[i]) * B_0;
							v_g = k_g_times_s_g * (static_cast<ParamType>(_a) + static_cast<ParamType>(_growthConstant) * pow(static_cast<ParamType>(_bins[i + 1]), static_cast<ParamType>(_p)));
							resCrystal[i] -= factor / static_cast<ParamType>(_binSizes[i]) * v_g * (0.5 * yCrystal[i] + 0.5 * yCrystal[i + 1]) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i + 1] - yCrystal[i]) / static_cast<ParamType>(_binCenterDists[i]);
						}
						else
						{
							// first order approximation
							resCrystal[i] += factor / static_cast<ParamType>(_binSizes[i]) * v_g * (0.5 * yCrystal[i - 1] + 0.5 * yCrystal[i]) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i] - yCrystal[i - 1]) / static_cast<ParamType>(_binCenterDists[i - 1]);
							// no flux BC
						}
					}

				}
				// compact weno schemes, still complicated, but much better than weno5
				else if (static_cast<ParamType>(_growthSchemeOrder) == 5)
				{
					std::vector<active> P_C_R_coeff_1(_nBins - 1);
					std::vector<active> P_C_R_coeff_2(_nBins - 1);  // negative 
					std::vector<active> P_L_coeff(_nBins - 1);
					std::vector<active> P_R_coeff(_nBins - 1);
					std::vector<active> IS_L_coeff(_nBins - 1);
					std::vector<active> IS_R_coeff(_nBins - 1);
					std::vector<active> IS_C_coeff_1(_nBins - 1);
					std::vector<active> IS_C_coeff_2(_nBins - 1);
					std::vector<active> IS_C_coeff_3(_nBins - 1);
					ParamType delta_sum = 0.0;
					ParamType delta_0_p1 = 0.0;
					ParamType delta_0_n1 = 0.0;
					// calculate the coefficients and store them, the first entry is 0
					for (int i = 1; i + 1 < _nBins; ++i)
					{
						delta_sum = static_cast<ParamType>(_binSizes[i]) + static_cast<ParamType>(_binSizes[i - 1]) + static_cast<ParamType>(_binSizes[i + 1]);
						delta_0_p1 = static_cast<ParamType>(_binSizes[i]) + static_cast<ParamType>(_binSizes[i + 1]);
						delta_0_n1 = static_cast<ParamType>(_binSizes[i]) + static_cast<ParamType>(_binSizes[i - 1]);
						P_C_R_coeff_1[i] = static_cast<ParamType>(_binSizes[i]) * (4.0 * delta_0_n1 - delta_sum) / 2.0 / delta_sum / delta_0_p1;
						P_C_R_coeff_2[i] = -static_cast<ParamType>(_binSizes[i]) * (delta_sum - 4.0 * static_cast<ParamType>(_binSizes[i + 1])) / 2.0 / delta_sum / delta_0_n1;
						P_L_coeff[i] = static_cast<ParamType>(_binSizes[i]) / delta_0_n1;
						P_R_coeff[i] = static_cast<ParamType>(_binSizes[i]) / delta_0_p1;
						IS_L_coeff[i] = 4.0 * P_L_coeff[i] * P_L_coeff[i];
						IS_R_coeff[i] = 4.0 * P_R_coeff[i] * P_R_coeff[i];
						IS_C_coeff_1[i] = (static_cast<ParamType>(_binSizes[i]) * static_cast<ParamType>(_binSizes[i]) * (static_cast<ParamType>(_binSizes[i]) + 3.0 * static_cast<ParamType>(_binSizes[i - 1]) - static_cast<ParamType>(_binSizes[i + 1])) * (static_cast<ParamType>(_binSizes[i]) + 3.0 * static_cast<ParamType>(_binSizes[i - 1]) - static_cast<ParamType>(_binSizes[i + 1])) + 156.0 * static_cast<ParamType>(_binSizes[i]) * static_cast<ParamType>(_binSizes[i]) * static_cast<ParamType>(_binSizes[i]) * static_cast<ParamType>(_binSizes[i])) / delta_sum / delta_sum / delta_0_p1 / delta_0_p1;
						IS_C_coeff_2[i] = (static_cast<ParamType>(_binSizes[i]) * static_cast<ParamType>(_binSizes[i]) * (static_cast<ParamType>(_binSizes[i]) + 3.0 * static_cast<ParamType>(_binSizes[i + 1]) - static_cast<ParamType>(_binSizes[i - 1])) * (static_cast<ParamType>(_binSizes[i]) + 3.0 * static_cast<ParamType>(_binSizes[i + 1]) - static_cast<ParamType>(_binSizes[i - 1])) + 156.0 * static_cast<ParamType>(_binSizes[i]) * static_cast<ParamType>(_binSizes[i]) * static_cast<ParamType>(_binSizes[i]) * static_cast<ParamType>(_binSizes[i])) / delta_sum / delta_sum / delta_0_n1 / delta_0_n1;
						IS_C_coeff_3[i] = (2.0 * static_cast<ParamType>(_binSizes[i]) * static_cast<ParamType>(_binSizes[i]) * (static_cast<ParamType>(_binSizes[i]) + 3.0 * static_cast<ParamType>(_binSizes[i + 1]) - static_cast<ParamType>(_binSizes[i - 1])) * (static_cast<ParamType>(_binSizes[i]) + 3.0 * static_cast<ParamType>(_binSizes[i - 1]) - static_cast<ParamType>(_binSizes[i + 1])) + 312.0 * static_cast<ParamType>(_binSizes[i]) * static_cast<ParamType>(_binSizes[i]) * static_cast<ParamType>(_binSizes[i]) * static_cast<ParamType>(_binSizes[i])) / delta_sum / delta_sum / delta_0_n1 / delta_0_p1;
					}
					// ode input
					StateParam IS_L_square = 0.0;
					StateParam IS_R_square = 0.0;
					StateParam IS_C_square = 0.0;
					StateParam alpha_L = 0.0;
					StateParam alpha_R = 0.0;
					StateParam alpha_C = 0.0;
					StateParam W_L = 0.0;
					StateParam W_R = 0.0;
					StateParam W_C = 0.0;
					StateParam P_L = 0.0;
					StateParam P_R = 0.0;
					StateParam P_C = 0.0;
					for (int i = 0; i < _nBins; ++i)
					{
						if (cadet_likely((i > 1) && (i + 1 < _nBins)))
						{
							// flux through the left face 
							IS_L_square = static_cast<ParamType>(IS_L_coeff[i - 1]) * (yCrystal[i - 1] - yCrystal[i - 2]) * (yCrystal[i - 1] - yCrystal[i - 2]);
							IS_R_square = static_cast<ParamType>(IS_R_coeff[i - 1]) * (yCrystal[i] - yCrystal[i - 1]) * (yCrystal[i] - yCrystal[i - 1]);
							IS_C_square = static_cast<ParamType>(IS_C_coeff_1[i - 1]) * (yCrystal[i] - yCrystal[i - 1]) * (yCrystal[i] - yCrystal[i - 1]) + static_cast<ParamType>(IS_C_coeff_2[i - 1]) * (yCrystal[i - 1] - yCrystal[i - 2]) * (yCrystal[i - 1] - yCrystal[i - 2]) + static_cast<ParamType>(IS_C_coeff_3[i - 1]) * (yCrystal[i - 1] - yCrystal[i - 2]) * (yCrystal[i] - yCrystal[i - 1]);
							alpha_L = 0.25 / (1e-10 + IS_L_square) / (1e-10 + IS_L_square);
							alpha_R = 0.25 / (1e-10 + IS_R_square) / (1e-10 + IS_R_square);
							alpha_C = 0.5 / (1e-10 + IS_C_square) / (1e-10 + IS_C_square);
							W_L = alpha_L / (alpha_L + alpha_R + alpha_C);
							W_R = alpha_R / (alpha_L + alpha_R + alpha_C);
							W_C = 1.0 - W_L - W_R;
							P_L = yCrystal[i - 1] + static_cast<ParamType>(P_L_coeff[i - 1]) * (yCrystal[i - 1] - yCrystal[i - 2]);
							P_R = yCrystal[i - 1] + static_cast<ParamType>(P_R_coeff[i - 1]) * (yCrystal[i] - yCrystal[i - 1]);
							P_C = yCrystal[i - 1] + static_cast<ParamType>(P_C_R_coeff_1[i]) * (yCrystal[i] - yCrystal[i - 1]) + static_cast<ParamType>(P_C_R_coeff_2[i]) * (yCrystal[i - 1] - yCrystal[i - 2]);
							resCrystal[i] += factor / static_cast<ParamType>(_binSizes[i]) * v_g * (W_L * P_L + W_C * P_C + W_R * P_R) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i] - yCrystal[i - 1]) / static_cast<ParamType>(_binCenterDists[i - 1]);
							// flux through the right face
							IS_L_square = static_cast<ParamType>(IS_L_coeff[i]) * (yCrystal[i] - yCrystal[i - 1]) * (yCrystal[i] - yCrystal[i - 1]);
							IS_R_square = static_cast<ParamType>(IS_R_coeff[i]) * (yCrystal[i + 1] - yCrystal[i]) * (yCrystal[i + 1] - yCrystal[i]);
							IS_C_square = static_cast<ParamType>(IS_C_coeff_1[i]) * (yCrystal[i + 1] - yCrystal[i]) * (yCrystal[i + 1] - yCrystal[i]) + static_cast<ParamType>(IS_C_coeff_2[i]) * (yCrystal[i] - yCrystal[i - 1]) * (yCrystal[i] - yCrystal[i - 1]) + static_cast<ParamType>(IS_C_coeff_3[i]) * (yCrystal[i] - yCrystal[i - 1]) * (yCrystal[i + 1] - yCrystal[i]);
							alpha_L = 0.25 / (1e-10 + IS_L_square) / (1e-10 + IS_L_square);
							alpha_R = 0.25 / (1e-10 + IS_R_square) / (1e-10 + IS_R_square);
							alpha_C = 0.5 / (1e-10 + IS_C_square) / (1e-10 + IS_C_square);
							W_L = alpha_L / (alpha_L + alpha_R + alpha_C);
							W_R = alpha_R / (alpha_L + alpha_R + alpha_C);
							W_C = 1.0 - W_L - W_R;
							P_L = yCrystal[i] + static_cast<ParamType>(P_L_coeff[i]) * (yCrystal[i] - yCrystal[i - 1]);
							P_R = yCrystal[i] + static_cast<ParamType>(P_R_coeff[i]) * (yCrystal[i + 1] - yCrystal[i]);
							P_C = yCrystal[i] + static_cast<ParamType>(P_C_R_coeff_1[i]) * (yCrystal[i + 1] - yCrystal[i]) + static_cast<ParamType>(P_C_R_coeff_2[i]) * (yCrystal[i] - yCrystal[i - 1]);
							v_g = k_g_times_s_g * (static_cast<ParamType>(_a) + static_cast<ParamType>(_growthConstant) * pow(static_cast<ParamType>(_bins[i + 1]), static_cast<ParamType>(_p)));
							resCrystal[i] -= factor / static_cast<ParamType>(_binSizes[i]) * v_g * (W_L * P_L + W_C * P_C + W_R * P_R) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i + 1] - yCrystal[i]) / static_cast<ParamType>(_binCenterDists[i]);
						}
						else if (i + 1 == _nBins)
						{
							// first order approximation
							resCrystal[i] += factor / static_cast<ParamType>(_binSizes[i]) * v_g * (0.5 * yCrystal[i - 1] + 0.5 * yCrystal[i]) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i] - yCrystal[i - 1]) / static_cast<ParamType>(_binCenterDists[i - 1]);
							// no flux BC
						}
						else if (i == 1)
						{
							// ideally cweno scheme should be applied to F_{ 1 + 1 / 2 }, but this saves some code and logical if statement
							resCrystal[i] += factor / static_cast<ParamType>(_binSizes[i]) * v_g * (0.5 * yCrystal[i - 1] + 0.5 * yCrystal[i]) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i] - yCrystal[i - 1]) / static_cast<ParamType>(_binCenterDists[i - 1]);
							v_g = k_g_times_s_g * (static_cast<ParamType>(_a) + static_cast<ParamType>(_growthConstant) * pow(static_cast<ParamType>(_bins[i + 1]), static_cast<ParamType>(_p)));
							resCrystal[i] -= factor / static_cast<ParamType>(_binSizes[i]) * v_g * (0.5 * yCrystal[i] + 0.5 * yCrystal[i + 1]) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i + 1] - yCrystal[i]) / static_cast<ParamType>(_binCenterDists[i]);
						}
						else
						{
							// boundary condition
							resCrystal[i] += factor / static_cast<ParamType>(_binSizes[i]) * B_0;
							// first order approximation
							v_g = k_g_times_s_g * (static_cast<ParamType>(_a) + static_cast<ParamType>(_growthConstant) * pow(static_cast<ParamType>(_bins[i + 1]), static_cast<ParamType>(_p)));
							resCrystal[i] -= factor / static_cast<ParamType>(_binSizes[i]) * v_g * (0.5 * yCrystal[i] + 0.5 * yCrystal[i + 1]) - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i + 1] - yCrystal[i]) / static_cast<ParamType>(_binCenterDists[i]);
						}
					}
				}
				return 0;
			}

			template <typename StateType, typename ResidualType, typename ParamType>
			int residualCombinedImpl(double t, unsigned int secIdx, const ColumnPosition& colPos,
				StateType const* yLiquid, StateType const* ySolid, ResidualType* resLiquid, ResidualType* resSolid, double factor, LinearBufferAllocator workSpace) const
			{
				return residualLiquidImpl<StateType, ResidualType, ParamType, double>(t, secIdx, colPos, yLiquid, resLiquid, factor, workSpace);
			}

			template <typename RowIterator>
			void jacobianLiquidImpl(double t, unsigned int secIdx, const ColumnPosition& colPos, double const* y, double factor, const RowIterator& jac, LinearBufferAllocator workSpace) const
			{
			}

			template <typename RowIteratorLiquid, typename RowIteratorSolid>
			void jacobianCombinedImpl(double t, unsigned int secIdx, const ColumnPosition& colPos, double const* yLiquid, double const* ySolid, double factor, const RowIteratorLiquid& jacLiquid, const RowIteratorSolid& jacSolid, LinearBufferAllocator workSpace) const
			{
			}
		};

		namespace reaction
		{
			void registerCrystallizationReaction(std::unordered_map<std::string, std::function<model::IDynamicReactionModel* ()>>& reactions)
			{
				reactions[CrystallizationReaction::identifier()] = []() { return new CrystallizationReaction(); };
			}
		}  // namespace reaction

	}  // namespace model

}  // namespace cadet
