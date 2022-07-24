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

//#include "LoggingUtils.hpp"
//#include "Logging.hpp"

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

				_breakageRateConstant = paramProvider.getDouble("CRY_BREAKAGE_RATE_CONSTANT");
				_parameters[makeParamId(hashString("CRY_BREAKAGE_RATE_CONSTANT"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_breakageRateConstant;

				_breakageKernelGamma = paramProvider.getDouble("CRY_BREAKAGE_KERNEL_GAMMA");
				_parameters[makeParamId(hashString("CRY_BREAKAGE_KERNEL_GAMMA"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_breakageKernelGamma;

				_breakageSelectionFunctionAlpha = paramProvider.getDouble("CRY_BREAKAGE_SELECTION_FUNCTION_ALPHA");
				_parameters[makeParamId(hashString("CRY_BREAKAGE_SELECTION_FUNCTION_ALPHA"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_breakageSelectionFunctionAlpha;

				_aggregationRateConstant = paramProvider.getDouble("CRY_AGGREGATION_RATE_CONSTANT");
				_parameters[makeParamId(hashString("CRY_AGGREGATION_RATE_CONSTANT"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_aggregationRateConstant;

				_aggregationIndex = paramProvider.getInt("CRY_AGGREGATION_INDEX");

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
			active _breakageRateConstant; // constant breakage rate constant
			active _breakageKernelGamma; // gamma in the breakage kernel
			active _breakageSelectionFunctionAlpha; // alpha in the selection function
			active _aggregationIndex; // determines which kernel to use
			active _aggregationRateConstant; // constant aggregation rate constant
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
					substrateConversion += yCrystal[i] * static_cast<ParamType>(_binCenters[i]) * static_cast<ParamType>(_binCenters[i]) * pow(1.0 + static_cast<ParamType>(_growthConstant) * static_cast<ParamType>(_binCenters[i]), static_cast<ParamType>(_p)) * static_cast<ParamType>(_binSizes[i]);
				}
				M *= massDensityShapeFactor;
				substrateConversion *= 3.0 * k_g_times_s_g;

				// B_0 = primary + secondary nucleation rate
				const StateParam B_0 = static_cast<ParamType>(_primaryNucleationRate) * pow(s, static_cast<ParamType>(_u)) + static_cast<ParamType>(_secondaryNucleationRate) * pow(s, static_cast<ParamType>(_b)) * pow(M, static_cast<ParamType>(_k));
				const ParamType x_c_3 = static_cast<ParamType>(_bins[0]) * static_cast<ParamType>(_bins[0]) * static_cast<ParamType>(_bins[0]);

				// disable the mass balance equation first
				//res[0] -= factor * massDensityShapeFactor * (B_0 * x_c_3 + substrateConversion);

				StateParam v_g = 0.0;

				// define breakage-related local parameters
				std::vector<active> Upsilon_source(_nBins);
				std::vector<active> Upsilon_sink(_nBins);

				// these parameters should be ParamType as they are independent on the StateType
				// StateType are reserved for y or res
				// StateParam are reserved for a combination of ParamType and StateType
				ParamType Upsilon_birth_sum = 0.0;
				ParamType Upsilon_death_sum = 0.0;
				const ParamType N_j = static_cast<ParamType>(_breakageKernelGamma) / (static_cast<ParamType>(_breakageKernelGamma) - 1);
				ParamType b_integral_birth = 0.0;
				// calculate upsilon and store it
				for (int i = 0; i < _nBins; ++i)
				{
					// reset the sum for each i
					Upsilon_birth_sum = 0.0;
					Upsilon_death_sum = 0.0;
					for (int j = 0; j < i + 1; ++j)
					{
						if (cadet_likely(i != j))
						{
							b_integral_birth = N_j * (pow(static_cast<ParamType>(_bins[j + 1]), 3.0 * static_cast<ParamType>(_breakageKernelGamma) - 3.0) - pow(static_cast<ParamType>(_bins[j]), 3.0 * static_cast<ParamType>(_breakageKernelGamma) - 3.0)) / pow(static_cast<ParamType>(_binCenters[i]), 3.0 * static_cast<ParamType>(_breakageKernelGamma) - 3.0);
							Upsilon_birth_sum += b_integral_birth * (static_cast<ParamType>(_binCenters[i]) * static_cast<ParamType>(_binCenters[i]) * static_cast<ParamType>(_binCenters[i]) - static_cast<ParamType>(_binCenters[j]) * static_cast<ParamType>(_binCenters[j]) * static_cast<ParamType>(_binCenters[j]));
							Upsilon_death_sum += b_integral_birth * static_cast<ParamType>(_binCenters[j]) * static_cast<ParamType>(_binCenters[j]) * static_cast<ParamType>(_binCenters[j]);
						}
						else
						{
							b_integral_birth = N_j * (pow(static_cast<ParamType>(_binCenters[j]), 3.0 * static_cast<ParamType>(_breakageKernelGamma) - 3.0) - pow(static_cast<ParamType>(_bins[j]), 3.0 * static_cast<ParamType>(_breakageKernelGamma) - 3.0)) / pow(static_cast<ParamType>(_binCenters[i]), 3.0 * static_cast<ParamType>(_breakageKernelGamma) - 3.0);
							Upsilon_death_sum += b_integral_birth * static_cast<ParamType>(_binCenters[j]) * static_cast<ParamType>(_binCenters[j]) * static_cast<ParamType>(_binCenters[j]);
						}
					}
					if (cadet_likely(i > 0))
					{
						Upsilon_source[i] = (N_j - 1.0) * (static_cast<ParamType>(_binCenters[i]) * static_cast<ParamType>(_binCenters[i]) * static_cast<ParamType>(_binCenters[i])) / Upsilon_birth_sum;
					}
					else
					{
						Upsilon_source[i] = (N_j - 1.0) * (static_cast<ParamType>(_binCenters[i]) * static_cast<ParamType>(_binCenters[i]) * static_cast<ParamType>(_binCenters[i]));
					}
					Upsilon_sink[i] = Upsilon_death_sum * Upsilon_source[i] / static_cast<ParamType>(_binCenters[i]) / static_cast<ParamType>(_binCenters[i]) / static_cast<ParamType>(_binCenters[i]);
				}
				// ode input
				// define breakage-related local parameters
				ParamType b_integral = 0.0;    
				ParamType selection_function = 0.0;  
				StateParam breakage_source = 0.0;
				StateParam breakage_sink = 0.0;
				// define aggregation-related local parameters
				StateParam sum_volume = 0.0;
				StateParam volume_k = 0.0;
				StateParam volume_j = 0.0;
				StateParam aggregation_correction_factor = 0.0;
				StateParam aggregation_source = 0.0;
				StateParam aggregation_sink = 0.0;
				for (int i = 0; i < _nBins; ++i)
				{
					// Flux through left face, reset the source term
					breakage_source = 0.0;
					for (int j = i; j < _nBins; ++j)
					{
						// selection function
						selection_function = static_cast<ParamType>(_breakageRateConstant) * pow(static_cast<ParamType>(_binCenters[j]), 3.0 * static_cast<ParamType>(_breakageSelectionFunctionAlpha));
						if (cadet_likely(i != j))
						{
							b_integral = N_j * (pow(static_cast<ParamType>(_bins[i + 1]), 3.0 * static_cast<ParamType>(_breakageKernelGamma) - 3.0) - pow(static_cast<ParamType>(_bins[i]), 3.0 * static_cast<ParamType>(_breakageKernelGamma) - 3.0)) / pow(static_cast<ParamType>(_binCenters[j]), 3.0 * static_cast<ParamType>(_breakageKernelGamma) - 3.0);
						}
						else
						{
							b_integral = N_j * (pow(static_cast<ParamType>(_binCenters[i]), 3.0 * static_cast<ParamType>(_breakageKernelGamma) - 3.0) - pow(static_cast<ParamType>(_bins[i]), 3.0 * static_cast<ParamType>(_breakageKernelGamma) - 3.0)) / pow(static_cast<ParamType>(_binCenters[j]), 3.0 * static_cast<ParamType>(_breakageKernelGamma) - 3.0);
						}
						breakage_source += selection_function * yCrystal[j] * b_integral * static_cast<ParamType>(Upsilon_source[j]) * static_cast<ParamType>(_binSizes[j]) / static_cast<ParamType>(_binSizes[i]);
					}
					aggregation_source = 0.0;
					aggregation_sink = 0.0;
					for (int j = 0; j < i + 1; ++j)
					{
						for (int k = j; k < i + 1; ++k)
						{
							sum_volume = pow(static_cast<ParamType>(_binCenters[k]) * static_cast<ParamType>(_binCenters[k]) * static_cast<ParamType>(_binCenters[k]) + static_cast<ParamType>(_binCenters[j]) * static_cast<ParamType>(_binCenters[j]) * static_cast<ParamType>(_binCenters[j]), 1.0 / 3.0);
							if (cadet_likely((sum_volume > static_cast<ParamType>(_bins[i])) && (sum_volume <= static_cast<ParamType>(_bins[i + 1]))))
							{
								//calculate the correction factor
								aggregation_correction_factor = (static_cast<ParamType>(_binCenters[j]) * static_cast<ParamType>(_binCenters[j]) * static_cast<ParamType>(_binCenters[j]) + static_cast<ParamType>(_binCenters[k]) * static_cast<ParamType>(_binCenters[k]) * static_cast<ParamType>(_binCenters[k])) / static_cast<ParamType>(_binCenters[i]) / static_cast<ParamType>(_binCenters[i]) / static_cast<ParamType>(_binCenters[i]);
								// 0.5 in the original scheme becomes 1 because the j <= k condition removes all duplicates, except for the center j=k
								// constant kernel 0, brownian kernel 1, smoluchowski kernel 2, golovin kernel 3, differential force kernel 4
								if (static_cast<ParamType>(_aggregationIndex) == 0)
								{
									if (cadet_likely(j != k))
									{
										aggregation_source += aggregation_correction_factor * static_cast<ParamType>(_aggregationRateConstant) * yCrystal[j] * yCrystal[k] * static_cast<ParamType>(_binSizes[j]) * static_cast<ParamType>(_binSizes[k]) / static_cast<ParamType>(_binSizes[i]);
									}
									else
									{
										aggregation_source += 0.5 * aggregation_correction_factor * static_cast<ParamType>(_aggregationRateConstant) * yCrystal[j] * yCrystal[k] * static_cast<ParamType>(_binSizes[j]) * static_cast<ParamType>(_binSizes[k]) / static_cast<ParamType>(_binSizes[i]);
									}
								}
								else if (static_cast<ParamType>(_aggregationIndex) == 1)
								{
									if (cadet_likely(j != k))
									{
										aggregation_source += aggregation_correction_factor * static_cast<ParamType>(_aggregationRateConstant) * (static_cast<ParamType>(_binCenters[k]) + static_cast<ParamType>(_binCenters[j])) * (static_cast<ParamType>(_binCenters[k]) + static_cast<ParamType>(_binCenters[j])) / static_cast<ParamType>(_binCenters[k]) / static_cast<ParamType>(_binCenters[j]) * yCrystal[j] * yCrystal[k] * static_cast<ParamType>(_binSizes[j]) * static_cast<ParamType>(_binSizes[k]) / static_cast<ParamType>(_binSizes[i]);
									}
									else
									{
										aggregation_source += 0.5 * aggregation_correction_factor * static_cast<ParamType>(_aggregationRateConstant) * (static_cast<ParamType>(_binCenters[k]) + static_cast<ParamType>(_binCenters[j])) * (static_cast<ParamType>(_binCenters[k]) + static_cast<ParamType>(_binCenters[j])) / static_cast<ParamType>(_binCenters[k]) / static_cast<ParamType>(_binCenters[j]) * yCrystal[j] * yCrystal[k] * static_cast<ParamType>(_binSizes[j]) * static_cast<ParamType>(_binSizes[k]) / static_cast<ParamType>(_binSizes[i]);
									}
								}
								else if (static_cast<ParamType>(_aggregationIndex) == 2)
								{
									if (cadet_likely(j != k))
									{
										aggregation_source += aggregation_correction_factor * static_cast<ParamType>(_aggregationRateConstant) * (static_cast<ParamType>(_binCenters[k]) + static_cast<ParamType>(_binCenters[j])) * (static_cast<ParamType>(_binCenters[k]) + static_cast<ParamType>(_binCenters[j])) * (static_cast<ParamType>(_binCenters[k]) + static_cast<ParamType>(_binCenters[j])) * yCrystal[j] * yCrystal[k] * static_cast<ParamType>(_binSizes[j]) * static_cast<ParamType>(_binSizes[k]) / static_cast<ParamType>(_binSizes[i]);
									}
									else
									{
										aggregation_source += 0.5 * aggregation_correction_factor * static_cast<ParamType>(_aggregationRateConstant) * (static_cast<ParamType>(_binCenters[k]) + static_cast<ParamType>(_binCenters[j])) * (static_cast<ParamType>(_binCenters[k]) + static_cast<ParamType>(_binCenters[j])) * (static_cast<ParamType>(_binCenters[k]) + static_cast<ParamType>(_binCenters[j])) * yCrystal[j] * yCrystal[k] * static_cast<ParamType>(_binSizes[j]) * static_cast<ParamType>(_binSizes[k]) / static_cast<ParamType>(_binSizes[i]);
									}
								}
								else if (static_cast<ParamType>(_aggregationIndex) == 3)
								{
									if (cadet_likely(j != k))
									{
										aggregation_source += aggregation_correction_factor * static_cast<ParamType>(_aggregationRateConstant) * (static_cast<ParamType>(_binCenters[k]) * static_cast<ParamType>(_binCenters[k]) * static_cast<ParamType>(_binCenters[k]) + static_cast<ParamType>(_binCenters[j]) * static_cast<ParamType>(_binCenters[j]) * static_cast<ParamType>(_binCenters[j])) * yCrystal[j] * yCrystal[k] * static_cast<ParamType>(_binSizes[j]) * static_cast<ParamType>(_binSizes[k]) / static_cast<ParamType>(_binSizes[i]);
									}
									else
									{
										aggregation_source += 0.5 * aggregation_correction_factor * static_cast<ParamType>(_aggregationRateConstant) * (static_cast<ParamType>(_binCenters[k]) * static_cast<ParamType>(_binCenters[k]) * static_cast<ParamType>(_binCenters[k]) + static_cast<ParamType>(_binCenters[j]) * static_cast<ParamType>(_binCenters[j]) * static_cast<ParamType>(_binCenters[j])) * yCrystal[j] * yCrystal[k] * static_cast<ParamType>(_binSizes[j]) * static_cast<ParamType>(_binSizes[k]) / static_cast<ParamType>(_binSizes[i]);
									}
								}
								else if (static_cast<ParamType>(_aggregationIndex) == 4)
								{
									if (cadet_likely(j != k))
									{
										aggregation_source += aggregation_correction_factor * static_cast<ParamType>(_aggregationRateConstant) * (static_cast<ParamType>(_binCenters[k]) + static_cast<ParamType>(_binCenters[j])) * (static_cast<ParamType>(_binCenters[k]) + static_cast<ParamType>(_binCenters[j])) * (static_cast<ParamType>(_binCenters[k]) * static_cast<ParamType>(_binCenters[k]) - static_cast<ParamType>(_binCenters[j]) * static_cast<ParamType>(_binCenters[j])) * yCrystal[j] * yCrystal[k] * static_cast<ParamType>(_binSizes[j]) * static_cast<ParamType>(_binSizes[k]) / static_cast<ParamType>(_binSizes[i]);
									}
									else
									{
										aggregation_source += 0.5 * aggregation_correction_factor * static_cast<ParamType>(_aggregationRateConstant) * (static_cast<ParamType>(_binCenters[k]) + static_cast<ParamType>(_binCenters[j])) * (static_cast<ParamType>(_binCenters[k]) + static_cast<ParamType>(_binCenters[j])) * (static_cast<ParamType>(_binCenters[k]) * static_cast<ParamType>(_binCenters[k]) - static_cast<ParamType>(_binCenters[j]) * static_cast<ParamType>(_binCenters[j])) * yCrystal[j] * yCrystal[k] * static_cast<ParamType>(_binSizes[j]) * static_cast<ParamType>(_binSizes[k]) / static_cast<ParamType>(_binSizes[i]);
									}
								}
							}
							else
							{
								// do nothing
							}
						}
					}
					if (cadet_likely(i > 0))
					{
						// Use v_g from previous iteration, this results in the same flux => conservative
						resCrystal[i] += factor / static_cast<ParamType>(_binSizes[i]) * (v_g * yCrystal[i - 1]) + factor * breakage_source + factor * aggregation_source - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i] - yCrystal[i - 1]) / static_cast<ParamType>(_binCenterDists[i - 1]);
					}
					else
					{
						// Left boundary condition, no aggregation source term is present in the first bin
						resCrystal[i] += factor / static_cast<ParamType>(_binSizes[i]) * B_0 + factor * breakage_source + factor * aggregation_source;
					}
					// Flux through right face
					// Calculate breakage sink term
					selection_function = static_cast<ParamType>(_breakageRateConstant) * pow(static_cast<ParamType>(_binCenters[i]), 3.0 * static_cast<ParamType>(_breakageSelectionFunctionAlpha));
					breakage_sink = yCrystal[i] * selection_function * static_cast<ParamType>(Upsilon_sink[i]);
					for (int j = 0; j < _nBins; ++j)
					{
						// constant kernel 0, brownian kernel 1, smoluchowski kernel 2, golovin kernel 3, differential force kernel 4
						if (static_cast<ParamType>(_aggregationIndex) == 0)
						{
							aggregation_sink += static_cast<ParamType>(_aggregationRateConstant) * yCrystal[j] * yCrystal[i] * static_cast<ParamType>(_binSizes[j]);
						}
						else if (static_cast<ParamType>(_aggregationIndex) == 1)
						{
							aggregation_sink += static_cast<ParamType>(_aggregationRateConstant) * (static_cast<ParamType>(_binCenters[i]) + static_cast<ParamType>(_binCenters[j])) * (static_cast<ParamType>(_binCenters[i]) + static_cast<ParamType>(_binCenters[j])) / static_cast<ParamType>(_binCenters[i]) / static_cast<ParamType>(_binCenters[j]) * yCrystal[j] * yCrystal[i] * static_cast<ParamType>(_binSizes[j]);
						}
						else if (static_cast<ParamType>(_aggregationIndex) == 2)
						{
							aggregation_sink += static_cast<ParamType>(_aggregationRateConstant) * (static_cast<ParamType>(_binCenters[i]) + static_cast<ParamType>(_binCenters[j])) * (static_cast<ParamType>(_binCenters[i]) + static_cast<ParamType>(_binCenters[j])) * (static_cast<ParamType>(_binCenters[i]) + static_cast<ParamType>(_binCenters[j])) * yCrystal[j] * yCrystal[i] * static_cast<ParamType>(_binSizes[j]);
						}
						else if (static_cast<ParamType>(_aggregationIndex) == 3)
						{
							aggregation_sink += static_cast<ParamType>(_aggregationRateConstant) * (static_cast<ParamType>(_binCenters[i]) * static_cast<ParamType>(_binCenters[i]) * static_cast<ParamType>(_binCenters[i]) + static_cast<ParamType>(_binCenters[j]) * static_cast<ParamType>(_binCenters[j]) * static_cast<ParamType>(_binCenters[j])) * yCrystal[j] * yCrystal[i] * static_cast<ParamType>(_binSizes[j]);
						}
						else if (static_cast<ParamType>(_aggregationIndex) == 4)
						{
							// absolute value
							if (i > j)
							{
								aggregation_sink += static_cast<ParamType>(_aggregationRateConstant) * (static_cast<ParamType>(_binCenters[i]) + static_cast<ParamType>(_binCenters[j])) * (static_cast<ParamType>(_binCenters[i]) + static_cast<ParamType>(_binCenters[j])) * (static_cast<ParamType>(_binCenters[i]) * static_cast<ParamType>(_binCenters[i]) - static_cast<ParamType>(_binCenters[j]) * static_cast<ParamType>(_binCenters[j])) * yCrystal[j] * yCrystal[i] * static_cast<ParamType>(_binSizes[j]);
							}
							else
							{
								aggregation_sink -= static_cast<ParamType>(_aggregationRateConstant) * (static_cast<ParamType>(_binCenters[i]) + static_cast<ParamType>(_binCenters[j])) * (static_cast<ParamType>(_binCenters[i]) + static_cast<ParamType>(_binCenters[j])) * (static_cast<ParamType>(_binCenters[i]) * static_cast<ParamType>(_binCenters[i]) - static_cast<ParamType>(_binCenters[j]) * static_cast<ParamType>(_binCenters[j])) * yCrystal[j] * yCrystal[i] * static_cast<ParamType>(_binSizes[j]);
							}
						}
					}
					v_g = k_g_times_s_g * pow(1.0 + static_cast<ParamType>(_growthConstant) * static_cast<ParamType>(_bins[i + 1]), static_cast<ParamType>(_p));
					if (cadet_likely(i + 1 < _nBins))
					{
						resCrystal[i] -= factor / static_cast<ParamType>(_binSizes[i]) * (v_g * yCrystal[i]) + factor * breakage_sink + factor * aggregation_sink - factor * static_cast<ParamType>(_growthDispersionRate) / static_cast<ParamType>(_binSizes[i]) * (yCrystal[i + 1] - yCrystal[i]) / static_cast<ParamType>(_binCenterDists[i]);
					}
					else
					{
						// sink term still exists for the last bin
						resCrystal[i] -= factor * breakage_sink + factor * aggregation_sink;
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
