// =============================================================================
//  CADET
//  
//  Copyright © 2008-present: The CADET-Core Authors
//            Please see the AUTHORS.md file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include <catch.hpp>

#include "ColumnTests.hpp"
#include "ParticleHelper.hpp"
#include "ReactionModelTests.hpp"
#include "Utils.hpp"
#include "JsonTestModels.hpp"


TEST_CASE("LRMP2D inlet DOF Jacobian", "[LRMP2D],[DG],[DG2D],[UnitOp],[Jacobian],[Inlet],[CILRMP2D]")
{
	cadet::test::column::testInletDofJacobian("LUMPED_RATE_MODEL_WITH_PORES_2D", "DG");
}

TEST_CASE("LRMP2D time derivative Jacobian vs FD", "[LRMP2D],[DG],[DG2D],[UnitOp],[Residual],[Jacobian],[CILRMP2D],[FD]")
{
	cadet::test::column::testTimeDerivativeJacobianFD("LUMPED_RATE_MODEL_WITH_PORES_2D", "DG", 1e-6, 0.0, 9e-4);
}

TEST_CASE("LRMP2D transport Jacobian", "[LRMP2D],[DG],[DG2D],[UnitOp],[Jacobian],[CILRMP2D]")
{
	const std::string relModelFilePath = std::string("/data/lrmp2d_bulkTransport_1comp_debug.json");
	cadet::JsonParameterProvider jpp = cadet::test::column::getReferenceFile(relModelFilePath);

	// get the number of radial ports
	jpp.pushScope("model");
	const int nUnits = jpp.getInt("NUNITS"); // there is one column and (nUnits-1) inlets, one per radial port
	const int columnIdx = 0;
	const std::string unitID = "000";
	const int nRad = nUnits - 1; // number of radial points or ports

	// we need to set flowRates for 2D models for the JacobianAD test, since velocity is only set with a call to the setFlowRate function,
	// which in turn is only called if we specify flow rates. We get the flow rates from the connections matrix
	jpp.pushScope("connections");
	jpp.pushScope("switch_000");

	const std::vector<double> connections = jpp.getDoubleArray("CONNECTIONS");
	std::vector<cadet::active> flowRate;
	for (int i = 1; i <= nRad; i++)
		flowRate.push_back(connections[i * 7 - 1]);

	jpp.popScope();
	jpp.popScope();
	jpp.pushScope("unit_" + unitID);

	// This test might run out of memory due to the required AD directions:
	// (axPolyDeg + 1) * axNElem * (radPolyDeg + 1) * radNElem * (nComp + nParType * (nComp + nBound))
	// result here is 588 pure dofs (21axPoints*14radPoints) * (1 + 1)
	for (int zElem = 1; zElem < 8; zElem++)
	{
		for (int rElem = 1; rElem < 8; rElem++)
		{
			jpp.pushScope("discretization");
			jpp.set("AX_NELEM", zElem);
			jpp.set("RAD_NELEM", rElem);

			jpp.popScope();

			cadet::test::column::testJacobianAD(jpp, 1e10, std::numeric_limits<float>::epsilon() * 100.0, &flowRate[0]); // @todo figure out why FD Jacobian pattern comparison doesnt work but AD Jacobian comparison does
		}
	}
}

TEST_CASE("LRMP2D with two component linear binding Jacobian", "[LRMP2D],[DG],[DG2D],[UnitOp],[Jacobian],[CILRMP2D]")
{
	const std::string relModelFilePath = std::string("/data/lrmp2d_dynLin_2comp_debug.json");
	cadet::JsonParameterProvider jpp = cadet::test::column::getReferenceFile(relModelFilePath);

	// get the number of radial ports
	jpp.pushScope("model");
	const int nUnits = jpp.getInt("NUNITS"); // there is one column and (nUnits-1) inlets, one per radial port
	const int columnIdx = 0;
	const std::string unitID = "000";
	const int nRad = nUnits - 1; // number of radial points or ports

	// we need to set flowRates for 2D models for the JacobianAD test, since velocity is only set with a call to the setFlowRate function,
	// which in turn is only called if we specify flow rates. We get the flow rates from the connections matrix
	jpp.pushScope("connections");
	jpp.pushScope("switch_000");

	const std::vector<double> connections = jpp.getDoubleArray("CONNECTIONS");
	std::vector<cadet::active> flowRate;
	for (int i = 1; i <= nRad; i++)
		flowRate.push_back(connections[i * 7 - 1]);

	jpp.popScope();
	jpp.popScope();
	jpp.pushScope("unit_" + unitID);

	// This test might run out of memory due to the required AD directions:
	// (axPolyDeg + 1) * axNElem * (radPolyDeg + 1) * radNElem * (nComp + nParType * (nComp + nBound))
	// result here is 1296 pure dofs (18axPoints*12radPoints) * (2 + 4)
	for (int zElem = 1; zElem < 7; zElem++)
	{
		for (int rElem = 1; rElem < 7; rElem++)
		{
			jpp.pushScope("discretization");
			jpp.set("AX_NELEM", zElem);
			jpp.set("RAD_NELEM", rElem);

			jpp.popScope();

			cadet::test::column::testJacobianAD(jpp, 1e10, std::numeric_limits<float>::epsilon() * 100.0, &flowRate[0]); // @todo figure out why FD Jacobian pattern comparison doesnt work but AD Jacobian comparison does
		}
	}
}

TEST_CASE("LRMP2D sensitivity Jacobians", "[LRMP2D],[UnitOp],[Sensitivity],[CILRMP2D]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("LUMPED_RATE_MODEL_WITH_PORES_2D", "DG");

	cadet::test::column::testFwdSensJacobians(jpp, 1e-4, 6e-7);
}

TEST_CASE("LRMP2D consistent initialization with linear binding", "[LRMP2D],[ConsistentInit],[CILRMP2D]")
{
	cadet::test::column::testConsistentInitializationLinearBinding("LUMPED_RATE_MODEL_WITH_PORES_2D", "DG", 1e-12, 1e-12, 0, 0);
	cadet::test::column::testConsistentInitializationLinearBinding("LUMPED_RATE_MODEL_WITH_PORES_2D", "DG", 1e-12, 1e-12, 1, 0);
	cadet::test::column::testConsistentInitializationLinearBinding("LUMPED_RATE_MODEL_WITH_PORES_2D", "DG", 1e-12, 1e-12, 0, 1);
	//cadet::test::column::testConsistentInitializationLinearBinding("LUMPED_RATE_MODEL_WITH_PORES_2D", "DG", 1e-12, 1e-12, 1, 1); // @todo AD with req binding does not work
}

//TEST_CASE("LRMP2D consistent initialization with SMA binding", "[LRMP2D],[ConsistentInit],[todo]")  // todo fix (also doesnt work for other models)
//{
//	std::vector<double> y(4 * 3 + 4 * 8 * 3 + 8 * 3 * 3 * (4 + 4) + 4 * 8 * 3, 0.0);
//	// Optimal values:
//	//	const double bindingCell[] = {1.2, 2.0, 1.0, 1.5, 858.034, 66.7896, 3.53273, 2.53153, 
//	//		1.0, 1.8, 1.5, 1.6, 856.173, 64.457, 5.73227, 2.85286};
//	const double bindingCell[] = { 1.2, 2.0, 1.0, 1.5, 840.0, 63.0, 3.0, 3.0,
//		1.0, 1.8, 1.5, 1.6, 840.0, 63.0, 6.0, 3.0 };
//	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4 * 3 + 4 * 8 * 3);
//	cadet::test::util::repeat(y.data() + 4 * 3 + 4 * 8 * 3, bindingCell, 16, 3 * 8 * 3 / 2);
//	cadet::test::util::populate(y.data() + 4 * 3 + 4 * 8 * 3 + 8 * 3 * 3 * (4 + 4), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4 * 8 * 3);
//
//	cadet::test::column::testConsistentInitializationSMABinding("LUMPED_RATE_MODEL_WITH_PORES_2D", "DG", y.data(), 1e-14, 1e-5, 0, 0);
//	cadet::test::column::testConsistentInitializationSMABinding("LUMPED_RATE_MODEL_WITH_PORES_2D", "DG", y.data(), 1e-14, 1e-5, 1, 0);
//	//cadet::test::column::testConsistentInitializationSMABinding("LUMPED_RATE_MODEL_WITH_PORES_2D", "DG", y.data(), 1e-14, 1e-5, 0, 1);
//	//cadet::test::column::testConsistentInitializationSMABinding("LUMPED_RATE_MODEL_WITH_PORES_2D", "DG", y.data(), 1e-14, 1e-5, 1, 1);
//}

TEST_CASE("LRMP2D consistent sensitivity initialization with linear binding", "[LRMP2D],[ConsistentInit],[Sensitivity],[CILRMP2D]")
{
	// Fill state vector with given initial values
	const unsigned int numDofs = 2 * 3 + 2 * 8 * 3 + 8 * 3 * 3 * (2 + 2) + 2 * 8 * 3;
	std::vector<double> y(numDofs, 0.0);
	std::vector<double> yDot(numDofs, 0.0);
	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, numDofs);
	cadet::test::util::populate(yDot.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.9)) + 1e-4; }, numDofs);

	// todo: kinetic binding doesnt work here
	//cadet::test::column::testConsistentInitializationSensitivity("LUMPED_RATE_MODEL_WITH_PORES_2D", "DG", y.data(), yDot.data(), true, 1e-14, 0, 0);
	cadet::test::column::testConsistentInitializationSensitivity("LUMPED_RATE_MODEL_WITH_PORES_2D", "DG", y.data(), yDot.data(), true, 1e-14, 1, 0);
	//cadet::test::column::testConsistentInitializationSensitivity("LUMPED_RATE_MODEL_WITH_PORES_2D", "DG", y.data(), yDot.data(), true, 1e-14, 0, 1);
	cadet::test::column::testConsistentInitializationSensitivity("LUMPED_RATE_MODEL_WITH_PORES_2D", "DG", y.data(), yDot.data(), true, 1e-14, 1, 1);
}

TEST_CASE("LRMP2D consistent sensitivity initialization with SMA binding", "[LRMP2D],[ConsistentInit],[Sensitivity],[CILRMP2D]")
{
	// Fill state vector with given initial values
	const unsigned int numDofs = 4 * 3 + 4 * 8 * 3 + 8 * 3 * 3 * (4 + 4) + 4 * 8 * 3;
	std::vector<double> y(numDofs, 0.0);
	std::vector<double> yDot(numDofs, 0.0);

	const double bindingCell[] = { 1.0, 1.8, 1.5, 1.6, 840.0, 63.0, 6.0, 3.0 };
	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4 * 3 + 4 * 8 * 3);
	cadet::test::util::repeat(y.data() + 4 * 3 + 4 * 8 * 3, bindingCell, 8, 3 * 8 * 3);
	cadet::test::util::populate(y.data() + 4 * 3 + 4 * 8 * 3 + 8 * 3 * 3 * (4 + 4), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4 * 8 * 3);

	cadet::test::util::populate(yDot.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.9)) + 1e-4; }, numDofs);

	// todo: kinetic binding doesnt work here
	//cadet::test::column::testConsistentInitializationSensitivity("LUMPED_RATE_MODEL_WITH_PORES_2D", "DG", y.data(), yDot.data(), false, 1e-9, 0, 0);
	cadet::test::column::testConsistentInitializationSensitivity("LUMPED_RATE_MODEL_WITH_PORES_2D", "DG", y.data(), yDot.data(), false, 1e-9, 1, 0);
	//cadet::test::column::testConsistentInitializationSensitivity("LUMPED_RATE_MODEL_WITH_PORES_2D", "DG", y.data(), yDot.data(), false, 1e-9, 0, 1);
	cadet::test::column::testConsistentInitializationSensitivity("LUMPED_RATE_MODEL_WITH_PORES_2D", "DG", y.data(), yDot.data(), false, 1e-9, 1, 1);
}

//TEST_CASE("LRMP2D LWE one vs two identical particle types match", "[LRMP2D],[Simulation],[ParticleType],[todo]") // todo
//{
//	cadet::test::particle::testOneVsTwoIdenticalParticleTypes("LUMPED_RATE_MODEL_WITH_PORES_2D", "DG", 1e-7, 5e-5);
//}
//
//TEST_CASE("LRMP2D LWE separate identical particle types match", "[LRMP2D],[Simulation],[ParticleType],[todo]") // todo
//{
//	cadet::test::particle::testSeparateIdenticalParticleTypes("LUMPED_RATE_MODEL_WITH_PORES_2D", "DG", 1e-15, 1e-15);
//}

TEST_CASE("LRMP2D linear binding single particle matches particle distribution", "[LRMP2D],[Simulation],[ParticleType],[todo]") // todo
{
	cadet::test::particle::testLinearMixedParticleTypes("LUMPED_RATE_MODEL_WITH_PORES_2D", "DG", 5e-8, 5e-5);
}

TEST_CASE("LRMP2D multiple particle types Jacobian analytic vs AD", "[LRMP2D],[Jacobian],[AD],[ParticleType],[todo]") // todo
{
	cadet::test::particle::testJacobianMixedParticleTypes("LUMPED_RATE_MODEL_WITH_PORES_2D", "DG");
}

TEST_CASE("LRMP2D multiple particle types time derivative Jacobian vs FD", "[LRMP2D],[UnitOp],[Residual],[Jacobian],[ParticleType],[CILRMP2D]")
{
	cadet::test::particle::testTimeDerivativeJacobianMixedParticleTypesFD("LUMPED_RATE_MODEL_WITH_PORES_2D", "DG", 1e-6, 0.0, 5e-3);
}

TEST_CASE("LRMP2D linear binding single particle matches spatially dependent particle distribution", "[LRMP2D],[Simulation],[ParticleType],[todo]") // todo
{
	cadet::test::particle::testLinearSpatiallyMixedParticleTypes("LUMPED_RATE_MODEL_WITH_PORES_2D", "DG", 5e-8, 5e-5);
}

TEST_CASE("LRMP2D dynamic reactions time derivative Jacobian vs FD bulk", "[LRMP2D],[Jacobian],[Residual],[ReactionModel],[CILRMP2D]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("LUMPED_RATE_MODEL_WITH_PORES_2D", "DG", true, false, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("LRMP2D dynamic reactions time derivative Jacobian vs FD particle", "[LRMP2D],[Jacobian],[Residual],[ReactionModel],[CILRMP2D]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("LUMPED_RATE_MODEL_WITH_PORES_2D", "DG", false, true, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("LRMP2D dynamic reactions time derivative Jacobian vs FD modified particle", "[LRMP2D],[Jacobian],[Residual],[ReactionModel],[CILRMP2D]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("LUMPED_RATE_MODEL_WITH_PORES_2D", "DG", false, true, true, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("LRMP2D dynamic reactions time derivative Jacobian vs FD bulk and particle", "[LRMP2D],[Jacobian],[Residual],[ReactionModel],[CILRMP2D]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("LUMPED_RATE_MODEL_WITH_PORES_2D", "DG", true, true, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("LRMP2D dynamic reactions time derivative Jacobian vs FD bulk and modified particle", "[LRMP2D],[Jacobian],[Residual],[ReactionModel],[CILRMP2D]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("LUMPED_RATE_MODEL_WITH_PORES_2D", "DG", true, true, true, 1e-6, 1e-14, 8e-4);
}

inline cadet::JsonParameterProvider createColumnWithTwoCompLinearBindingThreeParticleTypes()
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("LUMPED_RATE_MODEL_WITH_PORES_2D", "DG");

	const double parVolFrac[] = { 0.3, 0.6, 0.1 };
	const double parFactor[] = { 0.9, 0.8 };
	cadet::test::particle::extendModelToManyParticleTypes(jpp, 3, parFactor, parVolFrac);

	return jpp;
}

TEST_CASE("LRMP2D multi particle types dynamic reactions time derivative Jacobian vs FD bulk", "[LRMP2D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CILRMP2D]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, false, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("LRMP2D multi particle types dynamic reactions time derivative Jacobian vs FD particle", "[LRMP2D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CILRMP2D]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, false, true, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("LRMP2D multi particle types dynamic reactions time derivative Jacobian vs FD modified particle", "[LRMP2D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CILRMP2D]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, false, true, true, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("LRMP2D multi particle types dynamic reactions time derivative Jacobian vs FD bulk and particle", "[LRMP2D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CILRMP2D]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, true, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("LRMP2D multi particle types dynamic reactions time derivative Jacobian vs FD bulk and modified particle", "[LRMP2D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CILRMP2D]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, true, true, 1e-6, 1e-14, 8e-4);
}