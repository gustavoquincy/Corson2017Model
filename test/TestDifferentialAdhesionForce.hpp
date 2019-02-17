//
// Created by gq on 2/10/19.
//

#ifndef CHASTE_TESTDIFFERENTIALADHESIONFORCE_HPP
#define CHASTE_TESTDIFFERENTIALADHESIONFORCE_HPP

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "CellLabel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "CellAgesWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "FakePetscSetup.hpp"

#include "NagaiHondaDifferentialAdhesionForce.hpp"

#include "SimpleTargetAreaModifier.hpp"

class TestRunningDifferentialAdhesionSimulationsTutorial : public AbstractCellBasedTestSuite
{
public:

    void TestVertexBasedDifferentialAdhesionSimulation()
    {
        HoneycombVertexMeshGenerator generator(20, 20);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        p_mesh->SetCellRearrangementThreshold(0.1);

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(), std::vector<unsigned>(), p_diff_type);

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        cell_population.AddCellPopulationCountWriter<CellProliferativePhasesCountWriter>();
        cell_population.AddCellWriter<CellProliferativePhasesWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();

        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            if (RandomNumberGenerator::Instance()->ranf() < 0.5)
            {
                cell_iter->AddCellProperty(p_label);
            }
        }

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestVertexBasedDifferentialAdhesionSimulation");
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(5.0);

        MAKE_PTR(NagaiHondaDifferentialAdhesionForce<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(55.0);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(0.0);
        p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(1.0);
        p_force->SetNagaiHondaLabelledCellCellAdhesionEnergyParameter(6.0);
        p_force->SetNagaiHondaLabelledCellLabelledCellAdhesionEnergyParameter(3.0);
        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(12.0);
        p_force->SetNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter(40.0);
        simulator.AddForce(p_force);

        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        simulator.Solve();
    }

};

#endif //CHASTE_TESTDIFFERENTIALADHESIONFORCE_HPP
