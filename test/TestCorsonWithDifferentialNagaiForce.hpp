//
// Created by gq on 1/14/19.
//

#ifndef TESTCORSONWITHDIFFERENTIALNAGAIFORCE_HPP_
#define TESTCORSONWITHDIFFERENTIALNAGAIFORCE_HPP_

#include <cxxtest/TestSuite.h>
#include <properties/proliferative_types/DifferentiatedCellProliferativeType.hpp>
#include "CheckpointArchiveTypes.hpp"
//#include "AbstractCellBasedTestSuite.hpp"

#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "CellAgesWriter.hpp"
#include "CellCorsonWriter.hpp"
#include "CellIdWriter.hpp"
#include "CellLabel.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellsGenerator.hpp"
#include "CorsonSrnModel.hpp"
#include "CorsonTrackingModifier.hpp"

//#include "Exception.hpp"

#include "CylindricalHoneycombVertexMeshGenerator.hpp"

#include "NagaiHondaDifferentialAdhesionForce.hpp"
#include "WelikyOsterForce.hpp"
#include "FarhadifarForce.hpp"

#include "RandomMotionForce.hpp"

#include "OffLatticeSimulation.hpp"

#include "SmartPointers.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "WildTypeCellMutationState.hpp"

#include "DeltaNotchSrnModel.hpp"
/*
 * The next header defines the simulation class modifier corresponding to the Delta-Notch SRN model.
 * This modifier leads to the {{{CellData}}} cell property being updated at each timestep to deal with Delta-Notch signalling.
 */
#include "DeltaNotchTrackingModifier.hpp"


#include "CellBasedEventHandler.hpp"

#include "PetscSetupAndFinalize.hpp"


static const double M_TISSUE_RADIUS = 10; //18
static const double M_TIME_FOR_SIMULATION = 5; //5


class TestCorsonWithDifferentialNagaiForce : public AbstractCellBasedWithTimingsTestSuite
{

public:
    void TestCorsonWithDifferentialNagaiHondaForce()
    {
        SimulationTime::Instance()->Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);

        HoneycombVertexMeshGenerator generator(M_TISSUE_RADIUS, M_TISSUE_RADIUS);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        p_mesh->SetCellRearrangementThreshold(0.1);

        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellProperty> p_cell_type(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());
        CellsGenerator<UniformG1GenerationalCellCycleModel,2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_cell_type);

        //vector for setting cell_state_u to 0 in the following loop

        for (unsigned i=0; i<cells.size(); i++)
        {
            std::vector<double> initial_conditions;
            initial_conditions.push_back(0);
            initial_conditions.push_back(0);

            CorsonSrnModel* p_srn_model = new CorsonSrnModel();
            p_srn_model->SetInitialConditions(initial_conditions);

            cells[i]->SetSrnModel(p_srn_model);

            double birth_time = 0.0;
            cells[i]->SetBirthTime(birth_time);

            cells[i]->GetCellData()->SetItem("target area", 1.0);
        }
        //GenerateCells(p_mesh->GetNumElements(), cells);

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellCorsonWriter>();

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestCorsonOdeWithFarhadifarForce");

        simulator.SetDt(1.0/200.0); //1/200
        simulator.SetSamplingTimestepMultiple(10); //200, 50
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        MAKE_PTR(CorsonTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        /*p_force->SetNagaiHondaDeformationEnergyParameter(50.0);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(1.0);
        p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(1.0);
        p_force->SetNagaiHondaLabelledCellCellAdhesionEnergyParameter(2.0);
        p_force->SetNagaiHondaLabelledCellLabelledCellAdhesionEnergyParameter(1.0);
        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(10.0);
        p_force->SetNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter(20.0);*/
        MAKE_PTR(FarhadifarForce<2>, p_force);
        /*p_force->SetLineTensionParameter(0.12); //0.12
        p_force->SetPerimeterContractilityParameter(0.04); //0.04
        p_force->OutputForceParameters();*/
        simulator.AddForce(p_force);

        // Add some noise (to avoid local minimum)
        MAKE_PTR(RandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(0.005); //0.01 is the initial value
        simulator.AddForce(p_random_force);


        //CellBasedEventHandler::Reset();
        simulator.Solve();
        //CellBasedEventHandler::Headings();
        //CellBasedEventHandler::Report();

        /*//Iterate cell population to label those whose u is larger than 0.5
        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            if (cell_iter->GetCellData()->GetItem("cellstate u") > 0.3)
            {
                cell_iter->AddCellProperty(p_label);
            }
        }

        simulator.SetEndTime(M_TIME_TO_BIFURCATION + M_TIME_FOR_SIMULATION);
        simulator.Solve();*/


    }
};


#endif //TESTCORSONWITHDIFFERENTIALNAGAIFORCE_HPP_
