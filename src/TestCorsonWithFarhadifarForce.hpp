//
// Created by gq on 3/26/19.
//

#ifndef CHASTE_TESTCORSONWITHFARHADIFARFORCE_HPP
#define CHASTE_TESTCORSONWITHFARHADIFARFORCE_HPP

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

#include "CellBasedEventHandler.hpp"

#include "PetscSetupAndFinalize.hpp"


static const double M_TISSUE_RADIUS = 10; //18
static const double M_TIME_FOR_SIMULATION = 5;


class TestCorsonWithDifferentialNagaiForce : public AbstractCellBasedWithTimingsTestSuite
{
private:
    void GenerateCells(unsigned num_cells, std::vector<CellPtr>& rCells)
    {
        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_prolif_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());

        for (unsigned i = 0; i < num_cells; ++i)
        {
            std::vector<double> initial_conditions;
            //RandomNumberGenerator::Instance()->ranf()
            initial_conditions.push_back(0.01);

            UniformG1GenerationalCellCycleModel* p_cc_model = new UniformG1GenerationalCellCycleModel();
            p_cc_model->SetDimension(2);

            CorsonSrnModel* p_srn_model = new CorsonSrnModel();
            //p_srn_model->Initialise();
            p_srn_model->SetInitialConditions(initial_conditions);

            CellPtr p_cell(new Cell(p_state, p_cc_model, p_srn_model));
            p_cell->SetCellProliferativeType(p_prolif_type);

            double birth_time = 0.0;
            p_cell->SetBirthTime(birth_time);

            //p_cell->InitialiseCellCycleModel();
            p_cell->GetCellData()->SetItem("target area", 1.0);
            rCells.push_back(p_cell);


        }
    }
public:
    void TestCorsonWithFarhadifarForce()
    {
        SimulationTime::Instance()->Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);

        HoneycombVertexMeshGenerator generator(M_TISSUE_RADIUS, M_TISSUE_RADIUS);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        p_mesh->SetCellRearrangementThreshold(0.1);

        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumElements(), cells);

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellCorsonWriter>();

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestCorsonWifthFarhadifarForce");

        simulator.SetDt(1.0/200.0); //1/200
        simulator.SetSamplingTimestepMultiple(10); //200, 50
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        MAKE_PTR(CorsonTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);


        MAKE_PTR(FarhadifarForce<2>, p_force);
        /*p_force->SetLineTensionParameter(-0.85);
        p_force->SetBoundaryLineTensionParameter(-0.85);
        p_force->SetPerimeterContractilityParameter(0.1);*/
        simulator.AddForce(p_force);

        // Add some noise (to avoid local minimum)
        /*MAKE_PTR(RandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(0.005);
        simulator.AddForce(p_random_force);*/


        //CellBasedEventHandler::Reset();
        simulator.Solve();
        //CellBasedEventHandler::Headings();
        //CellBasedEventHandler::Report();




    }
};

#endif //CHASTE_TESTCORSONWITHFARHADIFARFORCE_HPP
