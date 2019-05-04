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
#include "CellIdWriter.hpp"
#include "CellLabel.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellsGenerator.hpp"

#include "CorsonSrnModel.hpp"
#include "CorsonTrackingModifier.hpp"
#include "CellCorsonWriter.hpp"
//#include "Exception.hpp"

#include "Toroidal2dVertexMesh.hpp"
#include "ToroidalHoneycombVertexMeshGenerator.hpp"

//Force Collection
#include "FarhadifarForce.hpp"
#include "NagaiHondaForce.hpp"
#include "WelikyOsterForce.hpp"

#include "RandomMotionForce.hpp"

#include "OffLatticeSimulation.hpp"

#include "SmartPointers.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "WildTypeCellMutationState.hpp"

#include "PlaneBoundaryCondition.hpp"

#include "Debug.hpp"
#include "CellBasedEventHandler.hpp"

#include "PetscSetupAndFinalize.hpp"


static const double M_TISSUE_RADIUS = 14; //18 10

static const double M_TIME_FOR_SIMULATION = 20; //5 10 20


class CorsonWithFarhadifarF : public AbstractCellBasedWithTimingsTestSuite
{
private:
    void GenerateCells(unsigned num_cells, std::vector<CellPtr>& rCells)
    {
        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_prolif_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());

        for (unsigned i = 0; i < num_cells; ++i)
        {
            std::vector<double> initial_conditions;
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
            p_cell->GetCellData()->SetItem("target area", 0.5*sqrt(3)); //1
            rCells.push_back(p_cell);

        }
    }
public:
    void TestCorsonWithFarhadifarForce()
    {
        SimulationTime::Instance()->Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);

        ToroidalHoneycombVertexMeshGenerator generator(M_TISSUE_RADIUS, M_TISSUE_RADIUS);
        Toroidal2dVertexMesh* p_mesh = generator.GetToroidalMesh();
        //p_mesh->SetCellRearrangementThreshold(0.5);

        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumElements(), cells);

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellCorsonWriter>();

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestCorsonWifthFarhadifarForce0503");

        /*
         *1/200
         *1/100
         * 1/50
         * 5.0/100
         * */
        simulator.SetDt(5.0/100);
        /*
         *
         * 200, 50, 10*/
        simulator.SetSamplingTimestepMultiple(10);

        simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        /*MAKE_PTR(WelikyOsterForce<2>, p_force);
        simulator.AddForce(p_force);*/
        MAKE_PTR(FarhadifarForce<2>, p_force);
        p_force->SetAreaElasticityParameter(1);
        p_force->SetBoundaryLineTensionParameter(0.12);
        p_force->SetLineTensionParameter(0.12);
        p_force->SetPerimeterContractilityParameter(0);
        simulator.AddForce(p_force);

    /*    //And box boundary condition (0,0) -> (0,10) -> (10,10) -> (10,0) -> (0,0)
        c_vector<double,2> point = zero_vector<double>(2);
        c_vector<double,2> normal = zero_vector<double>(2);
        normal(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc_down, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc_down);

        point(0) = -1;
        normal(0) = -1;
        normal(1) = 0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc_left, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc_left);

        point(1) = M_TISSUE_RADIUS;
        normal(1) = 1;
        normal(0) = 0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc_up, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc_up);

        point(0) = M_TISSUE_RADIUS+1;
        normal(0) = 1;
        normal(1) = 0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc_right, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc_right);*/


        // Add some noise (to avoid local minimum)
        MAKE_PTR(RandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(0.01); //0.001
        simulator.AddForce(p_random_force);

        MAKE_PTR(CorsonTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        CellBasedEventHandler::Reset();
        simulator.Solve();
        CellBasedEventHandler::Headings();
        CellBasedEventHandler::Report();

        /*out_stream& ParamsFile;
        ParamsFile << "\t\t\t<Simulation Time>" << M_TISSUE_RADIUS << "</Tissue Radius>\n";
*//*
        MAKE_PTR(CorsonTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        simulator.SetEndTime(M_TIME_FOR_SIMULATION+10);
        simulator.Solve();*/

    }
};

#endif //CHASTE_TESTCORSONWITHFARHADIFARFORCE_HPP
