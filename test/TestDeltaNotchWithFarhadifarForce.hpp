//
// Created by gq on 4/12/19.
//

#ifndef CHASTE_TESTDELTANOTCHWITHFARHADIFARFORCE_HPP
#define CHASTE_TESTDELTANOTCHWITHFARHADIFARFORCE_HPP


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
#include "DeltaNotchSrnModel.hpp"
#include "DeltaNotchTrackingModifier.hpp"
#include "CellDeltaNotchWriter.hpp"

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

#include "PlaneBoundaryCondition.hpp"
#include "SphereGeometryBoundaryCondition.hpp"

#include "Debug.hpp"
#include "CellBasedEventHandler.hpp"

#include "PetscSetupAndFinalize.hpp"


static const double M_TISSUE_RADIUS = 10; //18
static const double M_TIME_FOR_SIMULATION = 20; //5


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
            //RandomNumberGenerator::Instance()->ranf()
            initial_conditions.push_back(0.01);
            initial_conditions.push_back(0.01);

            UniformG1GenerationalCellCycleModel* p_cc_model = new UniformG1GenerationalCellCycleModel();
            p_cc_model->SetDimension(2);

            DeltaNotchSrnModel* p_srn_model = new DeltaNotchSrnModel();
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
    void TestDeltaNotchWithFarhadifarForce()
    {
        SimulationTime::Instance()->Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);

        /*CylindricalHoneycombVertexMeshGenerator generator(M_TISSUE_RADIUS, M_TISSUE_RADIUS);    // Parameters are: cells across, cells up
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();*/

        HoneycombVertexMeshGenerator generator(M_TISSUE_RADIUS, M_TISSUE_RADIUS);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        p_mesh->SetCellRearrangementThreshold(0.1);

        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumElements(), cells);

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellDeltaNotchWriter>();

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestDeltaNotchWithFarhadifarForce0415");

        simulator.SetDt(1.0/50); //1/200, 1/100
        simulator.SetSamplingTimestepMultiple(10); //200, 50
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        MAKE_PTR(DeltaNotchTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        MAKE_PTR(FarhadifarForce<2>, p_force);
        p_force->SetAreaElasticityParameter(1);
        p_force->SetBoundaryLineTensionParameter(0.12);
        p_force->SetLineTensionParameter(0.12);
        p_force->SetPerimeterContractilityParameter(0);
        simulator.AddForce(p_force);

        //And box boundary condition (0,0) -> (0,10) -> (10,10) -> (10,0) -> (0,0)
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
        normal(1) = 1.0;
        normal(0) = 0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc_up, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc_up);

        point(0) = M_TISSUE_RADIUS+1;
        normal(0) = 1;
        normal(1) = 0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc_right, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc_right);


        // Add some noise (to avoid local minimum)
        /*MAKE_PTR(RandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(0.001);
        simulator.AddForce(p_random_force);*/

        /*c_vector<double,2> point = zero_vector<double>(2);
        c_vector<double,2> normal = zero_vector<double>(2);
        normal(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc);*/
        /*c_vector<double,2> center = zero_vector<double>(2);
        double radius = 5;
        MAKE_PTR_ARGS(SphereGeometryBoundaryCondition<2>, p_sbc, (&cell_population, center, radius));
        simulator.AddCellPopulationBoundaryCondition(p_sbc);*/
        //CellBasedEventHandler::Reset();
        simulator.Solve();
        //CellBasedEventHandler::Headings();
        //CellBasedEventHandler::Report();




    }
};

#endif //CHASTE_TESTDELTANOTCHWITHFARHADIFARFORCE_HPP
