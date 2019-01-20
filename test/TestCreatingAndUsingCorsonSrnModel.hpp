//
// Created by gq on 1/19/19.
//

#ifndef CHASTE_TESTCREATINGANDUSINGCORSONSRNMODEL_HPP
#define CHASTE_TESTCREATINGANDUSINGCORSONSRNMODEL_HPP
/*


/*
 * = An example showing how to create a new subcellular reaction network (SRN) model and use it in a cell-based simulation. =
 *
 * == Introduction ==
 *
 * In the previous cell-based Chaste tutorials, we used existing cell-cycle and SRN models to define how cells
 * proliferate and update and subcellular model. In this tutorial, we show how to create a new SRN model class, and how this
 * can be used in a cell-based simulation.
 *
 * == Including header files ==
 *
 * We begin by including the necessary header files. */
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

/* The next header includes the Boost shared_ptr smart pointer, and defines some useful
 * macros to save typing when using it. */
#include "SmartPointers.hpp"
/* The next header includes the NEVER_REACHED macro, used in one of the methods below. */
#include "Exception.hpp"

/* The next header defines a base class for ode-based SRN models.
 * Our new SRN model will inherit from this abstract class. */
#include "AbstractOdeSrnModel.hpp"

/* These headers specify the methods to solve the ODE system. */
#include "AbstractOdeSystem.hpp"
#include "OdeSystemInformation.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"

/* This header specifies the ODE solvers. */
#include "CellCycleModelOdeSolver.hpp"

/* The following headers are needed for checkpointing. */
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/shared_ptr.hpp>

/* The remaining header files define classes that will be used in the cell-based
 * simulation test. We have encountered each of these header files in previous cell-based Chaste
 * tutorials. */
#include "CheckReadyToDivideAndPhaseIsUpdated.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "NagaiHondaForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "OffLatticeSimulation.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"


//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

/*
 * == Defining the SRN model and ODE system classes ==
 *
 * As an example, let us consider a SRN model in which we solve a simple ODE
 * dx/dt = -0.25*y
 * dy/dt = x
 * This has exact solution x = A cos 0.5t + B sin 0.5t
 * where A and B are determined by the initial condions.
 *
 * To implement this model we define a new SRN model, {{{MySrnModel}}},
 * which inherits from {{{AbstractOdeSrnModel}}} and
 * contains a {{{MyOdeSystem}}}.
 *
 * Note that usually this code would be separated out into a separate declaration in
 * a .hpp file and definition in a .cpp file.
 */
class MyOdeSystem : public AbstractOdeSystem
{
private:
    friend class boost::serialization::access;
    /* We only need to include the next block of code if we wish to be able
     * to archive (save or load) the ODE system (and therefore the SRN model) object in a cell-based simulation.
     * The code consists of a serialize method, in which we archive the ODE system
     * using the serialization code defined in the base class
     * {{{AbstractOdeSystem}}}.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractOdeSystem>(*this);
    }

    const double mtau = 1/2;


public:
    MyOdeSystem(std::vector<double> stateVariables=std::vector<double>()) : AbstractOdeSystem(1)
    {
        mpSystemInfo = OdeSystemInformation<MyOdeSystem>::Instance();

        if (stateVariables != std::vector<double>())
        {
            SetStateVariables(stateVariables);
        }
    }

    double SigmoidalFunction(double x) const
    {
        return (1 + tanh(2*x))/2;
    }

    void EvaluateYDerivatives(double time, const std::vector<double>& rY,
                              std::vector<double>& rDY)
    {
        double width = rCellPopulation.GetWidth(0);
        double x = rCellPopulation.GetLocationOfCellCentre(*)
        SimulationTime *p_simulation_time = SimulationTime::Instance();
        double time = p_simulation_time->GetTime();
        double L = mL1 * width;
        double l = ml1 * width;
        //mS0, mtau_g
        double s = mS0* SigmoidalFunction(1-time/mtau_g)*(exp(-pow(x,2)/2/pow(L,2))+exp(-pow(width-x,2)/2/pow(L,2))) +
                   SigmoidalFunction(time/mtau_g-1)*(exp(-pow(x,2)/2/pow(l,2))+exp(-pow(width-x,2)/2/pow(l,2)));

        rDY[0] = (SigmoidalFunction(2*(rY[0]-s)) - rY[0])/mtau;
    }
};

/* As in the ODE tutorials we need to define the ODE system informtion.
 */
template<>
void OdeSystemInformation<MyOdeSystem>::Initialise()
{
    this->mVariableNames.push_back("u");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(1.0);
/*
    this->mVariableNames.push_back("y");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.0);
*/
    this->mInitialised = true;
}


class MySrnModel : public AbstractOdeSrnModel
{
private:

    /* We only need to include the next block of code if we wish to be able
     * to archive (save or load) the SRN model object in a cell-based simulation.
     * The code consists of a serialize method, in which we archive the SRN
     * model using the serialization code defined in the base class
     * {{{AbstractOdeSrnModel}}}.
     */
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractOdeSrnModel>(*this);
    }

protected:
    /**
     * We need to define a protected copy-constructor for use by CreateSrnModel.
     * The only way for external code to create a copy of a SRN model
     * is by calling that method, to ensure that a model of the correct subclass is created.
     * This copy-constructor helps subclasses to ensure that all member variables are correctly copied when this happens.
     *
     * Note that the parent SRN model will have had ResetForDivision() called just before CreateSrnModel() is called,
     * so performing an exact copy of the parent is suitable behaviour. Any daughter-cell-specific initialisation
     * can be done in InitialiseDaughterCell().
     */
    MySrnModel(const MySrnModel& rModel)
            : AbstractOdeSrnModel(rModel)
    {
        /*
         * These lines copy the ODE system.
         */
        assert(rModel.GetOdeSystem());
        SetOdeSystem(new MyOdeSystem(rModel.GetOdeSystem()->rGetStateVariables()));
    }


    /* The first public method is a constructor, which just calls the base
     * constructor.  Note you can include an optional argument to specify the ODE solver.*/
public:

    MySrnModel()
            : AbstractOdeSrnModel(2, boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
    {

        mpOdeSolver = CellCycleModelOdeSolver<MySrnModel, RungeKutta4IvpOdeSolver>::Instance();
        mpOdeSolver->Initialise();
        SetDt(0.1);

        assert(mpOdeSolver->IsSetUp());
    }

    /* The second public method overrides {{{CreateSrnModel()}}}. This is a
     * builder method to create new copies of the SRN model. We call
     * the (protected) copy constructor which creates a copy of the cell cycle model.
     *
     */
    AbstractSrnModel* CreateSrnModel()
    {
        return new MySrnModel(*this);
    }

    /* The third public method overrides {{{Initialise()}}}. */
    void Initialise()
    {
        AbstractOdeSrnModel::Initialise(new MyOdeSystem);
    }

    /* The fourth public method runs the ODEs at each timestep and saves some results to {{{CellData}}}. */
    void SimulateToCurrentTime()
    {
        // run the ODE simulation as needed
        AbstractOdeSrnModel::SimulateToCurrentTime();

        /* this line outputs the ODE system variable to {{{CellData}}}. */
        mpCell->GetCellData()->SetItem("u",mpOdeSystem->rGetStateVariables()[0]);
    }

    /* Finally we define a method to output any parameters in our model, this needs to be included in every SRN model.*/
    void OutputSrnModelParameters(out_stream& rParamsFile)
    {
        // No new parameters to output, so just call method on direct parent class
        AbstractOdeSrnModel::OutputSrnModelParameters(rParamsFile);
    }
};

/* We need to include the next block of code if you want to be able to archive (save or load)
 * the SRN model object in a cell-based simulation. It is also required for writing out
 * the parameters file describing the settings for a simulation - it provides the unique
 * identifier for our new SRN model. Thus every SRN model class must provide this,
 * or you'll get errors when running simulations. */
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(MyOdeSystem)
CHASTE_CLASS_EXPORT(MySrnModel)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(MySrnModel)

/* Since we're defining the new SRN model and ODEs within the test file, we need to include the
 * following stanza as well, to make the code work with newer versions of the Boost libraries.
 * Normally the above export declaration would occur in the SRN model's .hpp file, and
 * the following lines would appear in the .cpp file.  See ChasteGuides/BoostSerialization for
 * more information.
 */
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(MyOdeSystem)
CHASTE_CLASS_EXPORT(MySrnModel)

/*
 * Need to re-include this after {{{SerializationExportWrapperForCpp.hpp}}}. This is to export the
 * components that would normally be in a seperate cpp file.
 */
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(MySrnModel)

/*
 * This completes the code for {{{MySrnModel}}}. Note that usually this code would
 * be separated out into a separate declaration in a .hpp file and definition in a .cpp file.
 *
 * === The Tests ===
 *
 * We now define the test class, which inherits from {{{AbstractCellBasedTestSuite}}}.
 */
class TestCreatingAndUsingANewSrnModelTutorial : public AbstractCellBasedTestSuite
{
public:
    /*
     * == Using the SRN model in a cell-based simulation ==
     *
     * We conclude with a brief test demonstrating how {{{MySrnModel}}} can be used
     * in a cell-based simulation.
     */
    void TestOffLatticeSimulationWithMySrnModel()
    {
        /* We use the honeycomb vertex mesh generator to create a vertex mesh.
         */
        HoneycombVertexMeshGenerator generator(2, 2);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        /* Next, we create some cells. First, define the cells vector. */
        std::vector<CellPtr> cells;
        /* We must create a shared_ptr to a {{{CellMutationState}}} with which to bestow the cells.
         * We make use of the macro MAKE_PTR to do this: the first argument is the class and
         * the second argument is the name of the shared_ptr. */
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        /* Then we loop over the nodes. */
        for (unsigned i=0; i<p_mesh->GetNumElements(); i++)
        {
            /* For each node we create a cell with our SRN model and simple Stochastic cell cycle model. */
            UniformG1GenerationalCellCycleModel* p_cell_cycle_model = new UniformG1GenerationalCellCycleModel();
            MySrnModel* p_srn_model = new MySrnModel;

            /* We choose to initialise the concentrations to random levels in each cell. */
            std::vector<double> initial_conditions;
            initial_conditions.push_back(1.0-2.0*RandomNumberGenerator::Instance()->ranf());
            //initial_conditions.push_back(1.0-2.0*RandomNumberGenerator::Instance()->ranf());
            p_srn_model->SetInitialConditions(initial_conditions);

            CellPtr p_cell(new Cell(p_state, p_cell_cycle_model, p_srn_model));
            p_cell->SetCellProliferativeType(p_stem_type);


            /* Now, we define a random birth time, chosen from [-T,0], where
             * T is the typical cell cycle duration
             */
            double birth_time = - RandomNumberGenerator::Instance()->ranf() * p_cell_cycle_model->GetAverageStemCellCycleTime();
            /* We then set the birth time and push the cell back into the vector of cells. */
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        /* Now that we have defined the mesh and cells, we can define the cell population, forces, areas modifier, and simulation
         * in the same way as the other tutorials. */
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestOffLatticeSimulationWithMySrnModel");
        simulator.SetEndTime(10.0);
        simulator.SetSamplingTimestepMultiple(50);

        MAKE_PTR(NagaiHondaForce<2>, p_force);
        simulator.AddForce(p_force);

        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        /* Finally to run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();
    }
    /*
     * To visualize the results, use Paraview. See the UserTutorials/VisualizingWithParaview tutorial for more information
     *
     * Load the file {{{/tmp/$USER/testoutput/TestOffLatticeSimulationWithMySrnModel/results_from_time_0/results.pvd}}},
     * and color by {{{x}}}.
     */
};


#endif //CHASTE_TESTCREATINGANDUSINGCORSONSRNMODEL_HPP
