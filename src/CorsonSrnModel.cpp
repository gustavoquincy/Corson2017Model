//
// Created by gq on 1/13/19.
//

#include "CorsonSrnModel.hpp"

CorsonSrnModel::CorsonSrnModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver) : AbstractOdeSrnModel(1, pOdeSolver)
{
    if (mpOdeSolver == boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
    {
#ifdef CHASTE_CVODE
        mpOdeSolver = CellCycleModelOdeSolver<CorsonSrnModel, CvodeAdaptor>::Instance();
		mpOdeSolver->Initialise();
		mpOdeSolver->SetMaxSteps(10000);
#else
        mpOdeSolver = CellCycleModelOdeSolver<CorsonSrnModel, RungeKutta4IvpOdeSolver>::Instance();
        mpOdeSolver->Initialise();
        SetDt(0.001);
#endif
    }
    assert(mpOdeSolver->IsSetUp());
}

CorsonSrnModel::CorsonSrnModel(const CorsonSrnModel& rModel)
        : AbstractOdeSrnModel(rModel)
{
    assert(rModel.GetOdeSystem());
    SetOdeSystem(new CorsonOdeSystem(rModel.GetOdeSystem()->rGetStateVariables()));
}

AbstractSrnModel* CorsonSrnModel::CreateSrnModel()
{
    return new CorsonSrnModel(*this);
}

void CorsonSrnModel::SimulateToCurrentTime()
{
    UpdateCorson();
    AbstractOdeSrnModel::SimulateToCurrentTime();
}

void CorsonSrnModel::Initialise()
{
    AbstractOdeSrnModel::Initialise(new CorsonOdeSystem);
}

void CorsonSrnModel::UpdateCorson()
{
    assert(mpOdeSystem != nullptr);
    assert(mpCell != nullptr);

    double s = mpCell->GetCellData()->GetItem("corson signaling");
    mpOdeSystem->SetParameter("Corson Signaling", s);
}

double CorsonSrnModel::GetCellStateParameter()
{
    assert(mpOdeSystem != nullptr);
    double u = mpOdeSystem->rGetStateVariables()[0];
    return u;
}

double CorsonSrnModel::GetSignalingParameter()
{
    assert(mpOdeSystem != nullptr);
    double s = mpOdeSystem->GetParameter("Corson Signaling");
    return s;
}

void CorsonSrnModel::OutputSrnModelParameters(out_stream& rParamsFile)
{
    AbstractOdeSrnModel::OutputSrnModelParameters(rParamsFile);
}

#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(CorsonSrnModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(CorsonSrnModel)