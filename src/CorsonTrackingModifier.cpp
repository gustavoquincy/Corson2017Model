//
// Created by gq on 1/13/19.
//

#include "CorsonTrackingModifier.hpp"
#include "CorsonSrnModel.hpp"

static const int mN = 324;

static const double mlambda = sqrt(1/mN);

static const double ml1 = 1.75 * mlambda;

static const double mD = 5e-5;

static const double ma0 = 5e-2;

static const double ma1 = 1 - ma0;

static const int mS0 = 2;

static const double mL1 = .2;

static const int mtau_g = 1;

template<unsigned DIM>
CorsonTrackingModifier<DIM>::CorsonTrackingModifier()
        : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
CorsonTrackingModifier<DIM>::~CorsonTrackingModifier()
{
}

template<unsigned DIM>
void CorsonTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void CorsonTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
double CorsonTrackingModifier<DIM>::SigmoidalFunction(double x) const {
    return (1 + tanh(2*x))/2;
}

template<unsigned DIM>
double CorsonTrackingModifier<DIM>::AutoSignalingGradient(double x, double time, double width) const
{
    double L = mL1 * width;
    double l = ml1 * width;
    return mS0* SigmoidalFunction(1-time/mtau_g)*(exp(-pow(x,2)/2/pow(L,2))+exp(-pow(width-x,2)/2/pow(L,2))) +
           SigmoidalFunction(time/mtau_g-1)*(exp(-pow(x,2)/2/pow(l,2))+exp(-pow(width-x,2)/2/pow(l,2)));
}


//intermediate definition for SignalingProductionFunction
template<unsigned DIM>
double CorsonTrackingModifier<DIM>::LigandActivityFunction(double u)  const
{
    return ma0 + 3*pow(u,3)/(1+pow(u,2))*ma1;
}

template<unsigned DIM>
double CorsonTrackingModifier<DIM>::SignalProductionFunction(double u) const
{
    return u * LigandActivityFunction(u);
}

template<unsigned DIM>
void CorsonTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{

    rCellPopulation.Update();

    //Get cell population width
    double width = rCellPopulation.GetWidth(0);//1

    //Access simulation time by quirying SimulationTime Class
    SimulationTime* p_simulation_time = SimulationTime::Instance();
    double t = p_simulation_time->GetTime();

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        CorsonSrnModel* p_model = static_cast<CorsonSrnModel*>(cell_iter->GetSrnModel());
        double u = p_model->GetCellStateParameter();

        cell_iter->GetCellData()->SetItem("cellstate u", u);
    }


    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        c_vector<double, DIM> centroid = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

        //reinitialize signal received from other cells after every loop
        double signal_received = 0;

        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter_2 = rCellPopulation.Begin();
             cell_iter_2 != rCellPopulation.End();
             ++cell_iter_2)
        {
            c_vector<double, DIM> centroid_2 = rCellPopulation.GetLocationOfCellCentre(*cell_iter_2);

            double distance_squared = pow(centroid[0]-centroid_2[0], 2) + pow(centroid[1]-centroid_2[1], 2);

            double coefficient = exp(-distance_squared/2/pow(ml1*width, 2));

            double this_cell_state = cell_iter_2->GetCellData()->GetItem("cellstate u");

            signal_received += coefficient * SignalProductionFunction(this_cell_state);
        }

        double x = centroid[0];

        //incorporate autonomous varied signaling gradient and signal received from other cells
        double s = AutoSignalingGradient(x, t, width) + signal_received;

        cell_iter->GetCellData()->SetItem("signaling s", s);

    }


}

template<unsigned DIM>
void CorsonTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class CorsonTrackingModifier<1>;
template class CorsonTrackingModifier<2>;
template class CorsonTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CorsonTrackingModifier)
