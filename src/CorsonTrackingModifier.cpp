//
// Created by gq on 1/13/19.
//

#include "CorsonTrackingModifier.hpp"
#include "CorsonSrnModel.hpp"


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
void CorsonTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{

    rCellPopulation.Update();

    double width = rCellPopulation.GetWidth(DIM);

    //unsure if simulationtime information can be accessed here
    SimulationTime* p_simulation_time = SimulationTime::Instance();

    double t = p_simulation_time->GetTime();

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        CorsonSrnModel* p_model = static_cast<CorsonSrnModel*>(cell_iter->GetSrnModel());
        double u = p_model->GetCellStateParameter();
        double s = p_model->GetSignalingParameter();

        cell_iter->GetCellData()->SetItem("corson signaling", s);
        cell_iter->GetCellData()->SetItem("corson cell state", u);
    }

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {

        c_vector<double, DIM> centroid = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

        //reinitialize signal received from other cells in every loop
        double signal = 0;

        for (typename AbstractCellPopulation<DIM>::Iterator
                     cell_iter_2 = rCellPopulation.Begin();
             cell_iter_2 != rCellPopulation.End();
             ++cell_iter_2)
        {
            c_vector<double, DIM> centroid_2 = rCellPopulation.GetLocationOfCellCentre(*cell_iter_2);
            double distance_squared = pow(centroid[0]-centroid_2[0],2) + pow(centroid[1]-centroid_2[1],2);
            double coefficient = exp(-distance_squared/2/pow(l_over_one * width,2));
            double this_cell_state = cell_iter_2->GetCellData()->GetItem("corson cell state");
            signal += coefficient * SignalProductionFunction(this_cell_state);
        }

        double x = rCellPopulation.GetLocationOfCellCentre(*cell_iter)[0];
        double s0 = AutoSignalingGradient(x,t) + signal;
        cell_iter->GetCellData()->SetItem("corson signaling", s0);

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
