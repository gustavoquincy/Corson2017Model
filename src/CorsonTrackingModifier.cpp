//
// Created by gq on 1/13/19.
//

#include "CorsonTrackingModifier.hpp"
#include "CorsonSrnModel.hpp"
#include "MathsCustomFunctions.hpp"
//static const int mN = 324;

//static const double mlambda = sqrt(1/mN);

//static const double ml1 = 1.75 * mlambda;

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
double CorsonTrackingModifier<DIM>::AutoSignalingGradient(double x, double time, double l) const
{
    //double L = mL1 * width;
    //double l = l1 * width;
    return mS0* SigmoidalFunction(1-time/mtau_g)*(exp(-SmallPow(x,2)/2/SmallPow(mL1,2))+exp(-SmallPow(1-x,2)/2/SmallPow(mL1,2))) +
           SigmoidalFunction(time/mtau_g-1)*(exp(-SmallPow(x,2)/2/SmallPow(l,2))+exp(-SmallPow(1-x,2)/2/SmallPow(l,2)));
}


//intermediate definition for SignalingProductionFunction
template<unsigned DIM>
double CorsonTrackingModifier<DIM>::LigandActivityFunction(double u)  const
{
    return ma0 + 3*SmallPow(u,3)/(1+2*u)*ma1; //Different from Eq.11
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

    int N = rCellPopulation.GetNumAllCells();

    double lambda = sqrt(1/N);

    double l1 = 1.75 * lambda;

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

        assert(u >= 0);
        assert (u <= 1);

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

            double distance_squared = SmallPow((centroid[0]-centroid_2[0])/width, 2) + SmallPow((centroid[1]-centroid_2[1])/width, 2);

            double coefficient;

            if (distance_squared == 0)
            {
                coefficient = 0; //special case that cii=0
            }
            else
            {
                coefficient = exp(-distance_squared/2/SmallPow(l1, 2));
            }

            double the_other_cell_u = cell_iter_2->GetCellData()->GetItem("cellstate u");

            assert (the_other_cell_u>=0);
            assert (the_other_cell_u<=1);

            double D_aster = coefficient * SignalProductionFunction(the_other_cell_u);

            assert(D_aster>=0);
            assert(D_aster<=1);

            signal_received += D_aster;
        }

        double x = centroid[0]/width;

        //incorporate autonomous varied signaling gradient and signal received from other cells
        double s = AutoSignalingGradient(x, t, l1) + signal_received;

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
