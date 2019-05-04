//
// Created by gq on 1/13/19.
//

#include "CorsonTrackingModifier.hpp"
#include "CorsonSrnModel.hpp"
#include "MathsCustomFunctions.hpp"
#include "Debug.hpp"
#include "FarhadifarForce.hpp"
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
    return ma0 + 3*SmallPow(u,3)/(2+ SmallPow(u,2))*ma1;
}

template<unsigned DIM>
double CorsonTrackingModifier<DIM>::SignalProductionFunction(double u) const
{
    //double uc = 0.4;
    //int k = 50;
    return u * LigandActivityFunction(u); //SigmoidalFunction(k*(u-uc)) (1) D(u) = u -> D(u) = .1 u i.e. much weaker mutual inhibition
}

template<unsigned DIM>
double CorsonTrackingModifier<DIM>::TargetAreaControllingFunction(double u) const
{
    return 2*SmallPow(u,2)+2*u+1;
}
template<unsigned DIM>
void CorsonTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{

    rCellPopulation.Update();

    double N = rCellPopulation.GetNumAllCells();

    double lambda = sqrt(1/N);

    double l1 = 1.25 * lambda; //1.75, 1.25

    mSignalingRangeParameter = l1/lambda;

    //Get cell population width
    double width = rCellPopulation.GetWidth(0);

    double height = rCellPopulation.GetWidth(1);

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        CorsonSrnModel* p_model = static_cast<CorsonSrnModel*>(cell_iter->GetSrnModel());
        double u = p_model->GetCellStateParameter();
        cell_iter->GetCellData()->SetItem("cellstate u", u);
        cell_iter->GetCellData()->SetItem("target area", (2*u+1)*0.5*sqrt(3)); //2*SmallPow(u,2)+2*u+1

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
            double dx = abs(centroid[0]-centroid_2[0]) <= width/2.0 ? abs(centroid[0]-centroid_2[0]) : width-abs(centroid[0]-centroid_2[0]);
            double dy = abs(centroid[1]-centroid_2[1]) <= height/2.0 ? abs(centroid[1]-centroid_2[1]) : height-abs(centroid[1]-centroid_2[1]);
            double distance = sqrt(SmallPow(dx,2) + SmallPow(dy,2));
            double coefficient;
            if (distance != 0)
            {
                coefficient = exp(-SmallPow(distance/width, 2)/SmallPow(l1, 2)/2);
            }
            else
            {
                coefficient = 0; //special case that cii=0
            }
            double the_other_cell_u = cell_iter_2->GetCellData()->GetItem("cellstate u");
            signal_received += coefficient * SignalProductionFunction(the_other_cell_u);
        }

        //Access simulation time by quirying SimulationTime Class
        /*SimulationTime* p_simulation_time = SimulationTime::Instance();
        double t = p_simulation_time->GetTime();
        double x = centroid[0]/width;
        double s =  signal_received + AutoSignalingGradient(x, t, l1); //no autosignaling gradient*/

        double s = signal_received; // no autosignaling gradient version
        cell_iter->GetCellData()->SetItem("signaling s", s);

    }

   /* //Hook up cell state u with cell target area
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
        cell_iter != rCellPopulation.End();
        ++cell_iter)
    {

    }*/
}

template<unsigned DIM>
void CorsonTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<S0>" << mS0 << "</S0>\n";
    *rParamsFile << "\t\t\t<L1>" << mL1 << "</L1>\n";
    *rParamsFile << "\t\t\t<tau_g>" << mtau_g << "</tau_g>\n";
    *rParamsFile << "\t\t\t<Signaling Range l/lambda>" << mSignalingRangeParameter << "</Signaling Range l/lambda>\n";
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
