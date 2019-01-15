//
// Created by gq on 1/13/19.
//

#include "CorsonOdeSystem.hpp"
#include "CellwiseOdeSystemInformation.hpp"

CorsonOdeSystem::CorsonOdeSystem(std::vector<double> stateVariables)
        : AbstractOdeSystem(1)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<CorsonOdeSystem>);

    SetDefaultInitialCondition(0, 1.0); // soon overwritten

    this->mParameters.push_back(1.0);

    if (stateVariables != std::vector<double>())
    {
        SetStateVariables(stateVariables);
    }
}

CorsonOdeSystem::~CorsonOdeSystem()
{
}

void CorsonOdeSystem::Init()
{
    n = 324;
    lambda = sqrt(1/n);
    l_over_one = 1.75 * lambda;
    tau = 1/2;
    D = 5e-5;
    a0 = 5e-2;
    a1 = 1-a0;
    S0 = 2;
    L_over_one = .2;
    tau_g = 1;
}

double SigmoidalFunction(double x)
{
    return (1 + tanh(x))/2;
}

double AutoSignalingGradient(double time, double x, double width)
{
    double L = L_over_one * width;
    double l = l_over_one * width;
    return S0*SigmoidalFunction(1-time/tau_g)*(exp(-pow(x,2)/2/pow(L,2))+exp(-pow(width-x,2)/2/pow(L,2))) +
            SigmoidalFunction(time/tau_g-1)*(exp(-pow(x,2)/2/pow(l,2))+exp(-pow(width-x,2)/2/pow(l,2)));
}

double LigandLevelFunction(double CellState) {
    return CellState;
}

double LigandActivityFunction(double CellState){
    return a0 + 3*pow(CellState,3)/(1+pow(CellState,2))*a1;
}

double SignalProductionFunction(double u){
    return LigandLevelFunction(u)*LigandActivityFunction(u);
}

void CorsonOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    double u = rY[0];
    double s = this->mParameters[0]; // Shorthand for "this->mParameter("Mean Delta");"

    rDY[0] = (SigmoidalFunction(2*(u-s)) - u) / tau;
}

template<>
void CellwiseOdeSystemInformation<CorsonOdeSystem>::Initialise()
{
    this->mVariableNames.push_back("Corson Cell State");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0); // will be filled in later

    this->mVariableNames.push_back("Corson Signaling");
    this->mVariableUnits.push_back("non-dim");

    this->mInitialised = true;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(CorsonOdeSystem)