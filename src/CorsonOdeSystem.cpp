//
// Created by gq on 1/13/19.
//

#include "CorsonOdeSystem.hpp"
#include "CellwiseOdeSystemInformation.hpp"


CorsonOdeSystem::CorsonOdeSystem(std::vector<double> stateVariables) : AbstractOdeSystem(1)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<CorsonOdeSystem>);

    SetDefaultInitialCondition(0, 0.00); // soon overwritten

    this->mParameters.push_back(0.0); //

    if (stateVariables != std::vector<double>())
    {
        SetStateVariables(stateVariables);
    }

}

CorsonOdeSystem::~CorsonOdeSystem()
{
}

double CorsonOdeSystem::SigmoidalFunction(double x) const
{
    return (1 + tanh(2*x))/2;
}

void CorsonOdeSystem::SetTimeScaleParameter(double timeScaleParameter) {
    mTimeScaleParameter = timeScaleParameter;
}

double CorsonOdeSystem::GetTimeScaleParameter() {
    return mTimeScaleParameter;
}

void CorsonOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{

    double tau = GetTimeScaleParameter();
    double u = rY[0];
    double s = this->mParameters[0];


    rDY[0] = (SigmoidalFunction(2*(u-s)) - u) / tau ; //tau=.5 i.e. 1/.5=2, here set 1/tau = 10 for faster pattern resolution
}

template<>
void CellwiseOdeSystemInformation<CorsonOdeSystem>::Initialise()
{
    this->mVariableNames.push_back("Cell State");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.00);

    this->mParameterNames.push_back("Signal");
    this->mParameterUnits.push_back("non-dim");
    //Initial value is not needed since s will be calculated in CorsonTrackingModifier

    this->mInitialised = true;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(CorsonOdeSystem)