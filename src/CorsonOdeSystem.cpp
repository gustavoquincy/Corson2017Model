//
// Created by gq on 1/13/19.
//

#include "CorsonOdeSystem.hpp"
#include "CellwiseOdeSystemInformation.hpp"


CorsonOdeSystem::CorsonOdeSystem(std::vector<double> stateVariables) : AbstractOdeSystem(1)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<CorsonOdeSystem>);

    //this->mStateVariables.push_back(1.0);
    SetDefaultInitialCondition(0, 0.0); // soon overwritten

    this->mParameters.push_back(0.2); //this refers to AbstactOdeSystem super class

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

void CorsonOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{


    double u = rY[0];
    double s = this->mParameters[0]; // Shorthand for "this->mParameter("Mean Delta");"

    rDY[0] = (SigmoidalFunction(2*(s)) - u) * 2;
}

template<>
void CellwiseOdeSystemInformation<CorsonOdeSystem>::Initialise()
{
    this->mVariableNames.push_back("Cell State"); // StateVariable name
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0); // will be filled in later

    this->mParameterNames.push_back("Signal");
    this->mParameterUnits.push_back("non-dim");

    this->mInitialised = true;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(CorsonOdeSystem)