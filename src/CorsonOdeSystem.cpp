//
// Created by gq on 1/13/19.
//

#include "CorsonOdeSystem.hpp"
#include "CellwiseOdeSystemInformation.hpp"

static const double mtau = 1/2; //1/2; make the gradient flatter

CorsonOdeSystem::CorsonOdeSystem(std::vector<double> stateVariables)
        : AbstractOdeSystem(1)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<CorsonOdeSystem>);

    SetDefaultInitialCondition(0, 1.0); // soon overwritten

    this->mParameters.push_back(0.5); //this refers to AbstactOdeSystem super class

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

    rDY[0] = (SigmoidalFunction(2*(u-s)) - u) / mtau;
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