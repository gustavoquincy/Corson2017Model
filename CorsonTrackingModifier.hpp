//
// Created by gq on 1/13/19.
//

#ifndef CORSONTRACKINGMODIFIER_HPP_
#define CORSONTRACKINGMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCellBasedSimulationModifier.hpp"


template<unsigned DIM>
class CorsonTrackingModifier : public AbstractCellBasedSimulationModifier<DIM,DIM>
{
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM,DIM> >(*this);
    }
private:

    const int mN = 324;
    const double mlambda = sqrt(1/mN);
    const double ml1 = 1.75 * mlambda;
    const double mD = 5e-5;
    const double ma0 = 5e-2;
    const double ma1 = 1 - ma0;
    const int mS0 = 2;
    const double mL1 = .2;
    const int mtau_g = 1;

    double SigmoidalFunction(double x) const;

public:

    CorsonTrackingModifier();

    virtual ~CorsonTrackingModifier();

    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    virtual void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory);

    void UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    void OutputSimulationModifierParameters(out_stream& rParamsFile);

    double AutoSignalingGradient(double x, double time, double width) const;

    double LigandActivityFunction(double u) const;

    double SignalProductionFunction(double u) const;
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CorsonTrackingModifier)

#endif //CHASTE_CORSONTRACKINGMODIFIER_HPP
