//
// Created by gq on 1/13/19.
//

#ifndef CHASTE_CORSONTRACKINGMODIFIER_HPP
#define CHASTE_CORSONTRACKINGMODIFIER_HPP

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


public:

    CorsonTrackingModifier();

    virtual ~CorsonTrackingModifier();

    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    virtual void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory);

    void UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CorsonTrackingModifier)

#endif //CHASTE_CORSONTRACKINGMODIFIER_HPP
