//
// Created by gq on 1/13/19.
//

#ifndef CORSONTRACKINGMODIFIER_HPP_
#define CORSONTRACKINGMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

#include <cmath>
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

    //Handy function to represent one term
    double SigmoidalFunction(double x) const;

    // eq(8), s0(xi,t)
    double AutoSignalingGradient(double x, double time, double width) const;

    //eq(11) a(u), D*(u) = a(u)D(u) = a(u)u
    double LigandActivityFunction(double u) const;

    //eq(10), D*(u)
    double SignalProductionFunction(double u) const;

public:
    //Default constructor
    CorsonTrackingModifier();

    //Destructor
    virtual ~CorsonTrackingModifier();

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * Specifies what to do in the simulation at the end of each time step.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Overridden SetupSolve() method.
     *
     * Specifies what to do in the simulation before the start of the time loop.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory);


    //Helper method to compute the signal level of each cell, including the time dependent gradient
    // and signal received from all the cells in the sheet
    void UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
    * Overridden OutputSimulationModifierParameters() method.
    * Output any simulation modifier parameters to file.
    *
    * @param rParamsFile the file stream to which the parameters are output
    */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);

};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CorsonTrackingModifier)

#endif //CHASTE_CORSONTRACKINGMODIFIER_HPP
