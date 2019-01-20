//
// Created by gq on 1/13/19.
//

#ifndef CORSONSRNMODEL_HPP_
#define CORSONSRNMODEL_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "CorsonOdeSystem.hpp"
#include "AbstractOdeSrnModel.hpp"

class CorsonSrnModel : public AbstractOdeSrnModel
{
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractOdeSrnModel>(*this);
    }

protected:
    CorsonSrnModel(const CorsonSrnModel& rModel);

public:

    CorsonSrnModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver = boost::shared_ptr<AbstractCellCycleModelOdeSolver>());

    AbstractSrnModel* CreateSrnModel();

    void Initialise();

    void SimulateToCurrentTime();

    void UpdateCorson();

    double GetCellStateParameter();

    double GetSignalingParameter();

    void OutputSrnModelParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(CorsonSrnModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(CorsonSrnModel)

#endif //CORSONSRNMODEL_HPP_
