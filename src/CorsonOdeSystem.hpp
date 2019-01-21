//
// Created by gq on 1/13/19.
//

#ifndef CORSONODESYSTEM_HPP_
#define CORSONODESYSTEM_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include <cmath>
#include <iostream>

#include "AbstractOdeSystem.hpp"


class CorsonOdeSystem : public AbstractOdeSystem
{
private:


    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractOdeSystem>(*this);
    }

    double SigmoidalFunction(double x) const;

public:

    CorsonOdeSystem(std::vector<double> stateVariables=std::vector<double>());

    ~CorsonOdeSystem();

    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY);


};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(CorsonOdeSystem)

namespace boost
{
    namespace serialization
    {
/**
 * Serialize information required to construct a CorsonOdeSystem.
 */
        template<class Archive>
        inline void save_construct_data(
                Archive & ar, const CorsonOdeSystem * t, const unsigned int file_version)
        {
            const std::vector<double>& state_variables = t->rGetConstStateVariables();
            ar & state_variables;
        }

/**
 * De-serialize constructor parameters and initialise a CorsonOdeSystem.
 */
        template<class Archive>
        inline void load_construct_data(
                Archive & ar, CorsonOdeSystem * t, const unsigned int file_version)
        {
            std::vector<double> state_variables;
            ar & state_variables;

            // Invoke inplace constructor to initialise instance
            ::new(t)CorsonOdeSystem(state_variables);
        }
    }
} // namespace ...

#endif //CORSONODESYSTEM_HPP_
