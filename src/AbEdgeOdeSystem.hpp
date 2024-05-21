#ifndef ABEDGEODESYSTEM_HPP_
#define ABEDGEODESYSTEM_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include <cmath>
#include <iostream>

#include "AbstractOdeSystem.hpp"

class AbEdgeOdeSystem : public AbstractOdeSystem
{
private:

    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractOdeSystem>(*this);
    }
public:

    /**
     * Default constructor.
     *
     * @param stateVariables optional initial conditions for state variables 
     *     (only used in archiving)
     */
    AbEdgeOdeSystem(std::vector<double> stateVariables = std::vector<double>());

    /**
     * Destructor.
     */
    ~AbEdgeOdeSystem();

    /**
     * Overridden EvaluateYDerivatives() method.
     * 
     * @param time used to evaluate the RHS.
     * @param rY value of the solution vector used to evaluate the RHS.
     * @param rDY filled in with the resulting derivatives
     */
    void EvaluateYDerivatives(
        double time,
        const std::vector<double>& rY,
        std::vector<double>& rDY) override;
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(AbEdgeOdeSystem)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a AbEdgeOdeSystem.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const AbEdgeOdeSystem * t, const unsigned int file_version)
{
    const std::vector<double>& state_variables = t->rGetConstStateVariables();
    ar & state_variables;
}

/**
 * De-serialize constructor parameters and initialise a AbEdgeOdeSystem.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, AbEdgeOdeSystem * t, const unsigned int file_version)
{
    std::vector<double> state_variables;
    ar & state_variables;

    // Invoke inplace constructor to initialise instance
    ::new(t)AbEdgeOdeSystem(state_variables);
}
}
} // namespace ...

#endif /*ABEDGEODESYSTEM_HPP_*/
