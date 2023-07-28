#ifndef FTDSEDGEODESYSTEM_HPP_
#define FTDSEDGEODESYSTEM_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include <cmath>
#include <iostream>

#include "AbstractOdeSystem.hpp"

class FtDsEdgeOdeSystem : public AbstractOdeSystem
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
    FtDsEdgeOdeSystem(std::vector<double> stateVariables = std::vector<double>());

    /**
     * Destructor.
     */
    ~FtDsEdgeOdeSystem();

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
CHASTE_CLASS_EXPORT(FtDsEdgeOdeSystem)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a FtDsEdgeOdeSystem.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const FtDsEdgeOdeSystem * t, const unsigned int file_version)
{
    const std::vector<double>& state_variables = t->rGetConstStateVariables();
    ar & state_variables;
}

/**
 * De-serialize constructor parameters and initialise a FtDsEdgeOdeSystem.
 */
template<class Archive>
inline void load_construct_data(
    Archive & ar, FtDsEdgeOdeSystem * t, const unsigned int file_version)
{
    std::vector<double> state_variables;
    ar & state_variables;

    // Invoke inplace constructor to initialise instance
    ::new(t)FtDsEdgeOdeSystem(state_variables);
}
}
} // namespace ...

#endif /*FTDSEDGEODESYSTEM_HPP_*/
