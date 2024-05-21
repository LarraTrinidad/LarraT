#include "CellwiseOdeSystemInformation.hpp"
#include "AbEdgeOdeSystem.hpp"

AbEdgeOdeSystem::AbEdgeOdeSystem(std::vector<double> stateVariables)
    : AbstractOdeSystem(4)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<AbEdgeOdeSystem>);

    /**
     * The state variables are as follows:
     *
     * 0 - unbound A
     * 1 - unbound B
     * 2 - complex C (A - neigh B)
     * 3 - complex neigh_C (neigh A - B)
     * 
     * The following initial condition is soon overwritten.
     */
    for (unsigned i = 0; i < 4; ++i)
    {
        SetDefaultInitialCondition(i, 0.0);
    }

    for (unsigned i = 0; i < 2; ++i)
    {
        this->mParameters.push_back(0.0);
    }

    if (stateVariables != std::vector<double>())
    {
        SetStateVariables(stateVariables);
    }
}

AbEdgeOdeSystem::~AbEdgeOdeSystem()
{
}

void AbEdgeOdeSystem::EvaluateYDerivatives(
    double time,
    const std::vector<double>& rY,
    std::vector<double>& rDY)
{
    /*
     * The ODE system is given by:
     *
     * d[A]/dt =       -k_on*[A]*[neigh_B] + k_off*[neigh_C],
     * d[B]/dt =       -k_on*[neigh_A]*[B] + k_off*[C],
     * d[C]/dt =        k_on*[neigh_A]*[B] - k_off*[C],
     * d[neigh_C]/dt =  k_on*[A]*[neigh_B] - k_off*[neigh_C].
     */
    const double A = rY[0];
    const double B = rY[1];
    const double C = rY[2];
    const double neigh_C = rY[3];

    const double neigh_A = this->mParameters[0];
    const double neigh_B = this->mParameters[1];
    //const double k_on = this->mParameters[2];
    //const double k_off = this->mParameters[3];

    rDY[0] = -A*neigh_B + neigh_C;
    rDY[1] = -neigh_A*B + C;
    rDY[2] = -rDY[1];
    rDY[3] = -rDY[0];
}

template<>
void CellwiseOdeSystemInformation<AbEdgeOdeSystem>::Initialise()
{
    // this->mInitialConditions will be filled in later
    this->mVariableNames.push_back("A");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("B");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("C");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("neigh_C");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mParameterNames.push_back("neighbour A");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("neighbour B");
    this->mParameterUnits.push_back("non-dim");

    //this->mParameterNames.push_back("k_on");
    //this->mParameterUnits.push_back("non-dim");

    //this->mParameterNames.push_back("k_off");
    //this->mParameterUnits.push_back("non-dim");

    this->mInitialised = true;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(AbEdgeOdeSystem)
