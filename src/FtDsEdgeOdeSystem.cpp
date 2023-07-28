#include "CellwiseOdeSystemInformation.hpp"
#include "FtDsEdgeOdeSystem.hpp"

FtDsEdgeOdeSystem::FtDsEdgeOdeSystem(std::vector<double> stateVariables)
    : AbstractOdeSystem(12)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<FtDsEdgeOdeSystem>);

    /**
     * The state variables are as follows:
     *
     * 0 - Ds
     * 1 - Ft
     * 2 - DsP
     * 3 - FtP
     * 4 - complex A (FtP-neigh_Ds)
     * 5 - complex B (FtP-neigh_DsP)
     * 6 - complex C (Ft-neigh_Ds)
     * 7 - complex D (Ft-neigh_DsP)
     * 8 - complex neigh_A (neigh_FtP-Ds)
     * 9 - complex neigh_B (DsP-neigh_FtP)
     * 10 - complex neigh_C (neigh_Ft-Ds)
     * 11 - complex neigh_D (neigh_Ft-DsP)
     *
     * We store the last state variable so that it can be written
     * to file at each time step alongside the others, and visualized.
     * 
     * The following initial condition is soon overwritten.
     */
    for (unsigned i = 0; i < 12; i++)
    {
        double value = (i < 3) ? 0.1 : 0.0;
        SetDefaultInitialCondition(i, value);
    }

    for (unsigned i = 0; i < 4; i++)
    {
        this->mParameters.push_back(0.0);
    }

    if (stateVariables != std::vector<double>())
    {
        SetStateVariables(stateVariables);
    }
}

FtDsEdgeOdeSystem::~FtDsEdgeOdeSystem()
{
}

void FtDsEdgeOdeSystem::EvaluateYDerivatives(
    double time,
    const std::vector<double>& rY,
    std::vector<double>& rDY)
{
    /*
     * Using the model in Hale et al. (2015), there are four different proteins:
     * unphosphorylated Ft (Ft); phosphorylated Ft (FtP); unphosphorylated Ds 
     * (Ds); and phosphorylated Ds (DsP). These can then form four different 
     * complexes: FtP-Ds (A); FtP-DsP (B); Ft-Ds (C); and Ft-DsP (D). Since Fj 
     * appears to have a more dominant effect to Ft, the binding strengths are 
     * A > B > C = D so the association constants for kon/koff can be given as 
     * 1, 1/2, 1/4, and 1/4, respectively.
     */

    const double Ds = rY[0];
    const double Ft = rY[1];
    const double DsP = rY[2];
    const double FtP = rY[3];
    const double A = rY[4];
    const double B = rY[5];
    const double C = rY[6];
    const double D = rY[7];
    const double neigh_A = rY[8];
    const double neigh_B = rY[9];
    const double neigh_C = rY[10];
    const double neigh_D = rY[11];

    const double neigh_Ds = this->mParameters[0];
    const double neigh_Ft = this->mParameters[1];
    const double neigh_DsP = this->mParameters[2];
    const double neigh_FtP = this->mParameters[3];

    // d[Ds]/dt
    rDY[0] = -neigh_FtP*Ds - neigh_Ft*Ds + neigh_A + 4*neigh_C; 

    // d[Ft]/dt
    rDY[1] = -neigh_Ds*Ft - neigh_DsP*Ft + 4*C + 4*D;

    // d[DsP]/dt
    rDY[2] = -neigh_FtP*DsP - neigh_Ft*DsP + 2*neigh_B + 4*neigh_D;

    // d[FtP]/dt
    rDY[3] =  -neigh_Ds*FtP - neigh_DsP*FtP + A + 2*B;

    // d[A]/dt
    rDY[4] = neigh_Ds*FtP - A;

    // d[B]/dt
    rDY[5] = neigh_DsP*FtP - 2*B;
    
    // d[C]/dt
    rDY[6]= neigh_Ds*Ft - 4*C; 
    
    // d[D]/dt
    rDY[7] = neigh_DsP*Ft - 4*D; 

    // d[neigh_A]/dt
    rDY[8] = neigh_FtP*Ds - neigh_A; 
    
    // d[neigh_B]/dt
    rDY[9] = neigh_FtP*DsP - 2*neigh_B; 
    
    // d[neigh_C]/dt
    rDY[10]= neigh_Ft*Ds - 4*neigh_C; 
    
    // d[neigh_D]/dt
    rDY[11] = neigh_Ft*DsP - 4*neigh_D;
}

template<>
void CellwiseOdeSystemInformation<FtDsEdgeOdeSystem>::Initialise()
{
    // this->mInitialConditions will be filled in later
    this->mVariableNames.push_back("Ds");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("Ft");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("DsP");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("FtP");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("A");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("B");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("C");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("D");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("neigh_A");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("neigh_B");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("neigh_C");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("neigh_D");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mParameterNames.push_back("neighbour Ds");
    this->mParameterUnits.push_back("non-dim");
    this->mParameterNames.push_back("neighbour Ft");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("neighbour DsP");
    this->mParameterUnits.push_back("non-dim");
    this->mParameterNames.push_back("neighbour FtP");
    this->mParameterUnits.push_back("non-dim");

    this->mInitialised = true;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(FtDsEdgeOdeSystem)
