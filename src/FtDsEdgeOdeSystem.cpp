#include "CellwiseOdeSystemInformation.hpp"
#include "FtDsEdgeOdeSystem.hpp"

FtDsEdgeOdeSystem::FtDsEdgeOdeSystem(std::vector<double> stateVariables)
    : AbstractOdeSystem(12)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<FtDsEdgeOdeSystem>);

    /**
     * The state variables are as follows:
     *
     * 0 - Notch concentration for this cell edge (Ds)
     * 1 - Delta concentration for this cell edge (Ft)
     * 2 - DsP
     * 3 - FtP
     * 4 - complex A (FtP-neigh_Ds)
     * 5 - complex B (FtP-neigh_DsP)
     * 6 - complex C (Ft-neigh_Ds)
     * 7 - complex D (Ft-neigh_DsP)
     * 8 - complex A_ (neigh_FtP-Ds)
     * 9 - complex B_ (DsP-neigh_FtP)
     * 10 - complex C_ (neigh_Ft-Ds)
     * 11 - complex D_ (neigh_Ft-DsP)
     *
     * We store the last state variable so that it can be written
     * to file at each time step alongside the others, and visualized.
     */

    SetDefaultInitialCondition(0, 0.1); // soon overwritten
    SetDefaultInitialCondition(1, 0.1); // soon overwritten
    SetDefaultInitialCondition(2, 0.1); // soon overwritten
    SetDefaultInitialCondition(3, 0.0); // soon overwritten
    SetDefaultInitialCondition(4, 0.0); // soon overwritten
    SetDefaultInitialCondition(5, 0.0); // soon overwritten
    SetDefaultInitialCondition(6, 0.0); // soon overwritten
    SetDefaultInitialCondition(7, 0.0); // soon overwritten
    SetDefaultInitialCondition(8, 0.0); // soon overwritten
    SetDefaultInitialCondition(9, 0.0); // soon overwritten
    SetDefaultInitialCondition(10, 0.0); // soon overwritten
    SetDefaultInitialCondition(11, 0.0); // soon overwritten

    this->mParameters.push_back(0.0);
    //By default zero. If no interior SRN model is specified, interior delta/notch is zero
    this->mParameters.push_back(0.0);
    this->mParameters.push_back(0.0);
    this->mParameters.push_back(0.0);


    if (stateVariables != std::vector<double>())
    {
        SetStateVariables(stateVariables);
    }
}

FtDsEdgeOdeSystem::~FtDsEdgeOdeSystem()
{

}



void FtDsEdgeOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{


    /** Using the model in Hale et al. (2015) - there are four different proteins:
   unphosphorylated Ft (Ft), phosphorylated Ds (Ds), unphosphorylated Ds (Ds), and phosphorylated Ds (DsP)
   which could then form four differnt complexes:
   FtP-Ds (A), FtP-DsP (B), Ft-Ds (C), and Ft-DsP (D)
   Since Fj appears to have a more dominant effect to Ft, the binding strengths are A > B > C = D
   so the association constants for kon/koff can be given as 1, 1/2, 1/4, and 1/4, respectively.
   */


    const double notch = rY[0];
    const double delta = rY[1];
    const double DsP = rY[2];
    const double FtP = rY[3];
    const double A = rY[4];
    const double B = rY[5];
    const double C = rY[6];
    const double D = rY[7];
    const double A_ = rY[8];
    const double B_ = rY[9];
    const double C_ = rY[10];
    const double D_ = rY[11];
    

    const double neigh_delta = this->mParameters[0]; // Shorthand for "this->mParameter("neighbor delta");"
    const double neigh_notch = this->mParameters[1];
    const double neigh_DsP = this->mParameters[2];
    const double neigh_FtP = this->mParameters[3];

    //const double interior_delta = this->mParameters[1];
    //const double interior_notch = this->mParameters[2];

    //const double ka_on = 1.0;
    //const double kb_on = 1.0;
    //const double kc_on = 1.0;
    //const double kd_on = 1.0;
    //const double ka_off = 1.0;
    //const double kb_off = 2.0;
    //const double kc_off = 4.0;
    //const double kd_off = 4.0;

    rDY[0] = -neigh_FtP * notch -  neigh_delta +  A_ + 2 * C_; // d[Ds]/dt

    rDY[1] = -delta * neigh_notch - delta * neigh_DsP + 4 * C + 4 * D; // d[Ft]/dt

    rDY[2] = -neigh_FtP * DsP -  neigh_delta * DsP + 2 * B_ + 4 * D_; // d[DsP]/dt

    rDY[3] =  FtP * neigh_notch -  FtP * neigh_DsP +  A + 4 * B; // d[FtP]/dt



    rDY[4] =  FtP * neigh_notch -  A; // d[A]/dt

    rDY[5] =  FtP * neigh_DsP - 2 * B; // d[B]/dt

    rDY[6]=  delta * neigh_notch - 4 * C; // d[C]/dt

    rDY[7] =  delta * neigh_DsP - 4 * D; // d[D]/dt



    rDY[8] =  notch * neigh_FtP -  A_; // d[A_]/dt

    rDY[9] =  neigh_FtP * DsP - 2 * B_; // d[B_]/dt

    rDY[10]=  neigh_delta * notch - 4 * C_; // d[C_]/dt

    rDY[11] =  neigh_delta * DsP - 4 * D_; // d[D_]/dt

}

template<>
void CellwiseOdeSystemInformation<FtDsEdgeOdeSystem>::Initialise()
{
    this->mVariableNames.push_back("Notch");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0); // will be filled in later

    this->mVariableNames.push_back("Delta");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0); // will be filled in later

    this->mVariableNames.push_back("DsP");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0); // will be filled in later

    this->mVariableNames.push_back("FtP");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0); // will be filled in later

    this->mVariableNames.push_back("A");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0); // will be filled in later

    this->mVariableNames.push_back("B");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0); // will be filled in later

    this->mVariableNames.push_back("C");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0); // will be filled in later

    this->mVariableNames.push_back("D");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0); // will be filled in later

    this->mVariableNames.push_back("A_");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0); // will be filled in later

    this->mVariableNames.push_back("B_");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0); // will be filled in later

    this->mVariableNames.push_back("C_");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0); // will be filled in later

    this->mVariableNames.push_back("D_");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0); // will be filled in later



    this->mParameterNames.push_back("neighbour delta");
    this->mParameterUnits.push_back("non-dim");
    this->mParameterNames.push_back("neighbour notch");
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
