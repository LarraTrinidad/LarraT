#include "FtDsEdgeSrnModel.hpp"

FtDsEdgeSrnModel::FtDsEdgeSrnModel(
    boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : AbstractOdeSrnModel(12, pOdeSolver)
{
    if (mpOdeSolver == boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
    {
#ifdef CHASTE_CVODE
        mpOdeSolver = CellCycleModelOdeSolver<FtDsEdgeSrnModel, CvodeAdaptor>::Instance();
        mpOdeSolver->Initialise();
        mpOdeSolver->SetMaxSteps(10000);
#else
        mpOdeSolver = CellCycleModelOdeSolver<FtDsEdgeSrnModel, RungeKutta4IvpOdeSolver>::Instance();
        mpOdeSolver->Initialise();
        SetDt(0.001);
#endif //CHASTE_CVODE
    }
    assert(mpOdeSolver->IsSetUp());
}

FtDsEdgeSrnModel::FtDsEdgeSrnModel(const FtDsEdgeSrnModel& rModel)
    : AbstractOdeSrnModel(rModel)
{
    /*
     * Set each member variable of the new SRN model that inherits its value 
     * from the parent.
     *
     * Note 1: some of the new SRN model's member variables will already have 
     * been correctly initialized in its constructor.
     *
     * Note 2: one or more of the new SRN model's member variables may be set/
     * overwritten as soon as InitialiseDaughterCell() is called on the new SRN 
     * model.
     *
     * Note 3: Only set the variables defined in this class. Variables defined
     * in parent classes will be defined there.
     */
    assert(rModel.GetOdeSystem());
    AbstractOdeSystem* p_parent_system(rModel.GetOdeSystem());
    SetOdeSystem(new FtDsEdgeOdeSystem(p_parent_system->rGetStateVariables()));
    for (unsigned i = 0; i < p_parent_system->GetNumberOfParameters(); ++i)
    {
        mpOdeSystem->SetParameter(i, p_parent_system->GetParameter(i));
    }
}

AbstractSrnModel* FtDsEdgeSrnModel::CreateSrnModel()
{
    return new FtDsEdgeSrnModel(*this);
}

void FtDsEdgeSrnModel::SimulateToCurrentTime()
{
    // Update information before running simulation
    UpdateFtDs();

    // Run the ODE simulation as needed
    AbstractOdeSrnModel::SimulateToCurrentTime();
}

void FtDsEdgeSrnModel::Initialise()
{
    AbstractOdeSrnModel::Initialise(new FtDsEdgeOdeSystem);
}

void FtDsEdgeSrnModel::InitialiseDaughterCell()
{
    assert(mpOdeSystem != nullptr);
    assert(mpCell != nullptr);

    // A new edge is initialised with zero concentrations
    mpOdeSystem->SetStateVariable("Ds", 0.0);
    mpOdeSystem->SetStateVariable("Ft", 0.0);
    mpOdeSystem->SetStateVariable("DsP", 0.0);
    mpOdeSystem->SetStateVariable("FtP", 0.0);
    mpOdeSystem->SetStateVariable("A", 0.0);
    mpOdeSystem->SetStateVariable("B", 0.0);
    mpOdeSystem->SetStateVariable("C", 0.0);
    mpOdeSystem->SetStateVariable("D", 0.0);
    mpOdeSystem->SetStateVariable("neigh_A", 0.0);
    mpOdeSystem->SetStateVariable("neigh_B", 0.0);
    mpOdeSystem->SetStateVariable("neigh_C", 0.0);
    mpOdeSystem->SetStateVariable("neigh_D", 0.0);

    mpOdeSystem->SetParameter("neighbour Ds", 0.0);
    mpOdeSystem->SetParameter("neighbour Ft", 0.0);
    mpOdeSystem->SetParameter("neighbour DsP", 0.0);
    mpOdeSystem->SetParameter("neighbour FtP", 0.0);
}

void FtDsEdgeSrnModel::UpdateFtDs()
{
    assert(mpOdeSystem != nullptr);
    assert(mpCell != nullptr);

    auto p_data = mpCell->GetCellEdgeData();

    double neigh_Ds = p_data->GetItem("neighbour Ds")[this->GetEdgeLocalIndex()];
    mpOdeSystem->SetParameter("neighbour Ds", neigh_Ds);

    double neigh_Ft = p_data->GetItem("neighbour Ft")[this->GetEdgeLocalIndex()];
    mpOdeSystem->SetParameter("neighbour Ft", neigh_Ft);

    double neigh_DsP = p_data->GetItem("neighbour DsP")[this->GetEdgeLocalIndex()];
    mpOdeSystem->SetParameter("neighbour DsP", neigh_DsP);

    double neigh_FtP = p_data->GetItem("neighbour FtP")[this->GetEdgeLocalIndex()];
    mpOdeSystem->SetParameter("neighbour FtP", neigh_FtP);
}

double FtDsEdgeSrnModel::GetDs()
{
    assert(mpOdeSystem != nullptr);
    double Ds = mpOdeSystem->rGetStateVariables()[0];
    return Ds;
}

void FtDsEdgeSrnModel::SetDs(double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[0] = value;
}

double FtDsEdgeSrnModel::GetFt()
{
    assert(mpOdeSystem != nullptr);
    double Ft = mpOdeSystem->rGetStateVariables()[1];
    return Ft;
}

void FtDsEdgeSrnModel::SetFt(double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[1] = value;
}

double FtDsEdgeSrnModel::GetDsP()
{
    assert(mpOdeSystem != nullptr);
    double DsP = mpOdeSystem->rGetStateVariables()[2];
    return DsP;
}

void FtDsEdgeSrnModel::SetDsP(double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[2] = value;
}

double FtDsEdgeSrnModel::GetFtP()
{
    assert(mpOdeSystem != nullptr);
    double FtP = mpOdeSystem->rGetStateVariables()[3];
    return FtP;
}

void FtDsEdgeSrnModel::SetFtP(double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[3] = value;
}

double FtDsEdgeSrnModel::GetA()
{
    assert(mpOdeSystem != nullptr);
    double A = mpOdeSystem->rGetStateVariables()[4];
    return A;
}

void FtDsEdgeSrnModel::SetA(double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[4] = value;
}

double FtDsEdgeSrnModel::GetB()
{
    assert(mpOdeSystem != nullptr);
    double B = mpOdeSystem->rGetStateVariables()[5];
    return B;
}

void FtDsEdgeSrnModel::SetB(double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[5] = value;
}

double FtDsEdgeSrnModel::GetC()
{
    assert(mpOdeSystem != nullptr);
    double C = mpOdeSystem->rGetStateVariables()[6];
    return C;
}

void FtDsEdgeSrnModel::SetC(double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[6] = value;
}

double FtDsEdgeSrnModel::GetD()
{
    assert(mpOdeSystem != nullptr);
    double D = mpOdeSystem->rGetStateVariables()[7];
    return D;
}

void FtDsEdgeSrnModel::SetD(double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[7] = value;
}

double FtDsEdgeSrnModel::GetNeighA()
{
    assert(mpOdeSystem != nullptr);
    double neigh_A = mpOdeSystem->rGetStateVariables()[8];
    return neigh_A;
}

void FtDsEdgeSrnModel::SetNeighA(double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[8] = value;
}

double FtDsEdgeSrnModel::GetNeighB()
{
    assert(mpOdeSystem != nullptr);
    double neigh_B = mpOdeSystem->rGetStateVariables()[9];
    return neigh_B;
}

void FtDsEdgeSrnModel::SetNeighB(double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[9] = value;
}

double FtDsEdgeSrnModel::GetNeighC()
{
    assert(mpOdeSystem != nullptr);
    double neigh_C = mpOdeSystem->rGetStateVariables()[10];
    return neigh_C;
}

void FtDsEdgeSrnModel::SetNeighC(double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[10] = value;
}

double FtDsEdgeSrnModel::GetNeighD()
{
    assert(mpOdeSystem != nullptr);
    double neigh_D = mpOdeSystem->rGetStateVariables()[11];
    return neigh_D;
}

void FtDsEdgeSrnModel::SetNeighD(double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[11] = value;
}

void FtDsEdgeSrnModel::OutputSrnModelParameters(out_stream& rParamsFile)
{
    AbstractOdeSrnModel::OutputSrnModelParameters(rParamsFile);
}

void FtDsEdgeSrnModel::AddSrnQuantities(
    AbstractSrnModel* pOtherSrn,
    const double scale)
{
    auto p_other_srn = static_cast<FtDsEdgeSrnModel*>(pOtherSrn);

    const double other_Ds = p_other_srn->GetDs();
    const double other_Ft = p_other_srn->GetFt();
    const double other_DsP = p_other_srn->GetDsP();
    const double other_FtP = p_other_srn->GetFtP();
    const double other_A = p_other_srn->GetA();
    const double other_B = p_other_srn->GetB();
    const double other_C = p_other_srn->GetC();
    const double other_D = p_other_srn->GetD();
    const double other_neigh_A = p_other_srn->GetNeighA();
    const double other_neigh_B = p_other_srn->GetNeighB();
    const double other_neigh_C = p_other_srn->GetNeighC();
    const double other_neigh_D = p_other_srn->GetNeighD();

    const double this_Ds = GetDs();
    const double this_Ft = GetFt();
    const double this_DsP = GetDsP();
    const double this_FtP = GetFtP();
    const double this_A = GetA();
    const double this_B = GetB();
    const double this_C = GetC();
    const double this_D = GetD();
    const double this_neigh_A = GetNeighA();
    const double this_neigh_B = GetNeighB();
    const double this_neigh_C = GetNeighC();
    const double this_neigh_D = GetNeighD();

    ///\todo This is unclear
    SetDs(this_Ds + scale*other_Ds);
    SetFt(this_Ft + scale*other_Ft);
    SetDsP(this_DsP + scale*other_DsP);
    SetFtP(this_FtP + scale*other_FtP);
    SetA(this_A + scale*other_A);
    SetB(this_B + scale*other_B);
    SetC(this_C + scale*other_C);
    SetD(this_D + scale*other_D);
    SetA(this_neigh_A + scale*other_neigh_A);
    SetB(this_neigh_B + scale*other_neigh_B);
    SetC(this_neigh_C + scale*other_neigh_C);
    SetD(this_neigh_D + scale*other_neigh_D);
}

void FtDsEdgeSrnModel::AddShrunkEdgeSrn(AbstractSrnModel* pShrunkEdgeSrn)
{
    /*
     * Here we assume that one half of SRN quantities are endocytosed and the 
     * remaining half are split between neighbouring junctions. Hence we add 
     * 1/4 of SRN variables.
     */
    AddSrnQuantities(pShrunkEdgeSrn, 0.25);
}

void FtDsEdgeSrnModel::AddMergedEdgeSrn(AbstractSrnModel* pMergedEdgeSrn)
{
    // Add all srn variables to this edge srn
    AddSrnQuantities(pMergedEdgeSrn);
}

void FtDsEdgeSrnModel::SplitEdgeSrn(const double relativePosition)
{
    // Edges with longer relative lengths after split have higher concentration
    ScaleSrnVariables(relativePosition);
}

// Declare identifier for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(FtDsEdgeSrnModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(FtDsEdgeSrnModel)
