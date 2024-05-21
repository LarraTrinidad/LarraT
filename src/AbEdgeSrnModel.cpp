#include "AbEdgeSrnModel.hpp"
#include "AbEdgeOdeSystem.hpp"

AbEdgeSrnModel::AbEdgeSrnModel(
    boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : AbstractOdeSrnModel(4, pOdeSolver)
{
    if (mpOdeSolver == boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
    {
#ifdef CHASTE_CVODE
        mpOdeSolver = CellCycleModelOdeSolver<AbEdgeSrnModel, CvodeAdaptor>::Instance();
        mpOdeSolver->Initialise();
        mpOdeSolver->SetMaxSteps(10000);
#else
        mpOdeSolver = CellCycleModelOdeSolver<AbEdgeSrnModel, RungeKutta4IvpOdeSolver>::Instance();
        mpOdeSolver->Initialise();
        SetDt(0.001);
#endif //CHASTE_CVODE
    }
    assert(mpOdeSolver->IsSetUp());
}

AbEdgeSrnModel::AbEdgeSrnModel(const AbEdgeSrnModel& rModel)
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
    SetOdeSystem(new AbEdgeOdeSystem(p_parent_system->rGetStateVariables()));
    for (unsigned i = 0; i < p_parent_system->GetNumberOfParameters(); ++i)
    {
        mpOdeSystem->SetParameter(i, p_parent_system->GetParameter(i));
    }
}

AbstractSrnModel* AbEdgeSrnModel::CreateSrnModel()
{
    return new AbEdgeSrnModel(*this);
}

void AbEdgeSrnModel::SimulateToCurrentTime()
{
    // Update information before running simulation
    UpdateAb();

    // Run the ODE simulation as needed
    AbstractOdeSrnModel::SimulateToCurrentTime();
}

void AbEdgeSrnModel::Initialise()
{
    AbstractOdeSrnModel::Initialise(new AbEdgeOdeSystem);
}

void AbEdgeSrnModel::InitialiseDaughterCell()
{
    assert(mpOdeSystem != nullptr);
    assert(mpCell != nullptr);

    // A new edge is initialised with zero concentrations
    ///\todo is this correct?
    mpOdeSystem->SetStateVariable("A", 1.0);
    mpOdeSystem->SetStateVariable("B", 1.0);
    mpOdeSystem->SetStateVariable("C", 0.0);
    mpOdeSystem->SetStateVariable("neigh_C", 0.0);

    mpOdeSystem->SetParameter("neighbour A", 1.0);
    mpOdeSystem->SetParameter("neighbour B", 1.0);
}

void AbEdgeSrnModel::UpdateAb()
{
    assert(mpOdeSystem != nullptr);
    assert(mpCell != nullptr);

    auto p_data = mpCell->GetCellEdgeData();

    double neigh_A = p_data->GetItem("neighbour A")[this->GetEdgeLocalIndex()];
    mpOdeSystem->SetParameter("neighbour A", neigh_A);

    double neigh_B = p_data->GetItem("neighbour B")[this->GetEdgeLocalIndex()];
    mpOdeSystem->SetParameter("neighbour B", neigh_B);
}

double AbEdgeSrnModel::GetA()
{
    assert(mpOdeSystem != nullptr);
    double A = mpOdeSystem->rGetStateVariables()[0];
    return A;
}

void AbEdgeSrnModel::SetA(double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[0] = value;
}

double AbEdgeSrnModel::GetB()
{
    assert(mpOdeSystem != nullptr);
    double B = mpOdeSystem->rGetStateVariables()[1];
    return B;
}

void AbEdgeSrnModel::SetB(double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[1] = value;
}

double AbEdgeSrnModel::GetC()
{
    assert(mpOdeSystem != nullptr);
    double C = mpOdeSystem->rGetStateVariables()[2];
    return C;
}

void AbEdgeSrnModel::SetC(double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[2] = value;
}

double AbEdgeSrnModel::GetNeighC()
{
    assert(mpOdeSystem != nullptr);
    double neigh_C = mpOdeSystem->rGetStateVariables()[3];
    return neigh_C;
}

void AbEdgeSrnModel::SetNeighC(double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[3] = value;
}

void AbEdgeSrnModel::OutputSrnModelParameters(out_stream& rParamsFile)
{
    AbstractOdeSrnModel::OutputSrnModelParameters(rParamsFile);
}

void AbEdgeSrnModel::AddSrnQuantities(
    AbstractSrnModel* pOtherSrn,
    const double scale)
{
    const double this_A = GetA();
    const double this_B = GetB();
    const double this_C = GetC();
    const double this_neigh_C = GetNeighC();

    auto p_other_srn = static_cast<AbEdgeSrnModel*>(pOtherSrn);
    const double other_A = p_other_srn->GetA();
    const double other_B = p_other_srn->GetB();
    const double other_C = p_other_srn->GetC();
    const double other_neigh_C = p_other_srn->GetNeighC();

    ///\todo This is unclear
    SetA(this_A + scale*other_A);
    SetB(this_B + scale*other_B);
    SetC(this_C + scale*other_C);
    SetNeighC(this_neigh_C + scale*other_neigh_C);
}

void AbEdgeSrnModel::AddShrunkEdgeSrn(AbstractSrnModel* pShrunkEdgeSrn)
{
    /*
     * Here we assume that one half of SRN quantities are endocytosed and the 
     * remaining half are split between neighbouring junctions. Hence we add 
     * 1/4 of SRN variables.
     */
    AddSrnQuantities(pShrunkEdgeSrn, 0.25); //0.25
}

void AbEdgeSrnModel::AddMergedEdgeSrn(AbstractSrnModel* pMergedEdgeSrn)
{
    // Add all srn variables to this edge srn
    AddSrnQuantities(pMergedEdgeSrn);
}

void AbEdgeSrnModel::SplitEdgeSrn(const double relativePosition)
{
    // Edges with longer relative lengths after split have higher concentration
    ScaleSrnVariables(relativePosition);
}

// Declare identifier for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(AbEdgeSrnModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(AbEdgeSrnModel)