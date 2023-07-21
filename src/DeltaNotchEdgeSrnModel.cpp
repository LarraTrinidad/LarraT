/*

Copyright (c) 2005-2021, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "DeltaNotchEdgeSrnModel.hpp"

DeltaNotchEdgeSrnModel::DeltaNotchEdgeSrnModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : AbstractOdeSrnModel(12, pOdeSolver)
{
    if (mpOdeSolver == boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
    {
#ifdef CHASTE_CVODE
        mpOdeSolver = CellCycleModelOdeSolver<DeltaNotchEdgeSrnModel, CvodeAdaptor>::Instance();
        mpOdeSolver->Initialise();
        mpOdeSolver->SetMaxSteps(10000);
#else
        mpOdeSolver = CellCycleModelOdeSolver<DeltaNotchEdgeSrnModel, RungeKutta4IvpOdeSolver>::Instance();
        mpOdeSolver->Initialise();
        SetDt(0.001);
#endif //CHASTE_CVODE
    }
    assert(mpOdeSolver->IsSetUp());
}

DeltaNotchEdgeSrnModel::DeltaNotchEdgeSrnModel(const DeltaNotchEdgeSrnModel& rModel)
    : AbstractOdeSrnModel(rModel)
{
    /*
     * Set each member variable of the new SRN model that inherits
     * its value from the parent.
     *
     * Note 1: some of the new SRN model's member variables
     * will already have been correctly initialized in its constructor.
     *
     * Note 2: one or more of the new SRN model's member variables
     * may be set/overwritten as soon as InitialiseDaughterCell() is called on
     * the new SRN model.
     *
     * Note 3: Only set the variables defined in this class. Variables defined
     * in parent classes will be defined there.
     */
    assert(rModel.GetOdeSystem());
    AbstractOdeSystem* p_parent_system(rModel.GetOdeSystem());
    SetOdeSystem(new DeltaNotchEdgeOdeSystem(p_parent_system->rGetStateVariables()));
    for (unsigned int i=0; i < p_parent_system->GetNumberOfParameters(); ++i)
        mpOdeSystem->SetParameter(i, p_parent_system->GetParameter(i));
}

AbstractSrnModel* DeltaNotchEdgeSrnModel::CreateSrnModel()
{
    return new DeltaNotchEdgeSrnModel(*this);
}

void DeltaNotchEdgeSrnModel::SimulateToCurrentTime()
{
    // Update information before running simulation
    UpdateDeltaNotch();
    // Run the ODE simulation as needed
    AbstractOdeSrnModel::SimulateToCurrentTime();
}

void DeltaNotchEdgeSrnModel::Initialise()
{
    AbstractOdeSrnModel::Initialise(new DeltaNotchEdgeOdeSystem);
}

void DeltaNotchEdgeSrnModel::InitialiseDaughterCell()
{
    assert(mpOdeSystem != nullptr);
    assert(mpCell != nullptr);


    //A new edge is initialised with zero concentrations
    mpOdeSystem->SetStateVariable("Notch",0.0);
    mpOdeSystem->SetStateVariable("Delta",0.0);
    mpOdeSystem->SetStateVariable("DsP",0.0);
    mpOdeSystem->SetStateVariable("FtP",0.0);
    mpOdeSystem->SetStateVariable("A",0.0);
    mpOdeSystem->SetStateVariable("B",0.0);
    mpOdeSystem->SetStateVariable("C",0.0);
    mpOdeSystem->SetStateVariable("D",0.0);
    mpOdeSystem->SetStateVariable("A_",0.0);
    mpOdeSystem->SetStateVariable("B_",0.0);
    mpOdeSystem->SetStateVariable("C_",0.0);
    mpOdeSystem->SetStateVariable("D_",0.0);
   

    mpOdeSystem->SetParameter("neighbour notch",0.0);
    mpOdeSystem->SetParameter("neighbour delta",0.0);
    mpOdeSystem->SetParameter("neighbour DsP",0.0);
    mpOdeSystem->SetParameter("neighbour FtP",0.0);

}


void DeltaNotchEdgeSrnModel::UpdateDeltaNotch()
{
    assert(mpOdeSystem != nullptr);
    assert(mpCell != nullptr);

    double neigh_delta
    = mpCell->GetCellEdgeData()->GetItem("neighbour delta")[this->GetEdgeLocalIndex()];
    mpOdeSystem->SetParameter("neighbour delta", neigh_delta);

   // double interior_delta = mpCell->GetCellData()->GetItem("interior delta");
   // mpOdeSystem->SetParameter("interior delta", interior_delta);

   // double interior_notch = mpCell->GetCellData()->GetItem("interior notch");
    //mpOdeSystem->SetParameter("interior notch", interior_notch);

}


double DeltaNotchEdgeSrnModel::GetNotch()
{
    assert(mpOdeSystem != nullptr);
    double notch = mpOdeSystem->rGetStateVariables()[0];
    return notch;
}

void DeltaNotchEdgeSrnModel::SetNotch(double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[0] = value;
}


double DeltaNotchEdgeSrnModel::GetDelta()
{
    assert(mpOdeSystem != nullptr);
    double delta = mpOdeSystem->rGetStateVariables()[1];
    return delta;
}

void DeltaNotchEdgeSrnModel::SetDelta(double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[1] = value;
}



double DeltaNotchEdgeSrnModel::GetDsP()
{
    assert(mpOdeSystem != nullptr);
    double DsP = mpOdeSystem->rGetStateVariables()[2];
    return DsP;
}

void DeltaNotchEdgeSrnModel::SetDsP(double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[2] = value;
}

double DeltaNotchEdgeSrnModel::GetFtP()
{
    assert(mpOdeSystem != nullptr);
    double FtP = mpOdeSystem->rGetStateVariables()[3];
    return FtP;
}

void DeltaNotchEdgeSrnModel::SetFtP(double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[3] = value;
}

double DeltaNotchEdgeSrnModel::GetA()
{
    assert(mpOdeSystem != nullptr);
    double A = mpOdeSystem->rGetStateVariables()[4];
    return A;
}

void DeltaNotchEdgeSrnModel::SetA(double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[4] = value;
}

double DeltaNotchEdgeSrnModel::GetB()
{
    assert(mpOdeSystem != nullptr);
    double B = mpOdeSystem->rGetStateVariables()[5];
    return B;
}

void DeltaNotchEdgeSrnModel::SetB(double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[5] = value;
}

double DeltaNotchEdgeSrnModel::GetC()
{
    assert(mpOdeSystem != nullptr);
    double C = mpOdeSystem->rGetStateVariables()[6];
    return C;
}

void DeltaNotchEdgeSrnModel::SetC(double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[6] = value;
}

double DeltaNotchEdgeSrnModel::GetD()
{
    assert(mpOdeSystem != nullptr);
    double D = mpOdeSystem->rGetStateVariables()[7];
    return D;
}

void DeltaNotchEdgeSrnModel::SetD(double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[7] = value;
}


double DeltaNotchEdgeSrnModel::GetA_()
{
    assert(mpOdeSystem != nullptr);
    double A_ = mpOdeSystem->rGetStateVariables()[8];
    return A_;
}

void DeltaNotchEdgeSrnModel::SetA_(double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[8] = value;
}

double DeltaNotchEdgeSrnModel::GetB_()
{
    assert(mpOdeSystem != nullptr);
    double B_ = mpOdeSystem->rGetStateVariables()[9];
    return B_;
}

void DeltaNotchEdgeSrnModel::SetB_(double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[9] = value;
}

double DeltaNotchEdgeSrnModel::GetC_()
{
    assert(mpOdeSystem != nullptr);
    double C_ = mpOdeSystem->rGetStateVariables()[10];
    return C_;
}

void DeltaNotchEdgeSrnModel::SetC_(double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[10] = value;
}

double DeltaNotchEdgeSrnModel::GetD_()
{
    assert(mpOdeSystem != nullptr);
    double D_ = mpOdeSystem->rGetStateVariables()[11];
    return D_;
}

void DeltaNotchEdgeSrnModel::SetD_(double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[11] = value;
}


double DeltaNotchEdgeSrnModel::GetNeighbouringDelta() const
{
    assert(mpOdeSystem != nullptr);
    return mpOdeSystem->GetParameter("neighbour delta");
}




double DeltaNotchEdgeSrnModel::GetNeighbouringDsP() const
{
    assert(mpOdeSystem != nullptr);
    return mpOdeSystem->GetParameter("neighbour DsP");
}

double DeltaNotchEdgeSrnModel::GetNeighbouringFtP() const
{
    assert(mpOdeSystem != nullptr);
    return mpOdeSystem->GetParameter("neighbour FtP");
}


double DeltaNotchEdgeSrnModel::GetInteriorDelta() const
{
    assert(mpOdeSystem != nullptr);
    return mpOdeSystem->GetParameter("interior delta");
}

double DeltaNotchEdgeSrnModel::GetInteriorNotch() const
{
    assert(mpOdeSystem != nullptr);
    return mpOdeSystem->GetParameter("interior notch");
}


void DeltaNotchEdgeSrnModel::OutputSrnModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    AbstractOdeSrnModel::OutputSrnModelParameters(rParamsFile);
}

void DeltaNotchEdgeSrnModel::AddSrnQuantities(AbstractSrnModel *p_other_srn,
                                              const double scale)
{
    auto other_srn = static_cast<DeltaNotchEdgeSrnModel*>(p_other_srn);
    const double other_delta = other_srn->GetDelta();
    const double other_notch = other_srn->GetNotch();
    const double other_DsP = other_srn->GetDsP();
    const double other_FtP = other_srn->GetFtP();
    const double other_A = other_srn->GetA();
    const double other_B = other_srn->GetB();
    const double other_C = other_srn->GetC();
    const double other_D = other_srn->GetD();
    const double other_A_ = other_srn->GetA_();
    const double other_B_ = other_srn->GetB_();
    const double other_C_ = other_srn->GetC_();
    const double other_D_ = other_srn->GetD_();


    const double this_delta = GetDelta();
    const double this_notch = GetNotch();
    const double this_DsP = GetDsP();
    const double this_FtP = GetFtP();
    const double this_A = GetA();
    const double this_B = GetB();
    const double this_C = GetC();
    const double this_D = GetD();
    const double this_A_ = GetA_();
    const double this_B_ = GetB_();
    const double this_C_ = GetC_();
    const double this_D_ = GetD_();

    SetDelta(this_delta+scale*other_delta);
    SetNotch(this_notch+scale*other_notch);
    SetDsP(this_DsP+scale*other_DsP);
    SetFtP(this_FtP+scale*other_FtP);
    SetA(this_A+scale*other_A);
    SetB(this_B+scale*other_B);
    SetC(this_C+scale*other_C);
    SetD(this_D+scale*other_D);
    SetA(this_A_+scale*other_A_);
    SetB(this_B_+scale*other_B_);
    SetC(this_C_+scale*other_C_);
    SetD(this_D_+scale*other_D_);
}


void DeltaNotchEdgeSrnModel::AddShrunkEdgeSrn(AbstractSrnModel *p_shrunk_edge_srn)
{
    // Here we assume that one half of srn quantities are endocytosed and the remaining
    // half are split between neighbouring junctions. Hence we add 1/4 of srn variables
    AddSrnQuantities(p_shrunk_edge_srn, 0.25);
}

void DeltaNotchEdgeSrnModel::AddMergedEdgeSrn(AbstractSrnModel* p_merged_edge_srn)
{
    // Add all srn variables to this edge srn
    AddSrnQuantities(p_merged_edge_srn);
}

void DeltaNotchEdgeSrnModel::SplitEdgeSrn(const double relative_position)
{
    //Edges with longer relative lengths after split have higher concentration
    ScaleSrnVariables(relative_position);
}


// Declare identifier for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(DeltaNotchEdgeSrnModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(DeltaNotchEdgeSrnModel)
