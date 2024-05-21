#include "FtDsEdgeTrackingModifier.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "CellSrnModel.hpp"
#include "FtDsEdgeSrnModel.hpp"
#include "Debug.hpp"

template<unsigned DIM>
FtDsEdgeTrackingModifier<DIM>::FtDsEdgeTrackingModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
      mUnboundProteinDiffusionCoefficient(0.03)
{
}

template<unsigned DIM>
FtDsEdgeTrackingModifier<DIM>::~FtDsEdgeTrackingModifier()
{
}

template<unsigned DIM>
void FtDsEdgeTrackingModifier<DIM>::UpdateAtEndOfTimeStep(
    AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Update the cell
    this->UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void FtDsEdgeTrackingModifier<DIM>::SetupSolve(
    AbstractCellPopulation<DIM,DIM>& rCellPopulation,
    std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void FtDsEdgeTrackingModifier<DIM>::UpdateCellData(
    AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{
    double D = mUnboundProteinDiffusionCoefficient;

    // For ease, store a static cast of the vertex-based cell population
    assert(dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation));
    auto p_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);

    // Iterate over cells
    for (auto cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // Get this cell's cell-level SRN model and number of edges
        assert(dynamic_cast<CellSrnModel*>(cell_iter->GetSrnModel()));
        auto p_cell_srn = static_cast<CellSrnModel*>(cell_iter->GetSrnModel());
        unsigned num_edges = p_cell_srn->GetNumEdgeSrn();

        // Cell's edge data
        std::vector<double> Ds_old;
        std::vector<double> Ft_old;
        std::vector<double> DsP_old;
        std::vector<double> FtP_old;
        std::vector<double> A_old;
        std::vector<double> B_old;
        std::vector<double> C_old;
        std::vector<double> D_old;
        std::vector<double> neighA_old;
        std::vector<double> neighB_old;
        std::vector<double> neighC_old;
        std::vector<double> neighD_old;

        std::vector<double> bound_Ft_old;
        std::vector<double> bound_FtP_old;
        std::vector<double> bound_Ds_old;
        std::vector<double> bound_DsP_old;
        std::vector<double> bound_totalFt_old;
        std::vector<double> bound_totalDs_old;

        std::vector<double> bound_totalneighFt_old;
        std::vector<double> bound_totalneighDs_old;

        std::vector<double> bound_neigh_Ft_old;
        std::vector<double> bound_neigh_FtP_old;
        std::vector<double> bound_neigh_Ds_old;
        std::vector<double> bound_neigh_DsP_old;

        std::vector<double> edge_lengths;

        for (unsigned edge_index = 0 ; edge_index  < num_edges; ++edge_index)
        {  
            // Store the current unbound protein concentrations on this edge
            auto p_edge_srn = boost::static_pointer_cast<FtDsEdgeSrnModel>(p_cell_srn->GetEdgeSrn(edge_index));
            Ds_old.push_back(p_edge_srn->GetDs());
            Ft_old.push_back(p_edge_srn->GetFt());
            DsP_old.push_back(p_edge_srn->GetDsP());
            FtP_old.push_back(p_edge_srn->GetFtP());
            A_old.push_back(p_edge_srn->GetA());
            B_old.push_back(p_edge_srn->GetB());
            C_old.push_back(p_edge_srn->GetC());
            D_old.push_back(p_edge_srn->GetD());
            neighA_old.push_back(p_edge_srn->GetNeighA());
            neighB_old.push_back(p_edge_srn->GetNeighB());
            neighC_old.push_back(p_edge_srn->GetNeighC());
            neighD_old.push_back(p_edge_srn->GetNeighD());

            bound_Ft_old.push_back((p_edge_srn->GetC()) + (p_edge_srn->GetD()));
            bound_FtP_old.push_back((p_edge_srn->GetA()) + (p_edge_srn->GetB()));
            bound_Ds_old.push_back((p_edge_srn->GetA()) + (p_edge_srn->GetC()));
            bound_DsP_old.push_back((p_edge_srn->GetB()) + (p_edge_srn->GetD()));
            bound_totalFt_old.push_back((p_edge_srn->GetA()) + (p_edge_srn->GetB())+(p_edge_srn->GetC()) + (p_edge_srn->GetD()));
            bound_totalDs_old.push_back((p_edge_srn->GetA()) + (p_edge_srn->GetB())+(p_edge_srn->GetC()) + (p_edge_srn->GetD()));

            bound_totalneighFt_old.push_back((p_edge_srn->GetNeighA()) + (p_edge_srn->GetNeighB())+(p_edge_srn->GetNeighC()) + (p_edge_srn->GetNeighD()));
            bound_totalneighDs_old.push_back((p_edge_srn->GetNeighA()) + (p_edge_srn->GetNeighB())+(p_edge_srn->GetNeighC()) + (p_edge_srn->GetNeighD()));

            bound_neigh_Ft_old.push_back((p_edge_srn->GetNeighC()) + (p_edge_srn->GetNeighD()));
            bound_neigh_FtP_old.push_back((p_edge_srn->GetNeighA()) + (p_edge_srn->GetNeighB()));
            bound_neigh_Ds_old.push_back((p_edge_srn->GetNeighA()) + (p_edge_srn->GetNeighC()));
            bound_neigh_DsP_old.push_back((p_edge_srn->GetNeighB()) + (p_edge_srn->GetNeighD()));

            // Store this edge's length
            auto p_element = p_population->GetElementCorrespondingToCell(*cell_iter);
            double edge_length = p_element->GetEdge(edge_index)->rGetLength();
            edge_lengths.push_back(edge_length);
        }

        /*
         * Update unbound protein concentrations based on a linear diffusive 
         * flux between neighbouring edges.
         */
        std::vector<double> Ds_new(num_edges);
        std::vector<double> Ft_new(num_edges);
        std::vector<double> DsP_new(num_edges);
        std::vector<double> FtP_new(num_edges);
        std::vector<double> A_new(num_edges);
        std::vector<double> B_new(num_edges);
        std::vector<double> C_new(num_edges);
        std::vector<double> D_new(num_edges);
        std::vector<double> neighA_new(num_edges);
        std::vector<double> neighB_new(num_edges);
        std::vector<double> neighC_new(num_edges);
        std::vector<double> neighD_new(num_edges);

        std::vector<double> bound_Ft_new(num_edges);
        std::vector<double> bound_FtP_new(num_edges);
        std::vector<double> bound_Ds_new(num_edges);
        std::vector<double> bound_DsP_new(num_edges);
        std::vector<double> bound_totalFt_new(num_edges);
        std::vector<double> bound_totalDs_new(num_edges);
        std::vector<double> bound_totalneighFt_new(num_edges);
        std::vector<double> bound_totalneighDs_new(num_edges);

        std::vector<double> bound_neigh_Ft_new(num_edges);
        std::vector<double> bound_neigh_FtP_new(num_edges);
        std::vector<double> bound_neigh_Ds_new(num_edges);
        std::vector<double> bound_neigh_DsP_new(num_edges);

           for (unsigned edge_index = 0 ; edge_index  < num_edges; ++edge_index)
        {
           auto p_edge_srn = boost::static_pointer_cast<FtDsEdgeSrnModel>(p_cell_srn->GetEdgeSrn(edge_index));
           unsigned prev_index = (edge_index == 0) ? num_edges - 1 : edge_index - 1;
           unsigned next_index = (edge_index == num_edges - 1) ? 0 : edge_index + 1;

           ///\todo consider validity of diffusive flux expression
           double dx = 0.5 * (edge_lengths[prev_index] + edge_lengths[edge_index]);
           double dt = SimulationTime::Instance()->GetTimeStep();
           double time_elapsed = SimulationTime::Instance()->GetTimeStepsElapsed();
           /* This if statement must be used as UpdateCellData is called for setup solve
             and at the end of each time step. Thus if no time has elapsed diffusion should
             not be occuring */
             if (time_elapsed == 0){
              Ds_new[edge_index] = Ds_old[edge_index];
              Ft_new[edge_index] = Ft_old[edge_index];
              DsP_new[edge_index] = DsP_old[edge_index];
              FtP_new[edge_index] = FtP_old[edge_index];
              A_new[edge_index] = A_old[edge_index];
              B_new[edge_index] = B_old[edge_index];
              C_new[edge_index] = C_old[edge_index];
              D_new[edge_index] = D_old[edge_index];
              neighA_new[edge_index] = neighA_old[edge_index];
              neighB_new[edge_index] = neighB_old[edge_index];
              neighC_new[edge_index] = neighC_old[edge_index];
              neighD_new[edge_index] = neighD_old[edge_index];
              
              bound_Ft_new[edge_index] = C_old[edge_index] + D_old[edge_index];
              bound_FtP_new[edge_index] = A_old[edge_index] + B_old[edge_index];
              bound_Ds_new[edge_index] = A_old[edge_index] + C_old[edge_index];
              bound_DsP_new[edge_index] = B_old[edge_index] + D_old[edge_index];
              bound_totalFt_new[edge_index] = A_old[edge_index] + B_old[edge_index] + C_old[edge_index] + D_old[edge_index];
              bound_totalDs_new[edge_index] = A_old[edge_index] + B_old[edge_index] + C_old[edge_index] + D_old[edge_index];

              bound_totalneighFt_new[edge_index] = neighA_old[edge_index] + neighB_old[edge_index] + neighC_old[edge_index] + neighD_old[edge_index];
              bound_totalneighDs_new[edge_index] = neighA_old[edge_index] + neighB_old[edge_index] + neighC_old[edge_index] + neighD_old[edge_index];

              bound_neigh_Ft_new[edge_index] = neighC_old[edge_index] + neighD_old[edge_index]; 
              bound_neigh_FtP_new[edge_index] = neighA_old[edge_index] + neighB_old[edge_index]; 
              bound_neigh_Ds_new[edge_index] = neighA_old[edge_index] + neighC_old[edge_index]; 
              bound_neigh_DsP_new[edge_index] = neighB_old[edge_index] + neighD_old[edge_index];

             } else {
            Ds_new[edge_index] = Ds_old[edge_index] + D*(Ds_old[prev_index] - 2*Ds_old[edge_index] + Ds_old[next_index])*dt/(dx*dx);
            Ft_new[edge_index] = Ft_old[edge_index] + D*(Ft_old[prev_index] - 2*Ft_old[edge_index] + Ft_old[next_index])*dt/(dx*dx);
            DsP_new[edge_index] = DsP_old[edge_index] + D*(DsP_old[prev_index] - 2*DsP_old[edge_index] + DsP_old[next_index])*dt/(dx*dx);
            FtP_new[edge_index] = FtP_old[edge_index] + D*(FtP_old[prev_index] - 2*FtP_old[edge_index] + FtP_old[next_index])*dt/(dx*dx);
            A_new[edge_index] = A_old[edge_index];
            B_new[edge_index] = B_old[edge_index];
            C_new[edge_index] = C_old[edge_index];
            D_new[edge_index] = D_old[edge_index];
            neighA_new[edge_index] = neighA_old[edge_index];
            neighB_new[edge_index] = neighB_old[edge_index];
            neighC_new[edge_index] = neighC_old[edge_index];
            neighD_new[edge_index] = neighD_old[edge_index];

            bound_Ft_new[edge_index] = C_old[edge_index] + D_old[edge_index];
            bound_FtP_new[edge_index] = A_old[edge_index] + B_old[edge_index];
            bound_Ds_new[edge_index] = A_old[edge_index] + C_old[edge_index];
            bound_DsP_new[edge_index] = B_old[edge_index] + D_old[edge_index];
            bound_totalFt_new[edge_index] = A_old[edge_index] + B_old[edge_index] + C_old[edge_index] + D_old[edge_index];
            bound_totalDs_new[edge_index] = A_old[edge_index] + B_old[edge_index] + C_old[edge_index] + D_old[edge_index];

            bound_totalneighFt_new[edge_index] = neighA_old[edge_index] + neighB_old[edge_index] + neighC_old[edge_index] + neighD_old[edge_index];
            bound_totalneighDs_new[edge_index] = neighA_old[edge_index] + neighB_old[edge_index] + neighC_old[edge_index] + neighD_old[edge_index];

            bound_neigh_Ft_new[edge_index] = neighC_old[edge_index] + neighD_old[edge_index]; 
            bound_neigh_FtP_new[edge_index] = neighA_old[edge_index] + neighB_old[edge_index]; 
            bound_neigh_Ds_new[edge_index] = neighA_old[edge_index] + neighC_old[edge_index]; 
            bound_neigh_DsP_new[edge_index] = neighB_old[edge_index] + neighD_old[edge_index]; 

             }
             
             p_edge_srn->SetDs(Ds_new[edge_index]);
             p_edge_srn->SetFt(Ft_new[edge_index]);
             p_edge_srn->SetDsP(DsP_new[edge_index]);
             p_edge_srn->SetFtP(FtP_new[edge_index]);
             p_edge_srn->SetA(A_new[edge_index]);
             p_edge_srn->SetB(B_new[edge_index]);
             p_edge_srn->SetC(C_new[edge_index]);
             p_edge_srn->SetD(D_new[edge_index]);
             p_edge_srn->SetNeighA(neighA_new[edge_index]);
             p_edge_srn->SetNeighB(neighB_new[edge_index]);
             p_edge_srn->SetNeighC(neighC_new[edge_index]);
             p_edge_srn->SetNeighD(neighD_new[edge_index]);


        }

        // Note: state variables must be in the same order as in FtDsOdeSystem
        cell_iter->GetCellEdgeData()->SetItem("edge Ds", Ds_new);
        cell_iter->GetCellEdgeData()->SetItem("edge Ft", Ft_new);
        cell_iter->GetCellEdgeData()->SetItem("edge DsP", DsP_new);
        cell_iter->GetCellEdgeData()->SetItem("edge FtP", FtP_new);
        cell_iter->GetCellEdgeData()->SetItem("edge A", A_new);
        cell_iter->GetCellEdgeData()->SetItem("edge B", B_new);
        cell_iter->GetCellEdgeData()->SetItem("edge C", C_new);
        cell_iter->GetCellEdgeData()->SetItem("edge D", D_new);
        cell_iter->GetCellEdgeData()->SetItem("neighbour A", neighA_new);
        cell_iter->GetCellEdgeData()->SetItem("neighbour B", neighB_new);
        cell_iter->GetCellEdgeData()->SetItem("neighbour C", neighC_new);
        cell_iter->GetCellEdgeData()->SetItem("neighbour D", neighD_new);
        
        cell_iter->GetCellEdgeData()->SetItem("bound Ft", bound_Ft_new);
        cell_iter->GetCellEdgeData()->SetItem("bound FtP", bound_FtP_new);
        cell_iter->GetCellEdgeData()->SetItem("bound Ds", bound_Ds_new);
        cell_iter->GetCellEdgeData()->SetItem("bound DsP", bound_DsP_new);
        cell_iter->GetCellEdgeData()->SetItem("bound total Ft", bound_totalFt_new);
        cell_iter->GetCellEdgeData()->SetItem("bound total Ds", bound_totalDs_new);

        cell_iter->GetCellEdgeData()->SetItem("bound total neighbour Ft", bound_totalneighFt_new);
        cell_iter->GetCellEdgeData()->SetItem("bound total neighbour Ds", bound_totalneighDs_new);

        cell_iter->GetCellEdgeData()->SetItem("bound neighbour Ft", bound_neigh_Ft_new);
        cell_iter->GetCellEdgeData()->SetItem("bound neighbour FtP", bound_neigh_FtP_new);
        cell_iter->GetCellEdgeData()->SetItem("bound neighbour Ds", bound_neigh_Ds_new);
        cell_iter->GetCellEdgeData()->SetItem("bound neighbour DsP", bound_neigh_DsP_new);

    }

    // After the edge data is filled, fill the edge neighbour data

    // Iterate over cells again
    for (auto cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // Get this cell's cell-level SRN model and number of edges
        auto p_cell_srn = static_cast<CellSrnModel*>(cell_iter->GetSrnModel());
        unsigned num_edges = p_cell_srn->GetNumEdgeSrn();

        std::vector<double> neigh_mean_Ds(num_edges);
        std::vector<double> neigh_mean_Ft(num_edges);
        std::vector<double> neigh_mean_DsP(num_edges);
        std::vector<double> neigh_mean_FtP(num_edges);

        for (unsigned edge_index = 0; edge_index < num_edges; ++edge_index)
        {
            // Get neighbouring cell's values
            ///\todo is this correct?
            auto elem_neighbours = p_population->GetNeighbouringEdgeIndices(*cell_iter, edge_index);
            for (auto neighbour : elem_neighbours)
            {
                auto p_cell = p_population->GetCellUsingLocationIndex(neighbour.first);
                auto p_data = p_cell->GetCellEdgeData();
                std::vector<double> neighbour_Ds_vec = p_data->GetItem("edge Ds");
                std::vector<double> neighbour_Ft_vec = p_data->GetItem("edge Ft");
                std::vector<double> neighbour_DsP_vec = p_data->GetItem("edge DsP");
                std::vector<double> neighbour_FtP_vec = p_data->GetItem("edge FtP");
                
                neigh_mean_Ds[edge_index] += neighbour_Ds_vec[neighbour.second] / elem_neighbours.size();
                neigh_mean_Ft[edge_index] += neighbour_Ft_vec[neighbour.second] / elem_neighbours.size();
                neigh_mean_DsP[edge_index] += neighbour_DsP_vec[neighbour.second] / elem_neighbours.size();
                neigh_mean_FtP[edge_index] += neighbour_FtP_vec[neighbour.second] / elem_neighbours.size();
                
            }
        }

        cell_iter->GetCellEdgeData()->SetItem("neighbour Ds", neigh_mean_Ds);
        cell_iter->GetCellEdgeData()->SetItem("neighbour Ft", neigh_mean_Ft);
        cell_iter->GetCellEdgeData()->SetItem("neighbour DsP", neigh_mean_DsP);
        cell_iter->GetCellEdgeData()->SetItem("neighbour FtP", neigh_mean_FtP);
          }
}

template<unsigned DIM>
double FtDsEdgeTrackingModifier<DIM>::GetUnboundProteinDiffusionCoefficient()
{
    return mUnboundProteinDiffusionCoefficient;
}

template<unsigned DIM>
void FtDsEdgeTrackingModifier<DIM>::SetUnboundProteinDiffusionCoefficient(
    double unboundProteinDiffusionCoefficient)
{
    mUnboundProteinDiffusionCoefficient = unboundProteinDiffusionCoefficient;
}

template<unsigned DIM>
void FtDsEdgeTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class FtDsEdgeTrackingModifier<1>;
template class FtDsEdgeTrackingModifier<2>;
template class FtDsEdgeTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"

EXPORT_TEMPLATE_CLASS_SAME_DIMS(FtDsEdgeTrackingModifier)