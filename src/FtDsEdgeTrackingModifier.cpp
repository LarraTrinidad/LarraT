#include "FtDsEdgeTrackingModifier.hpp"
#include "CellSrnModel.hpp"
#include "FtDsEdgeSrnModel.hpp"

template<unsigned DIM>
FtDsEdgeTrackingModifier<DIM>::FtDsEdgeTrackingModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
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
    AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    /*
     * Recovers each cell's edge levels proteins, and those of its neighbour's, 
     * then saves them.
     */
    for (auto cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        auto p_cell_srn_model = static_cast<CellSrnModel*>(cell_iter->GetSrnModel());

        // Cell's edge data
        std::vector<double> Ds_vec;
        std::vector<double> Ft_vec;
        std::vector<double> DsP_vec;
        std::vector<double> FtP_vec;
        std::vector<double> A_vec;
        std::vector<double> B_vec;
        std::vector<double> C_vec;
        std::vector<double> D_vec;
        std::vector<double> neigh_A_vec;
        std::vector<double> neigh_B_vec;
        std::vector<double> neigh_C_vec;
        std::vector<double> neigh_D_vec;

        std::vector<double> x_vec;
        std::vector<double> xm1_vec;
        std::vector<double> xp1_vec;
        for (unsigned i = 0 ; i  < p_cell_srn_model->GetNumEdgeSrn(); ++i)
        {
            auto p_model = boost::static_pointer_cast<FtDsEdgeSrnModel>(p_cell_srn_model->GetEdgeSrn(i));
            double this_Ds = p_model->GetDs();
            double this_Ft = p_model->GetFt();
            double this_DsP = p_model->GetDsP();
            double this_FtP = p_model->GetFtP();
            double this_A = p_model->GetA();
            double this_B = p_model->GetB();
            double this_C = p_model->GetC();
            double this_D = p_model->GetD();
            double this_neigh_A = p_model->GetNeighA();
            double this_neigh_B = p_model->GetNeighB();
            double this_neigh_C = p_model->GetNeighC();
            double this_neigh_D = p_model->GetNeighD();

            if (i == 0)
            {
                double this_xm1 = rCellPopulation.rGetMesh().GetDistanceBetweenNodes(p_cell_srn_model->GetNumEdgeSrn() -1, i);
                double this_x = rCellPopulation.rGetMesh().GetDistanceBetweenNodes(i, i + 1);
                double this_xp1 = rCellPopulation.rGetMesh().GetDistanceBetweenNodes(i + 1, i + 2);
                x_vec.push_back(this_x);
                xm1_vec.push_back(this_xm1);
                xp1_vec.push_back(this_xp1);
            }
            else if (i > 0 && i < p_cell_srn_model->GetNumEdgeSrn())
            {
                double this_xm1 = rCellPopulation.rGetMesh().GetDistanceBetweenNodes(i - 1, i);
                double this_x = rCellPopulation.rGetMesh().GetDistanceBetweenNodes(i, i + 1);
                double this_xp1 = rCellPopulation.rGetMesh().GetDistanceBetweenNodes(i + 1, i + 2);
                x_vec.push_back(this_x);
                xm1_vec.push_back(this_xm1);
                xp1_vec.push_back(this_xp1);
            }
  
            Ds_vec.push_back(this_Ds);
            Ft_vec.push_back(this_Ft);
            DsP_vec.push_back(this_DsP);
            FtP_vec.push_back(this_FtP);
            A_vec.push_back(this_A);
            B_vec.push_back(this_B);
            C_vec.push_back(this_C);
            D_vec.push_back(this_D);
            neigh_A_vec.push_back(this_neigh_A);
            neigh_B_vec.push_back(this_neigh_B);
            neigh_C_vec.push_back(this_neigh_C);
            neigh_D_vec.push_back(this_neigh_D);
        }
        
        for (unsigned i = 0 ; i  < p_cell_srn_model->GetNumEdgeSrn(); ++i)
        {
           double D = 0.03;
           auto p_model = boost::static_pointer_cast<FtDsEdgeSrnModel>(p_cell_srn_model->GetEdgeSrn(i));

           double x1 = 0.5 * (xm1_vec[i] + x_vec[i]);

           double ym1 = i - 1;
           double y1 = i;
           double yp1 = i + 1;
           
           if (i == 0)
           {
               ym1 = p_cell_srn_model->GetNumEdgeSrn() - 1;
           }
           else if (i == p_cell_srn_model->GetNumEdgeSrn() - 1)
           {
               yp1 = 0;
           }
           else
           {
                ym1 = i - 1;
                yp1 = i + 1;
           } 

           double this_Ds = Ds_vec[y1] + D*((Ds_vec[ym1] - 2*Ds_vec[y1] + Ds_vec[yp1])/(x1*x1));  
           Ds_vec[i] = this_Ds;

           double this_Ft = Ft_vec[y1] + D*((Ft_vec[ym1] - 2*Ft_vec[y1] + Ft_vec[yp1])/(x1*x1));  
           Ft_vec[i] = this_Ft;

           double this_DsP = DsP_vec[y1] + D*((DsP_vec[ym1] - 2*DsP_vec[y1] + DsP_vec[yp1])/(x1*x1));  
           DsP_vec[i] = this_DsP;

           double this_FtP = FtP_vec[y1] + D*((FtP_vec[ym1] - 2*FtP_vec[y1] + FtP_vec[yp1])/(x1*x1));  
           FtP_vec[i] = this_FtP;

           double this_A = A_vec[y1];
           A_vec[i] = this_A;

           double this_B = A_vec[y1];
           B_vec[i] = this_B;

           double this_C = C_vec[y1];
           C_vec[i] = this_C;

           double this_D = D_vec[y1];
           D_vec[i] = this_D;

           double this_neigh_A = neigh_A_vec[y1];
           neigh_A_vec[i] = this_neigh_A;

           double this_neigh_B = neigh_B_vec[y1];
           neigh_B_vec[i] = this_neigh_B;

           double this_neigh_C = neigh_C_vec[y1];
           neigh_C_vec[i] = this_neigh_C;

           double this_neigh_D = neigh_D_vec[y1];
           neigh_D_vec[i] = this_neigh_D;
        }

        // Note that the state variables must be in the same order as listed in FtDsOdeSystem
        cell_iter->GetCellEdgeData()->SetItem("edge Ds", Ds_vec);
        cell_iter->GetCellEdgeData()->SetItem("edge Ft", Ft_vec);
        cell_iter->GetCellEdgeData()->SetItem("edge DsP", DsP_vec);
        cell_iter->GetCellEdgeData()->SetItem("edge FtP", FtP_vec);
        cell_iter->GetCellEdgeData()->SetItem("edge A", A_vec);
        cell_iter->GetCellEdgeData()->SetItem("edge B", B_vec);
        cell_iter->GetCellEdgeData()->SetItem("edge C", C_vec);
        cell_iter->GetCellEdgeData()->SetItem("edge D", D_vec);
        cell_iter->GetCellEdgeData()->SetItem("edge neigh A", neigh_A_vec);
        cell_iter->GetCellEdgeData()->SetItem("edge neigh B", neigh_B_vec);
        cell_iter->GetCellEdgeData()->SetItem("edge neigh C", neigh_C_vec);
        cell_iter->GetCellEdgeData()->SetItem("edge neigh D", neigh_D_vec);
    }

    // After the edge data is filled, fill the edge neighbour data
    for (auto cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        auto p_cell_srn_model = static_cast<CellSrnModel*>(cell_iter->GetSrnModel());

        const unsigned int num_edges = p_cell_srn_model->GetNumEdgeSrn();
        std::vector<double> neigh_mean_Ds(num_edges);
        std::vector<double> neigh_mean_Ft(num_edges);
        std::vector<double> neigh_mean_DsP(num_edges);
        std::vector<double> neigh_mean_FtP(num_edges);

        std::vector<double> Ds_vec;
        std::vector<double> Ft_vec;
        std::vector<double> DsP_vec;
        std::vector<double> FtP_vec;
        for (unsigned i = 0 ; i  < p_cell_srn_model->GetNumEdgeSrn(); ++i)
        {
            auto p_model = boost::static_pointer_cast<FtDsEdgeSrnModel>(p_cell_srn_model->GetEdgeSrn(i));
            Ds_vec.push_back(p_model->GetDs());
            Ft_vec.push_back(p_model->GetFt());
            DsP_vec.push_back(p_model->GetDsP());
            FtP_vec.push_back(p_model->GetFtP());
        }
         
        for (unsigned i = 0; i < num_edges; ++i)
        {
            // Get neighbouring cell's values
            double mean_Ds = 0;
            double mean_Ft = 0;
            double mean_DsP = 0;
            double mean_FtP = 0;
            auto elem_neighbours = rCellPopulation.GetNeighbouringEdgeIndices(*cell_iter, i);
            for (auto neighbour : elem_neighbours)
            {
                auto p_cell = rCellPopulation.GetCellUsingLocationIndex(neighbour.first);
                auto p_data = p_cell->GetCellEdgeData();
                std::vector<double> neighbour_Ds_vec = p_data->GetItem("edge Ds");
                std::vector<double> neighbour_Ft_vec = p_data->GetItem("edge Ft");
                std::vector<double> neighbour_DsP_vec = p_data->GetItem("edge DsP");
                std::vector<double> neighbour_FtP_vec = p_data->GetItem("edge FtP");
                
                mean_Ds += neighbour_Ds_vec[neighbour.second];
                mean_Ft += neighbour_Ft_vec[neighbour.second];
                mean_DsP += neighbour_DsP_vec[neighbour.second];
                mean_FtP += neighbour_FtP_vec[neighbour.second];
            }

            ///\todo this seems unclear
            if (elem_neighbours.size() > 0)
            {
                neigh_mean_Ds[i] = mean_Ds;
                neigh_mean_Ft[i] = mean_Ft;
                neigh_mean_DsP[i] = mean_DsP;
                neigh_mean_FtP[i] = mean_FtP;
            }
        }
        cell_iter->GetCellEdgeData()->SetItem("neighbour Ds", neigh_mean_Ds);
        cell_iter->GetCellEdgeData()->SetItem("neighbour Ft", neigh_mean_Ft);
        cell_iter->GetCellEdgeData()->SetItem("neighbour DsP", neigh_mean_DsP);
        cell_iter->GetCellEdgeData()->SetItem("neighbour FtP", neigh_mean_FtP);
    }
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