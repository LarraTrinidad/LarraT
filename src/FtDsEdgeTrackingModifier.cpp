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
void FtDsEdgeTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Update the cell
    this->UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void FtDsEdgeTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void FtDsEdgeTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Recovers each cell's edge levels proteins, and those of its neighbor's
    // Then saves them
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        auto p_cell_edge_model = static_cast<CellSrnModel*>(cell_iter->GetSrnModel());

        /* Cells edge data */
        std::vector<double> notch_vec;
        std::vector<double> delta_vec;
        std::vector<double> DsP_vec;
        std::vector<double> FtP_vec;
        std::vector<double> A_vec;
        std::vector<double> B_vec;
        std::vector<double> C_vec;
        std::vector<double> D_vec;
        std::vector<double> A__vec;
        std::vector<double> B__vec;
        std::vector<double> C__vec;
        std::vector<double> D__vec;


        std::vector<double> x_vec;
        std::vector<double> xm1_vec;
        std::vector<double> xp1_vec;
        for (unsigned i = 0 ; i  < p_cell_edge_model->GetNumEdgeSrn(); i++)
        {
            boost::shared_ptr<FtDsEdgeSrnModel> p_model
            = boost::static_pointer_cast<FtDsEdgeSrnModel>(p_cell_edge_model->GetEdgeSrn(i));
            double this_notch = p_model->GetNotch();
            double this_delta = p_model->GetDelta();
            double this_DsP = p_model->GetDsP();
            double this_FtP = p_model->GetFtP();
            double this_A = p_model->GetA();
            double this_B = p_model->GetB();
            double this_C = p_model->GetC();
            double this_D = p_model->GetD();
            double this_A_ = p_model->GetA_();
            double this_B_ = p_model->GetB_();
            double this_C_ = p_model->GetC_();
            double this_D_ = p_model->GetD_();


            if ( i == 0){
            double this_xm1 = rCellPopulation.rGetMesh().GetDistanceBetweenNodes(p_cell_edge_model->GetNumEdgeSrn() -1, i);
            double this_x = rCellPopulation.rGetMesh().GetDistanceBetweenNodes(i, i + 1);
            double this_xp1 = rCellPopulation.rGetMesh().GetDistanceBetweenNodes(i + 1, i + 2);
            x_vec.push_back(this_x);
            xm1_vec.push_back(this_xm1);
            xp1_vec.push_back(this_xp1);
            }else if (i > 0 && i < p_cell_edge_model->GetNumEdgeSrn()){
            double this_xm1 = rCellPopulation.rGetMesh().GetDistanceBetweenNodes(i - 1, i);
            double this_x = rCellPopulation.rGetMesh().GetDistanceBetweenNodes(i, i + 1);
            double this_xp1 = rCellPopulation.rGetMesh().GetDistanceBetweenNodes(i + 1, i + 2);
            x_vec.push_back(this_x);
            xm1_vec.push_back(this_xm1);
            xp1_vec.push_back(this_xp1);
            }
  
            delta_vec.push_back(this_notch);
            notch_vec.push_back(this_delta);
            DsP_vec.push_back(this_DsP);
            FtP_vec.push_back(this_FtP);
            A_vec.push_back(this_A);
            B_vec.push_back(this_B);
            C_vec.push_back(this_C);
            D_vec.push_back(this_D);
            A__vec.push_back(this_A_);
            B__vec.push_back(this_B_);
            C__vec.push_back(this_C_);
            D__vec.push_back(this_D_);
        }
        
        for (unsigned i = 0 ; i  < p_cell_edge_model->GetNumEdgeSrn(); i++)
        {
           double D = 0.03; //pow(10, -2);
           boost::shared_ptr<FtDsEdgeSrnModel> p_model
           = boost::static_pointer_cast<FtDsEdgeSrnModel>(p_cell_edge_model->GetEdgeSrn(i));

           double x1 = (xm1_vec[i] + x_vec[i])/2.0;
           //double x2 = (xp1_vec[i] + x_vec[i])/2.0;
           //double x3 = (x1 + x2)/2.0;

           double ym1 = i - 1;
           double y1 = i;
           double yp1 = i + 1;
           
           if ( i == 0){
             ym1 = p_cell_edge_model->GetNumEdgeSrn() -1;
           } else if ( i == p_cell_edge_model->GetNumEdgeSrn() -1){
             yp1 = 0;
           } else {
             ym1 = i - 1;
             yp1 = i + 1;
           } 

           double this_notch = notch_vec[y1] + D*((notch_vec[ym1] - 2*notch_vec[y1] + notch_vec[yp1])/(x1*x1));  
           notch_vec[i] = this_notch;

           double this_delta = delta_vec[y1] + D*((delta_vec[ym1] - 2*delta_vec[y1] + delta_vec[yp1])/(x1*x1));  
           delta_vec[i] = this_delta;

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

           double this_A_ = A__vec[y1];
           A__vec[i] = this_A_;

           double this_B_ = A__vec[y1];
           B__vec[i] = this_B_;

           double this_C_ = C__vec[y1];
           C__vec[i] = this_C_;

           double this__D = D__vec[y1];
           D__vec[i] = this__D;
          //boundfrizzled_vec[i] = combined1_vec[i] + combined3_vec[i] + combined4_vec[i];
          //boundstrabismus_vec[i] = combined2_vec[i] + combinedm3_vec[i] + combinedm4_vec[i];

        }
        // Note that the state variables must be in the same order as listed in FtDsOdeSystem
        cell_iter->GetCellEdgeData()->SetItem("edge notch", notch_vec);
        cell_iter->GetCellEdgeData()->SetItem("edge delta", delta_vec);
        cell_iter->GetCellEdgeData()->SetItem("edge DsP", DsP_vec);
        cell_iter->GetCellEdgeData()->SetItem("edge FtP", FtP_vec);
        cell_iter->GetCellEdgeData()->SetItem("edge A", A_vec);
        cell_iter->GetCellEdgeData()->SetItem("edge B", B_vec);
        cell_iter->GetCellEdgeData()->SetItem("edge C", C_vec);
        cell_iter->GetCellEdgeData()->SetItem("edge D", D_vec);
        cell_iter->GetCellEdgeData()->SetItem("edge A_", A__vec);
        cell_iter->GetCellEdgeData()->SetItem("edge B_", B__vec);
        cell_iter->GetCellEdgeData()->SetItem("edge C_", C__vec);
        cell_iter->GetCellEdgeData()->SetItem("edge D_", D__vec);


    }

    //After the edge data is filled, fill the edge neighbour data
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
            cell_iter != rCellPopulation.End();
            ++cell_iter)
    {
        auto p_cell_edge_model = static_cast<CellSrnModel*>(cell_iter->GetSrnModel());
        const unsigned int n_cell_edges = p_cell_edge_model->GetNumEdgeSrn();
        std::vector<double> neigh_mean_delta(n_cell_edges);
        std::vector<double> neigh_mean_notch(n_cell_edges);
        std::vector<double> neigh_mean_DsP(n_cell_edges);
        std::vector<double> neigh_mean_FtP(n_cell_edges);
        std::vector<double> neigh_mean_A(n_cell_edges);
        std::vector<double> neigh_mean_B(n_cell_edges);
        std::vector<double> neigh_mean_C(n_cell_edges);
        std::vector<double> neigh_mean_D(n_cell_edges);
        std::vector<double> neigh_mean_A_(n_cell_edges);
        std::vector<double> neigh_mean_B_(n_cell_edges);
        std::vector<double> neigh_mean_C_(n_cell_edges);
        std::vector<double> neigh_mean_D_(n_cell_edges);


        std::vector<double> notch_vec;
        std::vector<double> delta_vec;
        std::vector<double> DsP_vec;
        std::vector<double> FtP_vec;
        std::vector<double> A_vec;
        std::vector<double> B_vec;
        std::vector<double> C_vec;
        std::vector<double> D_vec;
        std::vector<double> A__vec;
        std::vector<double> B__vec;
        std::vector<double> C__vec;
        std::vector<double> D__vec;
        //std::vector<double> frizzled_vec;

        
        for (unsigned i = 0 ; i  < p_cell_edge_model->GetNumEdgeSrn(); i++)
        {
            boost::shared_ptr<FtDsEdgeSrnModel> p_model
            = boost::static_pointer_cast<FtDsEdgeSrnModel>(p_cell_edge_model->GetEdgeSrn(i));
            double this_notch = p_model->GetNotch();
            double this_delta = p_model->GetDelta();
            double this_DsP = p_model->GetDsP();
            double this_FtP = p_model->GetFtP();
            double this_A = p_model->GetA();
            double this_B = p_model->GetB();
            double this_C = p_model->GetC();
            double this_D = p_model->GetD();
            double this_A_ = p_model->GetA_();
            double this_B_ = p_model->GetB_();
            double this_C_ = p_model->GetC_();
            double this_D_ = p_model->GetD_();


            notch_vec.push_back(this_notch);
            delta_vec.push_back(this_delta);
            DsP_vec.push_back(this_DsP);
            FtP_vec.push_back(this_FtP);
            A_vec.push_back(this_A);
            B_vec.push_back(this_B);
            C_vec.push_back(this_C);
            D_vec.push_back(this_D);
            A__vec.push_back(this_A_);
            B__vec.push_back(this_B_);
            C__vec.push_back(this_C_);
            D__vec.push_back(this_D_);
        }
         
        for (unsigned int i=0; i<n_cell_edges; ++i)
        {      
            double this_notch = notch_vec[i];
            double this_delta = delta_vec[i];
            double this_DsP = DsP_vec[i];
            double this_FtP = FtP_vec[i];
            double this_A = A_vec[i];
            double this_B = B_vec[i];
            double this_C = C_vec[i];
            double this_D = D_vec[i];
            double this_A_ = A__vec[i];
            double this_B_ = B__vec[i];
            double this_C_ = C__vec[i];
            double this_D_ = D__vec[i];


            //Get neighbouring cell's values of delta on this
            auto elemNeighbours = rCellPopulation.GetNeighbouringEdgeIndices(*cell_iter, i);
            double mean_delta = 0;
            double mean_notch = 0;
            double mean_DsP = 0;
            double mean_FtP = 0;
            double mean_A = 0;
            double mean_B = 0;
            double mean_C = 0;
            double mean_D = 0;
            double mean_A_ = 0;
            double mean_B_ = 0;
            double mean_C_ = 0;
            double mean_D_ = 0;



            for (auto neighbourIndex: elemNeighbours)
            {
                auto neighbourCell = rCellPopulation.GetCellUsingLocationIndex(neighbourIndex.first);
                std::vector<double> neighbour_delta_vec = neighbourCell->GetCellEdgeData()->GetItem("edge delta");
                std::vector<double> neighbour_notch_vec = neighbourCell->GetCellEdgeData()->GetItem("edge notch");
                std::vector<double> neighbour_DsP_vec = neighbourCell->GetCellEdgeData()->GetItem("edge DsP");
                std::vector<double> neighbour_FtP_vec = neighbourCell->GetCellEdgeData()->GetItem("edge FtP");
                std::vector<double> neighbour_A_vec = neighbourCell->GetCellEdgeData()->GetItem("edge A");
                std::vector<double> neighbour_B_vec = neighbourCell->GetCellEdgeData()->GetItem("edge B");
                std::vector<double> neighbour_C_vec = neighbourCell->GetCellEdgeData()->GetItem("edge C");
                std::vector<double> neighbour_D_vec = neighbourCell->GetCellEdgeData()->GetItem("edge D");
                std::vector<double> neighbour_A__vec = neighbourCell->GetCellEdgeData()->GetItem("edge A_");
                std::vector<double> neighbour_B__vec = neighbourCell->GetCellEdgeData()->GetItem("edge B_");
                std::vector<double> neighbour_C__vec = neighbourCell->GetCellEdgeData()->GetItem("edge C_");
                std::vector<double> neighbour_D__vec = neighbourCell->GetCellEdgeData()->GetItem("edge D_");
                
                mean_delta += neighbour_delta_vec[neighbourIndex.second];
                mean_notch += neighbour_notch_vec[neighbourIndex.second];
                mean_DsP += neighbour_DsP_vec[neighbourIndex.second];
                mean_FtP += neighbour_FtP_vec[neighbourIndex.second];
                mean_A += neighbour_A_vec[neighbourIndex.second];
                mean_B += neighbour_B_vec[neighbourIndex.second];
                mean_C += neighbour_C_vec[neighbourIndex.second];
                mean_D += neighbour_D_vec[neighbourIndex.second];
                mean_A_ += neighbour_A__vec[neighbourIndex.second];
                mean_B_ += neighbour_B__vec[neighbourIndex.second];
                mean_C_ += neighbour_C__vec[neighbourIndex.second];
                mean_D_ += neighbour_D__vec[neighbourIndex.second];
           }
            if (elemNeighbours.size()>0){
                mean_notch = mean_notch;
                mean_delta = mean_delta;
                mean_DsP = mean_DsP;
                mean_FtP = mean_FtP;
                mean_A = mean_A;
                mean_B = mean_B;
                mean_C = mean_C;
                mean_D = mean_D;
                mean_A_ = mean_A_;
                mean_B_ = mean_B_;
                mean_C_ = mean_C_;
                mean_D_ = mean_D_;


            //double account_notch = ((this_notch)+mean_notch) ;
            neigh_mean_notch[i] = 0*(this_notch)+ mean_notch;
            neigh_mean_delta[i] = 0*(this_delta)+ mean_delta;
            neigh_mean_DsP[i] = 0*(this_DsP)+ mean_DsP;
            neigh_mean_FtP[i] = 0*(this_FtP)+ mean_FtP;
            neigh_mean_A[i] = 0*(this_A)+ mean_A;
            neigh_mean_B[i] = 0*(this_B)+ mean_B;
            neigh_mean_C[i] = 0*(this_C)+ mean_C;
            neigh_mean_D[i] = 0*(this_D)+ mean_D;
            neigh_mean_A_[i] = 0*(this_A_)+ mean_A_;
            neigh_mean_B_[i] = 0*(this_B_)+ mean_B_;
            neigh_mean_C_[i] = 0*(this_C_)+ mean_C_;
            neigh_mean_D_[i] = 0*(this_D_)+ mean_D_;

            }
        }
        cell_iter->GetCellEdgeData()->SetItem("neighbour notch", neigh_mean_notch);
        cell_iter->GetCellEdgeData()->SetItem("neighbour delta", neigh_mean_delta);
        cell_iter->GetCellEdgeData()->SetItem("neighbour DsP", neigh_mean_DsP);
        cell_iter->GetCellEdgeData()->SetItem("neighbour FtP", neigh_mean_FtP);
        cell_iter->GetCellEdgeData()->SetItem("neighbour A", neigh_mean_A);
        cell_iter->GetCellEdgeData()->SetItem("neighbour B", neigh_mean_B);
        cell_iter->GetCellEdgeData()->SetItem("neighbour C", neigh_mean_C);
        cell_iter->GetCellEdgeData()->SetItem("neighbour D", neigh_mean_D);
        cell_iter->GetCellEdgeData()->SetItem("neighbour A_", neigh_mean_A_);
        cell_iter->GetCellEdgeData()->SetItem("neighbour B_", neigh_mean_B_);
        cell_iter->GetCellEdgeData()->SetItem("neighbour C_", neigh_mean_C_);
        cell_iter->GetCellEdgeData()->SetItem("neighbour D_", neigh_mean_D_);

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