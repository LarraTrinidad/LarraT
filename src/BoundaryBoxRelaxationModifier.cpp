#include "BoundaryBoxRelaxationModifier.hpp"
#include "Toroidal2dVertexMeshWithMutableSize.hpp"
#include "StressTensor.hpp"

///\todo Had to turn off LHS relaxation because sdk bounday condition can't account for it.

BoundaryBoxRelaxationModifier::BoundaryBoxRelaxationModifier()
    : AbstractCellBasedSimulationModifier<2>(),
      mStiffness(1000),
      mpForce(nullptr)
{
}

BoundaryBoxRelaxationModifier::~BoundaryBoxRelaxationModifier()
{
}

void BoundaryBoxRelaxationModifier::UpdateAtEndOfTimeStep(AbstractCellPopulation<2,2>& rCellPopulation)
{
    double dt = SimulationTime::Instance()->GetTimeStep();

    // Pointer to mesh
    AbstractMesh<2, 2>& r_mesh = rCellPopulation.rGetMesh();
    Toroidal2dVertexMeshWithMutableSize* p_mesh = static_cast<Toroidal2dVertexMeshWithMutableSize*>(&r_mesh);
    
    // Coords of box
    double current_x_lower = p_mesh->GetBoxCoords(0);
    double current_x_upper = p_mesh->GetBoxCoords(1);
    double current_y_lower = p_mesh->GetBoxCoords(2);
    double current_y_upper = p_mesh->GetBoxCoords(3);
    
    // Stress tensor
    c_matrix<double, 2,2> stress_tensor_2d = GetTissueStressTensor(rCellPopulation, mpForce.get());
    
    // Force in x
    double x_force = - stress_tensor_2d(0, 0);
    
    // Force in y
    double y_force = - stress_tensor_2d(1, 1);
    
    // Displacement of box is force/stiffness
    double deltaX = x_force / mStiffness;
    double deltaY = y_force / mStiffness;
    
    // Centroid of box
    double box_centroid_x = 0.5 * (current_x_lower + current_x_upper);
    double box_centroid_y = 0.5 * (current_y_lower + current_y_upper);

    // Angles to corners of box. Used to see which side boundary vertices are on
    double anglex0y0 = std::atan2(current_y_lower - box_centroid_y,
                                  current_x_lower - box_centroid_x);
    double anglex1y0 = std::atan2(current_y_lower - box_centroid_y,
                                  current_x_upper - box_centroid_x);
    
    // Find the boundary nodes
    std::set<unsigned> boundary_nodes = p_mesh->GetBoundaryNodes();
    
    double current_height = current_y_upper - current_y_lower;
    double new_height = current_height*(1-dt*deltaY);
    double height_diff = -new_height + current_height;
    
    // If it was a boundary node, check where it is and pull it appropriately
    for (auto n_index : boundary_nodes)
    {
      // Get the node
      Node<2>* p_node = p_mesh->GetNode(n_index);
      double node_x = p_node->rGetLocation()[0];
      double node_y = p_node->rGetLocation()[1];
      
      // Calculate angle to node
      double angleNode = std::atan2(node_y - box_centroid_y,
                                    node_x - box_centroid_x);

      /*
       * Stretch/compress at the correct boundaries. Note, don't need to stretch 
       * upper and RHS because they are deermined by height and width of box, so 
       * just stretch the box after.
       */

      // Bottom
      if (anglex0y0 <= angleNode && anglex1y0 >= angleNode)
      {
        p_node->rGetModifiableLocation()[1] += 0.5*height_diff;
      }
    }

    // Reset the size of the box
    p_mesh->SetBoxCoords(1, current_x_upper*(1 - dt*deltaX));
    p_mesh->SetBoxCoords(2, current_y_lower + 0.5*height_diff);
    p_mesh->SetBoxCoords(3, current_y_upper - 0.5*height_diff);
}

void BoundaryBoxRelaxationModifier::SetupSolve(
    AbstractCellPopulation<2, 2>& rCellPopulation,
    std::string outputDirectory)
{
}

void BoundaryBoxRelaxationModifier::SetStiffness(double stiffness)
{
  mStiffness = stiffness;
}

void BoundaryBoxRelaxationModifier::SetForce(
    boost::shared_ptr<FarhadifarForce<2> > pForce)
{
  mpForce = pForce;
}

void BoundaryBoxRelaxationModifier::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
  *rParamsFile << "\t\t\t<Stiffness>" << mStiffness << "</Stiffness>\n";
  AbstractCellBasedSimulationModifier<2>::OutputSimulationModifierParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(BoundaryBoxRelaxationModifier)