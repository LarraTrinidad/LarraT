#include "ExtrinsicPullModifier.hpp"

ExtrinsicPullModifier::ExtrinsicPullModifier()
    : AbstractCellBasedSimulationModifier<2>(),
      mApplyExtrinsicPullToAllNodes(true),
      mSpeed(1.0) 
{
}

ExtrinsicPullModifier::~ExtrinsicPullModifier()
{
}

void ExtrinsicPullModifier::UpdateAtEndOfTimeStep(
    AbstractCellPopulation<2,2>& rCellPopulation)
{
    double epsilon = 0.8;
    
    double dt = SimulationTime::Instance()->GetTimeStep();
    unsigned num_nodes = rCellPopulation.GetNumNodes();
    ChasteCuboid<2> bounds = rCellPopulation.rGetMesh().CalculateBoundingBox();
    double x_min = bounds.rGetLowerCorner()[0];
    double x_max = bounds.rGetUpperCorner()[0];
    
    if (mApplyExtrinsicPullToAllNodes)
    {
        // Pull on all nodes, with a constant strain rate
        double width = x_max - x_min;
        for (unsigned node_index = 0; node_index < num_nodes; ++node_index)
        {
            Node<2>* p_node = rCellPopulation.GetNode(node_index);
            double speed = mSpeed * (p_node->rGetLocation()[0] - x_min) / width;

            // Respect the SidekickBoundaryCondition...
            if (p_node->rGetLocation()[0] > x_min + epsilon)
            {
                p_node->rGetModifiableLocation()[0] -= speed*dt;
            }
        }
    }
    else
    {
        // Get the centroid of the boundary nodes
        // Define variables for centroid position
        double centroidX = 0;
        double centroidY = 0;
        double counter = 0;

        // Sum the X,Y coords for all points
        for (unsigned node_index = 0; node_index < num_nodes; ++node_index)
        {
            Node<2>* p_node = rCellPopulation.GetNode(node_index);
            
            std::set<unsigned> neighbourNodes = rCellPopulation.GetNeighbouringNodeIndices(node_index);
            
            bool is_boundary_node = false;
            if (neighbourNodes.size() < 3)
            {
               is_boundary_node = true;
            }
            else
            {
                for (auto neighbour : neighbourNodes)
                {
                    std::set<unsigned> neighbourNodes = rCellPopulation.GetNeighbouringNodeIndices(neighbour);
                    
                    if (neighbourNodes.size() < 3)
                    {
                        is_boundary_node = true;
                    }
                }
            }

            if (is_boundary_node)
            {
                centroidX += p_node->rGetLocation()[0];
                centroidY += p_node->rGetLocation()[1];
                counter += 1;
            }
        }

        // Divide by total to get centroid.
        centroidX /= counter;
        centroidY /= counter;
        
        // Make a vector of node indices:
        std::vector< int > boundaryNodes;
        for (unsigned node_index=0; node_index<num_nodes; node_index++)
        {
            // Get the node at this index
            Node<2>* p_node = rCellPopulation.GetNode(node_index);
            
            // Initialise as not a boundary node
            p_node->SetAsBoundaryNode(true); ///\todo should this be false?
            bool is_boundary_node = false;

            std::set<unsigned> neighbourNodes = rCellPopulation.GetNeighbouringNodeIndices(node_index);
            if (neighbourNodes.size() < 3)
            {
                // Update it as a boundary node
                is_boundary_node = true;
                p_node->SetAsBoundaryNode(true);
            }
            else
            {
                for (auto neighbour : neighbourNodes)
                {
                    // Node<2>* neighbour = rCellPopulation.(neighbourNodes[idx]);
                    std::set<unsigned> neighbourNodes = rCellPopulation.GetNeighbouringNodeIndices(neighbour);
                    
                    if (neighbourNodes.size() < 3)
                    {
                        // Update it as a boundary node
                        is_boundary_node = true;
                        p_node->SetAsBoundaryNode(true);
                    }
                }
            }

          if (is_boundary_node)
          {
              // Get the X and Y coords
              double currentX = p_node->rGetLocation()[0];
              double currentY = p_node->rGetLocation()[1];
              
              // If this is the first boundary node, just add it to the list
              if (boundaryNodes.size() == 0)
              {
                  boundaryNodes.push_back(node_index);
              }
              else
              {
                  /*
                   * Otherwise, compare it to the current nodes; put it in front 
                   * of the first node that makes a larger angle, relative to 
                   * centroid.
                   */
                  auto i = 0u; // counter
                  
                  bool done = false; // Breaks when we have placed the vertex
                  while (not done)
                  {
                      // If we have reached the end of the list, break
                      if (i == boundaryNodes.size())
                      {
                        boundaryNodes.push_back(node_index);
                        done = true;
                      }
                      
                      // Get the x and y coords from the next node in
                      // boundaryNodes to compare against.
                      Node<2>* p_node = rCellPopulation.GetNode(boundaryNodes[i]);

                      // Get the coords of the node at the current index in
                      // the ordered list of boundary nodes.
                      double previousX = p_node->rGetLocation()[0];
                      double previousY = p_node->rGetLocation()[1];
                      
                      // Calculate the angle that the two nodes make relative
                      // to the centroid
                      double prevAngle;
                      double dY = previousY - centroidY;
                      double dX = previousX - centroidX;
                      prevAngle = std::atan2( dY, dX );
                      double currentAngle;
                      dX = currentX - centroidX;
                      dY = currentY - centroidY;
                      currentAngle = std::atan2(dY, dX);
                      // If the angle of new node is smaller, it comes first
                      // so append it at this location.
                      if (currentAngle < prevAngle && not done)
                      {
                          boundaryNodes.insert(boundaryNodes.begin() + i, node_index);
                          done = true;
                      }
                      
                      // update count
                      i += 1;
                  }
                }
            }
        }
      
        // Get the top and bottom rightmost corners of bounding box
        double lower_x = bounds.rGetLowerCorner()[0];
        double lower_y = bounds.rGetLowerCorner()[1]; 
        
        double upperX =  lower_x; // bounds.rGetUpperCorner()[0];
        double upper_y = bounds.rGetUpperCorner()[1];
        
        // Make variables to hold indices of nodes closest to bounding box.
        unsigned upper_node_index = 0;
        unsigned lower_node_index = 0;

        // Euclid distance variables
        double dist;
        double upper_dist = 10000000.; // Initialise as huge
        double lower_dist = 10000000.; // Initialise as huge
        double temp_dist;
        double x_diff;
        double y_diff;
        for (unsigned i = 0; i < boundaryNodes.size(); ++i)
        {
            // Get the node pointer
            Node<2>* p_node = rCellPopulation.GetNode(boundaryNodes[i]);
            
            // Calculate distance to upper corner
            x_diff = upperX - p_node->rGetLocation()[0];
            y_diff = upper_y - p_node->rGetLocation()[1];
            dist = x_diff*x_diff + y_diff*y_diff;
            temp_dist = std::sqrt(dist);

            // If smaller, mark this as the closest node
            if (temp_dist < upper_dist)
            {
                upper_dist = temp_dist;
                upper_node_index = i;
            }
            
            // calc dist to lower corner
            x_diff = lower_x - p_node->rGetLocation()[0];
            y_diff = lower_y - p_node->rGetLocation()[1];
            dist = x_diff*x_diff + y_diff*y_diff;
            temp_dist = std::sqrt(dist);

            // If smaller, mark this as the closest node
            if (temp_dist < lower_dist)
            {
                lower_dist = temp_dist;
                lower_node_index = i;
            }
        }
        
        ///\todo Need to change this to leftmost nodes???
        for (unsigned i = 0; i < boundaryNodes.size(); ++i)
        {
            // If it lies between the top and bottom:
            if (i >= lower_node_index && i <= upper_node_index)
            {
                Node<2>* p_node = rCellPopulation.GetNode(boundaryNodes[i]);
                p_node->rGetModifiableLocation()[0] -= mSpeed*dt; // was += but changed to -= to move to the left
            }
            else if (i > upper_node_index)
            {
                break;
            }
        }
    }
}

void ExtrinsicPullModifier::SetupSolve(
    AbstractCellPopulation<2,2>& rCellPopulation,
    std::string outputDirectory)
{
}

void ExtrinsicPullModifier::ApplyExtrinsicPullToAllNodes(
    bool applyExtrinsicPullToAllNodes)
{
  mApplyExtrinsicPullToAllNodes = applyExtrinsicPullToAllNodes;
}

void ExtrinsicPullModifier::SetSpeed(double speed)
{
  mSpeed = speed;
}

void ExtrinsicPullModifier::OutputSimulationModifierParameters(
    out_stream& rParamsFile)
{
  *rParamsFile << "\t\t\t<ApplyExtrinsicPullToAllNodes>" << mApplyExtrinsicPullToAllNodes << "</ApplyExtrinsicPullToAllNodes>\n";
  *rParamsFile << "\t\t\t<Speed>" << mSpeed << "</Speed>\n";
  AbstractCellBasedSimulationModifier<2>::OutputSimulationModifierParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(ExtrinsicPullModifier)