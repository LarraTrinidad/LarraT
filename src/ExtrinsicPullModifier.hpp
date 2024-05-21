#ifndef EXTRINSICPULLMODIFIER_HPP_
#define EXTRINSICPULLMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCellBasedSimulationModifier.hpp"

class ExtrinsicPullModifier : public AbstractCellBasedSimulationModifier<2,2>
{
private:
  
  /** Needed for serialization. */
  friend class boost::serialization::access;
  /**
   * Boost Serialization method for archiving/checkpointing.
   * Archives the object and its member variables.
   *
   * @param archive  The boost archive.
   * @param version  The current version of this class.
   */
  template<class Archive>
  void serialize(Archive & archive, const unsigned int version)
  {
    archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<2,2> >(*this);
    archive & mApplyExtrinsicPullToAllNodes;
    archive & mPinAnteriorMostCells;
    archive & mSpeed;
  }
  
  bool mApplyExtrinsicPullToAllNodes;
  bool mPinAnteriorMostCells;
  double mSpeed;
  
public:
  
  /**
   * Default constructor.
   */
  ExtrinsicPullModifier();
  
  /**
   * Destructor.
   */
  virtual ~ExtrinsicPullModifier();
  
  /**
   * Overridden UpdateAtEndOfTimeStep() method.
   *
   * @param rCellPopulation reference to the cell population
   */
  virtual void UpdateAtEndOfTimeStep(
      AbstractCellPopulation<2,2>& rCellPopulation) override;
  
  /**
   * Overridden SetupSolve() method.
   * 
   * @param rCellPopulation reference to the cell population
   * @param outputDirectory the output directory, relative to where Chaste 
   *     output is stored
   */
  virtual void SetupSolve(
      AbstractCellPopulation<2,2>& rCellPopulation,
      std::string outputDirectory) override;
  
  /**
   * @param applyExtrinsicPullToAllNodes whether to apply the extrinsic pull to 
   *     all nodes in the tissue
   */
  void ApplyExtrinsicPullToAllNodes(bool applyExtrinsicPullToAllNodes);
  
  /**
   * @param speed the speed of the extrinsic pull
   */
  void SetSpeed(double speed);
  
  /**
   * Overridden OutputSimulationModifierParameters() method.
   *
   * @param rParamsFile the file stream to which the parameters are output
   */
  void OutputSimulationModifierParameters(out_stream& rParamsFile) override;
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(ExtrinsicPullModifier)
  
#endif /*EXTRINSICPULLMODIFIER_HPP_*/