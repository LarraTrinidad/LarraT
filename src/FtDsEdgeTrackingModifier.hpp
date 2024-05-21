#ifndef FTDSEDGERACKINGMODIFIER_HPP_
#define FTDSEDGERACKINGMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>


#include "AbstractCellBasedSimulationModifier.hpp"

template<unsigned DIM>
class FtDsEdgeTrackingModifier : public AbstractCellBasedSimulationModifier<DIM,DIM>
{
private:

    /**
     * Diffusion coefficient of unbound protein concentrations within the plasma 
     * membrane. Stored in this class, rather than in FtDsEdgeSrnModel or in 
     * FtDsEdgeOdeSystem, since it is only used in this class. Initialised to 
     * 0.03 in the constructor.
     */
    double mUnboundProteinDiffusionCoefficient;

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
        archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM,DIM> >(*this);
        archive & mUnboundProteinDiffusionCoefficient;
    }

public:

    /**
     * Default constructor.
     */
    FtDsEdgeTrackingModifier();

    /**
     * Destructor.
     */
    virtual ~FtDsEdgeTrackingModifier();

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * Specifies what to do in the simulation at the end of each time step.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Overridden SetupSolve() method.
     *
     * Specifies what to do in the simulation before the start of the time loop.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory);

    /**
     * Helper method.
     *
     * @param rCellPopulation reference to the cell population
     */
    void UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * @return mUnboundProteinDiffusionCoefficient.
     */
    double GetUnboundProteinDiffusionCoefficient();

    /**
     * Set mUnboundProteinDiffusionCoefficient.
     * 
     * @param unboundProteinDiffusionCoefficient the new value of 
     *     mUnboundProteinDiffusionCoefficient.
     */
    void SetUnboundProteinDiffusionCoefficient(
        double unboundProteinDiffusionCoefficient);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(FtDsEdgeTrackingModifier)

#endif //FTDSEDGERACKINGMODIFIER_HPP_
