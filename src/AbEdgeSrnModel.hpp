#ifndef ABEDGESRNMODEL_HPP_
#define ABEDGESRNMODEL_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "FtDsEdgeOdeSystem.hpp"
#include "AbstractOdeSrnModel.hpp"

class AbEdgeSrnModel : public AbstractOdeSrnModel
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the SRN model and member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractOdeSrnModel>(*this);
    }

protected:

    /**
     * Protected copy-constructor for use by CreateSrnModel(). The only way for 
     * external code to create a copy of a SRN model is by calling that method, 
     * to ensure that a model of the correct subclass is created. This 
     * copy-constructor helps subclasses to ensure that all member variables are 
     * correctly copied when this happens.
     *
     * This method is called by child classes to set member variables for a 
     * daughter cell upon cell division. Note that the parent SRN model will 
     * have had ResetForDivision() called just before CreateSrnModel() is 
     * called, so performing an exact copy of the parent is suitable behaviour. 
     * Any daughter-cell-specific initialisation can be done in 
     * InitialiseDaughterCell().
     *
     * @param rModel  the SRN model to copy.
     */
    AbEdgeSrnModel(const AbEdgeSrnModel& rModel);

public:

    /**
     * Default constructor calls base class.
     *
     * @param pOdeSolver An optional pointer to a cell-cycle model ODE solver 
     *     object (allows the use of different ODE solvers)
     */
    AbEdgeSrnModel(
        boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver = boost::shared_ptr<AbstractCellCycleModelOdeSolver>());

    /**
     * Overridden builder method to create new copies of this SRN model.
     *
     * @return a copy of the current SRN model.
     */
    virtual AbstractSrnModel* CreateSrnModel() override;

    /**
     * Overridden Initialise() method to initialise the SRN model at the start 
     * of a simulation.
     */
    virtual void Initialise() override;

    /**
     * This method is called when a new edge is created (e.g. after cell 
     * division or T1 swap)
     */
    virtual void InitialiseDaughterCell() override;

    /**
     * Overridden SimulateToTime() method for custom behaviour.
     */
    virtual void SimulateToCurrentTime() override;

    /**
     * Update the protein concentrations at neighbouring edge sensed by this 
     * edge; that is, fetch neighbour values from CellEdgeData object, storing 
     * the sensed information, into this model.
     */
    void UpdateAb();

    /** @return the unbound concentration of A at this edge. */
    double GetA();

    /**
     * Set the unbound concentration of A at this edge.
     * 
     * @param value
     */
    void SetA(double value);

    /** @return the unbound concentration of B at this edge. */
    double GetB();

    /**
     * Set the unbound concentration of B at this edge.
     * 
     * @param value
     */
    void SetB(double value);

    /** @return the complex concentration C at this edge. */
    double GetC();

    /**
     * Set the complex concentration C at this edge.
     * 
     * @param value
     */
    void SetC(double value);

    /**
     * @return the current complex C concentration at the neighbouring cell's 
     *     edge.
     */
    double GetNeighC();

    /**
     * Set the complex C concentration at the neighbouring cell's edge.
     * 
     * @param value
     */
    void SetNeighC(double value);

    /**
     * Overridden OutputSrnModelParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputSrnModelParameters(out_stream& rParamsFile) override;

    /**
     * Overridden AddSrnQuantities() method.
     *
     * @param pOtherSrn
     * @param scale
     */
    virtual void AddSrnQuantities(AbstractSrnModel* pOtherSrn,
                                  const double scale = 1.0) override;

    /**
     * Overridden AddShrunkEdgeSrn() method.
     * 
     * @param pShrunkEdgeSrn
     */
    virtual void AddShrunkEdgeSrn(AbstractSrnModel* pShrunkEdgeSrn) override;

    /**
     * Overridden AddMergedEdgeSrn() method.
     * 
     * @param pMergedEdgeSrn
     */
    virtual void AddMergedEdgeSrn(AbstractSrnModel* pMergedEdgeSrn) override;

    /**
     * Overridden SplitEdgeSrn() method.
     * 
     * @param relativePosition
     */
    virtual void SplitEdgeSrn(const double relativePosition) override;
};

typedef boost::shared_ptr<AbEdgeSrnModel> AbEdgeSrnModelPtr;

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(AbEdgeSrnModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(AbEdgeSrnModel)

#endif  /* ABEDGESRNMODEL_HPP_ */
