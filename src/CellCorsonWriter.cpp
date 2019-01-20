//
// Created by gq on 1/15/19.
//

#include "CellCorsonWriter.hpp"
#include "AbstractCellPopulation.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellCorsonWriter<ELEMENT_DIM, SPACE_DIM>::CellCorsonWriter()
        : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("CellCorsonOdeParameter.dat")
{
    this->mVtkCellDataName = "Cell state variable u";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellCorsonWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    double u = pCell->GetCellData()->GetItem("corson cell state");
    return u;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellCorsonWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    // Output the location index corresponding to this cell
    unsigned location_index = pCellPopulation->GetLocationIndexUsingCell(pCell);
    *this->mpOutStream << location_index << " ";

    // Output this cell's ID
    unsigned cell_id = pCell->GetCellId();
    *this->mpOutStream << cell_id << " ";

    // Output the position of this cell's centre
    c_vector<double, SPACE_DIM> centre_location = pCellPopulation->GetLocationOfCellCentre(pCell);
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        *this->mpOutStream << centre_location[i] << " ";
    }

    double cell_state_u = pCell->GetCellData()->GetItem("corson cell state");
    *this->mpOutStream << cell_state_u << " ";


    double cell_signaling_s = pCell->GetCellData()->GetItem("corson signaling");
    *this->mpOutStream << cell_signaling_s << " ";
}

// Explicit instantiation
template class CellCorsonWriter<1,1>;
template class CellCorsonWriter<1,2>;
template class CellCorsonWriter<2,2>;
template class CellCorsonWriter<1,3>;
template class CellCorsonWriter<2,3>;
template class CellCorsonWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellCorsonWriter)
