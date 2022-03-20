"""
This module provides a functionality to split the nodes in the
output layer of a Model.
"""

import numpy as np

from n3fit.backends import MetaLayer
from n3fit.backends import MetaModel
from n3fit.backends import operations as op


class SplitPDFLayers(MetaLayer):
    """Custom `Metalayer` that splits the output of a given Model
    into several parts. This enables us to split the output PDFs
    into different nucleus.

    Parameters
    ----------
        indices: list(int)
            indices to be extracted
        size: int
            size of the expected pdf model
    """

    def __init__(self, indices, size, **kwargs):
        self._indices = indices
        mask = np.zeros(size, dtype=bool)
        mask[indices] = True
        self._mask = mask
        super().__init__(**kwargs)

    def call(self, pdf):
        # Transpose in order to apply boolean mask
        pdfT = op.transpose(pdf)
        pdf_masked = op.boolean_mask(pdfT, self._mask)
        # Revert transpose of the masked PDF
        return op.transpose(pdf_masked)


def split_pdfmodels(pdf_model, map_pdfs, nfl=14):
    """Takes a trained PDF model and split the output layer into a
    sub-proton and/or nuclear PDFs, each with 14 nodes/flavours.
    
    Parameters
    ----------
        pdf_model: n3fit.MetaModel
            trained model with all the outputs
        ma_pdfs: list(dict)
            list of dictionary containing information on the
            nuclear targets

    Returns
    -------
        list(n3fit.MetaModel):
            containg all the PDF metamodel for each nuclear/proton
    """
    # Fetch number of outputs
    nopt = nfl * len(map_pdfs)
    # Construct the split lists
    split = [[i + j for j in range(nfl)] for i in range(0, nopt, nfl)]
    # Get the input and apply the model as a layer
    full_input, pdf_layer = pdf_model.apply_as_layer(pdf_model.input)
    list_pdf_models = []
    for slice in split:
        # Extract the output PDFs for a given slice
        splitted_layer = SplitPDFLayers(slice, nopt)(pdf_layer)
        # Generate the Model for the sliced PDFs
        splitted_model = MetaModel(full_input, splitted_layer)
        list_pdf_models.append(splitted_model)
    return list_pdf_models
