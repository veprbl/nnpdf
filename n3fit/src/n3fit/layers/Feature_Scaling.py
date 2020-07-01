import numpy as np

from n3fit.backends import MetaLayer
from tensorflow.keras import backend as K
import tensorflow as tf

from validphys.fkparser import load_fktable
from validphys.loader import FallbackLoader as Loader

fk_dis_datasets=[
    'NMCPD_D',
    # 'NMCPD_P',
    'NMC',
    'SLACP',
    'SLACD',
    'BCDMSP',
    'BCDMSD',
    'CHORUSNU',
    'CHORUSNB',
    'NTVNUDMN',
    'NTVNBDMN',
    'HERACOMBNCEM',
    'HERACOMBNCEP460',
    'HERACOMBNCEP575',
    'HERACOMBNCEP820',
    'HERACOMBNCEP920',
    'HERACOMBCCEM',
    'HERACOMBCCEP',
    'HERAF2CHARM',
    'H1HERAF2B',
    'ZEUSHERAF2B',
    ]

class Feature_Scaling(MetaLayer):
    """
        Applies a Normalisation of the x-grid distribution.
    """

    def __init__(
        self, **kwargs
    ):
        #def load_fk_tables(datasets):
        #    l = Loader()
        #    fk_xgrids = np.array([])
        #    for fk_dataset in datasets:
        #        print(f'loading {fk_dataset}')
        #        fk = l.check_fktable(setname=fk_dataset, theoryID=52, cfac=[])
        #        res = load_fktable(fk)
        #        fk_xgrids = np.concatenate([fk_xgrids, res.xgrid])
        #    return fk_xgrids

        #fk_xgrids = load_fk_tables(fk_dis_datasets)
        #scaled_xgrids = fk_xgrids
        #scaled_xgrids = (
        #2 * scaled_xgrids
        #+ scaled_xgrids ** 0.4
        #+ 1.2 * scaled_xgrids ** 0.3
        #+ 0.7 * scaled_xgrids ** 0.2
        #+ 0.5 * scaled_xgrids ** 0.1
        #)
        #self.max_ = scaled_xgrids.max()
        #self.min_ = scaled_xgrids.min()

        super().__init__(**kwargs)

    def call(self, x_raw):

        self.max_ = 5.330364801340823
        self.min_ = 0.1794173872000228

        def log10(x):
            numerator = K.log(x)
            denominator = K.log( tf.constant(10, dtype=numerator.dtype))
            return numerator/denominator

        flattened_xgrids = (
        2 * x_raw
        + x_raw ** 0.4
        + 1.2 * x_raw ** 0.3
        + 0.7 * x_raw ** 0.2
        + 0.5 * x_raw ** 0.1
        )

        feature_range_min = -1
        feature_range_max = 1

        scaled_xgrids = (flattened_xgrids - self.min_) / (
            self.max_ - self.min_
        ) * (feature_range_max - feature_range_min) + feature_range_min

        return scaled_xgrids
