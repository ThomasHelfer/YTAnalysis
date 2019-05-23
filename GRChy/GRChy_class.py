""""
This is the birth of GRChy

"""

import h5py

class GRChy():

    def __init__(self, path):
        h5 = h5py.load(path)    # Check

        #TODO
        # try:
        #     h5.....  (chombo_file)
        #
        # except Exception as e:
        #     # Error message

        attrs = np.array(h5.attrs, dtype=str)
        max_level = h5.attrs['max_level']
        fields = attrs[ 'component_' in attrs]



    def load_data(self, field):

        return 0

    def load_slice(self, field):

        return 0 