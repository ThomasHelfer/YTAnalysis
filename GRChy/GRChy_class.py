""""
This is the birth of GRChy

"""

import h5py

class GRChy():

    def __init__(self, path):
        h5 = h5py.load(path)    # Check

        try:
            h5.....  (chombo_file)

        except Exception as e:
            # Error message


        num_levels = h5.....

        fields = h5.... attributes


    def load_data(self, field):

        return 0

    def load_slice(self, field):

        return 0 