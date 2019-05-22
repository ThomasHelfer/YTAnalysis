"""
Mesh and Other Params
"""
# def base attributes
base_attrb = dict()
base_attrb['time'] = 0.0
base_attrb['iteration'] = 0
base_attrb['max_level'] = 0
base_attrb['num_components'] = len(component_names)
base_attrb['num_levels'] = 1
base_attrb['regrid_interval_0'] = 1
base_attrb['steps_since_regrid_0'] = 0
for comp,  name in enumerate(component_names):
    key = 'component_' + str(comp)
    tt = 'S' + str(len(name))
    base_attrb[key] = np.array(name, dtype=tt)


# def Chombo_global attributes
chombogloba_attrb = dict()
chombogloba_attrb['testReal'] = 0.0
chombogloba_attrb['SpaceDim'] = 3

# def level0 attributes
level_attrb = dict()
level_attrb['dt'] = 250.0
level_attrb['dx'] = float(L)/N
level_attrb['time'] = 0.0
level_attrb['is_periodic_0'] = 1
level_attrb['is_periodic_1'] = 1
level_attrb['is_periodic_2'] = 1
level_attrb['ref_ratio']= 2
level_attrb['tag_buffer_size'] = 3
prob_dom = (0, 0, 0, N-1, N-1, N-1)
prob_dt = np.dtype([('lo_i', '<i4'), ('lo_j', '<i4'), ('lo_k', '<i4'),
                    ('hi_i', '<i4'), ('hi_j', '<i4'), ('hi_k', '<i4')])
level_attrb['prob_domain'] = np.array(prob_dom, dtype=prob_dt)
boxes = np.array([(0, 0, 0, N-1, N-1, N-1)],
      dtype=[('lo_i', '<i4'), ('lo_j', '<i4'), ('lo_k', '<i4'), ('hi_i', '<i4'), ('hi_j', '<i4'), ('hi_k', '<i4')])


""""
CREATE HDF5
"""

if os.path.exists(filename):
    os.remove(filename)

h5file = h5.File(filename, 'w')  # New hdf5 file I want to create

# base attributes
for key in base_attrb.keys():
    h5file.attrs[key] = base_attrb[key]

# group: Chombo_global
chg = h5file.create_group('Chombo_global')
for key in chombogloba_attrb.keys():
    chg.attrs[key] = chombogloba_attrb[key]

# group: levels
l0 = h5file.create_group('level_0')
for key in level_attrb.keys():
    l0.attrs[key] = level_attrb[key]
sl0 = l0.create_group('data_attributes')
dadt = np.dtype([('intvecti', '<i4'), ('intvectj', '<i4'), ('intvectk', '<i4')])
sl0.attrs['ghost'] = np.array((3, 3, 3),  dtype=dadt)
sl0.attrs['outputGhost'] = np.array( (0, 0, 0),  dtype=dadt)
sl0.attrs['comps'] = base_attrb['num_components']
sl0.attrs['objectType'] = np.array('FArrayBox', dtype='S10')

# level datasets
dataset = np.zeros((base_attrb['num_components'], N, N, N))
for i, comp in enumerate(component_names):
    if comp in dset.keys():
        dataset[i] = dset[comp].T
fdset = []
for c in range(base_attrb['num_components']):
    fc = dataset[c].T.flatten()
    fdset.extend(fc)
fdset = np.array(fdset)

l0.create_dataset("Processors", data=np.array([0]))
l0.create_dataset("boxes",  data=boxes)
l0.create_dataset("data:offsets=0",  data=np.array([0, (base_attrb['num_components'])*N**3]))
l0.create_dataset("data:datatype=0",  data=fdset)

h5file.close()