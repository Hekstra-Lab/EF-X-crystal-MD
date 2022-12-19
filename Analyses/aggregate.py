import reciprocalspaceship as rs
from mdtools.utils import *

n_chains = 4
n_replicates = 5
last_epoch = 150

dataset = rs.read_mtz(f'processed_data/0/{last_epoch}/EF_10MV_cm_partial_{last_epoch}_chainwise_zero_subtraj_0_avg.mtz')
n_reflections = dataset.shape[0]

# need to do this so for phenix refinemenet, a bit of a hack...
def reformat_dataset(dataset, sigma=0.05):
    dataset['FOBS'] = dataset['FMODEL']
    dataset['SIGFOBS'] = (sigma*dataset['FOBS']).astype("Q")
    dataset = dataset[['FOBS', 'SIGFOBS']]
    
for phase_ in ['pos', 'neg', 'zero']:
    complex_reflections = np.zeros(n_reflections, dtype='complex128')
    for chain in range(n_chains):
        if phase_ == 'pos':
            if chain == 0 or chain == 2: 
                phase = 'pos'
            else:
                phase = 'neg'
        elif phase_ == 'neg':
            if chain == 0 or chain == 2: 
                phase = 'neg'
            else:
                phase = 'pos'
        else:
            phase = 'zero'
        for replicate in range(n_replicates):
            dataset = rs.read_mtz(f'processed_data/{replicate}/{last_epoch}/EF_10MV_cm_partial_{last_epoch}_chainwise_{phase}_subtraj_{chain}_avg.mtz')
            complex_reflections += mtz_to_cartesian_arr(dataset)
    complex_reflections /= (n_chains * n_replicates)
    dataset[:] = cartesian_arr_to_polar(complex_reflections)
    dataset.infer_mtz_dtypes(inplace = True)
    reformat_dataset(dataset)
    dataset.write_mtz(f'processed_data/aggregate/EF_10MV_cm_{phase}.mtz')

for epoch in [30, 60, 90, 120, 150]:
    for phase in ['pos', 'neg', 'zero']:
        for chain in range(n_chains):
            complex_reflections = np.zeros(n_reflections, dtype='complex128')
            for replicate in range(n_replicates):
                dataset = rs.read_mtz(f'processed_data/{replicate}/{epoch}/EF_10MV_cm_partial_{epoch}_chainwise_{phase}_subtraj_{chain}_avg.mtz')
                complex_reflections += mtz_to_cartesian_arr(dataset)
            complex_reflections /= n_replicates
            dataset[:] = cartesian_arr_to_polar(complex_reflections)
            dataset.infer_mtz_dtypes(inplace = True)
            reformat_dataset(dataset)
            dataset.write_mtz(f'processed_data/aggregate/{epoch}/EF_10MV_cm_chainwise_{phase}_subtraj_{chain}_avg.mtz')

    compute_all_difference_maps(f'processed_data/aggregate/{epoch}/EF_10MV_cm')