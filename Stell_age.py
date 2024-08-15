from astroquery.gaia import Gaia
import pandas as pd 

hsc2686 = pd.read_csv('/Users/xxz/LSR-labintern/HSC_2686.csv')
lynga3 = pd.read_csv('/Users/xxz/LSR-labintern/Lynga3.csv')

retrieval_type = 'ALL'          # Options are: 'EPOCH_PHOTOMETRY', 'MCMC_GSPPHOT', 'MCMC_MSC', 'XP_SAMPLED', 'XP_CONTINUOUS', 'RVS', 'ALL'
data_structure = 'INDIVIDUAL'     # Options are: 'INDIVIDUAL' or 'RAW'
data_release   = 'Gaia DR3'     # Options are: 'Gaia DR3' (default), 'Gaia DR2'

print(Gaia.load_data.__code__.co_argcount,': ',Gaia.load_data.__code__.co_varnames)

datalink = Gaia.load_data(ids=[hsc2686['SOURCE_ID']] and [lynga3['SOURCE_ID']],# <============ replace the array [sourceID] with an array with your sourceIDs
                          #data_release=data_release,
                          retrieval_type=retrieval_type,
                          #data_structure=data_structure,
                          verbose=True
                          )

print('datalink = ',datalink)

dl_keys  = [inp for inp in datalink.keys()]
dl_keys.sort()
print(f'The following Datalink products have been downloaded:')

for dl_key in dl_keys:
    print(f' * {dl_key}')
