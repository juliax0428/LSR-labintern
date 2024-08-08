import numpy as np
import pandas as pd
from astropy.io import ascii

# Read the data
lynga3 = pd.read_csv('/Users/xxz/LSR-labintern/Lynga3.csv')
lynga3_astro_param = pd.read_csv('/Users/xxz/LSR-labintern/astro_param_Lynga3.csv')
hsc2686 = pd.read_csv('/Users/xxz/LSR-labintern/HSC_2686.csv')
hsc2686_astro_param = pd.read_csv('/Users/xxz/LSR-labintern/astro_param_HSC_2686.csv')

# Extract coordinates
hsc2686_coord_ra = hsc2686['ra']
hsc2686_coord_dec = hsc2686['dec']
lynga3_coord_ra = lynga3['ra']
lynga3_coord_dec = lynga3['dec']

# Function to read CDS Table
def readCDSTable(tableName, ReadMeName):
    table = ascii.read(tableName, readme=ReadMeName)#, guess=False, fast_reader=False)
    return table

# Read LAMOST data
tableName = '/Users/xxz/Desktop/LAMOST_DR9/LAMA_stars.dat'
ReadMeName = '/Users/xxz/Desktop/LAMOST_DR9/LAMA_stars_ReadMe.txt'
lama = readCDSTable(tableName=tableName, ReadMeName=ReadMeName)

# Extract LAMOST coordinates
lama_ra = lama['RAdeg']
lama_dec = lama['DEdeg']
print('lama_ra = ',lama_ra)

# Function to calculate angular distance
# all angles must be in degrees
#@return: angular distance in arcsec
def angular_distance(ra1, dec1, ra2, dec2):
    from PyAstronomy import pyasl
    return pyasl.getAngDist(ra1,dec1,ra2,dec2)*3600.

# Calculate distances for HSC_2686
smallest_distances_hsc2686 = []
for j in range(len(hsc2686_coord_ra)):
    distances_hsc2686 = []
    for i in range(len(lama_ra)):
        distance = angular_distance(lama_ra[i], lama_dec[i], hsc2686_coord_ra[j], hsc2686_coord_dec[j])
        distances_hsc2686.append({'lama_id': lama['uid'][i], 'hsc2686_id': hsc2686['SOURCE_ID'][j], 'distance': distance})
    sorted_distances_hsc2686 = sorted(distances_hsc2686,key=lambda d: d['distance'])
    print('j = ',j,': smallest distance: ',sorted_distances_hsc2686[0])
    smallest_distances_hsc2686.append(sorted_distances_hsc2686[0])

for i in range(len(smallest_distances_hsc2686)):
    print('hsc2686_id = ',smallest_distances_hsc2686[i]['hsc2686_id'],': distance = ',smallest_distances_hsc2686[i]['distance'])

# Calculate distances for Lynga3
smallest_distances_lynga3 = []
for j in range(len(lynga3_coord_ra)):
    distances_lynga3 = []
    for i in range(len(lama_ra)):
        distance = angular_distance(lama_ra[i], lama_dec[i], lynga3_coord_ra[j], lynga3_coord_dec[j])
        distances_lynga3.append({'lama_id': lama['uid'][i], 'lynga3_id': lynga3['SOURCE_ID'][j], 'distance': distance})
    sorted_distances_lynga3 = sorted(distances_lynga3,key=lambda d: d['distance'])
    print('j = ',j,': smallest distance: ',sorted_distances_lynga3[0])
    smallest_distances_lynga3.append(sorted_distances_lynga3[0])

for i in range(len(smallest_distances_lynga3)):
    print('lynga3_id = ',smallest_distances_lynga3[i]['lynga3_id'],': distance = ',smallest_distances_lynga3[i]['distance'])

