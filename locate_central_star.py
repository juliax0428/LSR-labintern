import pandas as pd 
from gaia_query import hmsToDeg
from gaia_query import dmsToDeg
from PyAstronomy import pyasl

hsc2686 = pd.read_csv('/Users/xxz/Desktop/LSR-labintern/HSC_2686.csv')
lynga3 = pd.read_csv('/Users/xxz/Desktop/LSR-labintern/Lynga3.csv')

#Position of PN ([PN RA/DEC]: 15:16:41.00 -58:22:26.00)
pos_PN_hms = ('15:16:41.00')
pos_PN_dms = ('-58:22:26.00')
pos_PN_ra_deg = hmsToDeg(pos_PN_hms)
pos_PN_dec_deg = dmsToDeg(pos_PN_dms)
print ('Central Star of PN_ra_deg =', pos_PN_ra_deg, 'Central Star of PN_dec_deg', pos_PN_dec_deg)

#Find the Central Star
minDist = 1000.
central_star = None
for idx,pn in hsc2686.iterrows():
    print('pn = ',pn)
    print('pn[ra] = ',type(pn['ra']),': ',pn['ra'])
    dist = pyasl.getAngDist(pos_PN_ra_deg, pos_PN_dec_deg, pn['ra'], pn['dec'])
    if dist < minDist:
        minDist = dist
        central_star = pn
    if dist < 1/3600.:
        print ('Source ID of central star is', pn['SOURCE_ID'])
print('minDist = ',minDist)
print('Central star = ',central_star)       


#Central Star SOURCE_ID: 5877088151005347072
