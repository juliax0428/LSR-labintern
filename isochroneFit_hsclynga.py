'''
Python code for fitting isochrones for hsc 2686 and lynga 3. 
'''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.io import ascii

def readCDSTable(tableName, ReadMeName=None):
    if ReadMeName:
        table = ascii.read(tableName, readme=ReadMeName)
    else:
        table = ascii.read(tableName)  # Reading without a ReadMe file
    return table

def apparant_Gmag(Gmag, dist):
    return 5 * np.log10(dist) - 5 + Gmag

def G_mag_uncert(FG, e_FG):
    sigmaG_0 = 0.0027553202
    e_Gmag = np.sqrt((-2.5/np.log(10)*e_FG/FG)**2 + sigmaG_0**2)
    return e_Gmag
def G_BP_uncert(FGBP, e_FGBP):
    sigmaGBP_0 = 0.0027901700
    e_GBPmag = np.sqrt((-2.5/np.log(10)*e_FGBP/FGBP)**2 + sigmaGBP_0**2)
    return e_GBPmag
def G_RP_uncert(FGRP, e_FGRP):
    sigmaGRP_0 = 0.0037793818
    e_GRPmag= np.sqrt((-2.5/np.log(10)*e_FGRP/FGRP)**2 +sigmaGRP_0**2)
    return e_GRPmag


'''tableName = '/Users/xxz/Desktop/LSR-labintern/J_A+A_673_A114/clusters.dat'
memberName = '/Users/xxz/Desktop/LSR-labintern/J_A+A_673_A114/members.dat'
ReadMeName = '/Users/xxz/Desktop/LSR-labintern/J_A+A_673_A114/ReadMe'
clusters = readCDSTable(tableName=tableName,ReadMeName=ReadMeName)
members = readCDSTable(tableName=memberName,ReadMeName=ReadMeName)

cluster_data_hsc = []
cluster_data_lynga = []
for cluster in clusters:
    if cluster['Name'] =='HSC_2686':
        print('found cluster. name = ',cluster['Name'])
        cluster_data_hsc.append(cluster)
    elif cluster['Name'] =='Lynga_3':
        print('found cluster. name = ',cluster['Name'])
        cluster_data_lynga.append(cluster)

print('HSC2686 type = ', cluster_data_hsc[0]['Type'])
print('HSC2686 has logAge16 = ', cluster_data_hsc[0]['logAge16'])
print('HSC2686 has AV16 = ', cluster_data_hsc[0]['AV16'])
print('HSC2686 has dist16 = ', cluster_data_hsc[0]['dist16'])

print('Lynga3 type = ', cluster_data_lynga[0]['Type'])
print('Lynga3 has logAge16 = ', cluster_data_lynga[0]['logAge16'])
print('Lynga3 has AV16 = ', cluster_data_lynga[0]['AV16'])
print('Lynga3 has dist16 = ', cluster_data_lynga[0]['dist16'])'''


'''
HSC2686 has logAge16 =  6.4626441
HSC2686 has AV16 =  5.84
HSC2686 has dist16 =  4995.7033353

Lynga3 has logAge16 =  6.40001631
Lynga3 has AV16 =  5.925
Lynga3 has dist16 =  4928.85577221
'''

hsc2686 = pd.read_csv('/Users/xxz/Desktop/LSR-labintern/HSC_2686.csv')
lynga3 = pd.read_csv('/Users/xxz/Desktop/LSR-labintern/Lynga3.csv')

#isochrone: logAge 6.4-6.5, step=0.01
iso_hsc = '/Users/xxz/Desktop/LSR-labintern/iso_hsc2686.dat'
iso_lynga = '/Users/xxz/Desktop/LSR-labintern/iso_lynga3.dat'
T_iso_hsc = readCDSTable(tableName=iso_hsc)
T_iso_lynga = readCDSTable(tableName=iso_lynga)

def isochrone_fit_hsc(hsc):
    fig, ax = plt.subplots()
    unique_logAges = np.unique(hsc['logAge'])
    for age in unique_logAges:
        iso_df = hsc[hsc['logAge']==age]
        g_mag = apparant_Gmag(iso_df['mbolmag'], 4996.7033351)
        bp_rp = iso_df['G_BPmag'] - iso_df['G_RPmag']
        ax.scatter(bp_rp, g_mag, label=f'HSC 2686 logAge={age}', s=0.5, marker=',')
    ax.scatter(hsc2686['bp_rp'], hsc2686['phot_g_mean_mag'], label='HSC 2686', s=8)
    ax.errorbar(hsc2686['bp_rp'], hsc2686['phot_g_mean_mag'],
                xerr=abs(G_BP_uncert(hsc2686['phot_bp_mean_flux'],hsc2686['phot_bp_mean_flux_error']) + G_RP_uncert(hsc2686['phot_rp_mean_flux'],hsc2686['phot_rp_mean_flux_error'])), 
                yerr=abs(G_mag_uncert(hsc2686['phot_g_mean_flux'], hsc2686['phot_g_mean_flux_error'])), 
                fmt='o', capsize=2, elinewidth=0.5, markersize=3)    
    ax.invert_yaxis()
    ax.set_xlabel('G_BP - G_RP',fontsize=14)
    ax.set_ylabel('G Magnitude',fontsize=14)
    ax.set_title(f'isochrone fit for HSC 2686',fontsize=16)
    ax.legend()
    plt.savefig(fname='iso_hsc2686.png', dpi=300)
    plt.show()

def isochrone_fit_lynga(lynga):
    fig, ax = plt.subplots()
    unique_logAges = np.unique(lynga['logAge'])
    for age in unique_logAges:
        iso_df = lynga[lynga['logAge']==age]
        g_mag = apparant_Gmag(iso_df['mbolmag'], 4928.85577221)
        bp_rp = iso_df['G_BPmag'] - iso_df['G_RPmag']
        ax.scatter(bp_rp, g_mag, label=f'Lynga 3 logAge={age}', s=0.5, marker=',')
    ax.scatter(lynga3['bp_rp'], lynga3['phot_g_mean_mag'], label='Lynga 3', s=8)
    ax.errorbar(lynga3['bp_rp'], lynga3['phot_g_mean_mag'],
                xerr=abs(G_BP_uncert(lynga3['phot_bp_mean_flux'],lynga3['phot_bp_mean_flux_error']) + G_RP_uncert(lynga3['phot_rp_mean_flux'],lynga3['phot_rp_mean_flux_error'])), 
                yerr=abs(G_mag_uncert(lynga3['phot_g_mean_flux'], lynga3['phot_g_mean_flux_error'])),  
                fmt='o', capsize=2, elinewidth=0.5, markersize=3)
    ax.invert_yaxis()
    ax.set_xlabel('G_BP - G_RP',fontsize=14)
    ax.set_ylabel('G Magnitude',fontsize=14)
    ax.set_title(f'isochrone fit for Lynga 3',fontsize=16)
    ax.legend()
    plt.savefig(fname='iso_lynga3.png', dpi=300)
    plt.show()

isochrone_fit_hsc(T_iso_hsc)
isochrone_fit_lynga(T_iso_lynga)