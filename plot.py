'''
An attempt to plot similar plots as Fig.4 in 
Testing Cluster Membership of Planetary Nebulae with High-Precision Proper Motions. I. HST Observations of JaFu 1 Near the Globular Cluster Palomar 6
'''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 
from matplotlib.patches import Circle
from gaia_query import readCDSTable
from gaia_query import hmsToDeg
from gaia_query import dmsToDeg
from PyAstronomy import pyasl


# Read the data
hsc2686 = pd.read_csv('/Users/xxz/Desktop/LSR-labintern/HSC_2686.csv')
lynga3 = pd.read_csv('/Users/xxz/Desktop/LSR-labintern/Lynga3.csv')
#pleiades = pd.read_csv('/Users/xxz/Desktop/LSR-labintern/Pleiades.csv', sep=';', skiprows=2)
gaia_stars = pd.read_csv('/Users/xxz/Desktop/LSR-labintern/all_GAIA_stars_in_area.csv')
#czernik = pd.read_csv('/Users/xxz/Desktop/LSR-labintern/czernik20_cleaned.csv', sep=',')

#print('czernik', czernik)
#print('czernik.keys', czernik.columns)
#print('czernik datatype', czernik.dtypes)


#print('pleiades data is : \n',pleiades)
#print('pleiades.keys() \n', pleiades.keys())
#print('type of pleiades', pleiades.dtypes)
#print("czernik['Gmag']", czernik['Gmag'], 'datatype', czernik['Gmag'].dtype)
#print(czernik['BP-RP'].astype(float))

#Position of PN ([PN RA/DEC]: 15:16:41.00 -58:22:26.00)
pos_PN_hms = ('15:16:41.00')
pos_PN_dms = ('-58:22:26.00')
pos_PN_ra_deg = hmsToDeg(pos_PN_hms)
pos_PN_dec_deg = dmsToDeg(pos_PN_dms)
print ('Position of PN PN_ra_deg =', pos_PN_ra_deg, 'Position of PN PN_dec_deg', pos_PN_dec_deg)

#Find the Central Star
minDist = 1000.
central_star = None
for idx,pn in hsc2686.iterrows():
    #print('pn = ',pn)
    #print('pn[ra] = ',type(pn['ra']),': ',pn['ra'])
    dist = pyasl.getAngDist(pos_PN_ra_deg, pos_PN_dec_deg, pn['ra'], pn['dec'])
    if dist < minDist:
        minDist = dist
        central_star = pn
    if dist < 1/3600.:
        #print ('Source ID of central star is', pn['SOURCE_ID'])
        central_star_id = pn['SOURCE_ID']
print('minDist from PN 15:16:41.00 -58:22:26.00 = ',minDist, '\n Central star = ',central_star) 
#Central Star SOURCE_ID: 5877088151005347072

def apparant_Gmag(Gmag, d):
    return 5*np.log(d) - 5 + Gmag
def distance(plx):
    return plx*(10**3)

# Panel (a): 
# Positions diagram of our cluster stars.
def GaiaDR3_map():
    tableName = '/Users/xxz/Desktop/LSR-labintern/J_A+A_673_A114/clusters.dat'
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
    circle_hsc2686 = Circle((cluster_data_hsc[0]['RAdeg'], cluster_data_hsc[0]['DEdeg']), cluster_data_hsc[0]['rtot'],
                                color='blue', fill=False, linestyle='dashed')
    circle_lynga3 = Circle((cluster_data_lynga[0]['RAdeg'], cluster_data_lynga[0]['DEdeg']), cluster_data_lynga[0]['rtot'],
                               color='red', fill=False, linestyle='dashed')

    fig, ax = plt.subplots()
    ax.scatter(gaia_stars['ra'],  gaia_stars['dec'], color='grey', label='all stars', s=0.1)
    ax.scatter(hsc2686['ra'], hsc2686['dec'], color='blue', label='HSC_2686', s=4)
    ax.scatter(lynga3['ra'], lynga3['dec'], color='red', label='Lynga_3', s=4)
    ax.add_patch(circle_hsc2686)
    ax.add_patch(circle_lynga3)
    ax.plot(cluster_data_hsc[0]['RAdeg'],cluster_data_hsc[0]['DEdeg'], 
            'x', markersize=5, color = 'blue', label = 'Position of cluster HSC 2686')
    ax.plot(cluster_data_lynga[0]['RAdeg'],cluster_data_lynga[0]['DEdeg'], 
            'x', markersize=5, color = 'red', label = 'Position of cluster Lynga 3')
    ax.plot(np.mean(hsc2686['ra']), np.mean(hsc2686['dec']), '*', 
            color='blue', label='Mean position of stars of HSC 2686')
    ax.plot(np.mean(lynga3['ra']), np.mean(lynga3['dec']), '*',
            color = 'red', label = 'Mean position of stars of Lynga 3')
    ax.plot(pos_PN_ra_deg, pos_PN_dec_deg, 'D', markersize=4, color = 'green', label = 'Central Star of PN')

    ax.set_xlabel('Right Ascension (RA)',fontsize=18)
    ax.set_ylabel('Declination (DEC)',fontsize=18)
    ax.set_title('Gaia DR3 Map',fontsize=20)
    ax.legend()
    plt.savefig(fname='Gaia_map.png', dpi=300)
    plt.show()

# Panel(b):
# Color-magnitude diagram of all the stars. G magnitude versus GBP − GRP color.
def GaiaDR3_CMD():
    '''Color Magnitude diagram'''
    fig, ax = plt.subplots()
    ax.scatter(gaia_stars['bp_rp'], apparant_Gmag(gaia_stars['phot_g_mean_mag'],distance(gaia_stars['parallax'])), 
               color='grey', label='All stars in DR3', s=0.1)
    ax.scatter(hsc2686['bp_rp'], apparant_Gmag(hsc2686['phot_g_mean_mag'],distance(hsc2686['parallax'])), color='blue', label='HSC_2686', s=4)
    ax.scatter(lynga3['bp_rp'], apparant_Gmag(lynga3['phot_g_mean_mag'], distance(lynga3['parallax'])), color='red', label='Lynga_3', s=4)
    #ax.scatter(czernik['BP-RP'], apparant_Gmag(czernik['Gmag'], distance(czernik['Plx'])), label='czernik20', s=3)
    #ax.scatter(pleiades['BP-RP'].astype(float), pleiades['Gmag'].astype(float), color='purple', label='pleiades', s=4)
    ax.plot(central_star['bp_rp'], apparant_Gmag(central_star['phot_g_mean_mag'],distance(central_star['parallax'])), marker='D', markersize='3', color='green', label='Central Star')

    ax.set_xlabel('G_BP - G_RP')
    ax.set_ylabel('G Magnitude')
    ax.set_title('Gaia DR3 CMD')
    ax.invert_yaxis()
    ax.legend()
    plt.show()

# Panel (c):
# Vector-point diagram
def GaiaDR3_PM():
    '''Proper Motion Diagram'''
    fig, ax = plt.subplots()
    ax.scatter(gaia_stars['pmra'], gaia_stars['pmdec'], color='grey', 
               label='All stars in DR3', s=0.1)
    ax.scatter(hsc2686['pmra'], hsc2686['pmdec'], color='blue', label='HSC_2686', s=4)
    ax.scatter(lynga3['pmra'], lynga3['pmdec'], color='red', label='Lynga_3', s=4)
    #ax.scatter(pleiades['pmRA'].astype(float), pleiades['pmDE'].astype(float), color='purple', label='pleiades', s=3)
    #ax.scatter(czernik['pmRA'], czernik['pmDE'], label='czernik20', s=3)
    ax.plot(np.mean(hsc2686['pmra']), np.mean(hsc2686['pmdec']), '+',  label='Mean Proper Motion of hsc2686')
    ax.plot(np.mean(lynga3['pmra']), np.mean(lynga3['pmdec']), '+', label = 'Mean Proper Motion of Lynga 3')
    ax.plot(central_star['pmra'], central_star['pmdec'], marker='D', markersize='3', color='green', label='Central Star')
    #ax.set_xlim(-6,-4.4)
    #ax.set_ylim(-4.4,-2.8)
    ax.set_xlabel('Proper motion in RA', fontsize=14)
    ax.set_ylabel('Proper Motion in DE', fontsize=14)
    ax.set_title('Gaia DR3 Proper Motion', fontsize=16)
    ax.legend()
    plt.savefig('Gaia_pm.png', dpi=1000)
    plt.show()

#Panel (f):
#Vector-point diagram for the absolute PMs
def GaiaDR3_PM_uncertainties():
    fig, ax = plt.subplots()

    # Apply filters
    filter_hsc2686 = (hsc2686['pmra_error'] < 2) & (hsc2686['pmdec_error'] < 2)
    filter_lynga3 = (lynga3['pmra_error'] < 2) & (lynga3['pmdec_error'] < 2)
    filter_all_gaia = (gaia_stars['pmra_error'] < 2) & (gaia_stars['pmdec_error'] < 2)

    # Scatter plots
    ax.scatter(hsc2686[filter_hsc2686]['pmra'], hsc2686[filter_hsc2686]['pmdec'], 
               color='blue', label='HSC 2686', s=5)
    ax.scatter(lynga3[filter_lynga3]['pmra'], lynga3[filter_lynga3]['pmdec'], 
               color='red', label='Lynga 3', s=5)
    ax.scatter(gaia_stars[filter_all_gaia]['pmra'], gaia_stars[filter_all_gaia]['pmdec'], 
               color='grey', label='All stars in DR3', s=0.5)
    
    # Error bars
    ax.errorbar(hsc2686[filter_hsc2686]['pmra'], hsc2686[filter_hsc2686]['pmdec'],
                xerr=abs(hsc2686[filter_hsc2686]['pmra_error']), 
                yerr=abs(hsc2686[filter_hsc2686]['pmdec_error']), 
                fmt='o', color='blue', capsize=2, elinewidth=0.3, markersize=2)
    ax.errorbar(lynga3[filter_lynga3]['pmra'], lynga3[filter_lynga3]['pmdec'],
                xerr=abs(lynga3[filter_lynga3]['pmra_error']), 
                yerr=abs(lynga3[filter_lynga3]['pmdec_error']), 
                fmt='o', color='red', capsize=2, elinewidth=0.3, markersize=2)
    
    ax.plot(central_star['pmra'], central_star['pmdec'], marker='D', markersize=8,
            color='green', label='Central Star')
    ax.plot(np.mean(hsc2686['pmra']), np.mean(hsc2686['pmdec']), 
            marker='+', color='red', markersize=20, label='Mean PM of HSC 2686')
    ax.plot(np.mean(lynga3['pmra']), np.mean(lynga3['pmdec']), 
            marker='+', color='blue', markersize=20, label='Mean PM of Lynga 3')
    # Dynamic limits with padding
    ax.set_xlim(hsc2686['pmra'].min()-0.2, hsc2686['pmra'].max()+0.4)
    ax.set_ylim(hsc2686['pmdec'].min()-0.8, hsc2686['pmdec'].max()+0.5)
    
    # Labels and Title
    ax.set_xlabel('Proper motion in RA (mas/yr)', fontsize=14)
    ax.set_ylabel('Proper Motion in DEC (mas/yr)', fontsize=14)
    ax.set_title('Gaia DR3 Proper Motions with Uncertainty', fontsize=16)
    ax.legend()
    plt.savefig(fname=f'Gaia_pm_uncert.png', dpi=300)
    plt.show()

#Parallax vs Radial Velocity
def GaiaDR3_Plx():
    '''Parallax vs Radial Velocity'''
    fig, ax = plt.subplots()

    ax.scatter(gaia_stars['radial_velocity'], gaia_stars['parallax'], color='grey', 
               label='All stars in DR3', s=0.1)
    ax.scatter(hsc2686['radial_velocity'], hsc2686['parallax'], color='blue', label='HSC_2686', s=10)
    ax.scatter(lynga3['radial_velocity'], lynga3['parallax'], color='red', label='Lynga_3', s=10)
    ax.plot(central_star['radial_velocity'], central_star['parallax'], marker='D', markersize='8', color='green', label='Central Star')
    ax.set_xlim(-70,30)
    ax.set_ylim(0,0.6)
    ax.set_xlabel('Radial Velocity', fontsize=14)
    ax.set_ylabel('Parallax', fontsize=14)
    ax.set_title('Gaia DR3 Parallax against Radial Velocity', fontsize=16)
    ax.legend()
    plt.savefig(fname='plx_rv.png', dpi=300)
    plt.show()



GaiaDR3_map()
#GaiaDR3_CMD()
#GaiaDR3_PM()
#GaiaDR3_PM_uncertainties()
#GaiaDR3_Plx()