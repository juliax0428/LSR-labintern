import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 
from astropy.io import ascii
from astropy.table import Table

def readCDSTable(tableName, ReadMeName=None):
    if ReadMeName:
        table = ascii.read(tableName, readme=ReadMeName)
    else:
        table = ascii.read(tableName)  # Reading without a ReadMe file
    return table

def hms2deg(h, m, s):
    return (15. * s / 3600.) + (15. * m / 60.) + (h * 15.)

def dms2deg(sign, d, m, s):
    dec_deg = d + m / 60.0 + s / 3600.0
    # Apply the sign vectorized
    dec_deg = np.where(sign == '-', -dec_deg, dec_deg)
    return dec_deg

def distance(plx):
    return 1000. / plx  # mas to pc

def apparant_Gmag(Gmag, plx):
    return 5 * np.log10(distance(plx)) - 5 + Gmag

ClusterTableName = '/Users/xxz/Desktop/LSR-labintern/J_A+A_685_A83/clusters.dat'
MemberTableName_1 = '/Users/xxz/Desktop/LSR-labintern/J_A+A_685_A83/tablec1.dat'
MemberTableName_2 = '/Users/xxz/Desktop/LSR-labintern/J_A+A_685_A83/tablec2.dat'
MemberTableName_3 = '/Users/xxz/Desktop/LSR-labintern/J_A+A_685_A83/tablec3.dat'
ReadMeName = '/Users/xxz/Desktop/LSR-labintern/J_A+A_685_A83/ReadMe'
clusters = readCDSTable(tableName=ClusterTableName, ReadMeName=ReadMeName)
ges_param = readCDSTable(tableName=MemberTableName_1, ReadMeName=ReadMeName)
gaia_param = readCDSTable(tableName=MemberTableName_2, ReadMeName=ReadMeName)

hsc2686 = pd.read_csv('/Users/xxz/Desktop/LSR-labintern/HSC_2686.csv')
lynga3 = pd.read_csv('/Users/xxz/Desktop/LSR-labintern/Lynga3.csv')
gaia_stars = pd.read_csv('/Users/xxz/Desktop/LSR-labintern/all_GAIA_stars_in_area.csv')

i = 0
cluster_df = []
clusters_df = []
for cluster in clusters:
    if i >= 5:
        break
    print('cluster number is', cluster['Cluster'], ', cluster name is', cluster['Name'])
    ges_df = ges_param[ges_param['No'] == cluster['Cluster']]
    gaia_df = gaia_param[gaia_param['No'] == cluster['Cluster']]

    RA_deg = hms2deg(ges_df['RAh'], ges_df['RAm'], ges_df['RAs'])
    DE_deg = dms2deg(ges_df['DE-'], ges_df['DEd'], ges_df['DEm'], ges_df['DEs'])
    ges_df['RA_deg'] = RA_deg
    ges_df['DE_deg'] = DE_deg
    ges_df['NAME'] = cluster['Name']

    ges_df = ges_df.to_pandas()
    gaia_df = gaia_df.to_pandas()
    cluster_df = pd.merge(ges_df, gaia_df, on='CNAME')
    clusters_df.append(cluster_df)
    i += 1
clusters = pd.concat(clusters_df, ignore_index=True)

def GaiaDR3_CMD():
    '''Color Magnitude Diagram'''
    fig, ax = plt.subplots()
    ax.scatter(gaia_stars['bp_rp'], gaia_stars['phot_g_mean_mag'], 
               color='grey', label='All stars in DR3', s=0.1)
    ax.scatter(hsc2686['bp_rp'], hsc2686['phot_g_mean_mag'], color='blue', label='HSC_2686', s=4)
    ax.scatter(lynga3['bp_rp'], lynga3['phot_g_mean_mag'], color='red', label='Lynga_3', s=4)
    for name in clusters['NAME'].unique():
        cluster_data = clusters[clusters['NAME'] == name]
        ax.scatter(cluster_data['BP-RP'], cluster_data['Gmag'], label=name, s=1)
    
    ax.set_xlabel('G_BP - G_RP')
    ax.set_ylabel('G Magnitude')
    ax.set_title('Gaia DR3 CMD')
    ax.invert_yaxis()
    ax.legend()
    plt.show()

def GaiaDR3_PM():
    '''Proper Motion Diagram'''
    fig, ax = plt.subplots()
    ax.scatter(gaia_stars['pmra'], gaia_stars['pmdec'], color='grey', 
               label='All stars in DR3', s=0.1)
    ax.scatter(hsc2686['pmra'], hsc2686['pmdec'], color='blue', label='HSC_2686', s=3)
    ax.scatter(lynga3['pmra'], lynga3['pmdec'], color='red', label='Lynga_3', s=3)

    ax.plot(np.mean(hsc2686['pmra']), np.mean(hsc2686['pmdec']), '*', 
            color='orange', label='Mean Proper Motion of HSC_2686')
    ax.plot(np.mean(lynga3['pmra']), np.mean(lynga3['pmdec']), '*',
            color='yellow', label='Mean Proper Motion of Lynga 3')

    for name in clusters['NAME'].unique():
        cluster_data = clusters[clusters['NAME'] == name]
        ax.scatter(cluster_data['pmRA'], cluster_data['pmDE'], label=name, s=1)

    ax.set_xlabel('Proper motion in RA')
    ax.set_ylabel('Proper Motion in DE')
    ax.set_title('Gaia DR3 Proper Motion')
    ax.legend()
    plt.show()

def GaiaDR3_Plx():
    '''Parallax vs Radial Velocity'''
    fig, ax = plt.subplots()
    
    ax.scatter(gaia_stars['radial_velocity'], gaia_stars['parallax'], color='grey', 
               label='All stars in DR3', s=0.1)
    ax.scatter(hsc2686['radial_velocity'], hsc2686['parallax'], color='blue', label='HSC_2686', s=3)
    ax.scatter(lynga3['radial_velocity'], lynga3['parallax'], color='red', label='Lynga_3', s=3)

    for name in clusters['NAME'].unique():
        cluster_data = clusters[clusters['NAME'] == name]
        ax.scatter(cluster_data['RV'], cluster_data['plx'], label=name, s=1)

    ax.set_xlabel('Radial Velocity')
    ax.set_ylabel('Parallax')
    ax.set_title('Gaia DR3 Parallax against Radial Velocity')
    ax.legend()
    plt.show()

def plot_isochrones_log(ax, iso_file, cluster_file):
    '''
    Plot a range of isochrones for one Open Cluster
    '''
    if cluster_file.empty:
        print(f"Cluster file is empty for cluster: {iso_file['NAME'].values[0] if not iso_file.empty else 'Unknown'}")
        return

    unique_logAges = np.unique(iso_file['logAge'])
    for age in unique_logAges:
#        print('cluster_file[plx] = ',cluster_file['plx'])
#        STOP
        iso_df = iso_file[iso_file['logAge']==age]
        g_mag = apparant_Gmag(iso_df['mbolmag'], np.mean(cluster_file['plx']))
        bp_rp = iso_df['G_BPmag'] - iso_df['G_RPmag']
        ax.scatter(bp_rp, g_mag, label=f'{cluster_file["NAME"].values[0]} logAge={age}', s=1, marker=',')


def CMD_iso():
    iso_names = ['NGC6530', 'rho_oph', 'Trumpler_14', 'Cha_1', 'NGC2244']
    iso_files = {iso: f'/Users/xxz/Desktop/LSR-labintern/J_A+A_685_A83/{iso}.dat' for iso in iso_names}
    isochrone_data = {key: readCDSTable(file) for key, file in iso_files.items()}

    ash_1 = '/Users/xxz/Desktop/LSR-labintern/isochrone/ash1(logAge=8.8).dat'
    Majaess_190 = '/Users/xxz/Desktop/LSR-labintern/isochrone/majaess190(logAge=9.6).dat'
    T_ash = readCDSTable(tableName=ash_1)
    T_majaess = readCDSTable(tableName=Majaess_190)

    name_mapping = {
        'NGC 6530': 'NGC6530',
        'Rho Oph': 'rho_oph',
        'Trumpler 14': 'Trumpler_14',
        'Cha I': 'Cha_1',
        'NGC 2244': 'NGC2244'
    }

    # Loop through each open cluster and corresponding isochrone
    for name, iso_name in name_mapping.items():
        iso_file = isochrone_data[iso_name]
        cluster_file = clusters[clusters['NAME'] == name]

        fig, ax = plt.subplots()
        ax.scatter(hsc2686['bp_rp'], hsc2686['phot_g_mean_mag'], color='blue', label='HSC_2686', s=5)
        ax.scatter(lynga3['bp_rp'], lynga3['phot_g_mean_mag'], color='red', label='Lynga_3', s=5)

        for name in clusters['NAME'].unique():
            cluster_data = clusters[clusters['NAME'] == name]
            ax.scatter(cluster_data['BP-RP'], cluster_data['Gmag'], label=name, s=2)

        plot_isochrones_log(ax, iso_file, cluster_file)

        ax.set_xlabel('G_BP-G_RP')
        ax.set_ylabel('G Magnitude')
        ax.set_title(f'Gaia DR3 CMD with Isochrones for {iso_name}')
        ax.invert_yaxis()
        ax.legend()
        plt.show()

# Uncomment the desired function to run
# GaiaDR3_PM()
# GaiaDR3_Plx()
# GaiaDR3_CMD()
CMD_iso()