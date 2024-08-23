import pandas as pd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import ascii

hsc2686 = pd.read_csv('/Users/xxz/Desktop/LSR-labintern/HSC_2686.csv')
lynga3 = pd.read_csv('/Users/xxz/Desktop/LSR-labintern/Lynga3.csv')
all_gaia_stars = pd.read_csv('/Users/xxz/Desktop/LSR-labintern/all_GAIA_stars_in_area.csv')

def readCDSTable(tableName):
    table = ascii.read(tableName)#, guess=False, fast_reader=False)
    return table

iso_000 = '/Users/xxz/Desktop/iso_omega_0.0.dat'
iso_030 = '/Users/xxz/Desktop/iso_omega_0.3.dat'
iso_060 = '/Users/xxz/Desktop/iso_omega_0.6.dat'
iso_080 = '/Users/xxz/Desktop/iso_omega_0.8.dat'
iso_090 = '/Users/xxz/Desktop/iso_omega_0.9.dat'
iso_095 = '/Users/xxz/Desktop/iso_omega_0.95.dat'
iso_099 = '/Users/xxz/Desktop/iso_omega_0.99.dat'
iso_logage = '/Users/xxz/Desktop/log_age_gaia_bands.dat'
iso_logage_2 = '/Users/xxz/Desktop/log_age_8_10.dat'

T_iso_000 = readCDSTable(tableName=iso_000)
#T_iso_030 = readCDSTable(tableName=iso_030)
#T_iso_060 = readCDSTable(tableName=iso_060)
#T_iso_080 = readCDSTable(tableName=iso_080)
#T_iso_090 = readCDSTable(tableName=iso_090)
#T_iso_095 = readCDSTable(tableName=iso_095)
#T_iso_099 = readCDSTable(tableName=iso_099)
T_iso_logage = readCDSTable(tableName=iso_logage)
T_iso_logage_2 = readCDSTable(tableName=iso_logage_2)

print (T_iso_000.keys())
print (T_iso_logage.keys())

'''
Johnson-Cousins relationships
'''
def G_mag(V, I):
    abs_G = V -0.01746 + 0.008092*(V-I) - 0.2810*(V-I)**2 + 0.3655*(V-I)**3
    return abs_G

def G_BP_RP(V,I):
    bp_rp = -0.04212 + 1.286*(V-I) - 0.09494*(V-I)**2
    return bp_rp


'''Try with T_iso_000'''
def plot_isochrones(data_file):

    unique_logAges = np.unique(data_file['logAge'])
    num_isochrones_per_plot = 10
    num_subplots = int(np.ceil(len(unique_logAges) / num_isochrones_per_plot))

    for i in range(0, num_subplots):
        fig, ax =plt.subplots()
        start_idx = i * num_isochrones_per_plot
        end_idx = min((i + 1) * num_isochrones_per_plot, len(unique_logAges))
        for age in unique_logAges[start_idx:end_idx]:
            new_df = data_file[data_file['logAge'] == age]
            g_mag = G_mag(new_df['Vmag'], new_df['Imag'])
            bp_rp = G_BP_RP(new_df['Vmag'], new_df['Imag'])
            ax.scatter(bp_rp, g_mag, label=f'logAge={age}', s=3)
        ax.scatter(hsc2686['bp_rp'], hsc2686['phot_g_mean_mag'], color='blue', label='HSC_2686', s=4)
        ax.scatter(lynga3['bp_rp'], lynga3['phot_g_mean_mag'], color='red', label='Lynga_3', s=4)
        #ax.scatter(all_gaia_stars['bp_rp'], all_gaia_stars['phot_g_mean_mag'], 
               #color='grey', label='All stars in DR3', s=0.1)
        ax.invert_yaxis()
        ax.set_xlabel('G_BP - G_RP')
        ax.set_ylabel('G Magnitude')
        ax.set_title(f'Gaia DR3 CMD subplots{i+1}')
        ax.legend()
        plt.show()

def plot_isochrones_2 (data_file):
    unique_logAges = np.unique(data_file['logAge'])
    num_isochrones_per_plot = 10
    num_subplots = int(np.ceil(len(unique_logAges) / num_isochrones_per_plot))

    for i in range(0, num_subplots):
        fig, ax =plt.subplots()
        start_idx = i * num_isochrones_per_plot
        end_idx = min((i + 1) * num_isochrones_per_plot, len(unique_logAges))
        for age in unique_logAges[start_idx:end_idx]:
            new_df = data_file[data_file['logAge'] == age]
            g_mag = new_df['Gmag']
            bp_rp = new_df['G_BPmag']-new_df['G_RPmag']
            ax.scatter(bp_rp, g_mag, label=f'logAge={age}', s=3)
        ax.scatter(hsc2686['bp_rp'], hsc2686['phot_g_mean_mag'], color='blue', label='HSC_2686', s=4)
        ax.scatter(lynga3['bp_rp'], lynga3['phot_g_mean_mag'], color='red', label='Lynga_3', s=4)
        #ax.scatter(all_gaia_stars['bp_rp'], all_gaia_stars['phot_g_mean_mag'], 
               #color='grey', label='All stars in DR3', s=0.1)
        ax.invert_yaxis()
        ax.set_xlabel('G_BP - G_RP')
        ax.set_ylabel('G Magnitude')
        ax.set_title(f'Gaia DR3 CMD subplots{i+1}')
        ax.legend()
        plt.show()

plot_isochrones(T_iso_000)
plot_isochrones_2(T_iso_logage)
plot_isochrones_2(T_iso_logage_2)