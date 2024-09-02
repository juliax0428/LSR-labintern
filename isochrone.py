import pandas as pd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import ascii

hsc2686 = pd.read_csv('/Users/xxz/Desktop/LSR-labintern/HSC_2686.csv')
lynga3 = pd.read_csv('/Users/xxz/Desktop/LSR-labintern/Lynga3.csv')
all_gaia_stars = pd.read_csv('/Users/xxz/Desktop/LSR-labintern/all_GAIA_stars_in_area.csv')
czernik = pd.read_csv('/Users/xxz/Desktop/LSR-labintern/czernik20_cleaned.csv', sep=',')

def readCDSTable(tableName):
    table = ascii.read(tableName)#, guess=False, fast_reader=False)
    return table
def apparant_Gmag(Gmag, d):
    return 5*np.log(d) - 5 + Gmag
def dist(plx):
    return 1000./plx

iso_000 = '/Users/xxz/Desktop/LSR-labintern/isochrone/iso_000.dat'
iso_030 = '/Users/xxz/Desktop/LSR-labintern/isochrone/iso_omega_0.3.dat'
iso_060 = '/Users/xxz/Desktop/LSR-labintern/isochrone/iso_omega_0.6.dat'
iso_080 = '/Users/xxz/Desktop/LSR-labintern/isochrone/iso_omega_0.8.dat'
iso_090 = '/Users/xxz/Desktop/LSR-labintern/isochrone/iso_omega_0.9.dat'
iso_095 = '/Users/xxz/Desktop/LSR-labintern/isochrone/iso_omega_0.95.dat'
iso_099 = '/Users/xxz/Desktop/LSR-labintern/isochrone/iso_omega_0.99.dat'
iso_logage = '/Users/xxz/Desktop/LSR-labintern/isochrone/log_age_gaia_bands.dat'
iso_logage_2 = '/Users/xxz/Desktop/LSR-labintern/isochrone/log_age_8_10.dat'
ash_1 = '/Users/xxz/Desktop/LSR-labintern/isochrone/ash1(logAge=8.8).dat'
Majaess_190 = '/Users/xxz/Desktop/LSR-labintern/isochrone/majaess190(logAge=9.6).dat'

T_iso_000 = readCDSTable(tableName=iso_000)
#T_iso_030 = readCDSTable(tableName=iso_030)
#T_iso_060 = readCDSTable(tableName=iso_060)
#T_iso_080 = readCDSTable(tableName=iso_080)
#T_iso_090 = readCDSTable(tableName=iso_090)
#T_iso_095 = readCDSTable(tableName=iso_095)
#T_iso_099 = readCDSTable(tableName=iso_099)
T_iso_logage = readCDSTable(tableName=iso_logage)
T_iso_logage_2 = readCDSTable(tableName=iso_logage_2)
T_ash = readCDSTable(tableName=ash_1)
T_majaess = readCDSTable(tableName=Majaess_190)

#print (T_ash.keys())
#print (T_majaess.keys())
#print(T_ash.dtype)

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
def plot_isochrones_linear(data_file):

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
            ax.scatter(bp_rp, g_mag, label=f'logAge={age}={10**age/(10**9)}Gyrs', s=2, marker=',')
        ax.scatter(hsc2686['bp_rp'], hsc2686['phot_g_mean_mag'], color='blue', label='HSC_2686', s=4)
        ax.scatter(lynga3['bp_rp'], lynga3['phot_g_mean_mag'], color='red', label='Lynga_3', s=4)
        ax.scatter(czernik['BP-RP'], czernik['Gmag'], color='purple', label='czernik20', s=4)
        #ax.scatter(all_gaia_stars['bp_rp'], all_gaia_stars['phot_g_mean_mag'], 
               #color='grey', label='All stars in DR3', s=0.1)
        ax.invert_yaxis()
        ax.set_xlabel('G_BP - G_RP')
        ax.set_ylabel('G Magnitude')
        ax.set_title(f'Gaia DR3 CMD subplots{i+1}')
        ax.legend()
        plt.show()

def plot_isochrones_log (data_file):
    unique_logAges = np.unique(data_file['logAge'])
    num_isochrones_per_plot = 10
    num_subplots = int(np.ceil(len(unique_logAges) / num_isochrones_per_plot))

    for i in range(0, num_subplots):
        fig, ax =plt.subplots()
        start_idx = i * num_isochrones_per_plot
        end_idx = min((i + 1) * num_isochrones_per_plot, len(unique_logAges))
        for age in unique_logAges[start_idx:end_idx]:
            new_df = data_file[np.where(data_file['logAge'] == age)[0]]
            g_mag = apparant_Gmag(new_df['Gmag'], dist(new_df['Plx']))
            print('g_mag = ',g_mag)
            print('np.mean(g_mag) = ',np.mean(np.array(g_mag)))
            bp_rp = new_df['G_BPmag']-new_df['G_RPmag']
            ax.scatter(bp_rp, g_mag, label=f'logAge={age}={10**age/(10**9)}Gyrs', s=2, marker=',')
        ax.scatter(hsc2686['bp_rp'], hsc2686['phot_g_mean_mag'], color='blue', label='HSC_2686', s=4)
        ax.scatter(lynga3['bp_rp'], lynga3['phot_g_mean_mag'], color='red', label='Lynga_3', s=4)
        ax.scatter(czernik['BP-RP'], czernik['Gmag'], color='purple', label='czernik20', s=4)
        #ax.scatter(all_gaia_stars['bp_rp'], all_gaia_stars['phot_g_mean_mag'], 
               #color='grey', label='All stars in DR3', s=0.1)
        ax.invert_yaxis()
        ax.set_xlabel('G_BP - G_RP')
        ax.set_ylabel('G Magnitude')
        ax.set_title(f'Gaia DR3 CMD subplots{i+1}')
        ax.legend()
        plt.show()

def isochrone_changing_Av():

    Av_values = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0]
    Av_files = {Av: f'/Users/xxz/Desktop/LSR-labintern/iso_changing_av/iso_{Av}.dat' for Av in Av_values}
    isochrone_data = {key: readCDSTable(file) for key, file in Av_files.items()}

    unique_logAges = np.unique(isochrone_data[0]['logAge'])

    for age in unique_logAges:
        fig, ax = plt.subplots()
        
        for Av, df in isochrone_data.items():
            # Filter data by age
            new_df = df[np.isclose(df['logAge'], age)]
            if len(new_df) == 0:
                continue
            g_mag = new_df['Gmag']
            bp_rp = new_df['G_BPmag'] - new_df['G_RPmag']
            ax.scatter(bp_rp, g_mag, label=f'Av={Av}', s=2, marker=',')
        
        # Plotting settings
        ax.invert_yaxis()
        ax.set_xlabel('G_BP - G_RP')
        ax.set_ylabel('G Magnitude')
        ax.set_title(f'Isochrones for logAge={age}')
        ax.legend(loc='best')
        plt.show()

def iso_changing_metallicities():
    y




isochrone_changing_Av()

#plot_isochrones_linear(T_iso_000)
#plot_isochrones_log(T_iso_logage)
#plot_isochrones_log(T_iso_logage_2)
#plot_isochrones_log(T_majaess)
#plot_isochrones_log(T_ash)