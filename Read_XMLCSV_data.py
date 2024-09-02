import pandas as pd
from astropy.io.votable import parse
from astropy.io.votable import parse_single_table

# Example usage
file_path = '/Users/xxz/Desktop/LSR-labintern/Majaess.tsv'  # Update this with the correct path

# Manual processing of the file to handle complex formatting issues
cleaned_data = []
start_of_data = 0

# Read the file line by line and manually split by the correct delimiter, handling any inconsistencies
with open(file_path, 'r') as file:
    for line in file.readlines()[start_of_data:]:
        # Handle possible quotation and delimiter issues by manually processing each line
        if line.strip() and not line.startswith(('<!--', '<', '--')):  # Skip XML-like or comment lines
            # Split the line by semicolon, the presumed delimiter
            split_line = line.replace('\n', '').split(';')
            # Further cleaning if necessary, such as stripping quotes or additional whitespace
            cleaned_line = [item.strip('" ').replace('"', '') for item in split_line]
            cleaned_data.append(cleaned_line)

# Convert list of data into DataFrame
czernik_data = pd.DataFrame(cleaned_data, columns=[
    'RA_ICRS', 'DE_ICRS', 'Source', 'e_RA_ICRS', 'e_DE_ICRS', 'Plx', 'e_Plx', 'PM', 'pmRA', 'e_pmRA', 'pmDE', 'e_pmDE', 
    'RUWE', 'FG', 'e_FG', 'Gmag', 'FBP', 'e_FBP', 'BPmag', 'FRP', 'e_FRP', 'RPmag', 'BP-RP', 'RV', 'e_RV', 'Vbroad', 
    'GRVSmag', 'QSO', 'Gal', 'NSS', 'XPcont', 'XPsamp', 'RVS', 'EpochPh', 'EpochRV', 'MCMCGSP', 'MCMCMSC', 'And', 
    'Teff', 'logg', '[Fe/H]', 'Dist', 'A0', 'HIP', 'PS1', 'SDSS13', 'SKYM2', 'TYC2', 'URAT1', 'AllWISE', 'APASS9', 
    'GSC23', 'RAVE5', '2MASS', 'RAVE6', 'RAJ2000', 'DEJ2000'])

# Try to convert applicable columns to float, while handling errors by keeping original strings if conversion fails
for column in czernik_data.columns:
    try:
        czernik_data[column] = pd.to_numeric(czernik_data[column], errors='ignore')
    except Exception as e:
        print(f"Failed to convert {column}: {e}")

# Display the cleaned DataFrame and data types
czernik_data.head(), czernik_data.dtypes
