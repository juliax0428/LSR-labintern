import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia

#np.loadtxt('/Users/xxz/Desktop/J_A+A_673_A114/clusters.dat', delimiter='')

# string = xx:yy:zz.zzz
def hmsToDeg(string):
    try:
        h, m, s = [float(i) for i in string.split(':')]
        return (15. * s / 3600.) + (15. * m / 60.) + (h * 15.)
    except:
        print('hmsToDeg: string = <'+string+'>')
        STOP

def degToHMS(degrees):
    h = int(degrees / 15.)
    m = int((degrees - (h * 15.)) * 4.)
    s = (degrees - (m/4.) - (h*15.)) * 240.
    sStr = '%.3f' % (s)
    sStr = sStr.zfill(6)
    return '%02d:%02d:%s' % (h,m,sStr)

# string = xx:yy:zz.zzz
def dmsToDeg(string):
    try:
        d, m, s = [float(i) for i in string.split(':')]
    #    print('dmsToDeg: string = <'+string+'>: d = ',d,', m = ',m,', s = ',s)
        if string[0] == '-':
            d = 0. - d
            return 0. - (s / 3600. + m / 60. + d)
        return s / 3600. + m / 60. + d
    except:
        return float(string)

def degToDMS(degrees):
    d = int(degrees)
    m = int((degrees - d) * 60)
    s = (degrees - d - (m/60.)) * 3600.
    sStr = '%.3f' % abs(s)
    sStr = sStr.zfill(6)
    return '%02d:%02d:%s' % (d,abs(m),sStr)


def readGaiaMainTable():
    Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"  # Reselect Data Release 3, default

    Gaia.ROW_LIMIT = 3  # Ensure the default row limit.
    coord = SkyCoord(ra=hmsToDeg('13:15:18.74'), dec=dmsToDeg('-65:55:01.2'), unit=(u.degree, u.degree), frame='icrs')
    j = Gaia.cone_search_async(coord, radius=u.Quantity(1.0, u.deg))
    r = j.get_results()
    r.pprint(max_lines=3, max_width=530)

    gaiadr3_table = Gaia.load_table('gaiadr3.gaia_source')
    for column in gaiadr3_table.columns:
        print(column.name)
    return gaiadr3_table

def showGaiaTables():
    tables = Gaia.load_tables(only_names=True)
    for table in tables:
        print(table.get_qualified_name())


def readGaia_astrophysical_parameters():
    """ THIS ONE """
    table = Gaia.load_table('gaiadr3.astrophysical_parameters')
#    for column in table.columns:
#        print(column.name)
    return table

def readCDSTable(tableName,ReadMeName):
    from astropy.io import ascii
    table = ascii.read(tableName,
            readme=ReadMeName)
    return table

tableName = '/Users/xxz/Desktop/J_A+A_673_A114/clusters.dat'
ReadMeName = '/Users/xxz/Desktop/J_A+A_673_A114/ReadMe'

clusters = readCDSTable(tableName=tableName,ReadMeName=ReadMeName)
print('clusters = ',clusters)
print('dir(clusters) = ',dir(clusters))

