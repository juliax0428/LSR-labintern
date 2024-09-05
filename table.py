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

def save_csvfiles():
    ClusterTableName = '/Users/xxz/Desktop/LSR-labintern/J_A+A_685_A83/clusters.dat'
    MemberTableName_1 = '/Users/xxz/Desktop/LSR-labintern/J_A+A_685_A83/tablec1.dat'
    MemberTableName_2 = '/Users/xxz/Desktop/LSR-labintern/J_A+A_685_A83/tablec2.dat'
    MemberTableName_3 = '/Users/xxz/Desktop/LSR-labintern/J_A+A_685_A83/tablec3.dat'
    ReadMeName = '/Users/xxz/Desktop/LSR-labintern/J_A+A_685_A83/ReadMe'
    clusters = readCDSTable(tableName=ClusterTableName, ReadMeName=ReadMeName)

    i=0
    for cluster in clusters:
        if i >= 5:
            break
        print ('clusters name = ',cluster['Name'])
        print ('RAh', cluster['RAh'], 'RAm', cluster['RAm'],'RAs',cluster['RAs'])
        print ('DE-', cluster['DE-'], 'DEd', cluster['DEd'],'DEm',cluster['DEm'],'DEs',cluster['DEs'])
        print('Number of sources', cluster['N'])
        i += 1
    
save_csvfiles()