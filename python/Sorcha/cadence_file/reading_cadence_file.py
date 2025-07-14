'''
Idnetifikovati sve FOVs u kojima se nalazi centralna tačka DDF 
i dati njihovu raspodelu po kadenci 
(proteklom vremenu izmedju uzastopnih posmatranja centralne tačke datog DDF)

cadence fajlovi: http://astro-lsst-01.astro.washington.edu:8080/
'''

import numpy as np
import pandas as pd
import sqlite3
import matplotlib.pyplot as plt



def spherical_distance(lon1, lat1, lon2, lat2):
    # Convert degrees to radians
    lon1 = np.radians(lon1)
    lat1 = np.radians(lat1)
    lon2 = np.radians(lon2)
    lat2 = np.radians(lat2)
    
    # Spherical law of cosines formula
    delta_lon = lon2 - lon1
    distance = np.arccos(np.sin(lat1) * np.sin(lat2) + np.cos(lat1) * np.cos(lat2) * np.cos(delta_lon))
    
    return np.rad2deg(distance)
    
lsst_r = 1.5 # poluprecnik LSST FOV (deg) zapravo je 1.75 ali sam uzeo samo one FOVs koji su najdelje 1.5 stepen udaljeni od DDF

'''
lista DDFs
https://www.lsst.org/scientists/survey-design/ddf
'''

ddf_names = ['ELAIS S1', 'XMM-LSS', 'Extended Chandra Deep Field-South', 'COSMOS']
ddf_ra = [9.45, 35.25, 53.12, 150.1] # rektascenzije DDFs (deg)
ddf_dec = [-44, -4.5, -27.81, 2.18] # deklinacije DDFs (deg)


'''
BASELINE
'''

dbfile = 'baseline_v4.3.5_10yrs.db'

con = sqlite3.connect(dbfile)
# creating cursor
cur = con.cursor()


#observations = pd.read_sql_query("SELECT * FROM observations", con)
observations = pd.read_sql_query("SELECT observationStartMJD, fieldRA, fieldDec FROM observations", con)
info = pd.read_sql_query("SELECT * FROM info", con)

#

#observationId = np.array(observations['observationId'])
#note = np.array(observations['note'])
#ID_bl = np.array(observations['observationId'])
observationStartMJD_bl= np.array(observations['observationStartMJD'])
#visitTime_bl = np.array(observations['visitTime'])
#visitExposureTime_bl = np.array(observations['visitExposureTime'])
#filter_type_bl = np.array(observations['filter'])
#seeingFwhmGeom_bl = np.array(observations['seeingFwhmGeom'])
#seeingFwhmEff_bl = np.array(observations['seeingFwhmEff'])
#fiveSigmaDepth_bl = np.array(observations['fiveSigmaDepth'])
fieldRA_bl = np.array(observations['fieldRA'])
fieldDec_bl = np.array(observations['fieldDec'])
#rotSkyPos_bl = np.array(observations['rotSkyPos'])


con.close()

'''
OCEAN
'''
dbfile = 'ddf_ocean_ocean6_v4.3.5_10yrs.db'
con = sqlite3.connect(dbfile)
# creating cursor
cur = con.cursor()


#observations = pd.read_sql_query("SELECT * FROM observations", con)
observations = pd.read_sql_query("SELECT observationStartMJD, fieldRA, fieldDec FROM observations", con)
info = pd.read_sql_query("SELECT * FROM info", con)

#

#observationId = np.array(observations['observationId'])
#note = np.array(observations['note'])
#ID_oc = np.array(observations['observationId'])
observationStartMJD_oc= np.array(observations['observationStartMJD'])
#visitTime_oc = np.array(observations['visitTime'])
#visitExposureTime_oc = np.array(observations['visitExposureTime'])
#filter_type_oc = np.array(observations['filter'])
#seeingFwhmGeom_oc = np.array(observations['seeingFwhmGeom'])
#seeingFwhmEff_oc = np.array(observations['seeingFwhmEff'])
#fiveSigmaDepth_oc = np.array(observations['fiveSigmaDepth'])
fieldRA_oc = np.array(observations['fieldRA'])
fieldDec_oc = np.array(observations['fieldDec'])
#rotSkyPos_oc = np.array(observations['rotSkyPos'])

con.close()
#cadence = np.diff(observationStartMJD)

'''
Odredjivanje FOVs u koje upada centar DDF
'''


# LSST FOVs u koje upada centar ddf
ddf_fovs_bl = [[]] * len(ddf_names)
ddf_fovs_oc = [[]] * len(ddf_names)

for i in range(len(ddf_names)):
    
    # odredjivanje uglovnog rastojanja centra ddf od centra FOV
    dist_bl = spherical_distance(ddf_ra[i], ddf_dec[i], fieldRA_bl, fieldDec_bl)
    # selektovanje samo onih FOV u koje upada centar ddf
    ddr_ind_bl = list(np.where(dist_bl < lsst_r)[0])
    ddf_fovs_bl[i] = ddr_ind_bl
    
    dist_oc = spherical_distance(ddf_ra[i], ddf_dec[i], fieldRA_oc, fieldDec_oc)
    # selektovanje samo onih FOV u koje upada centar ddf
    ddr_ind_oc = list(np.where(dist_oc < lsst_r)[0])
    ddf_fovs_oc[i] = ddr_ind_oc
   

'''
Raspodele rektascenzije i deklinacije centralne tacke FOV za svaki DDF
'''    

for i in range(len(ddf_names)):
    
    plt.figure(figsize=(25, 6))  # Adjust figure size for more squared subplots
    plt.suptitle('Raspodela posmatranja po rektascenziji i deklinaciji za \n polje {}'.format(ddf_names[i]), 
                 fontsize=22, fontweight='bold', color='red')  # Increase distance between title and plots
    
    plt.subplots_adjust(wspace=0.8)  # Fine-tune horizontal space between subplots
    
    plt.subplot(131)
    plt.hist(fieldRA_bl[ddf_fovs_bl[i]], 100, alpha = 0.5, label = 'Baseline')
    plt.hist(fieldRA_oc[ddf_fovs_oc[i]], 100, alpha = 0.5, label = 'Ocean')
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.xlabel('rektascenzija (deg)', fontsize=20, fontweight='bold')
    plt.ylabel('broj posmatranja', fontsize=20, fontweight='bold')
    y_limits = plt.ylim()
    plt.plot([ddf_ra[i], ddf_ra[i]],y_limits, 'r', linewidth = 2)
    text_x = ddf_ra[i]-0.1
    text_y = (y_limits[0] + y_limits[1]) / 2  # Position the text in the middle of the y-axis
    plt.text(text_x, text_y, 'Centar DDF', color='red', fontsize=12, ha='center', va='center', rotation=90, fontweight = 'bold')
    plt.legend(fontsize = 16)
    plt.grid()
    
    plt.subplot(132)
    plt.hist(fieldDec_bl[ddf_fovs_bl[i]], 100, alpha = 0.5, label = 'Baseline')
    plt.hist(fieldDec_oc[ddf_fovs_oc[i]], 100, alpha = 0.5, label = 'Ocean')
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.xlabel('deklinacija (deg)', fontsize=20, fontweight='bold')
    plt.ylabel('broj posmatranja', fontsize=20, fontweight='bold')
    y_limits = plt.ylim()
    plt.plot([ddf_dec[i], ddf_dec[i]],y_limits, 'r', linewidth = 2)
    text_x = ddf_dec[i]-0.1
    text_y = (y_limits[0] + y_limits[1]) / 2  # Position the text in the middle of the y-axis
    plt.text(text_x, text_y, 'Centar DDF', color='red', fontsize=12, ha='center', va='center', rotation=90, fontweight = 'bold')
    plt.legend(fontsize = 16)
    plt.grid()
    
#    ## Plotting the bivariate histogram
#    plt.subplot(133)
#    plt.hist2d(fieldRA[ddf_fovs[i]], fieldDec[ddf_fovs[i]], bins=20, cmap='Blues')
#    plt.tick_params(axis='both', which='major', labelsize=14)
#    plt.xlabel('rektascenzija (deg)', fontsize=20, fontweight='bold')
#    plt.ylabel('deklinacija (deg)', fontsize=20, fontweight='bold')
#    plt.plot(ddf_ra[i], ddf_dec[i], 'or', markersize = 10)
#    
#    text_y = ddf_dec[i]+0.2
#    text_x = ddf_ra[i]
#    plt.text(text_x, text_y, 'Centar DDF', color='red', fontsize=12, ha='center', va='center', fontweight = 'bold')
    
#    cbar = plt.colorbar()  # Add colorbar to the last subplot
#    cbar.set_label('broj posmatranja', fontsize=20, fontweight='bold')  # Set label with desired font properties
#    cbar.ax.tick_params(labelsize=14)  # Increase colorbar ticks font size
#    plt.grid()
#    
#    plt.tight_layout(rect=[0, 0, 1, 0.95], pad=3.0)  # Adjust the layout with padding to include space for the title
    
    
    
'''
Odredjivanje vremena pocetka ekspozicije za ddf polja
'''

ddf_mjd_bl = [[]] * len(ddf_names)
ddf_mjd_oc = [[]] * len(ddf_names)
for i in range(len(ddf_names)):
    
    ddf_mjd_bl[i] = observationStartMJD_bl[ddf_fovs_bl[i]]
    ddf_mjd_oc[i] = observationStartMJD_oc[ddf_fovs_oc[i]]
    
    

'''
plotovanje raspodela kadenci
'''
for i in range(len(ddf_names)):
    
    plt.figure()  # Adjust figure size for more squared subplots
    plt.title('raspodela kadence za \n polje {}'.format(ddf_names[i]), 
              fontsize=22, fontweight='bold', color='red')

    plt.hist(np.diff(np.array(ddf_mjd_bl[i])*1440), 
             np.linspace(0, 100, 200), 
             label = 'Baseline', alpha = 0.5)  # Set bar color to black
    
    plt.hist(np.diff(np.array(ddf_mjd_oc[i])*1440), 
             np.linspace(0, 100, 200), 
             label = 'Ocean', alpha = 0.5)  # Set bar color to black
    
    
    plt.yscale('log')  # Set y-axis to log scale
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.xlabel('kadenca (min)', fontsize=20, fontweight='bold')
    plt.ylabel('broj kadenci', fontsize=20, fontweight='bold')
    plt.ylim(bottom=1)
    plt.legend(fontsize = 20)
    plt.grid()


   


