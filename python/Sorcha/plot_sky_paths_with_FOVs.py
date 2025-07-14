import rebound
import numpy as np
import matplotlib.pyplot as plt
from auxiliary_functions import kepler, ecc2true
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u

"""
Plots the apparent geocentric trajectory of an object obtained via integration with Rebound,
along with CSS detections and survey fields.

The user must specify:
- the name of the folder containing the results,
- the index of the object in the 'detected_objects.txt' file,
- and the dynamical model:
    full_SS = 1 includes all planets,
    full_SS = 0 includes only the Sun, Earth, and Moon.
"""

folder = 'results_1000'
choice = 866 # object we want to visualize
full_SS = 0


# ---------------------------------------------------------------------------
# constants
au = 149597870700.
mu = 1.32712440018e20
godina = 2*np.pi * np.sqrt(au**3/mu)


detections_file = folder + '/detections' + '.txt'
detected_objects_file = folder + '/detected_objects' + '.txt'
population_file = folder + '/detected_objects' + '.txt'

objectID = np.loadtxt(detections_file, usecols = (0), unpack = True, skiprows = 3, dtype = str)
chosen_object_ID = np.loadtxt(detected_objects_file, usecols = (0), unpack = True, skiprows = 3, dtype = str)

fovID_det, ra_det_1, ra_det_2, dec_det_1, dec_det_2  = np.loadtxt(detections_file, usecols = (1, 4, 6, 8, 10), unpack = True, skiprows = 3)
fovID_det = fovID_det.astype(int)


ISO = chosen_object_ID[choice]
fovID_det = fovID_det[objectID == ISO]
ra_det_1 = ra_det_1[objectID == ISO]
dec_det_1 = dec_det_1[objectID == ISO]

ra_det_2 = ra_det_2[objectID == ISO]
dec_det_2 = dec_det_2[objectID == ISO]

fovID, mjd1, mjd2, ra_fov, dec_fov = np.loadtxt('cadence_file/css.txt', usecols = (0, 1, 2, 4, 5), unpack = True)
fovID = fovID.astype(int)

indices = np.where(np.isin(fovID, fovID_det))[0]

ra_fov = ra_fov[indices] 
dec_fov = dec_fov[indices] 

run = ISO.split('_')[1]
        
OIF_file = folder + '/output' + '_' + run + '.csv'

oifID = np.loadtxt(OIF_file, skiprows = 1, usecols = (0), unpack = True, dtype = str, delimiter = ',')
oifMJD = np.loadtxt(OIF_file, skiprows = 1, usecols = (2), unpack = True, delimiter = ',')
oifMJD_TDB = np.loadtxt(OIF_file, skiprows = 1, usecols = (3), unpack = True, delimiter = ',')

oifMJD = oifMJD[oifID == ISO]
oifMJD_TDB = oifMJD_TDB[oifID == ISO]

delta_MJD = np.mean(oifMJD_TDB - oifMJD) - 2400000.

populationID = np.loadtxt(population_file, skiprows = 3, usecols = (0), unpack = True, dtype = str)
q, ecc, inc, node, argperi, t_p = np.loadtxt(population_file, skiprows = 3, usecols = (1, 2, 3, 4, 5, 6), unpack = True)
q = q[populationID == ISO][0]
ecc = ecc[populationID == ISO][0]
inc = inc[populationID == ISO][0]
node = node[populationID == ISO][0]
argperi = argperi[populationID == ISO][0]
t_p = t_p[populationID == ISO][0]

mjd_initial = min(mjd1)

a_iso = q/(1-ecc) * au

mean_motion = np.sqrt(-mu/a_iso**3)

mean_anomaly = mean_motion * (mjd_initial - t_p) * 86400

E = kepler(ecc, mean_anomaly, 1e-8)

f = ecc2true(E, ecc)

T = (mjd_initial - t_p) * 86400 / godina * 2*np.pi

# Ulazni podaci
korak = 0.1/365 # korak u godinama
total_time= (np.max(mjd2) - np.min(mjd1))/365.25 # ukupno vreme u godinama
#epoha = "JD2456200.5000" # pocetna epoha
epoha = "JD24" + str(mjd_initial + delta_MJD)

#---------------------------------------

sim=rebound.Simulation()
sim.add("Sun", date=epoha)
sim.add("399", date=epoha)
sim.add("301", date=epoha)
if full_SS == 1:
    sim.add("venus", date=epoha)
    sim.add("mercury", date=epoha)
    sim.add("mars", date=epoha)
    sim.add("jupiter", date=epoha)
    sim.add("saturn", date=epoha)
    sim.add("uran", date=epoha)
    sim.add("neptun", date=epoha)

sim.move_to_com()        # We always move to the center of momentum frame before an integration

sim.integrator = "ias15"

sim.add(primary=sim.particles[0],a = q / (1 - ecc), e = ecc, inc = np.deg2rad(inc), 
        omega = np.deg2rad(argperi), Omega = np.deg2rad(node), f = f)


year = 2.*np.pi # One year in units where G=1
times = np.arange(0.,total_time*year, korak*year)
timesMJD = times /2/np.pi * godina / 86400 + mjd_initial

x_gc, y_gc, z_gc = [], [], [] # cartesian geocentric coordinates

for i, t in enumerate(times):
    sim.integrate(t)
    try:
        sim.integrate(t)

    except rebound.Encounter as error:
        print(t/np.pi/2, error)
        
    x_gc.append(sim.particles[-1].x - sim.particles[1].x)
    y_gc.append(sim.particles[-1].y - sim.particles[1].y)
    z_gc.append(sim.particles[-1].z - sim.particles[1].z)
    

ra_reb, dec_reb = [], [] # ra, dec from rebound integration

x_gc = np.array(x_gc)
y_gc = np.array(y_gc)
z_gc = np.array(z_gc)

lon = np.rad2deg(np.arctan2(y_gc, x_gc))
lat = np.rad2deg(np.arctan(z_gc / np.sqrt(x_gc**2 + y_gc**2)))

for i in range(len(x_gc)):
    
    observation_time = Time(timesMJD[i], format='mjd')
    
    gc = SkyCoord(lon[i] * u.degree, lat[i] * u.degree, frame='geocentricmeanecliptic', obstime = observation_time)
    icrs = gc.icrs

    ra_reb.append(icrs.ra.value)
    dec_reb.append(icrs.dec.value)
    
    
def wrap_ra(ra_deg):
    """
    Wraps RA values intelligently:
    - If data spans around 0°, uses [-180°, 180°]
    - Otherwise uses [0°, 360°]
    """
    ra_deg = np.asarray(ra_deg)
    ra_mod = ra_deg % 360

    # If values span the 0/360 boundary (e.g. 359, 0, 1), use [-180, 180]
    if np.ptp(ra_mod) > 180:
        return (ra_mod + 180) % 360 - 180  # wrap to [-180, 180]
    else:
        return ra_mod  # wrap to [0, 360]

ra_det_1 = wrap_ra(ra_det_1)
ra_det_2 = wrap_ra(ra_det_2)
ra_fov = wrap_ra(ra_fov)
ra_reb = wrap_ra(ra_reb)
#for i in range(len(ra_det_1)):
#    if ra_det_1[i]>180:
#        ra_det_1[i]=ra_det_1[i]-360
#        
#
#        
#for i in range(len(ra_det_2)):
#    if ra_det_2[i]>180:
#        ra_det_2[i]=ra_det_2[i]-360
#        
#ra_det_1 = (ra_det_1 + 180) % 360 - 180
#ra_det_2 = (ra_det_2 + 180) % 360 - 180
#        
#        
#for i in range(len(ra_fov)):
#    if ra_fov[i]>180:
#        ra_fov[i]=ra_fov[i]-360
#        
#for i in range(len(ra_reb)):
#    if ra_reb[i]>180:
#        ra_reb[i]=ra_reb[i]-360
        
plt.figure()
plt.plot(ra_reb, dec_reb, 'b', linewidth = 3, label = 'apparent path (rebound)')
plt.plot(ra_det_1, dec_det_1, 'or', markersize = 6, label = 'CSS observations')
plt.plot(ra_det_2, dec_det_2, 'or', markersize = 6)
plt.xlabel('RA (deg)', fontsize = 20)
plt.ylabel('DEC (deg)', fontsize = 20)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)

plt.grid()

d_fov = 1.43

for i in range(len(ra_fov)):
    
    left_x = []
    left_y = []
    right_x = []
    right_y = []
    for j in range(-100,101):
         
        left_x.append(ra_fov[i] - d_fov  / np.cos(np.deg2rad(dec_fov[i]-d_fov /100*j)))
        left_y.append(dec_fov[i] - d_fov /100*j)
        
        right_x.append(ra_fov[i] + d_fov  / np.cos(np.deg2rad(dec_fov[i]-d_fov /100*j)))
        right_y.append(dec_fov[i] - d_fov /100*j)
    
    if i == 0:
        plt.plot(left_x, left_y, '--k', label = 'CSS FOVs')
        plt.plot(right_x, right_y, '--k')
        plt.plot([left_x[0], right_x[0]], [left_y[0], right_y[0]], '--k')
        plt.plot([left_x[-1], right_x[-1]], [left_y[-1], right_y[-1]], '--k')
    else:
        plt.plot(left_x, left_y, '--k')
        plt.plot(right_x, right_y, '--k')
        plt.plot([left_x[0], right_x[0]], [left_y[0], right_y[0]], '--k')
        plt.plot([left_x[-1], right_x[-1]], [left_y[-1], right_y[-1]], '--k')
plt.legend(fontsize = 20)