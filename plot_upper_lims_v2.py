import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time
from astropy.cosmology import WMAP9 as cosmo
from astropy.io import ascii as ascii
from astropy.table import Table, Row, Column
#from frb import frb
#from frb.galaxies import frbgalaxy
from scipy.interpolate import interp1d
import math
from scipy.interpolate import UnivariateSpline

plt.rcParams.update({'font.size': 30})
fig = plt.figure(figsize=(18,10))
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.105, bottom=0.11, right=0.96, top=0.85, wspace=0.2, hspace=0.2)

t0 = +1.

# LAST EPOCH AS TEMPLATE
days_180924_last = [249+1]
days_181112_last = [200+1]
days_190102_last = [149+1]
days_190608_last = [0.20+1,60.7+1]
days_190611_last = [0.49+1,58.5+1]
days_190711_last = [0.07+1,29.0+1]
days_190714_last = [0.48+1,25.7+1]
days_191001_last = [4.14+1,55.3+1]

# apparent magnitude with last epoch as template
apmag_180924_last = [22.5]
apmag_181112_last = [21.4]
apmag_190102_last = [21.8]
apmag_190608_last = [21.7,20.6]
apmag_190611_last = [21.1,21.4]
apmag_190711_last = [22.1,20.6]
apmag_190714_last = [21,20.4]
apmag_191001_last = [23.7,21.8]

# FIRST EPOCH AS TEMPLATE
days_180924_1st = [644+1]
days_181112_1st = [595+1]
days_190102_1st = [545+1]
days_190608_1st = [60.7+1,387+1]
days_190611_1st = [58.5+1,385+1]
days_190711_1st = [29.0+1,355+1]
days_190714_1st = [25.7+1,352+1]
days_191001_1st = [55.3+1,272+1]

# apparent magnitude with first epoch as template
apmag_180924_1st = [21.9]
apmag_181112_1st = [21.2]
apmag_190102_1st = [22.3]
apmag_190608_1st = [21.6,21.7]
apmag_190611_1st = [21.9,22]
apmag_190711_1st = [20.5,21.6]
apmag_190714_1st = [20.8,22.1]
apmag_191001_1st = [22,21.2]



# Distance modulus
distmod_180924 = 41.154939
distmod_181112 = 42.161412
distmod_190102 = 40.908352
distmod_190608 = 38.719884
distmod_190611 = 41.567514
distmod_190711 = 42.405241
distmod_190714 = 40.391713
distmod_191001 = 40.365562

# Galactic extinction from extinction calculator ned ipac
gext_180924 = 0.044
gext_181112 = 0.048
gext_190102 = 0.521
gext_190608 = 0.106
gext_190611 = 0.523
gext_190711 = 0.323
gext_190714 = 0.141
gext_191001 = 0.066

# Host extinction                                                            ######################## EBV=SDSSr?
#hext_180924 = frbgalaxy.FRBHost.by_name('180924').derived['EBV_photom']
#hext_181112 = frbgalaxy.FRBHost.by_name('181112').derived['EBV_photom']
#hext_190102 = frbgalaxy.FRBHost.by_name('190102').derived['EBV_photom']
#hext_190608 = frbgalaxy.FRBHost.by_name('190608').derived['EBV_photom']
#hext_190611 = frbgalaxy.FRBHost.by_name('190611').derived['EBV_photom']
#hext_190711 = frbgalaxy.FRBHost.by_name('190711').derived['EBV_photom']
#hext_190714 = frbgalaxy.FRBHost.by_name('190714').derived['EBV_photom']
#hext_191001 = frbgalaxy.FRBHost.by_name('191001').derived['EBV_photom']


# absolute magnitudes last
abmag_180924_last = np.array(apmag_180924_last) - distmod_180924 - gext_180924 
abmag_181112_last = np.array(apmag_181112_last) - distmod_181112 - gext_181112 
abmag_190102_last = np.array(apmag_190102_last) - distmod_190102 - gext_190102 
abmag_190608_last = np.array(apmag_190608_last) - distmod_190608 - gext_190608 
abmag_190611_last = np.array(apmag_190611_last) - distmod_190611 - gext_190611 
abmag_190711_last = np.array(apmag_190711_last) - distmod_190711 - gext_190711  
abmag_190714_last = np.array(apmag_190714_last) - distmod_190714 - gext_190714 
abmag_191001_last = np.array(apmag_191001_last) - distmod_191001 - gext_191001 

abmag_180924_1st = np.array(apmag_180924_1st) - distmod_180924 - gext_180924 
abmag_181112_1st = np.array(apmag_181112_1st) - distmod_181112 - gext_181112 
abmag_190102_1st = np.array(apmag_190102_1st) - distmod_190102 - gext_190102 
abmag_190608_1st = np.array(apmag_190608_1st) - distmod_190608 - gext_190608 
abmag_190611_1st = np.array(apmag_190611_1st) - distmod_190611 - gext_190611 
abmag_190711_1st = np.array(apmag_190711_1st) - distmod_190711 - gext_190711 
abmag_190714_1st = np.array(apmag_190714_1st) - distmod_190714 - gext_190714 
abmag_191001_1st = np.array(apmag_191001_1st) - distmod_191001 - gext_191001 





# SUPERNOVAE

####################
# SLSN = sn2015bn      explosion was 92 days before maximum light       # paper mail giuliano
# SN2010md            ????
# PTF12dam            day of explosion 56017
# https://watermark.silverchair.com/stv1360.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAqYwggKiBgkqhkiG9w0BBwagggKTMIICjwIBADCCAogGCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMzEWVGssUdJxaha5oAgEQgIICWa-1Mw994iI0MkG6iZj2dxz29jsEa6098hieU-ya7mzYXZMb2F3yp-KtmQnCbcESdPspfgMaGgL77eN7zl2vAzQsGF57g1Yxky9wSTyaPStMTFzJI-HphsqwbTB28yU37nVwmEiLmG7dM_Fr8dBnZjIn3O9MNJ8Uzbt6Y4WEa3-j3dKralWABSWQG-uxIyuJbGa3MjbMa7hcZLkOI5MwSncjj82lVvXaBWe17T_IDqea05fbyUjq2VQASzKgjgW3YMTf5Q1o2TV-PfXImltPDDXCO8pRaAnCj2FtgdND3_w7v7g9PmJXaDFS_Pf3h0CyNtM01Gstj4sZQbOkG8_ANvA1DNskZD-sGDbpJgYNrJdYK-guFy3dO33c8ruytv3M7GCWp64scVfLvbdcPHvZGhDHzHyOOnwlyUZEgw25T3bH-qjuccB5-W7JsD6sq49CzvDwKE4BXyoR3epTrpMyidnvOTzjbDfX1OGARaXPK5BLQpypcOwxDeeQxKCK0SRyQWLzg-z-AArfCdvb77VsHnHlUGR8J_jl8CmQCTfocUQfO8643haR1v95GLdYRJ3DyrJA8H8cIGuViiVSdXNaxX1BfY6kTd4L11cIOTNTIQ39ley9ofFjymKnTYK4I8BHDe_YWmXVkUkXiK4png1GdYV9CNPFN79XjeZ2DPOxgBMzJckWSMxW6XO4sTtAGvRqkCQlXeHXQj5U75gyZ8jhjcJhT_FYMt1R4Xt6pLLtoM8Ab9wNzfw_C0wgbHmfMk34Rfhxs331A7xrdaY7KeQ5uUvYgmTwB1lzBbY

# kilonova  from NGC4993 galaxy               exploto 17 de abril, primer dato 18 de abril

# Ia type SN = SN2014J    ~21 de enero explotó => primer dato 23 de enero => dia 2 primer dato 
# Ia type SN = SN2013gy        day mjd explosion 56629.45
# LSQ12gdj Ia type       https://watermark.silverchair.com/stu1723.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAqMwggKfBgkqhkiG9w0BBwagggKQMIICjAIBADCCAoUGCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMJmwHb5thEEPJwnjTAgEQgIICVoENLvilEAPmQ4QgMw9ed_0471NI8SrmaaMV99dtMUTxHlGQabYOrTWkuz5OIoeHwiHnHyn4gKZLQF7VAYHx4xJrhztMApV4h0-KMpgVGRZET4JM5gqRPRub5jgMOKvYbNIhAyD0-d_GrDNeGyf2m-9ACnXKKAHOyV_EjZknmpF57jfJom1a0EsJPb-nHpskZroO7TUrJuLy5hfTFAo8o-L3tcFYwC7RzH0rbVXQXmY5Mm6XypRVld81H_9c0b2ELa2guhYQS1q8t6ub4--Hs39FXYb0X84c24eeiN5ekjsblx7IzI1LDPQoKrD45N7iC6PnTHMkEuPGN7AXGGzUEWCjrMKRWBIzlqDF2bG7F0St2Gl6x7nUCiDkybiXtsvZUgPrtO-oTQh2u5RK0Ex7Kt_upKSv7061LCZgOhhfzjMmwXiwRMHfmmADU3-VMvK45xDu572NGTL21cmAoQZxG6Ld6Imhhy02_uLEQxHMRDU1NuVHa5X3sW2mrMGdqDLKWy0VrppSYtEUrKsUGLhVnagVfIO9iY2kSEUPF-qrGLXe1lFqOw53SOG0opgtROzO6oEwHy6wgQBgGesAZSSvibRuBGlBDKFbj_alF_MNa7yejhjpm5ilmaPh3p5S1TrQAgfhCp_D9wWqhnO-xQ1IfhZ64ybBmzwM_TOtvSGdJbru8BJwt-JCWlpxc6CUhrhsV8KTuCCaPmZkbb-VhR2UrtGpkYW1O39wKDQD0LSkcg2127JNLXrAhUOJ5_oJS0lqhsVZXm6iUyWrxKlHQ_ZtQUVbSj7PtSM

# https://www.nasa.gov/chandra/multimedia/supernova-sn2014j.html
# IIP type SN = SN2016X      day of explosion mjd=57406.0 (19 de enero)   =>   primer dato 21 de enero => dia 2 primer dato
#https://arxiv.org/pdf/1806.03855.pdf

# GRB SN = sn1998bw         bien centrada         https://arxiv.org/pdf/astro-ph/9808086.pdf
# GRB SN_2 = SN2010bh         day of explosion 55271.5 



####################

# el día más luminoso es el 49 (no centrado en la explosion) => a eso le resto lo 92 dias 
# 49 - 92 = -43
SLSN = Table()
time_mjd_SLSN = Time([57053.5,57072.59,57073.88,57077.17,57077.62,57080.08,57081.5,57082.94,57085.22,57086.01,57090.97,57093.02,57094.87,57095.93,57098.61,57102.52,57106.44,57110.16,57116.93,57120.15,57122.85,57123.87,57127.72,57129.94,57131.71,57132.93,57135.44,57135.94,57138.86,57138.88,57144.71,57144.89,57148.31,57150.03,57150.96,57155.76,57161.7,57162.73,57162.9,57166.92,57173.38,57176.95,57179.11,57184.96,57192.71,57199.38,57206.34,57212.39,57217.73,57230.35,57236.71,57373.33], format='mjd')
SLSN['time_mjd_SLSN'] = time_mjd_SLSN
SLSN['time_iso_SLSN'] = time_mjd_SLSN.iso
SLSN['days_SLSN_nocentered'] = time_mjd_SLSN - SLSN['time_mjd_SLSN'][0]
SLSN['days_SLSN_centered'] = SLSN['days_SLSN_nocentered'] + 43
SLSN['app_mag_SLSN'] = [17.33,16.98,17.02,16.95,16.92,16.94,17,16.92,16.95,16.9,16.82,16.78,16.75,16.77,16.8,16.72,16.84,16.85,16.91,16.99,17.01,17.12,17.21,17.16,17.30,17.20,17.28,17.26,17.29,17.33,17.5,17.41,17.29,17.43,17.39,17.44,17.66,17.66,17.56,17.68,17.75,17.9,17.88,18.06,18.12,18.23,18.17,18.31,18.37,18.42,18.53,20.5]
gext_SLSN = 0.057
z_SLSN = 0.1136
#K_SLSN = 2.5 * math.log10(1+z_SLSN)
SLSN['ab_mag_SLSN'] = SLSN['app_mag_SLSN'] - 38.63 - gext_SLSN
ab_mag_SLSN = SLSN['ab_mag_SLSN']
days_SLSN_centered = SLSN['days_SLSN_centered'].value + t0

# SLSN
PTF_ascii = ascii.read('/home/consuelo/projects/FRBs/LCOgt/SN/PTF12dam.csv')
cond = np.where(PTF_ascii['band']=='r')
PTF_ascii = PTF_ascii[cond]
time_mjd_PTF = PTF_ascii['time']
time_mjd_PTF = Time(time_mjd_PTF, format='mjd')
day_of_explosion_PTF = Time(56017, format='mjd')
PTF = Table()
PTF['time_mjd_PTF'] = time_mjd_PTF
PTF['time_iso_PTF'] = time_iso = time_mjd_PTF.iso
#PTF['days_PTF_nocentered'] = time_mjd_PTF - PTF['time_mjd_PTF'][0]
PTF['days_PTF_centered'] = PTF['time_mjd_PTF'] - day_of_explosion_PTF
z_PTF = 0.107
#K_PTF = 2.5 * math.log10(1+z_PTF)
dist_mod_PTF = cosmo.distmod(z_PTF).value
PTF['ab_magnitude'] = PTF_ascii['magnitude'] - dist_mod_PTF 
days_PTF_centered = PTF['days_PTF_centered'].value + t0
ab_mag_PTF = PTF['ab_magnitude']


# SLSN_II
SN2010md_ascii = ascii.read('/home/consuelo/projects/FRBs/LCOgt/SN/SN2010md.csv')
cond = np.where(SN2010md_ascii['band']=='r')
SN2010md_ascii = SN2010md_ascii[cond]
time_mjd_SN2010md = SN2010md_ascii['time']
time_mjd_SN2010md = Time(time_mjd_SN2010md, format='mjd')
SN2010md = Table()
SN2010md['time_mjd_SN2010md'] = time_mjd_SN2010md
SN2010md['time_iso_SN2010md'] = time_iso_SN2010md = time_mjd_SN2010md.iso
SN2010md['days_SN2010md_nocentered'] = time_mjd_SN2010md - SN2010md['time_mjd_SN2010md'][0]
SN2010md['days_SN2010md_centered'] = SN2010md['days_SN2010md_nocentered'] 
z_SN2010md = 0.0982
#K_SN2010md = 2.5 * math.log10(1+z_SN2010md)
dist_mod_SN2010md = cosmo.distmod(z_SN2010md).value
gext_SN2010md = 0.193
SN2010md['ab_magnitude'] = SN2010md_ascii['magnitude'] - dist_mod_SN2010md - gext_SN2010md 
days_SN2010md_centered = SN2010md['days_SN2010md_centered'].value + t0
ab_mag_SN2010md = SN2010md['ab_magnitude'] 


GRB_SN = Table()
time_mjd_GRB_SN = Time([50930.13,50930.15,50930.2,50930.30,50931.91,50932.18,50932.21,50934.82,50936.83,50937.90,50939.82,50940.27,50941.8,50942.26,50942.80,50943.28,50944.75,50946.89,50948.89,50950.71,50950.88,50951.84,50952.8,50953.73,50954.09,50954.84,50954.87,50955.19,50955.85,50955.87,50956.79,50959.16,50960.12,50961.65,50962.07,50964.11,50965.79,50966.86,50967.81,50968.80,50969.05,50970.21,50970.86,50971.88,50975.05,50981.03,50983.57,50984.62,50985.68,50986.03,51255.4,51274.4,51319.4,51346.4,51438.4],format='mjd')
GRB_SN['time_mjd_GRB_SN'] = time_mjd_GRB_SN.value
GRB_SN['time_iso_GRB_SN'] = time_mjd_GRB_SN.iso
GRB_SN['days_GRB_SN_centered'] = GRB_SN['time_mjd_GRB_SN'] - GRB_SN['time_mjd_GRB_SN'][0]
GRB_SN['app_mag_GRB_SN'] = [15.76,15.79,15.96,15.76,15.71,15.61,15.61,14.86,14.43,14.24,14.04,14,13.92,13.87,13.88,13.84,13.81,13.75,13.79,13.83,13.83,13.85,13.89,13.92,13.96,13.97,13.93,13.96,14.02,14.01,14.06,14.18,14.23,14.39,14.32,14.41,14.53,14.59,14.62,14.69,14.71,14.76,14.79,14.83,14.94,15.13,15.22,15.27,15.28,15.21,19.742,20.090,20.830,20.87,21.94]
gext_GRB = 0.134
z_GRB = 0.0085
#K_GRB = 2.5 * math.log10(1+z_GRB)
GRB_SN['ab_mag_GRB_SN'] = GRB_SN['app_mag_GRB_SN'] - 33.4 - gext_GRB 
days_GRB_SN_centered = GRB_SN['days_GRB_SN_centered'] + t0
ab_mag_GRB_SN = GRB_SN['ab_mag_GRB_SN']

GRB_SN_2 = Table()
GRB_SN_2['days_GRB_SN_2_centered'] = [1.49, 1.58, 2.48, 2.7, 3.5, 4.5, 5.48, 6.47, 8.48, 11.48, 12.47, 15.48, 16.48, 20.53, 22.49, 26.5, 33.49, 41.48, 84.44]
GRB_SN_2['app_mag_GRB_SN_2'] = [20.9, 20.91, 20.9, 20.77, 20.58, 20.4, 20.21, 20.11, 19.96, 19.91, 20.01, 20.07, 20.07, 20.36, 20.48, 20.78, 21.42, 21.83, 23.66]
gext_GRB_2 = 0.266
z_GRB_2 = 0.0593
#K_GRB_2 = 2.5 * math.log10(1+z_GRB_2)
GRB_SN_2['ab_mag_GRB_SN_2'] = GRB_SN_2['app_mag_GRB_SN_2'] - 37.02 - gext_GRB_2 
days_GRB_SN_2_centered = GRB_SN_2['days_GRB_SN_2_centered'] + t0
ab_mag_GRB_SN_2 = GRB_SN_2['ab_mag_GRB_SN_2']


kilonova = Table()
time_mjd_kilonova = Time([57983.75833,57983.96875,57984.76111,57984.96892,57985.77639,57985.97433,57986.97426], format='mjd')
time_explosion = Time([57982.52852], format='mjd')
kilonova['time_mjd_kilonova'] = time_mjd_kilonova
kilonova['time_iso_kilonova'] = time_mjd_kilonova.iso
kilonova['days_kilonova_centered'] = time_mjd_kilonova - time_explosion 
#kilonova['days_kilonova_centered'] = kilonova['days_kilonova_nocentered'] + 1

kilonova['ap_mag_kilonova'] = [17.89,17.99,18.8,19.13,19.52,19.81,20.53]
z_kilonova = 0.009843
#K_kilonova = 2.5 * math.log10(1+z_kilonova)
kilonova['ab_mag_kilonova'] = kilonova['ap_mag_kilonova'] - 33.01 - 0.282 
days_kilonova_centered = kilonova['days_kilonova_centered'].value + t0
ab_mag_kilonova = kilonova['ab_mag_kilonova']



# LSQ12gdj Ia type
LSQ_ascii = ascii.read('/home/consuelo/projects/FRBs/LCOgt/SN/LSQ12gdj.csv')
cond = np.where(LSQ_ascii['band']=='r')
LSQ_ascii = LSQ_ascii[cond]
time_mjd_LSQ = LSQ_ascii['time']
time_mjd_LSQ = Time(time_mjd_LSQ, format='mjd')
day_of_explosion_LSQ = Time(56235.9, format='mjd')
LSQ = Table()
LSQ['time_mjd_LSQ'] = time_mjd_LSQ
LSQ['time_iso_LSQ'] = time_iso = time_mjd_LSQ.iso
#LSQ['days_LSQ_nocentered'] = time_mjd_LSQ - LSQ['time_mjd_LSQ'][0]
LSQ['days_LSQ_centered'] = LSQ['time_mjd_LSQ'] - day_of_explosion_LSQ
z_LSQ = 0.03
#K_LSQ = 2.5 * math.log10(1+z_LSQ)
dist_mod_LSQ = cosmo.distmod(z_LSQ).value
gext_LSQ = 0.053
LSQ['ab_magnitude'] = LSQ_ascii['magnitude'] - dist_mod_LSQ - gext_LSQ 
days_LSQ_centered = LSQ['days_LSQ_centered'].value + t0
ab_mag_LSQ = LSQ['ab_magnitude']


# SN2009F Ia type           explosion 17 dias antes del maximo => 54828.15937
SN2009F_ascii = ascii.read('/home/consuelo/projects/FRBs/LCOgt/SN/SN2009f.csv')
cond = np.where(SN2009F_ascii['band']=='r')
SN2009F_ascii = SN2009F_ascii[cond]
time_mjd_SN2009F = SN2009F_ascii['time']
time_mjd_SN2009F = Time(time_mjd_SN2009F, format='mjd')
SN2009F = Table()
SN2009F['time_mjd_SN2009F'] = time_mjd_SN2009F
day_of_explosion_SN2009F = Time(54828.15, format='mjd')
SN2009F['days_SN2009F_centered'] = SN2009F['time_mjd_SN2009F'] - day_of_explosion_SN2009F
z_SN2009F = 0.012909
#K_SN2009F = 2.5 * math.log10(1+z_SN2009F)
dist_mod_SN2009F = cosmo.distmod(z_SN2009F).value
gext_SN2009F = 0.238
SN2009F['ab_magnitude'] = SN2009F_ascii['magnitude'] - dist_mod_SN2009F - gext_SN2009F 
days_SN2009F_centered = SN2009F['days_SN2009F_centered'].value + t0
ab_mag_SN2009F = SN2009F['ab_magnitude']





ms=24
msS=8



frb180924_last, = ax.plot(days_180924_last, abmag_180924_last, 'v',color='#ae3492', markersize=ms, label='FRB180924')
frb180924_1st, = ax.plot(days_180924_1st, abmag_180924_1st, 'v',fillstyle='none',markeredgewidth=2.2, color='#ae3492', markersize=ms)
frb181112_last, = ax.plot(days_181112_last, abmag_181112_last, 'v',color='#74b900', markersize=ms, label='FRB181112')
frb181112_1st, = ax.plot(days_181112_1st, abmag_181112_1st, 'v',fillstyle='none',markeredgewidth=2.2, color='#74b900', markersize=ms)
frb190102_last, = ax.plot(days_190102_last, abmag_190102_last, 'v',color='#ffb500', markersize=ms, label='FRB190102')
frb190102_1st, = ax.plot(days_190102_1st, abmag_190102_1st, 'v',fillstyle='none',markeredgewidth=2.2, color='#ffb500', markersize=ms)
frb190608_last, = ax.plot(days_190608_last, abmag_190608_last, 'v',color='#ff6700', markersize=ms, label='FRB190608')
frb190608_1st, = ax.plot(days_190608_1st, abmag_190608_1st, 'v',fillstyle='none',markeredgewidth=2.2, color='#ff6700', markersize=ms)
frb190611_last, = ax.plot(days_190611_last, abmag_190611_last, 'v',color='#3cbca1', markersize=ms, label='FRB190611')
frb190611_1st, = ax.plot(days_190611_1st, abmag_190611_1st, 'v',fillstyle='none',markeredgewidth=2.2, color='#3cbca1', markersize=ms)
frb190711_last, = ax.plot(days_190711_last, abmag_190711_last, 'v',color='#D47D8E', markersize=ms, label='FRB190711')
frb190711_1st, = ax.plot(days_190711_1st, abmag_190711_1st, 'v',fillstyle='none',markeredgewidth=2.2, color='#D47D8E', markersize=ms)
frb190714_last, = ax.plot(days_190714_last, abmag_190714_last, 'v',color='#df2a3c', markersize=ms, label='FRB190714')
frb190714_1st, = ax.plot(days_190714_1st, abmag_190714_1st, 'v',fillstyle='none',markeredgewidth=2.2, color='#df2a3c', markersize=ms)
frb191001_last, = ax.plot(days_191001_last, abmag_191001_last, 'v',color='#4838b6', markersize=ms, label='FRB191001')
frb191001_1st, = ax.plot(days_191001_1st, abmag_191001_1st, 'v',fillstyle='none',markeredgewidth=2.2, color='#4838b6', markersize=ms)

spl_GRB = UnivariateSpline(days_GRB_SN_centered, ab_mag_GRB_SN)
spl_GRB.set_smoothing_factor(0.2)
xs_GRB = np.linspace(1.65, 80, 100) + t0
#ab_GRB_SN, = ax.plot(days_GRB_SN_centered, ab_mag_GRB_SN, 'o',color='#ffb500', markersize=msS, label='GRB SN')	
ax.plot(xs_GRB,spl_GRB(xs_GRB),color='#ffb500',lw=1)

spl_GRB_2 = UnivariateSpline(days_GRB_SN_2_centered, ab_mag_GRB_SN_2)
spl_GRB_2.set_smoothing_factor(0.05)
xs_GRB_2 = np.linspace(1.65, 80, 100) + t0
#ab_GRB_SN_2, = ax.plot(days_GRB_SN_2_centered, ab_mag_GRB_SN_2, 'o',color='#ffb500', markersize=msS, label='GRB SN')	
ax.plot(xs_GRB_2,spl_GRB_2(xs_GRB_2),color='#ffb500',lw=1)

days_GRB = np.linspace(1.65, 80, 100) + t0
ax.fill_between(days_GRB, spl_GRB(xs_GRB), spl_GRB_2(xs_GRB_2), color='#ffb500', alpha=0.1)


spl_kilo = UnivariateSpline(days_kilonova_centered, ab_mag_kilonova)
spl_kilo.set_smoothing_factor(0.1)
xs_kilo = np.linspace(1, 5, 100) + t0
ab_kilonova, = ax.plot(days_kilonova_centered, ab_mag_kilonova, 'o',color='#ae3492', markersize=msS, label='Kilonova')
ax.plot(xs_kilo,spl_kilo(xs_kilo),color='#ae3492',lw=1)


# SLSN
spl_SLSN = UnivariateSpline(days_SLSN_centered, ab_mag_SLSN)
spl_SLSN.set_smoothing_factor(0.111)
xs_SLSN = np.linspace(6, 400, 100) + t0
#ab_SLSN, = ax.plot(days_SLSN_centered, ab_mag_SLSN, 'o',color='#df2a3c', markersize=msS, label='SL SN')
ax.plot(xs_SLSN,spl_SLSN(xs_SLSN),color='#df2a3c',lw=1)

# spl_PTF = UnivariateSpline(days_PTF_centered, ab_mag_PTF)
# spl_PTF.set_smoothing_factor(0.5)
# xs_PTF = np.linspace(1, 400, 100)
# ab_PTF, = ax.plot(days_PTF_centered, ab_mag_PTF, 'o',color='#df2a3c', markersize=msS, label='SL SN')
# ax.plot(xs_PTF,spl_PTF(xs_PTF),color='#df2a3c',lw=1)

spl_SN2010md = UnivariateSpline(days_SN2010md_centered, ab_mag_SN2010md)
spl_SN2010md.set_smoothing_factor(0.4)
xs_SN2010md = np.linspace(6, 400, 100) + t0
#ab_SN2010md, = ax.plot(days_SN2010md_centered, ab_mag_SN2010md, 'o',color='#df2a3c', markersize=msS, label='SL SN')
ax.plot(xs_SN2010md,spl_SN2010md(xs_SN2010md),color='#df2a3c',lw=1)

days_SLSN = np.linspace(6, 400, 100) + t0
ax.fill_between(days_SLSN, spl_SLSN(xs_SLSN), spl_SN2010md(xs_SN2010md), color='#df2a3c', alpha=0.1) 



# # II type
# spl = UnivariateSpline(days_SN2016X_centered, ab_mag_SN2016X)
# spl.set_smoothing_factor(0.2)
# xs = np.linspace(1.8, 180, 100)
# #ab_SN2016X, = ax.plot(days_SN2016X_centered, ab_mag_SN2016X, 'o',color='#2685b9', markersize=msS, label='IIP type SN')
# ax.plot(xs,spl(xs),color='#2685b9',lw=1)


# I type
spl1 = UnivariateSpline(days_LSQ_centered, ab_mag_LSQ)
spl1.set_smoothing_factor(0.1)
xs1 = np.linspace(10, 70, 42) + t0
#ab_LSQ, = ax.plot(days_LSQ_centered, ab_mag_LSQ, 'o',color='#ffa600', markersize=msS, label='Ia type SN')
ax.plot(xs1,spl1(xs1),color='#74b900',lw=1)

# spl0 = UnivariateSpline(days_SN2004ey_centered, ab_mag_SN2004ey)
# spl0.set_smoothing_factor(0.1)
# xs0 = np.linspace(2, 142, 42)
# ab_SN2004ey, = ax.plot(days_SN2004ey_centered, ab_mag_SN2004ey, 'o',color='#ffa600', markersize=msS, label='Ia type SN')
# ax.plot(xs0,spl1(xs0),color='#ffa600',lw=1)

spl2 = UnivariateSpline(days_SN2009F_centered, ab_mag_SN2009F)
spl2.set_smoothing_factor(0.02)
xs2 = np.linspace(10, 70, 42) + t0
#ab_SN2009F, = ax.plot(days_SN2009F_centered, ab_mag_SN2009F, 'o',color='#ffb500', markersize=msS, label='Ia type SN')
ax.plot(xs2,spl2(xs2),color='#74b900',lw=1)

days_Ia = np.linspace(10, 70, 42) + t0
ax.fill_between(days_Ia, spl2(xs2), spl1(xs1), color='#74b900', alpha=0.2) 

ax.set_xlabel('Time + 1 (days)')
ax.set_ylabel('Absolute Magnitude (mag)')
#ax.set_ylim(-21,-16)
ax.legend(handles=[frb180924_last,frb181112_last,frb190102_last,frb190608_last,frb190611_last,frb190711_last,frb190714_last,frb191001_last], bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left', ncol=4, mode="expand", borderaxespad=0., fontsize='small')
ax.annotate("GRB SN", (4,-17.6), size='small', color='#ffb500')
ax.annotate("SL SN", (82,-20.6), size='small', color='#df2a3c')
#ax.annotate("II type SN", (148,-15), size='small', color='#2685b9')
ax.annotate("Ia type SN", (12.5,-17.98), size='small', color='#74b900')
ax.annotate("Kilonova", (5.8,-13.3), size='small', color='#ae3492')
ax.invert_yaxis()
ax.set_xscale('log')
ax.set_xlim(-200,750)
ax.set_ylim(-12,-22.7)
ax.minorticks_on()
fig.savefig("superplot_v2.pdf")
plt.show()

