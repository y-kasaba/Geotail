"""
    Geotail PWI KPRM lib -- 2025/8/15
"""
import datetime
import math
import numpy as np

# ---------------------------------------------------------
# Geotail PWI KPRM FX (for orbit)
# ---------------------------------------------------------
def read_orbit_multi(Epoch0_min, Epoch0_max, mode_gap, mode_check, data_dir):
    """
        orbit.orbit_SM[epoch][3]   X/Y/Z in SM coordinate (Re)  
        orbit.[epoch]              R (Re)
        orbit.MLAT[epoch]          MLAT (deg)
        orbit.MLT[epoch]           MLT (h)
        orbit.tilt[epoch]          Tilt angle of Magnetic pole (deg)
        orbit.epoch[epoch]         Epoch
    """
    class struct:
        pass
    orbit     = struct()
    orbit.num = 0
    Epoch_min = datetime.datetime.strptime(Epoch0_min[0:10], "%Y-%m-%d")     # Start time
    Epoch_max = datetime.datetime.strptime(Epoch0_max[0:10], "%Y-%m-%d")     # End   time
    Epoch_len = Epoch_max - Epoch_min;   print('[Total]', Epoch_len.days+1, 'days      (', Epoch0_min[0:10], '-', Epoch0_max[0:10],  ')')
    Epoch     = Epoch_min
    for i in range (Epoch_len.days+1):
        str_Epoch = Epoch.strftime('%Y-%m-%d ')
        orbit1 = read_orbit(str_Epoch, data_dir)
        if orbit1.num>0:
            if orbit.num==0 or mode_check == 1:
                orbit = orbit1
                orbit.epoch = np.array(orbit1.epoch)
            else:
                orbit = orbit_add(orbit, orbit1)
        Epoch  = Epoch + datetime.timedelta(days=1)
        print(i, str_Epoch, orbit.num)

    # Fill NAN into GAP
    if mode_gap == 1:
        for i in range(orbit.num-1):
            dt = orbit.epoch[i+1] - orbit.epoch[i]
            if dt.total_seconds() > 60:
                orbit_nan(orbit, i);   orbit_nan(orbit, i+1)
    """
    print("[orbit_SM] ",    orbit.orbit_SM[0], orbit.orbit_SM[-1],    orbit.orbit_SM.shape,  orbit.orbit_SM.dtype)
    print("[R] ",           orbit.R[0],        orbit.R[-1],           orbit.R.shape,         orbit.R.dtype)
    print("[MLAT] ",        orbit.MLAT[0],     orbit.MLAT[-1],        orbit.MLAT.shape,      orbit.MLAT.dtype)
    print("[MLT] ",         orbit.MLT[0],      orbit.MLT[-1],         orbit.MLT.shape,       orbit.MLT.dtype)
    print("[tilt] ",        orbit.tilt[0],     orbit.tilt[-1],        orbit.tilt.shape,      orbit.tilt.dtype)
    print("[epoch] ",       orbit.epoch[0],    orbit.epoch[-1],       orbit.epoch.shape,     orbit.epoch.dtype)
    """
    return  orbit


def read_orbit(Epoch, data_dir):
    class struct:
        pass
    orbit1     = struct()
    orbit1.num = 0
    Epoch_YY = Epoch[2:4];  Epoch_MM = Epoch[5:7];  Epoch_DD = Epoch[8:10]
    name_FX_file  = name_FX_data(Epoch_YY, Epoch_MM, Epoch_DD, data_dir)
    try:
        with open(name_FX_file, 'r') as f:
            f.close()
    except FileNotFoundError:
        print("***ERROR*** ORB file - not found:", name_FX_file)
        return orbit1

    # Decode ORBIT
    f = open(name_FX_file, 'r')
    line_FX = f.readlines()
    f.close()
    line_FX = line_FX[1:len(line_FX)]
    #
    orbit1.num      = len(line_FX)
    orbit1.epoch    = []
    orbit1.orbit_SM = np.zeros((orbit1.num, 3))
    orbit1.R        = np.zeros((orbit1.num))
    orbit1.MLAT     = np.zeros((orbit1.num))
    orbit1.MLT      = np.zeros((orbit1.num))
    orbit1.tilt     = np.zeros((orbit1.num))

    for i in range(orbit1.num):
        # Epoch
        time_str = Epoch[0:11] + line_FX[i][0:12]
        time     = datetime.datetime.strptime(time_str, '%Y-%m-%d %H:%M:%S.%f')
        orbit1.epoch.append(time)
        #
        if len(line_FX[i]) != 76:
            print(i, time_str, ': Length error (default:76) - ', len(line_FX[i]))
            orbit_nan(orbit1, i)
        else:
            # SM & tilt angle
            orbit1.orbit_SM[i][0] = line_FX[i][39:47];  orbit1.orbit_SM[i][1] = line_FX[i][47:55];  orbit1.orbit_SM[i][2] = line_FX[i][55:63]
            orbit1.tilt[i] = line_FX[i][63:71]

            # MLAT
            if orbit1.orbit_SM[i][2] == 0.0:
                orbit1.MLAT[i] = 0.0
            else:
                orbit1.MLAT[i] = math.atan(orbit1.orbit_SM[i][2]/(orbit1.orbit_SM[i][0]**2 + orbit1.orbit_SM[i][1]**2)**0.5)*180./math.pi 

            # MLT
            if orbit1.orbit_SM[i][0] == 0.0:
                if orbit1.orbit_SM[i][1] < 0.0: orbit1.MLT[i] =  6.0
                else:                           orbit1.MLT[i] = 18.0
            else:
                orbit1.MLT[i] = math.atan(orbit1.orbit_SM[i][1]/orbit1.orbit_SM[i][0])*180./math.pi*24./360.
                if   orbit1.orbit_SM[i][0] > 0:   orbit1.MLT[i] = orbit1.MLT[i] + 12.
                elif orbit1.orbit_SM[i][1] > 0:   orbit1.MLT[i] = orbit1.MLT[i] + 24.

            # Distance
            orbit1.R[i] = (orbit1.orbit_SM[i][0]**2 + orbit1.orbit_SM[i][1]**2 + orbit1.orbit_SM[i][2]**2)**0.5
    return orbit1


def orbit_add(orbit, orbit1):
    orbit.orbit_SM = np.r_["0", orbit.orbit_SM, orbit1.orbit_SM]
    orbit.R        = np.r_["0", orbit.R,        orbit1.R]
    orbit.MLAT     = np.r_["0", orbit.MLAT,     orbit1.MLAT]
    orbit.MLT      = np.r_["0", orbit.MLT,      orbit1.MLT]
    orbit.tilt     = np.r_["0", orbit.tilt,     orbit1.tilt]
    orbit.epoch    = np.r_["0", orbit.epoch,    orbit1.epoch]
    orbit.num      = orbit.num + orbit1.num
    return orbit


def orbit_nan(orbit, i):
    orbit.orbit_SM[i] = np.nan
    orbit.R[i]        = np.nan
    orbit.MLAT[i]     = np.nan
    orbit.MLT[i]      = np.nan
    orbit.tilt[i]     = np.nan
    return


def name_FX_data(YY, MM, DD, data_dir):
    name_dir  = data_dir + '/KPRM/FX/'+YY+'/'
    name_file = YY + MM + DD + '_gtl_fx_SM_DSN.txt'
    name_FX_file = name_dir + name_file
    return  name_FX_file