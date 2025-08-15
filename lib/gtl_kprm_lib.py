"""
    Geotail KPRM MOMFXP lib -- 2025/8/15
"""
import datetime
import math
import numpy as np

# ---------------------------------------------------------
# Geotail MOMFXP
# ---------------------------------------------------------
def momfxp_read_multi(Epoch0_min, Epoch0_max, data_dir):
    """
        momfxp.X/Y/Z[epoch]     X/Y/Z in SM coordinate (Re)  
        momfxp.R[epoch]         R (Re)
        momfxp.MLAT[epoch]      MLAT (deg)
        momfxp.MLT[epoch]       MLT (h)
        momfxp.tilt[epoch]      Tilt angle of Magnetic pole (deg)
        momfxp.epoch[epoch]     Epoch
    """
    class struct:
        pass
    momfxp     = struct()
    momfxp.num = 0

    Epoch_min = datetime.datetime.strptime(Epoch0_min[0:10], "%Y-%m-%d")    # Start time
    Epoch_max = datetime.datetime.strptime(Epoch0_max[0:10], "%Y-%m-%d")    # End   time
    Epoch_len = Epoch_max - Epoch_min;  print('[Total]', Epoch_len.days+1, 'days      (', Epoch0_min[0:10], '-', Epoch0_max[0:10],  ')')
    Epoch     = Epoch_min
    for i in range (Epoch_len.days+1):
        str_Epoch = Epoch.strftime('%Y-%m-%d ')
        momfxp1 = momfxp_read(str_Epoch, data_dir)
        if momfxp1.num > 0:
            if momfxp.num == 0:
                momfxp = momfxp1
                momfxp.epoch = np.array(momfxp1.epoch)
            else:
                momfxp = momfxp_add(momfxp, momfxp1)
        Epoch = Epoch + datetime.timedelta(days=1)
        print(i, str_Epoch, momfxp.num)

    # Fill NAN into GAP
    """
    if mode_gap == 1:
        for i in range(momfxp.num-1):
            dt = momfxp.epoch[i+1] - momfxp.epoch[i]
            if dt.total_seconds() > 60:
                momfxp_nan(momfxp, i);   momfxp_nan(momfxp, i+1)
    """
    return  momfxp


def momfxp_read(Epoch, data_dir):
    # initialize
    class struct:
        pass
    momfxp1       = struct()
    momfxp1.num   = 0

    # file open
    Epoch_YY    = Epoch[2:4];  Epoch_MM = Epoch[5:7];  Epoch_DD = Epoch[8:10]
    name_MOMFXP_file = name_MOMFXP_data(Epoch_YY, Epoch_MM, Epoch_DD, data_dir)
    try:
        with open(name_MOMFXP_file, 'r') as f:
            f.close()
    except FileNotFoundError:
        print("***ERROR*** MOMFXP file - not found:", name_MOMFXP_file)
        return momfxp1

    # Decode momfxp
    f = open(name_MOMFXP_file, 'r')
    line_MOMFXP = f.readlines()
    f.close()
    line_MOMFXP = line_MOMFXP[1:len(line_MOMFXP)]
    num = len(line_MOMFXP)
    if num == 0:
        return momfxp1

    # variable set
    momfxp1.num = 60*24
    momfxp1.epoch = []
    momfxp1.N     = np.zeros(momfxp1.num)
    momfxp1.T     = np.zeros(momfxp1.num)
    momfxp1.I     = np.zeros(momfxp1.num)
    momfxp1.Vsc   = np.zeros(momfxp1.num)

    # loop: 00:00 - 23:59
    hh = np.int32(line_MOMFXP[0][0:2]);  mm = np.int32(line_MOMFXP[0][3:5]);  i_momfxp = hh*60 + mm
    j = 0;  k = 0;  N = 0.0;  T = 0.0;  Vsc = 0.0
    for i in range(60*24):
        hh = i//60;  mm = i%60;  time_str = Epoch[0:11] + f"{hh:02d}:{mm:02d}:00"
        time = datetime.datetime.strptime(time_str, '%Y-%m-%d %H:%M:%S')
        momfxp1.epoch.append(time)
        if i < i_momfxp:
            momfxp_nan(momfxp1, i)
        else:
            hh = np.int32(line_MOMFXP[j][0:2]);      mm = np.int32(line_MOMFXP[j][3:5]);  i_momfxp = hh*60 + mm
            while j < num:
                hh = np.int32(line_MOMFXP[j][0:2]);  mm = np.int32(line_MOMFXP[j][3:5]);  i_momfxp = hh*60 + mm
                if i < i_momfxp:
                    break
                # if ( (len(line_MOMFXP[j]) == 109 and line_MOMFXP[j][105:107] == 'Of') or (len(line_MOMFXP[j]) == 108 and line_MOMFXP[j][105:107] == 'On')):
                # if (len(line_MOMFXP[j]) == 108 and line_MOMFXP[j][105:107] == 'On'):
                if (len(line_MOMFXP[j]) == 109 and line_MOMFXP[j][105:107] == 'Of'):
                    if np.float32(line_MOMFXP[j][92:98]) == 0.0:
                        print(j, time_str, ': Vsc disable = 0.0')
                    N   +=  np.float32(line_MOMFXP[j][12:19])
                    T   += (np.float32(line_MOMFXP[j][37:43]) + np.float32(line_MOMFXP[j][43:49]))/2.0
                    Vsc +=  np.float32(line_MOMFXP[j][92:98])
                    k   += 1
                else:
                    print(j, time_str, ': Length error (default:109 or 108) - ', len(line_MOMFXP[j]))
                j += 1
            if k==0:
                momfxp_nan(momfxp1, i)
            else:
                momfxp1.N[i]   = N / k
                momfxp1.T[i]   = T / k
                momfxp1.I[i]   = lep_I(momfxp1.N[i], momfxp1.T[i])
                momfxp1.Vsc[i] = Vsc / k
                # print(i, "\t", momfxp1.epoch[i], f"{momfxp1.N[i]:5.2f}\t{momfxp1.T[i]:5.2f}\t{momfxp1.Vsc[i]:5.2f}\t{j:4d}\t{k:4d}\t{num:4d}")
            k = 0;  N = 0.0;  T = 0.0;  Vsc = 0.0
            if j >= num:
                break
    for j in range(i+1, 60*24):
        hh = j//60;  mm = j%60;  time_str = Epoch[0:11] + f"{hh:02d}:{mm:02d}:00"
        time = datetime.datetime.strptime(time_str, '%Y-%m-%d %H:%M:%S')
        momfxp1.epoch.append(time)
        momfxp_nan(momfxp1, j)

    momfxp1.epoch = np.array(momfxp1.epoch)
    return momfxp1


def momfxp_add(momfxp, momfxp1):
    momfxp.N     = np.r_["0", momfxp.N,     momfxp1.N]
    momfxp.T     = np.r_["0", momfxp.T,     momfxp1.T]
    momfxp.I     = np.r_["0", momfxp.I,     momfxp1.I]
    momfxp.Vsc   = np.r_["0", momfxp.Vsc,   momfxp1.Vsc]
    momfxp.epoch = np.r_["0", momfxp.epoch, momfxp1.epoch]
    momfxp.num   = momfxp.num + momfxp1.num
    return momfxp


def momfxp_nan(momfxp, i):
    momfxp.N[i]   = np.nan
    momfxp.T[i]   = np.nan
    momfxp.I[i]   = np.nan
    momfxp.Vsc[i] = np.nan
    return


def name_MOMFXP_data(YY, MM, DD, data_dir):
    name_dir  = data_dir + '/KPRM/MOMFXP/'+YY+'/'
    name_file = YY + MM + DD + '_gtl_momfxp_SM_DSN.txt'
    name_MOMFXP_file = name_dir + name_file
    return  name_MOMFXP_file




# ---------------------------------------------------------
# FX
# ---------------------------------------------------------
def fx_read_multi(Epoch0_min, Epoch0_max, data_dir):
    """
        fx.X/Y/Z[epoch]     X/Y/Z in SM coordinate (Re)  
        fx.R[epoch]         R (Re)
        fx.MLAT[epoch]      MLAT (deg)
        fx.MLT[epoch]       MLT (h)
        fx.tilt[epoch]      Tilt angle of Magnetic pole (deg)
        fx.epoch[epoch]     Epoch
    """
    class struct:
        pass
    fx     = struct()
    fx.num = 0

    Epoch_min = datetime.datetime.strptime(Epoch0_min[0:10], "%Y-%m-%d")    # Start time
    Epoch_max = datetime.datetime.strptime(Epoch0_max[0:10], "%Y-%m-%d")    # End   time
    Epoch_len = Epoch_max - Epoch_min;  print('[Total]', Epoch_len.days+1, 'days      (', Epoch0_min[0:10], '-', Epoch0_max[0:10],  ')')
    Epoch     = Epoch_min
    for i in range (Epoch_len.days+1):
        str_Epoch = Epoch.strftime('%Y-%m-%d ')
        fx1 = fx_read(str_Epoch, data_dir)
        if fx1.num > 0:
            if fx.num == 0:
                fx = fx1
                fx.epoch = np.array(fx1.epoch)
            else:
                fx = fx_add(fx, fx1)
        Epoch = Epoch + datetime.timedelta(days=1)
        print(i, str_Epoch, fx.num)

    # Fill NAN into GAP
    """
    if mode_gap == 1:
        for i in range(fx.num-1):
            dt = fx.epoch[i+1] - fx.epoch[i]
            if dt.total_seconds() > 60:
                fx_nan(fx, i);   fx_nan(fx, i+1)
    """
    return  fx


def fx_read(Epoch, data_dir):
    # initialize
    class struct:
        pass
    fx1       = struct()
    fx1.num   = 0

    # file open
    Epoch_YY    = Epoch[2:4];  Epoch_MM = Epoch[5:7];  Epoch_DD = Epoch[8:10]
    name_fx_file = name_fx_data(Epoch_YY, Epoch_MM, Epoch_DD, data_dir)
    try:
        with open(name_fx_file, 'r') as f:
            f.close()
    except FileNotFoundError:
        print("***ERROR*** fx file - not found:", name_fx_file)
        return fx1

    # Decode fx
    f = open(name_fx_file, 'r')
    line_fx = f.readlines()
    f.close()
    line_fx = line_fx[1:len(line_fx)]
    num = len(line_fx)
    if num == 0:
        return fx1

    # variable set
    fx1.num = 60*24
    fx1.epoch = []
    fx1.X     = np.zeros(fx1.num)
    fx1.Y     = np.zeros(fx1.num)
    fx1.Z     = np.zeros(fx1.num)
    fx1.R     = np.zeros(fx1.num)
    fx1.MLT   = np.zeros(fx1.num)
    fx1.MLAT  = np.zeros(fx1.num)
    fx1.tilt  = np.zeros(fx1.num)

    # loop: 00:00 - 23:59
    hh = np.int32(line_fx[0][0:2]);  mm = np.int32(line_fx[0][3:5]);  i_fx = hh*60 + mm
    j = 0;  k = 0;  X = 0.0;  Y = 0.0;  Z = 0.0;  tilt = 0.0
    for i in range(60*24):
        hh = i//60;  mm = i%60;  time_str = Epoch[0:11] + f"{hh:02d}:{mm:02d}:00"
        time = datetime.datetime.strptime(time_str, '%Y-%m-%d %H:%M:%S')
        fx1.epoch.append(time)
        if i < i_fx:
            fx_nan(fx1, i)
        else:
            hh = np.int32(line_fx[j][0:2]);  mm = np.int32(line_fx[j][3:5]);  i_fx = hh*60 + mm
            while j < num:
                hh = np.int32(line_fx[j][0:2]);  mm = np.int32(line_fx[j][3:5]);  i_fx = hh*60 + mm
                if i < i_fx:
                    break
                if len(line_fx[j]) != 76:
                    print(j, time_str, ': Length error (default:76) - ', len(line_fx[j]))
                else:
                    X    += np.float32(line_fx[j][39:47])
                    Y    += np.float32(line_fx[j][47:55])
                    Z    += np.float32(line_fx[j][55:63])
                    tilt += np.float32(line_fx[j][63:71])
                    k   += 1
                j += 1
            if k==0:
                fx_nan(fx1, i)
            else:
                fx1.X[i]    =   X / k
                fx1.Y[i]    =   Y / k
                fx1.Z[i]    =   Z / k
                fx1.tilt[i] = tilt / k
                fx1.R[i], fx1.MLAT[i], fx1.MLT[i] = orbit_XYZ_to_R_MLAT_MLT(fx1.X[i], fx1.Y[i], fx1.Z[i])
                # print(i, "\t", fx1.epoch[i], f"{fx1.X[i]:5.2f}\t{fx1.Y[i]:5.2f}\t{fx1.Z[i]:5.2f}\t{fx1.tilt[i]:5.2f}\t{j:4d}\t{k:4d}\t{num:4d}")
            k = 0;  X = 0.0;  Y = 0.0;  Z = 0.0;  tilt = 0.0
            if j >= num:
                break
    for j in range(i+1, 60*24):
        hh = j//60;  mm = j%60;  time_str = Epoch[0:11] + f"{hh:02d}:{mm:02d}:00"
        time = datetime.datetime.strptime(time_str, '%Y-%m-%d %H:%M:%S')
        fx1.epoch.append(time)
        fx_nan(fx1, j)
    fx1.epoch = np.array(fx1.epoch)
    return fx1


def fx_add(fx, fx1):
    fx.X     = np.r_["0", fx.X,     fx1.X]
    fx.Y     = np.r_["0", fx.Y,     fx1.Y]
    fx.Z     = np.r_["0", fx.Z,     fx1.Z]
    fx.tilt  = np.r_["0", fx.tilt,  fx1.tilt]
    fx.R     = np.r_["0", fx.R,     fx1.R]
    fx.MLAT  = np.r_["0", fx.MLAT,  fx1.MLAT]
    fx.MLT   = np.r_["0", fx.MLT,   fx1.MLT]
    fx.epoch = np.r_["0", fx.epoch, fx1.epoch]
    fx.num   = fx.num + fx1.num
    return fx


def fx_nan(fx, i):
    fx.X[i]    = np.nan
    fx.Y[i]    = np.nan
    fx.Z[i]    = np.nan
    fx.tilt[i] = np.nan
    fx.R[i]    = np.nan
    fx.MLAT[i] = np.nan
    fx.MLT[i]  = np.nan
    return


def name_fx_data(YY, MM, DD, data_dir):
    name_dir  = data_dir + '/KPRM/fx/'+YY+'/'
    name_file = YY + MM + DD + '_gtl_fx_SM_DSN.txt'
    name_fx_file = name_dir + name_file
    return  name_fx_file


def orbit_XYZ_to_R_MLAT_MLT(X, Y, Z):
    # R
    R = (X**2 + Y**2 + Z**2)**0.5
    # MLAT
    if X == 0.0 and Y == 0.0:
        if Z > 0.0: MLAT =  90.0
        else:       MLAT = -90.0
    else:
        MLAT = math.atan(Z/(X**2 + Y**2)**0.5)*180./math.pi
    # MLT
    if X == 0.0:
        if Y > 0.0: MLT = 18.0
        else:       MLT =  6.0
    else:
        MLT = math.atan(Y/X)*180./math.pi*24./360.
        if   X > 0:   MLT += 12.
        elif Y > 0:   MLT += 24.
    #
    return R, MLAT, MLT


# LEP I (nA/cm2)
def lep_I(N, T):
    # T = 10eV   --> V_e [m/s]    = ( 2E-19 [C] * 10 [eV] / 9E-31 [kg] )**0.5    = 1.5E6 m/s      [光速の0.5%程度]
    V = T * 1000.
    # N = 100/cc --> I_e (nA/cm2) =   2E-19 [C] * 100 [/cc] * 1E2 * 1.5E6 * 1E9  = 3 nA
    I = 2.0E-19 * N * 1.0E2 * V * 1.0E9
    return I