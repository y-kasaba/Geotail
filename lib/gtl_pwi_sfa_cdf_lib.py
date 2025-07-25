"""
    Geotail PWI SFA lib -- 2025/7/25
"""
import datetime
import math
import numpy as np
from spacepy import pycdf

# ---------------------------------------------------------
# Geotail PWI SFA
# ---------------------------------------------------------
def read_sfa_multi(Epoch0_min, Epoch0_max, mode_gap, mode_check, data_dir):
    """
            sfa.freq_e[0:384]             Frequency in E-field (Hz)   (1575 - 795363    Hz)
            sfa.freq_b[0:128]             Frequency in B-field (Hz)   (1575 - 12427.547 Hz)
            sfa.sfae_high[epoch][0:384]   Spectra in E-field   (dB V^2/m^2/Hz)
            sfa.sfab_high[epoch][0:128]   Spectra in B-field   (dB nT^2/Hz)
            sfa.antstat[epoch]            Antenna Status              (10)          00-11: ANT of E(1:U/WANT 0:V/PANT) & B(1:alpha 0:BZ)
            sfa.e_ampstat[epoch]          Receiver Gain               (11111)       00000-11111: gain of BAND-1/2/3/4/5 (0:L 1:H)
            sfa.epoch[epoch]              Epoch
    """
    class struct:
        pass
    sfa       = struct()

    Epoch_min = datetime.datetime.strptime(Epoch0_min[0:10], "%Y-%m-%d")     # Start epoch
    Epoch_max = datetime.datetime.strptime(Epoch0_max[0:10], "%Y-%m-%d")     # End   epoch
    Epoch_len = Epoch_max - Epoch_min;   print('[CDF]', Epoch_len.days+1, 'days      (', Epoch0_min[0:10], '-', Epoch0_max[0:10],  ')')
    Epoch     = Epoch_min
    for i in range (Epoch_len.days+1):
        str_Epoch = Epoch.strftime('%Y-%m-%d ')
        sfa1 = read_sfa(str_Epoch, data_dir)
        if i==0 or mode_check == 1:  
            sfa = sfa1
            sfa.freq_e = sfa1.freq_e
            sfa.freq_b = sfa1.freq_b
            sfa.df_e   = np.zeros((128*3))
            sfa.df_e[0:128]   = ( 12500 - 1562.5)/128.
            sfa.df_e[128:256] = (100000 - 12500 )/128.
            sfa.df_e[256:384] = (800000 - 100000)/128.
        else:
            sfa = sfa_add(sfa, sfa1)
        Epoch  = Epoch + datetime.timedelta(days=1)
        print(i+1, str_Epoch, sfa.num)

    # Fill NAN into GAP
    if mode_gap == 1:
        for i in range(sfa.num-1):
            dt = sfa.epoch[i+1] - sfa.epoch[i]
            if dt.total_seconds() > 60:
                sfa.sfae_high[i] = math.nan;  sfa.sfae_high[i+1] = math.nan

    print("   [freq_e] ", sfa.freq_e[0],       "\t",    sfa.freq_e[-1],        "\t",    sfa.freq_e.shape,    "\t",    sfa.freq_e.dtype)
    print("   [freq_b] ", sfa.freq_b[0],       "\t",    sfa.freq_b[-1],        "\t",    sfa.freq_b.shape,    "\t",    sfa.freq_b.dtype)
    print("[sfae_high] ", sfa.sfae_high[0][0], "\t",    sfa.sfae_high[-1][-1], "\t",    sfa.sfae_high.shape, "\t",    sfa.sfae_high.dtype)
    print("[sfab_high] ", sfa.sfab_high[0][0], "\t",    sfa.sfab_high[-1][-1], "\t",    sfa.sfab_high.shape, "\t",    sfa.sfab_high.dtype)
    print("  [antstat] ", sfa.antstat[0],      "\t",    sfa.antstat[-1],       "\t\t",  sfa.antstat.shape,   "\t",    sfa.antstat.dtype)
    print("[e_ampstat] ", sfa.e_ampstat[0],    "\t",    sfa.e_ampstat[-1],     "\t\t",  sfa.e_ampstat.shape, "\t",    sfa.e_ampstat.dtype)
    print("    [epoch] ", sfa.epoch[0],        "\t",    sfa.epoch[-1],         "\t",    sfa.epoch.shape,     "\t",    sfa.epoch.dtype)

    print("")
    # ANT & AMP selection
    sfa.antstat_E  = sfa.antstat    // 10
    sfa.e_ampstat5 =  sfa.e_ampstat         % 10
    sfa.e_ampstat4 = (sfa.e_ampstat //  10) % 10
    sfa.e_ampstat3 = (sfa.e_ampstat // 100) % 10
    print("[antstat(1:U/WANT 0:V/PANT)] ", sfa.antstat_E[0],   sfa.antstat_E[-1],  "\t", sfa.antstat_E.shape,  "\t",  sfa.antstat_E.dtype)
    print("[ampstat_B5  (1:High 0:Low)] ", sfa.e_ampstat5[0],  sfa.e_ampstat5[-1], "\t", sfa.e_ampstat5.shape, "\t",  sfa.e_ampstat5.dtype)
    print("[ampstat_B4  (1:High 0:Low)] ", sfa.e_ampstat4[0],  sfa.e_ampstat4[-1], "\t", sfa.e_ampstat4.shape, "\t",  sfa.e_ampstat4.dtype)
    print("[ampstat_B3  (1:High 0:Low)] ", sfa.e_ampstat3[0],  sfa.e_ampstat3[-1], "\t", sfa.e_ampstat3.shape, "\t",  sfa.e_ampstat3.dtype)

    return  sfa


def read_sfa(Epoch, data_dir):
    class struct:
        pass
    sfa1       = struct()
    Epoch_YYYY = Epoch[0:4];  Epoch_YY = Epoch[2:4];  Epoch_MM = Epoch[5:7];  Epoch_DD = Epoch[8:10]
    name_sfa_file = name_SFA_data(Epoch_YYYY, Epoch_YY, Epoch_MM, Epoch_DD, data_dir)
    print(name_sfa_file)
    try:
        with open(name_sfa_file, 'r') as f:
            f.close()
    except FileNotFoundError:
        print("***ERROR*** SFA file - not found:", name_sfa_file)
        return sfa1

    # Decode SFA
    with open(name_sfa_file, 'r') as f:
        cdf = pycdf.CDF(name_sfa_file)
        # print(cdf)
        sfa1.freq_e    = cdf['freq_e'][...]
        sfa1.freq_b    = cdf['freq_b'][...]
        sfa1.sfae_high = cdf['sfae_high'][...]
        sfa1.sfab_high = cdf['sfab_high'][...]
        sfa1.antstat   = np.uint16(cdf['antstat'][...])
        sfa1.e_ampstat = np.uint16(cdf['e_ampstat'][...])
        sfa1.epoch     = cdf['epoch'][...]
        sfa1.num       = sfa1.epoch.shape[0]
    return sfa1


def sfa_add(sfa, sfa1):
    sfa.sfae_high = np.r_["0", sfa.sfae_high, sfa1.sfae_high]
    sfa.sfab_high = np.r_["0", sfa.sfab_high, sfa1.sfab_high]
    sfa.antstat   = np.r_["0", sfa.antstat,   sfa1.antstat]
    sfa.e_ampstat = np.r_["0", sfa.e_ampstat, sfa1.e_ampstat]
    sfa.epoch     = np.r_["0", sfa.epoch,     sfa1.epoch]
    sfa.num       = sfa.num + sfa1.num
    return sfa


def name_SFA_data(YYYY, YY, MM, DD, data_dir):
    name_dir  = data_dir + '/PWI/cdf/'+YY+'/'+MM+'/'
    name_file = 'GE_H0_SFA_H_' + YYYY + MM + DD + '_v00.cdf'
    name_SFA_file = name_dir + name_file
    return  name_SFA_file