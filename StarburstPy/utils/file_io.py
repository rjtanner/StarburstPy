# -*- coding: utf-8 -*-
"""
    Functions for controlling file input and output.
    
    Author: Ryan Tanner
    email: ryan.tanner@nasa.gov
"""

import StarburstPy as sb
import os
import glob
import numpy as np

__all__ = ['_file_paths','write_input_txt','write_input_original','read_output_data']

class _file_paths(object):
    """
    Object for holding file paths to output directory, current directory etc.
    
    Default for everything is the current directory.
    
    """
    
    def __init__(self):
        
        self.cur_dir = os.path.curdir
        self.output_dir = os.path.curdir
        self.SB99_dir = os.path.curdir
    
    def _get_output_dir(self):
        return self.__output_dir
    
    def _set_output_dir(self, value):
        out_dir = value
        if os.path.isdir(out_dir):
            self.__output_dir = out_dir
        else:
            os.mkdir(out_dir)
            self.__output_dir = out_dir
            
    output_dir = property(_get_output_dir, _set_output_dir, None, None)
    
    def _get_SB99_dir(self):
        return self.__SB99_dir
    
    def _set_SB99_dir(self, value):
        self.__SB99_dir = value

    SB99_dir = property(_get_SB99_dir, _set_SB99_dir, None, None)


def read_output_data(model_name = None, output_dir = None):
    """
    Function for reading in all original Starburst99 output files.
    
    """
    
    if output_dir is None:
        output_dir = sb.indata.output_dir
    
    if not os.path.exists(output_dir):
        sb.sb_mess.error('{0} does not exist.'.format(output_dir), function = 'read_output_data')
        
    sb.sb_mess.message('Looking for output files in {0}.'.format(output_dir), function = 'read_output_data')
    
    file_list = []
    file_endings = ['input','output','yield','power','snr','quanta','color','ewidth','hires',
                        'ifaspec','irfeature','ovi','spectrum','sptyp1','sptyp2','uvline','wrlines']
    
    if model_name is None:
        sb.sb_mess.message('No model name given. Attempting to find one.', function = 'read_output_data')
        model_found = False
        
        for fe in file_endings:
            possible_file = glob.glob(os.path.join(output_dir,'*.'+fe))[0]
            if os.path.isfile(possible_file):
                model_name = os.path.splitext(os.path.basename(possible_file))[0]
                model_found = True
                break
        if model_found:   
            sb.sb_mess.message('Model {0} found.'.format(model_name), function = 'read_output_data')
        else:
            sb.sb_mess.error('Could not find Starburst99 model data in {0}.'.format(output_dir), function = 'read_output_data')
    
    for fe in file_endings:
        possible_file = os.path.join(output_dir,'{0}.'.format(model_name)+fe)
        if os.path.isfile(possible_file):
            file_list.append(possible_file)
            
    if len(file_list) == 0:
        sb.sb_mess.error('{0} contains no usable output files for {1} model.'.format(output_dir,model_name), function = 'read_output_data')
    
    data = {}
    headers = {}
    
    fdict = {'color' : ['color', 'Time', '130-V', '210-V', 'U-B', 'B-V', 'V-R', 'V-I', 
                           'V-J', 'V-H', 'V-K', 'V-L', 'Mag_V', 'Mag_B', 'Mag_Bol'],
                'ewidth' : ['recombination_lines', 'Time', 'Continuum_H_A', 'Luminosity_H_A', 'Equ_width_H_A', 
                            'Continuum_H_B', 'Luminosity_H_B', 'Equ_width_H_B', 
                            'Continuum_P_B', 'Luminosity_P_B', 'Equ_width_P_B', 
                            'Continuum_B_G', 'Luminosity_B_G', 'Equ_width_B_G'],
                'hires' : ['hires_spectra', 'Time', 'Wavelength', 'Log_Luminosity', 'Normalized'],
                'ifaspec' : ['infrared_spectra', 'Time', 'Wavelength', 'Log_Luminosity', 'Normalized'],
                'ovi' : ['uv_spectra1', 'Time', 'Wavelength', 'Log_Luminosity', 'Normalized'],
                'power' : ['power', 'Time', 'Power_All', 'Power_OB', 'Power_RSG','Power_LBV','Power_WR',
                           'Energy',
                           'Momentum_All', 'Momentum_OB', 'Momentum_RSG','Momentum_LBV','Momentum_WR'],
                'quanta' : ['quanta', 'Time', 'HI_rate', 'HI%', 
                            'HeI_rate', 'HeI%',
                            'HeII_rate', 'HeII%', 'Log_Luminosity'],
                'snr' : ['snr', 'Time', 'All_Rate', 'All_Power', 'All_Energy', 
                         'Type_IB_Rate', 'Type_IB_Power', 'Type_IB_Energy', 
                         'Typical_Mass', 'Lowest_Mass', 
                         'Stars_SN_Power', 'Stars_SN_Energy'],
                'spectrum' : ['spectrum', 'Time', 'Wavelength', 'Total', 'Stellar','Nebular'],
                'sptyp1' : ['spectral_type1', 'Time', 'All','O','WR','WN','WC','WR/O','WC/WN'],
                'uvline' : ['uv_spectra2', 'Time', 'Wavelength', 'Log_Luminosity', 'Normalized'],
                'yield' : ['chemical_yields','Time','H','He','C','N','O','Mg','Si','S','Fe',
                           'All_winds','All_SN','All_winds_SN','Total_mass']}

    for file in file_list:
        dkey = os.path.splitext(os.path.basename(file))[1][1:]
        if not (dkey == 'input' or dkey == 'output' or dkey == 'irfeature' or 
                dkey == 'sptyp2' or dkey == 'wrlines'):
            data[fdict[dkey][0]], headers[fdict[dkey][0]] = parse_file(file,fdict[dkey][1:])
    
    return data,headers


def parse_file(file_path,dkeys):
    """
    General function for parsing **almost** all Starburst99 output files.
    """

    with open(file_path) as fdata:
        fc = fdata.readlines()
    
    nrows = len(fc)
    is_header = True
    hline = 0
    header = []
    data = {}
    
    while is_header:
        header.append(fc[hline].replace('\n',''))
        sline = header[hline].split()
        if len(sline) > 0:
            if sline[0] == 'TIME':
                is_header = False
        hline += 1
        
    sline = fc[hline].replace('\n','').split()
    if sline[0] == '[s^-1]':
        hline += 1
    
    nrows = nrows - hline
    
    for i in range(len(dkeys)):
        data[dkeys[i]] = np.zeros([nrows])
        
    for line in range(nrows):
        sline = fc[line+hline].replace('\n','').split()
        res = [float(i) for i in sline]
        for i in range(len(dkeys)):
            data[dkeys[i]][line] = res[i]
            
    if data['Time'][0] == data['Time'][1]:
        jj = len(data['Time'])
        nx = 1
        while data['Time'][nx] == data['Time'][nx-1]:
            nx += 1
        ny = int(jj/nx)
        for i in range(len(dkeys)):
            data[dkeys[i]] = np.array(data[dkeys[i]]).reshape(ny,nx)
        data['Time'] = data['Time'][:,0]
        data['Wavelength'] = data['Wavelength'][0,:]
            
    return data,header


def write_input_txt(inputs, full_name):
    """
    Writes input parameters to simple textfile.
        - full_name includes the absolute or relative path to the file.
    """
    f = open(full_name, 'w')
            
    if inputs.fixed_mass:
        mass = inputs.total_mass
        sfr = 'N/A'
    else:
        mass = 'N/A'
        sfr = inputs.SFR
        
    if inputs.log_scale:
        dt = 'N/A'
        steps = inputs.nsteps
    else:
        dt = inputs.dt
        steps = 'N/A'
        
    if hasattr(inputs, 'IMF'):
        IMF = inputs.IMF
    else:
        IMF = ''
    
    lines = ['{:.<35}'.format('Model Name:')+'{0}\n\n'.format(inputs.model_name),
             '{:.<35}'.format('Fixed Mass:')+'{0}\n'.format(inputs.fixed_mass),
             '{:.<35}'.format('    Total Mass:')+'{0} M_sun\n'.format(mass),
             '{:.<35}'.format('Continuous Star Formation:')+'{0}\n'.format(not inputs.fixed_mass),
             '{:.<35}'.format('    SFR:')+'{0} M_sun/yr\n\n'.format(sfr),
             '{:.<35}'.format('IMF:')+'{0}\n'.format(IMF),
             '{:.<35}'.format('IMF Boundaries:')+'{0} M_sun\n'.format(inputs.IMF_bound),
             '{:.<35}'.format('IMF Exponents:')+'{0}\n'.format(inputs.IMF_exp),
             '{:.<35}'.format('    IMF Mass Limits:')+'{0}, {1} M_sun\n'.format(inputs.mass_lo,inputs.mass_hi),
             '{:.<35}'.format('    SN Mass Limits:')+'{0}, {1} M_sun\n\n'.format(inputs.SN_lo_mass,inputs.SN_hi_mass),
             '{:.<35}'.format('Evolutionary Track:')+'{0}\n'.format(inputs.evolve_track),
             '{:.<35}'.format('Metallicity')+'{0}\n'.format(inputs.metallicity),
             '{:.<35}'.format('Wind Model:')+'{0}\n\n'.format(inputs.wind_model),
             '{:.<35}'.format('Time Limits:')+'{0}-{1} yrs\n'.format(inputs.start_time,inputs.end_time),
             '{:.<35}'.format('Log Scale:')+'{0}\n'.format(inputs.log_scale),
             '{:.<35}'.format('    dt:')+'{0} yrs\n'.format(dt),
             '{:.<35}'.format('    Nsteps:')+'{0}\n\n'.format(steps),
             '{:.<35}'.format('Mass Grid:')+'{0}\n'.format(inputs.mass_grid),
             '{:.<35}'.format('Mass Track:')+'{0}\n'.format(inputs.mass_track),
             '{:.<35}'.format('Output Spectra dt:')+'{0} yrs\n'.format(inputs.spec_dt),
             '{:.<35}'.format('Low-res Spectra Output:')+'{0}\n'.format(inputs.low_res_spec),
             '{:.<35}'.format('High-res Spectra Library:')+'{0}\n'.format(inputs.high_res_spec_Z),
             '{:.<35}'.format('UV spectral library:')+'{0}\n'.format(inputs.UV_line_lib),
             '{:.<35}'.format('RGS microturbulance:')+'{0} km/s\n'.format(inputs.RGS_turb_vel),
             '{:.<35}'.format('RGS Abundance:')+'Solar: {0}\n\n'.format(inputs.RGS_abund),
             '{:*^35}'.format('Output')+'\n\n',
             '{:.<35}'.format('QUANTA:')+'{0}\n'.format(inputs.outputs['QUANTA']),
             '{:.<35}'.format('SNR:')+'{0}\n'.format(inputs.outputs['SNR']),
             '{:.<35}'.format('HRD:')+'{0}\n'.format(inputs.outputs['HRD']),
             '{:.<35}'.format('POWER:')+'{0}\n'.format(inputs.outputs['POWER']),
             '{:.<35}'.format('SP:')+'{0}\n'.format(inputs.outputs['SP']),
             '{:.<35}'.format('YIELDS:')+'{0}\n'.format(inputs.outputs['YIELDS']),
             '{:.<35}'.format('SPECTRUM:')+'{0}\n'.format(inputs.outputs['SPECTRUM']),
             '{:.<35}'.format('LINE:')+'{0}\n'.format(inputs.outputs['LINE']),
             '{:.<35}'.format('COLOR:')+'{0}\n'.format(inputs.outputs['COLOR']),
             '{:.<35}'.format('WIDTH:')+'{0}\n'.format(inputs.outputs['WIDTH']),
             '{:.<35}'.format('FEATURES:')+'{0}\n'.format(inputs.outputs['FEATURES']),
             '{:.<35}'.format('OVI:')+'{0}\n'.format(inputs.outputs['OVI']),
             '{:.<35}'.format('HIRES:')+'{0}\n'.format(inputs.outputs['HIRES']),
             '{:.<35}'.format('WRLINES:')+'{0}\n'.format(inputs.outputs['WRLINES']),
             '{:.<35}'.format('IFASPEC:')+'{0}\n'.format(inputs.outputs['IFASPEC'])]
    
    f.writelines(lines)
    
    f.close()
    
      
def write_input_original(inputs, full_name):
    """
    Creates an original Starburst99 input file.
    """
    
    f = open(full_name, 'w')
    
    NAME = inputs.model_name
    if inputs.fixed_mass:
        ISF = -1
    else:
        ISF = 1
    TOMA = inputs.total_mass/1e6
    SFR = inputs.SFR
    NINTERV = len(inputs.IMF_exp)
    XPONENT = ''
    for i in inputs.IMF_exp:
        XPONENT = XPONENT+',{0}'.format(i)
    XPONENT = XPONENT[1:]
    XMASLIM = ''
    for i in inputs.IMF_bound:
        XMASLIM = XMASLIM+',{0}'.format(i)
    XMASLIM = XMASLIM[1:]
    SNCUT = inputs.SN_lo_mass
    BHCUT = inputs.SN_hi_mass
    
    tracks = ['PadovaOrig', 'PadovaAGB', 'GenevaStd', 'GenevaHigh', 'Genevav00', 'Genevav40']
    PZ = [0.0004, 0.004, 0.008, 0.02, 0.05]
    GSH = [0.001, 0.004, 0.008, 0.02, 0.04]
    Gv = [0.002, 0.014]
    track_dict = {
                  tracks[0] : 
                      {PZ[0] : 31,
                       PZ[1] : 32,
                       PZ[2] : 33,
                       PZ[3] : 34,
                       PZ[4] : 35},
                  tracks[1] : 
                      {PZ[0] : 41,
                       PZ[1] : 42,
                       PZ[2] : 43,
                       PZ[3] : 44,
                       PZ[4] : 45},
                  tracks[2] : 
                      {GSH[0] : 11,
                       GSH[1] : 12,
                       GSH[2] : 13,
                       GSH[3] : 14,
                       GSH[4] : 15},
                  tracks[3] : 
                      {GSH[0] : 21,
                       GSH[1] : 22,
                       GSH[2] : 23,
                       GSH[3] : 24,
                       GSH[4] : 25},
                  tracks[4] : 
                      {Gv[0] : 52,
                       Gv[1] : 54},
                  tracks[5] : 
                      {Gv[0] : 62,
                       Gv[1] : 64},
                 }
    IZ = track_dict[inputs.evolve_track][inputs.metallicity]
    
    models = {'Evolution' : 0, 'Emperical' : 1, 'Theoretical' : 2, 'Elson' : 3}
    IWIND = models[inputs.wind_model]
    
    if inputs.start_time == 0:
        inputs.start_time = 1e4
    TIME1 = inputs.start_time/1e6
    if inputs.log_scale: # if time scale is log
        JTIME = 1
    else:
        JTIME = 0
        
    TBIV = inputs.dt/1e6
    ITBIV = inputs.nsteps
    TMAX = inputs.end_time/1e6
    
    grids = {'Full_Isochrone' : 3, 'Isochrone_Large' : 2, 'Large' : 1, 'Small' : 0}
    JMG = grids[inputs.mass_grid]
    
    if inputs.mass_track[0] == 0:
        LMIN = 0
    else:
        LMIN = '{0},{1}'.format(inputs.mass_track[0],inputs.mass_track[1])
        
    TDEL = inputs.spec_dt/1e6
    
    models = {'Planck' : 1, 'Lejeune' : 2, 'Lejuene-Schmutz' : 3, 
              'Lejeune-Hiller' : 4, 'Pauldrach-Hiller' : 5}
    IATMOS = models[inputs.low_res_spec]
    
    Zs = {0.001 : 1, 0.008 : 2, 0.02 : 3, 0.04 : 4}
    ILIB = Zs[inputs.high_res_spec_Z]
    
    if inputs.UV_line_lib == 'Solar':
        ILINE = 1
    else:
        ILINE = 2
    
    if inputs.RGS_abund:
        moop = 0
    else:
        moop = 1
    IVT = '{0},{1}'.format(inputs.RGS_turb_vel,moop)
    
    IO1 = 1 if inputs.outputs['QUANTA'] else -1
    IO2 = 1 if inputs.outputs['SNR'] else -1
    IO3 = 1 if inputs.outputs['HRD'] else -1
    IO4 = 1 if inputs.outputs['POWER'] else -1
    IO5 = 1 if inputs.outputs['SP'] else -1
    IO6 = 1 if inputs.outputs['YIELDS'] else -1
    IO7 = 1 if inputs.outputs['SPECTRUM'] else -1
    IO8 = 1 if inputs.outputs['LINE'] else -1
    IO9 = 1 if inputs.outputs['COLOR'] else -1
    IO10 = 1 if inputs.outputs['WIDTH'] else -1
    IO11 = 1 if inputs.outputs['FEATURES'] else -1
    IO12 = 1 if inputs.outputs['OVI'] else -1
    IO13 = 1 if inputs.outputs['HIRES'] else -1
    IO14 = 1 if inputs.outputs['WRLINES'] else -1
    IO15 = 1 if inputs.outputs['IFASPEC'] else -1
    
    lines = ['MODEL DESIGNATION:                                           [NAME]\n',
             '{0}\n'.format(NAME),
             'CONTINUOUS STAR FORMATION (>0) OR FIXED MASS (<=0):          [ISF]\n',
             '{0}\n'.format(ISF),
             'TOTAL STELLAR MASS [1E6 M_SOL] IF \'FIXED MASS\' IS CHOSEN:    [TOMA]\n',
             '{0}\n'.format(TOMA),
             'SFR [SOLAR MASSES PER YEAR] IF \'CONT. SF\' IS CHOSEN:         [SFR]\n',
             '{0}\n'.format(SFR),
             'NUMBER OF INTERVALS FOR THE IMF (KROUPA=2):                  [NINTERV]\n',
             '{0}\n'.format(NINTERV),
             'IMF EXPONENTS (KROUPA=1.3,2.3):                              [XPONENT]\n',
             '{0}\n'.format(XPONENT),
             'MASS BOUNDARIES FOR IMF (KROUPA=0.1,0.5,100) [SOLAR MASSES]: [XMASLIM]\n',
             '{0}\n'.format(XMASLIM),
             'SUPERNOVA CUT-OFF MASS [SOLAR MASSES]:                       [SNCUT]\n',
             '{0}\n'.format(SNCUT),
             'BLACK HOLE CUT-OFF MASS [SOLAR MASSES]:                      [BHCUT]\n',
             '{0}\n'.format(BHCUT),
             'METALLICITY + TRACKS:                                        [IZ]\n',
             'GENEVA STD: 11=0.001;  12=0.004; 13=0.008; 14=0.020; 15=0.040\n',
             'GENEVA HIGH:21=0.001;  22=0.004; 23=0.008; 24=0.020; 25=0.040\n',
             'PADOVA STD: 31=0.0004; 32=0.004; 33=0.008; 34=0.020; 35=0.050\n',
             'PADOVA AGB: 41=0.0004; 42=0.004; 43=0.008; 44=0.020; 45=0.050\n',
             'GENEVA v00: 51=0.001;  52=0.002; 53=0.008; 54=0.014; 55=0.040\n',
             'GENEVA v40: 61=0.001;  62=0.002; 63=0.008; 64=0.014; 65=0.040\n',
             '{0}\n'.format(IZ),
             'WIND MODEL (0: MAEDER; 1: EMP.; 2: THEOR.; 3: ELSON):        [IWIND]\n',
             '{0}\n'.format(IWIND),
             'INITIAL TIME [1.E6 YEARS]:                                   [TIME1]\n',
             '{0}\n'.format(TIME1),
             'TIME SCALE: LINEAR (=0) OR LOGARITHMIC (=1)                  [JTIME]\n',
             '{0}\n'.format(JTIME),
             'TIME STEP [1.e6 YEARS] (ONLY USED IF JTIME=0):               [TBIV]\n',
             '{0}\n'.format(TBIV),
             'NUMBER OF STEPS        (ONLY USED IF JTIME=1):               [ITBIV]\n',
             '{0}\n'.format(ITBIV),
             'LAST GRID POINT [1.e6 YEARS]:                                [TMAX]\n',
             '{0}\n'.format(TMAX),
             'SMALL (=0) OR LARGE (=1) MASS GRID;\n',
             'ISOCHRONE ON  LARGE GRID (=2) OR FULL ISOCHRONE (=3):        [JMG]\n',
             '{0}\n'.format(JMG),
             'LMIN, LMAX (ALL=0):                                          [LMIN,LMAX]\n',
             '{0}\n'.format(LMIN),
             'TIME STEP FOR PRINTING OUT THE SYNTHETIC SPECTRA [1.e6YR]:   [TDEL]\n',
             '{0}\n'.format(TDEL),
             'ATMOSPHERE: 1=PLA, 2=LEJ, 3=LEJ+SCH, 4=LEJ+SMI, 5=PAU+SMI:   [IATMOS]\n',
             '{0}\n'.format(IATMOS),
             'METALLICITY OF THE HIGH RESOLUTION MODELS                    [ILIB]\n',
             '(1=0.001, 2= 0.008, 3=0.020, 4=0.040):\n',
             '{0}\n'.format(ILIB),
             'LIBRARY FOR THE UV LINE SPECTRUM: (1=SOLAR, 2=LMC/SMC)       [ILINE]\n',
             '{0}\n'.format(ILINE),
             'RSG FEATURE: MICROTURB. VEL (1-6), SOL/NON-SOL ABUND (0,1)   [IVT,IRSG]\n',
             '{0}\n'.format(IVT),
             'OUTPUT FILES (NO<0, YES>=0)                                  [IO1,...]\n',
             '{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14}\n'.format(IO1,
              IO2,IO3,IO4,IO5,IO6,IO7,IO8,IO9,IO10,IO11,IO12,IO13,IO14,IO15),
             '******************************************************************\n',
             '  OUTPUT FILES:         1    SYNTHESIS.QUANTA\n',
             '                        2    SYNTHESIS.SNR\n',
             '                        3    SYNTHESIS.HRD\n',
             '                        4    SYNTHESIS.POWER\n',
             '                        5    SYNTHESIS.SP\n',
             '                        6    SYNTHESIS.YIELDS\n',
             '                        7    SYNTHESIS.SPECTRUM\n',
             '                        8    SYNTHESIS.LINE\n',
             '                        9    SYNTHESIS.COLOR\n',
             '                       10    SYNTHESIS.WIDTH\n',
             '                       11    SYNTHESIS.FEATURES\n',
             '                       12    SYNTHESIS.OVI\n',
             '                       13    SYNTHESIS.HIRES\n',
             '                       14    SYNTHESIS.WRLINES\n',
             '                       15    SYNTHESIS.IFASPEC\n',
             '\n',
             '  CORRESPONDENCE I VS. MASS:\n',
             '  M   120 115 110 105 100  95  90  85  80  75  70  65\n',
             '  I     1   2   3   4   5   6   7   8   9  10  11  12\n',
             '\n',
             '  M    60  55  50  45  40  35  30  25  20  15  10   5\n',
             '  I    13  14  15  16  17  18  19  20  21  22  23  24\n',
             '\n',
             '\n',
             '  M   120 119 118 117 116 115 114 113 112 111 110 109 108 107 106\n',
             '  I     1   2   3   4   5   6   7   8   9  10  11  12  13  14  15\n',
             '\n',
             '  M   105 104 103 102 101 100  99  98  97  96  95  94  93  92  91\n',
             '  I    16  17  18  19  20  21  22  23  24  25  26  27  28  29  30\n',
             '\n',
             '  M    90  89  88  87  86  85  84  83  82  81  80  79  78  77  76\n',
             '  I    31  32  33  34  35  36  37  38  39  40  41  42  43  44  45\n',
             '\n',
             '  M    75  74  73  72  71  70  69  68  67  66  65  64  63  62  61\n',
             '  I    46  47  48  49  50  51  52  53  54  55  56  57  58  59  60\n',
             '\n',
             '  M    60  59  58  57  56  55  54  53  52  51  50  49  48  47  46\n',
             '  I    61  62  63  64  65  66  67  68  69  70  71  72  73  74  75\n',
             '\n',
             '  M    45  44  43  42  41  40  39  38  37  36  35  34  33  32  31\n',
             '  I    76  77  78  79  80  81  82  83  84  85  86  87  88  89  90\n',
             '\n',
             '  M    30  29  28  27  26  25  24  23  22  21  20  19  18  17  16\n',
             '  I    91  92  93  94  95  96  97  98  99 100 101 102 103 104 105\n',
             '\n',
             '  M    15  14  13  12  11  10   9   8   7   6   5   4   3   2   1\n',
             '  I   106 107 108 109 110 111 112 113 114 115 116 117 118 119 120\n',
             '********************************************************************\n']
    f.writelines(lines)
    
    f.close()
    
    
def read_all_output(model_name):
    """
    Function to call all the other read output functions.
    """
    
    output_dir = sb.indata.output_dir
    
    headers = {}
    data = {}
    
    fdict = {'color' : ['color', 'Time', '130-V', '210-V', 'U-B', 'B-V', 'V-R', 'V-I', 
                           'V-J', 'V-H', 'V-K', 'V-L', 'Mag_V', 'Mag_B', 'Mag_Bol'],
                'ewidth' : ['recombination_lines', 'Time', 'Continuum_H_A', 'Luminosity_H_A', 'Equ_width_H_A', 
                            'Continuum_H_B', 'Luminosity_H_B', 'Equ_width_H_B', 
                            'Continuum_P_B', 'Luminosity_P_B', 'Equ_width_P_B', 
                            'Continuum_B_G', 'Luminosity_B_G', 'Equ_width_B_G'],
                'hires' : ['hires_spectra', 'Time', 'Wavelength', 'Log_Luminosity', 'Normalized'],
                'ifaspec' : ['infrared_spectra', 'Time', 'Wavelength', 'Log_Luminosity', 'Normalized'],
                'ovi' : ['uv_spectra1', 'Time', 'Wavelength', 'Log_Luminosity', 'Normalized'],
                'power' : ['power', 'Time', 'Power_All', 'Power_OB', 'Power_RSG','Power_LBV','Power_WR',
                           'Energy',
                           'Momentum_All', 'Momentum_OB', 'Momentum_RSG','Momentum_LBV','Momentum_WR'],
                'quanta' : ['quanta', 'Time', 'HI_rate', 'HI%', 
                            'HeI_rate', 'HeI%',
                            'HeII_rate', 'HeII%', 'Log_Luminosity'],
                'snr' : ['snr', 'Time', 'All_Rate', 'All_Power', 'All_Energy', 
                         'Type_IB_Rate', 'Type_IB_Power', 'Type_IB_Energy', 
                         'Typical_Mass', 'Lowest_Mass', 
                         'Stars_SN_Power', 'Stars_SN_Energy'],
                'spectrum' : ['spectrum', 'Time', 'Wavelength', 'Total', 'Stellar','Nebular'],
                'sptyp1' : ['spectral_type1', 'Time', 'All','O','WR','WN','WC','WR/O','WC/WN'],
                'uvline' : ['uv_spectra2', 'Time', 'Wavelength', 'Log_Luminosity', 'Normalized'],
                'yield' : ['chemical_yields','Time','H','He','C','N','O','Mg','Si','S','Fe',
                           'All_winds','All_SN','All_winds_SN','Total_mass']}
    
    fkeys = list(fdict.keys())

    for file in fkeys:
        headers[fdict[file][0]], data[fdict[file][0]] = parse_file(model_name,output_dir,fdict[file][1:],file)
    
    return headers,data


def parse_file(model_name,output_dir,dkeys,fsuffix):
    """
    General function for parsing almost all Starburst99 output files.
    """
    
    file_name = model_name+'.'+fsuffix
    full_name = os.path.join(output_dir,file_name)
    with open(full_name) as fdata:
        fc = fdata.readlines()
    
    nrows = len(fc)
    is_header = True
    hline = 0
    header = []
    data = {}
    
    while is_header:
        header.append(fc[hline].replace('\n',''))
        sline = header[hline].split()
        if len(sline) > 0:
            if sline[0] == 'TIME':
                is_header = False
        hline += 1
        
    sline = fc[hline].replace('\n','').split()
    if sline[0] == '[s^-1]':
        hline += 1
    
    nrows = nrows - hline
    
    for i in range(len(dkeys)):
        data[dkeys[i]] = np.zeros([nrows])
        
    for line in range(nrows):
        sline = fc[line+hline].replace('\n','').split()
        res = [float(i) for i in sline]
        for i in range(len(dkeys)):
            data[dkeys[i]][line] = res[i]
            
    if data['Time'][0] == data['Time'][1]:
        jj = len(data['Time'])
        nx = 1
        while data['Time'][nx] == data['Time'][nx-1]:
            nx += 1
        ny = int(jj/nx)
        for i in range(len(dkeys)):
            data[dkeys[i]] = np.array(data[dkeys[i]]).reshape(ny,nx)
        data['Time'] = data['Time'][:,0]
        data['Wavelength'] = data['Wavelength'][0,:]
            
    return header,data
            






















