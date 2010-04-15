"""
-----------------------------------------------------------------------
move-utils
Copyright (C) William Dickson, 2008.
  
wbd@caltech.edu
www.willdickson.com

Released under the LGPL Licence, Version 3

This file is part of move-motor.

move-utils is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
    
move-motor is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with move-utils.  If not, see
<http://www.gnu.org/licenses/>.

------------------------------------------------------------------------
libmove-motor.py

Purpose: Provides a ctypes wrapper for the move-motor library.

Author: William Dickson 
------------------------------------------------------------------------
"""
import os
import os.path
import scipy
import scipy.interpolate
import ConfigParser

BORFRC_DIR = os.path.join(os.environ['HOME'],'.borfrc')
DFLT_CAL_DIR = os.path.join(BORFRC_DIR,'motor_cal')
DFLT_MAP_DIR = os.path.join(BORFRC_DIR,'motor_map')

def get_ramp_moves(p0,p1,vmax,a,dt):
    """
    Generate ramp trajectories from vector p0 to vector p1. 
    
    Arguments:
      p0 = vector of starting points
      p1 = vector of ending points
      vmax = maximum allowed velocity
      a = constant acceleration
      
    Output:
      ramp_array = array of ramp trajectories

    """
    ramp_list = []
    for x0,x1 in zip(p0,p1):
        ramp = get_ramp(x0,x1,vmax,a,dt)
        ramp_list.append(ramp)
        
    # Make all ramps the same length by padding the short ones
    # Reshape at the same time so that all ramps are (maxlen,1)
    max_len = max(r.shape[0] for r in ramp_list)        
    for i,r in enumerate(ramp_list):
        if r.shape[0] < max_len:
            pad = r[-1]*scipy.ones((max_len-r.shape[0],))
            r = scipy.hstack((r,pad))
        ramp_list[i] = r.reshape((max_len,1))

    ramp_array = scipy.hstack(ramp_list)
    return ramp_array
        
def get_ramp(x0,x1,vmax,a,dt, output='ramp only'):
    """
    Generate a ramp trajectory from x0 to x1 with constant
    acceleration, a, to maximum velocity v_max. 

    Note, the main purlpose of this routine is to generate a
    trajectory from x0 to x1. For this reason v_max and a are adjusted
    slightly to work with the given time step.

    Arguments:
     x0 = starting position
     x1 = ending position
     vmax = maximum allowed velocity
     a = constant acceleration
     
     Keywords:
       output = 'ramp only' or 'full'
       when ramp only is selected then only the velocity ramp is returned. 
       If 'full' is selected the adjusted acceleration and maximum velocity 
       are also returned.
       
    Ouput:
      ramp = ramp trajectory form x0 to x1


    """
    # Insure we are dealing with floating point numbers
    x0, x1 = float(x0), float(x1)
    vmax, a = float(vmax), float(a)
    dt = float(dt)
    vmax, a = abs(vmax), abs(a) # Make sure that v_max and a are positive

    # Check to see if there is anything to do
    if x0==x1:
        return scipy.array([x0])

    # Get distance and sign indicating direction
    dist = abs(x1-x0)
    sign = scipy.sign(x1-x0)

    # Determine if we will reach v_max
    t2vmax = vmax/a
    t2halfdist = scipy.sqrt(0.5*dist/a)
    
    if t2vmax > t2halfdist:
        # Trajectory w/o constant velocity segment  
        T = scipy.sqrt(dist/a)
        n = int(scipy.round_((1.0/dt)*T))
         
        # Adjust accel and duration for rounding of n (discrete time steps)
        a = dist/(n*dt)**2
        T = scipy.sqrt(dist/a)
        
        # Generate trajectory
        t = scipy.linspace(0.0,2.0*T,2*n+1)
        def f1(t):
            return 0.5*sign*a*(t**2)
        def f2(t):
            s = t-T
            return f1(T)+ sign*a*T*s - 0.5*sign*a*s**2
        func_list = [f1,f2]
        cond_list = [t<=T, t>T]
        ramp = x0+scipy.piecewise(t,cond_list,func_list)
          
    else:
        # Trajectory w/ constant velocity segment
        # Compute acceleration time and adjust acceleration
        T1 = vmax/a 
        n = int(scipy.round_(T1/dt))
        a = vmax/(n*dt) # Adjusted acceleration 
        T1 = vmax/a # Adjusted acceleration time  

        # Compute and adjust constant velocity time
        T2 = dist/vmax - T1  
        m = int(scipy.round_(T2/dt))
        vmax = dist/(dt*(n+m)) # Adjusted max velocity  
        T2 = dist/vmax - T1 # Adjusted constant velocity time

        # Generate trajectory
        t = scipy.linspace(0.0,2.0*T1+T2,2*n+m+1)
        def f1(t):
            return 0.5*sign*a*(t**2)
        def f2(t):
            s = t-T1
            return f1(T1) + sign*vmax*s
        def f3(t):
            s = t-T1-T2
            return f2(T1+T2)+sign*vmax*s-0.5*sign*a*s**2 
        func_list = [f1,f2,f3]
        cond_list = [t<=T1, scipy.logical_and(t>T1,t<=T1+T2), t>T1+T2]
        ramp = x0+scipy.piecewise(t,cond_list,func_list)

    if output=='ramp only':
        return ramp
    elif output=='full':
        return ramp, vmax, a
    else:
        raise ValueError, 'unknown keyword option output=%s'%(output,)

def convert2int(move):
    """
    Converts continuous move to indices
    """
    move_index = scipy.around(move)
    return move_index.astype(scipy.integer)

def read_motor_maps(filename,mapdir=DFLT_MAP_DIR,caldir=DFLT_CAL_DIR): 
    """
    Reads motor mapping from configuration file. Read calibration files
    for RC motors if there are any.
    """
    base_filename = filename
    config = ConfigParser.ConfigParser()
    if os.path.exists(filename):
        config.read(filename)
    else:
        filename = os.path.join(mapdir,filename)
        if os.path.exists(filename):
            config.read(filename)
        else:
            raise RuntimeError, 'motor map file not found %s'%(base_filename,)
    
    motor_maps = {}
    for motor in config.sections():
        map = {}
        for opt in config.options(motor):
            val = config.get(motor,opt)
            map[opt] = val
        motor_maps[motor] = map
    
    # Convert values based on type of motor
    for motor,map in motor_maps.iteritems():
        map['clk'] = int(map['clk'])
        map['dir'] = int(map['dir'])
        map['number'] = int(map['number'])
        if map['type'].lower() == 'stepper':
            map['deg_per_ind'] = eval(map['deg_per_ind'])
        elif map['type'].lower() == 'rc':
            map['pulse_inc'] = float(map['pulse_inc'])
            map['pulse_dfl'] = float(map['pulse_dfl'])
            map['pulse_max'] = float(map['pulse_max'])
            map['pulse_min'] = float(map['pulse_min'])
            calfile = map['calfile']
            try:
                #map['caldata'] = scipy.loadtxt(calfile)
                map['caldata'] = read_motor_cal(calfile)
            except:
                calfile = os.path.join(caldir,calfile)
                #map['caldata'] = scipy.loadtxt(calfile)
                map['caldata'] = read_motor_cal(calfile)
        else:
            raise ValueError, 'unkown motor type'
    return motor_maps

def read_motor_cal(filename):
    """
    Read RC servo-motor calibration file. 1st entry in file should be the 
    bais angle subsequent entries should be usec deg pairs. The calibration 
    is corrected for bais and returned.
    """
    fd = open(filename,'r')
    cal = []
    for i, line in enumerate(fd.readlines()):
        line = line.split()
        if i == 0:
            bias = float(line[0])
        else:
            us = float(line[0])
            deg = float(line[1])
            cal.append([us,deg])
    cal = scipy.array(cal)
    cal = cal + bias
    fd.close()
    return cal

def _convert_deg2ind(kine_deg,map):
    """
    Converts kinematics from degrees to motor indices for single motor 
    with specified motor map.
    """
    kine_deg = scipy.array(kine_deg)
    if map['type'].lower() == 'stepper':
        kine_ind = kine_deg*(1.0/map['deg_per_ind'])
    elif map['type'].lower() == 'rc':
        caldata = sort_caldata(map['caldata'],1)
        cal_deg = caldata[:,1]
        cal_us = caldata[:,0]
        interp_func = scipy.interpolate.interp1d(cal_deg,cal_us,kind='linear')
        kine_us = interp_func(kine_deg)
        kine_us = kine_us - map['pulse_dfl']
        kine_ind = kine_us*(1.0/map['pulse_inc'])
    else:
        raise ValueError, 'uknown motor type %s'%(map['type'])
    return kine_ind

def _convert_ind2deg(kine_ind,map):
    """
    Converts kinematics from motor indices to degrees for single motor 
    with specified motor map.
    """
    kine_ind = scipy.array(kine_ind)
    if map['type'].lower() == 'stepper':
        kine_deg = kine_ind*map['deg_per_ind']
    elif map['type'].lower() == 'rc':
        caldata = sort_caldata(map['caldata'],0)
        cal_deg = caldata[:,1]
        cal_us = caldata[:,0]
        interp_func = scipy.interpolate.interp1d(cal_us,cal_deg,kind='linear')
        kine_us = kine_ind*map['pulse_inc'] + map['pulse_dfl']
        kine_deg = interp_func(kine_us)
    else:
        raise ValueError, 'unknown motor type %s'%(map['type'])
    return kine_deg

def deg2ind(kine_deg, motor_maps):
    """
    Converts array of kinematics from degrees to motor indices based motor maps.

    Note, kinematics may be of shape (N,K) or (K,) where N is the number of pts in each
    trajectory and K is the number of motors. Each column j of the array is assumed to 
    correspond to the kinematics for motor number j in the dictionaty of motor maps. 
    """

    kine_ind = scipy.zeros(kine_deg.shape)
    for motor, map in motor_maps.iteritems():
        n = map['number']
        if len(kine_deg.shape) == 1:
            kine_ind = _convert_deg2ind(kine_deg, map)
        elif len(kine_deg.shape) == 2:
            kine_ind[:,n] = _convert_deg2ind(kine_deg[:,n], map)
        else:
            raise ValueError, 'kine shape incompatible'
    return kine_ind

def ind2deg(kine_ind, motor_maps):
    """
    Converts array of kinematics from motor indices to degrees based motor maps.

    Note, kinematics may be of shape (N,K) or (K,) where N is the number of pts in each
    trajectory and K is the number of motors. Each column j of the array is assumed to 
    correspond to the kinematics for motor number j in the dictionaty of motor maps. 
    """
    kine_deg = scipy.zeros(kine_ind.shape)
    for motor, map in motor_maps.iteritems():
        n = map['number']
        if len(kine_ind.shape) == 1:
            kine_deg = _convert_ind2deg(kine_ind, map)
        elif len(kine_ind.shape) == 2:
            kine_deg[:,n] = _convert_ind2deg(kine_ind[:,n], map)
        else:
            raise ValueError, 'kine shape incompatible'
    return kine_deg

def sort_caldata(caldata,col):
    """
    Sorts calibration data into acending order for interpolation along
    the specified column.
    """
    def cmp_func(x,y):
        if x[col] < y[col]:
            return -1
        elif x[col] > y[col]:
            return 1
        else:
            return 0

    caldata_list = zip(list(caldata[:,0]), list(caldata[:,1]))
    caldata_list.sort(cmp=cmp_func)
    return  scipy.array(caldata_list)
    
def _zero_degpos_ind(map):
    """
    Returns the absolute index (relative to the default starting position)
    of the zero degree position for the motor specified by the map.
    """
    zero_degpos_ind = _convert_deg2ind(0.0,map)
    return zero_degpos_ind

def get_zero_degpos_ind(motor_maps):
    """
    Returns an array of the zero degree positions in indices of all motors specified 
    by the dictionary of motor maps.
    """
    motor_num_list = get_motor_num_list(motor_maps)
    num_motors = len(motor_num_list)
    zero_deg = scipy.zeros((num_motors,))
    zero_degpos_ind = deg2ind(zero_deg,motor_maps)
    return zero_degpos_ind

def _zero_indpos_deg(map):
    """
    Returns the position in degrees of the zero absolute index position.
    """
    zero_indpos_deg = _convert_ind2deg(0.0,map)
    return zero_indpos_deg

def get_zero_indpos_deg(motor_maps):
    """
    Returns an array of the zero absolute index positions in degrees for all motors
    specified by the dictionary of motor maps.
    """
    motor_num_list = get_motor_num_list(motor_maps)
    num_motors = len(motor_num_list)
    zero_ind = scipy.zeros((num_motors,))
    zero_indpos_deg = ind2deg(zero_ind,motor_maps)
    return zero_indpos_deg
    
def get_clkdir_pins(motor_maps):
    """
    Get tuples of clock and direction io pins in motor number order.
    """
    def cmp_func(x,y):
        if x['number'] < y['number']:
            return -1
        elif x['number'] > y['number']:
            return 1
        else:
            return 0

    map_list = [v for k,v in motor_maps.iteritems()]
    map_list.sort(cmp=cmp_func)
    clk_pins = [map['clk'] for map in map_list]
    dir_pins = [map['dir'] for map in map_list]
    return tuple(clk_pins), tuple(dir_pins)

def get_motor_num_list(motor_maps):
    """
    Returns list of motor numbers
    """
    num_list = [v['number'] for k,v in motor_maps.iteritems()]
    num_list.sort()
    return num_list

def get_num2name_map(maps):
    """
    Returns motor number to motor name dictionary
    """
    num2name = {}
    for k,v in maps.iteritems():
        num2name[maps[k]['number']] = k
    return num2name

def get_motor_cal(motor_maps,table_size=50):
    """
    Returns list of motors calibrations sorted by motor number.
    """
    cal_list = []

    # Sort motor calibrations by motor number
    num_key_pairs = [(v['number'],k) for k,v in motor_maps.iteritems()]
    num_key_pairs.sort()

    # Loop over motors and get motor calibration data
    for n,k in num_key_pairs:
        cal ={}
        if motor_maps[k]['type'] == 'RC':
            cal['type'] = 'table' 
            # Get bounds and produce lookup table
            max_pos =  motor_maps[k]['caldata'][:,1].max()
            min_pos =  motor_maps[k]['caldata'][:,1].min()
            deg_data = scipy.linspace(min_pos,max_pos,table_size)
            ind_data = deg2ind(deg_data, {k:motor_maps[k]})
            cal['deg_data'] = deg_data
            cal['ind_data'] = ind_data
        else:
            cal['type'] = 'mult'
            cal['deg_per_ind'] = motor_maps[k]['deg_per_ind']
        cal_list.append(cal)

    return cal_list

                          

