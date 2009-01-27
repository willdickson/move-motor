
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
import scipy


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
