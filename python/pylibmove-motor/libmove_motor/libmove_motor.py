#!/usr/bin/env python
import ctypes
import scipy
from signal import SIGINT

lib = ctypes.cdll.LoadLibrary("libmove-motor.so.1")

# Constants
MAX_MOTOR = lib.get_max_motor()
MIN_DT_NS = lib.get_min_dt()
MAX_DT_NS = lib.get_max_dt()
CLOCK_HI_NS = lib.get_clock_hi_ns()
S2NS = 1.0e9

# Structures
class motors_t(ctypes.Structure):
    _fields_ = [
        ('num', ctypes.c_int),
        ('dio_clk', MAX_MOTOR*ctypes.c_int),
        ('dio_dir', MAX_MOTOR*ctypes.c_int),
        ('dev_name',ctypes.c_char_p),
        ('subdev', ctypes.c_int),
        ]
class kine_t(ctypes.Structure):
    _fields_ = [
        ('data', ctypes.c_void_p),
        ('nrow', ctypes.c_int),
        ('ncol', ctypes.c_int),
        ('s0', ctypes.c_int),
        ('s1', ctypes.c_int),
        ]

# Functions
lib.outscan_kine.restype = ctypes.c_int
lib.outscan_kine.argstype = [
    kine_t,
    motors_t,
    ctypes.c_int,
    ctypes.c_void_p,
    ]
    
def outscan_kine(kine,motor_config,dt):
    """
    Outscan kinematics given the motor configuration and the time step.

    Inputs:
      kine = Nx(num motors) array of motor steps
      motor_config = motor configuration dictionary
      dt = time step in secs
    """
    
    dt_ns = ctypes.c_int()
    dt_ns = int(S2NS*dt)

    # Create and populate motors structure
    motors_struct = motors_t()
    motors_struct.num = motor_config['num']
    motors_struct.dio_clk = motor_config['dio_clk']
    motors_struct.dio_dir = motor_config['dio_dir']
    motors_struct.dev_name = motor_config['device']
    motors_struct.subdev = motor_config['subdevice']

    # Create kinematics structure
    kine_int = kine.astype(scipy.int_)    
    kine_struct = kine_t()
    kine_struct.data = kine_int.ctypes.data_as(ctypes.c_void_p)
    kine_struct.nrow = kine_int.ctypes.shape[0]
    kine_struct.ncol = kine_int.ctypes.shape[1]
    kine_struct.s0 = kine_int.ctypes.strides[0]
    kine_struct.s1 = kine_int.ctypes.strides[1]
    
    # Return array
    end_pos = (ctypes.c_int*motor_config['num'])() 
    for i in range(motor_config['num']):
        end_pos[i] = 10
    
    #Outscan data
    ret_val = lib.outscan_kine(kine_struct, motors_struct, dt_ns, end_pos)
    end_pos = scipy.array(end_pos)
    
    return end_pos, ret_val



# ------------------------------------------------------------------------
if __name__ == "__main__":

    # Test
    import time    
    dt = 1/1000.0
    
    motor_config = {
        'num' : 2,
        'dio_clk' : (0,2),
        'dio_dir' : (1,3),
        'device' : '/dev/comedi0',
        'subdevice' : 2
        }
    
    x = scipy.arange(5000,0,-1)
    y = scipy.arange(0,5000)/2
    kine = scipy.zeros((x.shape[0],motor_config['num']))
    kine[:,0] = x
    kine[:,1] = y

    end_pos, ret_val = outscan_kine(kine,motor_config,dt)
    if ret_val == SIGINT:
        print 'run was interrupted'
    else:
        print 'normal exit'
    print 'end position', end_pos



