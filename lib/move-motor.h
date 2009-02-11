/*---------------------------------------------------------------------
  move-motor
  Copyright (C) William Dickson, 2008.
  
  wbd@caltech.edu
  www.willdickson.com

  Released under the LGPL Licence, Version 3
  
  This file is part of move-motor.

  move-motor is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as
  published by the Free Software Foundation, either version 3 of the
  License, or (at your option) any later version.
    
  move-motor is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with move-motor.  If not, see
  <http://www.gnu.org/licenses/>.

----------------------------------------------------------------------
  move-motor.h

  Purpose: The move-motor library provides a simple mechanism for
  outscanning kinematics to motors in hard realtime from userspace 
  (RTAI lxrt).  The motors are controlled using clock and direction 
  commands via dio board using the comedi API. 
 
  Author: Will Dickson 
---------------------------------------------------------------------- */
#define SUCCESS 0
#define FAIL -1
#define DIO_HI 1
#define DIO_LO 0
#define MAX_MOTOR 12
#define MAX_DIO 24
#define ERR_SZ 200
#define MAX_DT_NS 10000000 // 100 Hz 
#define MIN_DT_NS 40000    // 20 kHz
#define CLOCK_HI_NS 20000  // 
#define TASK_NAME "move-motor"

// Structure for motor configuration
typedef struct {
  int num;
  int dio_clk[MAX_MOTOR];
  int dio_dir[MAX_MOTOR];
  char *dev_name;
  int subdev;
} motors_t;

// Structure for kinematics - outscan buffer
typedef struct {
  void *data;
  int nrow;
  int ncol; 
  int s0; 
  int s1;
} kine_t;


// Shared function prototypes
int outscan_kine(kine_t kine, motors_t motors, int dt_ns, int *end_pos);
