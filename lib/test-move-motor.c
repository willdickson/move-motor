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
  move-motor.c

  Purpose: Temporary file for testing and dvelopment of move-motor
  library.
 
  Author: Will Dickson 
---------------------------------------------------------------------- */
#include <stdio.h>
#include <stdlib.h>
#include "move-motor.h"

#define NUM_MOTOR 4
#define NROW 20000
#define DT_NS 300000
#define DEVICE "/dev/comedi0"
#define SUBDEV 2

int main(int argc, char* argv[])
{
  int i,j;
  motors_t motors;
  kine_t kine;
  int temp_data[NROW][NUM_MOTOR];
  int end_pos[NUM_MOTOR] = {0,0,0};
  
  // Sample set of kine
  kine.ncol = NUM_MOTOR;
  kine.data = (void*)temp_data;
  for (i=0;i<NROW; i++){
    for (j=0; j<NUM_MOTOR; j++) {
      temp_data[i][j] = i/(j+1);
    }
  }
  kine.ncol = NUM_MOTOR;
  kine.nrow = NROW;
  kine.s0 = NUM_MOTOR*sizeof(int);
  kine.s1 = sizeof(int);

  // Sample motor configuration
  motors.num = NUM_MOTOR;
  for (j=0; j<NUM_MOTOR; j++) {
    motors.dio_clk[j] = 2*j;
    motors.dio_dir[j] = 2*j+1;
  }
  motors.dev_name = DEVICE;
  motors.subdev = SUBDEV;


  outscan_kine(kine, motors, DT_NS, end_pos);
  for (i=0;i<NUM_MOTOR;i++) {
    printf("end_pos[%d] = %d\n", i, end_pos[i]);
  }

  return SUCCESS;
}
