// test-move-motor.c
//
// Simple test code for move-motor library
//
// ---------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include "move-motor.h"

#define NUM_MOTOR 4
#define NROW 20000
#define DT_NS 1000000
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
