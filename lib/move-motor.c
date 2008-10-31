// -------------------------------------------------------------------
// move-motor.c
//
// Purpose: The move-motor library provides a simple mechanism for
// outscanning kinematics to motors in hard realtime from userspace 
// (RTAI lxrt).  The motors are controlled using clock and direction 
// commands via dio board using the comedi API. 
// 
// Author: Will Dickson 09/30/2008
// -------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/mman.h>
#include <rtai_lxrt.h>
#include <rtai_comedi.h>
#include <signal.h>
#include "move-motor.h"

// RT-task parameters
const int PRIORITY=1;
const int STACK_SIZE=4096;
const int MSG_SIZE=0;

// Function prototypes
void sigint_func(int sig);
void print_err_msg(char *file, int line, char *err_msg);
int get_kine_pos(kine_t *kine, int motor_n, int index, int *pos);
int check_motors(motors_t motors);
int check_kine(kine_t kine);
int check_dt(int dt_ns);
int get_max_motor(void);
int get_max_dt(void);
int get_min_dt(void);
int get_clock_hi_ns(void);

// Global variables
int sigint_flag = 0;

// ------------------------------------------------------------------
// function: make_kine
//
// Purpose: Outscan kine to motors given a array of motor positions
// vs time (the kine) and the motor configuration.
// 
// Inputs:
//   kine  = structure containing motor kine array.
//   motors = motor configuration structure.
//   dt_ns  = realtime loop period
//
// Output: return code FAIL or SUCCESS.
//    
// ------------------------------------------------------------------
int outscan_kine(kine_t kine, motors_t motors, int dt_ns, int *end_pos)
{
  int i,j;
  int pos[MAX_MOTOR];   // Array of motor positions
  int pos_last;         // Last motor position
  int pos_diff;         // Difference between current and last position
  RTIME now_ns;         // The current time in ns
  RT_TASK *TASK;        // Real-time task
  void *device;         // Comedi device
  char err_msg[ERR_SZ]; // Error messages
  int return_flag = SUCCESS;
  void *sighandler;

  // Check real-time loop period
  if (check_dt(dt_ns) == FAIL) {
    print_err_msg(__FILE__, __LINE__, "invalid real-time loop period");
    return FAIL;
  }

  // Check motor config
  if (check_motors(motors) == FAIL) {
    print_err_msg(__FILE__, __LINE__, "motor configuration invalid");
    return FAIL;
  }

  // Check kinematics
  if (check_kine(kine) == FAIL) {
    print_err_msg(__FILE__, __LINE__, "kinematics invalid");
    return FAIL;
  }
  
  // Check that kine and motor config are compatible
  if (kine.ncol != motors.num) {
    snprintf(err_msg, ERR_SZ, "kine and motor configuration incompatible");
    print_err_msg(__FILE__, __LINE__, err_msg);
    return FAIL;
  }

  // Initialize pos
  for (j=0; j<motors.num; j++){
    if (get_kine_pos(&kine,j,0,&pos[j]) == FAIL) {
      print_err_msg(__FILE__, __LINE__, "initializing pos");
      return FAIL;
    }
  }

  // Setup SIGINT handler
  sighandler = signal(SIGINT,sigint_func);
  if (sighandler == SIG_ERR) {
    print_err_msg(__FILE__,__LINE__,"assigning SIGINT handler failed");
    return FAIL;
  }

  // Initialize realtime task
  rt_allow_nonroot_hrt();
  TASK = rt_task_init(nam2num(TASK_NAME), PRIORITY, STACK_SIZE, MSG_SIZE);
  if (TASK==0) {
    sighandler = signal(SIGINT,sighandler);
    print_err_msg(__FILE__, __LINE__, "rt_init_task failed\n");
    return FAIL;
  }
  rt_task_use_fpu(TASK,1);
  
  // Open comedi device 
  device = comedi_open(motors.dev_name);
  if (device==0) {
    rt_task_delete(TASK);
    sighandler = signal(SIGINT,sighandler);
    print_err_msg(__FILE__, __LINE__, "comedi_open failed\n");
    return FAIL;
  }

  // Set clk/dir pins to output 
  for (i=0; i<motors.num; i++) {
    comedi_dio_config(device, motors.subdev, motors.dio_clk[i], COMEDI_OUTPUT);
    comedi_dio_config(device, motors.subdev, motors.dio_dir[i], COMEDI_OUTPUT);
  }

  // Set all pins (ckl+dir) to low
  for (i=0; i<motors.num; i++) {
    comedi_dio_write(device, motors.subdev, motors.dio_clk[i], DIO_LO);
    comedi_dio_write(device, motors.subdev, motors.dio_dir[i], DIO_LO);
  }

  printf("Starting real-time outscan \n");
  fflush(stdout);

  // Setup periodic timer
  rt_set_oneshot_mode();
  start_rt_timer(0);

  // Go to real-time
  mlockall(MCL_CURRENT|MCL_FUTURE);
  rt_make_hard_real_time();

  // Outscan loop
  for (i=0; i<kine.nrow; i++) {

    now_ns = rt_get_time_ns();

    for (j=0; j<motors.num; j++) {
    
      // Get motor positions
      pos_last = pos[j];
      if (get_kine_pos(&kine,j,i,&pos[j]) == FAIL) {
	print_err_msg(__FILE__, __LINE__, "exiting realtime loop");
	return_flag = FAIL;
	goto RT_LOOP_EXIT; // Exit from loops and goto clean up code.
      }
      pos_diff = pos[j] - pos_last;

      // Set direction dio
      if (pos_diff > 0){
	comedi_dio_write(device, motors.subdev, motors.dio_dir[j], DIO_HI);
      }
      else {
	comedi_dio_write(device, motors.subdev, motors.dio_dir[j], DIO_LO);
      }

      // Set clock dio
      if (abs(pos_diff) > 0) {
	comedi_dio_write(device, motors.subdev, motors.dio_clk[j], DIO_HI);
      }

    } // end for j

    // Sleep for CLOCK_HI_NS and then set clocks lines low
    rt_sleep_until(nano2count(now_ns + CLOCK_HI_NS));
    for (j=0; j<motors.num; j++) {
      comedi_dio_write(device, motors.subdev, motors.dio_clk[j], DIO_LO);
    }

    // Check sigint flag
    if (sigint_flag == 1) {
      fprintf(stdout, "SIGINT - exiting real-time\n");
      return_flag = SIGINT;
      goto RT_LOOP_EXIT;
    }

    // Sleep until next period
    rt_sleep_until(nano2count(now_ns+dt_ns));

  } // end for i
 
RT_LOOP_EXIT:
  // Leave realtime
  rt_make_soft_real_time();
  munlockall();
  stop_rt_timer();

  printf("outscan done\n");

  // Set all pins (ckl+dir) to low
  for (i=0; i<motors.num; i++) {
    comedi_dio_write(device, motors.subdev, motors.dio_clk[i], DIO_LO);
    comedi_dio_write(device, motors.subdev, motors.dio_dir[i], DIO_LO);
  }
  
  // Clean up
  comedi_close(device);
  rt_task_delete(TASK);

  // Restore old SIGINT handler
  sighandler = signal(SIGINT,sighandler);
  if (sighandler == SIG_ERR) {
    print_err_msg(__FILE__,__LINE__,"restoring signal handler failed");
    return_flag = FAIL;
  }

  // Set final positon
  for (j=0; j<motors.num; j++) {
    end_pos[j] = pos[j];
  }
 
  return return_flag;
}

// ---------------------------------------------------------------
// sigint_func
//
// ---------------------------------------------------------------
void sigint_func(int sig) {
  sigint_flag = 1;
  return;
}


// ----------------------------------------------------------------
// function: get_kine_pos
//
// Purpose: Returns the kine position for motor # (motor_n) at 
// outscan index (index).
//
// Inputs:
//   kine   = pointer to kine structure
//   motor_n = motor nunber
//   index   = index at which to get position
//
// Output: kine position 
//   
// ----------------------------------------------------------------
int get_kine_pos(kine_t *kine, int motor_n, int index, int *pos)
{
  int s0, s1;
  int *ptr;

  // Check motor_n and index ranges
  if ((motor_n < 0) || (motor_n >= kine->ncol)) {
    print_err_msg(__FILE__,__LINE__,"(Fatal) motor_n out of range");
    return FAIL;
  }
  if ((index <0) || (index > kine->nrow)) {
    print_err_msg(__FILE__,__LINE__,"(Fatal) index out of range");
    return FAIL;
  }
  // Get kine position
  s0 = kine -> s0;
  s1 = kine -> s1;
  ptr = (int*)((kine -> data) + index*s0 + motor_n*s1);
  *pos = *ptr;
  return SUCCESS;
}

// ----------------------------------------------------------------
// function: check_kine
//
// Purpose: Check that kinematics are valid
//
// ----------------------------------------------------------------
int check_kine(kine_t kine)
{
  int i,j;
  int flag = SUCCESS;
  int p0,p1;
  int err_flag;

  for (i=0;i<(kine.nrow-1); i++) {
    for (j=0; j<kine.ncol; j++) {
      err_flag = get_kine_pos(&kine,j,i, &p0);
      err_flag = get_kine_pos(&kine,j,i+1, &p1);
      if (abs(p1-p0) > 1) {
	flag = FAIL;
      }
    }
  }
  return flag;
}

// ----------------------------------------------------------------
// function: check_dt
//
// Purpose: Check that the realtime loop period is valid
//
// ----------------------------------------------------------------
int check_dt(int dt_ns)
{
  int flag = SUCCESS;
  if (dt_ns > MAX_DT_NS) {
    flag = FAIL;
  }
  else if (dt_ns < MIN_DT_NS) {
    flag = FAIL;
  }
  else if (dt_ns < CLOCK_HI_NS) {
    flag = FAIL;
  }
  return flag;
}

// ----------------------------------------------------------------
// function: check_motors
//
// Purpose: Checks that motor configuration is valid
//
// ----------------------------------------------------------------
int check_motors(motors_t motors)
{
  int i,j;
  int flag = SUCCESS;

  // Check number of motors
  if ((motors.num <= 0) || (motors.num>MAX_MOTOR)) {
    return FAIL;
  }
 
  for (i=0; i<motors.num; i++) {   
    // Check clk and dir range
    if ((motors.dio_clk[i] < 0) || (motors.dio_clk[i] > MAX_DIO)) {
      flag = FAIL;
    }
    if ((motors.dio_dir[i] < 0) || (motors.dio_dir[i] > MAX_DIO)) {
      flag = FAIL;
    }
    // Check uniqness
    if (motors.dio_clk[i] == motors.dio_dir[i]) {
      flag = FAIL;
    }
    if (i<motors.num) {
      for (j=i+1; j < motors.num; j++) {
	if (motors.dio_clk[i] == motors.dio_clk[j]) {
	  flag = FAIL;
	}
	if (motors.dio_clk[i] == motors.dio_dir[j]) {
	  flag = FAIL;
	}
	if (motors.dio_dir[i] == motors.dio_clk[j]) {
	  flag = FAIL;
	}
	if (motors.dio_dir[i] == motors.dio_dir[j]) {
	  flag = FAIL;
	}
      } // End for j
    }
  } // End for i
  return flag;
}

// -----------------------------------------------------------------
// function: print_err_msg
//
// 
// Print simple error message
// -----------------------------------------------------------------
void print_err_msg(char *file, int line, char *err_msg)
{
  fprintf(stderr, "%s:%d Error, %s\n",file, line, err_msg);
  return;
}

// Simple functions for getting constants to python ctypes interface 
int get_max_motor(void) {return MAX_MOTOR;};
int get_max_dt(void) {return MAX_DT_NS;};
int get_min_dt(void) {return MIN_DT_NS;};
int get_clock_hi_ns(void) {return CLOCK_HI_NS;};




