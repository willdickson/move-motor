// ----------------------------------------------------------
// move-motor.h
//
// Purpose: The move-motor library provides a simple mechanism for
// outscanning kinematics to motors in hard realtime from userspace 
// (RTAI lxrt).  The motors are controlled using clock and direction 
// commands via dio board using the comedi API. 
// 
// Author: Will Dickson 09/30/2008
// ----------------------------------------------------------
#define SUCCESS 0
#define FAIL -1
#define DIO_HI 1
#define DIO_LO 0
#define MAX_MOTOR 12
#define MAX_DIO 24
#define ERR_SZ 200
#define MAX_DT_NS 10000000 // 100 Hz 
#define MIN_DT_NS 40000    // 25 kHz
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
