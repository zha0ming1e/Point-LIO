#pragma once
// #ifndef GNSS_PROCESS
// #define GNSS_PROCESS
#include "GNSS_Initialization.h"
#include "GNSS_Assignment.h"

using namespace gnss_comm;

#define WINDOW_SIZE (10) // should be 0

class GNSSProcess
{
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  GNSSProcess();
  ~GNSSProcess();
  
  void Reset();
  void inputIonoParams(double ts, const std::vector<double> &iono_params);
  void inputGNSSTimeDiff(const double t_diff);
  void processGNSS(const std::vector<ObsPtr> &gnss_meas, state_input &state, Eigen::Vector3d &omg);
  void processGNSS(const std::vector<ObsPtr> &gnss_meas, state_output &state);
  bool GNSSLIAlign();
  void updateGNSSStatistics(Eigen::Vector3d &pos);
  void inputpvt(double ts, double lat, double lon, double alt);
  void inputlla(double ts, double lat, double lon, double alt);
  void processIMUOutput(double dt, const Vector3d &linear_acceleration, const Vector3d &angular_velocity);
  void processIMU(double dt, const Vector3d &linear_acceleration, const Vector3d &angular_velocity);
  Eigen::Vector3d local2enu(Eigen::Matrix3d enu_rot, Eigen::Vector3d anc, Eigen::Vector3d &pos);
  void SetInit();
  bool AddFactor(gtsam::Rot3 rel_rot_, gtsam::Point3 rel_pos_, gtsam::Vector3 rel_v_, Eigen::Vector3d state_gravity, double delta_t, double time_current,
                Eigen::Vector3d ba, Eigen::Vector3d bg, Eigen::Vector3d omg, Eigen::Matrix3d rot);
  std::map<std::pair<double, int>, std::map<uint32_t, double[6]>> sat2cp; // 
  std::vector<ObsPtr> gnss_meas_buf[WINDOW_SIZE+1]; //
  std::vector<EphemBasePtr> gnss_ephem_buf[WINDOW_SIZE+1]; // 
  Eigen::Matrix3d rot_window[WINDOW_SIZE+1]; //
  Eigen::Vector3d pos_window[WINDOW_SIZE+1]; //
  Eigen::Vector3d vel_window[WINDOW_SIZE+1]; //
  Eigen::Vector3d Tex_imu_r;
  Eigen::Matrix3d Rex_imu_r;
  std::vector<double> pvt_time;
  std::vector<double> lla_time;
  std::vector<Eigen::Vector3d> pvt_holder;
  std::vector<Eigen::Vector3d> lla_holder;
  std::queue<std::vector<ObsPtr>> gnss_msg;

  bool gnss_online_init; // no use
  bool invalid_lidar = false;
  // double dt[4];
  // double ddt; 
  size_t id_accumulate; // 
  size_t frame_delete = 0; // 

  int frame_num = 0; // 
  double last_gnss_time = 0.0; //
  double gnss_sample_period = 0.1;

  double diff_t_gnss_local;
  double gnss_cp_std_threshold;
  double gnss_cp_time_threshold;
  Eigen::Vector3d ecef_pos, first_xyz_ecef_pvt, first_xyz_ecef_lla, first_lla_pvt, first_lla_lla;
  Eigen::Matrix3d Rot_gnss_init = Eigen::Matrix3d::Identity();
  bool gnss_ready;
  int frame_count; //
  int delete_thred;
  int wind_size = WINDOW_SIZE;
  bool nolidar = false;
  bool nolidar_cur;
  std::vector<Eigen::Vector3d> norm_vec_holder;
  // double para_yaw_enu_local[1];
  double para_rcv_dt[(WINDOW_SIZE+1)*4] = {0}; //
  double para_rcv_ddt[WINDOW_SIZE+1] = {0}; //
  Eigen::Vector3d anc_ecef, gravity_init;
  Eigen::Vector3d anc_local = Eigen::Vector3d::Zero();
  Eigen::Matrix3d R_ecef_enu;
  double yaw_enu_local = 0.0;
  
  void runISAM2opt(void);
  void GnssPsrDoppMeas(const ObsPtr &obs_, const EphemBasePtr &ephem_);
  void SvPosCals(const ObsPtr &obs_, const EphemBasePtr &ephem_);
  bool Evaluate(state_input &state, Eigen::Vector3d &omg);
  bool Evaluate(state_output &state);
  state_input state_;
  state_output state_const_;
  state_input state_last;
  state_output state_const_last;
  double relative_sqrt_info = 10;
  double cp_weight = 1.0;
  // double odo_weight = 1.0;
  IntegrationBase* pre_integration = new IntegrationBase{Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero()};
  GNSSAssignment* p_assign = new GNSSAssignment();
  private:
    const ObsPtr obs;
    const EphemBasePtr ephem;
    const std::vector<double> iono_paras;
    int freq_idx;
    double freq;
    Eigen::Vector3d sv_pos;
    Eigen::Vector3d sv_vel;
    double svdt, svddt, tgd;
    double pr_uura, dp_uura;
    Eigen::Matrix3d rot_pos;
};

// # endif