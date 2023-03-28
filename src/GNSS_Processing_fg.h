#pragma once
// #ifndef GNSS_PROCESS
// #define GNSS_PROCESS
#include <Eigen/Dense>
#include <Eigen/Core>
#include <gnss_comm/gnss_constant.hpp>
#include <gnss_comm/gnss_ros.hpp>
#include "GNSS_Initialization.h"
#include <common_lib.h>
#include <numeric>

#include <opencv2/core/eigen.hpp>

#include <gtsam/inference/Symbol.h>
#include <gtsam/nonlinear/Values.h>
#include <gtsam/nonlinear/Marginals.h>
#include <gtsam/geometry/Rot3.h>
#include <gtsam/geometry/Pose3.h>
#include <gtsam/geometry/Rot2.h>
#include <gtsam/geometry/Pose2.h>
#include <gtsam/slam/PriorFactor.h>
#include <gtsam/slam/BetweenFactor.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/NonlinearEquality.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <gtsam/nonlinear/ISAM2.h>

#include <gnss_factor/gnss_cp_factor.hpp>
#include <gnss_factor/gnss_lio_hard_factor.hpp>
#include <gnss_factor/gnss_lio_gravity_factor.hpp>
#include <gnss_factor/gnss_cp_factor_nolidar.hpp>
#include <gnss_factor/gnss_ddt_smooth_factor.hpp>
#include <gnss_factor/gnss_dt_ddt_factor.hpp>
#include <gnss_factor/gnss_lio_factor.hpp>
#include <gnss_factor/gnss_lio_factor_nolidar.hpp>
#include <gnss_factor/gnss_prior_factor.hpp>
#include <gnss_factor/gnss_psr_dopp_factor.hpp>
#include <gnss_factor/gnss_psr_dopp_factor_nolidar.hpp>

using gtsam::symbol_shorthand::R; // Pose3 ()
using gtsam::symbol_shorthand::P; // Pose3 (x,y,z,r,p,y)
// using gtsam::symbol_shorthand::V; // Vel   (xdot,ydot,zdot)
using gtsam::symbol_shorthand::B; // clock drift (dt_g,dt_r,dt_e,dt_c)
using gtsam::symbol_shorthand::C; // rate of clock drift  (ddt)
using gtsam::symbol_shorthand::E; // ext_p
using gtsam::symbol_shorthand::F; // pos, vel and imu bias (vel ba bg)
using gtsam::symbol_shorthand::A; // pos, vel
// using gtsam::symbol_shorthand::G; // ext_R
// using gtsam::symbol_shorthand::Y; // local enu (yaw)
// using gtsam::symbol_shorthand::A; // anchor point (anc) total = 18 dimensions

using namespace gnss_comm;

#define WINDOW_SIZE (10) // should be 0

class GNSSProcess
{
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  GNSSProcess();
  ~GNSSProcess();
  
  void Reset();
  void inputEphem(EphemBasePtr ephem_ptr);
  void inputIonoParams(double ts, const std::vector<double> &iono_params);
  void inputGNSSTimeDiff(const double t_diff);
  void Ephemfromrinex(const std::string &rinex_filepath);
  double str2double(const std::string &num_str);
  EphemPtr rinex_line2ephem(const std::vector<std::string> &ephem_lines);
  GloEphemPtr rinex_line2glo_ephem(const std::vector<std::string> &ephem_lines, const uint32_t gpst_leap_seconds);
  void rinex2ephems(const std::string &rinex_filepath, std::map<uint32_t, std::vector<EphemBasePtr>> &sat2ephem_);
  void rinex2iono_params(const std::string &rinex_filepath, std::vector<double> &iono_params);
  void processGNSS(const std::vector<ObsPtr> &gnss_meas, state_input &state, Eigen::Vector3d &omg);
  void processGNSS(const std::vector<ObsPtr> &gnss_meas, state_output &state);
  bool GNSSLIAlign();
  void updateGNSSStatistics(Eigen::Vector3d &pos);
  void inputpvt(double ts, double lat, double lon, double alt);
  void inputlla(double ts, double lat, double lon, double alt);
  void processIMUOutput(double dt, const Vector3d &linear_acceleration, const Vector3d &angular_velocity);
  void processIMU(double dt, const Vector3d &linear_acceleration, const Vector3d &angular_velocity);
  Eigen::Vector3d local2enu(state_input &state);
  Eigen::Vector3d local2enu(Eigen::Matrix3d enu_rot, Eigen::Vector3d anc, Eigen::Vector3d &pos);

  std::map<uint32_t, std::vector<EphemBasePtr>> sat2ephem_rnx;
  std::map<uint32_t, std::vector<EphemBasePtr>> sat2ephem;
  std::vector<double> latest_gnss_iono_params;
  std::map<uint32_t, std::map<double, size_t>> sat2time_index;
  std::map<uint32_t, std::map<double, size_t>> sat2time_index_rnx;
  std::map<std::pair<double, int>, std::map<uint32_t, double[6]>> sat2cp; // 
  std::vector<ObsPtr> gnss_meas_buf[WINDOW_SIZE+1]; //
  std::vector<EphemBasePtr> gnss_ephem_buf[WINDOW_SIZE+1]; // 
  std::map<uint32_t, uint32_t> sat_track_status; //
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

  // use GTSAM
  gtsam::NonlinearFactorGraph gtSAMgraph; // store factors //
  gtsam::Values initialEstimate; // store initial values of the node //

  gtsam::Values isamCurrentEstimate; // 
  gtsam::ISAM2 isam;
  bool gnss_online_init; // no use
  bool invalid_lidar = false;

  // size_t index_delete = 0; // 
  // size_t E_num = 0; // 
  std::deque<std::vector<size_t>> factor_id_frame; // 
  // std::deque<std::pair<double, int>> time_frame;
  size_t id_accumulate; // 
  size_t frame_delete = 0; // 
  bool outlier_rej = false;

  gtsam::noiseModel::Diagonal::shared_ptr margrotNoise;
  gtsam::noiseModel::Diagonal::shared_ptr margposNoise;
  // gtsam::noiseModel::Diagonal::shared_ptr priorvelNoise;
  gtsam::noiseModel::Diagonal::shared_ptr margNoise;
  gtsam::noiseModel::Diagonal::shared_ptr margdtNoise;
  gtsam::noiseModel::Diagonal::shared_ptr margddtNoise;

  gtsam::noiseModel::Diagonal::shared_ptr priorrotNoise;
  gtsam::noiseModel::Diagonal::shared_ptr priorposNoise;
  // gtsam::noiseModel::Diagonal::shared_ptr priorvelNoise;
  gtsam::noiseModel::Diagonal::shared_ptr priorNoise;
  gtsam::noiseModel::Diagonal::shared_ptr priordtNoise;
  // gtsam::noiseModel::Diagonal::shared_ptr margExtNoise;
  gtsam::noiseModel::Diagonal::shared_ptr dtNoise;
  gtsam::noiseModel::Diagonal::shared_ptr priorddtNoise;
  gtsam::noiseModel::Diagonal::shared_ptr ddtNoise;
  gtsam::noiseModel::Diagonal::shared_ptr odomNoise;
  // gtsam::noiseModel::Diagonal::shared_ptr odomNoiseIMU;
  gtsam::noiseModel::Base::shared_ptr odomNoiseIMU;
  gtsam::noiseModel::Base::shared_ptr robustpsrdoppNoise;
  gtsam::noiseModel::Base::shared_ptr robustcpNoise;
  // gtsam::noiseModel::Gaussian::shared_ptr testNoise;

  int frame_num = 0; // 
  double last_gnss_time = 0.0; //
  double gnss_sample_period = 0.1;

  double diff_t_gnss_local;
  double gnss_psr_std_threshold;
  double gnss_dopp_std_threshold;
  double gnss_cp_std_threshold;
  double gnss_cp_time_threshold;
  int gnss_track_num_threshold;
  Eigen::Vector3d ecef_pos, first_xyz_ecef_pvt, first_xyz_ecef_lla, first_lla_pvt, first_lla_lla;
  Eigen::Matrix3d Rot_gnss_init = Eigen::Matrix3d::Identity();
  double gnss_elevation_threshold;
  bool gnss_ready;
  int frame_count; //
  int delete_thred;
  int marg_thred;
  int wind_size = WINDOW_SIZE;
  bool quick_init = false;
  bool ephem_from_rinex = false;
  bool nolidar = false;
  bool nolidar_cur;
  // double para_yaw_enu_local[1];
  double para_rcv_dt[(WINDOW_SIZE+1)*4] = {0}; //
  double para_rcv_ddt[WINDOW_SIZE+1] = {0}; //
  Eigen::Vector3d anc_ecef, gravity_init;
  Eigen::Vector3d anc_local = Eigen::Vector3d::Zero();
  Eigen::Matrix3d R_ecef_enu;
  double yaw_enu_local = 0.0;
  
  void initNoises(void); 
  void runISAM2opt(void);
  void GnssPsrDoppMeas(const ObsPtr &obs_, const EphemBasePtr &ephem_);
  void SvPosCals(const ObsPtr &obs_, const EphemBasePtr &ephem_);
  bool Evaluate(state_input &state, Eigen::Vector3d &omg);
  bool Evaluate(state_output &state);
  double prior_noise = 0.01;
  double marg_noise = 0.01;
  state_input state_;
  state_output state_const_;
  state_input state_last;
  state_output state_const_last;
  double ddt_noise = 0.01;
  double dt_noise = 0.01;
  double odo_noise = 0.01;
  double psr_dopp_noise = 0.01;
  double cp_noise = 0.01;
  double relative_sqrt_info = 10;
  double cp_weight = 1.0;
  int change_ext = 1;
  // double odo_weight = 1.0;
  double outlier_thres = 0.1;
  IntegrationBase* pre_integration = new IntegrationBase{Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero()};

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
};

// # endif