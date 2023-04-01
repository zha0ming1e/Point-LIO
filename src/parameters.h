// #ifndef PARAM_H
// #define PARAM_H
#pragma once
#include <ros/ros.h>
#include <Eigen/Eigen>
#include <Eigen/Core>
#include <cstring>
#include "preprocess.h"
#include "GNSS_Processing_fg.h"
#include "IMU_Processing.h"
#include "LI_init/LI_init.h"
#include <sensor_msgs/NavSatFix.h>
#include <livox_ros_driver/CustomMsg.h>
#include <sensor_msgs/PointCloud2.h>
#include <mutex>
#include <omp.h>
#include <math.h>
#include <thread>
#include <fstream>
#include <csignal>
#include <unistd.h>
#include <Python.h>
#include <condition_variable>
#include <sensor_msgs/Imu.h>
#include <pcl/common/transforms.h>
#include <geometry_msgs/Vector3.h>

extern bool is_first_frame;
extern double lidar_end_time, first_lidar_time, time_con;
extern double last_timestamp_lidar, last_timestamp_imu;
extern int pcd_index;

extern state_input state_in;
extern state_output state_out;
extern std::string lid_topic, imu_topic;
extern bool prop_at_freq_of_imu, check_satu, con_frame, cut_frame;
extern bool use_imu_as_input, space_down_sample;
extern bool extrinsic_est_en, publish_odometry_without_downsample;
extern int  init_map_size, con_frame_num;
extern double match_s, satu_acc, satu_gyro, cut_frame_time_interval;
extern float  plane_thr;
extern double filter_size_surf_min, filter_size_map_min, fov_deg;
extern double cube_len; 
extern float  DET_RANGE;
extern bool   imu_en;
extern double imu_time_inte;
extern double laser_point_cov, acc_norm;
extern double acc_cov_input, gyr_cov_input, vel_cov;
extern double gyr_cov_output, acc_cov_output, b_gyr_cov, b_acc_cov;
extern double imu_meas_acc_cov, imu_meas_omg_cov; 
extern int    lidar_type, pcd_save_interval;
extern std::vector<double> gravity_init, gravity;
extern std::vector<double> extrinT, extrinT_gnss, offline_init_vec;
extern std::vector<double> extrinR, extrinR_gnss;
extern bool   runtime_pos_log, pcd_save_en, path_en;
extern bool   scan_pub_en, scan_body_pub_en;
extern shared_ptr<Preprocess> p_pre;
extern shared_ptr<LI_Init> Init_LI;
extern shared_ptr<ImuProcess> p_imu;
extern shared_ptr<GNSSProcess> p_gnss;
extern bool is_first_frame;

extern double time_diff_lidar_to_imu;
extern std::string gnss_ephem_topic, gnss_glo_ephem_topic, gnss_meas_topic, gnss_iono_params_topic;
extern std::string gnss_tp_info_topic, local_trigger_info_topic, rtk_pvt_topic, rtk_lla_topic;
extern std::vector<double> default_gnss_iono_params;
extern double gnss_local_time_diff, gnss_ekf_noise;
extern bool next_pulse_time_valid, update_gnss;
extern bool time_diff_valid, is_first_gnss;
extern double latest_gnss_time, next_pulse_time; 
extern double time_diff_gnss_local;
extern bool gnss_local_online_sync, nolidar; 
extern double li_init_gyr_cov, li_init_acc_cov, lidar_time_inte, first_imu_time;
extern int cut_frame_num, orig_odom_freq;
extern double online_refine_time; //unit: s
extern bool cut_frame_init;
extern bool GNSS_ENABLE;
extern double time_update_last, time_current, time_predict_last_const, t_last;
extern Eigen::Matrix3d Rot_gnss_init;

extern MeasureGroup Measures;

extern ofstream fout_out, fout_imu_pbp, fout_rtk;
void readParameters(ros::NodeHandle &n);
void open_file();
vect3 SO3ToEuler(const SO3 &orient);
void set_gnss_offline_init(bool nolidar_);
void cout_state_to_file();