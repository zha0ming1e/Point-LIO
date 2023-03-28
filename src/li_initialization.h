#pragma once

#include <common_lib.h>
#include "Estimator.h"
#define MAXN                (720000)

extern bool data_accum_finished, data_accum_start, online_calib_finish, refine_print;
extern int frame_num_init;
extern double time_lag_IMU_wtr_lidar, move_start_time, online_calib_starts_time; //, mean_acc_norm = 9.81;

extern double timediff_imu_wrt_lidar;
extern bool timediff_set_flg;
extern V3D gravity_lio;
extern mutex mtx_buffer;
extern condition_variable sig_buffer;
extern int scan_count;
extern int frame_ct, wait_num;
extern std::deque<PointCloudXYZI::Ptr>  lidar_buffer;
extern std::deque<double>               time_buffer;
extern std::deque<sensor_msgs::Imu::Ptr> imu_deque;
extern std::queue<std::vector<ObsPtr>> gnss_meas_buf;
extern std::mutex m_time;
extern bool lidar_pushed, imu_pushed;
extern double imu_first_time;
extern bool lose_lid;
extern sensor_msgs::Imu imu_last, imu_next;
extern PointCloudXYZI::Ptr  ptr_con;
extern double T1[MAXN], s_plot[MAXN], s_plot2[MAXN], s_plot3[MAXN], s_plot11[MAXN];

// extern sensor_msgs::Imu::ConstPtr imu_last_ptr;

void gnss_ephem_callback(const GnssEphemMsgConstPtr &ephem_msg);
void gnss_glo_ephem_callback(const GnssGloEphemMsgConstPtr &glo_ephem_msg);
void gnss_iono_params_callback(const StampedFloat64ArrayConstPtr &iono_msg);
void rtk_pvt_callback(const GnssPVTSolnMsgConstPtr &groundt_pvt);
void rtk_lla_callback(const sensor_msgs::NavSatFixConstPtr &lla_msg);
void gnss_meas_callback(const GnssMeasMsgConstPtr &meas_msg);
void local_trigger_info_callback(const ligo::LocalSensorExternalTriggerConstPtr &trigger_msg);
void gnss_tp_info_callback(const GnssTimePulseInfoMsgConstPtr &tp_msg);
void standard_pcl_cbk(const sensor_msgs::PointCloud2::ConstPtr &msg); 
void livox_pcl_cbk(const livox_ros_driver::CustomMsg::ConstPtr &msg); 
void imu_cbk(const sensor_msgs::Imu::ConstPtr &msg_in); 
void LI_Init_set();
bool sync_packages(MeasureGroup &meas, queue<std::vector<ObsPtr>> &gnss_msg);

// #endif