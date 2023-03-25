
#ifndef GNSS_PROCESS
#define GNSS_PROCESS
#include <opencv2/core/eigen.hpp>

#include <Eigen/Dense>
#include <gnss_comm/gnss_constant.hpp>
#include <gnss_comm/gnss_ros.hpp>
#include "GNSS_Initialization.hpp"
#include <common_lib_gnss.h>
#include <numeric>

using namespace gnss_comm;

#define WINDOW_SIZE (10)

#define PSR_TO_DOPP_RATIO (5)

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
  void processGNSS(const std::vector<ObsPtr> &gnss_meas, StatesGroupwithGNSS &state);
  void processGNSS(const std::vector<ObsPtr> &gnss_meas, StatesGroupwithGNSS2 &state);
  bool GNSSLIAlign();
  void updateGNSSStatistics(Eigen::Vector3d &pos, double yaw_cur, Eigen::Vector3d &anc_cur);
  void inputpvt(double ts, double lat, double lon, double alt);
  void inputlla(double ts, double lat, double lon, double alt);
  Eigen::Vector3d local2enu(StatesGroupwithGNSS &state);
  Eigen::Vector3d local2enu(double yaw, Eigen::Vector3d &anc, Eigen::Vector3d &pos);

  std::map<uint32_t, std::vector<EphemBasePtr>> sat2ephem_rnx;
  std::map<uint32_t, std::vector<EphemBasePtr>> sat2ephem;
  std::vector<double> latest_gnss_iono_params;
  std::map<uint32_t, std::map<double, size_t>> sat2time_index;
  std::map<uint32_t, std::map<double, size_t>> sat2time_index_rnx;
  std::map<double, std::map<uint32_t, double[3]>> sat2cp;
  std::vector<ObsPtr> gnss_meas_buf[WINDOW_SIZE+1];
  std::vector<EphemBasePtr> gnss_ephem_buf[WINDOW_SIZE+1];
  std::map<uint32_t, uint32_t> sat_track_status;
  Eigen::Vector3d pos_window[WINDOW_SIZE+1];
  Eigen::Vector3d vel_window[WINDOW_SIZE+1];
  Eigen::Vector3d Tex_imu_r;
  Eigen::Matrix3d Rex_imu_r;
  std::vector<double> pvt_time;
  std::vector<double> lla_time;
  std::vector<Eigen::Vector3d> pvt_holder;
  std::vector<Eigen::Vector3d> lla_holder;
  std::queue<std::vector<ObsPtr>> gnss_msg;

  double diff_t_gnss_local;
  double gnss_psr_std_threshold;
  double gnss_dopp_std_threshold;
  double gnss_cp_std_threshold;
  double gnss_cp_time_threshold;
  int gnss_track_num_threshold;
  Eigen::Vector3d ecef_pos, first_xyz_ecef_pvt, first_xyz_ecef_lla, first_lla_pvt, first_lla_lla;
  double gnss_elevation_threshold;
  bool gnss_ready;
  int frame_count;
  bool quick_init;
  bool ephem_from_rinex;
  bool nolidar;
  double para_yaw_enu_local[1];
  double para_rcv_dt[(WINDOW_SIZE+1)*4] = {0};
  double para_rcv_ddt[WINDOW_SIZE+1] = {0};
  Eigen::Vector3d anc_ecef;
  Eigen::Vector3d anc_local = Eigen::Vector3d::Zero();
  Eigen::Matrix3d R_ecef_enu;
  double yaw_enu_local = 0.0;

  std::map<uint32_t, std::map<double, size_t>> empty_map_i;
  std::map<uint32_t, std::vector<EphemBasePtr>> empty_map_e;
  std::vector<ObsPtr> empty_vec_o;
  std::vector<EphemBasePtr> empty_vec_e;
  std::vector<double> empty_vec_d;
  std::map<uint32_t, uint32_t> empty_map_t;

  void GnssPsrDoppMeas(const ObsPtr &obs_, const EphemBasePtr &ephem_);
  void SvPosCals(const ObsPtr &obs_, const EphemBasePtr &ephem_);
  bool Evaluate(StatesGroupwithGNSS &state, Eigen::VectorXd &residuals, Eigen::MatrixXd &jacobians, int &meas_size);
  bool Evaluate(StatesGroupwithGNSS2 &state, Eigen::VectorXd &residuals, Eigen::MatrixXd &jacobians, int &meas_size);
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
    double relative_sqrt_info;
};

GNSSProcess::GNSSProcess()
    : diff_t_gnss_local(0.0), gnss_track_num_threshold(20)
{
  Reset();
}

GNSSProcess::~GNSSProcess() {}

void GNSSProcess::Reset() 
{
  ROS_WARN("Reset GNSSProcess");
  sat2time_index.swap(empty_map_i);
  sat2ephem.swap(empty_map_e);
  latest_gnss_iono_params.swap(empty_vec_d);
  sat_track_status.swap(empty_map_t);
  // gnss_meas_buf.clear();
  // gnss_ephem_buf.clear();
  gnss_ready = false;
  frame_count = 0;
  quick_init = false;
  ephem_from_rinex = false;
  nolidar = false;
}

void GNSSProcess::Ephemfromrinex(const std::string &rinex_filepath)
{
  rinex2ephems(rinex_filepath, sat2ephem_rnx);
  std::map<uint32_t, std::vector<EphemBasePtr>>::iterator it;
  for (it = sat2ephem_rnx.begin(); it != sat2ephem_rnx.end(); it++)
  {
    for (int j = 0; j < it->second.size(); j++)
    {
      sat2time_index_rnx[it->first].emplace(time2sec(it->second[j]->toe), j);
    }
  }
}

void GNSSProcess::inputEphem(EphemBasePtr ephem_ptr) // 
{
    double toe = time2sec(ephem_ptr->toe);
    // if a new ephemeris comes
    if (sat2time_index.count(ephem_ptr->sat) == 0 || sat2time_index.at(ephem_ptr->sat).count(toe) == 0)
    {
        sat2ephem[ephem_ptr->sat].emplace_back(ephem_ptr);
        sat2time_index[ephem_ptr->sat].emplace(toe, sat2ephem.at(ephem_ptr->sat).size()-1);
    }
}

void GNSSProcess::inputIonoParams(double ts, const std::vector<double> &iono_params) // 
{
    if (iono_params.size() != 8)    return;

    // update ionosphere parameters
    latest_gnss_iono_params.swap(empty_vec_d);
    std::copy(iono_params.begin(), iono_params.end(), std::back_inserter(latest_gnss_iono_params));
}

void GNSSProcess::inputpvt(double ts, double lat, double lon, double alt) // 
{
  Eigen::Vector3d lla;
  lla << lat, lon, alt;
  if (pvt_time.empty())
  {
    first_lla_pvt = lla;
    first_xyz_ecef_pvt = geo2ecef(lla);
    cout << "first ecef xyz:" << first_xyz_ecef_pvt.transpose() << endl;
  }
  Eigen::Vector3d xyz_ecef = geo2ecef(lla);
  Eigen::Vector3d xyz_enu = ecef2enu(first_lla_pvt, xyz_ecef - first_xyz_ecef_pvt);
  pvt_time.push_back(ts);
  pvt_holder.push_back(xyz_enu);
}

void GNSSProcess::inputlla(double ts, double lat, double lon, double alt) // 
{
  Eigen::Vector3d lla;
  lla << lat, lon, alt;
  if (lla_time.empty())
  {
    first_lla_lla = lla;
    first_xyz_ecef_lla = geo2ecef(lla);
  }
  Eigen::Vector3d xyz_ecef = geo2ecef(lla);
  Eigen::Vector3d xyz_enu = ecef2enu(first_lla_lla, xyz_ecef - first_xyz_ecef_lla);
  lla_time.push_back(ts);
  lla_holder.push_back(xyz_enu);
}

Eigen::Vector3d GNSSProcess::local2enu(StatesGroupwithGNSS &state)
{
  Eigen::Vector3d enu_pos;
  if (!nolidar)
  {
    Eigen::Matrix3d R_enu_local_;
    R_enu_local_ = Eigen::AngleAxisd(state.yaw_enu_local, Eigen::Vector3d::UnitZ());
    enu_pos = R_enu_local_ * (state.rot_end * Tex_imu_r + state.pos_end - anc_local); // 

    Eigen::Matrix3d R_ecef_enu_ = ecef2rotation(state.anc);
    Eigen::Vector3d ecef_pos_ = state.anc + R_ecef_enu_ * enu_pos;
    // Eigen::Vector3d lla_pos = ecef2geo(first_xyz_enu_pvt);
    enu_pos = ecef2enu(first_lla_pvt, ecef_pos_ - first_xyz_ecef_pvt);
  }
  else
  {
    Eigen::Vector3d pos_r = state.rot_end * Tex_imu_r + state.pos_end;
    // Eigen::Vector3d lla_pos = ecef2geo(first_xyz_enu_pvt);
    enu_pos = ecef2enu(first_lla_pvt, pos_r - first_xyz_ecef_pvt);
  }

  return enu_pos;
}

Eigen::Vector3d GNSSProcess::local2enu(double yaw, Eigen::Vector3d &anc, Eigen::Vector3d &pos)
{
  Eigen::Vector3d enu_pos;
  if (!nolidar)
  {
    Eigen::Matrix3d R_enu_local_;
    R_enu_local_ = Eigen::AngleAxisd(yaw, Eigen::Vector3d::UnitZ());
    enu_pos = R_enu_local_ * (pos - anc_local); // 

    Eigen::Matrix3d R_ecef_enu_ = ecef2rotation(anc);
    Eigen::Vector3d ecef_pos_ = anc + R_ecef_enu_ * enu_pos;
    // Eigen::Vector3d lla_pos = ecef2geo(first_xyz_enu_pvt);
    enu_pos = ecef2enu(first_lla_pvt, ecef_pos_ - first_xyz_ecef_pvt);
  }
  else
  {
    Eigen::Vector3d pos_r = pos;
    // Eigen::Vector3d lla_pos = ecef2geo(first_xyz_enu_pvt);
    enu_pos = ecef2enu(first_lla_pvt, pos_r - first_xyz_ecef_pvt);
  }

  return enu_pos;
}

void GNSSProcess::inputGNSSTimeDiff(const double t_diff) // 
{
    diff_t_gnss_local = t_diff;
}

void GNSSProcess::processGNSS(const std::vector<ObsPtr> &gnss_meas, StatesGroupwithGNSS &state)
{
  std::vector<ObsPtr> valid_meas;
  std::vector<EphemBasePtr> valid_ephems;
  if (gnss_meas.empty())  
  {
    if (gnss_ready)
    {
      gnss_meas_buf[0].swap(empty_vec_o);
      gnss_ephem_buf[0].swap(empty_vec_e);
    }
    return;
  }
  for (auto obs : gnss_meas)
  {
      // filter according to system
      uint32_t sys = satsys(obs->sat, NULL);
      if (sys != SYS_GPS && sys != SYS_GLO && sys != SYS_GAL && sys != SYS_BDS)
          continue;
    
      size_t ephem_index = -1;
      EphemBasePtr best_ephem_cur;
    if (!ephem_from_rinex)
    {
      // if not got cooresponding ephemeris yet
      if (sat2ephem.count(obs->sat) == 0)
          continue;
      
      if (obs->freqs.empty())    continue;       // no valid signal measurement
      // int 
      freq_idx = -1;
      L1_freq(obs, &freq_idx);
      if (freq_idx < 0)   continue;              // no L1 observation
      
      double obs_time = time2sec(obs->time);
      std::map<double, size_t> time2index = sat2time_index.at(obs->sat);
      double ephem_time = EPH_VALID_SECONDS;
      for (auto ti : time2index)
      {
          if (std::abs(ti.first - obs_time) < ephem_time)
          {
              ephem_time = std::abs(ti.first - obs_time);
              ephem_index = ti.second;
          }
      }
      std::map<double, size_t>().swap(time2index);
      if (ephem_time >= EPH_VALID_SECONDS)
      {
          cerr << "ephemeris not valid anymore\n";
          continue;
      }
      best_ephem_cur = sat2ephem.at(obs->sat).at(ephem_index);
    }
    else
    {
      // cout << "gnss ready:" << gnss_ready << endl;
      if (sat2ephem_rnx.count(obs->sat) == 0)
          continue;
      if (obs->freqs.empty())    continue;       // no valid signal measurement
      
      freq_idx = -1;
      L1_freq(obs, &freq_idx);
      if (freq_idx < 0)   continue;              // no L1 observation
      
      double obs_time = time2sec(obs->time);
      std::map<double, size_t> time2index = sat2time_index_rnx.at(obs->sat);
      double ephem_time = EPH_VALID_SECONDS;
      for (auto ti : time2index)
      {
          if (std::abs(ti.first - obs_time) < ephem_time)
          {
              ephem_time = std::abs(ti.first - obs_time);
              ephem_index = ti.second;
          }
      }
      std::map<double, size_t>().swap(time2index);
      if (ephem_time >= EPH_VALID_SECONDS)
      {
          cerr << "ephemeris not valid anymore\n";
          continue;
      }
      best_ephem_cur = sat2ephem_rnx.at(obs->sat).at(ephem_index);
    }
      const EphemBasePtr &best_ephem = best_ephem_cur;
      // filter by tracking status
      LOG_IF(FATAL, freq_idx < 0) << "No L1 observation found.\n";
      if (obs->psr_std[freq_idx]  > gnss_psr_std_threshold ||
          obs->dopp_std[freq_idx] > gnss_dopp_std_threshold) //||
          // obs->cp_std[freq_idx] * 0.004 > gnss_cp_std_threshold)
      {
          sat_track_status[obs->sat] = 0;
          continue;
      }
      else
      {
          if (sat_track_status.count(obs->sat) == 0)
              sat_track_status[obs->sat] = 0;
          ++ sat_track_status[obs->sat];
      }
      if (sat_track_status[obs->sat] < gnss_track_num_threshold)
      {
          continue;           // not being tracked for enough epochs
      }
      // filter by elevation angle
      if (gnss_ready) // && !quick_it) // gnss initialization is completed, then filter the sat by elevation angle // need to be defined
      {
          Eigen::Vector3d sat_ecef;
          if (sys == SYS_GLO)
              sat_ecef = geph2pos(obs->time, std::dynamic_pointer_cast<GloEphem>(best_ephem), NULL);
          else
              sat_ecef = eph2pos(obs->time, std::dynamic_pointer_cast<Ephem>(best_ephem), NULL);
          double azel[2] = {0, M_PI/2.0};
          Eigen::Vector3d pos_gnss = state.pos_end + state.rot_end * Tex_imu_r;
          updateGNSSStatistics(pos_gnss, state.yaw_enu_local, state.anc);
          sat_azel(ecef_pos, sat_ecef, azel); // ecef_pos should be updated for this time step // coarse value is acceptable as well TODO
          if (azel[1] < gnss_elevation_threshold*M_PI/180.0)
              continue;
      }
      valid_meas.push_back(obs);
      valid_ephems.push_back(best_ephem);
  }
  if (!gnss_ready)
  {
    if (valid_meas.empty() || valid_meas.size() < 5) return; // right or not?
    pos_window[frame_count] = state.pos_end + state.rot_end * Tex_imu_r;
    vel_window[frame_count] = state.vel_end;
    gnss_meas_buf[frame_count] = valid_meas; 
    gnss_ephem_buf[frame_count] = valid_ephems;
    frame_count ++;
    gnss_ready = GNSSLIAlign();
  }
  else
  {  
    gnss_meas_buf[0] = valid_meas; 
    gnss_ephem_buf[0] = valid_ephems;
    // optimization();
  }
}

void GNSSProcess::processGNSS(const std::vector<ObsPtr> &gnss_meas, StatesGroupwithGNSS2 &state)
{
  std::vector<ObsPtr> valid_meas;
  std::vector<EphemBasePtr> valid_ephems;
  if (gnss_meas.empty())  
  {
    if (gnss_ready)
    {
      gnss_meas_buf[0].swap(empty_vec_o);
      gnss_ephem_buf[0].swap(empty_vec_e);
    }
    return;
  }
  for (auto obs : gnss_meas)
  {
      // filter according to system
      uint32_t sys = satsys(obs->sat, NULL);
      if (sys != SYS_GPS && sys != SYS_GLO && sys != SYS_GAL && sys != SYS_BDS)
          continue;
    
      size_t ephem_index = -1;
      EphemBasePtr best_ephem_cur;
    if (!ephem_from_rinex)
    {
      // if not got cooresponding ephemeris yet
      if (sat2ephem.count(obs->sat) == 0)
          continue;
      
      if (obs->freqs.empty())    continue;       // no valid signal measurement
      // int 
      freq_idx = -1;
      L1_freq(obs, &freq_idx);
      if (freq_idx < 0)   continue;              // no L1 observation
      
      double obs_time = time2sec(obs->time);
      std::map<double, size_t> time2index = sat2time_index.at(obs->sat);
      double ephem_time = EPH_VALID_SECONDS;
      for (auto ti : time2index)
      {
          if (std::abs(ti.first - obs_time) < ephem_time)
          {
              ephem_time = std::abs(ti.first - obs_time);
              ephem_index = ti.second;
          }
      }
      std::map<double, size_t>().swap(time2index);
      if (ephem_time >= EPH_VALID_SECONDS)
      {
          cerr << "ephemeris not valid anymore\n";
          continue;
      }
      best_ephem_cur = sat2ephem.at(obs->sat).at(ephem_index);
    }
    else
    {
      // cout << "gnss ready:" << gnss_ready << endl;
      if (sat2ephem_rnx.count(obs->sat) == 0)
          continue;
      if (obs->freqs.empty())    continue;       // no valid signal measurement
      
      freq_idx = -1;
      L1_freq(obs, &freq_idx);
      if (freq_idx < 0)   continue;              // no L1 observation
      
      double obs_time = time2sec(obs->time);
      std::map<double, size_t> time2index = sat2time_index_rnx.at(obs->sat);
      double ephem_time = EPH_VALID_SECONDS;
      for (auto ti : time2index)
      {
          if (std::abs(ti.first - obs_time) < ephem_time)
          {
              ephem_time = std::abs(ti.first - obs_time);
              ephem_index = ti.second;
          }
      }
      std::map<double, size_t>().swap(time2index);
      if (ephem_time >= EPH_VALID_SECONDS)
      {
          cerr << "ephemeris not valid anymore\n";
          continue;
      }
      best_ephem_cur = sat2ephem_rnx.at(obs->sat).at(ephem_index);
    }
      const EphemBasePtr &best_ephem = best_ephem_cur;
      // filter by tracking status
      LOG_IF(FATAL, freq_idx < 0) << "No L1 observation found.\n";
      if (gnss_ready)
      {
      if (obs->psr_std[freq_idx]  > gnss_psr_std_threshold ||
          obs->dopp_std[freq_idx] > gnss_dopp_std_threshold) //||
          // obs->cp_std[freq_idx] * 0.004 > gnss_cp_std_threshold)
      {
          sat_track_status[obs->sat] = 0;
          continue;
      }
      else
      {
          if (sat_track_status.count(obs->sat) == 0)
              sat_track_status[obs->sat] = 0;
          ++ sat_track_status[obs->sat];
      }
      }
      else
      {
      if (obs->psr_std[freq_idx]  > gnss_psr_std_threshold ||
          obs->dopp_std[freq_idx] > gnss_dopp_std_threshold) //||
          // obs->cp_std[freq_idx] * 0.004 > gnss_cp_std_threshold)
      {
          sat_track_status[obs->sat] = 0;
          continue;
      }
      else
      {
          if (sat_track_status.count(obs->sat) == 0)
              sat_track_status[obs->sat] = 0;
          ++ sat_track_status[obs->sat];
      }
      }
      if (sat_track_status[obs->sat] < gnss_track_num_threshold)
      {
          continue;           // not being tracked for enough epochs
      }
      // filter by elevation angle
      if (gnss_ready) // && !quick_it) // gnss initialization is completed, then filter the sat by elevation angle // need to be defined
      {
          Eigen::Vector3d sat_ecef;
          if (sys == SYS_GLO)
              sat_ecef = geph2pos(obs->time, std::dynamic_pointer_cast<GloEphem>(best_ephem), NULL);
          else
              sat_ecef = eph2pos(obs->time, std::dynamic_pointer_cast<Ephem>(best_ephem), NULL);
          double azel[2] = {0, M_PI/2.0};
          Eigen::Vector3d pos_gnss = state.pos_end + state.rot_end * Tex_imu_r;
          updateGNSSStatistics(pos_gnss, state.yaw_enu_local, state.anc);
          sat_azel(ecef_pos, sat_ecef, azel); // ecef_pos should be updated for this time step // coarse value is acceptable as well TODO
          if (azel[1] < gnss_elevation_threshold*M_PI/180.0)
              continue;
      }
      valid_meas.push_back(obs);
      valid_ephems.push_back(best_ephem);
  }
  if (!gnss_ready)
  {
    if (valid_meas.empty() || valid_meas.size() < 5) return; // right or not?
    pos_window[frame_count] = state.pos_end + state.rot_end * Tex_imu_r;
    vel_window[frame_count] = state.vel_end;
    gnss_meas_buf[frame_count] = valid_meas; 
    gnss_ephem_buf[frame_count] = valid_ephems;
    frame_count ++;
    gnss_ready = GNSSLIAlign();
  }
  else
  {  
    gnss_meas_buf[0] = valid_meas; 
    gnss_ephem_buf[0] = valid_ephems;
    // optimization();
  }
}

// void GNSSProcess::processGNSS(const std::vector<ObsPtr> &gnss_meas, StatesGroupwithGNSS2 &state)
// {
//   std::vector<ObsPtr> valid_meas;
//   std::vector<ObsPtr> valid_cp_meas;
//   std::vector<EphemBasePtr> valid_ephems;
//   std::vector<EphemBasePtr> valid_cp_ephems;
//   if (gnss_meas.empty())  
//   {
//     if (gnss_ready)
//     {
//       std::vector<ObsPtr>().swap(gnss_meas_buf[0]);
//       std::vector<ObsPtr>().swap(gnss_meas_buf[1]);
//       std::vector<EphemBasePtr>().swap(gnss_ephem_buf[0]);
//       std::vector<EphemBasePtr>().swap(gnss_ephem_buf[1]);
//     }
//     return;
//   }
//   for (auto obs : gnss_meas)
//   {
//       // filter according to system
//       uint32_t sys = satsys(obs->sat, NULL);
//       if (sys != SYS_GPS && sys != SYS_GLO && sys != SYS_GAL && sys != SYS_BDS)
//           continue;
    
//       if (obs->freqs.empty())    continue;  // no valid signal measurement
//       freq_idx = -1;
//       L1_freq(obs, &freq_idx);
//       if (freq_idx < 0)   continue;              // no L1 observation

//       size_t ephem_index = -1;
//       EphemBasePtr best_ephem_cur;
//     if (!ephem_from_rinex)
//     {
//       // if not got cooresponding ephemeris yet
//       if (sat2ephem.count(obs->sat) == 0)
//           continue;
      
//       // if (obs->freqs.empty())    continue;       // no valid signal measurement
//       // int 
//       // freq_idx = -1;
//       // L1_freq(obs, &freq_idx);
//       // if (freq_idx < 0)   continue;              // no L1 observation
      
//       double obs_time = time2sec(obs->time);
//       std::map<double, size_t> time2index = sat2time_index.at(obs->sat);
//       double ephem_time = EPH_VALID_SECONDS;
//       for (auto ti : time2index)
//       {
//           if (std::abs(ti.first - obs_time) < ephem_time)
//           {
//               ephem_time = std::abs(ti.first - obs_time);
//               ephem_index = ti.second;
//           }
//       }
//       std::map<double, size_t>().swap(time2index);
//       if (ephem_time >= EPH_VALID_SECONDS)
//       {
//           cerr << "ephemeris not valid anymore\n";
//           continue;
//       }
//       best_ephem_cur = sat2ephem.at(obs->sat).at(ephem_index);
//     }
//     else
//     {
//       // cout << "gnss ready:" << gnss_ready << endl;
//       if (sat2ephem_rnx.count(obs->sat) == 0)
//           continue;
//       // if (obs->freqs.empty())    continue;       // no valid signal measurement
      
//       // freq_idx = -1;
//       // L1_freq(obs, &freq_idx);
//       // if (freq_idx < 0)   continue;              // no L1 observation
      
//       double obs_time = time2sec(obs->time);
//       std::map<double, size_t> time2index = sat2time_index_rnx.at(obs->sat);
//       double ephem_time = EPH_VALID_SECONDS;
//       for (auto ti : time2index)
//       {
//           if (std::abs(ti.first - obs_time) < ephem_time)
//           {
//               ephem_time = std::abs(ti.first - obs_time);
//               ephem_index = ti.second;
//           }
//       }
//       std::map<double, size_t>().swap(time2index);
//       if (ephem_time >= EPH_VALID_SECONDS)
//       {
//           cerr << "ephemeris not valid anymore\n";
//           continue;
//       }
//       best_ephem_cur = sat2ephem_rnx.at(obs->sat).at(ephem_index);
//     }
//       const EphemBasePtr &best_ephem = best_ephem_cur;
//       if (gnss_ready) // && !quick_it) // gnss initialization is completed, then filter the sat by elevation angle // need to be defined
//       {
//           Eigen::Vector3d sat_ecef;
//           if (sys == SYS_GLO)
//               sat_ecef = geph2pos(obs->time, std::dynamic_pointer_cast<GloEphem>(best_ephem), NULL);
//           else
//               sat_ecef = eph2pos(obs->time, std::dynamic_pointer_cast<Ephem>(best_ephem), NULL);
//           double azel[2] = {0, M_PI/2.0};
//           Eigen::Vector3d pos_gnss = state.pos_end + state.rot_end * Tex_imu_r;
//           updateGNSSStatistics(pos_gnss, state.yaw_enu_local, state.anc);
//           sat_azel(ecef_pos, sat_ecef, azel); // ecef_pos should be updated for this time step // coarse value is acceptable as well TODO
//           if (azel[1] >= gnss_elevation_threshold*M_PI/180.0)
//           {
//             if (obs->cp_std[freq_idx] * 0.004 < gnss_cp_std_threshold && obs->cp[freq_idx] > 0 && obs->psr_std[freq_idx] < gnss_psr_std_threshold)
//             {
//               valid_cp_meas.push_back(obs);
//               valid_cp_ephems.push_back(best_ephem);
//             }
//           }
//       }
//       // filter by tracking status
//       if (!gnss_ready)
//       {
//       LOG_IF(FATAL, freq_idx < 0) << "No L1 observation found.\n";
//       if (!quick_init)
//       {
//       if (obs->psr_std[freq_idx]  > gnss_psr_std_threshold ||
//           obs->dopp_std[freq_idx] > gnss_dopp_std_threshold)
//       {
//           sat_track_status[obs->sat] = 0;
//           continue;
//       }
//       else
//       {
//           if (sat_track_status.count(obs->sat) == 0)
//               sat_track_status[obs->sat] = 0;
//           ++ sat_track_status[obs->sat];
//       }
//       }
//       else
//       {
//         if (obs->psr_std[freq_idx]  > gnss_psr_std_threshold ||
//           obs->dopp_std[freq_idx] > gnss_dopp_std_threshold)
//         {
//           sat_track_status[obs->sat] = 0;
//           continue;
//         }
//         else
//         {
//           if (sat_track_status.count(obs->sat) == 0)
//               sat_track_status[obs->sat] = 0;
//           ++ sat_track_status[obs->sat];
//         }
//       }
//       }
//       else
//       {
//       LOG_IF(FATAL, freq_idx < 0) << "No L1 observation found.\n";
//       if (obs->psr_std[freq_idx] > gnss_psr_std_threshold / 2 ||
//           obs->dopp_std[freq_idx] > gnss_dopp_std_threshold)
//       {
//           sat_track_status[obs->sat] = 0;
//           continue;
//       }
//       else
//       {
//           if (sat_track_status.count(obs->sat) == 0)
//               sat_track_status[obs->sat] = 0;
//           ++ sat_track_status[obs->sat];
//       }
//       }

//       if (sat_track_status[obs->sat] < gnss_track_num_threshold)
//       {
//           continue;           // not being tracked for enough epochs
//       }
//       // filter by elevation angle
//       if (gnss_ready) // && !quick_it) // gnss initialization is completed, then filter the sat by elevation angle // need to be defined
//       {
//           Eigen::Vector3d sat_ecef;
//           if (sys == SYS_GLO)
//               sat_ecef = geph2pos(obs->time, std::dynamic_pointer_cast<GloEphem>(best_ephem), NULL);
//           else
//               sat_ecef = eph2pos(obs->time, std::dynamic_pointer_cast<Ephem>(best_ephem), NULL);
//           double azel[2] = {0, M_PI/2.0};
//           Eigen::Vector3d pos_gnss = state.pos_end + state.rot_end * Tex_imu_r;
//           updateGNSSStatistics(pos_gnss, state.yaw_enu_local, state.anc);
//           sat_azel(ecef_pos, sat_ecef, azel); // ecef_pos should be updated for this time step // coarse value is acceptable as well TODO
//           if (azel[1] < gnss_elevation_threshold*M_PI/180.0)
//               continue;
//       }
//       valid_meas.push_back(obs);
//       valid_ephems.push_back(best_ephem);
//   }
//   if (!gnss_ready)
//   {
//     if (valid_meas.empty() || valid_meas.size() < 5) return; // right or not?
//     pos_window[frame_count] = state.pos_end + state.rot_end * Tex_imu_r;
//     vel_window[frame_count] = state.vel_end;
//     gnss_meas_buf[frame_count] = valid_meas; 
//     gnss_ephem_buf[frame_count] = valid_ephems;
//     frame_count ++;
//     gnss_ready = GNSSLIAlign();
//   }
//   else
//   {  
//     gnss_meas_buf[0] = valid_meas; 
//     gnss_ephem_buf[0] = valid_ephems;
//     gnss_meas_buf[1] = valid_cp_meas; 
//     gnss_ephem_buf[1] = valid_cp_ephems;
//     std::vector<ObsPtr>().swap(valid_cp_meas);
//     std::vector<ObsPtr>().swap(valid_meas);
//     std::vector<EphemBasePtr>().swap(valid_cp_ephems);
//     std::vector<EphemBasePtr>().swap(valid_ephems);
//   }
// }


bool GNSSProcess::GNSSLIAlign()
{
  for (uint32_t i = 0; i < (WINDOW_SIZE+1); ++i)
  {
    if (gnss_meas_buf[i].empty() || gnss_meas_buf[i].size() < 5)
    {
      if (frame_count == WINDOW_SIZE + 1)
      {
        for (uint32_t i = 0; i < (WINDOW_SIZE); ++i)
        {
          gnss_meas_buf[i] = gnss_meas_buf[i+1];
          gnss_ephem_buf[i] = gnss_ephem_buf[i+1];
          pos_window[i] = pos_window[i+1];
          vel_window[i] = vel_window[i+1];
        }
        frame_count = WINDOW_SIZE;

        gnss_meas_buf[frame_count].swap(empty_vec_o);
        gnss_ephem_buf[frame_count].swap(empty_vec_e);              
      }
      return false;
    }
  }

  // check horizontal velocity excitation
  if (!quick_init)
  {
    Eigen::Vector2d avg_hor_vel(0.0, 0.0);
    for (uint32_t i = 0; i < (WINDOW_SIZE+1); ++i)
        avg_hor_vel += vel_window[i].head<2>().cwiseAbs();
    avg_hor_vel /= (WINDOW_SIZE+1);
    if (avg_hor_vel.norm() < 0.3)
    {
      std::cerr << "velocity excitation not enough for GNSS-LI alignment.\n";
      for (uint32_t i = 0; i < (WINDOW_SIZE); ++i)
      {
        gnss_meas_buf[i] = gnss_meas_buf[i+1];
        gnss_ephem_buf[i] = gnss_ephem_buf[i+1];
        pos_window[i] = pos_window[i+1];
        vel_window[i] = vel_window[i+1];
      }
      frame_count = WINDOW_SIZE;
      gnss_meas_buf[frame_count].swap(empty_vec_o);
      gnss_ephem_buf[frame_count].swap(empty_vec_e);    
      return false;
    }
  }
  std::vector<std::vector<ObsPtr>> curr_gnss_meas_buf;
  std::vector<std::vector<EphemBasePtr>> curr_gnss_ephem_buf;
  for (uint32_t i = 0; i < (WINDOW_SIZE+1); ++i)
  {
      curr_gnss_meas_buf.push_back(gnss_meas_buf[i]);
      curr_gnss_ephem_buf.push_back(gnss_ephem_buf[i]);
  }

  GNSSLIInitializer gnss_li_initializer(curr_gnss_meas_buf, curr_gnss_ephem_buf, latest_gnss_iono_params);

  // 1. get a rough global location
  Eigen::Matrix<double, 7, 1> rough_xyzt;
  rough_xyzt.setZero();
  if (!gnss_li_initializer.coarse_localization(rough_xyzt))
  {
      std::cerr << "Fail to obtain a coarse location.\n";
      for (uint32_t i = 0; i < (WINDOW_SIZE); ++i)
      {
        gnss_meas_buf[i] = gnss_meas_buf[i+1];
        gnss_ephem_buf[i] = gnss_ephem_buf[i+1];
        pos_window[i] = pos_window[i+1];
        vel_window[i] = vel_window[i+1];
      }
      frame_count = WINDOW_SIZE;
      gnss_meas_buf[frame_count].swap(empty_vec_o);
      gnss_ephem_buf[frame_count].swap(empty_vec_e);
      return false;
  }

  if (!quick_init)
  {
  // 2. perform yaw alignment
  std::vector<Eigen::Vector3d> local_vs;
  for (uint32_t i = 0; i < (WINDOW_SIZE+1); ++i)
      local_vs.push_back(vel_window[i]); // values at gnss measurement
  Eigen::Vector3d rough_anchor_ecef = rough_xyzt.head<3>();
  double aligned_yaw = 0;
  double aligned_rcv_ddt = 0;
  if (!gnss_li_initializer.yaw_alignment(local_vs, rough_anchor_ecef, aligned_yaw, aligned_rcv_ddt))
  {
      std::cerr << "Fail to align ENU and local frames.\n";
      for (uint32_t i = 0; i < (WINDOW_SIZE); ++i)
      {
        gnss_meas_buf[i] = gnss_meas_buf[i+1];
        gnss_ephem_buf[i] = gnss_ephem_buf[i+1];

        pos_window[i] = pos_window[i+1];
        vel_window[i] = vel_window[i+1];
      }
      frame_count = WINDOW_SIZE;
      gnss_meas_buf[frame_count].swap(empty_vec_o);
      gnss_ephem_buf[frame_count].swap(empty_vec_e);
      return false;
  }
  // std::cout << "aligned_yaw is " << aligned_yaw*180.0/M_PI << '\n';

  // 3. perform anchor refinement
  std::vector<Eigen::Vector3d> local_ps;
  for (uint32_t i = 0; i < (WINDOW_SIZE+1); ++i)
      local_ps.push_back(pos_window[i]); // values at gnss measurement
  Eigen::Matrix<double, 7, 1> refined_xyzt;
  refined_xyzt.setZero();
  if (!gnss_li_initializer.anchor_refinement(local_ps, aligned_yaw, 
      aligned_rcv_ddt, rough_xyzt, refined_xyzt))
  {
      std::cerr << "Fail to refine anchor point.\n";
      for (uint32_t i = 0; i < (WINDOW_SIZE); ++i)
      {
        gnss_meas_buf[i] = gnss_meas_buf[i+1];
        gnss_ephem_buf[i] = gnss_ephem_buf[i+1];

        pos_window[i] = pos_window[i+1];
        vel_window[i] = vel_window[i+1];
      }
      frame_count = WINDOW_SIZE;
      gnss_meas_buf[frame_count].swap(empty_vec_o);
      gnss_ephem_buf[frame_count].swap(empty_vec_e);
      return false;
  }
  // std::cout << "refined anchor point is " << std::setprecision(20) 
  //           << refined_xyzt.head<3>().transpose() << '\n';

  // restore GNSS states
  uint32_t one_observed_sys = static_cast<uint32_t>(-1);
  for (uint32_t k = 0; k < 4; ++k)
  {
      if (rough_xyzt(k+3) != 0) // why use rough value?
      {
          one_observed_sys = k;
          break;
      }
  }
  for (uint32_t i = 0; i < (WINDOW_SIZE+1); ++i)
  {
      para_rcv_ddt[i] = aligned_rcv_ddt;
      for (uint32_t k = 0; k < 4; ++k)
      {
          if (rough_xyzt(k+3) == 0)
              para_rcv_dt[i*4+k] = refined_xyzt(3+one_observed_sys) + aligned_rcv_ddt * (time2sec(gnss_meas_buf[i][0]->time) - time2sec(gnss_meas_buf[0][0]->time)); // why multiply i? should it be dt? how to get dt (TODO)
          else
              para_rcv_dt[i*4+k] = refined_xyzt(3+k) + aligned_rcv_ddt * (time2sec(gnss_meas_buf[i][0]->time) - time2sec(gnss_meas_buf[0][0]->time));
      }
  }
  anc_ecef = refined_xyzt.head<3>();
  anc_local = pos_window[0]; // [WINDOW_SIZE];
  yaw_enu_local = aligned_yaw;
  R_ecef_enu = ecef2rotation(anc_ecef);
  // if (!gnss_li_initializer.yaw_refinement(local_vs, anc_ecef, local_ps, aligned_yaw, aligned_rcv_ddt))
  // {
  //     std::cerr << "Fail to align ENU and local frames.\n";
  //     for (uint32_t i = 0; i < (WINDOW_SIZE); ++i)
  //     {
  //       gnss_meas_buf[i] = gnss_meas_buf[i+1];
  //       gnss_ephem_buf[i] = gnss_ephem_buf[i+1];

  //       pos_window[i] = pos_window[i+1];
  //       vel_window[i] = vel_window[i+1];
  //     }
  //     frame_count = WINDOW_SIZE;
  //     gnss_meas_buf[frame_count].clear();
  //     gnss_ephem_buf[frame_count].clear();
  //     return false;
  // }
  
  // yaw_enu_local = aligned_yaw;
  // for (uint32_t i = 0; i < (WINDOW_SIZE+1); ++i)
  // {
  //     para_rcv_ddt[i] = aligned_rcv_ddt;
  // }
  }
  else
  {
    for (uint32_t k = 0; k < 4; k++)
    {
      if (rough_xyzt(3+k) == 0)
      {
        std::cerr << "Fail to quick init anchor point.\n";
        for (uint32_t i = 0; i < (WINDOW_SIZE); ++i)
        {
          gnss_meas_buf[i] = gnss_meas_buf[i+1]; // change the strategy
          gnss_ephem_buf[i] = gnss_ephem_buf[i+1];

          pos_window[i] = pos_window[i+1];
          vel_window[i] = vel_window[i+1];
        }
        frame_count = WINDOW_SIZE;
        gnss_meas_buf[frame_count].swap(empty_vec_o);
        gnss_ephem_buf[frame_count].swap(empty_vec_e);
        return false;
      }
    }
    anc_ecef = rough_xyzt.head<3>();
    anc_local = pos_window[0]; // [WINDOW_SIZE];
    yaw_enu_local = 0.0;
    R_ecef_enu = ecef2rotation(anc_ecef);
    para_rcv_ddt[WINDOW_SIZE] = 128.0;
    for (uint32_t k_ = 0; k_ < 4; ++k_)
    {
      para_rcv_dt[WINDOW_SIZE*4+k_] = rough_xyzt(3+k_);
    }
    for (uint32_t k_ = 2; k_ < WINDOW_SIZE+1; k_++)
    {
      std::vector<ObsPtr>().swap(gnss_meas_buf[k_]);
      std::vector<EphemBasePtr>().swap(gnss_ephem_buf[k_]);
    }
    return true;
  }
    for (uint32_t k_ = 2; k_ < WINDOW_SIZE+1; k_++)
    {
      std::vector<ObsPtr>().swap(gnss_meas_buf[k_]);
      std::vector<EphemBasePtr>().swap(gnss_ephem_buf[k_]);
    }
  return true;
}

void GNSSProcess::updateGNSSStatistics(Eigen::Vector3d &pos, double yaw_cur, Eigen::Vector3d &anc_cur) // delete
{
  if (!nolidar)
  {
    Eigen::Matrix3d R_enu_local_;
    R_enu_local_ = Eigen::AngleAxisd(yaw_cur, Eigen::Vector3d::UnitZ());
    Eigen::Vector3d enu_pos = R_enu_local_ * (pos - anc_local);
    R_ecef_enu = ecef2rotation(anc_cur);
    ecef_pos = anc_cur + R_ecef_enu * enu_pos;
  }
  else
  {
    ecef_pos = pos;
  }
}

void GNSSProcess::GnssPsrDoppMeas(const ObsPtr &obs_, const EphemBasePtr &ephem_) 
{
  freq = L1_freq(obs_, &freq_idx);
  LOG_IF(FATAL, freq < 0) << "No L1 observation found.";

  uint32_t sys = satsys(obs_->sat, NULL);
  double tof = obs_->psr[freq_idx] / LIGHT_SPEED;
  gtime_t sv_tx = time_add(obs_->time, -tof);

  if (sys == SYS_GLO)
  {
      GloEphemPtr glo_ephem = std::dynamic_pointer_cast<GloEphem>(ephem_);
      svdt = geph2svdt(sv_tx, glo_ephem);
      sv_tx = time_add(sv_tx, -svdt);
      sv_pos = geph2pos(sv_tx, glo_ephem, &svdt);
      sv_vel = geph2vel(sv_tx, glo_ephem, &svddt);
      tgd = 0.0;
      pr_uura = 2.0 * (obs_->psr_std[freq_idx]/0.16);
      dp_uura = 2.0 * (obs_->dopp_std[freq_idx]/0.256);
  }
  else
  {
      EphemPtr eph = std::dynamic_pointer_cast<Ephem>(ephem_);
      svdt = eph2svdt(sv_tx, eph); // used in eva
      sv_tx = time_add(sv_tx, -svdt);
      sv_pos = eph2pos(sv_tx, eph, &svdt); // used in eva
      sv_vel = eph2vel(sv_tx, eph, &svddt); // used in eva
      tgd = eph->tgd[0];
      if (sys == SYS_GAL)
      {
          pr_uura = (eph->ura - 2.0) * (obs_->psr_std[freq_idx]/0.16);
          dp_uura = (eph->ura - 2.0) * (obs_->dopp_std[freq_idx]/0.256);
      }
      else
      {
          pr_uura = (eph->ura - 1.0) * (obs_->psr_std[freq_idx]/0.16);
          dp_uura = (eph->ura - 1.0) * (obs_->dopp_std[freq_idx]/0.256);
      }
  }
  LOG_IF(FATAL, pr_uura <= 0) << "pr_uura is " << pr_uura; // get those parameters mainly, both used in eva
  LOG_IF(FATAL, dp_uura <= 0) << "dp_uura is " << dp_uura;
  relative_sqrt_info = 1.0;
}

void GNSSProcess::SvPosCals(const ObsPtr &obs_, const EphemBasePtr &ephem_) 
{
  freq = L1_freq(obs_, &freq_idx);
  LOG_IF(FATAL, freq < 0) << "No L1 observation found.";

  uint32_t sys = satsys(obs_->sat, NULL);
  double tof = obs_->psr[freq_idx] / LIGHT_SPEED;
  gtime_t sv_tx = time_add(obs_->time, -tof);

  if (sys == SYS_GLO)
  {
      GloEphemPtr glo_ephem = std::dynamic_pointer_cast<GloEphem>(ephem_);
      svdt = geph2svdt(sv_tx, glo_ephem);
      sv_tx = time_add(sv_tx, -svdt);
      sv_pos = geph2pos(sv_tx, glo_ephem, &svdt);
      sv_vel = geph2vel(sv_tx, glo_ephem, &svddt);
  }
  else
  {
      EphemPtr eph = std::dynamic_pointer_cast<Ephem>(ephem_);
      svdt = eph2svdt(sv_tx, eph); // used in eva
      sv_tx = time_add(sv_tx, -svdt);
      sv_pos = eph2pos(sv_tx, eph, &svdt); // used in eva
      sv_vel = eph2vel(sv_tx, eph, &svddt); // used in eva
  }
}

bool GNSSProcess::Evaluate(StatesGroupwithGNSS &state, Eigen::VectorXd &residuals, Eigen::MatrixXd &jacobians, int &meas_size)
{
  if (gnss_meas_buf[0].empty()) // || gnss_meas_buf[0].size() < 4
  {
    residuals.setZero();
    jacobians.setZero();
    return false;
  }
  double rcv_dt[4];
  rcv_dt[0] = state.dt_g;
  rcv_dt[1] = state.dt_r;
  rcv_dt[2] = state.dt_e;
  rcv_dt[3] = state.dt_c;
  double rcv_ddt = state.ddt;
  double yaw_diff = state.yaw_enu_local;
  Eigen::Vector3d ref_ecef;
  if (nolidar)
  {
    ref_ecef = state.pos_end + state.rot_end * Tex_imu_r;
  }
  else
  {
   ref_ecef = state.anc;
  }
  double s1, s2, e2, a, ep, p, h, lat, ds1dx, ds2dx, ds1dy, ds2dy, ds1dz, ds2dz, sins, coss;
  e2 = EARTH_ECCE_2;
  a = EARTH_SEMI_MAJOR; // _glo?
  ep = ref_ecef(0)*ref_ecef(0) + ref_ecef(1) * ref_ecef(1);
  p = a*a*(1-e2);
  h = ref_ecef(2)*ref_ecef(2)*a*a;
  s1 = ref_ecef(2) + e2/(1-e2) * sqrt(p) * pow(ref_ecef(2)*a/sqrt(h+ep*p),3);
  s2 = sqrt(ep) - a * e2 * pow((ep*p)/(h+ep*p),1.5);
  lat = atan(s1/s2);
  sins = -s1/(s1*s1+s2*s2);
  coss = s2/(s1*s1+s2*s2);
  Eigen::Vector3d R1TE3, E1, vecP, vecV, vecLon, vecLat;
  R1TE3 << 0.0, -sin(lat-M_PI/2), cos(lat-M_PI/2);
  E1 << 1.0, 0.0, 0.0;

  ds1dx = e2/(1-e2) * sqrt(p) * a * ref_ecef(2) * h * (-3) * p * ref_ecef(0) / pow(h+ep*p,2.5);
  ds1dy = e2/(1-e2) * sqrt(p) * a * ref_ecef(2) * h * (-3) * p * ref_ecef(1) / pow(h+ep*p,2.5);
  ds1dz = 1 + e2/(1-e2) * sqrt(p) * 3 *sqrt(h/(h+ep*p)) * a * a * ref_ecef(2) * ep * p / pow(h+ep*p,2);

  ds2dx = ref_ecef(0) / sqrt(ep) - a * e2 * pow(p,1.5) * 3 * sqrt(ep)*ref_ecef(0)*h/pow(h+ep*p,2.5);
  ds2dy = ref_ecef(1) / sqrt(ep) - a * e2 * pow(p,1.5) * 3 * sqrt(ep)*ref_ecef(1)*h/pow(h+ep*p,2.5);
  ds2dz = a*e2*3 * pow(p,1.5) * a * a * ref_ecef(2) * pow(ep, 1.5) / pow(h+ep*p, 2.5);

  vecLon << -ref_ecef(1)/ep, ref_ecef(0)/ep, 0.0;
  vecLat << coss * ds1dx + sins * ds2dx, coss * ds1dy + sins * ds2dy, coss * ds1dz + sins * ds2dz;

  const Eigen::Vector3d local_pos = state.rot_end * Tex_imu_r + state.pos_end;
  const Eigen::Vector3d local_vel = state.vel_end; // check if it is right

  double sin_yaw_diff = std::sin(yaw_diff);
  double cos_yaw_diff = std::cos(yaw_diff);
  Eigen::Matrix3d R_enu_local;
  R_enu_local << cos_yaw_diff, -sin_yaw_diff, 0,
                  sin_yaw_diff,  cos_yaw_diff, 0,
                  0           ,  0           , 1;
  Eigen::Matrix3d R_ecef_enu_cur = ecef2rotation(ref_ecef); // provide anchor value
  Eigen::Matrix3d R_ecef_local = R_ecef_enu_cur * R_enu_local;

  Eigen::Vector3d P_ecef;
  if (nolidar)
  {
    P_ecef = local_pos;
  }
  else
  {
   P_ecef = R_ecef_local * (local_pos - anc_local) + ref_ecef;
  }
  Eigen::Vector3d V_ecef;
  if (nolidar)
  {
   V_ecef = local_vel;
  }
  else
  {
    V_ecef = R_ecef_local * local_vel;
  }

  vecP = R_enu_local * (local_pos - anc_local);
  vecV = R_enu_local * local_vel;
  Eigen::Matrix3d hatP, hatV;
  hatP << 0.0, -vecP(2), vecP(1),
          vecP(2), 0.0, -vecP(0),
          -vecP(1), vecP(0), 0.0;
  hatV << 0.0, -vecV(2), vecV(1),
          vecV(2), 0.0, -vecV(0),
          -vecV(1), vecV(0), 0.0;

  const std::vector<ObsPtr> &curr_obs = gnss_meas_buf[0];
  const std::vector<EphemBasePtr> &curr_ephem = gnss_ephem_buf[0];

  std::map<double, std::map<uint32_t, double[3]> >::iterator it;
  double time_list[curr_obs.size()];
  double time_min = time2sec(curr_obs[0]->time);
  for (uint32_t j = 0; j < curr_obs.size(); ++j)
  {
    time_list[j] = 0.0;
    for (it = sat2cp.begin(); it != sat2cp.end(); ++it)
    {
      std::map<uint32_t, double[3]>::iterator it_sat;
      it_sat = it->second.find(curr_obs[j]->sat);
      if (it_sat != it->second.end())
      {
        time_list[j] = it->first;
        if (time_list[j] < time_min)
        {
          time_min = time_list[j];
        }
        break;
      }
    }
  }

  int best_sat = -1;
  double min_cp_std = 100000000;
  for (uint32_t j = 0; j < curr_obs.size(); ++j)
  {
    if (time_list[j] == time_min)
    {
      freq = L1_freq(curr_obs[j], &freq_idx);
      LOG_IF(FATAL, freq < 0) << "No L1 observation found."; 
      if (curr_obs[j]->cp_std[freq_idx] * 0.004 < min_cp_std && curr_obs[j]->cp[freq_idx] > 0)
      {
        min_cp_std = curr_obs[j]->cp_std[freq_idx] * 0.004;
        best_sat = j;
      }
    }
  }

  std::deque<uint32_t> pair_sat_copy;

  std::deque<double> meas_sats;
  std::deque<double> meas_cov_sats;
  std::deque<double> meas_time_sats;
  std::deque<double> meas_sats_copy;
  std::deque<double> meas_cov_sats_copy;
  // std::deque<double> meas_time_sats_copy;
  
  if (best_sat > -1)
  {
    std::deque<uint32_t> pair_sat;

    for (uint32_t j = 0; j < curr_obs.size(); ++j)
    {
      if (j != best_sat)
      {
        uint32_t sys = satsys(curr_obs[j]->sat, NULL);
        freq = L1_freq(curr_obs[j], &freq_idx);

        if (sys == satsys(curr_obs[best_sat]->sat, NULL) && curr_obs[j]->cp_std[freq_idx] * 0.004 < gnss_cp_std_threshold)
        {
          pair_sat.push_back(j);
        }
      }
    }

    for (uint32_t j = 0; j < curr_obs.size(); ++j)
    {
      std::map<uint32_t, double[3]>::iterator it_old_best, it_old;
      std::vector<double> meas;
      std::vector<double> meas_cov;
      std::vector<double> meas_time;
      if (j == pair_sat.front())
      {
        for (it = sat2cp.begin(); it != sat2cp.end(); ++it)
        {
          it_old = it->second.find(curr_obs[j]->sat);
          it_old_best = it->second.find(curr_obs[best_sat]->sat);
          if (it_old != it->second.end() && it_old_best != it->second.end())
          {
            meas.push_back(it_old_best->second[0] - it_old->second[0] - it_old_best->second[1] + it_old->second[1]);
            meas_cov.push_back(it_old_best->second[2] * it_old_best->second[2] + it_old->second[2] * it_old->second[2]);
            meas_time.push_back(it->first);
          }
        }
        if (meas.size() > 0)
        {
          double sum = std::accumulate(meas.begin(), meas.end(), 0.0);
          double mean = sum / meas.size(); //  meas[meas.size()-1]; //
          std::vector<double> empty_vec_d;
          meas.swap(empty_vec_d);
          meas_sats.push_back(mean);
          meas_cov_sats.push_back(meas_cov[meas_cov.size()-1]);
          meas_time_sats.push_back(meas_time[meas_time.size()-1]);
          meas_cov.swap(empty_vec_d);
          meas_time.swap(empty_vec_d);
          pair_sat_copy.push_back(pair_sat.front());
        }
        pair_sat.pop_front();
      }
    }
  }

  Eigen::MatrixXd jacobians_cur(curr_obs.size() * 2, jacobians.cols());
  Eigen::VectorXd residuals_cur(curr_obs.size() * 2);
  residuals_cur.setZero();
  jacobians_cur.setZero();
  std::map<uint32_t, double[3]> curr_cp_map;

  std::vector<double> meas_cp;
  std::vector<double> esti_cp;
  std::vector<double> cov_cp;
  std::vector<double> delta_t_cp;
  std::vector<Eigen::Vector3d> jaco_cp;
  double meas_cp_best, esti_cp_best, cov_cp_best;
  Eigen::Vector3d jaco_cp_best;

  for (uint32_t j = 0; j < curr_obs.size(); ++j)
  {
    const uint32_t sys = satsys(curr_obs[j]->sat, NULL);
    const uint32_t sys_idx = gnss_comm::sys2idx.at(sys);
    GnssPsrDoppMeas(curr_obs[j], curr_ephem[j]); //, latest_gnss_iono_params);
  
    double ion_delay = 0, tro_delay = 0;
    double azel[2] = {0, M_PI/2.0};
    if (P_ecef.norm() > 0)
    {
        sat_azel(P_ecef, sv_pos, azel);
        Eigen::Vector3d rcv_lla = ecef2geo(P_ecef);
        tro_delay = calculate_trop_delay(curr_obs[j]->time, rcv_lla, azel);
        ion_delay = calculate_ion_delay(curr_obs[j]->time, latest_gnss_iono_params, rcv_lla, azel); // rely on local pose
    }
    double sin_el = sin(azel[1]);
    double sin_el_2 = sin_el*sin_el;
    double pr_weight = sin_el_2 / pr_uura * relative_sqrt_info; // not requisite
    double dp_weight = sin_el_2 / dp_uura * relative_sqrt_info * PSR_TO_DOPP_RATIO; // not requisite

    Eigen::Vector3d rcv2sat_ecef = sv_pos - P_ecef;
    Eigen::Vector3d rcv2sat_unit = rcv2sat_ecef.normalized();

    freq = L1_freq(curr_obs[j], &freq_idx);
    
    const double wavelength = LIGHT_SPEED / freq;

    if (curr_obs[j]->cp_std[freq_idx] * 0.004 < gnss_cp_std_threshold)
    {
      if (curr_obs[j]->cp[freq_idx] > 0)
      {
        curr_cp_map[curr_obs[j]->sat][0] = curr_obs[j]->cp[freq_idx] * wavelength;
        curr_cp_map[curr_obs[j]->sat][1] = rcv2sat_ecef.norm();
        curr_cp_map[curr_obs[j]->sat][2] = curr_obs[j]->cp_std[freq_idx] * 0.004;
      }
    }

    if (j == best_sat)
    {
      meas_cp_best = curr_obs[j]->cp[freq_idx] * wavelength;
      esti_cp_best = rcv2sat_ecef.norm();
      jaco_cp_best = rcv2sat_unit;
      cov_cp_best  = curr_obs[j]->cp_std[freq_idx] * 0.004 * curr_obs[j]->cp_std[freq_idx] * 0.004;
    }

    if (j == pair_sat_copy.front())
    {
      if (curr_obs[j]->cp[freq_idx] > 0)
      {
      meas_cp.push_back(curr_obs[j]->cp[freq_idx] * wavelength);
      cov_cp.push_back(curr_obs[j]->cp_std[freq_idx] * curr_obs[j]->cp_std[freq_idx] * 0.004 * 0.004);
      delta_t_cp.push_back(time2sec(curr_obs[j]->time) - meas_time_sats.front());
      esti_cp.push_back(rcv2sat_ecef.norm());
      jaco_cp.push_back(rcv2sat_unit);
      meas_sats_copy.push_back(meas_sats.front());
      meas_cov_sats_copy.push_back(meas_cov_sats.front());
      }
      pair_sat_copy.pop_front();
      meas_sats.pop_front();
      meas_cov_sats.pop_front();
      meas_time_sats.pop_front();
    }

    const double psr_sagnac = EARTH_OMG_GPS*(sv_pos(0)*P_ecef(1)-sv_pos(1)*P_ecef(0))/LIGHT_SPEED;
    double psr_estimated = rcv2sat_ecef.norm() + psr_sagnac + rcv_dt[sys_idx] - svdt*LIGHT_SPEED + // why not multiply light_speed?  
                              ion_delay + tro_delay + tgd*LIGHT_SPEED;
    const double dopp_sagnac = EARTH_OMG_GPS/LIGHT_SPEED*(sv_vel(0)*P_ecef(1)+
              sv_pos(0)*V_ecef(1) - sv_vel(1)*P_ecef(0) - sv_pos(1)*V_ecef(0));
    double dopp_estimated = (sv_vel - V_ecef).dot(rcv2sat_unit) + rcv_ddt + dopp_sagnac - svddt*LIGHT_SPEED;
    // cout << "CHECK VALUE:" << ion_delay << ";" << tro_delay << ";" << svdt*LIGHT_SPEED << ";" << rcv_dt[sys_idx] << ";" << tgd*LIGHT_SPEED << ";" << psr_sagnac << endl;
    // cout << "CHECK VALUE dp:" << rcv_ddt << ";" << svddt*LIGHT_SPEED << ";" << dopp_sagnac << ";" << state.anc.transpose() << ";" << state.yaw_enu_local << endl;
    Eigen::Matrix3d hat_T;
    hat_T << SKEW_SYM_MATRX(Tex_imu_r);
    if (!nolidar)
    {    
      residuals_cur[2*j] = -(psr_estimated - curr_obs[j]->psr[freq_idx]) * pr_weight;

      residuals_cur[2*j+1] = -(dopp_estimated + curr_obs[j]->dopp[freq_idx]*wavelength) * dp_weight; // / wavelength; // 
      // J_Pi

      jacobians_cur.block<1, 3>(2*j,3) = -rcv2sat_unit.transpose() * R_ecef_local * pr_weight; // why use normalized result?
      jacobians_cur.block<1, 3>(2*j,0) = rcv2sat_unit.transpose() * R_ecef_local * hat_T * pr_weight; // why use normalized result?
      jacobians_cur.block<1, 3>(2*j,24) = -rcv2sat_unit.transpose() * (Eye3d + R_ecef_enu_cur * hatP * R1TE3 * vecLon.transpose() - R_ecef_enu_cur * hatP * E1 * vecLat.transpose()) * pr_weight; // simplify

      const double norm3 = pow(rcv2sat_ecef.norm(), 3);
      const double norm2 = rcv2sat_ecef.squaredNorm();
      Eigen::Matrix3d unit2rcv_pos;
      for (size_t i = 0; i < 3; ++i)
      {
          for (size_t k = 0; k < 3; ++k)
          {
              if (i == k)
                  unit2rcv_pos(i, k) = (norm2-rcv2sat_ecef(i)*rcv2sat_ecef(i))/norm3;
              else
                  unit2rcv_pos(i, k) = (-rcv2sat_ecef(i)*rcv2sat_ecef(k))/norm3;
          }
      }
      unit2rcv_pos *= -1;
      jacobians_cur.block<1, 3>(2*j+1,3) = (sv_vel-V_ecef).transpose() * unit2rcv_pos * 
          R_ecef_local * dp_weight; //  / wavelength; //
      jacobians_cur.block<1, 3>(2*j+1,0) = -(sv_vel-V_ecef).transpose() * unit2rcv_pos * 
          R_ecef_local * hat_T * dp_weight; //  / wavelength; //
      jacobians_cur.block<1, 3>(2*j+1,24) = (sv_vel-V_ecef).transpose() * unit2rcv_pos * (Eye3d + R_ecef_enu_cur * hatP * R1TE3 * vecLon.transpose() - R_ecef_enu_cur * hatP * E1 * vecLat.transpose()) * dp_weight
                                      + rcv2sat_unit.transpose() * (-1.0) * (R_ecef_enu_cur * hatV * R1TE3 * vecLon.transpose() - R_ecef_enu_cur * hatV * E1 * vecLat.transpose()) * dp_weight;
    

      // J_Vi
      jacobians_cur.block<1, 3>(2*j+1,6) = rcv2sat_unit.transpose() * (-1.0) * 
          R_ecef_local * dp_weight; //  / wavelength; // 
      // J_rcv_dt
      jacobians_cur(2*j,19+sys_idx) = 1.0 * pr_weight;

      // J_rcv_ddt
      jacobians_cur(2*j+1,23) = 1.0 * dp_weight; //  / wavelength; //

      // J_yaw_diff
      Eigen::Matrix3d d_yaw;
      d_yaw << 0, -1, 0, 
                1, 0, 0, 
                0, 0, 0;
      jacobians_cur(2*j,18) = -rcv2sat_unit.dot( (R_ecef_local * d_yaw * (local_pos - anc_local))) * pr_weight;
      jacobians_cur(2*j+1,18) = -rcv2sat_unit.dot((R_ecef_local * d_yaw * local_vel)) * dp_weight + (sv_vel-V_ecef).dot( unit2rcv_pos * 
          R_ecef_local * d_yaw * (local_pos - anc_local)) * dp_weight; //  / wavelength; //
    }
    else
    {    
      residuals_cur[2*j] = -(psr_estimated - curr_obs[j]->psr[freq_idx]) * pr_weight;

      
      residuals_cur[2*j+1] = -(dopp_estimated + curr_obs[j]->dopp[freq_idx]*wavelength) * dp_weight; // / wavelength; // 
      // J_Pi

      jacobians_cur.block<1, 3>(2*j,3) = -rcv2sat_unit.transpose() * pr_weight; 
      jacobians_cur.block<1, 3>(2*j,0) = rcv2sat_unit.transpose() * hat_T * pr_weight; 

      const double norm3 = pow(rcv2sat_ecef.norm(), 3);
      const double norm2 = rcv2sat_ecef.squaredNorm();
      Eigen::Matrix3d unit2rcv_pos;
      for (size_t i = 0; i < 3; ++i)
      {
          for (size_t k = 0; k < 3; ++k)
          {
              if (i == k)
                  unit2rcv_pos(i, k) = (norm2-rcv2sat_ecef(i)*rcv2sat_ecef(i))/norm3;
              else
                  unit2rcv_pos(i, k) = (-rcv2sat_ecef(i)*rcv2sat_ecef(k))/norm3;
          }
      }
      unit2rcv_pos *= -1;
      jacobians_cur.block<1, 3>(2*j+1,3) = (sv_vel-V_ecef).transpose() * unit2rcv_pos * dp_weight; //  / wavelength; //
      jacobians_cur.block<1, 3>(2*j+1,0) = -(sv_vel-V_ecef).transpose() * unit2rcv_pos * hat_T * dp_weight; //  / wavelength; //

      // J_Vi
      jacobians_cur.block<1, 3>(2*j+1,6) = rcv2sat_unit.transpose() * (-1.0) * dp_weight; //  / wavelength; // 
      // J_rcv_dt
      jacobians_cur(2*j,19+sys_idx) = 1.0 * pr_weight;

      // J_rcv_ddt
      jacobians_cur(2*j+1,23) = 1.0 * dp_weight; //  / wavelength; //
    }
  }

  residuals.resize(curr_obs.size() * 2 + meas_sats_copy.size());
  jacobians.resize(curr_obs.size()* 2 + meas_sats_copy.size(), jacobians.cols());
  meas_size = curr_obs.size() * 2 + meas_sats_copy.size();
  residuals.setZero();
  jacobians.setZero();
  for (int j = 0; j < curr_obs.size() * 2; j++)
  {
    residuals(j) = residuals_cur(j);
    jacobians.block<1,27>(j,0) = jacobians_cur.block<1,27>(j,0);
  }

  double cp_weight = 0.10;
  for (uint32_t j = 0; j < meas_sats_copy.size(); ++j)
  {
    // cp_weight = 0.01 / (cov_cp_best + cov_cp[j] + meas_cov_sats_copy[j]) / (1+2*delta_t_cp[j]+delta_t_cp[j]*delta_t_cp[j]);
    // cout << "check cp std:" << cp_weight << endl;
    // cout << "check cp:" << meas_sats_copy[j] << ";" << meas_cp_best << ";" << meas_cp[j] << ";" << esti_cp_best << ";" << esti_cp[j] << endl;
    residuals[2*curr_obs.size() + j] = -((meas_sats_copy[j] - (meas_cp_best - meas_cp[j]) + (esti_cp_best - esti_cp[j]))) * cp_weight;
    Eigen::Matrix3d hat_T;
    hat_T << SKEW_SYM_MATRX(Tex_imu_r);
    if (!nolidar)
    {
      jacobians.block<1, 3>(2*curr_obs.size() + j,3) = (-jaco_cp_best.transpose() + jaco_cp[j].transpose()) * R_ecef_local * cp_weight; // why use normalized result?
      jacobians.block<1, 3>(2*curr_obs.size() + j,0) = (jaco_cp_best.transpose() - jaco_cp[j].transpose()) * R_ecef_local * hat_T * cp_weight; // why use normalized result?
      
      jacobians.block<1, 3>(2*curr_obs.size() + j,24) = (-jaco_cp_best.transpose() + jaco_cp[j].transpose()) * (Eye3d + R_ecef_enu_cur * hatP * R1TE3 * vecLon.transpose() - R_ecef_enu_cur * hatP * E1 * vecLat.transpose()) * cp_weight; // simplify

      Eigen::Matrix3d d_yaw;
      d_yaw << 0, -1, 0, 
                1, 0, 0, 
                0, 0, 0;
      jacobians(2*curr_obs.size() + j,18) = (-jaco_cp_best + jaco_cp[j]).dot( (R_ecef_local * d_yaw * (local_pos - anc_local))) * cp_weight;
    }
    else
    {
      jacobians.block<1, 3>(2*curr_obs.size() + j,3) = (-jaco_cp_best.transpose() + jaco_cp[j].transpose()) * cp_weight; // why use normalized result?
      jacobians.block<1, 3>(2*curr_obs.size() + j,0) = (jaco_cp_best.transpose() - jaco_cp[j].transpose()) * hat_T * cp_weight; // why use normalized result?
    }
  }
  sat2cp[time2sec(curr_obs[0]->time)] = curr_cp_map;
  std::map<double, std::map<uint32_t, double[3]>>::iterator it_old;
  if (!sat2cp.empty())
  {
    it_old = sat2cp.begin();
    if (time2sec(curr_obs[0]->time) - it_old->first > gnss_cp_time_threshold)
    {
      std::map<uint32_t, double[3]>().swap(it_old->second);
      size_t del_size = sat2cp.erase(it_old->first);
    }
  }
  return true;
}

bool GNSSProcess::Evaluate(StatesGroupwithGNSS2 &state, Eigen::VectorXd &residuals, Eigen::MatrixXd &jacobians, int &meas_size)
{
  if (gnss_meas_buf[0].empty()) // || gnss_meas_buf[0].size() < 4
  {
    residuals.setZero();
    jacobians.setZero();
    return false;
  }
  double rcv_dt[4];
  rcv_dt[0] = state.dt_g;
  rcv_dt[1] = state.dt_r;
  rcv_dt[2] = state.dt_e;
  rcv_dt[3] = state.dt_c;
  double rcv_ddt = state.ddt;
  double yaw_diff = state.yaw_enu_local;
  Eigen::Vector3d ref_ecef;
  if (nolidar)
  {
    ref_ecef = state.pos_end + state.rot_end * Tex_imu_r;
  }
  else
  {
   ref_ecef = state.anc;
  }
  double s1, s2, e2, a, ep, p, h, lat, ds1dx, ds2dx, ds1dy, ds2dy, ds1dz, ds2dz, sins, coss;
  e2 = EARTH_ECCE_2;
  a = EARTH_SEMI_MAJOR; // _glo?
  ep = ref_ecef(0)*ref_ecef(0) + ref_ecef(1) * ref_ecef(1);
  p = a*a*(1-e2);
  h = ref_ecef(2)*ref_ecef(2)*a*a;
  s1 = ref_ecef(2) + e2/(1-e2) * sqrt(p) * pow(ref_ecef(2)*a/sqrt(h+ep*p),3);
  s2 = sqrt(ep) - a * e2 * pow((ep*p)/(h+ep*p),1.5);
  lat = atan(s1/s2);
  sins = -s1/(s1*s1+s2*s2);
  coss = s2/(s1*s1+s2*s2);
  Eigen::Vector3d R1TE3, E1, vecP, vecV, vecLon, vecLat;
  R1TE3 << 0.0, -sin(lat-M_PI/2), cos(lat-M_PI/2);
  E1 << 1.0, 0.0, 0.0;

  ds1dx = e2/(1-e2) * sqrt(p) * a * ref_ecef(2) * h * (-3) * p * ref_ecef(0) / pow(h+ep*p,2.5);
  ds1dy = e2/(1-e2) * sqrt(p) * a * ref_ecef(2) * h * (-3) * p * ref_ecef(1) / pow(h+ep*p,2.5);
  ds1dz = 1 + e2/(1-e2) * sqrt(p) * 3 *sqrt(h/(h+ep*p)) * a * a * ref_ecef(2) * ep * p / pow(h+ep*p,2);

  ds2dx = ref_ecef(0) / sqrt(ep) - a * e2 * pow(p,1.5) * 3 * sqrt(ep)*ref_ecef(0)*h/pow(h+ep*p,2.5);
  ds2dy = ref_ecef(1) / sqrt(ep) - a * e2 * pow(p,1.5) * 3 * sqrt(ep)*ref_ecef(1)*h/pow(h+ep*p,2.5);
  ds2dz = a*e2*3 * pow(p,1.5) * a * a * ref_ecef(2) * pow(ep, 1.5) / pow(h+ep*p, 2.5);

  vecLon << -ref_ecef(1)/ep, ref_ecef(0)/ep, 0.0;
  vecLat << coss * ds1dx + sins * ds2dx, coss * ds1dy + sins * ds2dy, coss * ds1dz + sins * ds2dz;

  const Eigen::Vector3d local_pos = state.rot_end * Tex_imu_r + state.pos_end;
  const Eigen::Vector3d local_vel = state.vel_end; // check if it is right

  double sin_yaw_diff = std::sin(yaw_diff);
  double cos_yaw_diff = std::cos(yaw_diff);
  Eigen::Matrix3d R_enu_local;
  R_enu_local << cos_yaw_diff, -sin_yaw_diff, 0,
                  sin_yaw_diff,  cos_yaw_diff, 0,
                  0           ,  0           , 1;
  Eigen::Matrix3d R_ecef_enu_cur = ecef2rotation(ref_ecef); // provide anchor value
  Eigen::Matrix3d R_ecef_local = R_ecef_enu_cur * R_enu_local;

  Eigen::Vector3d P_ecef;
  if (nolidar)
  {
    P_ecef = local_pos;
  }
  else
  {
   P_ecef = R_ecef_local * (local_pos - anc_local) + ref_ecef;
  }
  Eigen::Vector3d V_ecef;
  if (nolidar)
  {
   V_ecef = local_vel;
  }
  else
  {
    V_ecef = R_ecef_local * local_vel;
  }

  vecP = R_enu_local * (local_pos - anc_local);
  vecV = R_enu_local * local_vel;
  Eigen::Matrix3d hatP, hatV;
  hatP << 0.0, -vecP(2), vecP(1),
          vecP(2), 0.0, -vecP(0),
          -vecP(1), vecP(0), 0.0;
  hatV << 0.0, -vecV(2), vecV(1),
          vecV(2), 0.0, -vecV(0),
          -vecV(1), vecV(0), 0.0;

  const std::vector<ObsPtr> &curr_obs = gnss_meas_buf[0];
  const std::vector<EphemBasePtr> &curr_ephem = gnss_ephem_buf[0];

  std::map<double, std::map<uint32_t, double[3]> >::iterator it;
  double time_list[curr_obs.size()];
  double time_min = time2sec(curr_obs[0]->time);
  for (uint32_t j = 0; j < curr_obs.size(); ++j)
  {
    time_list[j] = 0.0;
    for (it = sat2cp.begin(); it != sat2cp.end(); ++it)
    {
      std::map<uint32_t, double[3]>::iterator it_sat;
      it_sat = it->second.find(curr_obs[j]->sat);
      if (it_sat != it->second.end())
      {
        time_list[j] = it->first;
        if (time_list[j] < time_min)
        {
          time_min = time_list[j];
        }
        break;
      }
    }
  }

  int best_sat = -1;
  double min_cp_std = 100000000;
  for (uint32_t j = 0; j < curr_obs.size(); ++j)
  {
    if (time_list[j] == time_min)
    {
      freq = L1_freq(curr_obs[j], &freq_idx);
      LOG_IF(FATAL, freq < 0) << "No L1 observation found."; 
      if (curr_obs[j]->cp_std[freq_idx] * 0.004 < min_cp_std && curr_obs[j]->cp[freq_idx] > 0)
      {
        min_cp_std = curr_obs[j]->cp_std[freq_idx] * 0.004;
        best_sat = j;
      }
    }
  }

  std::deque<uint32_t> pair_sat_copy;

  std::deque<double> meas_sats;
  std::deque<double> meas_cov_sats;
  std::deque<double> meas_time_sats;
  std::deque<double> meas_sats_copy;
  std::deque<double> meas_cov_sats_copy;
  // std::deque<double> meas_time_sats_copy;
  
  if (best_sat > -1)
  {
    std::deque<uint32_t> pair_sat;

    for (uint32_t j = 0; j < curr_obs.size(); ++j)
    {
      if (j != best_sat)
      {
        uint32_t sys = satsys(curr_obs[j]->sat, NULL);
        freq = L1_freq(curr_obs[j], &freq_idx);

        if (sys == satsys(curr_obs[best_sat]->sat, NULL) && curr_obs[j]->cp_std[freq_idx] * 0.004 < gnss_cp_std_threshold)
        {
          pair_sat.push_back(j);
        }
      }
    }

    for (uint32_t j = 0; j < curr_obs.size(); ++j)
    {
      std::map<uint32_t, double[3]>::iterator it_old_best, it_old;
      std::vector<double> meas;
      std::vector<double> meas_cov;
      std::vector<double> meas_time;
      if (j == pair_sat.front())
      {
        for (it = sat2cp.begin(); it != sat2cp.end(); ++it)
        {
          it_old = it->second.find(curr_obs[j]->sat);
          it_old_best = it->second.find(curr_obs[best_sat]->sat);
          if (it_old != it->second.end() && it_old_best != it->second.end())
          {
            meas.push_back(it_old_best->second[0] - it_old->second[0] - it_old_best->second[1] + it_old->second[1]);
            meas_cov.push_back(it_old_best->second[2] * it_old_best->second[2] + it_old->second[2] * it_old->second[2]);
            meas_time.push_back(it->first);
          }
        }
        if (meas.size() > 0)
        {
          double sum = std::accumulate(meas.begin(), meas.end(), 0.0);
          double mean = sum / meas.size(); // meas[meas.size()-1]; // 
          std::vector<double> empty_vec_d;
          meas.swap(empty_vec_d);
          meas_sats.push_back(mean);
          meas_cov_sats.push_back(meas_cov[meas_cov.size()-1]);
          meas_time_sats.push_back(meas_time[meas_time.size()-1]);
          meas_cov.swap(empty_vec_d);
          meas_time.swap(empty_vec_d);
          pair_sat_copy.push_back(pair_sat.front());
        }
        pair_sat.pop_front();
      }
    }
  }

  Eigen::MatrixXd jacobians_cur(curr_obs.size() * 2, jacobians.cols());
  Eigen::VectorXd residuals_cur(curr_obs.size() * 2);
  residuals_cur.setZero();
  jacobians_cur.setZero();
  std::map<uint32_t, double[3]> curr_cp_map;

  std::vector<double> meas_cp;
  std::vector<double> esti_cp;
  std::vector<double> cov_cp;
  std::vector<double> delta_t_cp;
  std::vector<Eigen::Vector3d> jaco_cp;
  double meas_cp_best, esti_cp_best, cov_cp_best;
  Eigen::Vector3d jaco_cp_best;

  for (uint32_t j = 0; j < curr_obs.size(); ++j)
  {
    const uint32_t sys = satsys(curr_obs[j]->sat, NULL);
    const uint32_t sys_idx = gnss_comm::sys2idx.at(sys);
    GnssPsrDoppMeas(curr_obs[j], curr_ephem[j]); //, latest_gnss_iono_params);
  
    double ion_delay = 0, tro_delay = 0;
    double azel[2] = {0, M_PI/2.0};
    if (P_ecef.norm() > 0)
    {
        sat_azel(P_ecef, sv_pos, azel);
        Eigen::Vector3d rcv_lla = ecef2geo(P_ecef);
        tro_delay = calculate_trop_delay(curr_obs[j]->time, rcv_lla, azel);
        ion_delay = calculate_ion_delay(curr_obs[j]->time, latest_gnss_iono_params, rcv_lla, azel); // rely on local pose
    }
    double sin_el = sin(azel[1]);
    double sin_el_2 = sin_el*sin_el;
    double pr_weight = sin_el_2 / pr_uura * relative_sqrt_info; // not requisite
    double dp_weight = sin_el_2 / dp_uura * relative_sqrt_info * PSR_TO_DOPP_RATIO; // not requisite

    Eigen::Vector3d rcv2sat_ecef = sv_pos - P_ecef;
    Eigen::Vector3d rcv2sat_unit = rcv2sat_ecef.normalized();

    freq = L1_freq(curr_obs[j], &freq_idx);
    
    const double wavelength = LIGHT_SPEED / freq;

    if (curr_obs[j]->cp_std[freq_idx] * 0.004 < gnss_cp_std_threshold)
    {
      if (curr_obs[j]->cp[freq_idx] > 0)
      {
        curr_cp_map[curr_obs[j]->sat][0] = curr_obs[j]->cp[freq_idx] * wavelength;
        curr_cp_map[curr_obs[j]->sat][1] = rcv2sat_ecef.norm();
        curr_cp_map[curr_obs[j]->sat][2] = curr_obs[j]->cp_std[freq_idx] * 0.004;
      }
    }

    if (j == best_sat)
    {
      meas_cp_best = curr_obs[j]->cp[freq_idx] * wavelength;
      esti_cp_best = rcv2sat_ecef.norm();
      jaco_cp_best = rcv2sat_unit;
      cov_cp_best  = curr_obs[j]->cp_std[freq_idx] * 0.004 * curr_obs[j]->cp_std[freq_idx] * 0.004;
    }

    if (j == pair_sat_copy.front())
    {
      if (curr_obs[j]->cp[freq_idx] > 0)
      {
      meas_cp.push_back(curr_obs[j]->cp[freq_idx] * wavelength);
      cov_cp.push_back(curr_obs[j]->cp_std[freq_idx] * curr_obs[j]->cp_std[freq_idx] * 0.004 * 0.004);
      delta_t_cp.push_back(time2sec(curr_obs[j]->time) - meas_time_sats.front());
      esti_cp.push_back(rcv2sat_ecef.norm());
      jaco_cp.push_back(rcv2sat_unit);
      meas_sats_copy.push_back(meas_sats.front());
      meas_cov_sats_copy.push_back(meas_cov_sats.front());
      }
      pair_sat_copy.pop_front();
      meas_sats.pop_front();
      meas_cov_sats.pop_front();
      meas_time_sats.pop_front();
    }

    const double psr_sagnac = EARTH_OMG_GPS*(sv_pos(0)*P_ecef(1)-sv_pos(1)*P_ecef(0))/LIGHT_SPEED;
    double psr_estimated = rcv2sat_ecef.norm() + psr_sagnac + rcv_dt[sys_idx] - svdt*LIGHT_SPEED + // why not multiply light_speed?  
                              ion_delay + tro_delay + tgd*LIGHT_SPEED;
    const double dopp_sagnac = EARTH_OMG_GPS/LIGHT_SPEED*(sv_vel(0)*P_ecef(1)+
              sv_pos(0)*V_ecef(1) - sv_vel(1)*P_ecef(0) - sv_pos(1)*V_ecef(0));
    double dopp_estimated = (sv_vel - V_ecef).dot(rcv2sat_unit) + rcv_ddt + dopp_sagnac - svddt*LIGHT_SPEED;
    // cout << "CHECK VALUE:" << ion_delay << ";" << tro_delay << ";" << svdt*LIGHT_SPEED << ";" << rcv_dt[sys_idx] << ";" << tgd*LIGHT_SPEED << ";" << psr_sagnac << endl;
    // cout << "CHECK VALUE dp:" << rcv_ddt << ";" << svddt*LIGHT_SPEED << ";" << dopp_sagnac << ";" << state.anc.transpose() << ";" << state.yaw_enu_local << endl;
    cout << "check weight:" << dp_weight << ";" << pr_weight << endl;
    Eigen::Matrix3d hat_T;
    hat_T << SKEW_SYM_MATRX(Tex_imu_r);
    if (!nolidar)
    {    
      residuals_cur[2*j] = -(psr_estimated - curr_obs[j]->psr[freq_idx]) * pr_weight;

      residuals_cur[2*j+1] = -(dopp_estimated + curr_obs[j]->dopp[freq_idx]*wavelength) * dp_weight; // / wavelength; // 
      // J_Pi

      jacobians_cur.block<1, 3>(2*j,3) = -rcv2sat_unit.transpose() * R_ecef_local * pr_weight; // why use normalized result?
      jacobians_cur.block<1, 3>(2*j,0) = rcv2sat_unit.transpose() * R_ecef_local * hat_T * pr_weight; // why use normalized result?
      jacobians_cur.block<1, 3>(2*j,30) = -rcv2sat_unit.transpose() * (Eye3d + R_ecef_enu_cur * hatP * R1TE3 * vecLon.transpose() - R_ecef_enu_cur * hatP * E1 * vecLat.transpose()) * pr_weight; // simplify

      const double norm3 = pow(rcv2sat_ecef.norm(), 3);
      const double norm2 = rcv2sat_ecef.squaredNorm();
      Eigen::Matrix3d unit2rcv_pos;
      for (size_t i = 0; i < 3; ++i)
      {
          for (size_t k = 0; k < 3; ++k)
          {
              if (i == k)
                  unit2rcv_pos(i, k) = (norm2-rcv2sat_ecef(i)*rcv2sat_ecef(i))/norm3;
              else
                  unit2rcv_pos(i, k) = (-rcv2sat_ecef(i)*rcv2sat_ecef(k))/norm3;
          }
      }
      unit2rcv_pos *= -1;
      jacobians_cur.block<1, 3>(2*j+1,3) = (sv_vel-V_ecef).transpose() * unit2rcv_pos * 
          R_ecef_local * dp_weight; //  / wavelength; //
      jacobians_cur.block<1, 3>(2*j+1,0) = -(sv_vel-V_ecef).transpose() * unit2rcv_pos * 
          R_ecef_local * hat_T * dp_weight; //  / wavelength; //
      jacobians_cur.block<1, 3>(2*j+1,30) = (sv_vel-V_ecef).transpose() * unit2rcv_pos * (Eye3d + R_ecef_enu_cur * hatP * R1TE3 * vecLon.transpose() - R_ecef_enu_cur * hatP * E1 * vecLat.transpose()) * dp_weight
                                      + rcv2sat_unit.transpose() * (-1.0) * (R_ecef_enu_cur * hatV * R1TE3 * vecLon.transpose() - R_ecef_enu_cur * hatV * E1 * vecLat.transpose()) * dp_weight;
    

      // J_Vi
      jacobians_cur.block<1, 3>(2*j+1,6) = rcv2sat_unit.transpose() * (-1.0) * 
          R_ecef_local * dp_weight; //  / wavelength; // 
      // J_rcv_dt
      jacobians_cur(2*j,25+sys_idx) = 1.0 * pr_weight;

      // J_rcv_ddt
      jacobians_cur(2*j+1,29) = 1.0 * dp_weight; //  / wavelength; //

      // J_yaw_diff
      Eigen::Matrix3d d_yaw;
      d_yaw << 0, -1, 0, 
                1, 0, 0, 
                0, 0, 0;
      jacobians_cur(2*j,24) = -rcv2sat_unit.dot( (R_ecef_local * d_yaw * (local_pos - anc_local))) * pr_weight;
      jacobians_cur(2*j+1,24) = -rcv2sat_unit.dot((R_ecef_local * d_yaw * local_vel)) * dp_weight + (sv_vel-V_ecef).dot( unit2rcv_pos * 
          R_ecef_local * d_yaw * (local_pos - anc_local)) * dp_weight; //  / wavelength; //
    }
    else
    {    
      residuals_cur[2*j] = -(psr_estimated - curr_obs[j]->psr[freq_idx]) * pr_weight;

      
      residuals_cur[2*j+1] = -(dopp_estimated + curr_obs[j]->dopp[freq_idx]*wavelength) * dp_weight; // / wavelength; // 
      // J_Pi

      jacobians_cur.block<1, 3>(2*j,3) = -rcv2sat_unit.transpose() * pr_weight; 
      jacobians_cur.block<1, 3>(2*j,0) = rcv2sat_unit.transpose() * hat_T * pr_weight; 

      const double norm3 = pow(rcv2sat_ecef.norm(), 3);
      const double norm2 = rcv2sat_ecef.squaredNorm();
      Eigen::Matrix3d unit2rcv_pos;
      for (size_t i = 0; i < 3; ++i)
      {
          for (size_t k = 0; k < 3; ++k)
          {
              if (i == k)
                  unit2rcv_pos(i, k) = (norm2-rcv2sat_ecef(i)*rcv2sat_ecef(i))/norm3;
              else
                  unit2rcv_pos(i, k) = (-rcv2sat_ecef(i)*rcv2sat_ecef(k))/norm3;
          }
      }
      unit2rcv_pos *= -1;
      jacobians_cur.block<1, 3>(2*j+1,3) = (sv_vel-V_ecef).transpose() * unit2rcv_pos * dp_weight; //  / wavelength; //
      jacobians_cur.block<1, 3>(2*j+1,0) = -(sv_vel-V_ecef).transpose() * unit2rcv_pos * hat_T * dp_weight; //  / wavelength; //

      // J_Vi
      jacobians_cur.block<1, 3>(2*j+1,6) = rcv2sat_unit.transpose() * (-1.0) * dp_weight; //  / wavelength; // 
      // J_rcv_dt
      jacobians_cur(2*j,25+sys_idx) = 1.0 * pr_weight;

      // J_rcv_ddt
      jacobians_cur(2*j+1,29) = 1.0 * dp_weight; //  / wavelength; //
    }
  }

  residuals.resize(curr_obs.size() * 2 + meas_sats_copy.size());
  jacobians.resize(curr_obs.size()* 2 + meas_sats_copy.size(), jacobians.cols());
  meas_size = curr_obs.size() * 2 + meas_sats_copy.size();
  residuals.setZero();
  jacobians.setZero();
  for (int j = 0; j < curr_obs.size() * 2; j++)
  {
    residuals(j) = residuals_cur(j);
    jacobians.block<1,33>(j,0) = jacobians_cur.block<1,33>(j,0);
  }

  double cp_weight = 0.10;
  for (uint32_t j = 0; j < meas_sats_copy.size(); ++j)
  {
    // cp_weight = 0.01 / (cov_cp_best + cov_cp[j] + meas_cov_sats_copy[j]) / (1+2*delta_t_cp[j]+delta_t_cp[j]*delta_t_cp[j]);
    // cout << "check cp std:" << cp_weight << endl;
    // cout << "check cp:" << meas_sats_copy[j] << ";" << meas_cp_best << ";" << meas_cp[j] << ";" << esti_cp_best << ";" << esti_cp[j] << endl;
    residuals[2*curr_obs.size() + j] = -((meas_sats_copy[j] - (meas_cp_best - meas_cp[j]) + (esti_cp_best - esti_cp[j]))) * cp_weight;
    Eigen::Matrix3d hat_T;
    hat_T << SKEW_SYM_MATRX(Tex_imu_r);
    if (!nolidar)
    {
      jacobians.block<1, 3>(2*curr_obs.size() + j,3) = (-jaco_cp_best.transpose() + jaco_cp[j].transpose()) * R_ecef_local * cp_weight; // why use normalized result?
      jacobians.block<1, 3>(2*curr_obs.size() + j,0) = (jaco_cp_best.transpose() - jaco_cp[j].transpose()) * R_ecef_local * hat_T * cp_weight; // why use normalized result?
      
      jacobians.block<1, 3>(2*curr_obs.size() + j,30) = (-jaco_cp_best.transpose() + jaco_cp[j].transpose()) * (Eye3d + R_ecef_enu_cur * hatP * R1TE3 * vecLon.transpose() - R_ecef_enu_cur * hatP * E1 * vecLat.transpose()) * cp_weight; // simplify

      Eigen::Matrix3d d_yaw;
      d_yaw << 0, -1, 0, 
                1, 0, 0, 
                0, 0, 0;
      jacobians(2*curr_obs.size() + j,24) = (-jaco_cp_best + jaco_cp[j]).dot( (R_ecef_local * d_yaw * (local_pos - anc_local))) * cp_weight;
    }
    else
    {
      jacobians.block<1, 3>(2*curr_obs.size() + j,3) = (-jaco_cp_best.transpose() + jaco_cp[j].transpose()) * cp_weight; // why use normalized result?
      jacobians.block<1, 3>(2*curr_obs.size() + j,0) = (jaco_cp_best.transpose() - jaco_cp[j].transpose()) * hat_T * cp_weight; // why use normalized result?
    }
  }
  sat2cp[time2sec(curr_obs[0]->time)] = curr_cp_map;
  std::map<double, std::map<uint32_t, double[3]>>::iterator it_old;
  if (!sat2cp.empty())
  {
    it_old = sat2cp.begin();
    if (time2sec(curr_obs[0]->time) - it_old->first > gnss_cp_time_threshold)
    {
      std::map<uint32_t, double[3]>().swap(it_old->second);
      size_t del_size = sat2cp.erase(it_old->first);
    }
  }
  return true;
}

// bool GNSSProcess::Evaluate(StatesGroupwithGNSS2 &state, Eigen::VectorXd &residuals, Eigen::MatrixXd &jacobians, int &meas_size)
// {
//   if (gnss_meas_buf[1].empty() && gnss_meas_buf[0].empty()) 
//   {
//     residuals.setZero();
//     jacobians.setZero();
//     return false;
//   }
//   double rcv_dt[4];
//   rcv_dt[0] = state.dt_g;
//   rcv_dt[1] = state.dt_r;
//   rcv_dt[2] = state.dt_e;
//   rcv_dt[3] = state.dt_c;
//   double rcv_ddt = state.ddt;
//   double yaw_diff = state.yaw_enu_local;
//   Eigen::Vector3d ref_ecef;
//   if (nolidar)
//   {
//    ref_ecef = state.pos_end;
//   }
//   else
//   {
//    ref_ecef = state.anc;
//   }
//   double s1, s2, e2, a, ep, p, h, lat, ds1dx, ds2dx, ds1dy, ds2dy, ds1dz, ds2dz, sins, coss;
//   e2 = EARTH_ECCE_2;
//   a = EARTH_SEMI_MAJOR; // _glo?
//   ep = ref_ecef(0)*ref_ecef(0) + ref_ecef(1) * ref_ecef(1);
//   p = a*a*(1-e2);
//   h = ref_ecef(2)*ref_ecef(2)*a*a;
//   s1 = ref_ecef(2) + e2/(1-e2) * sqrt(p) * pow(ref_ecef(2)*a/sqrt(h+ep*p),3);
//   s2 = sqrt(ep) - a * e2 * pow((ep*p)/(h+ep*p),1.5);
//   lat = atan(s1/s2);
//   sins = -s1/(s1*s1+s2*s2);
//   coss = s2/(s1*s1+s2*s2);
//   Eigen::Vector3d R1TE3, E1, vecP, vecV, vecLon, vecLat;
//   R1TE3 << 0.0, -sin(lat-M_PI/2), cos(lat-M_PI/2);
//   E1 << 1.0, 0.0, 0.0;

//   ds1dx = e2/(1-e2) * sqrt(p) * a * ref_ecef(2) * h * (-3) * p * ref_ecef(0) / pow(h+ep*p,2.5);
//   ds1dy = e2/(1-e2) * sqrt(p) * a * ref_ecef(2) * h * (-3) * p * ref_ecef(1) / pow(h+ep*p,2.5);
//   ds1dz = 1 + e2/(1-e2) * sqrt(p) * 3 *sqrt(h/(h+ep*p)) * a * a * ref_ecef(2) * ep * p / pow(h+ep*p,2);

//   ds2dx = ref_ecef(0) / sqrt(ep) - a * e2 * pow(p,1.5) * 3 * sqrt(ep)*ref_ecef(0)*h/pow(h+ep*p,2.5);
//   ds2dy = ref_ecef(1) / sqrt(ep) - a * e2 * pow(p,1.5) * 3 * sqrt(ep)*ref_ecef(1)*h/pow(h+ep*p,2.5);
//   ds2dz = a*e2*3 * pow(p,1.5) * a * a * ref_ecef(2) * pow(ep, 1.5) / pow(h+ep*p, 2.5);

//   vecLon << -ref_ecef(1)/ep, ref_ecef(0)/ep, 0.0;
//   vecLat << coss * ds1dx + sins * ds2dx, coss * ds1dy + sins * ds2dy, coss * ds1dz + sins * ds2dz;

//   const Eigen::Vector3d local_pos = state.rot_end * Tex_imu_r + state.pos_end;
//   const Eigen::Vector3d local_vel = state.vel_end; // check if it is right

//   double sin_yaw_diff = std::sin(yaw_diff);
//   double cos_yaw_diff = std::cos(yaw_diff);
//   Eigen::Matrix3d R_enu_local;
//   R_enu_local << cos_yaw_diff, -sin_yaw_diff, 0,
//                   sin_yaw_diff,  cos_yaw_diff, 0,
//                   0           ,  0           , 1;
//   Eigen::Matrix3d R_ecef_enu_cur = ecef2rotation(ref_ecef); // provide anchor value
//   Eigen::Matrix3d R_ecef_local = R_ecef_enu_cur * R_enu_local;

//   Eigen::Vector3d P_ecef;
//   if (nolidar)
//   {
//     P_ecef = local_pos;
//   }
//   else
//   {
//    P_ecef = R_ecef_local * (local_pos - anc_local) + ref_ecef;
//   }
//   Eigen::Vector3d V_ecef;
//   if (nolidar)
//   {
//    V_ecef = local_vel;
//   }
//   else
//   {
//     V_ecef = R_ecef_local * local_vel;
//   }

//   vecP = R_enu_local * (local_pos - anc_local);
//   vecV = R_enu_local * local_vel;
//   Eigen::Matrix3d hatP, hatV;
//   hatP << 0.0, -vecP(2), vecP(1),
//           vecP(2), 0.0, -vecP(0),
//           -vecP(1), vecP(0), 0.0;
//   hatV << 0.0, -vecV(2), vecV(1),
//           vecV(2), 0.0, -vecV(0),
//           -vecV(1), vecV(0), 0.0;

//   const std::vector<ObsPtr> &curr_obs = gnss_meas_buf[0];
//   const std::vector<EphemBasePtr> &curr_ephem = gnss_ephem_buf[0];
//   const std::vector<ObsPtr> &curr_cp_obs = gnss_meas_buf[1];
//   const std::vector<EphemBasePtr> &curr_cp_ephem = gnss_ephem_buf[1];

//   std::map<double, std::map<uint32_t, double[2]> >::iterator it;
//   double time_list[curr_cp_obs.size()];
//   double time_min = time2sec(curr_cp_obs[0]->time);
//   for (uint32_t j = 0; j < curr_cp_obs.size(); ++j)
//   {
//     time_list[j] = 0.0;
//     for (it = sat2cp.begin(); it != sat2cp.end(); ++it)
//     {
//       std::map<uint32_t, double[2]>::iterator it_sat;
//       it_sat = it->second.find(curr_cp_obs[j]->sat);
//       if (it_sat != it->second.end())
//       {
//         time_list[j] = it->first;
//         if (time_list[j] < time_min)
//         {
//           time_min = time_list[j];
//         }
//         break;
//       }
//     }
//   }

//   int best_sat = -1;
//   double min_cp_std = 100000000;
//   for (uint32_t j = 0; j < curr_cp_obs.size(); ++j)
//   {
//     if (time_list[j] == time_min)
//     {
//       freq = L1_freq(curr_cp_obs[j], &freq_idx);
//       LOG_IF(FATAL, freq < 0) << "No L1 observation found."; 
//       if (curr_cp_obs[j]->cp_std[freq_idx] * 0.004 < min_cp_std && curr_cp_obs[j]->cp[freq_idx] > 0)
//       {
//         min_cp_std = curr_cp_obs[j]->cp_std[freq_idx] * 0.004;
//         best_sat = j;
//       }
//     }
//   }

//   std::deque<uint32_t> pair_sat_copy;

//   std::deque<double> meas_sats;
//   std::deque<double> meas_sats_copy;
  
//   if (best_sat > -1)
//   {
//     std::deque<uint32_t> pair_sat;

//     for (uint32_t j = 0; j < curr_cp_obs.size(); ++j)
//     {
//       if (j != best_sat)
//       {
//         uint32_t sys = satsys(curr_cp_obs[j]->sat, NULL);
//         freq = L1_freq(curr_cp_obs[j], &freq_idx);

//         if (sys == satsys(curr_cp_obs[best_sat]->sat, NULL) && curr_cp_obs[j]->cp_std[freq_idx] * 0.004 < gnss_cp_std_threshold)
//         {
//           pair_sat.push_back(j);
//         }
//       }
//     }

//     for (uint32_t j = 0; j < curr_cp_obs.size(); ++j)
//     {
//       std::map<uint32_t, double[2]>::iterator it_old_best, it_old;
//       std::vector<double> meas;
//       if (j == pair_sat.front())
//       {
//         for (it = sat2cp.begin(); it != sat2cp.end(); ++it)
//         {
//           it_old = it->second.find(curr_cp_obs[j]->sat);
//           it_old_best = it->second.find(curr_cp_obs[best_sat]->sat);
//           if (it_old != it->second.end() && it_old_best != it->second.end())
//           {
//             meas.push_back(it_old_best->second[0] - it_old->second[0] - it_old_best->second[1] + it_old->second[1]);
//           }
//         }
//         if (meas.size() > 0)
//         {
//           double sum = std::accumulate(meas.begin(), meas.end(), 0.0);
//           double mean = sum / meas.size();
//           std::vector<double> empty_vec_d;
//           meas.swap(empty_vec_d);
//           meas_sats.push_back(mean);
//           pair_sat_copy.push_back(pair_sat.front());
//         }
//         pair_sat.pop_front();
//       }
//     }
//   }

//   std::map<uint32_t, double[2]> curr_cp_map;

//   std::vector<double> meas_cp;
//   std::vector<double> esti_cp;
//   std::vector<Eigen::Vector3d> jaco_cp;
//   double meas_cp_best, esti_cp_best;
//   Eigen::Vector3d jaco_cp_best;

//   for (uint32_t j = 0; j < curr_cp_obs.size(); ++j)
//   {
//     SvPosCals(curr_cp_obs[j], curr_cp_ephem[j]); //, latest_gnss_iono_params);
//     Eigen::Vector3d rcv2sat_ecef = sv_pos - P_ecef;
//     Eigen::Vector3d rcv2sat_unit = rcv2sat_ecef.normalized();

//     freq = L1_freq(curr_cp_obs[j], &freq_idx);
    
//     const double wavelength = LIGHT_SPEED / freq;

//     if (curr_cp_obs[j]->cp_std[freq_idx] * 0.004 < gnss_cp_std_threshold)
//     {
//       if (curr_cp_obs[j]->cp[freq_idx] > 0)
//       {
//         curr_cp_map[curr_cp_obs[j]->sat][0] = curr_cp_obs[j]->cp[freq_idx] * wavelength;
//         curr_cp_map[curr_cp_obs[j]->sat][1] = rcv2sat_ecef.norm();
//       }
//     }

//     if (j == best_sat)
//     {
//       meas_cp_best = curr_cp_obs[j]->cp[freq_idx] * wavelength;
//       esti_cp_best = rcv2sat_ecef.norm();
//       jaco_cp_best = rcv2sat_unit;
//     }

//     if (j == pair_sat_copy.front())
//     {
//       if (curr_cp_obs[j]->cp[freq_idx] > 0)
//       {
//       meas_cp.push_back(curr_cp_obs[j]->cp[freq_idx] * wavelength);
//       esti_cp.push_back(rcv2sat_ecef.norm());
//       jaco_cp.push_back(rcv2sat_unit);
//       meas_sats_copy.push_back(meas_sats.front());
//       }
//       pair_sat_copy.pop_front();
//       meas_sats.pop_front();
//     }
//   }

//   // Eigen::MatrixXd jacobians_cur(curr_obs.size() * 2, jacobians.cols());
//   // Eigen::VectorXd residuals_cur(curr_obs.size() * 2);
//   // residuals_cur.setZero();
//   // jacobians_cur.setZero();

//   meas_size = curr_obs.size() * 2 + meas_sats_copy.size();
//   if (meas_size < 1) return false;
//   residuals.resize(meas_size);
//   jacobians.resize(meas_size, jacobians.cols());
//   residuals.setZero();
//   jacobians.setZero();
//   // std::map<uint32_t, double[2]> curr_cp_map;

//   // std::vector<double> meas_cp;
//   // std::vector<double> esti_cp;
//   // std::vector<Eigen::Vector3d> jaco_cp;
//   // double meas_cp_best, esti_cp_best;
//   // Eigen::Vector3d jaco_cp_best;

//   for (uint32_t j = 0; j < curr_obs.size(); ++j)
//   {
//     const uint32_t sys = satsys(curr_obs[j]->sat, NULL);
//     const uint32_t sys_idx = gnss_comm::sys2idx.at(sys);
//     GnssPsrDoppMeas(curr_obs[j], curr_ephem[j]); //, latest_gnss_iono_params);
  
//     double ion_delay = 0, tro_delay = 0;
//     double azel[2] = {0, M_PI/2.0};
//     if (P_ecef.norm() > 0)
//     {
//         sat_azel(P_ecef, sv_pos, azel);
//         Eigen::Vector3d rcv_lla = ecef2geo(P_ecef);
//         tro_delay = calculate_trop_delay(curr_obs[j]->time, rcv_lla, azel);
//         ion_delay = calculate_ion_delay(curr_obs[j]->time, latest_gnss_iono_params, rcv_lla, azel); // rely on local pose
//     }
//     double sin_el = sin(azel[1]);
//     double sin_el_2 = sin_el*sin_el;
//     double pr_weight = sin_el_2 / pr_uura * relative_sqrt_info; // not requisite
//     double dp_weight = sin_el_2 / dp_uura * relative_sqrt_info * PSR_TO_DOPP_RATIO; // not requisite

//     Eigen::Vector3d rcv2sat_ecef = sv_pos - P_ecef;
//     Eigen::Vector3d rcv2sat_unit = rcv2sat_ecef.normalized();

//     freq = L1_freq(curr_obs[j], &freq_idx);
    
//     const double wavelength = LIGHT_SPEED / freq;

//     // if (curr_obs[j]->cp_std[freq_idx] * 0.004 < gnss_cp_std_threshold)
//     // {
//     //   if (curr_obs[j]->cp[freq_idx] > 0)
//     //   {
//     //     curr_cp_map[curr_obs[j]->sat][0] = curr_obs[j]->cp[freq_idx] * wavelength;
//     //     curr_cp_map[curr_obs[j]->sat][1] = rcv2sat_ecef.norm();
//     //   }
//     // }

//     // if (j == best_sat)
//     // {
//     //   meas_cp_best = curr_obs[j]->cp[freq_idx] * wavelength;
//     //   // cout << "line 1346:" << curr_obs[j]->cp[freq_idx] << ";" << wavelength << endl;
//     //   esti_cp_best = rcv2sat_ecef.norm();
//     //   jaco_cp_best = rcv2sat_unit;
//     // }

//     // if (j == pair_sat_copy.front())
//     // {
//     //   if (curr_obs[j]->cp[freq_idx] > 0)
//     //   {
//     //   meas_cp.push_back(curr_obs[j]->cp[freq_idx] * wavelength);
//     //   // cout << "line 1348:" << curr_obs[j]->cp[freq_idx] << ";" << wavelength << endl;
//     //   esti_cp.push_back(rcv2sat_ecef.norm());
//     //   jaco_cp.push_back(rcv2sat_unit);
//     //   meas_sats_copy.push_back(meas_sats.front());
//     //   }
//     //   pair_sat_copy.pop_front();
//     //   meas_sats.pop_front();
//     // }

//     const double psr_sagnac = EARTH_OMG_GPS*(sv_pos(0)*P_ecef(1)-sv_pos(1)*P_ecef(0))/LIGHT_SPEED;
//     double psr_estimated = rcv2sat_ecef.norm() + psr_sagnac + rcv_dt[sys_idx] - svdt*LIGHT_SPEED + // why not multiply light_speed?  
//                               ion_delay + tro_delay + tgd*LIGHT_SPEED;
//     const double dopp_sagnac = EARTH_OMG_GPS/LIGHT_SPEED*(sv_vel(0)*P_ecef(1)+
//               sv_pos(0)*V_ecef(1) - sv_vel(1)*P_ecef(0) - sv_pos(1)*V_ecef(0));
//     double dopp_estimated = (sv_vel - V_ecef).dot(rcv2sat_unit) + rcv_ddt + dopp_sagnac - svddt*LIGHT_SPEED;
//     // cout << "CHECK VALUE:" << ion_delay << ";" << tro_delay << ";" << svdt*LIGHT_SPEED << ";" << rcv_dt[sys_idx] << ";" << tgd*LIGHT_SPEED << ";" << psr_sagnac << endl;
//     // cout << "CHECK VALUE dp:" << rcv_ddt << ";" << svddt*LIGHT_SPEED << ";" << dopp_sagnac << ";" << state.anc.transpose() << ";" << state.yaw_enu_local << endl;
      // Eigen::Matrix3d hat_T;
//       hat_T << SKEW_SYM_MATRX(Tex_imu_r);
//     if (!nolidar)
//     {    
//       residuals[2*j] = -(psr_estimated - curr_obs[j]->psr[freq_idx]) * pr_weight;

//       residuals[2*j+1] = -(dopp_estimated + curr_obs[j]->dopp[freq_idx]*wavelength) * dp_weight; // / wavelength; // 
//       // J_Pi

//       jacobians.block<1, 3>(2*j,3) = -rcv2sat_unit.transpose() * R_ecef_local * pr_weight; // why use normalized result?
//       jacobians.block<1, 3>(2*j,30) = -rcv2sat_unit.transpose() * (Eye3d + R_ecef_enu_cur * hatP * R1TE3 * vecLon.transpose() - R_ecef_enu_cur * hatP * E1 * vecLat.transpose()) * pr_weight; // simplify
//       
//       jacobians.block<1, 3>(2*j,0) = rcv2sat_unit.transpose() * R_ecef_local * state.rot_end * hat_T * pr_weight;
//       const double norm3 = pow(rcv2sat_ecef.norm(), 3);
//       const double norm2 = rcv2sat_ecef.squaredNorm();
//       Eigen::Matrix3d unit2rcv_pos;
//       for (size_t i = 0; i < 3; ++i)
//       {
//           for (size_t k = 0; k < 3; ++k)
//           {
//               if (i == k)
//                   unit2rcv_pos(i, k) = (norm2-rcv2sat_ecef(i)*rcv2sat_ecef(i))/norm3;
//               else
//                   unit2rcv_pos(i, k) = (-rcv2sat_ecef(i)*rcv2sat_ecef(k))/norm3;
//           }
//       }
//       unit2rcv_pos *= -1;
//       jacobians.block<1, 3>(2*j+1,3) = (sv_vel-V_ecef).transpose() * unit2rcv_pos * 
//           R_ecef_local * dp_weight; //  / wavelength; //
//       jacobians.block<1, 3>(2*j+1,0) = -(sv_vel-V_ecef).transpose() * unit2rcv_pos * 
//           R_ecef_local * state.rot_end * hat_T * dp_weight; //  / wavelength; //
//       jacobians.block<1, 3>(2*j+1,30) = (sv_vel-V_ecef).transpose() * unit2rcv_pos * (Eye3d + R_ecef_enu_cur * hatP * R1TE3 * vecLon.transpose() - R_ecef_enu_cur * hatP * E1 * vecLat.transpose()) * dp_weight
//                                       + rcv2sat_unit.transpose() * (-1.0) * (R_ecef_enu_cur * hatV * R1TE3 * vecLon.transpose() - R_ecef_enu_cur * hatV * E1 * vecLat.transpose()) * dp_weight;
    

//       // J_Vi
//       jacobians.block<1, 3>(2*j+1,6) = rcv2sat_unit.transpose() * (-1.0) * 
//           R_ecef_local * dp_weight; //  / wavelength; // 
//       // J_rcv_dt
//       jacobians(2*j,25+sys_idx) = 1.0 * pr_weight;

//       // J_rcv_ddt
//       jacobians(2*j+1,29) = 1.0 * dp_weight; //  / wavelength; //

//       // J_yaw_diff
//       Eigen::Matrix3d d_yaw;
//       d_yaw << 0, -1, 0, 
//                 1, 0, 0, 
//                 0, 0, 0;
//       jacobians(2*j,24) = -rcv2sat_unit.dot( (R_ecef_local * d_yaw * (local_pos - anc_local))) * pr_weight;
//       jacobians(2*j+1,24) = -rcv2sat_unit.dot((R_ecef_local * d_yaw * local_vel)) * dp_weight + (sv_vel-V_ecef).dot( unit2rcv_pos * 
//           R_ecef_local * d_yaw * (local_pos - anc_local)) * dp_weight; //  / wavelength; //
//     }
//     else
//     {    
//       residuals[2*j] = -(psr_estimated - curr_obs[j]->psr[freq_idx]) * pr_weight;

      
//       residuals[2*j+1] = -(dopp_estimated + curr_obs[j]->dopp[freq_idx]*wavelength) * dp_weight; // / wavelength; // 
//       // J_Pi

//       jacobians.block<1, 3>(2*j,3) = -rcv2sat_unit.transpose() * pr_weight; 
//       jacobians.block<1, 3>(2*j,0) = rcv2sat_unit.transpose() * hat_T * pr_weight; 

//       const double norm3 = pow(rcv2sat_ecef.norm(), 3);
//       const double norm2 = rcv2sat_ecef.squaredNorm();
//       Eigen::Matrix3d unit2rcv_pos;
//       for (size_t i = 0; i < 3; ++i)
//       {
//           for (size_t k = 0; k < 3; ++k)
//           {
//               if (i == k)
//                   unit2rcv_pos(i, k) = (norm2-rcv2sat_ecef(i)*rcv2sat_ecef(i))/norm3;
//               else
//                   unit2rcv_pos(i, k) = (-rcv2sat_ecef(i)*rcv2sat_ecef(k))/norm3;
//           }
//       }
//       unit2rcv_pos *= -1;
//       jacobians.block<1, 3>(2*j+1,3) = (sv_vel-V_ecef).transpose() * unit2rcv_pos * dp_weight; //  / wavelength; //
//       jacobians.block<1, 3>(2*j+1,0) = -(sv_vel-V_ecef).transpose() * unit2rcv_pos * hat_T * dp_weight; //  / wavelength; //

//       // J_Vi
//       jacobians.block<1, 3>(2*j+1,6) = rcv2sat_unit.transpose() * (-1.0) * dp_weight; //  / wavelength; // 
//       // J_rcv_dt
//       jacobians(2*j,25+sys_idx) = 1.0 * pr_weight;

//       // J_rcv_ddt
//       jacobians(2*j+1,29) = 1.0 * dp_weight; //  / wavelength; //
//     }
//   }

//   // residuals.resize(curr_obs.size() * 2 + meas_sats_copy.size());
//   // jacobians.resize(curr_obs.size()* 2 + meas_sats_copy.size(), jacobians.cols());
//   // meas_size = curr_obs.size() * 2 + meas_sats_copy.size();
//   // residuals.setZero();
//   // jacobians.setZero();
//   // for (int j = 0; j < curr_cp_obs.size() * 2; j++)
//   // {
//   //   residuals(j) = residuals_cur(j);
//   //   jacobians.block<1,33>(j,0) = jacobians_cur.block<1,33>(j,0);
//   // }

//   for (uint32_t j = 0; j < meas_sats_copy.size(); ++j)
//   {
//     // cout << "check cp:" << meas_sats_copy[j] << ";" << meas_cp_best << ";" << meas_cp[j] << ";" << esti_cp_best << ";" << esti_cp[j] << endl;
//     residuals[2*curr_obs.size() + j] = meas_sats_copy[j] - (meas_cp_best - meas_cp[j]) + (esti_cp_best - esti_cp[j]);
      //Eigen::Matrix3d hat_T;
//       hat_T << SKEW_SYM_MATRX(Tex_imu_r);
//     if (!nolidar)
//     {
//       jacobians.block<1, 3>(2*curr_obs.size() + j,3) = -(-jaco_cp_best.transpose() * R_ecef_local + jaco_cp[j].transpose() * R_ecef_local); // why use normalized result?
//       
//       jacobians.block<1, 3>(2*curr_obs.size() + j,0) = (-jaco_cp_best.transpose() * R_ecef_local * state.rot_end * hat_T + jaco_cp[j].transpose() * R_ecef_local * state.rot_end * hat_T); // why use normalized result?
//       jacobians.block<1, 3>(2*curr_obs.size() + j,30) = -(-jaco_cp_best.transpose() + jaco_cp[j].transpose()) * (Eye3d + R_ecef_enu_cur * hatP * R1TE3 * vecLon.transpose() - R_ecef_enu_cur * hatP * E1 * vecLat.transpose()); // simplify

//       Eigen::Matrix3d d_yaw;
//       d_yaw << 0, -1, 0, 
//                 1, 0, 0, 
//                 0, 0, 0;
//       jacobians(2*curr_obs.size() + j,24) = jaco_cp_best.dot( (R_ecef_local * d_yaw * (local_pos - anc_local))) - jaco_cp[j].dot( (R_ecef_local * d_yaw * (local_pos - anc_local)));
//     }
//     else
//     {
//       jacobians.block<1, 3>(2*curr_obs.size() + j,3) = -(-jaco_cp_best.transpose() + jaco_cp[j].transpose()); // why use normalized result?
//       jacobians.block<1, 3>(2*curr_obs.size() + j,0) = (-jaco_cp_best.transpose() + jaco_cp[j].transpose()) * hat_T; // why use normalized result?
//     }
//   }
//   sat2cp[time2sec(curr_cp_obs[0]->time)] = curr_cp_map;
//   if (!sat2cp.empty())
//   {
//     std::map<double, std::map<uint32_t, double[2]>>::iterator it_old;
//     it_old = sat2cp.begin();
//     while (it_old != sat2cp.end())
//     {
//     if (time2sec(curr_cp_obs[0]->time) - it_old->first > gnss_cp_time_threshold)
//     {
//       std::map<uint32_t, double[2]>().swap(it_old->second);
//       it_old = sat2cp.erase(it_old);
//     }
//     else
//     {
//       break;
//     }
//     // if (sat2cp.empty()) break;
//     // it_old = sat2cp.begin();
//     }
//   }
//   // if (meas_size < 1)
//   // {
//   //   return false;
//   // }
//   return true;
// }

double GNSSProcess::str2double(const std::string &num_str)
{
    size_t D_pos = num_str.find("D");
    std::string tmp_str = num_str;
    if (D_pos != std::string::npos)
        tmp_str = tmp_str.replace(D_pos, 1, "e");
    return std::stod(tmp_str);
}

EphemPtr GNSSProcess::rinex_line2ephem(const std::vector<std::string> &ephem_lines)
{
    LOG_IF(FATAL, ephem_lines.size() != 8) << "Ephemeris record should contain 8 lines";
    uint32_t sat_sys = SYS_NONE;
    if      (ephem_lines[0].at(0) == 'G')    sat_sys = SYS_GPS;
    else if (ephem_lines[0].at(0) == 'C')    sat_sys = SYS_BDS;
    else if (ephem_lines[0].at(0) == 'E')    sat_sys = SYS_GAL;
    LOG_IF(FATAL, sat_sys == SYS_NONE) << "Satellite system is not supported: " << ephem_lines[0].at(0);

    EphemPtr ephem(new Ephem());
    uint32_t prn = static_cast<uint32_t>(std::stoi(ephem_lines[0].substr(1, 2)));
    ephem->sat = sat_no(sat_sys, prn);
    double epoch[6];
    epoch[0] = static_cast<double>(std::stoi(ephem_lines[0].substr(4, 4)));
    epoch[1] = static_cast<double>(std::stoi(ephem_lines[0].substr(9, 2)));
    epoch[2] = static_cast<double>(std::stoi(ephem_lines[0].substr(12, 2)));
    epoch[3] = static_cast<double>(std::stoi(ephem_lines[0].substr(15, 2)));
    epoch[4] = static_cast<double>(std::stoi(ephem_lines[0].substr(18, 2)));
    epoch[5] = static_cast<double>(std::stoi(ephem_lines[0].substr(21, 2)));
    ephem->toc = epoch2time(epoch);
    if (sat_sys == SYS_BDS)     ephem->toc.time += 14;     // BDS-GPS time correction
    ephem->af0 = str2double(ephem_lines[0].substr(23, 19));
    ephem->af1 = str2double(ephem_lines[0].substr(42, 19));
    ephem->af2 = str2double(ephem_lines[0].substr(61, 19));

    // the second line
    if (sat_sys == SYS_GPS)
        ephem->iode  = str2double(ephem_lines[1].substr(4, 19));
    ephem->crs       = str2double(ephem_lines[1].substr(23, 19));
    ephem->delta_n   = str2double(ephem_lines[1].substr(42, 19));
    ephem->M0        = str2double(ephem_lines[1].substr(61, 19));

    // the third line
    ephem->cuc = str2double(ephem_lines[2].substr(4, 19));
    ephem->e = str2double(ephem_lines[2].substr(23, 19));
    ephem->cus = str2double(ephem_lines[2].substr(42, 19));
    double sqrt_A = str2double(ephem_lines[2].substr(61, 19));
    ephem->A = sqrt_A * sqrt_A;

    // the forth line
    ephem->toe_tow = str2double(ephem_lines[3].substr(4, 19));
    ephem->cic = str2double(ephem_lines[3].substr(23, 19));
    ephem->OMG0 = str2double(ephem_lines[3].substr(42, 19));
    ephem->cis = str2double(ephem_lines[3].substr(61, 19));

    // the fifth line
    ephem->i0 = str2double(ephem_lines[4].substr(4, 19));
    ephem->crc = str2double(ephem_lines[4].substr(23, 19));
    ephem->omg = str2double(ephem_lines[4].substr(42, 19));
    ephem->OMG_dot = str2double(ephem_lines[4].substr(61, 19));

    // the sixth line
    ephem->i_dot = str2double(ephem_lines[5].substr(4, 19));
    if  (sat_sys == SYS_GAL)
    {
        uint32_t ephe_source = static_cast<uint32_t>(str2double(ephem_lines[5].substr(23, 19)));
        if (!(ephe_source & 0x01))  
        {
            // LOG(ERROR) << "not contain I/NAV E1-b info, skip this ephemeris";
            return ephem;   // only parse I/NAV E1-b ephemeris
        }
    }
    ephem->week = static_cast<uint32_t>(str2double(ephem_lines[5].substr(42, 19)));
    if (sat_sys == SYS_GPS || sat_sys == SYS_GAL)     ephem->toe = gpst2time(ephem->week, ephem->toe_tow);
    else if (sat_sys == SYS_BDS)                      ephem->toe = bdt2time(ephem->week, ephem->toe_tow+14);
    // if (sat_sys == SYS_GAL)     ephem->toe = gst2time(ephem->week, ephem->toe_tow);

    // the seventh line
    ephem->ura = str2double(ephem_lines[6].substr(4, 19));
    ephem->health = static_cast<uint32_t>(str2double(ephem_lines[6].substr(23, 19)));
    ephem->tgd[0] = str2double(ephem_lines[6].substr(42, 19));
    if (sat_sys == SYS_BDS || sat_sys == SYS_GAL)
        ephem->tgd[1] = str2double(ephem_lines[6].substr(61, 19));
    if (sat_sys == SYS_GPS)     ephem->iodc = str2double(ephem_lines[6].substr(61, 19));

    // the eighth line
    double ttr_tow = str2double(ephem_lines[7].substr(4, 19));
    // GAL week = GST week + 1024 + rollover, already align with GPS week!!!
    if      (sat_sys == SYS_GPS || sat_sys == SYS_GAL)   ephem->ttr = gpst2time(ephem->week, ttr_tow);
    else if (sat_sys == SYS_BDS)   ephem->ttr = bdt2time(ephem->week, ttr_tow);

    // convert time system to parameter GPST
    if (sat_sys == SYS_BDS)
    {
        uint32_t week = 0;
        ephem->toe_tow = time2gpst(ephem->toe, &week);
        ephem->week = week;
    }

    return ephem;
}

GloEphemPtr GNSSProcess::rinex_line2glo_ephem(const std::vector<std::string> &ephem_lines, const uint32_t gpst_leap_seconds)
{
    LOG_IF(FATAL, ephem_lines.size() != 4) << "GLO ephemeris record should contain 8 lines";
    LOG_IF(FATAL, ephem_lines[0].at(0) != 'R') << "Not a valid GLO ephemeris record";
    GloEphemPtr glo_ephem(new GloEphem());

    uint32_t prn = static_cast<uint32_t>(std::stoi(ephem_lines[0].substr(1, 2)));
    glo_ephem->sat = sat_no(SYS_GLO, prn);
    double epoch[6];
    epoch[0] = static_cast<double>(std::stoi(ephem_lines[0].substr(4, 4)));
    epoch[1] = static_cast<double>(std::stoi(ephem_lines[0].substr(9, 2)));
    epoch[2] = static_cast<double>(std::stoi(ephem_lines[0].substr(12, 2)));
    epoch[3] = static_cast<double>(std::stoi(ephem_lines[0].substr(15, 2)));
    epoch[4] = static_cast<double>(std::stoi(ephem_lines[0].substr(18, 2)));
    epoch[5] = static_cast<double>(std::stoi(ephem_lines[0].substr(21, 2)));
    glo_ephem->toe = epoch2time(epoch);
    glo_ephem->toe.time += gpst_leap_seconds;
    glo_ephem->tau_n = -1.0 * str2double(ephem_lines[0].substr(23, 19));
    glo_ephem->gamma = str2double(ephem_lines[0].substr(42, 19));

    // the second line
    glo_ephem->pos[0] = str2double(ephem_lines[1].substr(4, 19)) * 1e3;
    glo_ephem->vel[0] = str2double(ephem_lines[1].substr(23, 19)) * 1e3;
    glo_ephem->acc[0] = str2double(ephem_lines[1].substr(42, 19)) * 1e3;
    glo_ephem->health = static_cast<uint32_t>(str2double(ephem_lines[1].substr(61, 19)));

    // the third line
    glo_ephem->pos[1] = str2double(ephem_lines[2].substr(4, 19)) * 1e3;
    glo_ephem->vel[1] = str2double(ephem_lines[2].substr(23, 19)) * 1e3;
    glo_ephem->acc[1] = str2double(ephem_lines[2].substr(42, 19)) * 1e3;
    glo_ephem->freqo  = static_cast<int>(str2double(ephem_lines[2].substr(61, 19)));

    // the forth line
    glo_ephem->pos[2] = str2double(ephem_lines[3].substr(4, 19)) * 1e3;
    glo_ephem->vel[2] = str2double(ephem_lines[3].substr(23, 19)) * 1e3;
    glo_ephem->acc[2] = str2double(ephem_lines[3].substr(42, 19)) * 1e3;
    glo_ephem->age  = static_cast<uint32_t>(str2double(ephem_lines[3].substr(61, 19)));

    return glo_ephem;
}

void GNSSProcess::rinex2ephems(const std::string &rinex_filepath, std::map<uint32_t, std::vector<EphemBasePtr>> &sat2ephem_)
{
    uint32_t gpst_leap_seconds = static_cast<uint32_t>(-1);
    std::ifstream ephem_file(rinex_filepath);
    std::string line;
    while(std::getline(ephem_file, line))
    {
        if (line.find("RINEX VERSION / TYPE") != std::string::npos && line.find("3.04") == std::string::npos)
        {
            LOG(ERROR) << "Only RINEX 3.04 is supported for observation file";
            return;
        }
        else if (line.find("LEAP SECONDS") != std::string::npos && line.find("BDS") == std::string::npos)
            gpst_leap_seconds = static_cast<uint32_t>(std::stoi(line.substr(4, 6)));
        else if (line.find("END OF HEADER") != std::string::npos)
            break;
    }
    LOG_IF(FATAL, gpst_leap_seconds == static_cast<uint32_t>(-1)) << "No leap second record found";

    while(std::getline(ephem_file, line))
    {
        if (line.at(0) == 'G' || line.at(0) == 'C' || line.at(0) == 'E')
        {
            std::vector<std::string> ephem_lines;
            ephem_lines.push_back(line);
            for (size_t i = 0; i < 7; ++i)
            {
                std::getline(ephem_file, line);
                ephem_lines.push_back(line);
            }
            EphemPtr ephem = rinex_line2ephem(ephem_lines);
            if (!ephem || ephem->ttr.time == 0)  continue;
            if (sat2ephem_.count(ephem->sat) == 0)
                sat2ephem_.emplace(ephem->sat, std::vector<EphemBasePtr>());
            sat2ephem_.at(ephem->sat).push_back(ephem);
        }
        else if (line.at(0) == 'R')
        {
            std::vector<std::string> ephem_lines;
            ephem_lines.push_back(line);
            for (size_t i = 0; i < 3; ++i)
            {
                std::getline(ephem_file, line);
                ephem_lines.push_back(line);
            }
            GloEphemPtr glo_ephem = rinex_line2glo_ephem(ephem_lines, gpst_leap_seconds);
            if (sat2ephem_.count(glo_ephem->sat) == 0)
                sat2ephem_.emplace(glo_ephem->sat, std::vector<EphemBasePtr>());
            sat2ephem_.at(glo_ephem->sat).push_back(glo_ephem);
        }
    }
}

#endif


