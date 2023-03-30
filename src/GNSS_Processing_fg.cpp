#include "GNSS_Processing_fg.h"

GNSSProcess::GNSSProcess()
    : diff_t_gnss_local(0.0), gnss_track_num_threshold(20)
{
  Reset();
  // initNoises();
}

GNSSProcess::~GNSSProcess() {}

void GNSSProcess::Reset() 
{
  ROS_WARN("Reset GNSSProcess");
  std::map<std::pair<double, int>, std::map<uint32_t, double[6]>> empty_map_c;
  sat2cp.swap(empty_map_c);
  // sat2time_index.swap(empty_map_i);
  // sat2ephem.swap(empty_map_e);
  // latest_gnss_iono_params.swap(empty_vec_d);
  for (size_t i = 0; i < WINDOW_SIZE+1; i++)
  {
    std::vector<ObsPtr> empty_vec_o;
    std::vector<EphemBasePtr> empty_vec_e;
    gnss_meas_buf[i].swap(empty_vec_o);
    gnss_ephem_buf[i].swap(empty_vec_e);
  }

  std::map<uint32_t, uint32_t> empty_map_t;
  sat_track_status.swap(empty_map_t);
  gtSAMgraph.resize(0); 
  initialEstimate.clear();
  isamCurrentEstimate.clear();
  // index_delete = 0;
  frame_delete = 0;
  // E_num = 0;
  factor_id_frame.clear();
  id_accumulate = 0;
  frame_num = 0;
  last_gnss_time = 0.0;
  frame_count = 0;
  invalid_lidar = false;
  Rot_gnss_init.setIdentity();

  gnss_ready = false;
  // if (nolidar)
  {
    pre_integration->repropagate(Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero());
  }

  gtsam::ISAM2Params parameters;
  parameters.relinearizeThreshold = 0.1;
  parameters.relinearizeSkip = 1; // may matter? improtant!
  isam = gtsam::ISAM2(parameters);
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
  rinex2iono_params(rinex_filepath, latest_gnss_iono_params);
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
  std::vector<double> empty_vec_d;
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

Eigen::Vector3d GNSSProcess::local2enu(Eigen::Matrix3d R_enu_local_, Eigen::Vector3d anc, Eigen::Vector3d &pos)
{
  Eigen::Vector3d enu_pos;
  if (!nolidar)
  {
    enu_pos = R_enu_local_ * (pos - anc_local); // 

    // Eigen::Matrix3d R_ecef_enu_ = ecef2rotation(anc);
    // Eigen::Vector3d ecef_pos_ = anc + R_ecef_enu_ * enu_pos;
    Eigen::Vector3d ecef_pos_ = anc + enu_pos;
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

void GNSSProcess::initNoises( void ) // maybe usable!
{
    gtsam::Vector priorrotNoiseVector3(3);
    priorrotNoiseVector3 << prior_noise, prior_noise, prior_noise;
    priorrotNoise = gtsam::noiseModel::Diagonal::Variances(priorrotNoiseVector3);

    gtsam::Vector priorposNoiseVector12(12);
    priorposNoiseVector12 << prior_noise, prior_noise, prior_noise, prior_noise, prior_noise, prior_noise,
                            prior_noise, prior_noise, prior_noise, prior_noise, prior_noise, prior_noise;
    priorposNoise = gtsam::noiseModel::Diagonal::Variances(priorposNoiseVector12);

    // gtsam::Vector priorvelNoiseVector3(3);
    // priorvelNoiseVector3 << prior_noise, prior_noise, prior_noise;
    // priorvelNoise = gtsam::noiseModel::Diagonal::Variances(priorvelNoiseVector3);

    gtsam::Vector priorNoiseVector6(6);
    priorNoiseVector6 << prior_noise * 1, prior_noise * 1, prior_noise * 1, prior_noise * 1, prior_noise * 1, prior_noise * 1; 
    //, prior_noise, prior_noise, prior_noise, prior_noise, prior_noise, prior_noise, prior_noise, prior_noise;
    priorNoise = gtsam::noiseModel::Diagonal::Variances(priorNoiseVector6);

    gtsam::Vector priordtNoiseVector4(4);
    priordtNoiseVector4 << prior_noise, prior_noise, prior_noise, prior_noise;
    priordtNoise = gtsam::noiseModel::Diagonal::Variances(priordtNoiseVector4);

    // gtsam::Vector margExtNoiseVector4(4);
    // margExtNoiseVector4 << 1e-6, 1e-6, 1e-6, 1e-6;
    // margExtNoise = gtsam::noiseModel::Diagonal::Variances(margExtNoiseVector4);

    gtsam::Vector priorddtNoiseVector1(1);
    priorddtNoiseVector1 << prior_noise;
    priorddtNoise = gtsam::noiseModel::Diagonal::Variances(priorddtNoiseVector1);

    gtsam::Vector margrotNoiseVector3(3);
    margrotNoiseVector3 << marg_noise, marg_noise, marg_noise;
    margrotNoise = gtsam::noiseModel::Diagonal::Variances(margrotNoiseVector3);

    gtsam::Vector margposNoiseVector9(9);
    margposNoiseVector9 << marg_noise, marg_noise, marg_noise, marg_noise, marg_noise, marg_noise,
                            marg_noise, marg_noise, marg_noise; //, marg_noise, marg_noise, marg_noise;
    margposNoise = gtsam::noiseModel::Diagonal::Variances(margposNoiseVector9);

    // gtsam::Vector priorvelNoiseVector3(3);
    // priorvelNoiseVector3 << prior_noise, prior_noise, prior_noise;
    // priorvelNoise = gtsam::noiseModel::Diagonal::Variances(priorvelNoiseVector3);

    gtsam::Vector margNoiseVector3(3);
    margNoiseVector3 << prior_noise * 1, prior_noise * 1, prior_noise * 1; //, marg_noise, marg_noise, marg_noise, marg_noise, marg_noise, marg_noise, 
                        // marg_noise, marg_noise;
    margNoise = gtsam::noiseModel::Diagonal::Variances(margNoiseVector3);

    gtsam::Vector margdtNoiseVector4(4);
    margdtNoiseVector4 << marg_noise, marg_noise, marg_noise, marg_noise;
    margdtNoise = gtsam::noiseModel::Diagonal::Variances(margdtNoiseVector4);

    // gtsam::Vector margExtNoiseVector4(4);
    // margExtNoiseVector4 << 1e-6, 1e-6, 1e-6, 1e-6;
    // margExtNoise = gtsam::noiseModel::Diagonal::Variances(margExtNoiseVector4);

    gtsam::Vector margddtNoiseVector1(1);
    margddtNoiseVector1 << prior_noise;
    margddtNoise = gtsam::noiseModel::Diagonal::Variances(margddtNoiseVector1);

    gtsam::Vector dtNoiseVector4(4);
    dtNoiseVector4 << dt_noise, dt_noise, dt_noise, dt_noise;
    dtNoise = gtsam::noiseModel::Diagonal::Variances(dtNoiseVector4);

    gtsam::Vector ddtNoiseVector1(1);
    ddtNoiseVector1 << ddt_noise;
    ddtNoise = gtsam::noiseModel::Diagonal::Variances(ddtNoiseVector1);

    // gtsam::Vector odomNoiseVector9(9);
    // odomNoiseVector9 << odo_noise, odo_noise, odo_noise, odo_noise * 100, odo_noise * 100, odo_noise * 100, odo_noise, odo_noise, odo_noise;
    // odomNoise = gtsam::noiseModel::Diagonal::Variances(odomNoiseVector9); // should be related to the time, maybe proportional

    gtsam::Vector odomNoiseVector3(3);
    odomNoiseVector3 << odo_noise, odo_noise, odo_noise; //, odo_noise, odo_noise, odo_noise, odo_noise, odo_noise, odo_noise;
    odomNoise = gtsam::noiseModel::Diagonal::Variances(odomNoiseVector3); // should be related to the imu noise
    gtsam::Vector odomNoiseVector15(15);
    odomNoiseVector15 << odo_noise, odo_noise, odo_noise, odo_noise, odo_noise, odo_noise, odo_noise, odo_noise, odo_noise,
                        odo_noise, odo_noise, odo_noise, odo_noise, odo_noise, odo_noise;
    odomNoiseIMU = gtsam::noiseModel::Diagonal::Variances(odomNoiseVector15); // should be related to the imu noise
    // odomNoiseIMU = gtsam::noiseModel::Robust::Create(
    //                 gtsam::noiseModel::mEstimator::Cauchy::Create(1), // optional: replacing Cauchy by DCS or GemanMcClure is okay but Cauchy is empirically good.
    //                 gtsam::noiseModel::Diagonal::Variances(odomNoiseVector15));

    double psrNoiseScore = psr_dopp_noise; // constant is ok...
    double doppNoiseScore = psr_dopp_noise; // constant is ok...
    gtsam::Vector robustpsrdoppNoiseVector2(2); // gtsam::Pose3 factor has 6 elements (6D)
    robustpsrdoppNoiseVector2 << psrNoiseScore, doppNoiseScore;
    double cpNoiseScore = cp_noise; // 1e9
    gtsam::Vector robustcpNoiseVector1(1); // gps factor has 3 elements (xyz)
    robustcpNoiseVector1 << cpNoiseScore; // means only caring altitude here. (because LOAM-like-methods tends to be asymptotically flyging)
    if (outlier_rej)
    {
      robustpsrdoppNoise = gtsam::noiseModel::Robust::Create(
                      gtsam::noiseModel::mEstimator::Cauchy::Create(outlier_thres), // optional: replacing Cauchy by DCS or GemanMcClure is okay but Cauchy is empirically good.
                      // gtsam::noiseModel::mEstimator::Huber::Create(outlier_thres), // optional: replacing Cauchy by DCS or GemanMcClure is okay but Cauchy is empirically good.
                      gtsam::noiseModel::Diagonal::Variances(robustpsrdoppNoiseVector2));

      robustcpNoise = gtsam::noiseModel::Robust::Create(
                      gtsam::noiseModel::mEstimator::Cauchy::Create(outlier_thres), // optional: replacing Cauchy by DCS or GemanMcClure is okay but Cauchy is empirically good.
                      // gtsam::noiseModel::mEstimator::Huber::Create(outlier_thres), // optional: replacing Cauchy by DCS or GemanMcClure is okay but Cauchy is empirically good.
                      gtsam::noiseModel::Diagonal::Variances(robustcpNoiseVector1));
    }
    else
    {
      robustpsrdoppNoise = gtsam::noiseModel::Diagonal::Variances(robustpsrdoppNoiseVector2);
      robustcpNoise = gtsam::noiseModel::Diagonal::Variances(robustcpNoiseVector1);
    }
    // testNoise = gtsam::noiseModel::Gaussian::Covariance()
    // robustpsrdoppNoise = gtsam::noiseModel::Diagonal::Variances(robustpsrdoppNoiseVector2);

} // initNoises

void GNSSProcess::processGNSS(const std::vector<ObsPtr> &gnss_meas, state_input &state, Eigen::Vector3d &omg)
{
  std::vector<ObsPtr> valid_meas;
  std::vector<EphemBasePtr> valid_ephems;
  if (gnss_meas.empty())  
  {
    if (gnss_ready)
    {
      std::vector<ObsPtr> empty_vec_o;
      std::vector<EphemBasePtr> empty_vec_e;
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
      if (obs->psr_std[freq_idx]  > gnss_psr_std_threshold / 3 ||
          obs->dopp_std[freq_idx] > gnss_dopp_std_threshold / 3) //||
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
          // if (frame_num < gnss_track_num_threshold && !gnss_online_init && nolidar)
          // {
          //   last_gnss_time = time2sec(obs->time);
          //   pre_integration->repropagate(state.ba, state.bg);
          // }
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
          Eigen::Vector3d pos_gnss = state.pos + state.rot.normalized().toRotationMatrix() * Tex_imu_r;
          updateGNSSStatistics(pos_gnss);
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
    // if (!nolidar)
    // {
    //   rot_window[frame_count] = state.rot;
    //   pos_window[frame_count] = state.pos + state.rot * Tex_imu_r;
    //   Eigen::Matrix3d omg_skew;
    //   omg_skew << SKEW_SYM_MATRX(omg);
    //   vel_window[frame_count] = state.vel + omg_skew * Tex_imu_r;
    // }
    // else
    {
      rot_window[frame_count] = state.rot.normalized().toRotationMatrix();
      pos_window[frame_count] = state.pos + state.rot.normalized().toRotationMatrix() * Tex_imu_r;
      Eigen::Matrix3d omg_skew;
      omg_skew << SKEW_SYM_MATRX(omg);
      vel_window[frame_count] = state.vel + state.rot.normalized().toRotationMatrix() * omg_skew * Tex_imu_r;
      // vel_window[frame_count] = state.vel;
    }
    gnss_meas_buf[frame_count] = valid_meas; 
    gnss_ephem_buf[frame_count] = valid_ephems;
    frame_count ++;
    gnss_ready = GNSSLIAlign();
    if (gnss_ready)
    {
      ROS_INFO("GNSS Initialization is done");
      state_ = state;
      state_last = state;
    }
  }
  else
  {  
    gnss_meas_buf[0] = valid_meas; 
    gnss_ephem_buf[0] = valid_ephems;
  }
}

void GNSSProcess::processGNSS(const std::vector<ObsPtr> &gnss_meas, state_output &state)
{
  std::vector<ObsPtr> valid_meas;
  std::vector<EphemBasePtr> valid_ephems;
  if (gnss_meas.empty())  
  {
    if (gnss_ready)
    {
      std::vector<ObsPtr> empty_vec_o;
      std::vector<EphemBasePtr> empty_vec_e;
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
      if (obs->psr_std[freq_idx]  > gnss_psr_std_threshold / 3 ||
          obs->dopp_std[freq_idx] > gnss_dopp_std_threshold / 3) //||
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
          // if (frame_num < gnss_track_num_threshold && !gnss_online_init && nolidar)
          // {
          //   last_gnss_time = time2sec(obs->time);
          //   pre_integration->repropagate(state.ba, state.bg);
          // }
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
          Eigen::Vector3d pos_gnss = state.pos + state.rot.normalized().toRotationMatrix() * Tex_imu_r;
          updateGNSSStatistics(pos_gnss);
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
    // if (!nolidar)
    // {
    //   rot_window[frame_count] = state.rot;
    //   pos_window[frame_count] = state.pos + state.rot * Tex_imu_r;
    //   Eigen::Matrix3d omg_skew;
    //   omg_skew << SKEW_SYM_MATRX(state.omg);
    //   vel_window[frame_count] = state.vel + omg_skew * Tex_imu_r;
    // }
    // else
    {
      rot_window[frame_count] = state.rot.normalized().toRotationMatrix();
      pos_window[frame_count] = state.pos + state.rot.normalized().toRotationMatrix() * Tex_imu_r;
      Eigen::Matrix3d omg_skew;
      omg_skew << SKEW_SYM_MATRX(state.omg);
      vel_window[frame_count] = state.vel + state.rot.normalized().toRotationMatrix() * omg_skew * Tex_imu_r;
      // vel_window[frame_count] = state.vel;
    }
    gnss_meas_buf[frame_count] = valid_meas; 
    gnss_ephem_buf[frame_count] = valid_ephems;
    frame_count ++;
    gnss_ready = GNSSLIAlign();
    if (gnss_ready)
    {
      state_const_ = state;
      state_const_last = state;
    }
  }
  else
  {  
    gnss_meas_buf[0] = valid_meas; 
    gnss_ephem_buf[0] = valid_ephems;
  }
}

void GNSSProcess::runISAM2opt(void) //
{
  gtsam::FactorIndices delete_factor;
  gtsam::FactorIndices().swap(delete_factor);

  if (gnss_ready)
  {
    bool delete_happen = false;
    if (frame_num - frame_delete > delete_thred) // (graph_whole1.size() - index_delete > 4000)
    {
      delete_happen = true;
    while (frame_num - frame_delete > delete_thred) // (graph_whole1.size() - index_delete > 3000)
    { 
      if (!factor_id_frame.empty())       
      {
        // if (frame_delete > 0)
        {
        for (size_t i = 0; i < factor_id_frame[0].size(); i++)
        {
          // if (factor_id_frame[0][i] != 0 && factor_id_frame[0][i] != 1 || nolidar)
          {
            delete_factor.push_back(factor_id_frame[0][i]);
          }
        }
        // index_delete += factor_id_frame[0].size();
        }
      
        factor_id_frame.pop_front();
        frame_delete ++;
      }
      if (factor_id_frame.empty()) break;
    }
    }

    if (delete_happen && !nolidar)
    {
      if (frame_delete > 0)
      {
        if (frame_delete >= change_ext)
        {
        gtsam::noiseModel::Gaussian::shared_ptr updatedERNoise = gtsam::noiseModel::Gaussian::Covariance(isam.marginalCovariance(P(0)) * 1); // important
        gtsam::noiseModel::Gaussian::shared_ptr updatedEPNoise = gtsam::noiseModel::Gaussian::Covariance(isam.marginalCovariance(E(0)) * 1); // important
        gtsam::PriorFactor<gtsam::Rot3> init_ER(P(0),isamCurrentEstimate.at<gtsam::Rot3>(P(0)), updatedERNoise); //  margrotNoise); //
        gtsam::PriorFactor<gtsam::Vector3> init_EP(E(0),isamCurrentEstimate.at<gtsam::Vector3>(E(0)), updatedEPNoise); // margrotNoise); // 
        gtSAMgraph.add(init_ER);
        gtSAMgraph.add(init_EP);
        // factor_id_frame[0].push_back(id_accumulate);
        // factor_id_frame[0].push_back(id_accumulate+1);
        factor_id_frame[frame_num - 1 - frame_delete].push_back(id_accumulate);
        factor_id_frame[frame_num - 1 - frame_delete].push_back(id_accumulate+1);
        id_accumulate += 2;
        change_ext = frame_num;
        }
      size_t j = 0;
      for (; j < marg_thred; j++)
      {
        // get updated noise before reset
        gtsam::noiseModel::Gaussian::shared_ptr updatedRotNoise = gtsam::noiseModel::Gaussian::Covariance(isam.marginalCovariance(R(frame_delete+j))); // important
        gtsam::noiseModel::Gaussian::shared_ptr updatedPosNoise = gtsam::noiseModel::Gaussian::Covariance(isam.marginalCovariance(A(frame_delete+j))); // important
        // gtsam::noiseModel::Gaussian::shared_ptr updatedPosNoise = gtsam::noiseModel::Gaussian::Covariance(isam.marginalCovariance(F(frame_delete+j))); // important
        gtsam::noiseModel::Gaussian::shared_ptr updatedDtNoise = gtsam::noiseModel::Gaussian::Covariance(isam.marginalCovariance(B(frame_delete+j))); // important
        gtsam::noiseModel::Gaussian::shared_ptr updatedDdtNoise = gtsam::noiseModel::Gaussian::Covariance(isam.marginalCovariance(C(frame_delete+j))); // important

        gtsam::PriorFactor<gtsam::Rot3> init_rot(R(frame_delete+j),isamCurrentEstimate.at<gtsam::Rot3>(R(frame_delete+j)), updatedRotNoise); // margrotNoise);
        // gtsam::PriorFactor<gtsam::Vector12> init_vel(F(frame_delete+j), isamCurrentEstimate.at<gtsam::Vector12>(F(frame_delete+j)), updatedPosNoise); // margposNoise);
        gtsam::PriorFactor<gtsam::Vector6> init_vel(A(frame_delete+j), isamCurrentEstimate.at<gtsam::Vector6>(A(frame_delete+j)), updatedPosNoise); // margposNoise);
        gtsam::PriorFactor<gtsam::Vector4> init_dt(B(frame_delete+j), isamCurrentEstimate.at<gtsam::Vector4>(B(frame_delete+j)), updatedDtNoise); // margdtNoise);
        gtsam::PriorFactor<gtsam::Vector1> init_ddt(C(frame_delete+j), isamCurrentEstimate.at<gtsam::Vector1>(C(frame_delete+j)), updatedDdtNoise); // margddtNoise);
        gtSAMgraph.add(init_rot);
        gtSAMgraph.add(init_vel);
        gtSAMgraph.add(init_dt);
        gtSAMgraph.add(init_ddt);
        factor_id_frame[0].push_back(id_accumulate+(j)*4);
        factor_id_frame[0].push_back(id_accumulate+1+(j)*4);
        factor_id_frame[0].push_back(id_accumulate+2+(j)*4);
        factor_id_frame[0].push_back(id_accumulate+3+(j)*4);
      }
      // id_accumulate += (j-1) * 4;
      id_accumulate += j * 4;
      }
      isam.update(gtSAMgraph, initialEstimate);
      gtSAMgraph.resize(0); // will the initialEstimate change?
      initialEstimate.clear();
      isam.update(gtSAMgraph, initialEstimate, delete_factor);   
    }
    else if (delete_happen && nolidar)
    {
      if (frame_delete > 0) // (frame_delete == 0)
      {
      size_t j = 0;
      // for (; j < 10; j++)
      for (; j < marg_thred; j++)
      {
        gtsam::noiseModel::Gaussian::shared_ptr updatedRotNoise = gtsam::noiseModel::Gaussian::Covariance(isam.marginalCovariance(R(frame_delete+j))); // important
        gtsam::noiseModel::Gaussian::shared_ptr updatedPosNoise = gtsam::noiseModel::Gaussian::Covariance(isam.marginalCovariance(F(frame_delete+j))); // important
        gtsam::noiseModel::Gaussian::shared_ptr updatedDtNoise = gtsam::noiseModel::Gaussian::Covariance(isam.marginalCovariance(B(frame_delete+j))); // important
        gtsam::noiseModel::Gaussian::shared_ptr updatedDdtNoise = gtsam::noiseModel::Gaussian::Covariance(isam.marginalCovariance(C(frame_delete+j))); // important

        gtsam::PriorFactor<gtsam::Rot3> init_rot(R(frame_delete+j),isamCurrentEstimate.at<gtsam::Rot3>(R(frame_delete+j)), updatedRotNoise); // margrotNoise);
        gtsam::PriorFactor<gtsam::Vector12> init_vel(F(frame_delete+j), isamCurrentEstimate.at<gtsam::Vector12>(F(frame_delete+j)), updatedPosNoise); // margposNoise);
        gtsam::PriorFactor<gtsam::Vector4> init_dt(B(frame_delete+j), isamCurrentEstimate.at<gtsam::Vector4>(B(frame_delete+j)), updatedDtNoise); // margdtNoise);
        gtsam::PriorFactor<gtsam::Vector1> init_ddt(C(frame_delete+j), isamCurrentEstimate.at<gtsam::Vector1>(C(frame_delete+j)), updatedDdtNoise); // margddtNoise); // could delete?
        gtSAMgraph.add(init_rot);
        gtSAMgraph.add(init_vel);
        gtSAMgraph.add(init_dt);
        gtSAMgraph.add(init_ddt);
        
        {
          factor_id_frame[0].push_back(id_accumulate+j*4);
          factor_id_frame[0].push_back(id_accumulate+1+j*4);
          factor_id_frame[0].push_back(id_accumulate+2+j*4);
          factor_id_frame[0].push_back(id_accumulate+3+j*4);
          // factor_id_frame[0].push_back(id_accumulate+4+j*4);
        }
      }
      id_accumulate += j * 4;
      }
      isam.update(gtSAMgraph, initialEstimate);
      gtSAMgraph.resize(0); // will the initialEstimate change?
      initialEstimate.clear();
      isam.update(gtSAMgraph, initialEstimate, delete_factor);
      // gtSAMgraph.resize(0); // will the initialEstimate change?
      // initialEstimate.clear();
      // isam.update();
    }
    else
    {
      isam.update(gtSAMgraph, initialEstimate); //, delete_factor);

      gtSAMgraph.resize(0); // will the initialEstimate change?
      initialEstimate.clear();
      isam.update();
    }
  }
  else
  {
    isam.update(gtSAMgraph, initialEstimate);

    gtSAMgraph.resize(0); // will the initialEstimate change?
    initialEstimate.clear();
    isam.update();
  }
  isamCurrentEstimate = isam.calculateEstimate();
  
  if (nolidar) // || invalid_lidar)
  {
    pre_integration->repropagate(isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(6),
                                isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(9));
  }
  else
  {
    pre_integration->repropagate(Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero());
  }
  std::map<std::pair<double, int>, std::map<uint32_t, double[6]>>::iterator it_old;
  if (!sat2cp.empty())
  {
    it_old = sat2cp.begin();
    while (it_old->first.second < frame_delete)
    {
      std::map<uint32_t, double[6]>().swap(it_old->second);
      size_t del_size = sat2cp.erase(it_old->first);
      if (sat2cp.empty()) break;
      it_old = sat2cp.begin();
    }
  }
}

bool GNSSProcess::GNSSLIAlign()
{
  for (uint32_t i = 0; i < (wind_size+1); i++)
  {
    if (gnss_meas_buf[i].empty() || gnss_meas_buf[i].size() < 5) // need IMU to prop
    {
      // if (frame_count == wind_size + 1)
      // {
        if (i >= frame_count)
        {
          return false;
        }
        uint32_t j = i;
        bool shift_happen = false;
        for (; j < (wind_size); j++)
        {
          shift_happen = true;
          gnss_meas_buf[j] = gnss_meas_buf[j+1];
          gnss_ephem_buf[j] = gnss_ephem_buf[j+1];
          rot_window[j] = rot_window[j+1];
          pos_window[j] = pos_window[j+1];
          vel_window[j] = vel_window[j+1];
        }
        if (shift_happen || j == wind_size)
        { 
          frame_count -= 1;
          std::vector<ObsPtr> empty_vec_o;
          std::vector<EphemBasePtr> empty_vec_e;
          gnss_meas_buf[wind_size].swap(empty_vec_o); // maybe problem
          gnss_ephem_buf[wind_size].swap(empty_vec_e);          
        }    
      // }
      return false;
    }
  }
  
  for (uint32_t i = 0; i < wind_size; i++)
  {
    if (time2sec(gnss_meas_buf[i+1][0]->time) - time2sec(gnss_meas_buf[i][0]->time) > 15 * gnss_sample_period) // need IMU to prop
    {
      // if (frame_count == wind_size + 1)
      // {
        for (uint32_t j = i+1; j < wind_size+1; ++j)
        {
          gnss_meas_buf[j-i-1] = gnss_meas_buf[j];
          gnss_ephem_buf[j-i-1] = gnss_ephem_buf[j];
          rot_window[j-i-1] = rot_window[j];
          pos_window[j-i-1] = pos_window[j];
          vel_window[j-i-1] = vel_window[j];
        }
        frame_count -= i+1;
        for (uint32_t j = wind_size-i; j < wind_size+1; ++j)
        {
          std::vector<ObsPtr> empty_vec_o;
          std::vector<EphemBasePtr> empty_vec_e;
          gnss_meas_buf[j].swap(empty_vec_o);
          gnss_ephem_buf[j].swap(empty_vec_e); 
        }             
      // }
      return false;
    }
  }

  // check horizontal velocity excitation
  if (!quick_init)
  {
    Eigen::Vector2d avg_hor_vel(0.0, 0.0);
    for (uint32_t i = 0; i < (wind_size+1); ++i)
        avg_hor_vel += vel_window[i].head<2>().cwiseAbs();
    avg_hor_vel /= (wind_size+1);
    if (avg_hor_vel.norm() < 0.3) // if (0) ?
    {
      std::cerr << "velocity excitation not enough for GNSS-LI alignment.\n";
      for (uint32_t i = 0; i < (wind_size); ++i)
      {
        gnss_meas_buf[i] = gnss_meas_buf[i+1];
        gnss_ephem_buf[i] = gnss_ephem_buf[i+1];
        rot_window[i] = rot_window[i+1];
        pos_window[i] = pos_window[i+1];
        vel_window[i] = vel_window[i+1];
      }
      frame_count = wind_size;
      std::vector<ObsPtr> empty_vec_o;
      std::vector<EphemBasePtr> empty_vec_e;
      gnss_meas_buf[frame_count].swap(empty_vec_o);
      gnss_ephem_buf[frame_count].swap(empty_vec_e);    
      return false;
    }
  }
  std::vector<std::vector<ObsPtr>> curr_gnss_meas_buf;
  std::vector<std::vector<EphemBasePtr>> curr_gnss_ephem_buf;
  for (uint32_t i = 0; i < (wind_size+1); ++i)
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
      for (uint32_t i = 0; i < (wind_size); ++i)
      {
        gnss_meas_buf[i] = gnss_meas_buf[i+1];
        gnss_ephem_buf[i] = gnss_ephem_buf[i+1];
        rot_window[i] = rot_window[i+1];
        pos_window[i] = pos_window[i+1];
        vel_window[i] = vel_window[i+1];
      }
      frame_count = wind_size;
      std::vector<ObsPtr> empty_vec_o;
      std::vector<EphemBasePtr> empty_vec_e;
      gnss_meas_buf[frame_count].swap(empty_vec_o);
      gnss_ephem_buf[frame_count].swap(empty_vec_e);
      return false;
  }

  if (!quick_init)
  {
    // 2. perform yaw alignment
    std::vector<Eigen::Vector3d> local_vs;
    for (uint32_t i = 0; i < (wind_size+1); ++i)
        local_vs.push_back(vel_window[i]); // values at gnss measurement
    Eigen::Vector3d rough_anchor_ecef = rough_xyzt.head<3>();
    double aligned_yaw = 0;
    double aligned_rcv_ddt = 0;
    if (!gnss_li_initializer.yaw_alignment(local_vs, rough_anchor_ecef, aligned_yaw, aligned_rcv_ddt))
    {
        std::cerr << "Fail to align ENU and local frames.\n";
        for (uint32_t i = 0; i < (wind_size); ++i)
        {
          gnss_meas_buf[i] = gnss_meas_buf[i+1];
          gnss_ephem_buf[i] = gnss_ephem_buf[i+1];

          rot_window[i] = rot_window[i+1];
          pos_window[i] = pos_window[i+1];
          vel_window[i] = vel_window[i+1];
        }
        frame_count = wind_size;
        std::vector<ObsPtr> empty_vec_o;
        std::vector<EphemBasePtr> empty_vec_e;
        gnss_meas_buf[frame_count].swap(empty_vec_o);
        gnss_ephem_buf[frame_count].swap(empty_vec_e);
        return false;
    }

    // 3. perform anchor refinement
    std::vector<Eigen::Vector3d> local_ps;
    for (uint32_t i = 0; i < (wind_size+1); ++i)
        local_ps.push_back(pos_window[i]); // values at gnss measurement
    Eigen::Matrix<double, 7, 1> refined_xyzt;
    refined_xyzt.setZero();
    if (!gnss_li_initializer.anchor_refinement(local_ps, aligned_yaw, 
        aligned_rcv_ddt, rough_xyzt, refined_xyzt))
    {
        std::cerr << "Fail to refine anchor point.\n";
        for (uint32_t i = 0; i < (wind_size); ++i)
        {
          gnss_meas_buf[i] = gnss_meas_buf[i+1];
          gnss_ephem_buf[i] = gnss_ephem_buf[i+1];

          rot_window[i] = rot_window[i+1];
          pos_window[i] = pos_window[i+1];
          vel_window[i] = vel_window[i+1];
        }
        frame_count = wind_size;
        std::vector<ObsPtr> empty_vec_o;
        std::vector<EphemBasePtr> empty_vec_e;
        gnss_meas_buf[frame_count].swap(empty_vec_o);
        gnss_ephem_buf[frame_count].swap(empty_vec_e);
        return false;
    }

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
    for (uint32_t i = 0; i < (wind_size+1); ++i)
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
    anc_local = pos_window[0]; // [WINDOW_SIZE]; // ?
    yaw_enu_local = aligned_yaw; // can we align the rot matrix?
    R_ecef_enu = ecef2rotation(anc_ecef);
    //
    if (!nolidar)
    {
      anc_local = Rot_gnss_init.transpose() * anc_local;
      Eigen::Matrix3d R_enu_local_;
      R_enu_local_ = R_ecef_enu * Eigen::AngleAxisd(yaw_enu_local, Eigen::Vector3d::UnitZ()) * Rot_gnss_init;
      // double init_value[8];
      // init_value[0] = para_rcv_dt[0]; init_value[1] = para_rcv_dt[1]; init_value[2] = para_rcv_dt[2]; init_value[3] = para_rcv_dt[3];
      // init_value[4] = para_rcv_ddt[0]; init_value[5] = anc_ecef[0]; init_value[6] = anc_ecef[1]; init_value[7] = anc_ecef[2];
      // mtxPosegraph.lock();
      {
        // prior factor 
        Eigen::Matrix<double, 6, 1> init_vel_bias_vector;
        init_vel_bias_vector.block<3,1>(0,0) = Rot_gnss_init.transpose() * pos_window[0];
        init_vel_bias_vector.block<3,1>(3,0) = Rot_gnss_init.transpose() * vel_window[0];
        // init_vel_bias_vector.block<6,1>(6,0) = Eigen::Matrix<double, 6, 1>::Zero();
        initialEstimate.insert(R(0), gtsam::Rot3(Rot_gnss_init.transpose() * rot_window[0]));
        // initialEstimate.insert(F(0), gtsam::Vector12(init_vel_bias_vector));
        initialEstimate.insert(A(0), gtsam::Vector6(init_vel_bias_vector));
        initialEstimate.insert(B(0), gtsam::Vector4(para_rcv_dt[0], para_rcv_dt[1], para_rcv_dt[2], para_rcv_dt[3]));
        initialEstimate.insert(C(0), gtsam::Vector1(para_rcv_ddt[0]));
        // initialEstimate.insert(Y(0), gtsam::Vector1(yaw_enu_local));
        initialEstimate.insert(E(0), gtsam::Vector3(anc_ecef[0], anc_ecef[1], anc_ecef[2]));
        initialEstimate.insert(P(0), gtsam::Rot3(R_enu_local_));
        // gtSAMgraph.add(glio::PriorFactor(B(0), C(0), E(E_num), P(E_num), init_value, R_enu_local_, priorNoise));

        gtsam::PriorFactor<gtsam::Rot3> init_rot_ext(P(0), gtsam::Rot3(gtsam::Rot3(R_enu_local_)), priorrotNoise);
        gtsam::PriorFactor<gtsam::Vector3> init_pos_ext(E(0), gtsam::Vector3(anc_ecef[0], anc_ecef[1], anc_ecef[2]), margNoise);
        gtsam::PriorFactor<gtsam::Vector4> init_dt(B(0), gtsam::Vector4(para_rcv_dt[0], para_rcv_dt[1], para_rcv_dt[2], para_rcv_dt[3]), priordtNoise);
        gtsam::PriorFactor<gtsam::Vector1> init_ddt(C(0), gtsam::Vector1(para_rcv_ddt[0]), priorddtNoise);
        gtsam::PriorFactor<gtsam::Rot3> init_rot_(R(0), gtsam::Rot3(Rot_gnss_init.transpose() * rot_window[0]), priorrotNoise);
        gtsam::PriorFactor<gtsam::Vector6> init_vel_(A(0), gtsam::Vector6(init_vel_bias_vector), priorNoise); // priorposNoise);
        // gtsam::PriorFactor<gtsam::Vector12> init_vel_(F(0), gtsam::Vector12(init_vel_bias_vector), priorposNoise);
        gtSAMgraph.add(init_rot_ext);
        gtSAMgraph.add(init_pos_ext);
        gtSAMgraph.add(init_dt);
        gtSAMgraph.add(init_ddt);
        gtSAMgraph.add(init_rot_);
        gtSAMgraph.add(init_vel_);
        factor_id_frame.push_back(std::vector<size_t>{0, 1, 2, 3, 4, 5});
        // time_frame.push_back(std::pair<double, int>(0.0, 0));

        for (size_t i = 1; i < wind_size+1; i++)
        {
          init_vel_bias_vector.block<3,1>(3,0) = Rot_gnss_init.transpose() * vel_window[i];
          init_vel_bias_vector.block<3,1>(0,0) = Rot_gnss_init.transpose() * pos_window[i];
          initialEstimate.insert(R(i), gtsam::Rot3(Rot_gnss_init.transpose() * rot_window[i]));
          // initialEstimate.insert(F(i), gtsam::Vector12(init_vel_bias_vector));
          initialEstimate.insert(A(i), gtsam::Vector6(init_vel_bias_vector));
          initialEstimate.insert(B(i), gtsam::Vector4(para_rcv_dt[i*4], para_rcv_dt[i*4+1], para_rcv_dt[i*4+2], para_rcv_dt[i*4+3]));
          initialEstimate.insert(C(i), gtsam::Vector1(para_rcv_ddt[i]));

          gtsam::PriorFactor<gtsam::Rot3> init_rot(R(i), gtsam::Rot3(Rot_gnss_init.transpose() * rot_window[i]), priorrotNoise);
          // gtsam::PriorFactor<gtsam::Vector12> init_vel(F(i), gtsam::Vector12(init_vel_bias_vector), priorposNoise);
          gtsam::PriorFactor<gtsam::Vector6> init_vel(A(i), gtsam::Vector6(init_vel_bias_vector), priorNoise); // priorposNoise);
          gtsam::PriorFactor<gtsam::Vector4> init_dt(B(i), gtsam::Vector4(para_rcv_dt[i*4], para_rcv_dt[i*4+1], para_rcv_dt[i*4+2], para_rcv_dt[i*4+3]), priordtNoise);
          gtsam::PriorFactor<gtsam::Vector1> init_ddt(C(i), gtsam::Vector1(para_rcv_ddt[i]), priorddtNoise);
          gtSAMgraph.add(init_rot);
          gtSAMgraph.add(init_vel);
          gtSAMgraph.add(init_dt);
          gtSAMgraph.add(init_ddt);
          factor_id_frame.push_back(std::vector<size_t>{6 + 4 * (i-1), 7 + 4 * (i-1), 8 + 4 * (i-1), 9 + 4 * (i-1)}); //, 7 + 4 * (i-1)});
          // time_frame.push_back(std::pair<double, int>(0.0, i));
        }
        id_accumulate = 4 * wind_size + 6;
      }   
    }
    else
    {
      Eigen::Matrix3d R_enu_local_;
      R_enu_local_ = Eigen::AngleAxisd(yaw_enu_local, Eigen::Vector3d::UnitZ());

        for (size_t i = 0; i < wind_size+1; i++)
        {
          gtsam::PriorFactor<gtsam::Rot3> init_rot(R(i), gtsam::Rot3(R_ecef_enu * R_enu_local_ * rot_window[i]), priorrotNoise);
          // Eigen::Vector3d init_vel = ;
          Eigen::Matrix<double, 12, 1> init_vel_bias_vector;
          init_vel_bias_vector.block<3,1>(0,0) = anc_ecef + R_ecef_enu * R_enu_local_ * (pos_window[i] - pos_window[0] - rot_window[i] * Tex_imu_r);
          init_vel_bias_vector.block<3,1>(3,0) = R_ecef_enu * R_enu_local_ * vel_window[i];
          init_vel_bias_vector.block<6,1>(6,0) = Eigen::Matrix<double, 6, 1>::Zero();
          gtsam::PriorFactor<gtsam::Vector12> init_vel_bias(F(i), gtsam::Vector12(init_vel_bias_vector), priorposNoise);
          gtsam::PriorFactor<gtsam::Vector4> init_dt(B(i), gtsam::Vector4(para_rcv_dt[i*4], para_rcv_dt[i*4+1], para_rcv_dt[i*4+2], para_rcv_dt[i*4+3]), priordtNoise);
          gtsam::PriorFactor<gtsam::Vector1> init_ddt(C(i), gtsam::Vector1(para_rcv_ddt[i]), priorddtNoise);
          gtSAMgraph.add(init_rot);
          gtSAMgraph.add(init_vel_bias);
          gtSAMgraph.add(init_dt);
          gtSAMgraph.add(init_ddt);
          factor_id_frame.push_back(std::vector<size_t>{i * 4, i * 4 + 1, i * 4  + 2, i * 4 + 3});
          // time_frame.push_back(std::pair<double, int>(0.0, i));
          initialEstimate.insert(R(i), gtsam::Rot3(R_ecef_enu * R_enu_local_ * rot_window[i]));
          initialEstimate.insert(F(i), gtsam::Vector12(init_vel_bias_vector));
          initialEstimate.insert(B(i), gtsam::Vector4(para_rcv_dt[i*4], para_rcv_dt[i*4+1], para_rcv_dt[i*4+2], para_rcv_dt[i*4+3]));
          initialEstimate.insert(C(i), gtsam::Vector1(para_rcv_ddt[i]));
        }
        id_accumulate = 4 + 4 * wind_size;
    }
    // mtxPosegraph.unlock();
    last_gnss_time = time2sec(gnss_meas_buf[wind_size][0]->time);
    frame_num = frame_count;
  }
  else
  {
    for (uint32_t k = 0; k < 4; k++)
    {
      if (rough_xyzt(3+k) == 0)
      {
        std::cerr << "Fail to quick init anchor point.\n";
        for (uint32_t i = 0; i < (wind_size); ++i)
        {
          gnss_meas_buf[i] = gnss_meas_buf[i+1]; // change the strategy
          gnss_ephem_buf[i] = gnss_ephem_buf[i+1];

          rot_window[i] = rot_window[i+1];
          pos_window[i] = pos_window[i+1];
          vel_window[i] = vel_window[i+1];
        }
        frame_count = wind_size;
        std::vector<ObsPtr> empty_vec_o;
        std::vector<EphemBasePtr> empty_vec_e;
        gnss_meas_buf[frame_count].swap(empty_vec_o);
        gnss_ephem_buf[frame_count].swap(empty_vec_e);
        return false;
      }
    }
    anc_ecef = rough_xyzt.head<3>();
    anc_local = pos_window[0]; // [WINDOW_SIZE]; // ?
    yaw_enu_local = 0.0;
    R_ecef_enu = ecef2rotation(anc_ecef);
    para_rcv_ddt[wind_size] = 128.0;
    for (uint32_t k_ = 0; k_ < 4; ++k_)
    {
      para_rcv_dt[wind_size*4+k_] = rough_xyzt(3+k_);
    }

    // use gtsam
    if (!nolidar)
    {
      anc_local = Rot_gnss_init.transpose() * anc_local;
      // mtxPosegraph.lock();
      {
        Eigen::Matrix3d R_enu_local_;
        R_enu_local_ = R_ecef_enu * Eigen::AngleAxisd(yaw_enu_local, Eigen::Vector3d::UnitZ()) * Rot_gnss_init;   
        Eigen::Matrix<double, 6, 1> init_vel_bias_vector;
        init_vel_bias_vector.block<3,1>(0,0) = Rot_gnss_init.transpose() * pos_window[0];
        init_vel_bias_vector.block<3,1>(3,0) = Rot_gnss_init.transpose() * vel_window[0];
        // init_vel_bias_vector.block<6,1>(6,0) = Eigen::Matrix<double, 6, 1>::Zero();
        gtsam::PriorFactor<gtsam::Vector3> init_anc(E(0), gtsam::Vector3(anc_ecef[0], anc_ecef[1], anc_ecef[2]),margNoise); // a bug
        gtsam::PriorFactor<gtsam::Rot3> init_ext_rot(P(0), gtsam::Rot3(R_enu_local_),priorrotNoise); // a bug
        gtSAMgraph.add(init_anc);
        gtSAMgraph.add(init_ext_rot);
        initialEstimate.insert(E(0), gtsam::Vector3(anc_ecef[0], anc_ecef[1], anc_ecef[2]));
        initialEstimate.insert(P(0), gtsam::Rot3(R_enu_local_));
        // gtsam::PriorFactor<gtsam::Vector12> init_pos(F(0), gtsam::Vector12(init_vel_bias_vector), priorposNoise);
        gtsam::PriorFactor<gtsam::Vector6> init_pos(A(0), gtsam::Vector6(init_vel_bias_vector), priorNoise); // priorposNoise);
        gtsam::PriorFactor<gtsam::Rot3> init_rot(R(0), gtsam::Rot3(Rot_gnss_init.transpose() * rot_window[0]), priorrotNoise);
        gtsam::PriorFactor<gtsam::Vector4> init_dt(B(0), gtsam::Vector4(para_rcv_dt[0], para_rcv_dt[1], para_rcv_dt[2], para_rcv_dt[3]), priordtNoise);
        factor_id_frame.push_back(std::vector<size_t>{0, 1, 2, 3, 4});

        initialEstimate.insert(R(0), gtsam::Rot3(Rot_gnss_init.transpose() * rot_window[0]));
        initialEstimate.insert(A(0), gtsam::Vector6(init_vel_bias_vector));
        // initialEstimate.insert(F(0), gtsam::Vector12(init_vel_bias_vector));
        initialEstimate.insert(B(0), gtsam::Vector4(para_rcv_dt[0], para_rcv_dt[1], para_rcv_dt[2], para_rcv_dt[3]));
        initialEstimate.insert(C(0), gtsam::Vector1(para_rcv_ddt[0]));

        for (size_t i = 1; i < wind_size+1; i++)
        {
          init_vel_bias_vector.block<3,1>(0,0) = Rot_gnss_init.transpose() * pos_window[i];
          init_vel_bias_vector.block<3,1>(3,0) = Rot_gnss_init.transpose() * vel_window[i];
          gtsam::PriorFactor<gtsam::Rot3> init_rot(R(i), gtsam::Rot3(Rot_gnss_init.transpose() * rot_window[i]), priorrotNoise);
          gtsam::PriorFactor<gtsam::Vector6> init_pos(A(i), gtsam::Vector6(init_vel_bias_vector), priorNoise); //priorposNoise);
          gtsam::PriorFactor<gtsam::Vector4> init_dt(B(i), gtsam::Vector4(para_rcv_dt[i*4], para_rcv_dt[i*4+1], para_rcv_dt[i*4+2], para_rcv_dt[i*4+3]), priordtNoise);
          
          gtSAMgraph.add(init_rot); // no rotation as well, so pos should be divided into rot and trans, or cov should be modified.
          gtSAMgraph.add(init_pos);
          gtSAMgraph.add(init_dt);
          factor_id_frame.push_back(std::vector<size_t>{2 + i * 3, i * 3 + 3, i * 3 + 4});
          initialEstimate.insert(R(i), gtsam::Rot3(Rot_gnss_init.transpose() * rot_window[i]));
          initialEstimate.insert(A(i), gtsam::Vector6(init_vel_bias_vector));
          // initialEstimate.insert(F(i), gtsam::Vector12(init_vel_bias_vector));
          initialEstimate.insert(B(i), gtsam::Vector4(para_rcv_dt[i*4], para_rcv_dt[i*4+1], para_rcv_dt[i*4+2], para_rcv_dt[i*4+3]));
          initialEstimate.insert(C(i), gtsam::Vector1(para_rcv_ddt[i]));
        }
        id_accumulate = 3 * wind_size + 5;
      }  
    }
    else
    {   
      Eigen::Matrix3d R_enu_local_;
      R_enu_local_ = Eigen::AngleAxisd(yaw_enu_local, Eigen::Vector3d::UnitZ());   
        for (size_t i = 0; i < wind_size+1; i++)
        {
          gtsam::PriorFactor<gtsam::Rot3> init_rot(R(i), gtsam::Rot3(R_ecef_enu * R_enu_local_ * rot_window[i]), priorrotNoise);
          gtsam::PriorFactor<gtsam::Vector4> init_dt(B(i), gtsam::Vector4(para_rcv_dt[i*4], para_rcv_dt[i*4+1], para_rcv_dt[i*4+2], para_rcv_dt[i*4+3]), priordtNoise);
          gtSAMgraph.add(init_rot);
          gtSAMgraph.add(init_dt);
          factor_id_frame.push_back(std::vector<size_t>{i * 3, i * 3 + 1, i * 3 + 2});
          // Eigen::Vector3d init_vel = ;
          Eigen::Matrix<double, 12, 1> init_vel_bias_vector;
          init_vel_bias_vector.block<3,1>(0,0) = anc_ecef + R_ecef_enu * R_enu_local_ * (pos_window[i] - pos_window[0] - rot_window[i] * Tex_imu_r);
          init_vel_bias_vector.block<3,1>(3,0) = R_ecef_enu * R_enu_local_ * vel_window[i];
          init_vel_bias_vector.block<6,1>(6,0) = Eigen::Matrix<double, 6, 1>::Zero();
          gtsam::PriorFactor<gtsam::Vector12> init_pos(F(i), gtsam::Vector12(init_vel_bias_vector), priorposNoise);
          gtSAMgraph.add(init_pos);
          
          initialEstimate.insert(R(i), gtsam::Rot3(R_ecef_enu * R_enu_local_ * rot_window[i]));
          initialEstimate.insert(F(i), gtsam::Vector12(init_vel_bias_vector));
          initialEstimate.insert(B(i), gtsam::Vector4(para_rcv_dt[i*4], para_rcv_dt[i*4+1], para_rcv_dt[i*4+2], para_rcv_dt[i*4+3]));
          initialEstimate.insert(C(i), gtsam::Vector1(para_rcv_ddt[i]));
        }
        id_accumulate = 3 * wind_size + 3;
    }
    frame_num = frame_count;
    last_gnss_time = time2sec(gnss_meas_buf[wind_size][0]->time);

    for (uint32_t k_ = 2; k_ < wind_size+1; k_++)
    {
      std::vector<ObsPtr>().swap(gnss_meas_buf[k_]);
      std::vector<EphemBasePtr>().swap(gnss_ephem_buf[k_]);
    }
  }

  for (uint32_t k_ = 1; k_ < wind_size+1; k_++)
  {
    std::vector<ObsPtr>().swap(gnss_meas_buf[k_]);
    std::vector<EphemBasePtr>().swap(gnss_ephem_buf[k_]);
  }
  runISAM2opt();
  return true;
}

void GNSSProcess::updateGNSSStatistics(Eigen::Vector3d &pos) // delete
{
  if (!nolidar)
  {
    Eigen::Vector3d anc_cur = isamCurrentEstimate.at<gtsam::Vector3>(E(0));
    Eigen::Matrix3d R_enu_local_ = isamCurrentEstimate.at<gtsam::Rot3>(P(0)).matrix();

    Eigen::Vector3d enu_pos = R_enu_local_ * (pos - anc_local);
    // R_ecef_enu = ecef2rotation(anc_cur);
    // ecef_pos = anc_cur + R_ecef_enu * enu_pos;
    ecef_pos = anc_cur + enu_pos;
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
  // relative_sqrt_info = 10;
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

bool GNSSProcess::Evaluate(state_input &state, Eigen::Vector3d &omg)
{
  if (gnss_meas_buf[0].empty() || gnss_meas_buf[0].size() < 4)
  {
    // cout << "no valid gnss" << endl;
    return false;
  }

  double rcv_dt[4];
  bool rcv_sys[4];
  rcv_sys[0] = false; rcv_sys[1] = false; rcv_sys[2] = false; rcv_sys[3] = false;
  double rcv_ddt;
  double time_current = time2sec(gnss_meas_buf[0][0]->time);
  invalid_lidar = false;
  if (!nolidar)
  {
    invalid_lidar = nolidar_cur;
    // for (size_t i = 0; i < 15; i++)
    // {
    //   if (state.cov(i, i) > 1.0) invalid_lidar = true;
    // }
  }
  if ((time_current - last_gnss_time > 15 * gnss_sample_period && nolidar) || (time_current - last_gnss_time > 15 * gnss_sample_period && invalid_lidar && !nolidar))
  {
    Reset();
    return false;
  }

  if (nolidar_cur)
  {
    nolidar_cur = false;
  }
  // cout << "check time diff:" << time_current - last_gnss_time << endl;
  {
    rcv_ddt = isamCurrentEstimate.at<gtsam::Vector1>(C(frame_num-1))[0];
    rcv_dt[0] = isamCurrentEstimate.at<gtsam::Vector4>(B(frame_num-1))[0] + rcv_ddt * (time_current - last_gnss_time);
    rcv_dt[1] = isamCurrentEstimate.at<gtsam::Vector4>(B(frame_num-1))[1] + rcv_ddt * (time_current - last_gnss_time);
    rcv_dt[2] = isamCurrentEstimate.at<gtsam::Vector4>(B(frame_num-1))[2] + rcv_ddt * (time_current - last_gnss_time);
    rcv_dt[3] = isamCurrentEstimate.at<gtsam::Vector4>(B(frame_num-1))[3] + rcv_ddt * (time_current - last_gnss_time);
  }

  const std::vector<ObsPtr> &curr_obs = gnss_meas_buf[0];
  const std::vector<EphemBasePtr> &curr_ephem = gnss_ephem_buf[0];

  // find best sat in the current gnss measurements
  std::map<std::pair<double, int>, std::map<uint32_t, double[6]> >::iterator it;
  double time_list_[curr_obs.size()];
  double time_min = time2sec(curr_obs[0]->time);
  for (uint32_t j = 0; j < curr_obs.size(); j++)
  {
    time_list_[j] = 0.0;
    for (it = sat2cp.begin(); it != sat2cp.end(); it++)
    {
      std::map<uint32_t, double[6]>::iterator it_sat;
      it_sat = it->second.find(curr_obs[j]->sat);
      if (it_sat != it->second.end())
      {
        time_list_[j] = it->first.first;
        if (time_list_[j] < time_min)
        {
          time_min = time_list_[j];
        }
        break;
      }
    }
  }

  int best_sat = -1;
  double min_cp_std = 100000000;
  for (uint32_t j = 0; j < curr_obs.size(); j++)
  {
    if (time_list_[j] == time_min)
    {
      freq = L1_freq(curr_obs[j], &freq_idx);
      LOG_IF(FATAL, freq < 0) << "No L1 observation found."; 
      if (curr_obs[j]->cp_std[freq_idx] * 0.004 < min_cp_std && curr_obs[j]->cp_std[freq_idx] * 0.004 < gnss_cp_std_threshold && curr_obs[j]->cp[freq_idx] > 10)
      {
        min_cp_std = curr_obs[j]->cp_std[freq_idx] * 0.004;
        best_sat = j;
      }
    }
  }

  std::deque<uint32_t> pair_sat_copy;
  std::deque<uint32_t>().swap(pair_sat_copy);

  std::deque<double> meas_sats;
  std::deque<double> meas_cov_sats;
  std::deque<double> meas_time_sats;
  std::deque<int> meas_index_sats;
  std::deque<Eigen::Vector3d> meas_svpos_sats;
  std::deque<Eigen::Vector3d> meas_svpos_best;
  
  // find pair sats to the best one in the current gnss measurements
  if (best_sat > -1)
  {
    std::deque<uint32_t> pair_sat;
    std::deque<uint32_t>().swap(pair_sat);

    for (uint32_t j = 0; j < curr_obs.size(); j++)
    {
      if (j != best_sat)
      {
        uint32_t sys = satsys(curr_obs[j]->sat, NULL);
        freq = L1_freq(curr_obs[j], &freq_idx);

        if (sys == satsys(curr_obs[best_sat]->sat, NULL) && curr_obs[j]->cp_std[freq_idx] * 0.004 < gnss_cp_std_threshold)
        {
          if (curr_obs[j]->cp[freq_idx] > 10)
          {
            pair_sat.push_back(j);
          }
        }
      }
    }

    for (uint32_t j = 0; j < curr_obs.size(); j++)
    {
      std::map<uint32_t, double[6]>::iterator it_old_best, it_old;
      double meas;
      double meas_cov;
      double meas_time;
      int meas_index;
      Eigen::Vector3d meas_svpos;
      Eigen::Vector3d best_svpos;
      if (pair_sat.size() > 0)
      {
        if (j == pair_sat.front())
        {
          bool cp_found = false;
          for (it = sat2cp.begin(); it != sat2cp.end(); it++)
          {
            it_old = it->second.find(curr_obs[j]->sat);
            it_old_best = it->second.find(curr_obs[best_sat]->sat);
            if (it_old != it->second.end() && it_old_best != it->second.end())
            {
              cp_found = true;
              meas = it_old_best->second[0] - it_old->second[0]; // - it_old_best->second[1] + it_old->second[1]);
              meas_cov = (it_old_best->second[2] * it_old_best->second[2] + it_old->second[2] * it_old->second[2]); // 
              meas_time = it->first.first;
              meas_index = it->first.second;
              meas_svpos << it_old->second[3], it_old->second[4], it_old->second[5];
              best_svpos << it_old_best->second[3], it_old_best->second[4], it_old_best->second[5];
            }
          }
        
          if (cp_found)
          {
            meas_sats.push_back(meas);
            meas_cov_sats.push_back(meas_cov);
            meas_time_sats.push_back(meas_time);
            meas_index_sats.push_back(meas_index);
            meas_svpos_sats.push_back(meas_svpos);
            meas_svpos_best.push_back(best_svpos);
            pair_sat_copy.push_back(j);
          }
          pair_sat.pop_front();
        }
      }
    }
    
  }

  std::map<uint32_t, double[6]> curr_cp_map;

  std::vector<double> meas_cp;
  std::vector<double> cov_cp;
  std::vector<Eigen::Vector3d> sv_pos_pair;
  double cov_cp_best, meas_cp_best; //, esti_cp_best, 
  Eigen::Vector3d sv_pos_best;
  std::vector<size_t> factor_id_cur;

  // if (gnss_online_init || !nolidar || frame_num >= 10)
  M3D omg_skew;

  omg_skew << SKEW_SYM_MATRX(omg);
  Eigen::Vector3d hat_omg_T = omg_skew * Tex_imu_r;
  {   
  for (uint32_t j = 0; j < curr_obs.size(); j++)
  {
    const uint32_t sys = satsys(curr_obs[j]->sat, NULL);
    const uint32_t sys_idx = gnss_comm::sys2idx.at(sys);
    GnssPsrDoppMeas(curr_obs[j], curr_ephem[j]); //, latest_gnss_iono_params);
    

    freq = L1_freq(curr_obs[j], &freq_idx); // save
    
    const double wavelength = LIGHT_SPEED / freq; // save
    if (curr_obs[j]->cp_std[freq_idx] * 0.004 < gnss_cp_std_threshold)
    {
      if (curr_obs[j]->cp[freq_idx] > 10)
      {
        curr_cp_map[curr_obs[j]->sat][0] = curr_obs[j]->cp[freq_idx] * wavelength;
        curr_cp_map[curr_obs[j]->sat][2] = curr_obs[j]->cp_std[freq_idx] * 0.004;
        curr_cp_map[curr_obs[j]->sat][3] = sv_pos[0];
        curr_cp_map[curr_obs[j]->sat][4] = sv_pos[1];
        curr_cp_map[curr_obs[j]->sat][5] = sv_pos[2];
      }
    }

    if (j == best_sat)
    {
      meas_cp_best = curr_obs[j]->cp[freq_idx] * wavelength;
      sv_pos_best = sv_pos;
      cov_cp_best  = curr_obs[j]->cp_std[freq_idx] * 0.004 * curr_obs[j]->cp_std[freq_idx] * 0.004;
    }

    if (pair_sat_copy.size() > 0)
    {
    if (j == pair_sat_copy.front())
    {
      meas_cp.push_back(curr_obs[j]->cp[freq_idx] * wavelength);
      cov_cp.push_back(curr_obs[j]->cp_std[freq_idx] * curr_obs[j]->cp_std[freq_idx] * 0.004 * 0.004);
      sv_pos_pair.push_back(sv_pos);
      pair_sat_copy.pop_front();
    }
    }
    /////////////////////////////////
    double values[30];
    values[0] = Tex_imu_r[0]; values[1] = Tex_imu_r[1]; values[2] = Tex_imu_r[2]; values[3] = anc_local[0]; values[4] = anc_local[1]; values[5] = anc_local[2];
    values[6] = sv_pos[0]; values[7] = sv_pos[1]; values[8] = sv_pos[2]; values[9] = sv_vel[0]; values[10] = sv_vel[1]; values[11] = sv_vel[2];
    values[12] = svdt; values[13] = tgd; values[14] = svddt; values[15] = pr_uura; values[16] = dp_uura; values[17] = relative_sqrt_info;
    values[18] = latest_gnss_iono_params[0]; values[19] = latest_gnss_iono_params[1]; values[20] = latest_gnss_iono_params[2]; values[21] = latest_gnss_iono_params[3]; values[22] = latest_gnss_iono_params[4]; values[23] = latest_gnss_iono_params[5];
    values[24] = latest_gnss_iono_params[6]; values[25] = latest_gnss_iono_params[7]; values[26] = time_current; values[27] = freq; values[28] = curr_obs[j]->psr[freq_idx]; values[29] = curr_obs[j]->dopp[freq_idx];

    {
    rcv_sys[sys_idx] = true;
    if (!nolidar)
    {   
      gtSAMgraph.add(glio::GnssPsrDoppFactor(R(frame_num), A(frame_num), B(frame_num), C(frame_num), E(0), P(0), values, sys_idx, hat_omg_T, robustpsrdoppNoise));
    }
    else
    {    
      gtSAMgraph.add(glio::GnssPsrDoppFactorNolidar(R(frame_num), F(frame_num), B(frame_num), C(frame_num), values, sys_idx, hat_omg_T, robustpsrdoppNoise)); // not work
    }
    {
      factor_id_cur.push_back(id_accumulate);
    }
    id_accumulate += 1;
    }
  }
  }

  Eigen::Matrix<double, 12, 1> init_vel_bias_vector;
  init_vel_bias_vector.block<3,1>(0,0) = state.pos;
  init_vel_bias_vector.block<3,1>(3,0) = state.vel;
  init_vel_bias_vector.block<3,1>(6,0) = state.ba;
  init_vel_bias_vector.block<3,1>(9,0) = state.bg;

  // if (gnss_online_init || !nolidar || frame_num >= 10)
  {
  gtSAMgraph.add(glio::DdtSmoothFactor(C(frame_num-1), C(frame_num), ddtNoise));
  gtSAMgraph.add(glio::DtDdtFactor(B(frame_num-1), B(frame_num), C(frame_num-1), C(frame_num), rcv_sys, time_current - last_gnss_time, dtNoise)); // not work
  factor_id_frame[frame_num-1-frame_delete].push_back(id_accumulate);
  factor_id_frame[frame_num-1-frame_delete].push_back(id_accumulate+1);
  id_accumulate += 2;


  if (!nolidar) // && !invalid_lidar)
  {
    // Eigen::Matrix<double, 15, 15> state_cov = Eigen::Matrix<double, 15, 15>::Identity();
    // state_cov.block<9, 9>(0, 0) = state.cov.block<9, 9>(0, 0) * 100;
    // gtsam::noiseModel::Gaussian::shared_ptr LioNoise = gtsam::noiseModel::Gaussian::Covariance(state_cov);
    Eigen::Matrix3d last_rot = state_last.rot.normalized().toRotationMatrix(); // isamCurrentEstimate.at<gtsam::Rot3>(R(frame_num-1)).matrix(); // state_const_.rot; // 
    // cout << "check time period" << pre_integration->sum_dt << ";" << time_current - last_gnss_time <<  endl;
    Eigen::Vector3d last_pos = state_last.pos; // isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(0); // state_.pos; // 
    Eigen::Vector3d last_vel = state_last.vel; // isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(3); // state_.vel; //
    state_last = state;
    double delta_t = time_current - last_gnss_time;
    // Eigen::Vector3d cur_grav = state.rot.transpose() * state.gravity; //  
    gtsam::Rot3 rel_rot_(last_rot.transpose() * state.rot.normalized().toRotationMatrix());
    gtsam::Point3 rel_pos_(last_rot.transpose() * (state.pos - last_pos - last_vel * delta_t - 0.5 * state.gravity * delta_t * delta_t)); 
    // (state.pos - last_pos);
    gtsam::Vector3 rel_v_(last_rot.transpose() * (state.vel - last_vel - state.gravity * delta_t)); // (state.vel - last_vel);
    gtSAMgraph.add(glio::GnssLioFactor(R(frame_num-1), A(frame_num-1), R(frame_num), A(frame_num), 1.0, rel_rot_, rel_pos_, rel_v_,
                  state.gravity, delta_t, margposNoise)); // odomNoiseIMU)); // LioNoise)); // 
    factor_id_frame[frame_num-1-frame_delete].push_back(id_accumulate);
    id_accumulate += 1;
    // gtSAMgraph.add(glio::GnssLioGravFactor(R(frame_num), cur_grav, state.gravity, odomNoise));
    // Eigen::Matrix<double, 15, 15> state_cov = Eigen::Matrix<double, 15, 15>::Identity(); 
    // state_cov.block<9, 9>(0, 0) = state.cov.block<9, 9>(0, 0) * 1000;
    // gtsam::noiseModel::Gaussian::shared_ptr LioNoise = gtsam::noiseModel::Gaussian::Covariance(state_cov);
    // gtSAMgraph.add(glio::GnssLioHardFactor(R(frame_num), F(frame_num), init_vel_bias_vector, state.rot, LioNoise)); // odomNoiseIMU));
    // factor_id_cur.push_back(id_accumulate);
    // id_accumulate += 1;
    // gtsam::Rot3 rel_rot(pre_integration->delta_q); // (last_rot.transpose() * state.rot);
    // // double delta_t = time_current - last_gnss_time;
    // gtsam::Point3 rel_pos(pre_integration->delta_p);
    // gtsam::Vector3 rel_v(pre_integration->delta_v); 
    // gtSAMgraph.add(glio::GnssLioFactorNolidar(R(frame_num-1), F(frame_num-1), R(frame_num), F(frame_num), rel_rot, rel_pos, rel_v, 
    //               state_.gravity, time_current - last_gnss_time, state_.ba, state_.bg, pre_integration, odomNoiseIMU));
    // factor_id_frame[frame_num-1-frame_delete].push_back(id_accumulate);
    // id_accumulate += 1;
  }
  else
  {
    // if (!nolidar)
    // {
    //   state = state_;
    // }
    gtsam::Rot3 rel_rot(pre_integration->delta_q); // (last_rot.transpose() * state.rot);
    // double delta_t = time_current - last_gnss_time;
    gtsam::Point3 rel_pos(pre_integration->delta_p);
    gtsam::Vector3 rel_v(pre_integration->delta_v); 
    gtSAMgraph.add(glio::GnssLioFactorNolidar(R(frame_num-1), F(frame_num-1), R(frame_num), F(frame_num), rel_rot, rel_pos, rel_v, 
                  state.gravity, time_current - last_gnss_time, state.ba, state.bg, pre_integration, odomNoiseIMU));

    factor_id_frame[frame_num-1-frame_delete].push_back(id_accumulate);
  
    id_accumulate += 1;
  }
  }

  if (!nolidar) // && invalid_lidar)
  {
    Eigen::Matrix<double, 6, 1> init_vel_bias_vector_imu;
    init_vel_bias_vector_imu.block<3,1>(0,0) = state.pos;
    init_vel_bias_vector_imu.block<3,1>(3,0) = state.vel;
    // init_vel_bias_vector_imu.block<3,1>(6,0) = state_.ba;
    // init_vel_bias_vector_imu.block<3,1>(9,0) = state_.bg;
    initialEstimate.insert(A(frame_num), gtsam::Vector6(init_vel_bias_vector_imu));
    initialEstimate.insert(R(frame_num), gtsam::Rot3(state.rot.normalized().toRotationMatrix()));  
  }
  else
  {
    initialEstimate.insert(F(frame_num), gtsam::Vector12(init_vel_bias_vector));
    initialEstimate.insert(R(frame_num), gtsam::Rot3(state.rot.normalized().toRotationMatrix())); 
  }                 
  initialEstimate.insert(C(frame_num), gtsam::Vector1(rcv_ddt));
  initialEstimate.insert(B(frame_num), gtsam::Vector4(rcv_dt[0], rcv_dt[1], rcv_dt[2], rcv_dt[3]));     
  
  // if (gnss_online_init || !nolidar || frame_num >= 10)
  {
  for (uint32_t j = 0; j < meas_index_sats.size(); j++)
  {
    double values[20];
    values[0] = Tex_imu_r[0]; values[1] = Tex_imu_r[1]; values[2] = Tex_imu_r[2]; values[3] = anc_local[0]; values[4] = anc_local[1]; values[5] = anc_local[2];
    values[6] = meas_svpos_best[j][0]; values[7] = meas_svpos_best[j][1]; values[8] = meas_svpos_best[j][2]; values[9] = meas_svpos_sats[j][0]; values[10] = meas_svpos_sats[j][1]; values[11] = meas_svpos_sats[j][2];
    values[12] = sv_pos_best[0]; values[13] = sv_pos_best[1]; values[14] = sv_pos_best[2]; values[15] = sv_pos_pair[j][0]; values[16] = sv_pos_pair[j][1]; values[17] = sv_pos_pair[j][2];
    values[18] = meas_cp_best - meas_cp[j] - meas_sats[j]; values[19] = cp_weight; 
    if (!nolidar)
    {
      {
        gtSAMgraph.add(glio::GnssCpFactor(E(0), P(0), R(meas_index_sats[j]), A(meas_index_sats[j]), R(frame_num), A(frame_num), values, robustcpNoise));
      }
    }
    else
    {
      gtSAMgraph.add(glio::GnssCpFactorNolidar(R(meas_index_sats[j]), F(meas_index_sats[j]), R(frame_num), F(frame_num), values, robustcpNoise)); // not work
    }
    // factor_id_cur.push_back(id_accumulate);
    factor_id_frame[meas_index_sats[j]-frame_delete].push_back(id_accumulate);
    id_accumulate += 1;
  }
  }
  // if (frame_num % delete_thred > 0) // || nolidar)
  {
    factor_id_frame.push_back(factor_id_cur);
    std::vector<size_t>().swap(factor_id_cur);
  }
  frame_num ++;
  // if (gnss_online_init || !nolidar || frame_num >= 10)
  {
    runISAM2opt();

    // state.cov.block<3,3>(0, 0) = isam.marginalCovariance(R(frame_num-1));
    // state.cov.block<6,6>(3, 3) = isam.marginalCovariance(F(frame_num-1)).block<6, 6>(0, 0);
    
    if (nolidar)
    {
      state.rot = Eigen::Quaterniond(isamCurrentEstimate.at<gtsam::Rot3>(R(frame_num-1)).matrix());
      cout << "check norm:" << state.rot.norm() << ";" << isamCurrentEstimate.at<gtsam::Rot3>(R(frame_num-1)).matrix().norm() << endl;
      state.pos = isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(0);
      state.vel = isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(3);
      state.ba = isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(6);
      state.bg = isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(9);
      state.gravity = ecef2rotation(state.pos) * gravity_init;
    }
    else
    {
      state_.rot = isamCurrentEstimate.at<gtsam::Rot3>(R(frame_num-1)).matrix();
      state_.pos = isamCurrentEstimate.at<gtsam::Vector6>(A(frame_num-1)).segment<3>(0);
      state_.vel = isamCurrentEstimate.at<gtsam::Vector6>(A(frame_num-1)).segment<3>(3);
      // state_.ba = isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(6);
      // state_.bg = isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(9);
      // state.gravity = isamCurrentEstimate.at<gtsam::Rot3>(P(0)).matrix().transpose() * gravity_init;
      state_.gravity = isamCurrentEstimate.at<gtsam::Rot3>(P(0)).matrix().transpose() * ecef2rotation(isamCurrentEstimate.at<gtsam::Vector3>(E(0))) * gravity_init;
      state.gravity = state_.gravity; // 
      // state.rot = isamCurrentEstimate.at<gtsam::Rot3>(R(frame_num-1)).matrix();
      // state.pos = isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(0);
      // state.vel = isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(3);
      // state.ba = isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(6);
      // state.bg = isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(9);
      // state.gravity = state_.gravity; // ecef2rotation(state.pos) * gravity_init;
    }    
    // state_ = state;
    if (!nolidar && invalid_lidar)
    {
      state = state_;
      // state.cov = MD(15, 15)::Identity() * 0.01;
      printf("invalid lidar!!!!!!!");
      // state.cov.setIdentity();
      // state.cov.block<3, 3>(0, 0) = isam.marginalCovariance(R(frame_num-1));
      // state.cov.block<6, 6>(3, 3) = isam.marginalCovariance(F(frame_num-1)).block<6, 6>(0, 0);
      // state.cov.block<3, 3>(9, 9) = isam.marginalCovariance(F(frame_num-1)).block<3, 3>(9, 9);
      // state.cov.block<3, 3>(12, 12) = isam.marginalCovariance(F(frame_num-1)).block<3, 3>(6, 6);
      // invalid_lidar = false;
    }
    last_gnss_time = time_current;
    // state_ = state;
  }
  sat2cp[std::pair<double, int>(time2sec(curr_obs[0]->time), frame_num-1)] = curr_cp_map;
  std::map<std::pair<double, int>, std::map<uint32_t, double[6]>>::iterator it_old;
  if (!sat2cp.empty())
  {
    it_old = sat2cp.begin();
    while (time2sec(curr_obs[0]->time) - it_old->first.first > gnss_cp_time_threshold)
    {
      std::map<uint32_t, double[6]>().swap(it_old->second);
      size_t del_size = sat2cp.erase(it_old->first);
      if (sat2cp.empty()) break;
      it_old = sat2cp.begin();
    }
  }
  return true;
}

bool GNSSProcess::Evaluate(state_output &state)
{
  if (gnss_meas_buf[0].empty() || gnss_meas_buf[0].size() < 4)
  {
    // cout << "no valid gnss" << endl;
    return false;
  }

  double rcv_dt[4];
  double rcv_ddt;
  bool rcv_sys[4];
  rcv_sys[0] = false; rcv_sys[1] = false; rcv_sys[2] = false; rcv_sys[3] = false;
  double time_current = time2sec(gnss_meas_buf[0][0]->time);
  invalid_lidar = false;
  if (!nolidar)
  {
    invalid_lidar = nolidar_cur;
    // for (size_t i = 0; i < 21; i++)
    // {
    //   if (state.cov(i, i) > 10.0) invalid_lidar = true;
    // }
  }
  if ((time_current - last_gnss_time > 15 * gnss_sample_period && nolidar) || (time_current - last_gnss_time > 15 * gnss_sample_period && invalid_lidar && !nolidar))
  {
    Reset();
    return false;
  }

  if (nolidar_cur)
  {
    nolidar_cur = false;
  }
  // cout << "check time diff:" << time_current - last_gnss_time << endl;
  {
    rcv_ddt = isamCurrentEstimate.at<gtsam::Vector1>(C(frame_num-1))[0]; // modify!
    rcv_dt[0] = isamCurrentEstimate.at<gtsam::Vector4>(B(frame_num-1))[0] + rcv_ddt * (time_current - last_gnss_time);
    rcv_dt[1] = isamCurrentEstimate.at<gtsam::Vector4>(B(frame_num-1))[1] + rcv_ddt * (time_current - last_gnss_time);
    rcv_dt[2] = isamCurrentEstimate.at<gtsam::Vector4>(B(frame_num-1))[2] + rcv_ddt * (time_current - last_gnss_time);
    rcv_dt[3] = isamCurrentEstimate.at<gtsam::Vector4>(B(frame_num-1))[3] + rcv_ddt * (time_current - last_gnss_time);
  }

  const std::vector<ObsPtr> &curr_obs = gnss_meas_buf[0];
  const std::vector<EphemBasePtr> &curr_ephem = gnss_ephem_buf[0];

  // find best sat in the current gnss measurements
  std::map<std::pair<double, int>, std::map<uint32_t, double[6]> >::iterator it;
  double time_list_[curr_obs.size()];
  double time_min = time2sec(curr_obs[0]->time);
  for (uint32_t j = 0; j < curr_obs.size(); j++)
  {
    time_list_[j] = 0.0;
    for (it = sat2cp.begin(); it != sat2cp.end(); it++)
    {
      std::map<uint32_t, double[6]>::iterator it_sat;
      it_sat = it->second.find(curr_obs[j]->sat);
      if (it_sat != it->second.end())
      {
        time_list_[j] = it->first.first;
        if (time_list_[j] < time_min)
        {
          time_min = time_list_[j];
        }
        break;
      }
    }
  }

  int best_sat = -1;
  double min_cp_std = 100000000;
  for (uint32_t j = 0; j < curr_obs.size(); j++)
  {
    if (time_list_[j] == time_min)
    {
      freq = L1_freq(curr_obs[j], &freq_idx);
      LOG_IF(FATAL, freq < 0) << "No L1 observation found."; 
      if (curr_obs[j]->cp_std[freq_idx] * 0.004 < min_cp_std && curr_obs[j]->cp_std[freq_idx] * 0.004 < gnss_cp_std_threshold && curr_obs[j]->cp[freq_idx] > 10)
      {
        min_cp_std = curr_obs[j]->cp_std[freq_idx] * 0.004;
        best_sat = j;
      }
    }
  }

  std::deque<uint32_t> pair_sat_copy;
  std::deque<uint32_t>().swap(pair_sat_copy);

  std::deque<double> meas_sats;
  std::deque<double> meas_cov_sats;
  std::deque<double> meas_time_sats;
  std::deque<int> meas_index_sats;
  std::deque<Eigen::Vector3d> meas_svpos_sats;
  std::deque<Eigen::Vector3d> meas_svpos_best;
  
  // find pair sats to the best one in the current gnss measurements
  if (best_sat > -1)
  {
    std::deque<uint32_t> pair_sat;
    std::deque<uint32_t>().swap(pair_sat);

    for (uint32_t j = 0; j < curr_obs.size(); j++)
    {
      if (j != best_sat)
      {
        uint32_t sys = satsys(curr_obs[j]->sat, NULL);
        freq = L1_freq(curr_obs[j], &freq_idx);

        if (sys == satsys(curr_obs[best_sat]->sat, NULL) && curr_obs[j]->cp_std[freq_idx] * 0.004 < gnss_cp_std_threshold)
        {
          if (curr_obs[j]->cp[freq_idx] > 10)
          {
            pair_sat.push_back(j);
          }
        }
      }
    }

    for (uint32_t j = 0; j < curr_obs.size(); j++)
    {
      std::map<uint32_t, double[6]>::iterator it_old_best, it_old;
      double meas;
      double meas_cov;
      double meas_time;
      int meas_index;
      Eigen::Vector3d meas_svpos;
      Eigen::Vector3d best_svpos;
      if (pair_sat.size() > 0)
      {
        if (j == pair_sat.front())
        {
          bool cp_found = false;
          for (it = sat2cp.begin(); it != sat2cp.end(); it++)
          {
            it_old = it->second.find(curr_obs[j]->sat);
            it_old_best = it->second.find(curr_obs[best_sat]->sat);
            if (it_old != it->second.end() && it_old_best != it->second.end())
            {
              cp_found = true;
              meas = it_old_best->second[0] - it_old->second[0]; // - it_old_best->second[1] + it_old->second[1]);
              meas_cov = (it_old_best->second[2] * it_old_best->second[2] + it_old->second[2] * it_old->second[2]); // 
              meas_time = it->first.first;
              meas_index = it->first.second;
              meas_svpos << it_old->second[3], it_old->second[4], it_old->second[5];
              best_svpos << it_old_best->second[3], it_old_best->second[4], it_old_best->second[5];
            }
          }
        
          if (cp_found)
          {
            meas_sats.push_back(meas);
            meas_cov_sats.push_back(meas_cov);
            meas_time_sats.push_back(meas_time);
            meas_index_sats.push_back(meas_index);
            meas_svpos_sats.push_back(meas_svpos);
            meas_svpos_best.push_back(best_svpos);
            pair_sat_copy.push_back(j);
          }
          pair_sat.pop_front();
        }
      }
    }
    
  }

  std::map<uint32_t, double[6]> curr_cp_map;

  std::vector<double> meas_cp;
  std::vector<double> cov_cp;
  std::vector<Eigen::Vector3d> sv_pos_pair;
  double cov_cp_best, meas_cp_best; //, esti_cp_best, 
  Eigen::Vector3d sv_pos_best;
  std::vector<size_t> factor_id_cur;

    M3D omg_skew;

    omg_skew << SKEW_SYM_MATRX(state.omg);
    Eigen::Vector3d hat_omg_T = omg_skew * Tex_imu_r;;
    
  // if (gnss_online_init || !nolidar || frame_num >= 10)
  {
  for (uint32_t j = 0; j < curr_obs.size(); j++)
  {
    const uint32_t sys = satsys(curr_obs[j]->sat, NULL);
    const uint32_t sys_idx = gnss_comm::sys2idx.at(sys);
    GnssPsrDoppMeas(curr_obs[j], curr_ephem[j]); //, latest_gnss_iono_params);
    

    freq = L1_freq(curr_obs[j], &freq_idx); // save
    
    const double wavelength = LIGHT_SPEED / freq; // save
    if (curr_obs[j]->cp_std[freq_idx] * 0.004 < gnss_cp_std_threshold)
    {
      if (curr_obs[j]->cp[freq_idx] > 10)
      {
        curr_cp_map[curr_obs[j]->sat][0] = curr_obs[j]->cp[freq_idx] * wavelength;
        curr_cp_map[curr_obs[j]->sat][2] = curr_obs[j]->cp_std[freq_idx] * 0.004;
        curr_cp_map[curr_obs[j]->sat][3] = sv_pos[0];
        curr_cp_map[curr_obs[j]->sat][4] = sv_pos[1];
        curr_cp_map[curr_obs[j]->sat][5] = sv_pos[2];
      }
    }

    if (j == best_sat)
    {
      meas_cp_best = curr_obs[j]->cp[freq_idx] * wavelength;
      sv_pos_best = sv_pos;
      cov_cp_best  = curr_obs[j]->cp_std[freq_idx] * 0.004 * curr_obs[j]->cp_std[freq_idx] * 0.004;
    }

    if (pair_sat_copy.size() > 0)
    {
    if (j == pair_sat_copy.front())
    {
      meas_cp.push_back(curr_obs[j]->cp[freq_idx] * wavelength);
      cov_cp.push_back(curr_obs[j]->cp_std[freq_idx] * curr_obs[j]->cp_std[freq_idx] * 0.004 * 0.004);
      sv_pos_pair.push_back(sv_pos);
      pair_sat_copy.pop_front();
    }
    }
    /////////////////////////////////
    double values[30];
    values[0] = Tex_imu_r[0]; values[1] = Tex_imu_r[1]; values[2] = Tex_imu_r[2]; values[3] = anc_local[0]; values[4] = anc_local[1]; values[5] = anc_local[2];
    values[6] = sv_pos[0]; values[7] = sv_pos[1]; values[8] = sv_pos[2]; values[9] = sv_vel[0]; values[10] = sv_vel[1]; values[11] = sv_vel[2];
    values[12] = svdt; values[13] = tgd; values[14] = svddt; values[15] = pr_uura; values[16] = dp_uura; values[17] = relative_sqrt_info;
    values[18] = latest_gnss_iono_params[0]; values[19] = latest_gnss_iono_params[1]; values[20] = latest_gnss_iono_params[2]; values[21] = latest_gnss_iono_params[3]; values[22] = latest_gnss_iono_params[4]; values[23] = latest_gnss_iono_params[5];
    values[24] = latest_gnss_iono_params[6]; values[25] = latest_gnss_iono_params[7]; values[26] = time_current; values[27] = freq; values[28] = curr_obs[j]->psr[freq_idx]; values[29] = curr_obs[j]->dopp[freq_idx];

    {
    rcv_sys[sys_idx] = true;
    if (!nolidar)
    {    
      gtSAMgraph.add(glio::GnssPsrDoppFactor(R(frame_num), A(frame_num), B(frame_num), C(frame_num), E(0), P(0), values, sys_idx, hat_omg_T, robustpsrdoppNoise));
    }
    else
    {    
      gtSAMgraph.add(glio::GnssPsrDoppFactorNolidar(R(frame_num), A(frame_num), B(frame_num), C(frame_num), values, sys_idx, hat_omg_T, robustpsrdoppNoise)); // not work
    }
    {
      factor_id_cur.push_back(id_accumulate);
    }
    id_accumulate += 1;
    }
  }
  }

  Eigen::Matrix<double, 12, 1> init_vel_bias_vector;
  init_vel_bias_vector.block<3,1>(0,0) = state.pos;
  init_vel_bias_vector.block<3,1>(3,0) = state.vel;
  init_vel_bias_vector.block<3,1>(6,0) = state.ba;
  init_vel_bias_vector.block<3,1>(9,0) = state.bg;

  // if (gnss_online_init || !nolidar || frame_num >= 10)
  {
  gtSAMgraph.add(glio::DdtSmoothFactor(C(frame_num-1), C(frame_num), ddtNoise));
  gtSAMgraph.add(glio::DtDdtFactor(B(frame_num-1), B(frame_num), C(frame_num-1), C(frame_num), rcv_sys, time_current - last_gnss_time, dtNoise)); // not work
  factor_id_frame[frame_num-1-frame_delete].push_back(id_accumulate);
  factor_id_frame[frame_num-1-frame_delete].push_back(id_accumulate+1);
  id_accumulate += 2;

  if (!nolidar) // && !invalid_lidar)
  {
    Eigen::Matrix3d last_rot = state_const_last.rot.normalized().toRotationMatrix(); // isamCurrentEstimate.at<gtsam::Rot3>(R(frame_num-1)); // state_const_.rot; // 
    Eigen::Vector3d last_pos = state_const_last.pos; // isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(0); // state_const_.pos; // 
    Eigen::Vector3d last_vel = state_const_last.vel; // isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(3); // state_const_.vel; // 
    state_const_last = state;
    double delta_t = time_current - last_gnss_time;
    gtsam::Rot3 rel_rot(last_rot.transpose() * state.rot.normalized().toRotationMatrix());
    gtsam::Point3 rel_pos(last_rot.transpose() * (state.pos - last_pos)); // - last_vel * delta_t - 0.5 * state.gravity * delta_t * delta_t)); 
    // Eigen::Matrix<double, 9, 1> init_vel_bias_vector;
    // init_vel_bias_vector.block<3,1>(0,0) = state.vel;
    // init_vel_bias_vector.block<3,1>(3,0) = state.ba;
    // init_vel_bias_vector.block<3,1>(6,0) = state.bg;
    Eigen::Vector3d cur_grav = state.rot.normalized().toRotationMatrix().transpose() * state.gravity;
    gtsam::Vector3 rel_v(last_rot.transpose() * (state.vel - last_vel)); // - state.gravity * delta_t)); // (state.vel - last_vel);
    // initialEstimate.insert(F(frame_num), gtsam::Vector9(init_vel_bias_vector));
    gtSAMgraph.add(glio::GnssLioFactor(R(frame_num-1), A(frame_num-1), R(frame_num), A(frame_num), 1.0, rel_rot, rel_pos, rel_v,
                state.gravity, delta_t, margposNoise)); // odomNoiseIMU));
    factor_id_frame[frame_num-1-frame_delete].push_back(id_accumulate);
    id_accumulate += 1;
    gtSAMgraph.add(glio::GnssLioGravFactor(R(frame_num), cur_grav, state.gravity, odomNoise));
    // Eigen::Matrix<double, 15, 15> state_cov = Eigen::Matrix<double, 15, 15>::Identity();
    // state_cov.block<9, 9>(0, 0) = state.cov.block<9, 9>(0, 0) * 100;
    // gtsam::noiseModel::Gaussian::shared_ptr LioNoise = gtsam::noiseModel::Gaussian::Covariance(state_cov);
    // gtSAMgraph.add(glio::GnssLioHardFactor(R(frame_num), F(frame_num), init_vel_bias_vector, state.rot, LioNoise)); // odomNoiseIMU));
    factor_id_cur.push_back(id_accumulate);
    id_accumulate += 1;
    // gtsam::Rot3 rel_rot(pre_integration->delta_q); //(last_rot.transpose() * state.rot); is able to be added;
    // gtsam::Point3 rel_pos(pre_integration->delta_p); //(state.pos - last_pos);
    // gtsam::Vector3 rel_v(pre_integration->delta_v); //(state.vel - last_vel);
    // gtSAMgraph.add(glio::GnssLioFactorNolidar(R(frame_num-1), F(frame_num-1), R(frame_num), F(frame_num), rel_rot, rel_pos, rel_v,
    //               state_const_.gravity, time_current - last_gnss_time, state_const_.ba, state_const_.bg, pre_integration, odomNoiseIMU));
    // factor_id_frame[frame_num-1-frame_delete].push_back(id_accumulate);
    // id_accumulate += 1;
  }
  else
  {
    // if (!nolidar)
    // {
    //   state = state_const_;
    // }
    gtsam::Rot3 rel_rot(pre_integration->delta_q); //(last_rot.transpose() * state.rot);
    gtsam::Point3 rel_pos(pre_integration->delta_p); //(state.pos - last_pos);
    gtsam::Vector3 rel_v(pre_integration->delta_v); //(state.vel - last_vel);
    gtSAMgraph.add(glio::GnssLioFactorNolidar(R(frame_num-1), F(frame_num-1), R(frame_num), F(frame_num), rel_rot, rel_pos, rel_v,
                  state.gravity, time_current - last_gnss_time, state.ba, state.bg, pre_integration, odomNoiseIMU));
    factor_id_frame[frame_num-1-frame_delete].push_back(id_accumulate);
    id_accumulate += 1;
  }
  }

  if (!nolidar) // && invalid_lidar)
  {
    Eigen::Matrix<double, 6, 1> init_vel_bias_vector_imu;
    init_vel_bias_vector_imu.block<3,1>(0,0) = state.pos;
    init_vel_bias_vector_imu.block<3,1>(3,0) = state.vel;
    // init_vel_bias_vector_imu.block<3,1>(6,0) = state_const_.bias_a;
    // init_vel_bias_vector_imu.block<3,1>(9,0) = state_const_.bias_g;
    initialEstimate.insert(A(frame_num), gtsam::Vector6(init_vel_bias_vector_imu));
    initialEstimate.insert(R(frame_num), gtsam::Rot3(state.rot.normalized().toRotationMatrix()));  
  }
  else
  {
    initialEstimate.insert(F(frame_num), gtsam::Vector12(init_vel_bias_vector));
    initialEstimate.insert(R(frame_num), gtsam::Rot3(state.rot.normalized().toRotationMatrix()));                  
  }
  initialEstimate.insert(C(frame_num), gtsam::Vector1(rcv_ddt));
  initialEstimate.insert(B(frame_num), gtsam::Vector4(rcv_dt[0], rcv_dt[1], rcv_dt[2], rcv_dt[3]));                

  // if (gnss_online_init || !nolidar || frame_num >= 10)
  {
  for (uint32_t j = 0; j < meas_index_sats.size(); j++)
  {
    double values[20];
    values[0] = Tex_imu_r[0]; values[1] = Tex_imu_r[1]; values[2] = Tex_imu_r[2]; values[3] = anc_local[0]; values[4] = anc_local[1]; values[5] = anc_local[2];
    values[6] = meas_svpos_best[j][0]; values[7] = meas_svpos_best[j][1]; values[8] = meas_svpos_best[j][2]; values[9] = meas_svpos_sats[j][0]; values[10] = meas_svpos_sats[j][1]; values[11] = meas_svpos_sats[j][2];
    values[12] = sv_pos_best[0]; values[13] = sv_pos_best[1]; values[14] = sv_pos_best[2]; values[15] = sv_pos_pair[j][0]; values[16] = sv_pos_pair[j][1]; values[17] = sv_pos_pair[j][2];
    values[18] = meas_cp_best - meas_cp[j] - meas_sats[j]; values[19] = cp_weight; 
    if (!nolidar)
    {
      {
        gtSAMgraph.add(glio::GnssCpFactor(E(0), P(0), R(meas_index_sats[j]), A(meas_index_sats[j]), R(frame_num), A(frame_num), values, robustcpNoise));
      }
    }
    else
    {
      gtSAMgraph.add(glio::GnssCpFactorNolidar(R(meas_index_sats[j]), F(meas_index_sats[j]), R(frame_num), F(frame_num), values, robustcpNoise)); // not work
    }
    // factor_id_cur.push_back(id_accumulate);
    factor_id_frame[meas_index_sats[j]-frame_delete].push_back(id_accumulate);
    id_accumulate += 1;
  }
  }
  // if (frame_num % delete_thred > 0) // || nolidar)
  {
    factor_id_frame.push_back(factor_id_cur);
    std::vector<size_t>().swap(factor_id_cur);
  }
  frame_num ++;
  // if (gnss_online_init || !nolidar || frame_num >= 10)
  {
    runISAM2opt();
    // state.cov.block<3,3>(0, 0) = isam.marginalCovariance(R(frame_num-1));
    // state.cov.block<6,6>(3, 3) = isam.marginalCovariance(F(frame_num-1)).block<6, 6>(0, 0);

    if (nolidar)
    {
      state.rot = isamCurrentEstimate.at<gtsam::Rot3>(R(frame_num-1)).matrix();
      state.pos = isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(0);
      state.vel = isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(3);
      state.ba = isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(6);
      state.bg = isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(9);
      state.gravity = ecef2rotation(state.pos) * gravity_init;
    }
    else
    {
      state_const_.rot = isamCurrentEstimate.at<gtsam::Rot3>(R(frame_num-1)).matrix();
      // gtsam::Point3 update_pos
      state_const_.pos = isamCurrentEstimate.at<gtsam::Vector6>(A(frame_num-1)).segment<3>(0);
      // gtsam::Vector3 update_vel;
      state_const_.vel = isamCurrentEstimate.at<gtsam::Vector6>(A(frame_num-1)).segment<3>(3);
      // state_const_.bias_a = isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(6);
      // state_const_.bias_g = isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(9);
      state_const_.gravity = isamCurrentEstimate.at<gtsam::Rot3>(P(0)).matrix().transpose() * ecef2rotation(isamCurrentEstimate.at<gtsam::Vector3>(E(0))) * gravity_init;
      state.gravity = state_const_.gravity;
    }
    // state_const_ = state;
    if (!nolidar && invalid_lidar)
    {
      state = state_const_;
      // state.cov = MD(21, 21)::Identity() * 0.01;
      printf("invalid lidar!!!!!!!!!");
      // state.cov.setIdentity();
      // state.cov.block<3, 3>(0, 0) = isam.marginalCovariance(R(frame_num-1));
      // state.cov.block<6, 6>(3, 3) = isam.marginalCovariance(F(frame_num-1)).block<6, 6>(0, 0);
      // state.cov.block<6, 6>(9, 9) = 0.001 * state.cov.block<6, 6>(9, 9);
      // state.cov.block<3, 3>(15, 15) = isam.marginalCovariance(F(frame_num-1)).block<3, 3>(9, 9);
      // state.cov.block<3, 3>(18, 18) = isam.marginalCovariance(F(frame_num-1)).block<3, 3>(6, 6);
      // invalid_lidar = false;
    }
    last_gnss_time = time_current;
    // state_const_ = state;
  }
  sat2cp[std::pair<double, int>(time2sec(curr_obs[0]->time), frame_num-1)] = curr_cp_map;
  std::map<std::pair<double, int>, std::map<uint32_t, double[6]>>::iterator it_old;
  if (!sat2cp.empty())
  {
    it_old = sat2cp.begin();
    while (time2sec(curr_obs[0]->time) - it_old->first.first > gnss_cp_time_threshold)
    {
      std::map<uint32_t, double[6]>().swap(it_old->second);
      size_t del_size = sat2cp.erase(it_old->first);
      if (sat2cp.empty()) break;
      it_old = sat2cp.begin();
    }
  }
  return true;
}

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

void GNSSProcess::rinex2iono_params(const std::string &rinex_filepath, std::vector<double> &iono_params)
{
    iono_params.resize(8);
    std::ifstream file(rinex_filepath);
    std::string line;

    // check first line, mainly RINEX version
    if (!(std::getline(file, line) && line.find("RINEX VERSION")  != std::string::npos 
            && line.find("3.04") != std::string::npos))
    {
        LOG(ERROR) << "Only RINEX 3.04 is supported";
        return;
    }

    bool find_alpha = false, find_beta = false;
    while(std::getline(file, line))
    {
        if (line.find("IONOSPHERIC CORR") != std::string::npos && line.find("GPSA") != std::string::npos)
        {

            // parse ion alpha value
            for (size_t i = 0; i < 4; ++i)
            {
                std::string value_str = line.substr(i*12+5, 12);
                iono_params[i] = str2double(value_str);
            }
            find_alpha = true;
        }
        else if (line.find("IONOSPHERIC CORR") != std::string::npos && line.find("GPSB") != std::string::npos)
        {
            // parse ion beta value
            for (size_t i = 0; i < 4; ++i)
            {
                std::string value_str = line.substr(i*12+5, 12);
                iono_params[i+4] = str2double(value_str);
            }
            find_beta = true;
        }

        if(find_alpha && find_beta)
            break;
    }
    file.close();
}

void GNSSProcess::processIMU(double dt, const Vector3d &linear_acceleration, const Vector3d &angular_velocity)
{
  Eigen::Vector3d dr = angular_velocity * dt;
  state_.rot.boxplus(dr);

  Eigen::Vector3d acc_imu = state_.rot * linear_acceleration + state_.gravity;

  state_.pos += state_.vel * dt + 0.5 * acc_imu * dt * dt;

  state_.vel += acc_imu * dt;
}

void GNSSProcess::processIMUOutput(double dt, const Vector3d &linear_acceleration, const Vector3d &angular_velocity)
{
  Eigen::Vector3d dr = angular_velocity * dt;
  state_const_.rot.boxplus(dr);

  Eigen::Vector3d acc_imu = state_const_.rot * linear_acceleration + state_const_.gravity;

  state_const_.pos += state_const_.vel * dt + 0.5 * acc_imu * dt * dt;

  state_const_.vel += acc_imu * dt;
}

// #endif


