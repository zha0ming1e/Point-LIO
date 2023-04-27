#include "GNSS_Processing_fg.h"

GNSSProcess::GNSSProcess()
    : diff_t_gnss_local(0.0)
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
  for (size_t i = 0; i < WINDOW_SIZE+1; i++)
  {
    std::vector<ObsPtr> empty_vec_o;
    std::vector<EphemBasePtr> empty_vec_e;
    gnss_meas_buf[i].swap(empty_vec_o);
    gnss_ephem_buf[i].swap(empty_vec_e);
  }
  p_assign->change_ext = 1;
  std::map<uint32_t, uint32_t> empty_map_t;
  p_assign->sat_track_status.swap(empty_map_t);
  p_assign->gtSAMgraph.resize(0); 
  p_assign->initialEstimate.clear();
  p_assign->isamCurrentEstimate.clear();
  // index_delete = 0;
  frame_delete = 0;
  // E_num = 0;
  p_assign->factor_id_frame.clear();
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
  parameters.relinearizeSkip = 5; // may matter? improtant!
  p_assign->isam = gtsam::ISAM2(parameters);
}

void GNSSProcess::inputIonoParams(double ts, const std::vector<double> &iono_params) // 
{
  if (iono_params.size() != 8)    return;

  // update ionosphere parameters
  std::vector<double> empty_vec_d;
  p_assign->latest_gnss_iono_params.swap(empty_vec_d);
  std::copy(iono_params.begin(), iono_params.end(), std::back_inserter(p_assign->latest_gnss_iono_params));
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

  if (gnss_ready)
  {
    Eigen::Vector3d pos_gnss = state.pos + state.rot * Tex_imu_r; // .normalized()
    updateGNSSStatistics(pos_gnss);
  }
  p_assign->processGNSSBase(gnss_meas, valid_meas, valid_ephems, gnss_ready, ecef_pos);
  
  if (!gnss_ready)
  {
    if (valid_meas.empty() || valid_meas.size() < 5) return; // right or not?
    {
      rot_window[frame_count] = state.rot; //.normalized().toRotationMatrix();
      pos_window[frame_count] = state.pos + state.rot * Tex_imu_r; // .normalized()
      Eigen::Matrix3d omg_skew;
      omg_skew << SKEW_SYM_MATRX(omg);
      vel_window[frame_count] = state.vel + state.rot * omg_skew * Tex_imu_r; // .normalized().toRotationMatrix()
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

  if (gnss_ready)
  {
    Eigen::Vector3d pos_gnss = state.pos + state.rot * Tex_imu_r; // .normalized()
    updateGNSSStatistics(pos_gnss);
  }
  p_assign->processGNSSBase(gnss_meas, valid_meas, valid_ephems, gnss_ready, ecef_pos);
  
  if (!gnss_ready)
  {
    if (valid_meas.empty() || valid_meas.size() < 5) return; // right or not?
    {
      rot_window[frame_count] = state.rot; //.normalized().toRotationMatrix();
      pos_window[frame_count] = state.pos + state.rot * Tex_imu_r; // .normalized()
      Eigen::Matrix3d omg_skew;
      omg_skew << SKEW_SYM_MATRX(state.omg);
      vel_window[frame_count] = state.vel + state.rot * omg_skew * Tex_imu_r; // .normalized().toRotationMatrix()
      // vel_window[frame_count] = state.vel;
    }
    gnss_meas_buf[frame_count] = valid_meas; 
    gnss_ephem_buf[frame_count] = valid_ephems;
    frame_count ++;
    gnss_ready = GNSSLIAlign();
    if (gnss_ready)
    {
      ROS_INFO("GNSS Initialization is done");
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
      if (!p_assign->factor_id_frame.empty())       
      {
        // if (frame_delete > 0)
        {
        for (size_t i = 0; i < p_assign->factor_id_frame[0].size(); i++)
        {
          // if (p_assign->factor_id_frame[0][i] != 0 && p_assign->factor_id_frame[0][i] != 1 || nolidar)
          {
            delete_factor.push_back(p_assign->factor_id_frame[0][i]);
          }
        }
        // index_delete += p_assign->factor_id_frame[0].size();
        }
      
        p_assign->factor_id_frame.pop_front();
        frame_delete ++;
      }
      if (p_assign->factor_id_frame.empty()) break;
    }
    }

    if (delete_happen)
    {
      p_assign->delete_variables(nolidar, frame_delete, frame_num, id_accumulate, delete_factor);
    }
    else
    {
      p_assign->isam.update(p_assign->gtSAMgraph, p_assign->initialEstimate);
      p_assign->gtSAMgraph.resize(0); // will the initialEstimate change?
      p_assign->initialEstimate.clear();
      p_assign->isam.update();
    }
  }
  else
  {
    p_assign->isam.update(p_assign->gtSAMgraph, p_assign->initialEstimate);
    p_assign->gtSAMgraph.resize(0); // will the initialEstimate change?
    p_assign->initialEstimate.clear();
    p_assign->isam.update();
  }
  p_assign->isamCurrentEstimate = p_assign->isam.calculateEstimate();
  
  if (nolidar) // || invalid_lidar)
  {
    pre_integration->repropagate(p_assign->isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(6),
                                p_assign->isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(9));
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
        for (uint32_t j = frame_count; j < wind_size+1; ++j) // wind_size-i
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

  std::vector<std::vector<ObsPtr>> curr_gnss_meas_buf;
  std::vector<std::vector<EphemBasePtr>> curr_gnss_ephem_buf;
  for (uint32_t i = 0; i < (wind_size+1); ++i)
  {
      curr_gnss_meas_buf.push_back(gnss_meas_buf[i]);
      curr_gnss_ephem_buf.push_back(gnss_ephem_buf[i]);
  }

  GNSSLIInitializer gnss_li_initializer(curr_gnss_meas_buf, curr_gnss_ephem_buf, p_assign->latest_gnss_iono_params);

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
  {
    bool dt_fail[4];
    memset(dt_fail, false, sizeof(dt_fail));
    size_t num_dt_fail = 0, dt_success = 0;
    for (uint32_t k_ = 0; k_ < 4; k_++)
    {
      if (rough_xyzt(3+k_) == 0)
      {
        num_dt_fail += 1;
        dt_fail[k_] = true;
      }
      else
      {
        dt_success = k_;
      }
    }
    if (num_dt_fail == 4)
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
    anc_ecef = rough_xyzt.head<3>();
    anc_local = Rot_gnss_init.transpose() * pos_window[0]; // [WINDOW_SIZE]; // ?
    yaw_enu_local = 0.0;
    R_ecef_enu = ecef2rotation(anc_ecef);
    para_rcv_ddt[wind_size] = 128.0;
    for (uint32_t k_ = 0; k_ < 4; ++k_)
    {
      if (dt_fail[k_])
      {
        para_rcv_dt[wind_size*4+k_] = rough_xyzt(3+dt_success);
      }
      else
      {
        para_rcv_dt[wind_size*4+k_] = rough_xyzt(3+k_);
      }
    }
    SetInit();
    frame_num = 1; // frame_count;
    last_gnss_time = time2sec(gnss_meas_buf[wind_size][0]->time);
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
    Eigen::Vector3d anc_cur;
    Eigen::Matrix3d R_enu_local_;
    // if (frame_num == 1)
    // {
    //   anc_cur = anc_ecef;
    //   R_enu_local_ = R_ecef_enu * Eigen::AngleAxisd(yaw_enu_local, Eigen::Vector3d::UnitZ()) * Rot_gnss_init;
    // }
    // else
    {
      anc_cur = p_assign->isamCurrentEstimate.at<gtsam::Vector3>(E(0));
      R_enu_local_ = p_assign->isamCurrentEstimate.at<gtsam::Rot3>(P(0)).matrix();
    }
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
  if (gnss_meas_buf[0].empty()) // || gnss_meas_buf[0].size() < 4)
  {
    // cout << "no valid gnss" << endl;
    return false;
  }

  double time_current = time2sec(gnss_meas_buf[0][0]->time);
  double delta_t = time_current - last_gnss_time;

  gtsam::Rot3 rel_rot = gtsam::Rot3(pre_integration->delta_q);
  gtsam::Point3 rel_pos, pos;
  gtsam::Vector3 rel_vel, vel; 
  Eigen::Matrix3d rot = Eigen::Matrix3d::Identity();
  rel_pos = pre_integration->delta_p;
  rel_vel = pre_integration->delta_v; 
  if (!nolidar) // && !invalid_lidar)
  {
    // Eigen::Matrix3d last_rot = state_last.rot.normalized().toRotationMatrix(); // isamCurrentEstimate.at<gtsam::Rot3>(R(frame_num-1)).matrix(); // state_const_.rot; // 
    // // cout << "check time period" << pre_integration->sum_dt << ";" << time_current - last_gnss_time <<  endl;
    // Eigen::Vector3d last_pos = state_last.pos; // isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(0); // state_.pos; // 
    // Eigen::Vector3d last_vel = state_last.vel; // isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(3); // state_.vel; //
    // Eigen::Vector3d cur_grav = state.rot.transpose() * state.gravity; //  
    // rel_rot = gtsam::Rot3(last_rot.transpose() * state.rot.normalized().toRotationMatrix());
    rot = state.rot; //.normalized().toRotationMatrix(); 
    // rel_pos = last_rot.transpose() * (state.pos - last_pos - last_vel * delta_t - 0.5 * state.gravity * delta_t * delta_t); 
    pos = state.pos; // (state.pos - last_pos);
    // rel_vel = last_rot.transpose() * (state.vel - last_vel - state.gravity * delta_t); // (state.vel - last_vel);
    vel = state.vel; // 
    // gtsam::Rot3 rel_rot(pre_integration->delta_q); // (last_rot.transpose() * state.rot);
    // // double delta_t = time_current - last_gnss_time;
    // gtsam::Point3 rel_pos(pre_integration->delta_p);
    // gtsam::Vector3 rel_v(pre_integration->delta_v); 
  }
  else
  {
    pos = state.ba;
    vel = state.bg;
  }
  
  if (!nolidar) // && invalid_lidar)
  {
    Eigen::Matrix<double, 6, 1> init_vel_bias_vector_imu;
    init_vel_bias_vector_imu.block<3,1>(0,0) = state.pos;
    init_vel_bias_vector_imu.block<3,1>(3,0) = state.vel;
    // init_vel_bias_vector_imu.block<3,1>(6,0) = state_.ba;
    // init_vel_bias_vector_imu.block<3,1>(9,0) = state_.bg;
    p_assign->initialEstimate.insert(A(frame_num), gtsam::Vector6(init_vel_bias_vector_imu));
    p_assign->initialEstimate.insert(R(frame_num), gtsam::Rot3(state.rot));  // .normalized().toRotationMatrix()
  }
  else
  {
    Eigen::Matrix<double, 12, 1> init_vel_bias_vector;
    init_vel_bias_vector.block<3,1>(0,0) = state.pos;
    init_vel_bias_vector.block<3,1>(3,0) = state.vel;
    init_vel_bias_vector.block<3,1>(6,0) = state.ba;
    init_vel_bias_vector.block<3,1>(9,0) = state.bg;
    p_assign->initialEstimate.insert(F(frame_num), gtsam::Vector12(init_vel_bias_vector));
    p_assign->initialEstimate.insert(R(frame_num), gtsam::Rot3(state.rot)); // .normalized().toRotationMatrix()
  }              
  rot_pos = state.rot; //.normalized().toRotationMatrix();
  if (AddFactor(rel_rot, rel_pos, rel_vel, state.gravity, delta_t, time_current, pos, vel, omg, rot))
  {
    frame_num ++;
    if (!nolidar) state_last = state;
    runISAM2opt();
  }
  else
  {
    return false;
  }

  // state.cov.block<3,3>(0, 0) = isam.marginalCovariance(R(frame_num-1));
  // state.cov.block<6,6>(3, 3) = isam.marginalCovariance(F(frame_num-1)).block<6, 6>(0, 0);
  
  if (nolidar)
  {
    state.rot = p_assign->isamCurrentEstimate.at<gtsam::Rot3>(R(frame_num-1)).matrix();
    // state.rot.normalize();
    state.pos = p_assign->isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(0);
    state.vel = p_assign->isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(3);
    state.ba = p_assign->isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(6);
    state.bg = p_assign->isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(9);
    state.gravity = ecef2rotation(state.pos) * gravity_init;
  }
  else
  {
    state_.rot = p_assign->isamCurrentEstimate.at<gtsam::Rot3>(R(frame_num-1)).matrix();
    // state_.rot.normalize();
    state_.pos = p_assign->isamCurrentEstimate.at<gtsam::Vector6>(A(frame_num-1)).segment<3>(0);
    state_.vel = p_assign->isamCurrentEstimate.at<gtsam::Vector6>(A(frame_num-1)).segment<3>(3);
    // state_.ba = p_assign->isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(6);
    // state_.bg = p_assign->isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(9);
    // state.gravity = p_assign->isamCurrentEstimate.at<gtsam::Rot3>(P(0)).matrix().transpose() * gravity_init;
    state_.gravity = p_assign->isamCurrentEstimate.at<gtsam::Rot3>(P(0)).matrix().transpose() * ecef2rotation(p_assign->isamCurrentEstimate.at<gtsam::Vector3>(E(0))) * gravity_init;
    // state.gravity = state_.gravity; // 
    // state.rot = p_assign->isamCurrentEstimate.at<gtsam::Rot3>(R(frame_num-1)).matrix();
    // state.pos = p_assign->isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(0);
    // state.vel = p_assign->isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(3);
    // state.ba = p_assign->isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(6);
    // state.bg = p_assign->isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(9);
    // state.gravity = state_.gravity; // ecef2rotation(state.pos) * gravity_init;
  }    
  
  last_gnss_time = time_current;
  // state_ = state;
  // const std::vector<ObsPtr> &curr_obs = gnss_meas_buf[0];
  std::map<std::pair<double, int>, std::map<uint32_t, double[6]>>::iterator it_old;
  if (!sat2cp.empty())
  {
    it_old = sat2cp.begin();
    while (time_current - it_old->first.first > gnss_cp_time_threshold)
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

  double time_current = time2sec(gnss_meas_buf[0][0]->time);
  double delta_t = time_current - last_gnss_time;

  gtsam::Rot3 rel_rot = gtsam::Rot3(pre_integration->delta_q);;
  gtsam::Point3 rel_pos, pos;
  gtsam::Vector3 rel_vel, vel; 
  Eigen::Matrix3d rot = Eigen::Matrix3d::Identity();
  rel_pos = pre_integration->delta_p;
  rel_vel = pre_integration->delta_v; 
  if (!nolidar) // && !invalid_lidar)
  {
    // Eigen::Matrix3d last_rot = state_const_last.rot.normalized().toRotationMatrix(); // isamCurrentEstimate.at<gtsam::Rot3>(R(frame_num-1)).matrix(); // state_const_.rot; // 
    // // cout << "check time period" << pre_integration->sum_dt << ";" << time_current - last_gnss_time <<  endl;
    // Eigen::Vector3d last_pos = state_const_last.pos; // isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(0); // state_.pos; // 
    // Eigen::Vector3d last_vel = state_const_last.vel; // isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(3); // state_.vel; //
    // Eigen::Vector3d cur_grav = state.rot.transpose() * state.gravity; //  
    rot = state.rot; //.normalized().toRotationMatrix();
    // rel_rot = gtsam::Rot3(last_rot.transpose() * state.rot.normalized().toRotationMatrix());
    pos = state.pos; // last_rot.transpose() * (state.pos - last_pos - last_vel * delta_t - 0.5 * state.gravity * delta_t * delta_t); 
    // (state.pos - last_pos);
    vel = state.vel; // last_rot.transpose() * (state.vel - last_vel - state.gravity * delta_t); // (state.vel - last_vel);
    // gtsam::Rot3 rel_rot(pre_integration->delta_q); // (last_rot.transpose() * state.rot);
    // // double delta_t = time_current - last_gnss_time;
    // gtsam::Point3 rel_pos(pre_integration->delta_p);
    // gtsam::Vector3 rel_v(pre_integration->delta_v); 
  }
  else
  {
    pos = state.ba;
    vel = state.bg;
  }
  
  if (!nolidar) // && invalid_lidar)
  {
    Eigen::Matrix<double, 6, 1> init_vel_bias_vector_imu;
    init_vel_bias_vector_imu.block<3,1>(0,0) = state.pos;
    init_vel_bias_vector_imu.block<3,1>(3,0) = state.vel;
    p_assign->initialEstimate.insert(A(frame_num), gtsam::Vector6(init_vel_bias_vector_imu));
    p_assign->initialEstimate.insert(R(frame_num), gtsam::Rot3(state.rot));  // .normalized().toRotationMatrix()
  }
  else
  {
    Eigen::Matrix<double, 12, 1> init_vel_bias_vector;
    init_vel_bias_vector.block<3,1>(0,0) = state.pos;
    init_vel_bias_vector.block<3,1>(3,0) = state.vel;
    init_vel_bias_vector.block<3,1>(6,0) = state.ba;
    init_vel_bias_vector.block<3,1>(9,0) = state.bg;
    p_assign->initialEstimate.insert(F(frame_num), gtsam::Vector12(init_vel_bias_vector));
    p_assign->initialEstimate.insert(R(frame_num), gtsam::Rot3(state.rot)); // .normalized().toRotationMatrix()
  }              
  rot_pos = state.rot; //.normalized().toRotationMatrix();
  if (AddFactor(rel_rot, rel_pos, rel_vel, state.gravity, delta_t, time_current, pos, vel, state.omg, rot))
  {
    frame_num ++;
    if (!nolidar) state_const_last = state;
    runISAM2opt();
  }
  else
  {
    return false;
  }

  // state.cov.block<3,3>(0, 0) = isam.marginalCovariance(R(frame_num-1));
  // state.cov.block<6,6>(3, 3) = isam.marginalCovariance(F(frame_num-1)).block<6, 6>(0, 0);
  
  if (nolidar)
  {
    state.rot = p_assign->isamCurrentEstimate.at<gtsam::Rot3>(R(frame_num-1)).matrix();
    // state.rot.normalize();
    state.pos = p_assign->isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(0);
    state.vel = p_assign->isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(3);
    state.ba = p_assign->isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(6);
    state.bg = p_assign->isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(9);
    state.gravity = ecef2rotation(state.pos) * gravity_init;
  }
  else
  {
    state_const_.rot = p_assign->isamCurrentEstimate.at<gtsam::Rot3>(R(frame_num-1)).matrix();
    // state_const_.rot.normalize();
    state_const_.pos = p_assign->isamCurrentEstimate.at<gtsam::Vector6>(A(frame_num-1)).segment<3>(0);
    state_const_.vel = p_assign->isamCurrentEstimate.at<gtsam::Vector6>(A(frame_num-1)).segment<3>(3);
    // state_.ba = p_assign->isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(6);
    // state_.bg = p_assign->isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(9);
    // state.gravity = p_assign->isamCurrentEstimate.at<gtsam::Rot3>(P(0)).matrix().transpose() * gravity_init;
    state_const_.gravity = p_assign->isamCurrentEstimate.at<gtsam::Rot3>(P(0)).matrix().transpose() * ecef2rotation(p_assign->isamCurrentEstimate.at<gtsam::Vector3>(E(0))) * gravity_init;
    // state.gravity = state_const_.gravity; // 
    // state.rot = p_assign->isamCurrentEstimate.at<gtsam::Rot3>(R(frame_num-1)).matrix();
    // state.pos = p_assign->isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(0);
    // state.vel = p_assign->isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(3);
    // state.ba = p_assign->isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(6);
    // state.bg = p_assign->isamCurrentEstimate.at<gtsam::Vector12>(F(frame_num-1)).segment<3>(9);
    // state.gravity = state_.gravity; // ecef2rotation(state.pos) * gravity_init;
  }    
  last_gnss_time = time_current;
  std::map<std::pair<double, int>, std::map<uint32_t, double[6]>>::iterator it_old;
  if (!sat2cp.empty())
  {
    it_old = sat2cp.begin();
    while (time_current - it_old->first.first > gnss_cp_time_threshold)
    {
      std::map<uint32_t, double[6]>().swap(it_old->second);
      size_t del_size = sat2cp.erase(it_old->first);
      if (sat2cp.empty()) break;
      it_old = sat2cp.begin();
    }
  }
  return true;
}

bool GNSSProcess::AddFactor(gtsam::Rot3 rel_rot, gtsam::Point3 rel_pos, gtsam::Vector3 rel_vel, Eigen::Vector3d state_gravity, double delta_t, double time_current,
                Eigen::Vector3d ba, Eigen::Vector3d bg, Eigen::Vector3d omg, Eigen::Matrix3d rot)
{
  double rcv_dt[4];
  bool rcv_sys[4];
  rcv_sys[0] = false; rcv_sys[1] = false; rcv_sys[2] = false; rcv_sys[3] = false;
  double rcv_ddt;
  invalid_lidar = false;
  if (!nolidar)
  {
    invalid_lidar = nolidar_cur;
    size_t num_norm = norm_vec_holder.size();
    if (num_norm > 2) // 10)
    {
    Eigen::MatrixXd A(num_norm, 3);
    for (size_t i = 0; i < num_norm; i++)
    {
      A.row(i) = norm_vec_holder[i].transpose();
    }
    std::vector<Eigen::Vector3d>().swap(norm_vec_holder);
    BDCSVD<Eigen::MatrixXd> svd(A, ComputeThinU | ComputeThinV);
    Eigen::VectorXd singular_values = svd.singularValues().cwiseAbs();
    std::sort(singular_values.data(), singular_values.data() + singular_values.size(), std::greater<double>());
    if (singular_values[1] / singular_values[0] < 0.1) invalid_lidar = true;
    }
    else
    {
      invalid_lidar = true;
      std::vector<Eigen::Vector3d>().swap(norm_vec_holder);
    }
  }
  if ((delta_t > 15 * gnss_sample_period && nolidar) || (delta_t > 15 * gnss_sample_period && invalid_lidar && !nolidar))
  {
    Reset();
    return false;
  }
  if (nolidar_cur) nolidar_cur = false;
  rcv_ddt = p_assign->isamCurrentEstimate.at<gtsam::Vector1>(C(frame_num-1))[0];
  rcv_dt[0] = p_assign->isamCurrentEstimate.at<gtsam::Vector4>(B(frame_num-1))[0] + rcv_ddt * delta_t;
  rcv_dt[1] = p_assign->isamCurrentEstimate.at<gtsam::Vector4>(B(frame_num-1))[1] + rcv_ddt * delta_t;
  rcv_dt[2] = p_assign->isamCurrentEstimate.at<gtsam::Vector4>(B(frame_num-1))[2] + rcv_ddt * delta_t;
  rcv_dt[3] = p_assign->isamCurrentEstimate.at<gtsam::Vector4>(B(frame_num-1))[3] + rcv_ddt * delta_t;

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
  omg_skew << SKEW_SYM_MATRX(omg);
  Eigen::Vector3d hat_omg_T = omg_skew * Tex_imu_r;
   
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
    values[18] = p_assign->latest_gnss_iono_params[0]; values[19] = p_assign->latest_gnss_iono_params[1]; values[20] = p_assign->latest_gnss_iono_params[2]; values[21] = p_assign->latest_gnss_iono_params[3]; 
    values[22] = p_assign->latest_gnss_iono_params[4]; values[23] = p_assign->latest_gnss_iono_params[5]; values[24] = p_assign->latest_gnss_iono_params[6]; values[25] = p_assign->latest_gnss_iono_params[7]; 
    values[26] = time_current; values[27] = freq; values[28] = curr_obs[j]->psr[freq_idx]; values[29] = curr_obs[j]->dopp[freq_idx];

    rcv_sys[sys_idx] = true;
    if (!nolidar)
    { 
      p_assign->gtSAMgraph.add(glio::GnssPsrDoppFactor(R(frame_num), A(frame_num), B(frame_num), C(frame_num), E(0), P(0), invalid_lidar, values, sys_idx, hat_omg_T, p_assign->robustpsrdoppNoise));
      // p_assign->gtSAMgraph.add(glio::GnssPsrDoppFactorPos(A(frame_num), B(frame_num), C(frame_num), E(0), P(0), invalid_lidar, values, sys_idx, rot_pos, hat_omg_T, p_assign->robustpsrdoppNoise));
    }
    else
    {    
      p_assign->gtSAMgraph.add(glio::GnssPsrDoppFactorNolidar(R(frame_num), F(frame_num), B(frame_num), C(frame_num), values, sys_idx, hat_omg_T, p_assign->robustpsrdoppNoise)); // not work
      // p_assign->gtSAMgraph.add(glio::GnssPsrDoppFactorNolidarPos(F(frame_num), B(frame_num), C(frame_num), values, sys_idx, hat_omg_T, rot_pos, p_assign->robustpsrdoppNoise)); // not work
    }
    factor_id_cur.push_back(id_accumulate);
    id_accumulate += 1;
  }
  sat2cp[std::pair<double, int>(time2sec(curr_obs[0]->time), frame_num)] = curr_cp_map;
  // if (frame_num == 1 && !nolidar)
  // {
  //   Eigen::Matrix3d R_enu_local_;
  //   R_enu_local_ = R_ecef_enu * Rot_gnss_init; // Eigen::AngleAxisd(yaw_enu_local, Eigen::Vector3d::UnitZ()) * 
  //   p_assign->initialEstimate.insert(E(0), gtsam::Vector3(anc_ecef[0], anc_ecef[1], anc_ecef[2]));
  //   p_assign->initialEstimate.insert(P(0), gtsam::Rot3(R_enu_local_));

  //   gtsam::PriorFactor<gtsam::Rot3> init_rot_ext(P(0), gtsam::Rot3(gtsam::Rot3(R_enu_local_)), p_assign->margrotNoise);
  //   gtsam::PriorFactor<gtsam::Vector3> init_pos_ext(E(0), gtsam::Vector3(anc_ecef[0], anc_ecef[1], anc_ecef[2]), p_assign->margNoise);
  //   p_assign->gtSAMgraph.add(init_rot_ext);
  //   p_assign->gtSAMgraph.add(init_pos_ext);
  //   p_assign->factor_id_frame[frame_num-1].push_back(id_accumulate);
  //   p_assign->factor_id_frame[frame_num-1].push_back(id_accumulate+1);
  //   id_accumulate += 2;
  // }
  p_assign->gtSAMgraph.add(glio::DdtSmoothFactor(C(frame_num-1), C(frame_num), p_assign->ddtNoise));
  // p_assign->gtSAMgraph.add(gtsam::PriorFactor<gtsam::Vector1>(C(frame_num), gtsam::Vector1(rcv_ddt), p_assign->ddtNoise));
  p_assign->gtSAMgraph.add(glio::DtDdtFactor(B(frame_num-1), B(frame_num), C(frame_num-1), C(frame_num), rcv_sys, delta_t, p_assign->dtNoise)); // not work
  p_assign->factor_id_frame[frame_num-1-frame_delete].push_back(id_accumulate);
  p_assign->factor_id_frame[frame_num-1-frame_delete].push_back(id_accumulate+1);
  id_accumulate += 2;
  if (!nolidar)
  {
    // Eigen::Matrix<double, 15, 15> state_cov = Eigen::Matrix<double, 15, 15>::Identity();
    // state_cov.block<9, 9>(0, 0) = state.cov.block<9, 9>(0, 0) * 100;
    // gtsam::noiseModel::Gaussian::shared_ptr LioNoise = gtsam::noiseModel::Gaussian::Covariance(state_cov); 
    if (invalid_lidar)
    {
      p_assign->gtSAMgraph.add(glio::GnssLioHardFactor(R(frame_num), A(frame_num), ba, bg, rot, sqrt_lidar, p_assign->odomNoise)); //LioNoise)); // odomNoiseIMU));
      factor_id_cur.push_back(id_accumulate);
      id_accumulate += 1;
      // odo_weight = 3.0;
      // p_assign->gtSAMgraph.add(glio::GnssLioFactor(R(frame_num-1), A(frame_num-1), R(frame_num), A(frame_num), 1.0, rel_rot, rel_pos, rel_vel,
                      // state_gravity, delta_t, pre_integration->covariance, p_assign->odomNoise)); // odomNoiseIMU)); // LioNoise)); // 
      // p_assign->factor_id_frame[frame_num-1-frame_delete].push_back(id_accumulate);
      // id_accumulate += 1;
      // Eigen::Vector3d anc_cur = p_assign->isamCurrentEstimate.at<gtsam::Vector3>(E(0));
      // Eigen::Matrix3d R_enu_local_ = p_assign->isamCurrentEstimate.at<gtsam::Rot3>(P(0)).matrix();
      // gtsam::noiseModel::Gaussian::shared_ptr updatedERNoise = gtsam::noiseModel::Gaussian::Covariance(p_assign->isam.marginalCovariance(P(0)) * 10); // important
      // gtsam::noiseModel::Gaussian::shared_ptr updatedEPNoise = gtsam::noiseModel::Gaussian::Covariance(p_assign->isam.marginalCovariance(E(0)) * 10); // important
      // gtsam::PriorFactor<gtsam::Rot3> init_ER(P(0),p_assign->isamCurrentEstimate.at<gtsam::Rot3>(P(0)), updatedERNoise); // 
      // gtsam::PriorFactor<gtsam::Vector3> init_EP(E(0),p_assign->isamCurrentEstimate.at<gtsam::Vector3>(E(0)), updatedEPNoise); // 
      // p_assign->gtSAMgraph.add(init_ER);
      // p_assign->gtSAMgraph.add(init_EP);
      // gtsam::PriorFactor<gtsam::Rot3> init_rot_ext(P(0), gtsam::Rot3(gtsam::Rot3(R_enu_local_)), p_assign->margdtNoise); // maybe not need, no, need and should be small.
      // gtsam::PriorFactor<gtsam::Vector3> init_pos_ext(E(0), gtsam::Vector3(anc_cur[0], anc_cur[1], anc_cur[2]), p_assign->margddtNoise);
      // p_assign->gtSAMgraph.add(init_rot_ext);
      // p_assign->gtSAMgraph.add(init_pos_ext);
      // factor_id_cur.push_back(id_accumulate);
      // factor_id_cur.push_back(id_accumulate+1);
      // p_assign->factor_id_frame[frame_num-1-frame_delete].push_back(id_accumulate);
      // p_assign->factor_id_frame[frame_num-1-frame_delete].push_back(id_accumulate+1);
      // id_accumulate += 2;
    }
    else
    {
      // gtSAMgraph.add(glio::GnssLioGravFactor(R(frame_num), cur_grav, state.gravity, odomNoise));
      // Eigen::Matrix<double, 15, 15> state_cov = Eigen::Matrix<double, 15, 15>::Identity(); 
      // state_cov.block<9, 9>(0, 0) = state.cov.block<9, 9>(0, 0) * 1000;
      // gtsam::noiseModel::Gaussian::shared_ptr LioNoise = gtsam::noiseModel::Gaussian::Covariance(state_cov);
      p_assign->gtSAMgraph.add(glio::GnssLioHardFactor(R(frame_num), A(frame_num), ba, bg, rot, sqrt_lidar, p_assign->margposNoise)); //LioNoise)); // odomNoiseIMU));
      factor_id_cur.push_back(id_accumulate);
      id_accumulate += 1;
      // odo_weight = 2.0;
    }
    // gtSAMgraph.add(glio::GnssLioFactorNolidar(R(frame_num-1), F(frame_num-1), R(frame_num), F(frame_num), rel_rot, rel_pos, rel_v, 
    //               state_.gravity, delta_t, state_.ba, state_.bg, pre_integration, odomNoiseIMU));
    // factor_id_frame[frame_num-1-frame_delete].push_back(id_accumulate);
    // id_accumulate += 1;
  }
  else
  {
    p_assign->gtSAMgraph.add(glio::GnssLioFactorNolidar(R(frame_num-1), F(frame_num-1), R(frame_num), F(frame_num), rel_rot, rel_pos, rel_vel, 
                  state_gravity, delta_t, ba, bg, pre_integration, p_assign->odomNoiseIMU));
    p_assign->factor_id_frame[frame_num-1-frame_delete].push_back(id_accumulate);
    id_accumulate += 1;
  }
  p_assign->initialEstimate.insert(C(frame_num), gtsam::Vector1(rcv_ddt));
  p_assign->initialEstimate.insert(B(frame_num), gtsam::Vector4(rcv_dt[0], rcv_dt[1], rcv_dt[2], rcv_dt[3]));  

  for (uint32_t j = 0; j < meas_index_sats.size(); j++)
  {
    double values[20];
    values[0] = Tex_imu_r[0]; values[1] = Tex_imu_r[1]; values[2] = Tex_imu_r[2]; values[3] = anc_local[0]; values[4] = anc_local[1]; values[5] = anc_local[2];
    values[6] = meas_svpos_best[j][0]; values[7] = meas_svpos_best[j][1]; values[8] = meas_svpos_best[j][2]; values[9] = meas_svpos_sats[j][0]; values[10] = meas_svpos_sats[j][1]; values[11] = meas_svpos_sats[j][2];
    values[12] = sv_pos_best[0]; values[13] = sv_pos_best[1]; values[14] = sv_pos_best[2]; values[15] = sv_pos_pair[j][0]; values[16] = sv_pos_pair[j][1]; values[17] = sv_pos_pair[j][2];
    values[18] = meas_cp_best - meas_cp[j] - meas_sats[j]; values[19] = cp_weight; 
    if (!nolidar)
    {
      p_assign->gtSAMgraph.add(glio::GnssCpFactor(E(0), P(0), R(meas_index_sats[j]), A(meas_index_sats[j]), R(frame_num), A(frame_num), invalid_lidar, values, p_assign->robustcpNoise));
      // Eigen::Matrix3d rot_before = p_assign->isamCurrentEstimate.at<gtsam::Rot3>(R(meas_index_sats[j])).matrix();
      // p_assign->gtSAMgraph.add(glio::GnssCpFactorPos(E(0), P(0), A(meas_index_sats[j]), A(frame_num), invalid_lidar, values, rot_before, rot_pos, p_assign->robustcpNoise));
    }
    else
    {
      p_assign->gtSAMgraph.add(glio::GnssCpFactorNolidar(R(meas_index_sats[j]), F(meas_index_sats[j]), R(frame_num), F(frame_num), values, p_assign->robustcpNoise)); // not work
      // Eigen::Matrix3d rot_before = p_assign->isamCurrentEstimate.at<gtsam::Rot3>(R(meas_index_sats[j])).matrix();
      // p_assign->gtSAMgraph.add(glio::GnssCpFactorNolidarPos(F(meas_index_sats[j]), F(frame_num), values, rot_before, rot_pos, p_assign->robustcpNoise)); // not work
    }
    // factor_id_cur.push_back(id_accumulate);
    p_assign->factor_id_frame[meas_index_sats[j]-frame_delete].push_back(id_accumulate);
    id_accumulate += 1;
  }

  {
    p_assign->factor_id_frame.push_back(factor_id_cur);
    std::vector<size_t>().swap(factor_id_cur);
  }
  return true;
}

void GNSSProcess::processIMU(double dt, const Vector3d &linear_acceleration, const Vector3d &angular_velocity)
{
  Eigen::Vector3d dr = angular_velocity * dt;
  state_.rot.boxplus(dr);

  Eigen::Vector3d acc_imu = state_.rot * linear_acceleration + state_.gravity; // .normalized()

  state_.pos += state_.vel * dt + 0.5 * acc_imu * dt * dt;

  state_.vel += acc_imu * dt;
}

void GNSSProcess::processIMUOutput(double dt, const Vector3d &linear_acceleration, const Vector3d &angular_velocity)
{
  Eigen::Vector3d dr = angular_velocity * dt;
  state_const_.rot.boxplus(dr);

  Eigen::Vector3d acc_imu = state_const_.rot * linear_acceleration + state_const_.gravity; // .normalized()

  state_const_.pos += state_const_.vel * dt + 0.5 * acc_imu * dt * dt;

  state_const_.vel += acc_imu * dt;
}

void GNSSProcess::SetInit()
{
  if (!nolidar)
  {
    Eigen::Matrix3d R_enu_local_;
    R_enu_local_ = R_ecef_enu * Rot_gnss_init; // * Eigen::AngleAxisd(yaw_enu_local, Eigen::Vector3d::UnitZ()) 
    // prior factor 
    Eigen::Matrix<double, 6, 1> init_vel_bias_vector;
    init_vel_bias_vector.block<3,1>(0,0) = Rot_gnss_init.transpose() * pos_window[wind_size];
    init_vel_bias_vector.block<3,1>(3,0) = Rot_gnss_init.transpose() * vel_window[wind_size];
    // dt[0] = para_rcv_dt[wind_size*4]; dt[1] = para_rcv_dt[wind_size*4+1], dt[2] = para_rcv_dt[wind_size*4+2], dt[3] = para_rcv_dt[wind_size*4+3];
    // ddt = para_rcv_ddt[wind_size];
    p_assign->initialEstimate.insert(R(0), gtsam::Rot3(Rot_gnss_init.transpose() * rot_window[wind_size]));
    // p_assign->initialEstimate.insert(F(0), gtsam::Vector12(init_vel_bias_vector));
    p_assign->initialEstimate.insert(A(0), gtsam::Vector6(init_vel_bias_vector));
    p_assign->initialEstimate.insert(B(0), gtsam::Vector4(para_rcv_dt[wind_size*4], para_rcv_dt[wind_size*4+1], para_rcv_dt[wind_size*4+2], para_rcv_dt[wind_size*4+3]));
    // p_assign->initialEstimate.insert(B(0), gtsam::Vector4(0.0, 0.0, 0.0, 0.0));
    // p_assign->initialEstimate.insert(C(0), gtsam::Vector1(para_rcv_ddt[wind_size]));
    p_assign->initialEstimate.insert(C(0), gtsam::Vector1(0.0));
    // p_assign->initialEstimate.insert(Y(0), gtsam::Vector1(yaw_enu_local));
    p_assign->initialEstimate.insert(E(0), gtsam::Vector3(anc_ecef[0], anc_ecef[1], anc_ecef[2]));
    p_assign->initialEstimate.insert(P(0), gtsam::Rot3(R_enu_local_));

    gtsam::PriorFactor<gtsam::Rot3> init_rot_ext(P(0), gtsam::Rot3(gtsam::Rot3(R_enu_local_)), p_assign->priorextrotNoise);
    gtsam::PriorFactor<gtsam::Vector3> init_pos_ext(E(0), gtsam::Vector3(anc_ecef[0], anc_ecef[1], anc_ecef[2]), p_assign->priorextposNoise);
    gtsam::PriorFactor<gtsam::Vector4> init_dt(B(0), gtsam::Vector4(para_rcv_dt[wind_size*4], para_rcv_dt[wind_size*4+1], para_rcv_dt[wind_size*4+2], para_rcv_dt[wind_size*4+3]), p_assign->priordtNoise);
    // gtsam::PriorFactor<gtsam::Vector4> init_dt(B(0), gtsam::Vector4(0.0, 0.0, 0.0, 0.0), p_assign->priordtNoise);
    // gtsam::PriorFactor<gtsam::Vector1> init_ddt(C(0), gtsam::Vector1(para_rcv_ddt[wind_size]), p_assign->priorddtNoise);
    gtsam::PriorFactor<gtsam::Vector1> init_ddt(C(0), gtsam::Vector1(0.0), p_assign->priorddtNoise);
    gtsam::PriorFactor<gtsam::Rot3> init_rot_(R(0), gtsam::Rot3(Rot_gnss_init.transpose() * rot_window[wind_size]), p_assign->priorrotNoise);
    gtsam::PriorFactor<gtsam::Vector6> init_vel_(A(0), gtsam::Vector6(init_vel_bias_vector), p_assign->priorNoise); // priorposNoise);
    // gtsam::PriorFactor<gtsam::Vector12> init_vel_(F(0), gtsam::Vector12(init_vel_bias_vector), priorposNoise);
    p_assign->gtSAMgraph.add(init_rot_ext);
    p_assign->gtSAMgraph.add(init_pos_ext);
    p_assign->gtSAMgraph.add(init_dt);
    p_assign->gtSAMgraph.add(init_ddt);
    p_assign->gtSAMgraph.add(init_rot_);
    p_assign->gtSAMgraph.add(init_vel_);
    p_assign->factor_id_frame.push_back(std::vector<size_t>{0, 1, 2, 3, 4, 5});
    id_accumulate += 6;
  }
  else
  {
  //   Eigen::Matrix3d R_enu_local_;
  //   R_enu_local_ = Eigen::AngleAxisd(yaw_enu_local, Eigen::Vector3d::UnitZ());
    // dt[0] = para_rcv_dt[wind_size*4], dt[1] = para_rcv_dt[wind_size*4+1], dt[2] = para_rcv_dt[wind_size*4+2], dt[3] = para_rcv_dt[wind_size*4+3];
    // ddt = para_rcv_ddt[wind_size];
    gtsam::PriorFactor<gtsam::Rot3> init_rot(R(0), gtsam::Rot3(R_ecef_enu * rot_window[wind_size]), p_assign->priorrotNoise); //  * R_enu_local_
    Eigen::Matrix<double, 12, 1> init_vel_bias_vector;
    init_vel_bias_vector.block<3,1>(0,0) = anc_ecef + R_ecef_enu * (pos_window[wind_size] - pos_window[0] - rot_window[wind_size] * Tex_imu_r); //  * R_enu_local_
    init_vel_bias_vector.block<3,1>(3,0) = R_ecef_enu * vel_window[wind_size]; // R_enu_local_ * 
    init_vel_bias_vector.block<6,1>(6,0) = Eigen::Matrix<double, 6, 1>::Zero();
    gtsam::PriorFactor<gtsam::Vector12> init_vel_bias(F(0), gtsam::Vector12(init_vel_bias_vector), p_assign->priorposNoise);
    gtsam::PriorFactor<gtsam::Vector4> init_dt(B(0), gtsam::Vector4(para_rcv_dt[wind_size*4], para_rcv_dt[wind_size*4+1], para_rcv_dt[wind_size*4+2], para_rcv_dt[wind_size*4+3]), p_assign->priordtNoise);
    gtsam::PriorFactor<gtsam::Vector1> init_ddt(C(0), gtsam::Vector1(0.0), p_assign->priorddtNoise); // para_rcv_ddt[wind_size]
    p_assign->gtSAMgraph.add(init_rot);
    p_assign->gtSAMgraph.add(init_vel_bias);
    p_assign->gtSAMgraph.add(init_dt);
    p_assign->gtSAMgraph.add(init_ddt);
    p_assign->factor_id_frame.push_back(std::vector<size_t>{0, 1, 2, 3}); //{i * 4, i * 4 + 1, i * 4  + 2, i * 4 + 3});
    p_assign->initialEstimate.insert(R(0), gtsam::Rot3(R_ecef_enu * rot_window[wind_size])); // R_enu_local_ * 
    p_assign->initialEstimate.insert(F(0), gtsam::Vector12(init_vel_bias_vector));
    p_assign->initialEstimate.insert(B(0), gtsam::Vector4(para_rcv_dt[wind_size*4], para_rcv_dt[wind_size*4+1], para_rcv_dt[wind_size*4+2], para_rcv_dt[wind_size*4+3]));
    p_assign->initialEstimate.insert(C(0), gtsam::Vector1(0.0)); // para_rcv_ddt[wind_size]
    id_accumulate += 4;
  }
}