#include "GNSS_Assignment.h"

GNSSAssignment::GNSSAssignment() : gnss_track_num_threshold(20) {}

void GNSSAssignment::initNoises( void ) // maybe usable!
{
    gtsam::Vector priorrotNoiseVector3(3);
    priorrotNoiseVector3 << prior_noise, prior_noise, prior_noise;
    // priorrotNoiseVector3 << prior_noise / 1000, prior_noise / 1000, prior_noise / 1000;
    priorrotNoise = gtsam::noiseModel::Diagonal::Variances(priorrotNoiseVector3);

    gtsam::Vector priorposNoiseVector12(12);
    // priorposNoiseVector12 << prior_noise / 1000, prior_noise / 1000, prior_noise / 1000, prior_noise / 1000, prior_noise / 1000, prior_noise / 1000,
    //                         prior_noise / 1000, prior_noise / 1000, prior_noise / 1000, prior_noise / 1000, prior_noise / 1000, prior_noise / 1000;
    priorposNoiseVector12 << prior_noise, prior_noise, prior_noise, prior_noise, prior_noise, prior_noise,
                            prior_noise, prior_noise, prior_noise, prior_noise, prior_noise, prior_noise;
    priorposNoise = gtsam::noiseModel::Diagonal::Variances(priorposNoiseVector12);

    // gtsam::Vector priorvelNoiseVector3(3);
    // priorvelNoiseVector3 << prior_noise, prior_noise, prior_noise;
    // priorvelNoise = gtsam::noiseModel::Diagonal::Variances(priorvelNoiseVector3);

    gtsam::Vector priorNoiseVector6(6);
    // priorNoiseVector6 << prior_noise / 1000, prior_noise / 1000, prior_noise / 1000, prior_noise / 1000, prior_noise / 1000, prior_noise / 1000; 
    priorNoiseVector6 << prior_noise, prior_noise, prior_noise, prior_noise, prior_noise, prior_noise; 
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
    margrotNoiseVector3 << prior_noise, prior_noise, prior_noise;
    margrotNoise = gtsam::noiseModel::Diagonal::Variances(margrotNoiseVector3);

    gtsam::Vector margposNoiseVector9(9);
    margposNoiseVector9 << marg_noise, marg_noise, marg_noise, marg_noise, marg_noise, marg_noise,
                            marg_noise, marg_noise, marg_noise; //, marg_noise, marg_noise, marg_noise;
    margposNoise = gtsam::noiseModel::Diagonal::Variances(margposNoiseVector9);

    // gtsam::Vector priorvelNoiseVector3(3);
    // priorvelNoiseVector3 << prior_noise, prior_noise, prior_noise;
    // priorvelNoise = gtsam::noiseModel::Diagonal::Variances(priorvelNoiseVector3);

    gtsam::Vector margNoiseVector3(3);
    margNoiseVector3 << prior_noise, prior_noise, prior_noise; //, marg_noise, marg_noise, marg_noise, marg_noise, marg_noise, marg_noise, 
                        // marg_noise, marg_noise;
    margNoise = gtsam::noiseModel::Diagonal::Variances(margNoiseVector3);

    gtsam::Vector margdtNoiseVector3(3);
    margdtNoiseVector3 << prior_noise, prior_noise, prior_noise; //, marg_noise;
    margdtNoise = gtsam::noiseModel::Diagonal::Variances(margdtNoiseVector3);

    // gtsam::Vector margExtNoiseVector4(4);
    // margExtNoiseVector4 << 1e-6, 1e-6, 1e-6, 1e-6;
    // margExtNoise = gtsam::noiseModel::Diagonal::Variances(margExtNoiseVector4);

    gtsam::Vector margddtNoiseVector3(3);
    margddtNoiseVector3 << prior_noise, prior_noise, prior_noise; //prior_noise;
    margddtNoise = gtsam::noiseModel::Diagonal::Variances(margddtNoiseVector3);

    gtsam::Vector dtNoiseVector4(4);
    dtNoiseVector4 << dt_noise, dt_noise, dt_noise, dt_noise;
    dtNoise = gtsam::noiseModel::Diagonal::Variances(dtNoiseVector4);

    gtsam::Vector ddtNoiseVector1(1);
    ddtNoiseVector1 << ddt_noise;
    ddtNoise = gtsam::noiseModel::Diagonal::Variances(ddtNoiseVector1);

    gtsam::Vector odomNoiseVector9(9);
    odomNoiseVector9 << odo_noise, odo_noise, odo_noise, odo_noise, odo_noise, odo_noise, odo_noise, odo_noise, odo_noise;
    odomNoise = gtsam::noiseModel::Diagonal::Variances(odomNoiseVector9); // should be related to the time, maybe proportional

    // gtsam::Vector odomNoiseVector3(3);
    // odomNoiseVector3 << odo_noise, odo_noise, odo_noise; //, odo_noise, odo_noise, odo_noise, odo_noise, odo_noise, odo_noise;
    // odomNoise = gtsam::noiseModel::Diagonal::Variances(odomNoiseVector3); // should be related to the imu noise
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

void GNSSAssignment::Ephemfromrinex(const std::string &rinex_filepath)
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

void GNSSAssignment::inputEphem(EphemBasePtr ephem_ptr) // 
{
    double toe = time2sec(ephem_ptr->toe);
    // if a new ephemeris comes
    if (sat2time_index.count(ephem_ptr->sat) == 0 || sat2time_index.at(ephem_ptr->sat).count(toe) == 0)
    {
        sat2ephem[ephem_ptr->sat].emplace_back(ephem_ptr);
        sat2time_index[ephem_ptr->sat].emplace(toe, sat2ephem.at(ephem_ptr->sat).size()-1);
    }
}

void GNSSAssignment::rinex2iono_params(const std::string &rinex_filepath, std::vector<double> &iono_params)
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

void GNSSAssignment::rinex2ephems(const std::string &rinex_filepath, std::map<uint32_t, std::vector<EphemBasePtr>> &sat2ephem_)
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

void GNSSAssignment::processGNSSBase(const std::vector<ObsPtr> &gnss_meas, std::vector<ObsPtr> &valid_meas, std::vector<EphemBasePtr> &valid_ephems, bool gnss_ready, Eigen::Vector3d ecef_pos)
{
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
      freq_idx_ = -1;
      L1_freq(obs, &freq_idx_);
      if (freq_idx_ < 0)   continue;              // no L1 observation
      
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
      
      freq_idx_ = -1;
      L1_freq(obs, &freq_idx_);
      if (freq_idx_ < 0)   continue;              // no L1 observation
      
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
      LOG_IF(FATAL, freq_idx_ < 0) << "No L1 observation found.\n";
      if (gnss_ready)
      {
      if (obs->psr_std[freq_idx_]  > gnss_psr_std_threshold ||
          obs->dopp_std[freq_idx_] > gnss_dopp_std_threshold) //||
          // obs->cp_std[freq_idx_] * 0.004 > gnss_cp_std_threshold)
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
      if (obs->psr_std[freq_idx_]  > gnss_psr_std_threshold / 3 ||
          obs->dopp_std[freq_idx_] > gnss_dopp_std_threshold / 3) //||
          // obs->cp_std[freq_idx_] * 0.004 > gnss_cp_std_threshold)
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
          sat_azel(ecef_pos, sat_ecef, azel); // ecef_pos should be updated for this time step // coarse value is acceptable as well TODO
          if (azel[1] < gnss_elevation_threshold*M_PI/180.0)
              continue;
      }
      valid_meas.push_back(obs);
      valid_ephems.push_back(best_ephem);
  }
}

void GNSSAssignment::delete_variables(bool nolidar, size_t frame_delete, int frame_num, size_t &id_accumulate, gtsam::FactorIndices delete_factor)
{
    if (!nolidar)
    {
      if (frame_delete > 0)
      {
        if (frame_delete >= change_ext)
        {
            // gtsam::noiseModel::Gaussian::shared_ptr updatedERNoise = gtsam::noiseModel::Gaussian::Covariance(isam.marginalCovariance(P(0)) * 1); // important
            // gtsam::noiseModel::Gaussian::shared_ptr updatedEPNoise = gtsam::noiseModel::Gaussian::Covariance(isam.marginalCovariance(E(0)) * 1); // important
            gtsam::PriorFactor<gtsam::Rot3> init_ER(P(0),isamCurrentEstimate.at<gtsam::Rot3>(P(0)), margrotNoise); // updatedERNoise); //  
            gtsam::PriorFactor<gtsam::Vector3> init_EP(E(0),isamCurrentEstimate.at<gtsam::Vector3>(E(0)), margrotNoise); // updatedEPNoise); //
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
    else
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
}  

double GNSSAssignment::str2double(const std::string &num_str)
{
    size_t D_pos = num_str.find("D");
    std::string tmp_str = num_str;
    if (D_pos != std::string::npos)
        tmp_str = tmp_str.replace(D_pos, 1, "e");
    return std::stod(tmp_str);
}

EphemPtr GNSSAssignment::rinex_line2ephem(const std::vector<std::string> &ephem_lines)
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

GloEphemPtr GNSSAssignment::rinex_line2glo_ephem(const std::vector<std::string> &ephem_lines, const uint32_t gpst_leap_seconds)
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
