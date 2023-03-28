#include "GNSS_Initialization.h"

GNSSLIInitializer::GNSSLIInitializer(const std::vector<std::vector<ObsPtr>> &gnss_meas_buf_, 
    const std::vector<std::vector<EphemBasePtr>> &gnss_ephem_buf_, const std::vector<double> &iono_params_)
        : gnss_meas_buf(gnss_meas_buf_), gnss_ephem_buf(gnss_ephem_buf_), iono_params(iono_params_)
{
    num_all_meas = 0;
    all_sat_states.clear();
    for (uint32_t i = 0; i < gnss_meas_buf.size(); ++i)
    {
        num_all_meas += gnss_meas_buf[i].size();
        all_sat_states.push_back(sat_states(gnss_meas_buf[i], gnss_ephem_buf[i]));
    }
}

bool GNSSLIInitializer::coarse_localization(Eigen::Matrix<double, 7, 1> &result)
{
    result.setZero();
    std::vector<ObsPtr> accum_obs;
    std::vector<EphemBasePtr> accum_ephems;
    for (uint32_t i = 0; i < gnss_meas_buf.size(); ++i)
    {
        std::copy(gnss_meas_buf[i].begin(), gnss_meas_buf[i].end(), std::back_inserter(accum_obs));
        std::copy(gnss_ephem_buf[i].begin(), gnss_ephem_buf[i].end(), std::back_inserter(accum_ephems));
    }
    Eigen::Matrix<double, 7, 1> xyzt = Psr_pos(accum_obs, accum_ephems, iono_params);
    // cout << "check 68:" << iono_params[0] << ";" << iono_params[7] << ";" << satsys(accum_obs[0]->sat, NULL) << ";" << xyzt.transpose() << endl;
    if (xyzt.topLeftCorner<3, 1>().norm() == 0)
    {
        std::cerr << "Failed to obtain a rough reference location.\n";
        return false;
    }

    for (uint32_t k = 0; k < 4; ++k)
    {
        if (fabs(xyzt(k+3)) < 1)
            xyzt(k+3) = 0;          // not observed yet
    }

    result = xyzt;
    return true;
}

bool GNSSLIInitializer::yaw_alignment(const std::vector<Eigen::Vector3d> &local_vs, 
    const Eigen::Vector3d &rough_anchor_ecef, double &aligned_yaw, double &rcv_ddt)
{
    aligned_yaw = 0;
    rcv_ddt = 0;

    double est_yaw = 0;
    double est_rcv_ddt = 0;

    Eigen::Matrix3d rough_R_ecef_enu = ecef2rotation(rough_anchor_ecef);
    uint32_t align_iter = 0;
    double align_dx_norm = 1.0;
    while (align_iter < MAX_ITERATION && align_dx_norm > CONVERGENCE_EPSILON)
    {
        Eigen::MatrixXd align_G(num_all_meas, 2);
        align_G.setZero();
        align_G.col(1).setOnes();
        Eigen::VectorXd align_b(num_all_meas);
        align_b.setZero();
        Eigen::Matrix3d align_R_enu_local(Eigen::AngleAxisd(est_yaw, Eigen::Vector3d::UnitZ()));
        Eigen::Matrix3d align_tmp_M;
        align_tmp_M << -sin(est_yaw), -cos(est_yaw), 0,
                        cos(est_yaw), -sin(est_yaw), 0,
                        0       , 0        , 0;
        
        uint32_t align_counter = 0;
        for (uint32_t i = 0; i < gnss_meas_buf.size(); ++i)
        {
            Eigen::Matrix<double, 4, 1> ecef_vel_ddt;
            ecef_vel_ddt.head<3>() = rough_R_ecef_enu * align_R_enu_local * local_vs[i];
            ecef_vel_ddt(3) = est_rcv_ddt;
            Eigen::VectorXd epoch_res;
            Eigen::MatrixXd epoch_J;
            dopp_res(ecef_vel_ddt, rough_anchor_ecef, gnss_meas_buf[i], all_sat_states[i], epoch_res, epoch_J);
            align_b.segment(align_counter, gnss_meas_buf[i].size()) = epoch_res;
            align_G.block(align_counter, 0, gnss_meas_buf[i].size(), 1) = 
                epoch_J.leftCols(3)*rough_R_ecef_enu*align_tmp_M*local_vs[i];
            align_counter += gnss_meas_buf[i].size();
        }
        Eigen::VectorXd dx = -(align_G.transpose()*align_G).inverse() * align_G.transpose() * align_b;
        est_yaw += dx(0);
        est_rcv_ddt += dx(1);
        align_dx_norm = dx.norm();
        ++ align_iter;
    }

    if (align_iter > MAX_ITERATION)
    {
        std::cerr << "Fail to initialize yaw offset.\n";
        return false;
    }

    aligned_yaw = est_yaw;
    if (aligned_yaw > M_PI)
        aligned_yaw -= floor(est_yaw/(2.0*M_PI) + 0.5) * (2.0*M_PI);
    else if (aligned_yaw < -M_PI)
        aligned_yaw -=  ceil(est_yaw/(2.0*M_PI) - 0.5) * (2.0*M_PI);

    rcv_ddt = est_rcv_ddt;

    return true;
}

bool GNSSLIInitializer::anchor_refinement(const std::vector<Eigen::Vector3d> &local_ps, 
    const double aligned_yaw, const double aligned_ddt, 
    const Eigen::Matrix<double, 7, 1> &rough_ecef_dt, Eigen::Matrix<double, 7, 1> &refined_ecef_dt)
{
    refined_ecef_dt.setZero();

    Eigen::Matrix3d aligned_R_enu_local(Eigen::AngleAxisd(aligned_yaw, Eigen::Vector3d::UnitZ()));

    // refine anchor point and receiver clock bias
    Eigen::Vector3d refine_anchor = rough_ecef_dt.head<3>();
    Eigen::Vector4d refine_dt = rough_ecef_dt.tail<4>();
    uint32_t refine_iter = 0;
    double refine_dx_norm = 1.0;
    Eigen::Vector3d refine_pos;
    refine_pos.setZero();
    for (uint32_t k = 0; k < local_ps.size(); k++)
    {
        refine_pos += local_ps[k];
    }
    refine_pos /= local_ps.size();
    std::vector<uint32_t> unobserved_sys;
    for (uint32_t k = 0; k < 4; ++k)
    {
        if (rough_ecef_dt(3+k) == 0)
            unobserved_sys.push_back(k);
    }

    while (refine_iter < MAX_ITERATION && refine_dx_norm > CONVERGENCE_EPSILON)
    {
        Eigen::MatrixXd refine_G(num_all_meas+unobserved_sys.size(), 7);
        Eigen::VectorXd refine_b(num_all_meas+unobserved_sys.size());
        refine_G.setZero();
        refine_b.setZero();
        uint32_t refine_counter = 0;
        Eigen::Matrix3d refine_R_ecef_enu = ecef2rotation(refine_anchor);
        Eigen::Matrix3d refine_R_ecef_local = refine_R_ecef_enu * aligned_R_enu_local;

        for (uint32_t i = 0; i < gnss_meas_buf.size(); ++i)
        {
            Eigen::Matrix<double, 7, 1> ecef_xyz_dt;
            ecef_xyz_dt.head<3>() = refine_R_ecef_local * (local_ps[i] - local_ps[0]) + refine_anchor;
            ecef_xyz_dt.tail<4>() = refine_dt + aligned_ddt * Eigen::Vector4d::Ones() * (time2sec(gnss_meas_buf[i][0]->time) - time2sec(gnss_meas_buf[0][0]->time)); // ? must be wrong !

            Eigen::VectorXd epoch_res;
            Eigen::MatrixXd epoch_J;
            std::vector<Eigen::Vector2d> tmp_atmos_delay, tmp_sv_azel;
            psr_res(ecef_xyz_dt, gnss_meas_buf[i], all_sat_states[i], iono_params, 
                epoch_res, epoch_J, tmp_atmos_delay, tmp_sv_azel);
            refine_b.segment(refine_counter, gnss_meas_buf[i].size()) = epoch_res;
            refine_G.middleRows(refine_counter, gnss_meas_buf[i].size()) = epoch_J;
            refine_counter += gnss_meas_buf[i].size();
        }
        for (uint32_t k : unobserved_sys)
        {
            refine_b(refine_counter) = 0;
            refine_G(refine_counter, k+3) = 1.0;
            ++ refine_counter;
        }

        Eigen::VectorXd dx = -(refine_G.transpose()*refine_G).inverse() * refine_G.transpose() * refine_b;
        refine_anchor += dx.head<3>();
        refine_dt += dx.tail<4>();
        refine_dx_norm = dx.norm();
        ++ refine_iter;
    }

    if (refine_iter > MAX_ITERATION)
    {
        std::cerr << "Fail to perform anchor refinement.\n";
        return false;
    }

    refined_ecef_dt.head<3>() = refine_anchor;
    refined_ecef_dt.tail<4>() = refine_dt;

    return true;
}

bool GNSSLIInitializer::yaw_refinement(const std::vector<Eigen::Vector3d> &local_vs, const Eigen::Vector3d &refined_anchor_ecef,
    const std::vector<Eigen::Vector3d> &local_ps, double &aligned_yaw, double &rcv_ddt)
{
    aligned_yaw = 0;
    rcv_ddt = 0;

    double est_yaw = 0;
    double est_rcv_ddt = 0;

    Eigen::Matrix3d refined_R_ecef_enu = ecef2rotation(refined_anchor_ecef);
    uint32_t align_iter = 0;
    double align_dx_norm = 1.0;
    while (align_iter < MAX_ITERATION && align_dx_norm > CONVERGENCE_EPSILON)
    {
        Eigen::MatrixXd align_G(num_all_meas, 2);
        align_G.setZero();
        align_G.col(1).setOnes();
        Eigen::VectorXd align_b(num_all_meas);
        align_b.setZero();
        Eigen::Matrix3d align_R_enu_local(Eigen::AngleAxisd(est_yaw, Eigen::Vector3d::UnitZ()));
        Eigen::Matrix3d align_tmp_M;
        align_tmp_M << -sin(est_yaw), -cos(est_yaw), 0,
                        cos(est_yaw), -sin(est_yaw), 0,
                        0       , 0        , 0;
        
        uint32_t align_counter = 0;
        for (uint32_t i = 0; i < gnss_meas_buf.size(); ++i)
        {
            Eigen::Matrix<double, 4, 1> ecef_vel_ddt;
            ecef_vel_ddt.head<3>() = refined_R_ecef_enu * align_R_enu_local * local_vs[i];
            ecef_vel_ddt(3) = est_rcv_ddt;
            Eigen::VectorXd epoch_res;
            Eigen::MatrixXd epoch_J;
            Eigen::Vector3d ecef_pos = refined_R_ecef_enu * align_R_enu_local * (local_ps[i] - local_ps[0]) + refined_anchor_ecef;
            dopp_res(ecef_vel_ddt, ecef_pos, gnss_meas_buf[i], all_sat_states[i], epoch_res, epoch_J);
            align_b.segment(align_counter, gnss_meas_buf[i].size()) = epoch_res;
            align_G.block(align_counter, 0, gnss_meas_buf[i].size(), 1) = 
                epoch_J.leftCols(3)*refined_R_ecef_enu*align_tmp_M*local_vs[i];
            align_counter += gnss_meas_buf[i].size();
        }
        Eigen::VectorXd dx = -(align_G.transpose()*align_G).inverse() * align_G.transpose() * align_b;
        est_yaw += dx(0);
        est_rcv_ddt += dx(1);
        align_dx_norm = dx.norm();
        ++ align_iter;
    }

    if (align_iter > MAX_ITERATION)
    {
        std::cerr << "Fail to initialize yaw offset.\n";
        return false;
    }

    aligned_yaw = est_yaw;
    if (aligned_yaw > M_PI)
        aligned_yaw -= floor(est_yaw/(2.0*M_PI) + 0.5) * (2.0*M_PI);
    else if (aligned_yaw < -M_PI)
        aligned_yaw -=  ceil(est_yaw/(2.0*M_PI) - 0.5) * (2.0*M_PI);

    rcv_ddt = est_rcv_ddt;

    return true;
}


  Eigen::Matrix<double, 7, 1> GNSSLIInitializer::Psr_pos(const std::vector<ObsPtr> &obs, 
        const std::vector<EphemBasePtr> &ephems, const std::vector<double> &iono_params)
    {
        Eigen::Matrix<double, 7, 1> result;
        result.setZero();

        std::vector<ObsPtr> valid_obs;
        std::vector<EphemBasePtr> valid_ephems;
        filter_L1(obs, ephems, valid_obs, valid_ephems);
        if (valid_obs.size() < 4)
        {
            LOG(ERROR) << "[gnss_comm::psr_pos] GNSS observation not enough.\n";
            return result;
        }

        std::vector<SatStatePtr> all_sat_states = sat_states(valid_obs, valid_ephems);

        Eigen::Matrix<double, 7, 1> xyzt;
        xyzt.setZero();
        double dx_norm = 1.0;
        uint32_t num_iter = 0;
        while(num_iter < MAX_ITER_PVT && dx_norm > EPSILON_PVT)
        {
            Eigen::VectorXd b;
            Eigen::MatrixXd G;
            std::vector<Eigen::Vector2d> atmos_delay;
            std::vector<Eigen::Vector2d> all_sv_azel;
            psr_res(xyzt, valid_obs, all_sat_states, iono_params, b, G, atmos_delay, all_sv_azel);

            std::vector<uint32_t> good_idx;
            for (uint32_t i = 0; i < valid_obs.size(); ++i)
            {
                if (G.row(i).norm() <= 0)   continue;       // res not computed
                // if (all_sv_azel[i].y() > 15.0/180.0*M_PI)  
                    good_idx.push_back(i);
            }
            int sys_mask[4] = {0, 0, 0, 0};
            // for (uint32_t i = 0; i < valid_obs.size(); ++i)
            for (uint32_t i = 0; i < good_idx.size(); ++i)
            {
                const uint32_t obs_sys = satsys(valid_obs[good_idx[i]]->sat, NULL);
                sys_mask[sys2idx.at(obs_sys)] = 1;
            }
            uint32_t num_extra_constraint = 4;
            for (uint32_t k = 0; k < 4; ++k)    num_extra_constraint -= sys_mask[k];
            LOG_IF(FATAL, num_extra_constraint >= 4) << "[gnss_comm::psr_pos] too many extra-clock constraints.\n";

            const uint32_t good_num = good_idx.size();
            LOG_IF(ERROR, good_num < 4) << "too few good obs: " << good_num;

            Eigen::MatrixXd good_G(good_num+num_extra_constraint, 7);
            Eigen::VectorXd good_b(good_num+num_extra_constraint);
            Eigen::MatrixXd good_W(good_num+num_extra_constraint, good_num+num_extra_constraint);
            good_W.setZero();
            for (uint32_t i = 0; i < good_num; ++i)
            {
                const uint32_t origin_idx = good_idx[i];
                const ObsPtr &this_obs = valid_obs[origin_idx];
                good_G.row(i) = G.row(origin_idx);
                good_b(i) = b(origin_idx);
                // compose weight
                const double sin_el = sin(all_sv_azel[origin_idx].y());
                double weight = sin_el*sin_el;
                // cout << "weight:" << weight << endl;
                int l1_idx = -1;
                L1_freq(this_obs, &l1_idx);
                LOG_IF(FATAL, l1_idx < 0) << "[gnss_comm::psr_pos] no L1 observation found.\n";

                if (this_obs->psr_std[l1_idx] > 0)
                    weight /= (this_obs->psr_std[l1_idx]/0.16);
                const uint32_t obs_sys = satsys(this_obs->sat, NULL);
                if (obs_sys == SYS_GPS || obs_sys == SYS_BDS)
                    weight /= valid_ephems[origin_idx]->ura-1;
                else if (obs_sys == SYS_GAL)
                    weight /= valid_ephems[origin_idx]->ura-2;
                else if (obs_sys == SYS_GLO)
                    weight /= 4;
                good_W(i, i) = weight;
            }
            uint32_t tmp_count = good_num;
            // add extra pseudo measurement to contraint unobservable clock bias
            for (size_t k = 0; k < 4; ++k)
            {
                if (!sys_mask[k])
                {
                    good_G.row(tmp_count).setZero();
                    good_G(tmp_count, k+3) = 1.0;
                    good_b(tmp_count) = 0;
                    good_W(tmp_count, tmp_count) = 1000;       // large weight
                    ++tmp_count;
                }
            }

            // ready for solving
            Eigen::VectorXd dx = -(good_G.transpose()*good_W*good_G).inverse() * good_G.transpose() * good_W * good_b;
            // cout << "check G and b:" << good_G << endl;
            // cout << "check w" << good_W << endl;
            dx_norm = dx.norm();
            // cout << "check" << dx.transpose() << ";" << xyzt.transpose() << endl;
            xyzt += dx;
            // LOG(INFO) << "cov is \n" << (G.transpose()*W*G).inverse();
            ++num_iter;
        }
        if (num_iter == MAX_ITER_PVT)
        {
            LOG(WARNING) << "[gnss_comm::psr_pos] XYZT solver reached maximum iterations.\n";
            return result;
        }

        result = xyzt;
        return result;
    }