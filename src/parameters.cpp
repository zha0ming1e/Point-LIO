#include "parameters.h"

bool is_first_frame = true;
double lidar_end_time = 0.0, first_lidar_time = 0.0, time_con = 0.0;
double last_timestamp_lidar = -1.0, last_timestamp_imu = -1.0;
int pcd_index = 0;

state_input state_in;
state_output state_out;
std::string lid_topic, imu_topic;
bool prop_at_freq_of_imu, check_satu, con_frame, cut_frame;
bool use_imu_as_input, space_down_sample, publish_odometry_without_downsample;
int  init_map_size, con_frame_num;
double match_s, satu_acc, satu_gyro, cut_frame_time_interval;
float  plane_thr;
double filter_size_surf_min, filter_size_map_min, fov_deg;
double cube_len; 
float  DET_RANGE;
bool   imu_en;
double imu_time_inte, gnss_ekf_noise = 0.01;
double laser_point_cov, acc_norm;
double vel_cov, acc_cov_input, gyr_cov_input;
double gyr_cov_output, acc_cov_output, b_gyr_cov, b_acc_cov;
double imu_meas_acc_cov, imu_meas_omg_cov; 
int    lidar_type, pcd_save_interval;
std::vector<double> gravity_init;
std::vector<double> extrinT, extrinT_gnss(3, 0.0), offline_init_vec(9, 0.0);
std::vector<double> extrinR, extrinR_gnss(9, 0.0);
bool   runtime_pos_log, pcd_save_en, path_en, extrinsic_est_en = true;
bool   scan_pub_en, scan_body_pub_en;
shared_ptr<Preprocess> p_pre;
shared_ptr<LI_Init> Init_LI;
shared_ptr<ImuProcess> p_imu;
shared_ptr<GNSSProcess> p_gnss;
double time_update_last = 0.0, time_current = 0.0, time_predict_last_const = 0.0, t_last = 0.0;
double time_diff_lidar_to_imu = 0.0;

std::string gnss_ephem_topic, gnss_glo_ephem_topic, gnss_meas_topic, gnss_iono_params_topic;
std::string gnss_tp_info_topic, local_trigger_info_topic, rtk_pvt_topic, rtk_lla_topic;
std::vector<double> default_gnss_iono_params(8, 0.0);
double gnss_local_time_diff = 18.0;
bool next_pulse_time_valid = false, update_gnss = false;
bool time_diff_valid = false, is_first_gnss = true;
double latest_gnss_time = -1, next_pulse_time = 0.0; 
double time_diff_gnss_local = 0.0;
bool gnss_local_online_sync = true, nolidar = false; 
double li_init_gyr_cov = 0.1, li_init_acc_cov = 0.1, lidar_time_inte = 0.1, first_imu_time = 0.0;
int cut_frame_num = 1, orig_odom_freq = 10;
double online_refine_time = 20.0; //unit: s
bool cut_frame_init = true;
bool GNSS_ENABLE = true;
Eigen::Matrix3d Rot_gnss_init(Eye3d);

MeasureGroup Measures;

ofstream fout_out, fout_imu_pbp, fout_rtk;

void readParameters(ros::NodeHandle &nh)
{
  p_pre.reset(new Preprocess());
  Init_LI.reset(new LI_Init());
  p_imu.reset(new ImuProcess());
  p_gnss.reset(new GNSSProcess());
  nh.param<bool>("prop_at_freq_of_imu", prop_at_freq_of_imu, 1);
  nh.param<bool>("use_imu_as_input", use_imu_as_input, 1);
  nh.param<bool>("check_satu", check_satu, 1);
  nh.param<int>("init_map_size", init_map_size, 100);
  nh.param<bool>("space_down_sample", space_down_sample, 1);
  nh.param<double>("mapping/satu_acc",satu_acc,3.0);
  nh.param<double>("mapping/satu_gyro",satu_gyro,35.0);
  nh.param<double>("mapping/acc_norm",acc_norm,1.0);
  nh.param<float>("mapping/plane_thr", plane_thr, 0.05f);
  nh.param<int>("point_filter_num", p_pre->point_filter_num, 2);
  nh.param<std::string>("common/lid_topic",lid_topic,"/livox/lidar");
  nh.param<std::string>("common/imu_topic", imu_topic,"/livox/imu");
  nh.param<bool>("common/con_frame",con_frame,false);
  nh.param<int>("common/con_frame_num",con_frame_num,1);
  nh.param<bool>("common/cut_frame",cut_frame,false);
  nh.param<double>("common/cut_frame_time_interval",cut_frame_time_interval,0.1);
  nh.param<double>("common/time_diff_lidar_to_imu",time_diff_lidar_to_imu,0.0);
  nh.param<double>("filter_size_surf",filter_size_surf_min,0.5);
  nh.param<double>("filter_size_map",filter_size_map_min,0.5);
  nh.param<double>("cube_side_length",cube_len,200);
  nh.param<float>("mapping/det_range",DET_RANGE,300.f);
  nh.param<double>("mapping/fov_degree",fov_deg,180);
  nh.param<bool>("mapping/imu_en",imu_en,true);
  nh.param<bool>("mapping/extrinsic_est_en",extrinsic_est_en,true);
  nh.param<double>("mapping/imu_time_inte",imu_time_inte,0.005);
  nh.param<double>("mapping/lidar_meas_cov",laser_point_cov,0.1);
  nh.param<double>("mapping/acc_cov_input",acc_cov_input,0.1);
  nh.param<double>("mapping/vel_cov",vel_cov,20);
  nh.param<double>("mapping/gyr_cov_input",gyr_cov_input,0.1);
  nh.param<double>("mapping/gyr_cov_output",gyr_cov_output,0.1);
  nh.param<double>("mapping/acc_cov_output",acc_cov_output,0.1);
  nh.param<double>("mapping/b_gyr_cov",b_gyr_cov,0.0001);
  nh.param<double>("mapping/b_acc_cov",b_acc_cov,0.0001);
  nh.param<double>("mapping/imu_meas_acc_cov",imu_meas_acc_cov,0.1);
  nh.param<double>("mapping/imu_meas_omg_cov",imu_meas_omg_cov,0.1);
  nh.param<double>("preprocess/blind", p_pre->blind, 1.0);
  nh.param<int>("preprocess/lidar_type", lidar_type, 1);
  nh.param<int>("preprocess/scan_line", p_pre->N_SCANS, 16);
  nh.param<int>("preprocess/scan_rate", p_pre->SCAN_RATE, 10);
  nh.param<int>("preprocess/timestamp_unit", p_pre->time_unit, 1);
  nh.param<double>("mapping/match_s", match_s, 81);
  nh.param<std::vector<double>>("mapping/gravity", gravity_init, std::vector<double>());
  nh.param<std::vector<double>>("mapping/extrinsic_T", extrinT, std::vector<double>());
  nh.param<std::vector<double>>("mapping/extrinsic_R", extrinR, std::vector<double>());
  nh.param<bool>("odometry/publish_odometry_without_downsample", publish_odometry_without_downsample, false);
  nh.param<bool>("publish/path_en",path_en, true);
  nh.param<bool>("publish/scan_publish_en",scan_pub_en,1);
  nh.param<bool>("publish/scan_bodyframe_pub_en",scan_body_pub_en,1);
  nh.param<bool>("runtime_pos_log_enable", runtime_pos_log, 0);
  nh.param<bool>("pcd_save/pcd_save_en", pcd_save_en, false);
  nh.param<int>("pcd_save/interval", pcd_save_interval, -1);

  nh.param<double>("mapping/lidar_time_inte",lidar_time_inte,0.1);
  nh.param<double>("mapping/lidar_meas_cov",laser_point_cov,0.1);
  nh.param<double>("gnss/psr_dopp_weight",p_gnss->relative_sqrt_info, 10);
  nh.param<double>("gnss/cp_weight",p_gnss->cp_weight, 0.1);
  // nh.param<double>("gnss/odo_weight",p_gnss->odo_weight, 0.1);
  nh.param<double>("gnss/gnss_ekf_noise",gnss_ekf_noise,0.01);
  nh.param<vector<double>>("gnss/gnss_extrinsic_T", extrinT_gnss, vector<double>());
  nh.param<vector<double>>("gnss/gnss_extrinsic_R", extrinR_gnss, vector<double>());
  nh.param<vector<double>>("gnss/offline_init_vec", offline_init_vec, vector<double>());

    nh.param<bool>("initialization/LIInit_en", p_imu->UseLIInit, false);
    if (p_imu->UseLIInit)
    {
        nh.param<double>("initialization/acc_cov",li_init_acc_cov,0.1);
        nh.param<double>("initialization/gyr_cov",li_init_gyr_cov,0.1);
        nh.param<bool>("initialization/cut_frame", cut_frame_init, true);
        nh.param<int>("initialization/cut_frame_num", cut_frame_num, 1);
        nh.param<int>("initialization/orig_odom_freq", orig_odom_freq, 10);
        nh.param<double>("initialization/online_refine_time", online_refine_time, 20.0);
        // nh.param<double>("initialization/mean_acc_norm", mean_acc_norm, 9.81);
        nh.param<double>("initialization/data_accum_length", Init_LI->data_accum_length, 300);
        p_imu->LI_init_done = false;
    }
    p_imu->gravity_ << VEC_FROM_ARRAY(gravity_init);
    nh.param<bool>("gnss/gnss_enable", GNSS_ENABLE, false);
    nh.param<bool>("gnss/outlier_rejection", p_gnss->p_assign->outlier_rej, false);
    cout << "gnss enable:" << GNSS_ENABLE << endl;
    if (GNSS_ENABLE)
    {
        nh.param<string>("gnss/gnss_ephem_topic",gnss_ephem_topic,"/ublox_driver/ephem");
        nh.param<string>("gnss/gnss_glo_ephem_topic",gnss_glo_ephem_topic,"/ublox_driver/glo_ephem");
        nh.param<string>("gnss/gnss_meas_topic",gnss_meas_topic,"/ublox_driver/range_meas");
        nh.param<string>("gnss/gnss_iono_params_topic",gnss_iono_params_topic,"/ublox_driver/iono_params");
        nh.param<string>("gnss/rtk_pvt_topic",rtk_pvt_topic,"/ublox_driver/receiver_pvt");
        nh.param<string>("gnss/rtk_lla_topic",rtk_lla_topic,"/ublox_driver/receiver_lla");
        nh.param<string>("gnss/gnss_tp_info_topic",gnss_tp_info_topic,"/ublox_driver/time_pulse_info");
        nh.param<vector<double>>("gnss/gnss_iono_default_parameters",default_gnss_iono_params,vector<double>());
        p_gnss->gravity_init << VEC_FROM_ARRAY(gravity_init);
        nh.param<bool>("gnss/gnss_local_online_sync",gnss_local_online_sync,1);
        if (gnss_local_online_sync)
        {
            nh.param<string>("gnss/local_trigger_info_topic",local_trigger_info_topic,"/external_trigger");
        }
        else
        {
            nh.param<double>("gnss/gnss_local_time_diff",gnss_local_time_diff, 18.0);
        }
        nh.param<double>("gnss/gnss_elevation_thres",p_gnss->p_assign->gnss_elevation_threshold, 30.0);
        nh.param<double>("gnss/prior_noise",p_gnss->p_assign->prior_noise, 0.010);
        nh.param<double>("gnss/marg_noise",p_gnss->p_assign->marg_noise, 0.010);
        nh.param<double>("gnss/b_acc_noise",p_gnss->pre_integration->acc_w, 0.010);
        nh.param<double>("gnss/b_omg_noise",p_gnss->pre_integration->gyr_w, 0.010);
        nh.param<double>("gnss/acc_noise",p_gnss->pre_integration->acc_n, 0.010);
        nh.param<double>("gnss/omg_noise",p_gnss->pre_integration->gyr_n, 0.010);
        nh.param<double>("gnss/ddt_noise",p_gnss->p_assign->ddt_noise, 0.010);
        nh.param<double>("gnss/dt_noise",p_gnss->p_assign->dt_noise, 0.010);
        nh.param<double>("gnss/psr_dopp_noise",p_gnss->p_assign->psr_dopp_noise,0.1);
        nh.param<double>("gnss/odo_noise",p_gnss->p_assign->odo_noise,0.1);
        nh.param<double>("gnss/cp_noise",p_gnss->p_assign->cp_noise,0.1);
        nh.param<double>("gnss/gnss_psr_std_thres",p_gnss->p_assign->gnss_psr_std_threshold, 2.0);
        nh.param<double>("gnss/gnss_dopp_std_thres",p_gnss->p_assign->gnss_dopp_std_threshold, 2.0);
        nh.param<double>("gnss/gnss_cp_std_thres",p_gnss->gnss_cp_std_threshold, 2.0);
        nh.param<double>("gnss/gnss_cp_time_thres",p_gnss->gnss_cp_time_threshold, 2.0);
        nh.param<int>("gnss/gnss_track_num_thres",p_gnss->p_assign->gnss_track_num_threshold, 20);
        nh.param<int>("gnss/gtsam_variable_thres",p_gnss->delete_thred, 200);
        nh.param<int>("gnss/gtsam_marg_variable_thres",p_gnss->p_assign->marg_thred, 1);
        nh.param<double>("gnss/outlier_thres",p_gnss->p_assign->outlier_thres, 0.1);
        nh.param<double>("gnss/gnss_sample_period",p_gnss->gnss_sample_period, 0.1);
        nh.param<bool>("gnss/nolidar",nolidar, false);
        nh.param<bool>("gnss/quick_init",p_gnss->quick_init, false);
        nh.param<bool>("gnss/online_init",p_gnss->gnss_online_init, false);
        nh.param<bool>("gnss/ephem_from_rinex",p_gnss->p_assign->ephem_from_rinex, false);
        p_gnss->p_assign->initNoises();
        if (nolidar)
        {
            nh.param<int>("gnss/window_size",p_gnss->wind_size, 1);
            p_imu->UseLIInit = false;
            p_imu->LI_init_done = true;
        }
    }
    else
    {
        nh.param<string>("gnss/rtk_pvt_topic",rtk_pvt_topic,"/ublox_driver/receiver_pvt");
    }
}

vect3 SO3ToEuler(const SO3 &orient) 
{
	Eigen::Matrix<double, 3, 1> _ang;
	Eigen::Vector4d q_data = orient.coeffs().transpose();
	//scalar w=orient.coeffs[3], x=orient.coeffs[0], y=orient.coeffs[1], z=orient.coeffs[2];
	double sqw = q_data[3]*q_data[3];
	double sqx = q_data[0]*q_data[0];
	double sqy = q_data[1]*q_data[1];
	double sqz = q_data[2]*q_data[2];
	double unit = sqx + sqy + sqz + sqw; // if normalized is one, otherwise is correction factor
	double test = q_data[3]*q_data[1] - q_data[2]*q_data[0];

	if (test > 0.49999*unit) { // singularity at north pole
	
		_ang << 2 * std::atan2(q_data[0], q_data[3]), M_PI/2, 0;
		double temp[3] = {_ang[0] * 57.3, _ang[1] * 57.3, _ang[2] * 57.3};
		vect3 euler_ang(temp, 3);
		return euler_ang;
	}
	if (test < -0.49999*unit) { // singularity at south pole
		_ang << -2 * std::atan2(q_data[0], q_data[3]), -M_PI/2, 0;
		double temp[3] = {_ang[0] * 57.3, _ang[1] * 57.3, _ang[2] * 57.3};
		vect3 euler_ang(temp, 3);
		return euler_ang;
	}
		
	_ang <<
			std::atan2(2*q_data[0]*q_data[3]+2*q_data[1]*q_data[2] , -sqx - sqy + sqz + sqw),
			std::asin (2*test/unit),
			std::atan2(2*q_data[2]*q_data[3]+2*q_data[1]*q_data[0] , sqx - sqy - sqz + sqw);
	double temp[3] = {_ang[0] * 57.3, _ang[1] * 57.3, _ang[2] * 57.3};
	vect3 euler_ang(temp, 3);
	return euler_ang;
}

void open_file()
{

    fout_out.open(DEBUG_FILE_DIR("mat_out.txt"),ios::out);
    fout_imu_pbp.open(DEBUG_FILE_DIR("imu_pbp.txt"),ios::out);
    if (GNSS_ENABLE)
    {fout_rtk.open(DEBUG_FILE_DIR("pos_rtk.txt"),ios::out);}
    if (fout_out && fout_imu_pbp)
        cout << "~~~~"<<ROOT_DIR<<" file opened" << endl;
    else
        cout << "~~~~"<<ROOT_DIR<<" doesn't exist" << endl;

}

void set_gnss_offline_init(bool nolidar_)
{
    if (!nolidar_)
    {
        p_gnss->gnss_ready = true;
        // double init_value[8];
        // init_value[0] = offline_init_vec[5]; init_value[1] = offline_init_vec[6]; init_value[2] = offline_init_vec[7]; init_value[3] = offline_init_vec[8];
        // init_value[4] = offline_init_vec[4]; init_value[5] = offline_init_vec[0]; init_value[6] = offline_init_vec[1]; init_value[7] = offline_init_vec[2];

        Eigen::Matrix3d R_enu_local_, rot_init;
        Eigen::Vector3d ecef_pos;
        ecef_pos << offline_init_vec[0], offline_init_vec[1], offline_init_vec[2];
        R_enu_local_ = ecef2rotation(ecef_pos) * Eigen::AngleAxisd(offline_init_vec[3], Eigen::Vector3d::UnitZ()) * Rot_gnss_init;  
        // p_gnss->gtSAMgraph.add(glio::PriorFactor(B(0), C(0), E(0), P(0), init_value, R_enu_local_, p_gnss->priorNoise));
        // p_gnss->time_frame.push_back(std::pair<double, int>(0.0, 0));

        Eigen::Matrix<double, 6, 1> init_vel_bias_vector;
        if (use_imu_as_input)
        {
            init_vel_bias_vector.block<3,1>(0,0) = kf_input.x_.pos;
            init_vel_bias_vector.block<3,1>(3,0) = kf_input.x_.vel;
            rot_init = kf_input.x_.rot.normalized().toRotationMatrix();
        }
        else
        {
            init_vel_bias_vector.block<3,1>(0,0) = kf_output.x_.pos;
            init_vel_bias_vector.block<3,1>(3,0) = kf_output.x_.vel;
            rot_init = kf_output.x_.rot.normalized().toRotationMatrix();
        }
        // init_vel_bias_vector.block<6,1>(6,0) = Eigen::Matrix<double, 6, 1>::Zero();
        gtsam::PriorFactor<gtsam::Rot3> init_rot_(R(0), gtsam::Rot3(rot_init), p_gnss->p_assign->priorrotNoise);
        gtsam::PriorFactor<gtsam::Vector6> init_vel_(A(0), gtsam::Vector6(init_vel_bias_vector), p_gnss->p_assign->priorNoise);
        // gtsam::PriorFactor<gtsam::Vector12> init_vel_(F(0), gtsam::Vector12(init_vel_bias_vector), p_gnss->priorposNoise);
        gtsam::PriorFactor<gtsam::Rot3> init_rot_ext(P(0), gtsam::Rot3(R_enu_local_), p_gnss->p_assign->priorrotNoise);
        gtsam::PriorFactor<gtsam::Vector3> init_pos_ext(E(0), gtsam::Vector3(offline_init_vec[0], offline_init_vec[1], offline_init_vec[2]), p_gnss->p_assign->margNoise);
        gtsam::PriorFactor<gtsam::Vector4> init_dt_(B(0), gtsam::Vector4(offline_init_vec[5], offline_init_vec[6], offline_init_vec[7], offline_init_vec[8]), p_gnss->p_assign->priordtNoise);
        gtsam::PriorFactor<gtsam::Vector1> init_ddt_(C(0), gtsam::Vector1(offline_init_vec[4]), p_gnss->p_assign->priorddtNoise);
        p_gnss->p_assign->gtSAMgraph.add(init_rot_);
        p_gnss->p_assign->gtSAMgraph.add(init_vel_);
        p_gnss->p_assign->gtSAMgraph.add(init_rot_ext);
        p_gnss->p_assign->gtSAMgraph.add(init_pos_ext);
        p_gnss->p_assign->gtSAMgraph.add(init_dt_);
        p_gnss->p_assign->gtSAMgraph.add(init_ddt_);
        p_gnss->p_assign->factor_id_frame.push_back(std::vector<size_t>{0, 1, 2, 3, 4, 5});

        p_gnss->p_assign->initialEstimate.insert(E(0), gtsam::Vector3(offline_init_vec[0], offline_init_vec[1], offline_init_vec[2]));
        p_gnss->p_assign->initialEstimate.insert(P(0), gtsam::Rot3(R_enu_local_));
        p_gnss->p_assign->initialEstimate.insert(R(0), gtsam::Rot3(rot_init));
        p_gnss->p_assign->initialEstimate.insert(A(0), gtsam::Vector6((init_vel_bias_vector)));
        // p_gnss->initialEstimate.insert(F(0), gtsam::Vector12((init_vel_bias_vector)));
        p_gnss->p_assign->initialEstimate.insert(B(0), gtsam::Vector4(offline_init_vec[5], offline_init_vec[6], offline_init_vec[7], offline_init_vec[8]));
        p_gnss->p_assign->initialEstimate.insert(C(0), gtsam::Vector1(offline_init_vec[4]));
        p_gnss->id_accumulate = 6; // 54;
    
        p_gnss->last_gnss_time = Measures.lidar_beg_time; // imu_first_time;
        p_gnss->state_ = kf_input.x_;
        p_gnss->state_const_ = kf_output.x_;
        p_gnss->frame_num = 1; // 11;
        p_gnss->runISAM2opt();
    }
    else
    {
        p_gnss->gnss_ready = true;
        is_first_frame = false;
        Eigen::Matrix<double, 12, 1> init_vel_bias_vector;

        if (use_imu_as_input)
        {
            state_in.pos(0) = offline_init_vec[0];
            state_in.pos(1) = offline_init_vec[1];
            state_in.pos(2) = offline_init_vec[2];
            Eigen::Matrix3d R_enu_local_, R_ecef_enu;
            R_ecef_enu = ecef2rotation(state_in.pos);
            R_enu_local_ = Eigen::AngleAxisd(offline_init_vec[3], Eigen::Vector3d::UnitZ());
            state_in.rot = R_ecef_enu * R_enu_local_;
            state_in.rot.normalize();
            state_in.pos -= state_in.rot.normalized().toRotationMatrix() * p_gnss->Tex_imu_r;
            state_in.gravity = R_ecef_enu * kf_input.x_.gravity; // * R_enu_local_ 
            kf_input.change_x(state_in);
        
            gtsam::PriorFactor<gtsam::Rot3> init_rot(R(0), gtsam::Rot3(state_in.rot.normalized().toRotationMatrix()), p_gnss->p_assign->priorposNoise);
            gtsam::PriorFactor<gtsam::Vector4> init_dt(B(0), gtsam::Vector4(offline_init_vec[5], offline_init_vec[6], offline_init_vec[7], offline_init_vec[8]), p_gnss->p_assign->priordtNoise);
            gtsam::PriorFactor<gtsam::Vector1> init_ddt(C(0), gtsam::Vector1(offline_init_vec[4]), p_gnss->p_assign->priorddtNoise);
            // Eigen::Matrix<double, 12, 1> init_vel_bias_vector;
            init_vel_bias_vector.block<3,1>(0,0) = state_in.pos;
            init_vel_bias_vector.block<3,1>(3,0) = state_in.vel;
            init_vel_bias_vector.block<3,1>(6,0) = state_in.ba;
            init_vel_bias_vector.block<3,1>(9,0) = state_in.bg;
            gtsam::PriorFactor<gtsam::Vector12> init_vel_bias(F(0), gtsam::Vector12(init_vel_bias_vector), p_gnss->p_assign->priorposNoise);
            p_gnss->p_assign->initialEstimate.insert(R(0), gtsam::Rot3(state_in.rot.normalized().toRotationMatrix()));
            p_gnss->state_ = state_in;

            p_gnss->p_assign->gtSAMgraph.add(init_rot);
            p_gnss->p_assign->gtSAMgraph.add(init_vel_bias);
            p_gnss->p_assign->gtSAMgraph.add(init_dt);
            p_gnss->p_assign->gtSAMgraph.add(init_ddt);
        }
        else
        {
            state_out.pos(0) = offline_init_vec[0];
            state_out.pos(1) = offline_init_vec[1];
            state_out.pos(2) = offline_init_vec[2];
            Eigen::Matrix3d R_enu_local_, R_ecef_enu;
            R_ecef_enu = ecef2rotation(state_out.pos);
            R_enu_local_ = Eigen::AngleAxisd(offline_init_vec[3], Eigen::Vector3d::UnitZ());
            state_out.rot = R_ecef_enu * R_enu_local_;
            state_out.rot.normalize();
            state_out.pos -= state_out.rot.normalized().toRotationMatrix() * p_gnss->Tex_imu_r;
            state_out.gravity = R_ecef_enu * kf_output.x_.gravity; // * R_enu_local_ 
            kf_output.change_x(state_out);
        
            gtsam::PriorFactor<gtsam::Rot3> init_rot(R(0), gtsam::Rot3(state_in.rot.normalized().toRotationMatrix()), p_gnss->p_assign->priorposNoise);
            gtsam::PriorFactor<gtsam::Vector4> init_dt(B(0), gtsam::Vector4(offline_init_vec[5], offline_init_vec[6], offline_init_vec[7], offline_init_vec[8]), p_gnss->p_assign->priordtNoise);
            gtsam::PriorFactor<gtsam::Vector1> init_ddt(C(0), gtsam::Vector1(offline_init_vec[4]), p_gnss->p_assign->priorddtNoise);
            init_vel_bias_vector.block<3,1>(0,0) = state_out.pos;
            init_vel_bias_vector.block<3,1>(3,0) = state_out.vel;
            init_vel_bias_vector.block<3,1>(6,0) = state_out.ba;
            init_vel_bias_vector.block<3,1>(9,0) = state_out.bg;
            gtsam::PriorFactor<gtsam::Vector12> init_vel_bias(F(0), gtsam::Vector12(init_vel_bias_vector), p_gnss->p_assign->priorposNoise);
            p_gnss->p_assign->initialEstimate.insert(R(0), gtsam::Rot3(state_out.rot.normalized().toRotationMatrix()));
            p_gnss->state_const_ = state_out;
        
            p_gnss->p_assign->gtSAMgraph.add(init_rot);
            p_gnss->p_assign->gtSAMgraph.add(init_vel_bias);
            p_gnss->p_assign->gtSAMgraph.add(init_dt);
            p_gnss->p_assign->gtSAMgraph.add(init_ddt);
        }
        p_gnss->p_assign->factor_id_frame.push_back(std::vector<size_t>{0, 1, 2, 3});
        // p_gnss->time_frame.push_back(std::pair<double, int>(0.0, 0));
        
        p_gnss->p_assign->initialEstimate.insert(F(0), gtsam::Vector12(init_vel_bias_vector));
        p_gnss->p_assign->initialEstimate.insert(B(0), gtsam::Vector4(offline_init_vec[5], offline_init_vec[6], offline_init_vec[7], offline_init_vec[8]));
        p_gnss->p_assign->initialEstimate.insert(C(0), gtsam::Vector1(offline_init_vec[4]));
        p_gnss->id_accumulate = 4; // 54;
        p_gnss->last_gnss_time = time_current;
        p_gnss->frame_num = 1; // 11;
        p_gnss->runISAM2opt();
    }
}

void cout_state_to_file()
{
    if (use_imu_as_input)
    {
        Eigen::Vector3d pos_enu;
        if (!nolidar)
        {
            Eigen::Vector3d pos_r = p_gnss->state_.rot * p_gnss->Tex_imu_r + p_gnss->state_.pos;
            Eigen::Matrix3d enu_rot = p_gnss->p_assign->isamCurrentEstimate.at<gtsam::Rot3>(P(0)).matrix();
            Eigen::Vector3d anc_cur = p_gnss->p_assign->isamCurrentEstimate.at<gtsam::Vector3>(E(0));
            pos_enu = p_gnss->local2enu(enu_rot, anc_cur, pos_r);
        }
        else
        {
            Eigen::Vector3d pos_r = kf_input.x_.rot.normalized().toRotationMatrix() * p_gnss->Tex_imu_r + kf_input.x_.pos;
            pos_enu = p_gnss->local2enu(Eigen::Matrix3d::Zero(), Eigen::Vector3d::Zero(), pos_r);
        }
        if (!p_gnss->gnss_online_init)
        {
            V3D euler_cur = SO3ToEuler(kf_input.x_.rot);
            fout_out << setw(20) << t_last - first_imu_time << " " << euler_cur.transpose() << " " << pos_enu.transpose() << " " << kf_input.x_.vel.transpose() \
                        <<" "<<kf_input.x_.gravity.transpose()<<" "<<p_gnss->p_assign->isamCurrentEstimate.at<gtsam::Vector1>(C(p_gnss->frame_num-1))[0] 
                        << " " <<p_gnss->p_assign->isamCurrentEstimate.at<gtsam::Vector1>(C(p_gnss->frame_num-1))[0]<< " " <<p_gnss->p_assign->isamCurrentEstimate.at<gtsam::Vector1>(C(p_gnss->frame_num-1))[0]<<" "<<p_gnss->p_assign->isamCurrentEstimate.at<gtsam::Vector4>(B(p_gnss->frame_num-1)).transpose()<<endl;
        }
        else
        {
            if (nolidar)
            {
                V3D euler_cur = SO3ToEuler(kf_input.x_.rot);
                fout_out << setw(20) << t_last - first_imu_time << " " << euler_cur.transpose() << " " << pos_enu.transpose() << " " << kf_input.x_.vel.transpose() \
                        <<" "<<kf_input.x_.gravity.transpose()<<" "<<p_gnss->p_assign->isamCurrentEstimate.at<gtsam::Vector1>(C(p_gnss->frame_num-1))[0]  
                        << " " << p_gnss->p_assign->isamCurrentEstimate.at<gtsam::Vector1>(C(p_gnss->frame_num-1))[0] << " " <<p_gnss->p_assign->isamCurrentEstimate.at<gtsam::Vector1>(C(p_gnss->frame_num-1))[0]<<" "<<p_gnss->p_assign->isamCurrentEstimate.at<gtsam::Vector4>(B(p_gnss->frame_num-1)).transpose()<<endl;
            }
            else
            {
                V3D euler_cur = SO3ToEuler(kf_input.x_.rot); // euler_cur.transpose()
                Eigen::Vector3d euler_ext = RotMtoEuler(p_gnss->p_assign->isamCurrentEstimate.at<gtsam::Rot3>(P(0)).matrix()); // euler_cur.transpose()
                fout_out << setw(20) << t_last - first_imu_time << " " << euler_cur.transpose() << " " << pos_enu.transpose() << " " << kf_input.x_.pos.transpose() \
                            <<" "<<p_gnss->p_assign->isamCurrentEstimate.at<gtsam::Vector3>(E(0)).transpose()<<" "<<  
                            euler_ext.segment<2>(0).transpose() << " " <<p_gnss->p_assign->isamCurrentEstimate.at<gtsam::Vector1>(C(p_gnss->frame_num-1))[0]<<" "
                            <<p_gnss->p_assign->isamCurrentEstimate.at<gtsam::Vector4>(B(p_gnss->frame_num-1)).transpose()<<endl;
            }
        }
    }
    else
    {
        Eigen::Vector3d pos_enu;
        if (!nolidar)
        {
            Eigen::Vector3d pos_r = p_gnss->state_const_.rot * p_gnss->Tex_imu_r + p_gnss->state_const_.pos; // maybe improper
            Eigen::Matrix3d enu_rot = p_gnss->p_assign->isamCurrentEstimate.at<gtsam::Rot3>(P(0)).matrix();
            Eigen::Vector3d anc_cur = p_gnss->p_assign->isamCurrentEstimate.at<gtsam::Vector3>(E(0));
            pos_enu = p_gnss->local2enu(enu_rot, anc_cur, pos_r);
        }
        else
        {
            Eigen::Vector3d pos_r = kf_output.x_.rot.normalized().toRotationMatrix() * p_gnss->Tex_imu_r + kf_output.x_.pos;
            pos_enu = p_gnss->local2enu(Eigen::Matrix3d::Zero(), Eigen::Vector3d::Zero(), pos_r);
        }
        if (!p_gnss->gnss_online_init)
        {
            Eigen::Vector3d euler_cur = SO3ToEuler(kf_output.x_.rot);
            fout_out << setw(20) << time_predict_last_const - first_imu_time << " " << euler_cur.transpose() << " " << pos_enu.transpose() << " " << kf_output.x_.vel.transpose() \
                        <<" "<<kf_output.x_.omg.transpose()<<" "<<kf_output.x_.acc.transpose()<<" "<<kf_output.x_.gravity.transpose()<<" "<<p_gnss->p_assign->isamCurrentEstimate.at<gtsam::Vector1>(C(p_gnss->frame_num-1))[0] \
                        << " " << p_gnss->p_assign->isamCurrentEstimate.at<gtsam::Vector1>(C(p_gnss->frame_num-1))[0] << " " <<p_gnss->p_assign->isamCurrentEstimate.at<gtsam::Vector1>(C(p_gnss->frame_num-1))[0]<<" "<<p_gnss->p_assign->isamCurrentEstimate.at<gtsam::Vector4>(B(p_gnss->frame_num-1)).transpose()<<endl;
        }
        else
        {
            if (nolidar)
            {
                V3D euler_cur = SO3ToEuler(kf_input.x_.rot);
                fout_out << setw(20) << time_predict_last_const - first_imu_time << " " << euler_cur.transpose() << " " << pos_enu.transpose() << " " << kf_output.x_.vel.transpose() \
                        <<" "<<kf_output.x_.omg.transpose()<<" "<<kf_output.x_.acc.transpose()<<" "<<kf_output.x_.gravity.transpose()<<" "<<p_gnss->p_assign->isamCurrentEstimate.at<gtsam::Vector1>(C(p_gnss->frame_num-1))[0] \
                        << " " << p_gnss->p_assign->isamCurrentEstimate.at<gtsam::Vector1>(C(p_gnss->frame_num-1))[0] << " " <<p_gnss->p_assign->isamCurrentEstimate.at<gtsam::Vector1>(C(p_gnss->frame_num-1))[0]<<" "<<p_gnss->p_assign->isamCurrentEstimate.at<gtsam::Vector4>(B(p_gnss->frame_num-1)).transpose()<<endl;
            }
            else
            {
                Eigen::Vector3d euler_cur = SO3ToEuler(kf_output.x_.rot); // euler_cur.transpose()
                Eigen::Vector3d euler_ext = SO3ToEuler(p_gnss->p_assign->isamCurrentEstimate.at<gtsam::Rot3>(P(0)).matrix());
                fout_out << setw(20) << time_predict_last_const - first_imu_time << " " << euler_cur.transpose() << " " << pos_enu.transpose() << " " << kf_output.x_.pos.transpose() \
                            <<" "<<kf_output.x_.omg.transpose()<<" "<<kf_output.x_.acc.transpose()<<" "<<p_gnss->p_assign->isamCurrentEstimate.at<gtsam::Vector3>(E(0)).transpose()<<" "
                            << euler_ext.segment<2>(0).transpose() << " " <<p_gnss->p_assign->isamCurrentEstimate.at<gtsam::Vector1>(C(p_gnss->frame_num-1))[0]<<" "<<p_gnss->p_assign->isamCurrentEstimate.at<gtsam::Vector4>(B(p_gnss->frame_num-1)).transpose()<<endl;
            }
        }
    }
}