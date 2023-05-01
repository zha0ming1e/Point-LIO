// #include <../include/IKFoM/IKFoM_toolkit/esekfom/esekfom.hpp>
#include "Estimator.h"

PointCloudXYZI::Ptr normvec(new PointCloudXYZI(100000, 1));
std::vector<int> time_seq;
PointCloudXYZI::Ptr feats_down_body(new PointCloudXYZI(10000, 1));
PointCloudXYZI::Ptr feats_down_world(new PointCloudXYZI(10000, 1));
std::vector<V3D> pbody_list;
std::vector<PointVector> Nearest_Points; 
std::shared_ptr<IVoxType> ivox_ = nullptr;                    // localmap in ivox
std::vector<float> pointSearchSqDis(NUM_MATCH_POINTS);
bool   point_selected_surf[100000] = {0};
std::vector<M3D> crossmat_list;
int effct_feat_num = 0;
int k = 0;
int idx = -1;
esekfom::esekf<state_input, 18, input_ikfom> kf_input;
esekfom::esekf<state_output, 24, input_ikfom> kf_output;
input_ikfom input_in;
V3D angvel_avr, acc_avr, acc_avr_norm;
int feats_down_size = 0;  
V3D Lidar_T_wrt_IMU(Zero3d);
M3D Lidar_R_wrt_IMU(Eye3d);

Eigen::Matrix<double, 18, 18> process_noise_cov_input()
{
	Eigen::Matrix<double, 18, 18> cov;
	cov.setZero();
	cov.block<3, 3>(3, 3).diagonal() << gyr_cov_input, gyr_cov_input, gyr_cov_input;
	cov.block<3, 3>(6, 6).diagonal() << acc_cov_input, acc_cov_input, acc_cov_input;
	cov.block<3, 3>(9, 9).diagonal() << b_gyr_cov, b_gyr_cov, b_gyr_cov;
	cov.block<3, 3>(12, 12).diagonal() << b_acc_cov, b_acc_cov, b_acc_cov;
	// MTK::get_cov<process_noise_input>::type cov = MTK::get_cov<process_noise_input>::type::Zero();
	// MTK::setDiagonal<process_noise_input, vect3, 0>(cov, &process_noise_input::ng, gyr_cov_input);// 0.03
	// MTK::setDiagonal<process_noise_input, vect3, 3>(cov, &process_noise_input::na, acc_cov_input); // *dt 0.01 0.01 * dt * dt 0.05
	// MTK::setDiagonal<process_noise_input, vect3, 6>(cov, &process_noise_input::nbg, b_gyr_cov); // *dt 0.00001 0.00001 * dt *dt 0.3 //0.001 0.0001 0.01
	// MTK::setDiagonal<process_noise_input, vect3, 9>(cov, &process_noise_input::nba, b_acc_cov);   //0.001 0.05 0.0001/out 0.01
	return cov;
}

Eigen::Matrix<double, 24, 24> process_noise_cov_output()
{
	Eigen::Matrix<double, 24, 24> cov;
	cov.setZero();
	cov.block<3, 3>(6, 6).diagonal() << vel_cov, vel_cov, vel_cov;
	cov.block<3, 3>(9, 9).diagonal() << gyr_cov_output, gyr_cov_output, gyr_cov_output;
	cov.block<3, 3>(12, 12).diagonal() << acc_cov_output, acc_cov_output, acc_cov_output;
	cov.block<3, 3>(18, 18).diagonal() << b_gyr_cov, b_gyr_cov, b_gyr_cov;
	cov.block<3, 3>(21, 21).diagonal() << b_acc_cov, b_acc_cov, b_acc_cov;
	return cov;
}

Eigen::Matrix<double, 18, 1> get_f_input(state_input &s, const input_ikfom &in)
{
	Eigen::Matrix<double, 18, 1> res = Eigen::Matrix<double, 18, 1>::Zero();
	vect3 omega;
	in.gyro.boxminus(omega, s.bg);
	vect3 a_inertial = s.rot * (in.acc-s.ba); // .normalized()
	for(int i = 0; i < 3; i++ ){
		res(i) = s.vel[i];
		res(i + 3) = omega[i]; 
		res(i + 6) = a_inertial[i] + s.gravity[i]; 
	}
	return res;
}

Eigen::Matrix<double, 24, 1> get_f_output(state_output &s, const input_ikfom &in)
{
	Eigen::Matrix<double, 24, 1> res = Eigen::Matrix<double, 24, 1>::Zero();
	vect3 a_inertial = s.rot * s.acc; // .normalized()
	for(int i = 0; i < 3; i++ ){
		res(i) = s.vel[i];
		res(i + 3) = s.omg[i]; 
		res(i + 6) = a_inertial[i] + s.gravity[i]; 
	}
	return res;
}

Eigen::Matrix<double, 18, 18> df_dx_input(state_input &s, const input_ikfom &in)
{
	Eigen::Matrix<double, 18, 18> cov = Eigen::Matrix<double, 18, 18>::Zero();
	cov.template block<3, 3>(0, 6) = Eigen::Matrix3d::Identity();
	vect3 acc_;
	in.acc.boxminus(acc_, s.ba);
	vect3 omega;
	in.gyro.boxminus(omega, s.bg);
	cov.template block<3, 3>(6, 3) = -s.rot*MTK::hat(acc_); // .normalized().toRotationMatrix()
	cov.template block<3, 3>(6, 12) = -s.rot; //.normalized().toRotationMatrix();
	// Eigen::Matrix<state_ikfom::scalar, 2, 1> vec = Eigen::Matrix<state_ikfom::scalar, 2, 1>::Zero();
	// Eigen::Matrix<state_ikfom::scalar, 3, 2> grav_matrix;
	// s.S2_Mx(grav_matrix, vec, 21);
	cov.template block<3, 3>(6, 15) = Eigen::Matrix3d::Identity(); // grav_matrix; 
	cov.template block<3, 3>(3, 9) = -Eigen::Matrix3d::Identity(); 
	return cov;
}

Eigen::Matrix<double, 24, 24> df_dx_output(state_output &s, const input_ikfom &in)
{
	Eigen::Matrix<double, 24, 24> cov = Eigen::Matrix<double, 24, 24>::Zero();
	cov.template block<3, 3>(0, 6) = Eigen::Matrix3d::Identity();
	cov.template block<3, 3>(6, 3) = -s.rot*MTK::hat(s.acc); // .normalized().toRotationMatrix()
	cov.template block<3, 3>(6, 12) = s.rot; //.normalized().toRotationMatrix();
	// Eigen::Matrix<state_ikfom::scalar, 2, 1> vec = Eigen::Matrix<state_ikfom::scalar, 2, 1>::Zero();
	// Eigen::Matrix<state_ikfom::scalar, 3, 2> grav_matrix;
	// s.S2_Mx(grav_matrix, vec, 21);
	cov.template block<3, 3>(6, 15) = Eigen::Matrix3d::Identity(); // grav_matrix; 
	cov.template block<3, 3>(3, 9) = Eigen::Matrix3d::Identity(); 
	return cov;
}

void h_model_input(state_input &s, esekfom::dyn_share_modified<double> &ekfom_data)
{
	bool match_in_map = false;
	VF(4) pabcd;
	pabcd.setZero();
	normvec->resize(time_seq[k]);
	int effect_num_k = 0;
	for (int j = 0; j < time_seq[k]; j++)
	{
		PointType &point_body_j  = feats_down_body->points[idx+j+1];
		PointType &point_world_j = feats_down_world->points[idx+j+1];
		pointBodyToWorld(&point_body_j, &point_world_j); 
		V3D p_body = pbody_list[idx+j+1];
		V3D p_world;
		p_world << point_world_j.x, point_world_j.y, point_world_j.z;
		{
			auto &points_near = Nearest_Points[idx+j+1];
            ivox_->GetClosestPoint(point_world_j, points_near, NUM_MATCH_POINTS); // 
			if ((points_near.size() < NUM_MATCH_POINTS) || pointSearchSqDis[NUM_MATCH_POINTS - 1] > 5) // 5)
			{
				point_selected_surf[idx+j+1] = false;
			}
			else
			{
				point_selected_surf[idx+j+1] = false;
				if (esti_plane(pabcd, points_near, plane_thr)) //(planeValid)
				{
					float pd2 = pabcd(0) * point_world_j.x + pabcd(1) * point_world_j.y + pabcd(2) * point_world_j.z + pabcd(3);
					
					if (p_body.norm() > match_s * pd2 * pd2)
					{
						point_selected_surf[idx+j+1] = true;
						normvec->points[j].x = pabcd(0);
						normvec->points[j].y = pabcd(1);
						normvec->points[j].z = pabcd(2);
						normvec->points[j].intensity = pabcd(3);
						effect_num_k ++;
					}
				}  
			}
		}
	}
	if (effect_num_k == 0) 
	{
		ekfom_data.valid = false;
		return;
	}
	ekfom_data.M_Noise = laser_point_cov;
	ekfom_data.h_x.resize(effect_num_k, 6);
	ekfom_data.h_x = Eigen::MatrixXd::Zero(effect_num_k, 6); //12);
	ekfom_data.z.resize(effect_num_k);
	int m = 0;
	// V3D last_norm_vec = V3D::Zero();
	// if (!p_gnss->norm_vec_holder.empty()) 
	// {
	// 	last_norm_vec = p_gnss->norm_vec_holder.back();
	// }
	for (int j = 0; j < time_seq[k]; j++)
	{
		// ekfom_data.converge = false;
		if(point_selected_surf[idx+j+1])
		{
			V3D norm_vec(normvec->points[j].x, normvec->points[j].y, normvec->points[j].z);
			p_gnss->norm_vec_holder.push_back(norm_vec);
			// if (fabs(last_norm_vec.transpose() * norm_vec) > 0.9)
			// {
			// 	ekfom_data.converge = true;
			// }
			// if (extrinsic_est_en)
			// {
			// 	V3D p_body = pbody_list[idx+j+1];
			// 	M3D p_crossmat, p_imu_crossmat;
			// 	p_crossmat << SKEW_SYM_MATRX(p_body);
			// 	V3D point_imu = s.offset_R_L_I * p_body + s.offset_T_L_I;
			// 	p_imu_crossmat << SKEW_SYM_MATRX(point_imu);
			// 	V3D C(s.rot.conjugate().normalized() * norm_vec);
			// 	V3D A(p_imu_crossmat * C);
			// 	V3D B(p_crossmat * s.offset_R_L_I.conjugate() * C);
			// 	ekfom_data.h_x.block<1, 12>(m, 0) << norm_vec(0), norm_vec(1), norm_vec(2), VEC_FROM_ARRAY(A), VEC_FROM_ARRAY(B), VEC_FROM_ARRAY(C);
			// }
			// else
			{   
				M3D point_crossmat = crossmat_list[idx+j+1];
				V3D C(s.rot.transpose() * norm_vec); // conjugate().normalized()
				V3D A(point_crossmat * C);
				ekfom_data.h_x.block<1, 6>(m, 0) << norm_vec(0), norm_vec(1), norm_vec(2), VEC_FROM_ARRAY(A); //, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
			}
			ekfom_data.z(m) = -norm_vec(0) * feats_down_world->points[idx+j+1].x -norm_vec(1) * feats_down_world->points[idx+j+1].y -norm_vec(2) * feats_down_world->points[idx+j+1].z-normvec->points[j].intensity;
			// if (ekfom_data.converge)
			// {
			// 	ekfom_data.h_x.block<1, 6>(m, 0) *= 2;
			// 	ekfom_data.z(m) *= 2;
			// }
			m++;
		}
	}
	effct_feat_num += effect_num_k;
}

void h_model_output(state_output &s, esekfom::dyn_share_modified<double> &ekfom_data)
{
	bool match_in_map = false;
	VF(4) pabcd;
	pabcd.setZero();
	normvec->resize(time_seq[k]);
	int effect_num_k = 0;
	for (int j = 0; j < time_seq[k]; j++)
	{
		PointType &point_body_j  = feats_down_body->points[idx+j+1];
		PointType &point_world_j = feats_down_world->points[idx+j+1];
		pointBodyToWorld(&point_body_j, &point_world_j); 
		V3D p_body = pbody_list[idx+j+1];
		V3D p_world;
		p_world << point_world_j.x, point_world_j.y, point_world_j.z;
		{
			auto &points_near = Nearest_Points[idx+j+1];
			
            ivox_->GetClosestPoint(point_world_j, points_near, NUM_MATCH_POINTS); // 
			
			if ((points_near.size() < NUM_MATCH_POINTS) || pointSearchSqDis[NUM_MATCH_POINTS - 1] > 5)
			{
				point_selected_surf[idx+j+1] = false;
			}
			else
			{
				point_selected_surf[idx+j+1] = false;
				if (esti_plane(pabcd, points_near, plane_thr)) //(planeValid)
				{
					float pd2 = pabcd(0) * point_world_j.x + pabcd(1) * point_world_j.y + pabcd(2) * point_world_j.z + pabcd(3);
					
					if (p_body.norm() > match_s * pd2 * pd2)
					{
						// point_selected_surf[i] = true;
						point_selected_surf[idx+j+1] = true;
						normvec->points[j].x = pabcd(0);
						normvec->points[j].y = pabcd(1);
						normvec->points[j].z = pabcd(2);
						normvec->points[j].intensity = pabcd(3);
						effect_num_k ++;
					}
				}  
			}
		}
	}
	if (effect_num_k == 0) 
	{
		ekfom_data.valid = false;
		return;
	}
	ekfom_data.M_Noise = laser_point_cov;
	ekfom_data.h_x.resize(effect_num_k, 6);
	ekfom_data.h_x = Eigen::MatrixXd::Zero(effect_num_k, 6); // 12);
	ekfom_data.z.resize(effect_num_k);
	int m = 0;
	// V3D last_norm_vec = V3D::Zero();
	// if (!p_gnss->norm_vec_holder.empty()) 
	// {
	// 	last_norm_vec = p_gnss->norm_vec_holder.back();
	// }
	for (int j = 0; j < time_seq[k]; j++)
	{
		// ekfom_data.converge = false;
		if(point_selected_surf[idx+j+1])
		{
			V3D norm_vec(normvec->points[j].x, normvec->points[j].y, normvec->points[j].z);
			p_gnss->norm_vec_holder.push_back(norm_vec);
			// if (fabs(last_norm_vec.transpose() * norm_vec) > 0.9)
			// {
			// 	ekfom_data.converge = true;
			// }
			// if (extrinsic_est_en)
			// {
			// 	V3D p_body = pbody_list[idx+j+1];
			// 	M3D p_crossmat, p_imu_crossmat;
			// 	p_crossmat << SKEW_SYM_MATRX(p_body);
			// 	V3D point_imu = s.offset_R_L_I * p_body + s.offset_T_L_I;
			// 	p_imu_crossmat << SKEW_SYM_MATRX(point_imu);
			// 	V3D C(s.rot.conjugate().normalized() * norm_vec);
			// 	V3D A(p_imu_crossmat * C);
			// 	V3D B(p_crossmat * s.offset_R_L_I.conjugate() * C);
			// 	ekfom_data.h_x.block<1, 12>(m, 0) << norm_vec(0), norm_vec(1), norm_vec(2), VEC_FROM_ARRAY(A), VEC_FROM_ARRAY(B), VEC_FROM_ARRAY(C);
			// }
			// else
			{   
				M3D point_crossmat = crossmat_list[idx+j+1];
				V3D C(s.rot.transpose() * norm_vec); // conjugate().normalized()
				V3D A(point_crossmat * C);
				ekfom_data.h_x.block<1, 6>(m, 0) << norm_vec(0), norm_vec(1), norm_vec(2), VEC_FROM_ARRAY(A); //, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
			}
			ekfom_data.z(m) = -norm_vec(0) * feats_down_world->points[idx+j+1].x -norm_vec(1) * feats_down_world->points[idx+j+1].y -norm_vec(2) * feats_down_world->points[idx+j+1].z-normvec->points[j].intensity;
			// if (ekfom_data.converge)
			// {
			// 	ekfom_data.h_x.block<1, 6>(m, 0) *= 2;
			// 	ekfom_data.z(m) *= 2;
			// }
			m++;
		}
	}
	effct_feat_num += effect_num_k;
}

void h_model_IMU_output(state_output &s, esekfom::dyn_share_modified<double> &ekfom_data)
{
    std::memset(ekfom_data.satu_check, false, 6);
	ekfom_data.z_IMU.block<3,1>(0, 0) = angvel_avr - s.omg - s.bg;
	ekfom_data.z_IMU.block<3,1>(3, 0) = acc_avr * G_m_s2 / acc_norm - s.acc - s.ba;
    ekfom_data.R_IMU << imu_meas_omg_cov, imu_meas_omg_cov, imu_meas_omg_cov, imu_meas_acc_cov, imu_meas_acc_cov, imu_meas_acc_cov;
	if(check_satu)
	{
		if(fabs(angvel_avr(0)) >= 0.99 * satu_gyro)
		{
			ekfom_data.satu_check[0] = true; 
			ekfom_data.z_IMU(0) = 0.0;
		}
		
		if(fabs(angvel_avr(1)) >= 0.99 * satu_gyro) 
		{
			ekfom_data.satu_check[1] = true;
			ekfom_data.z_IMU(1) = 0.0;
		}
		
		if(fabs(angvel_avr(2)) >= 0.99 * satu_gyro)
		{
			ekfom_data.satu_check[2] = true;
			ekfom_data.z_IMU(2) = 0.0;
		}
		
		if(fabs(acc_avr(0)) >= 0.99 * satu_acc)
		{
			ekfom_data.satu_check[3] = true;
			ekfom_data.z_IMU(3) = 0.0;
		}

		if(fabs(acc_avr(1)) >= 0.99 * satu_acc) 
		{
			ekfom_data.satu_check[4] = true;
			ekfom_data.z_IMU(4) = 0.0;
		}

		if(fabs(acc_avr(2)) >= 0.99 * satu_acc) 
		{
			ekfom_data.satu_check[5] = true;
			ekfom_data.z_IMU(5) = 0.0;
		}
	}
}

void h_model_GNSS_input(state_input &s, esekfom::dyn_share_modified<double> &ekfom_data)
{
	// Eigen::Vector3d res_r;
	// p_gnss->state_.rot.boxminus(res_r, s.rot);
	ekfom_data.h_GNSS.setIdentity();
	ekfom_data.h_GNSS *= p_gnss->odo_weight;
	ekfom_data.z_GNSS.setZero();
	// ekfom_data.h_GNSS.block<3, 3>(3, 3) = Eigen::Matrix3d::Zero(); // Jacob_right_inv<double>(res_r); // 
	// ekfom_data.h_GNSS.block<3, 3>(6, 6) = Eigen::Matrix3d::Zero(); // 
	// ekfom_data.z_GNSS.block<3, 1>(3, 0) = res_r;
	ekfom_data.z_GNSS.block<3, 1>(0, 0) = p_gnss->odo_weight * (p_gnss->state_.pos - s.pos); // 
	// ekfom_data.z_GNSS.block<3, 1>(6, 0) = p_gnss->state_.vel - s.vel;
	ekfom_data.M_Noise = gnss_ekf_noise;
}

void h_model_GNSS_output(state_output &s, esekfom::dyn_share_modified<double> &ekfom_data)
{
	// Eigen::Vector3d res_r;
	// p_gnss->state_const_.rot.boxminus(res_r, s.rot);
	ekfom_data.h_GNSS.setIdentity();
	ekfom_data.h_GNSS *= p_gnss->odo_weight;
	ekfom_data.z_GNSS.setZero();
	// ekfom_data.h_GNSS.block<3, 3>(3, 3) = Eigen::Matrix3d::Zero(); // Jacob_right_inv<double>(res_r); // 
	// ekfom_data.h_GNSS.block<3, 3>(6, 6) = Eigen::Matrix3d::Zero(); // Jacob_right_inv<double>(res_r); // 
	// ekfom_data.z_GNSS.block<3, 1>(3, 0) = res_r;
	ekfom_data.z_GNSS.block<3, 1>(0, 0) = p_gnss->odo_weight * (p_gnss->state_const_.pos - s.pos); // 
	// ekfom_data.z_GNSS.block<3, 1>(6, 0) = p_gnss->state_const_.vel - s.vel;
	ekfom_data.M_Noise = gnss_ekf_noise;
}

void pointBodyToWorld(PointType const * const pi, PointType * const po)
{    
    V3D p_body(pi->x, pi->y, pi->z);
    
    V3D p_global;
	// if (extrinsic_est_en)
	// {	
	// 	if (!use_imu_as_input)
	// 	{
	// 		p_global = kf_output.x_.rot.normalized() * (kf_output.x_.offset_R_L_I * p_body + kf_output.x_.offset_T_L_I) + kf_output.x_.pos;
	// 	}
	// 	else
	// 	{
	// 		p_global = kf_input.x_.rot.normalized() * (kf_input.x_.offset_R_L_I * p_body + kf_input.x_.offset_T_L_I) + kf_input.x_.pos;
	// 	}
	// }
	// else
	{
		if (!use_imu_as_input)
		{
			p_global = kf_output.x_.rot * (Lidar_R_wrt_IMU * p_body + Lidar_T_wrt_IMU) + kf_output.x_.pos; // .normalized()
		}
		else
		{
			p_global = kf_input.x_.rot * (Lidar_R_wrt_IMU * p_body + Lidar_T_wrt_IMU) + kf_input.x_.pos; // .normalized()
		}
	}

    po->x = p_global(0);
    po->y = p_global(1);
    po->z = p_global(2);
    po->intensity = pi->intensity;
}

void LI_Init_update()
{
	int rematch_num = 0;
	bool nearest_search_en = true;
	state_input state_propagat = p_imu->state_LI_Init;
	/*** iterated state estimation ***/

	for (int iterCount = 0; iterCount < 5; iterCount++) {

		/** closest surface search and residual computation **/
		for (int i = 0; i < feats_down_size; i++) {
			PointType &point_body = feats_down_body->points[i];
			PointType &point_world = feats_down_world->points[i];
			V3D p_body(point_body.x, point_body.y, point_body.z);
			/// transform to world frame
			p_imu->pointBodyToWorld_li_init(&point_body, &point_world);
			std::vector<float> pointSearchSqDis(NUM_MATCH_POINTS);
			auto &points_near = Nearest_Points[i];
			uint8_t search_flag = 0;

			if (nearest_search_en) {
				/** Find the closest surfaces in the map **/
                ivox_->GetClosestPoint(point_world, points_near, NUM_MATCH_POINTS); // 
				if (points_near.size() < NUM_MATCH_POINTS)
					point_selected_surf[i] = false;
				else
					point_selected_surf[i] = !(pointSearchSqDis[NUM_MATCH_POINTS - 1] > 5);
			}

			if (!point_selected_surf[i] || points_near.size() < NUM_MATCH_POINTS) {
				point_selected_surf[i] = false;
				continue;
			}

			point_selected_surf[i] = false;
			VD(4) pabcd;
			pabcd.setZero();
			if (esti_plane(pabcd, points_near, 0.1)) //(planeValid)
			{
				float pd2 = pabcd(0) * point_world.x + pabcd(1) * point_world.y + pabcd(2) * point_world.z +
							pabcd(3);
				float s = 1 - 0.9 * fabs(pd2) / sqrt(p_body.norm());

				if (s > 0.9) {
					point_selected_surf[i] = true;
					normvec->points[i].x = pabcd(0);
					normvec->points[i].y = pabcd(1);
					normvec->points[i].z = pabcd(2);
					normvec->points[i].intensity = pd2;
				}
			}
		}
		int effect_feat_num = 0;
		for (int i = 0; i < feats_down_size; i++) {
			if (point_selected_surf[i]) {
				// laserCloudOri->points[effect_feat_num] = feats_down_body->points[i];
				// corr_normvect->points[effect_feat_num] = normvec->points[i];
				effect_feat_num++;
			}
		}

		/*** Computation of Measurement Jacobian matrix H and measurents vector ***/

		MatrixXd Hsub(effect_feat_num, 6);
		MatrixXd Hsub_T_R_inv(6, effect_feat_num);
		VectorXd meas_vec(effect_feat_num);

		Hsub.setZero();
		Hsub_T_R_inv.setZero();
		meas_vec.setZero();
		int m_ = 0;

		for (int i = 0; i < feats_down_size; i++) {
			if (point_selected_surf[i])
			{
				const PointType &laser_p = feats_down_body->points[i];
				V3D point_this_L(laser_p.x, laser_p.y, laser_p.z);

				// V3D point_this = Lidar_R_wrt_IMU * point_this_L + Lidar_T_wrt_IMU;
				M3D point_crossmat;
				point_crossmat << SKEW_SYM_MATRX(point_this_L);

				/*** get the normal vector of closest surface/corner ***/
				const PointType &norm_p = normvec->points[i];
				V3D norm_vec(norm_p.x, norm_p.y, norm_p.z);

				/*** calculate the Measurement Jacobian matrix H ***/
			
				V3D G(point_crossmat * p_imu->state_LI_Init.rot.transpose() * norm_vec); // conjugate().normalized()
				Hsub.row(m_) << VEC_FROM_ARRAY(G), norm_p.x, norm_p.y, norm_p.z; //, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;                    
				Hsub_T_R_inv.col(m_) = Hsub.row(m_).transpose() * 100;
				/*** Measurement: distance to the closest surface/corner ***/
				meas_vec(m_) = -norm_p.intensity;
				m_ ++;
			}
		}

		MatrixXd K(12, effect_feat_num);

		bool EKF_stop_flg = false;
		bool flg_EKF_converged = false;

		/*** Iterative Kalman Filter Update ***/
		Eigen::Matrix<double, 12, 1> solution;
		Eigen::Matrix<double, 12, 12> H_T_H;
		solution.setZero();
		H_T_H.setZero();
		H_T_H.block<6, 6>(0, 0) = Hsub_T_R_inv * Hsub;
		MD(12, 12) &&K_1 = (H_T_H + p_imu->state_cov.inverse()).inverse();
		K = K_1.block<12, 6>(0, 0) * Hsub_T_R_inv;
		Eigen::Matrix<double, 12, 1> vec;
		vec.block<3,1>(3,0) = Log<double>(p_imu->state_LI_Init.rot.transpose() * state_propagat.rot); // .normalized().toRotationMatrix()
		vec.block<3,1>(0,0) = state_propagat.pos - p_imu->state_LI_Init.pos;
		vec.block<3,1>(6,0) = state_propagat.vel - p_imu->state_LI_Init.vel;
		vec.block<3,1>(9,0) = state_propagat.bg - p_imu->state_LI_Init.bg;
		solution = K * meas_vec + vec - K * Hsub * vec.block<6, 1>(0, 0);

		//state update
		p_imu->state_LI_Init.rot.boxplus(solution.block<3, 1>(3, 0)); // += solution;
		p_imu->state_LI_Init.pos.boxplus(solution.block<3, 1>(0, 0)) ; // += solution;
		p_imu->state_LI_Init.vel.boxplus(solution.block<3, 1>(6, 0)) ; // += solution;
		p_imu->state_LI_Init.bg.boxplus(solution.block<3, 1>(9, 0)) ; // += solution;

		if ((solution.block<3, 1>(0, 0).norm() * 57.3 < 0.01) && (solution.block<3, 1>(3, 0).norm() * 100 < 0.015))
			flg_EKF_converged = true;

		/*** Rematch Judgement ***/
		nearest_search_en = false;
		if (flg_EKF_converged || ((rematch_num == 0) && (iterCount == (5 - 2)))) {
			nearest_search_en = true;
			rematch_num++;
		}

		/*** Convergence Judgements and Covariance Update ***/
		if (!EKF_stop_flg && (rematch_num >= 2 || (iterCount == 5 - 1))) {
		// if (flg_EKF_inited) {
			/*** Covariance Update ***/
			Eigen::Matrix<double, 12, 12> G;
			G.setZero();
			G.block<12, 6>(0, 0) = K * Hsub;
			p_imu->state_cov = (Eigen::Matrix<double, 12, 12>::Identity() - G) * p_imu->state_cov;
		// }
			EKF_stop_flg = true;
		}

		if (EKF_stop_flg) break;
	}
}