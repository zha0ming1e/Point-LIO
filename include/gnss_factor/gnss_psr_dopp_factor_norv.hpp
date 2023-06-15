#ifndef GNSS_PSR_DOPP_FACTOR_NORV_H_
#define GNSS_PSR_DOPP_FACTOR_NORV_H_

#include <vector>
#include <Eigen/Dense>
#include <gtsam/nonlinear/Marginals.h>
#include <gtsam/geometry/Pose3.h>
#include <gtsam/inference/Symbol.h>
#include <gtsam/base/Vector.h>
using namespace gnss_comm;

/* 
**  parameters[0]: position and orientation at time k
**  parameters[1]: velocity and acc/gyro bias at time k
**  parameters[2]: position and orientation at time k+1
**  parameters[3]: velocity and acc/gyro bias at time k+1
**  parameters[4]: receiver clock bias in light travelling distance (m)
**  parameters[5]: receiver clock bias change rate in clock bias light travelling distance per second (m/s)
**  parameters[6]: yaw difference between ENU and local coordinate (rad)
**  parameters[7]: anchor point's ECEF coordinate
**  
 */
namespace glio {

#define PSR_TO_DOPP_RATIO (5)

class GnssPsrDoppFactorNoRV : public gtsam::NoiseModelFactor4<gtsam::Vector3, gtsam::Vector4, gtsam::Vector3, gtsam::Rot3> 
{
    public: 
        // EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        // GnssPsrDoppFactor() = delete;
        GnssPsrDoppFactorNoRV(gtsam::Key j1, gtsam::Key j2, gtsam::Key j3, gtsam::Key j4, bool invalid_lidar_, double values_[30], int sys_idx_, 
        Eigen::Vector3d hat_omg_T_, const gtsam::SharedNoiseModel& model) :
        hat_omg_T(hat_omg_T_), sys_idx(sys_idx_), invalid_lidar(invalid_lidar_),
        gtsam::NoiseModelFactor4<gtsam::Vector3, gtsam::Vector4, gtsam::Vector3, gtsam::Rot3>(model, j1, j2, j3, j4) {
            Tex_imu_r << values_[0], values_[1], values_[2];
            anc_local << values_[3], values_[4], values_[5];
            sv_pos << values_[6], values_[7], values_[8];
            sv_vel << values_[9], values_[10], values_[11];
            svdt = values_[12];
            tgd = values_[13];
            svddt = values_[14];
            pr_uura = values_[15];
            dp_uura = values_[16];
            relative_sqrt_info = values_[17];
            for (int i = 0; i < 8; i++)
            {
                latest_gnss_iono_params.push_back(values_[18+i]);
            }
            time_current = values_[26];
            freq = values_[27];
            psr_measured = values_[28];
            dopp_measured = values_[29];
        }

        virtual ~GnssPsrDoppFactorNoRV() {}

        gtsam::Vector evaluateError(const gtsam::Vector3 &pos_vel_bias, const gtsam::Vector4 &dt, const gtsam::Vector3 &ext_p, const gtsam::Rot3 &ext_R, //const gtsam::Vector3 &anc, 
            boost::optional<gtsam::Matrix&> H1 = boost::none, boost::optional<gtsam::Matrix&> H2 = boost::none, 
            boost::optional<gtsam::Matrix&> H3 = boost::none, boost::optional<gtsam::Matrix&> H4 = boost::none) const //, boost::optional<gtsam::Matrix&> H7 = boost::none) const
        {
            Eigen::Vector3d ref_ecef = ext_p;

            const Eigen::Vector3d local_pos = Tex_imu_r + pos_vel_bias.segment<3>(0);

            Eigen::Matrix3d R_ecef_local = ext_R.matrix(); // R_ecef_enu_cur * R_enu_local;

            Eigen::Vector3d P_ecef = R_ecef_local * (local_pos - anc_local) + ref_ecef;

            double ion_delay = 0, tro_delay = 0;
            double azel[2] = {0, M_PI/2.0};
            if (P_ecef.norm() > 0)
            {
                sat_azel(P_ecef, sv_pos, azel);
                Eigen::Vector3d rcv_lla = ecef2geo(P_ecef);
                tro_delay = calculate_trop_delay(sec2time(time_current), rcv_lla, azel);
                ion_delay = calculate_ion_delay(sec2time(time_current), latest_gnss_iono_params, rcv_lla, azel); // rely on local pose
            }
            double sin_el = sin(azel[1]);
            double sin_el_2 = sin_el*sin_el;
            double pr_weight = sin_el_2 / pr_uura * relative_sqrt_info; // not requisite

            Eigen::Vector3d rcv2sat_ecef = sv_pos - P_ecef;
            Eigen::Vector3d rcv2sat_unit = rcv2sat_ecef.normalized();
            

            const double psr_sagnac = EARTH_OMG_GPS*(sv_pos(0)*P_ecef(1)-sv_pos(1)*P_ecef(0))/LIGHT_SPEED;
            double psr_estimated = rcv2sat_ecef.norm() + psr_sagnac + dt[sys_idx] - svdt*LIGHT_SPEED + // why not multiply light_speed?  
                                    ion_delay + tro_delay + tgd*LIGHT_SPEED;

            gtsam::Vector1 residual;
            
            {    
                residual[0] = (psr_estimated - psr_measured) * pr_weight;

                if (H1)
                {
                    (*H1) = gtsam::Matrix::Zero(1,3);
                    (*H1).block<1,3>(0,0) = -rcv2sat_unit.transpose() * R_ecef_local * pr_weight;
                }

                if (H2)
                {
                    (*H2) = gtsam::Matrix::Zero(1,4);
                    (*H2)(0, sys_idx) = 1.0 * pr_weight;
                }

                if (H3)
                {
                    (*H3) = gtsam::Matrix::Zero(1,3);
                    // if (!invalid_lidar)
                    {
                    (*H3).block<1,3>(0,0) = -rcv2sat_unit.transpose() * pr_weight;
                    // (*H5).block<1,3>(0,0) = -rcv2sat_unit.transpose() * (Eye3d + R_ecef_enu_cur * hatP * R1TE3 * vecLon.transpose() - R_ecef_enu_cur * hatP * E1 * vecLat.transpose()) * pr_weight;
                    }
                }

                if (H4)
                {
                    (*H4) = gtsam::Matrix::Zero(1,3);
                    // if (!invalid_lidar)
                    {
                    Eigen::Vector3d pos_v = local_pos - anc_local;
                    Eigen::Matrix3d d_pos;
                    if (hat_omg_T.norm() > 0.3)
                    {
                        d_pos << 0.0, -pos_v[2], pos_v[1], 
                                    pos_v[2], 0.0, -pos_v[0], 
                                    -pos_v[1], pos_v[0], 0.0;
                        // d_pos << 0.0, 0.0, pos_v[1], 
                        //             0.0, 0.0, -pos_v[0], 
                        //             0.0, 0.0, 0.0;
                        // d_vel << 0.0, 0.0, local_vel[1], 
                        //             0.0, 0.0, -local_vel[0], 
                        //             0.0, 0.0, 0.0;
                    // }
                    // else
                    // {
                    //     d_pos << 0.0, -pos_v[2], 0.0, 
                    //                 pos_v[2], 0.0, 0.0, 
                    //                 -pos_v[1], pos_v[0], 0.0;
                    //     d_vel << 0.0, -local_vel[2], 0.0, 
                    //                 local_vel[2], 0.0, 0.0, 
                    //                 -local_vel[1], local_vel[0], 0.0;
                    // }
                    (*H4).block<1,3>(0,0) = rcv2sat_unit.transpose() * (R_ecef_local * d_pos) * pr_weight;
                    }
                    }
                }
                return residual;
            }
        }
    private:
        Eigen::Vector3d Tex_imu_r, anc_local, sv_pos, sv_vel, hat_omg_T;
        double svdt, tgd, svddt, pr_uura, dp_uura, relative_sqrt_info, time_current, freq, psr_measured, dopp_measured;
        std::vector<double> latest_gnss_iono_params;
        int sys_idx;
        bool invalid_lidar;
};
}

#endif