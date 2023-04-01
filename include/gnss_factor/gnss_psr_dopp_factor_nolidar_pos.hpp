#ifndef GNSS_PSR_DOPP_FACTOR_NOLIDAR_POS_H_
#define GNSS_PSR_DOPP_FACTOR_NOLIDAR_POS_H_

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

class GnssPsrDoppFactorNolidarPos : public gtsam::NoiseModelFactor3<gtsam::Vector12, gtsam::Vector4, gtsam::Vector1>
{
    public: 
        // EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        // GnssPsrDoppFactor() = delete;
        GnssPsrDoppFactorNolidarPos(gtsam::Key j1, gtsam::Key j2, gtsam::Key j3, double values_[30], int sys_idx_, Eigen::Vector3d hat_omg_T_, Eigen::Matrix3d rot_, const gtsam::SharedNoiseModel& model) :
        gtsam::NoiseModelFactor3<gtsam::Vector12, gtsam::Vector4, gtsam::Vector1>(model, j1, j2, j3), sys_idx(sys_idx_), hat_omg_T(hat_omg_T_), rot(rot_) {
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
        virtual ~GnssPsrDoppFactorNolidarPos() {}
        gtsam::Vector evaluateError(const gtsam::Vector12 &pos_vel_bias, const gtsam::Vector4 &dt, const gtsam::Vector1 &ddt,
            boost::optional<gtsam::Matrix&> H1 = boost::none, boost::optional<gtsam::Matrix&> H2 = boost::none, boost::optional<gtsam::Matrix&> H3 = boost::none) const
        {          
            const Eigen::Vector3d local_pos = rot * Tex_imu_r + pos_vel_bias.segment<3>(0);
            const Eigen::Vector3d local_vel = pos_vel_bias.segment<3>(3) + hat_omg_T;

            Eigen::Vector3d P_ecef = local_pos;
            
            Eigen::Vector3d V_ecef = local_vel;

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
            double dp_weight = sin_el_2 / dp_uura * relative_sqrt_info * PSR_TO_DOPP_RATIO; // not requisite

            Eigen::Vector3d rcv2sat_ecef = sv_pos - P_ecef;
            Eigen::Vector3d rcv2sat_unit = rcv2sat_ecef.normalized();
            
            const double wavelength = LIGHT_SPEED / freq;

            const double psr_sagnac = EARTH_OMG_GPS*(sv_pos(0)*P_ecef(1)-sv_pos(1)*P_ecef(0))/LIGHT_SPEED;
            double psr_estimated = rcv2sat_ecef.norm() + psr_sagnac + dt[sys_idx] - svdt*LIGHT_SPEED + // why not multiply light_speed?  
                                    ion_delay + tro_delay + tgd*LIGHT_SPEED;
            const double dopp_sagnac = EARTH_OMG_GPS/LIGHT_SPEED*(sv_vel(0)*P_ecef(1)+
                    sv_pos(0)*V_ecef(1) - sv_vel(1)*P_ecef(0) - sv_pos(1)*V_ecef(0));
            double dopp_estimated = (sv_vel - V_ecef).dot(rcv2sat_unit) + ddt[0] + dopp_sagnac - svddt*LIGHT_SPEED;

            gtsam::Vector2 residual;
            
            {    
                residual[0] = (psr_estimated - psr_measured) * pr_weight;
            
                residual[1] = (dopp_estimated + dopp_measured*wavelength) * dp_weight;
                // cout << "check residual:" << residual.transpose() << endl;
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

                if (H1)
                {
                    (*H1) = gtsam::Matrix::Zero(2,12);
                    (*H1).block<1,3>(0,0) = -rcv2sat_unit.transpose() * pr_weight; 
                    (*H1).block<1,3>(1,0) = (sv_vel-V_ecef).transpose() * unit2rcv_pos * dp_weight;
                    (*H1).block<1,3>(1,3) = rcv2sat_unit.transpose() * (-1.0) * dp_weight;
                }

                if (H2)
                {
                    (*H2) = gtsam::Matrix::Zero(2,4);
                    (*H2)(0, sys_idx) = 1.0 * pr_weight;
                }

                if (H3)
                {
                    (*H3) = gtsam::Matrix::Zero(2,1);
                    (*H3)(1,0) = 1.0 * dp_weight;
                }
                return residual;
            }
        }
    private:
        Eigen::Vector3d Tex_imu_r, anc_local, sv_pos, sv_vel, hat_omg_T;
        Eigen::Matrix3d rot;
        double svdt, tgd, svddt, pr_uura, dp_uura, relative_sqrt_info, time_current, freq, psr_measured, dopp_measured;
        std::vector<double> latest_gnss_iono_params;
        int sys_idx;        
};
}

#endif