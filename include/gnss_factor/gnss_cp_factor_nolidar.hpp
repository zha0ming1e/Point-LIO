#ifndef GNSS_CP_FACTOR_NOLIDAR_H_
#define GNSS_CP_FACTOR_NOLIDAR_H_

#include <vector>
#include <Eigen/Dense>
#include <gtsam/nonlinear/Marginals.h>
#include <gtsam/geometry/Pose3.h>
#include <gtsam/base/Vector.h>
#include <gtsam/inference/Symbol.h>
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

class GnssCpFactorNolidar : public gtsam::NoiseModelFactor4<gtsam::Rot3, gtsam::Vector12, gtsam::Rot3, gtsam::Vector12>
{
    public: 
        GnssCpFactorNolidar(gtsam::Key i1, gtsam::Key i2, gtsam::Key j1, gtsam::Key j2, double values_[20], const gtsam::SharedNoiseModel& model) :
        gtsam::NoiseModelFactor4<gtsam::Rot3, gtsam::Vector12, gtsam::Rot3, gtsam::Vector12>(model, i1, i2, j1, j2) {
            Tex_imu_r << values_[0], values_[1], values_[2];
            anc_local << values_[3], values_[4], values_[5];
            sv_pos_bi << values_[6], values_[7], values_[8];
            sv_pos_pi << values_[9], values_[10], values_[11];
            sv_pos_bj << values_[12], values_[13], values_[14];
            sv_pos_pj << values_[15], values_[16], values_[17];
            cp_measured = values_[18];
            cp_weight = values_[19];
        }

        virtual ~GnssCpFactorNolidar() {}
        
        gtsam::Vector evaluateError(const gtsam::Rot3 &rot1, const gtsam::Vector12 &pos1, const gtsam::Rot3 &rot2, const gtsam::Vector12 &pos2, 
            boost::optional<gtsam::Matrix&> H1 = boost::none, boost::optional<gtsam::Matrix&> H2 = boost::none, boost::optional<gtsam::Matrix&> H3 = boost::none, boost::optional<gtsam::Matrix&> H4 = boost::none) const
        {
            const Eigen::Vector3d local_pos1 = rot1 * Tex_imu_r + pos1.segment<3>(0);
            const Eigen::Vector3d local_pos2 = rot2 * Tex_imu_r + pos2.segment<3>(0);

            Eigen::Vector3d P_ecef1, P_ecef2;
            {
                P_ecef1 = local_pos1;
                P_ecef2 = local_pos2;
            }

            Eigen::Vector3d rcv2sat_ecef_bj = sv_pos_bj - P_ecef2;
            Eigen::Vector3d rcv2sat_ecef_bi = sv_pos_bi - P_ecef1;
            Eigen::Vector3d rcv2sat_ecef_pj = sv_pos_pj - P_ecef2;
            Eigen::Vector3d rcv2sat_ecef_pi = sv_pos_pi - P_ecef1;
            Eigen::Vector3d rcv2sat_unit_bj = rcv2sat_ecef_bj.normalized();
            Eigen::Vector3d rcv2sat_unit_bi = rcv2sat_ecef_bi.normalized();
            Eigen::Vector3d rcv2sat_unit_pj = rcv2sat_ecef_pj.normalized();
            Eigen::Vector3d rcv2sat_unit_pi = rcv2sat_ecef_pi.normalized();

            Eigen::Matrix3d hat_T;
            hat_T << SKEW_SYM_MATRX(Tex_imu_r);

            gtsam::Vector1 residual;

            {
                residual[0] = (rcv2sat_ecef_bj.norm() - rcv2sat_ecef_pj.norm() - rcv2sat_ecef_bi.norm() + rcv2sat_ecef_pi.norm() - cp_measured) * cp_weight;

                if (H1)
                {
                    (*H1) = gtsam::Matrix::Zero(1, 3);
                    (*H1).block<1,3>(0,0) = -(rcv2sat_unit_bi-rcv2sat_unit_pi).transpose() * hat_T * cp_weight;
                }

                if (H2)
                {
                    (*H2) = gtsam::Matrix::Zero(1, 12);
                    (*H2).block<1,3>(0,0) = (rcv2sat_unit_bi-rcv2sat_unit_pi).transpose() * cp_weight;
                }

                if (H3)
                {
                    (*H3) = gtsam::Matrix::Zero(1, 3);
                    (*H3).block<1,3>(0,0) = (rcv2sat_unit_bj-rcv2sat_unit_pj).transpose() * hat_T * cp_weight;
                }

                if (H4)
                {
                    (*H4) = gtsam::Matrix::Zero(1, 12);
                    (*H4).block<1,3>(0,0) = -(rcv2sat_unit_bj-rcv2sat_unit_pj).transpose() * cp_weight;
                }
                return residual;
            }
        }
    private:
        Eigen::Vector3d sv_pos_bi, sv_pos_bj, sv_pos_pi, sv_pos_pj, Tex_imu_r, anc_local;
        double cp_measured, cp_weight;
};
}

#endif