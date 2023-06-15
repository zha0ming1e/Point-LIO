#ifndef GNSS_CP_FACTOR_NORV_H_
#define GNSS_CP_FACTOR_NORV_H_

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

class GnssCpFactorNoRV : public gtsam::NoiseModelFactor4<gtsam::Vector3, gtsam::Rot3, gtsam::Vector3, gtsam::Vector3>
{
    public: 
        GnssCpFactorNoRV(gtsam::Key f1, gtsam::Key f2, gtsam::Key i1, gtsam::Key j1, bool invalid_lidar_, double values_[20], Eigen::Vector3d RTex2_, Eigen::Vector3d vel_, const gtsam::SharedNoiseModel& model) :
        gtsam::NoiseModelFactor4<gtsam::Vector3, gtsam::Rot3, gtsam::Vector3, gtsam::Vector3>(model, f1, f2, i1, j1) {
            Tex_imu_r << values_[0], values_[1], values_[2];
            anc_local << values_[3], values_[4], values_[5];
            sv_pos_bi << values_[6], values_[7], values_[8];
            sv_pos_pi << values_[9], values_[10], values_[11];
            sv_pos_bj << values_[12], values_[13], values_[14];
            sv_pos_pj << values_[15], values_[16], values_[17];
            cp_measured = values_[18];
            cp_weight = values_[19];
            invalid_lidar = invalid_lidar_;
            RTex2 = RTex2_;
            vel = vel_;
        }

        virtual ~GnssCpFactorNoRV() {}
        
        gtsam::Vector evaluateError(const gtsam::Vector3 &ext_p, const gtsam::Rot3 &ext_R, const gtsam::Vector3 &pos1, const gtsam::Vector3 &pos2, 
            boost::optional<gtsam::Matrix&> H1 = boost::none, boost::optional<gtsam::Matrix&> H2 = boost::none, 
            boost::optional<gtsam::Matrix&> H3 = boost::none, boost::optional<gtsam::Matrix&> H4 = boost::none) const
        {
            Eigen::Vector3d ref_ecef = ext_p;

            const Eigen::Vector3d local_pos1 = Tex_imu_r + pos1.segment<3>(0);
            const Eigen::Vector3d local_pos2 = RTex2 + pos2.segment<3>(0);

            // Eigen::Matrix3d R_enu_local = ext_R.matrix();
            // Eigen::Matrix3d R_ecef_enu_cur = ecef2rotation(ref_ecef); // provide anchor value
            Eigen::Matrix3d R_ecef_local = ext_R.matrix(); // R_ecef_enu_cur * R_enu_local;

            Eigen::Vector3d P_ecef1, P_ecef2;
            {
                P_ecef1 = R_ecef_local * (local_pos1 - anc_local) + ref_ecef;
                P_ecef2 = R_ecef_local * (local_pos2 - anc_local) + ref_ecef;
            }

            Eigen::Vector3d rcv2sat_ecef_bj = sv_pos_bj - P_ecef2;
            Eigen::Vector3d rcv2sat_ecef_bi = sv_pos_bi - P_ecef1;
            Eigen::Vector3d rcv2sat_ecef_pj = sv_pos_pj - P_ecef2;
            Eigen::Vector3d rcv2sat_ecef_pi = sv_pos_pi - P_ecef1;
            Eigen::Vector3d rcv2sat_unit_bj = rcv2sat_ecef_bj.normalized();
            Eigen::Vector3d rcv2sat_unit_bi = rcv2sat_ecef_bi.normalized();
            Eigen::Vector3d rcv2sat_unit_pj = rcv2sat_ecef_pj.normalized();
            Eigen::Vector3d rcv2sat_unit_pi = rcv2sat_ecef_pi.normalized();

            gtsam::Vector1 residual;
            {
                residual[0] = (rcv2sat_ecef_bj.norm() - rcv2sat_ecef_pj.norm() - rcv2sat_ecef_bi.norm() + rcv2sat_ecef_pi.norm() - cp_measured) * cp_weight;
                // cout << "check cp residual:" << residual << endl;
                if (H1)
                {
                    (*H1) = gtsam::Matrix::Zero(1,3);
                    
                    // if (!invalid_lidar)
                    {
                    (*H1).block<1,3>(0,0) = (-rcv2sat_unit_bj + rcv2sat_unit_pj).transpose() * cp_weight + (rcv2sat_unit_bi - rcv2sat_unit_pi).transpose() * cp_weight;
                    // (*H1).block<1,3>(0,0) = (-rcv2sat_unit_bj + rcv2sat_unit_pj).transpose() * (Eye3d + R_ecef_enu_cur * hatP2 * R1TE3 * vecLon.transpose() - R_ecef_enu_cur * hatP2 * E1 * vecLat.transpose()) * cp_weight
                    //          + (rcv2sat_unit_bi - rcv2sat_unit_pi).transpose() * (Eye3d + R_ecef_enu_cur * hatP1 * R1TE3 * vecLon.transpose() - R_ecef_enu_cur * hatP1 * E1 * vecLat.transpose()) * cp_weight;
                    }
                }

                if (H2)
                {
                    (*H2) = gtsam::Matrix::Zero(1,3);
                    // if (!invalid_lidar)
                    {
                    Eigen::Matrix3d d_pos1, d_pos2;
                    Eigen::Vector3d pos_v1 = local_pos1 - anc_local;
                    Eigen::Vector3d pos_v2 = local_pos2 - anc_local;
                    if (vel.norm() > 0.3)
                    {
                        d_pos1 << 0.0, -pos_v1[2], pos_v1[1], 
                                    pos_v1[2], 0.0, -pos_v1[0], 
                                    -pos_v1[1], pos_v1[0], 0.0;
                        d_pos2 << 0.0, -pos_v2[2], pos_v2[1], 
                                    pos_v2[2], 0.0, -pos_v2[0], 
                                    -pos_v2[1], pos_v2[0], 0.0;
                        // d_pos1 << 0.0, 0.0, pos_v1[1], 
                        //             0.0, 0.0, -pos_v1[0], 
                        //             0.0, 0.0, 0.0;
                        // d_pos2 << 0.0, 0.0, pos_v2[1], 
                        //             0.0, 0.0, -pos_v2[0], 
                        //             0.0, 0.0, 0.0;
                    // }
                    // else
                    // {
                    //     d_pos1 << 0.0, -pos_v1[2], 0.0, 
                    //                 pos_v1[2], 0.0, 0.0, 
                    //                 -pos_v1[1], pos_v1[0], 0.0;
                    //     d_pos2 << 0.0, -pos_v2[2], 0.0, 
                    //                 pos_v2[2], 0.0, 0.0, 
                    //                 -pos_v2[1], pos_v2[0], 0.0;
                    // }
                    (*H2).block<1,3>(0,0) = (rcv2sat_unit_bj - rcv2sat_unit_pj).transpose() * R_ecef_local * d_pos2 * cp_weight
                                -(rcv2sat_unit_bi - rcv2sat_unit_pi).transpose() * R_ecef_local * d_pos1 * cp_weight;
                    // printf("check hessian:%f, %f, %f\n", (*H2)(0, 0), (*H2)(0, 1), (*H2)(0, 2));
                    }
                    }
                }

                if (H3)
                {
                    (*H3) = gtsam::Matrix::Zero(1, 3);
                    (*H3).block<1,3>(0,0) = (rcv2sat_unit_bi-rcv2sat_unit_pi).transpose() * R_ecef_local * cp_weight;
                }

                if (H4)
                {
                    (*H4) = gtsam::Matrix::Zero(1, 3);
                    (*H4).block<1,3>(0,0) = -(rcv2sat_unit_bj-rcv2sat_unit_pj).transpose() * R_ecef_local * cp_weight;
                }
                return residual;
            }
        }
    private:
        Eigen::Vector3d sv_pos_bi, sv_pos_bj, sv_pos_pi, sv_pos_pj, Tex_imu_r, anc_local, RTex2, vel;
        double cp_measured, cp_weight;
        bool invalid_lidar;
};
}

#endif