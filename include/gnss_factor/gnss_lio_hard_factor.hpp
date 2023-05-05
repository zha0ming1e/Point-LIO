#ifndef GNSS_LIO_HARD_FACTOR_H_
#define GNSS_LIO_HARD_FACTOR_H_

#include <Eigen/Dense>

#include <gtsam/nonlinear/Marginals.h>
#include <gtsam/base/Vector.h>
#include <gtsam/geometry/Pose3.h>
#include <gtsam/inference/Symbol.h>
using namespace gnss_comm;

/* 
**  initial anc, yaw, ddt, dt;
 */
namespace glio {

class GnssLioHardFactor : public gtsam::NoiseModelFactor2<gtsam::Rot3, gtsam::Vector6> //, gtsam::Vector3>
{
    public: 
        GnssLioHardFactor(gtsam::Key j1, gtsam::Key j2, Eigen::Vector3d pos_lio_, Eigen::Vector3d vel_lio_, Eigen::Matrix3d rot_lio_, Eigen::Matrix<double, 9, 9>  sqrt_lidar_, bool no_weight_, const gtsam::SharedNoiseModel& model) : 
        rot_lio(rot_lio_), pos_lio(pos_lio_), vel_lio(vel_lio_), sqrt_lidar(sqrt_lidar_), no_weight(no_weight_), gtsam::NoiseModelFactor2<gtsam::Rot3, gtsam::Vector6>(model, j1, j2) {}
        
        virtual ~GnssLioHardFactor() {}
        gtsam::Vector evaluateError(const gtsam::Rot3 &rot, const gtsam::Vector6 &pos_vel,
        boost::optional<gtsam::Matrix&> H1 = boost::none, boost::optional<gtsam::Matrix&> H2 = boost::none) const
        {
            Eigen::Matrix3d res_R = rot_lio.transpose() * rot.matrix();
            Eigen::Vector3d res_r = gtsam::Rot3::Logmap(gtsam::Rot3(res_R));
            if (H1) 
            {
                (*H1) = gtsam::Matrix::Zero(9, 3);
                (*H1).block<3, 3>(3, 0) = Jacob_right_inv<double>(res_r); // gtsam::Matrix::Identity(4, 4);
                if (!no_weight) (*H1) = sqrt_lidar * (*H1);
            }
            if (H2)
            {
                (*H2) = gtsam::Matrix::Zero(9, 6);
                (*H2).block<3, 3>(0, 0) = gtsam::Matrix::Identity(3, 3);
                (*H2).block<3, 3>(6, 3) = gtsam::Matrix::Identity(3, 3);
                // (*H2).block<6, 6>(9, 6) = gtsam::Matrix::Zero(6, 6);
                if (!no_weight) (*H2) = sqrt_lidar * (*H2);
            }
            gtsam::Vector residual(9);
            residual.segment<3>(3) = res_r;
            residual.segment<3>(0) = pos_vel.block<3, 1>(0, 0) - pos_lio;
            residual.segment<3>(6) = pos_vel.block<3, 1>(3, 0) - vel_lio;
            if (!no_weight) residual = sqrt_lidar * residual;
            return residual;
        }
    private:
        Eigen::Matrix<double, 3, 1> pos_lio;
        Eigen::Matrix<double, 3, 1> vel_lio;
        Eigen::Matrix3d rot_lio;
        Eigen::Matrix<double, 9, 9> sqrt_lidar;
        bool no_weight;
};
}
#endif