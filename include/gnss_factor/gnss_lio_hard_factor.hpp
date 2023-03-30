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

class GnssLioHardFactor : public gtsam::NoiseModelFactor2<gtsam::Rot3, gtsam::Vector12> //, gtsam::Vector3>
{
    public: 
        GnssLioHardFactor(gtsam::Key j1, gtsam::Key j2, Eigen::Matrix<double, 12, 1> pos_vel_bias_lio_, Eigen::Matrix3d rot_lio_, const gtsam::SharedNoiseModel& model) : 
        rot_lio(rot_lio_), pos_vel_bias_lio(pos_vel_bias_lio_), gtsam::NoiseModelFactor2<gtsam::Rot3, gtsam::Vector12>(model, j1, j2) {}
        
        virtual ~GnssLioHardFactor() {}
        gtsam::Vector evaluateError(const gtsam::Rot3 &rot, const gtsam::Vector12 &pos_vel_bias,
        boost::optional<gtsam::Matrix&> H1 = boost::none, boost::optional<gtsam::Matrix&> H2 = boost::none) const
        {
            Eigen::Matrix3d res_R = rot_lio.transpose() * rot.matrix();
            Eigen::Vector3d res_r = gtsam::Rot3::Logmap(gtsam::Rot3(res_R));
            if (H1) 
            {
                (*H1) = gtsam::Matrix::Zero(15, 3);
                (*H1).block<3, 3>(0, 0) = Jacob_right_inv<double>(res_r); // gtsam::Matrix::Identity(4, 4);
            }
            if (H2)
            {
                (*H2) = gtsam::Matrix::Zero(15, 12);
                (*H2).block<12, 12>(3, 0) = gtsam::Matrix::Identity(12, 12);
                // (*H2).block<6, 6>(9, 6) = gtsam::Matrix::Zero(6, 6);
            }
            gtsam::Vector residual(15);
            residual.segment<3>(0) = res_r;
            residual.segment<12>(3) = pos_vel_bias - pos_vel_bias_lio;
            return residual;
        }
    private:
        Eigen::Matrix<double, 12, 1> pos_vel_bias_lio;
        Eigen::Matrix3d rot_lio;
};
}
#endif