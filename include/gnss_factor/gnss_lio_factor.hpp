#ifndef GNSS_LIO_FACTOR_H_
#define GNSS_LIO_FACTOR_H_

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

class GnssLioFactor : public gtsam::NoiseModelFactor4<gtsam::Rot3, gtsam::Vector6, gtsam::Rot3, gtsam::Vector6>
{
    public: 
        GnssLioFactor(gtsam::Key i1, gtsam::Key i2, gtsam::Key j1, gtsam::Key j2, double odo_weight_, gtsam::Rot3 rel_rot_, gtsam::Point3 rel_pos_, gtsam::Vector3 rel_vel_, 
        const gtsam::Vector3 &grav_, double dt_, const gtsam::SharedNoiseModel& model) :
        gtsam::NoiseModelFactor4<gtsam::Rot3, gtsam::Vector6, gtsam::Rot3, gtsam::Vector6>(model, i1, i2, j1, j2), rel_rot(rel_rot_), rel_pos(rel_pos_), 
        rel_vel(rel_vel_), grav(grav_), dt(dt_), odo_weight(odo_weight_) {
        }
        virtual ~GnssLioFactor() {}
        gtsam::Vector evaluateError(const gtsam::Rot3 &rot1, const gtsam::Vector6 &pos_vel_bias1, const gtsam::Rot3 &rot2, const gtsam::Vector6 &pos_vel_bias2, 
            boost::optional<gtsam::Matrix&> H1 = boost::none, boost::optional<gtsam::Matrix&> H2 = boost::none, boost::optional<gtsam::Matrix&> H3 = boost::none, boost::optional<gtsam::Matrix&> H4 = boost::none) const
        {

            Eigen::Matrix3d d = rot1.transpose() * rot2.matrix();
            Eigen::Vector3d delta_p = rot1.transpose() * (pos_vel_bias2.segment<3>(0) - pos_vel_bias1.segment<3>(0) - pos_vel_bias1.segment<3>(3) * dt - 0.5 * grav * dt * dt);
            Eigen::Vector3d delta_v = rot1.transpose() * (pos_vel_bias2.segment<3>(3) - pos_vel_bias1.segment<3>(3) - grav * dt);
            double sqrt_info = 0.1 / dt * odo_weight; // could change the rmse

            Eigen::Matrix3d res_R = rel_rot.transpose() * d;//.matrix();
            Eigen::Vector3d res_r = gtsam::Rot3::Logmap(gtsam::Rot3(res_R));

            if (H1) 
            {
                (*H1) = gtsam::Matrix::Zero(9,3);
                (*H1).block<3,3>(0,0) = -Jacob_right_inv<double>(res_r) * d.transpose();  //- d.transpose(); 
                (*H1).block<3,3>(3,0) = skew_sym_mat(delta_p); 
                (*H1).block<3,3>(6,0) = skew_sym_mat(delta_v); 
                (*H1) = sqrt_info * (*H1);
            }
            if (H2) 
            {
                (*H2) = gtsam::Matrix::Zero(9,6);
                (*H2).block<3,3>(3,0) = - sqrt_info * rot1.transpose(); // - gtsam::Matrix::Identity(3,3); 
                (*H2).block<3,3>(3,3) = - sqrt_info * rot1.transpose() * dt; // gtsam::Matrix::Identity(3,3);
                (*H2).block<3,3>(6,3) = - sqrt_info * rot1.transpose(); // - gtsam::Matrix::Identity(3,3);
                // (*H2).block<3,3>(9,6) = - 0.01 * gtsam::Matrix::Identity(3,3); // dt
                // (*H2).block<3,3>(12,9) = - 0.01 * gtsam::Matrix::Identity(3,3); // dt
                // (*H3) = sqrt_info * (*H3);
            }
            if (H3) 
            {
                (*H3) = gtsam::Matrix::Zero(9,3);
                (*H3).block<3,3>(0,0) = Jacob_right_inv<double>(res_r); // gtsam::Matrix::Identity(3,3); 
                (*H3) = sqrt_info * (*H3);
            }
            if (H4) 
            {
                (*H4) = gtsam::Matrix::Zero(9,6);
                (*H4).block<3,3>(3,0) = sqrt_info * rot1.transpose(); // gtsam::Matrix::Identity(3,3); 
                (*H4).block<3,3>(6,3) = sqrt_info * rot1.transpose(); // gtsam::Matrix::Identity(3,3);
                // (*H4).block<3,3>(9,6) = 0.1 * gtsam::Matrix::Identity(3,3); // dt
                // (*H4).block<3,3>(12,9) = 0.1 * gtsam::Matrix::Identity(3,3); // dt
                // (*H4).block<3,3>(9,6) = 0.1 * gtsam::Matrix::Identity(3,3); // dt
                // (*H4).block<3,3>(12,9) = 0.1 * gtsam::Matrix::Identity(3,3); // dt
                // (*H6) = sqrt_info * (*H6);
            }
            gtsam::Vector residual(9);
            // Eigen::Matrix3d res_R = rel_rot.transpose() * d.matrix();
            residual.segment<3>(0) = res_r; // gtsam::Rot3::Logmap(gtsam::Rot3(res_R));
            residual.segment<3>(3) = delta_p - rel_pos; // pos2- pos1 - rel_pos;
            residual.segment<3>(6) = delta_v - rel_vel; // vel2 - vel1 - rel_vel;
            // residual.segment<6>(9) = vel2.segment<6>(3) - vel1.segment<6>(3);
            // residual = sqrt_info * residual;
            // residual.segment<6>(9) = 0.1 * (pos_vel_bias2.segment<6>(6) - pos_vel_bias1.segment<6>(6)); // dt gtsam::Vector6::Zero(); //
            residual.segment<9>(0) = sqrt_info * residual.segment<9>(0);
            return residual;
        }
    private:
        gtsam::Rot3 rel_rot;
        gtsam::Point3 rel_pos;
        gtsam::Vector3 rel_vel;
        gtsam::Vector3 grav;
        double dt, odo_weight;
};
}

#endif