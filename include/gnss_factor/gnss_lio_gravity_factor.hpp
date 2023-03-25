#ifndef GNSS_LIO_GRAV_FACTOR_H_
#define GNSS_LIO_GRAV_FACTOR_H_

#include <Eigen/Dense>

#include <gtsam/nonlinear/Marginals.h>
#include <gtsam/base/Vector.h>
#include <gtsam/geometry/Pose3.h>
#include <gtsam/inference/Symbol.h>

/* 
**  initial anc, yaw, ddt, dt;
 */
namespace glio {

class GnssLioGravFactor : public gtsam::NoiseModelFactor1<gtsam::Rot3> //, gtsam::Vector3>
{
    public: 
        GnssLioGravFactor(gtsam::Key j1, Eigen::Matrix<double, 3, 1> grav_body_, Eigen::Vector3d grav_world_, const gtsam::SharedNoiseModel& model) : 
        grav_body(grav_body_), grav_world(grav_world_), gtsam::NoiseModelFactor1<gtsam::Rot3>(model, j1) {}
        
        virtual ~GnssLioGravFactor() {}
        gtsam::Vector evaluateError(const gtsam::Rot3 &rot,
        boost::optional<gtsam::Matrix&> H1 = boost::none) const
        {
            if (H1) 
            {
                (*H1) = gtsam::Matrix::Zero(3, 3);
                (*H1).block<3, 3>(0, 0) = -rot.matrix() * skew_sym_mat(grav_body); // gtsam::Matrix::Identity(4, 4);
            }
            
            gtsam::Vector residual(3);
            residual.segment<3>(0) = rot.matrix() * grav_body - grav_world;
            return residual;
        }
    private:
        Eigen::Matrix<double, 3, 1> grav_body, grav_world;
};
}
#endif