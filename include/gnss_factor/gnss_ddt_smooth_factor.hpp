#ifndef DDT_SMOOTH_FACTOR_H_
#define DDT_SMOOTH_FACTOR_H_

#include <Eigen/Dense>

// #include <gtsam/geometry/Rot3.h>
// #include <gtsam/geometry/Pose3.h>
// #include <gtsam/slam/PriorFactor.h>
// #include <gtsam/slam/BetweenFactor.h>
#include <gtsam/nonlinear/Marginals.h>
// #include <gtsam/nonlinear/Values.h>
#include <gtsam/inference/Symbol.h>
#include <gtsam/base/Vector.h>

/* 
**  parameters[0]: rev_ddt (t)     in light travelling distance (m)
**  parameters[1]: rev_ddt (t+1)   in light travelling distance (m)
 */
namespace glio {

class DdtSmoothFactor : public gtsam::NoiseModelFactor2<gtsam::Vector1, gtsam::Vector1>
{
    public: 
        DdtSmoothFactor(gtsam::Key i, gtsam::Key j, const gtsam::SharedNoiseModel& model) : 
        gtsam::NoiseModelFactor2<gtsam::Vector1, gtsam::Vector1>(model, i, j) {
        }
        virtual ~DdtSmoothFactor() {}
        gtsam::Vector evaluateError(const gtsam::Vector1 &ddti, const gtsam::Vector1 &ddtj, boost::optional<gtsam::Matrix&> H1 = boost::none, boost::optional<gtsam::Matrix&> H2 = boost::none) const
        {
            if (H1) //(*H1) = (gtsam::Matrix(1,1)<<-1.0).finished();
            {
                (*H1) = - gtsam::Matrix::Identity(1, 1);
            }
            if (H2) //(*H2) = (gtsam::Matrix(1,1)<<1.0).finished();
            {
                (*H2) = gtsam::Matrix::Identity(1, 1);
            }
            gtsam::Vector1 residual;
            residual(0) = ddtj[0] - ddti[0];
            return residual;
            // return (gtsam::Vector(1) << ddtj - ddti).finished();
        }
    // private:
        // double weight_;
};
} 

#endif