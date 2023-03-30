#ifndef PRIOR_FACTOR_H_
#define PRIOR_FACTOR_H_

#include <Eigen/Dense>

#include <gtsam/nonlinear/Marginals.h>
#include <gtsam/base/Vector.h>
#include <gtsam/inference/Symbol.h>
using namespace gnss_comm;

/* 
**  initial anc, yaw, ddt, dt;
 */
namespace glio {

class PriorFactor : public gtsam::NoiseModelFactor4<gtsam::Vector4, gtsam::Vector1, gtsam::Vector3, gtsam::Rot3> //, gtsam::Vector3>
{
    public: 
        PriorFactor(gtsam::Key j1, gtsam::Key j2, gtsam::Key f1, gtsam::Key f2, double init_value_[8], Eigen::Matrix3d rot_, const gtsam::SharedNoiseModel& model) : 
        rot(rot_), gtsam::NoiseModelFactor4<gtsam::Vector4, gtsam::Vector1, gtsam::Vector3, gtsam::Rot3>(model, j1, j2, f1, f2) {
            memcpy(init_value, init_value_, 8 * sizeof(double));
        }
        
        virtual ~PriorFactor() {}
        gtsam::Vector evaluateError(const gtsam::Vector4 &dtj, const gtsam::Vector1 &ddtj, const gtsam::Vector3 &ext_p, const gtsam::Rot3 &ext_R, 
        boost::optional<gtsam::Matrix&> H1 = boost::none, boost::optional<gtsam::Matrix&> H2 = boost::none, boost::optional<gtsam::Matrix&> H3 = boost::none, boost::optional<gtsam::Matrix&> H4 = boost::none) const
        {
            if (H1) 
            {
                (*H1) = gtsam::Matrix::Zero(11, 4);
                (*H1).block<4, 4>(0, 0) = gtsam::Matrix::Identity(4, 4);
            }
            if (H2)
            {
                (*H2) = gtsam::Matrix::Zero(11, 1);
                (*H2)(4, 0) = 1.0;
            }
            if (H3) 
            {
                (*H3) = gtsam::Matrix::Zero(11,3);
                (*H3).block<3,3>(5,0) = gtsam::Matrix::Identity(3,3);
            }
            if (H4) 
            {
                (*H4) = gtsam::Matrix::Zero(11,3);
                (*H4).block<3,3>(8,0) = gtsam::Matrix::Identity(3,3);
            }
            gtsam::Vector11 residual;
            residual << dtj(0) - init_value[0], dtj(1) - init_value[1], dtj(2) - init_value[2], dtj(3) - init_value[3], ddtj(0) - init_value[4], ext_p(0) - init_value[5], ext_p(1) - init_value[6], 
                                        ext_p(2) - init_value[7], 0.0, 0.0, 0.0;
            Eigen::Matrix3d res_R = rot.transpose() * ext_R.matrix();
            residual.segment<3>(8) = gtsam::Rot3::Logmap(gtsam::Rot3(res_R));
            return residual;
        }
    private:
        double init_value[8];
        Eigen::Matrix3d rot;
};
}

#endif