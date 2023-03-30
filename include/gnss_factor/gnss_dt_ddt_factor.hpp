
#ifndef DT_DDT_FACTOR_H_
#define DT_DDT_FACTOR_H_

#include <map>
#include <Eigen/Dense>
#include <gtsam/inference/Symbol.h>
#include <gtsam/base/Vector.h>
#include <gtsam/nonlinear/Marginals.h>
using namespace gnss_comm;

/* 
**  parameters[0]: rev_dt (t)   in light travelling distance (m)
**  parameters[1]: rev_dt (t+1) in light travelling distance (m)
**  parameters[2]: rev_ddt(t)   in light travelling distance per second (m/s)
**  parameters[3]: rev_ddt(t+1) in light travelling distance per second (m/s)
 */
namespace glio {

class DtDdtFactor : public gtsam::NoiseModelFactor4<gtsam::Vector4, gtsam::Vector4, gtsam::Vector1, gtsam::Vector1>
{
    public: 
        DtDdtFactor(gtsam::Key i1, gtsam::Key j1, gtsam::Key i2, gtsam::Key j2, bool rcv_sys_[4], const double delta_t_, const gtsam::SharedNoiseModel& model) : 
        gtsam::NoiseModelFactor4<gtsam::Vector4, gtsam::Vector4, gtsam::Vector1, gtsam::Vector1>(model, i1, j1, i2, j2), delta_t(delta_t_) 
        {
            // rcv_sys[0] = rcv_sys_[0]; rcv_sys[1] = rcv_sys_[1]; rcv_sys[2] = rcv_sys_[2]; rcv_sys[3] = rcv_sys_[3];
            rcv_sys[0] = true; rcv_sys[1] = true; rcv_sys[2] = true; rcv_sys[3] = true;
        }
        virtual ~DtDdtFactor() {}
        gtsam::Vector evaluateError(const gtsam::Vector4 &dti, const gtsam::Vector4 &dtj, const gtsam::Vector1 &ddti, const gtsam::Vector1 &ddtj, boost::optional<gtsam::Matrix&> H1 = boost::none, boost::optional<gtsam::Matrix&> H2 = boost::none, boost::optional<gtsam::Matrix&> H3 = boost::none, boost::optional<gtsam::Matrix&> H4 = boost::none) const
        {
            if (H1) 
            {
                (*H1) = gtsam::Matrix::Zero(4, 4); 
                if (rcv_sys[0]) (*H1)(0, 0) = -1.0;
                if (rcv_sys[1]) (*H1)(1, 1) = -1.0;
                if (rcv_sys[2]) (*H1)(2, 2) = -1.0;
                if (rcv_sys[3]) (*H1)(3, 3) = -1.0;
            }
            if (H2) 
            {
                (*H2) = gtsam::Matrix::Zero(4, 4); 
                if (rcv_sys[0]) (*H2)(0, 0) = 1.0;
                if (rcv_sys[1]) (*H2)(1, 1) = 1.0;
                if (rcv_sys[2]) (*H2)(2, 2) = 1.0;
                if (rcv_sys[3]) (*H2)(3, 3) = 1.0;
            }
            if (H3) 
            {
                (*H3) = gtsam::Matrix::Zero(4, 1);
                if (rcv_sys[0]) (*H3)(0, 0) = - 0.5 * delta_t;
                if (rcv_sys[1]) (*H3)(1, 0) = - 0.5 * delta_t;
                if (rcv_sys[2]) (*H3)(2, 0) = - 0.5 * delta_t;
                if (rcv_sys[3]) (*H3)(3, 0) = - 0.5 * delta_t;
            }
            if (H4) 
            {
                (*H4) = gtsam::Matrix::Zero(4, 1);
                if (rcv_sys[0]) (*H4)(0, 0) = - 0.5 * delta_t;
                if (rcv_sys[1]) (*H4)(1, 0) = - 0.5 * delta_t;
                if (rcv_sys[2]) (*H4)(2, 0) = - 0.5 * delta_t;
                if (rcv_sys[3]) (*H4)(3, 0) = - 0.5 * delta_t;
            }

            gtsam::Vector4 ddt_4;
            ddt_4 = gtsam::Matrix::Zero(4, 1);
            if (rcv_sys[0]) ddt_4(0, 0) = dtj[0] - dti[0] - 0.5 * (ddti[0] + ddtj[0]) * delta_t;
            if (rcv_sys[1]) ddt_4(1, 0) = dtj[1] - dti[1] - 0.5 * (ddti[0] + ddtj[0]) * delta_t;
            if (rcv_sys[2]) ddt_4(2, 0) = dtj[2] - dti[2] - 0.5 * (ddti[0] + ddtj[0]) * delta_t;
            if (rcv_sys[3]) ddt_4(3, 0) = dtj[3] - dti[3] - 0.5 * (ddti[0] + ddtj[0]) * delta_t;
            // ddt_4 = dtj - dti - ddt_4;
            return ddt_4;
        }
    
    private:
        double delta_t;
        bool rcv_sys[4];
};
}

#endif