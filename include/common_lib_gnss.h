#ifndef COMMON_LIB_GNSS_H
#define COMMON_LIB_GNSS_H

#include <common_lib.h>

using namespace std;
using namespace Eigen;

#define DIM_STATE_GNSS (24)

struct StatesGroupwithGNSS2
{
    StatesGroupwithGNSS2() {
		this->rot_end = M3D::Identity();
		this->pos_end = Zero3d;
        this->vel_end = Zero3d;
        this->omg     = Zero3d;
        this->acc     = Zero3d;
        this->gravity = Zero3d;
        this->bias_g  = Zero3d;
        this->bias_a  = Zero3d;
        this->yaw_enu_local = 0.0;
        this->dt_g = 0.0;
        this->dt_r = 0.0;
        this->dt_e = 0.0;
        this->dt_c = 0.0;
        this->ddt = 0.0;
        this->anc = Zero3d;
        this->cov     = MD(DIM_STATE_GNSS+9,DIM_STATE_GNSS+9)::Identity() * 0.01;
        // this->cov.block<9,9>(9,9) = MD(9,9)::Identity() * 0.01;
        this->Q       = MD(26, 26)::Identity();
	};

    StatesGroupwithGNSS2(const StatesGroupwithGNSS2& b) {
		this->rot_end = b.rot_end;
		this->pos_end = b.pos_end;
        this->vel_end = b.vel_end;
        this->omg     = b.omg;
        this->acc     = b.acc;
        this->gravity = b.gravity;

        this->yaw_enu_local = b.yaw_enu_local;
        this->dt_g = b.dt_g;
        this->dt_r = b.dt_r;
        this->dt_e = b.dt_e;
        this->dt_c = b.dt_c;
        this->ddt = b.ddt;

        this->anc = b.anc;

        this->bias_g  = b.bias_g;
        this->bias_a  = b.bias_a;
        this->cov     = b.cov;
        this->Q       = b.Q;
	};

    StatesGroupwithGNSS2& operator=(const StatesGroupwithGNSS2& b)
	{
        this->rot_end = b.rot_end;
		this->pos_end = b.pos_end;
        this->vel_end = b.vel_end;
        this->omg     = b.omg;
        this->acc     = b.acc;
        this->gravity = b.gravity;

        this->yaw_enu_local = b.yaw_enu_local;
        this->dt_g = b.dt_g;
        this->dt_r = b.dt_r;
        this->dt_e = b.dt_e;
        this->dt_c = b.dt_c;
        this->ddt = b.ddt;

        this->anc = b.anc;

        this->bias_g  = b.bias_g;
        this->bias_a  = b.bias_a;
        this->cov     = b.cov;
        this->Q       = b.Q;
        return *this;
	};

    
    StatesGroupwithGNSS2 operator+(const Matrix<double, DIM_STATE_GNSS+9, 1> &state_add)
	{
        StatesGroupwithGNSS2 a;
		a.rot_end = this->rot_end * Exp(state_add(0,0), state_add(1,0), state_add(2,0));
		a.pos_end = this->pos_end + state_add.block<3,1>(3,0);
        a.vel_end = this->vel_end + state_add.block<3,1>(6,0);
        a.omg     = this->omg     + state_add.block<3,1>(9,0);
        a.acc     = this->acc     + state_add.block<3,1>(12,0);
        a.gravity = this->gravity + state_add.block<3,1>(15,0);
        a.bias_g  = this->bias_g  + state_add.block<3,1>(18,0);
        a.bias_a  = this->bias_a  + state_add.block<3,1>(21,0);

        a.yaw_enu_local = this->yaw_enu_local + state_add(24);
        a.dt_g = this->dt_g + state_add(25);
        a.dt_r = this->dt_r + state_add(26);
        a.dt_e = this->dt_e + state_add(27);
        a.dt_c = this->dt_c + state_add(28);
        a.ddt = this->ddt + state_add(29);

        a.anc = this->anc + state_add.block<3,1>(30,0);
        // a.cov     = this->cov;
		return a;
	};

    StatesGroupwithGNSS2& operator+=(const Matrix<double, DIM_STATE_GNSS+9, 1> &state_add)
	{
        this->rot_end  = this->rot_end * Exp(state_add(0,0), state_add(1,0), state_add(2,0));
		this->pos_end += state_add.block<3,1>(3,0);
        this->vel_end += state_add.block<3,1>(6,0);
        this->omg     += state_add.block<3,1>(9,0);
        this->acc     += state_add.block<3,1>(12,0);
        this->gravity += state_add.block<3,1>(15,0);
        this->bias_g  += state_add.block<3,1>(18,0);
        this->bias_a  += state_add.block<3,1>(21,0);

        this->yaw_enu_local += state_add(24);
        this->dt_g += state_add(25);
        this->dt_r += state_add(26);
        this->dt_e += state_add(27);
        this->dt_c += state_add(28);
        this->ddt += state_add(29);

        this->anc += state_add.block<3,1>(30,0);
		return *this;
	};

    void update(const Matrix<double, DIM_STATE_GNSS + 9, 1> &state_add)
	{
        this->rot_end  = this->rot_end * Exp(state_add(0,0), state_add(1,0), state_add(2,0));
		this->pos_end += state_add.block<3,1>(3,0);
        // std::cout << "why stuck here?" << state_add.block<3,1>(3,0) << std::endl;
        this->vel_end += state_add.block<3,1>(6,0);
        this->omg     += state_add.block<3,1>(9,0);
        this->acc     += state_add.block<3,1>(12,0);
        this->gravity += state_add.block<3,1>(15,0);
        // bias could be updated as well due to the existence of covariance.

        this->bias_g  += state_add.block<3,1>(18,0);
        this->bias_a  += state_add.block<3,1>(21,0);

        this->yaw_enu_local += state_add(24);
        this->dt_g += state_add(25);
        this->dt_r += state_add(26);
        this->dt_e += state_add(27);
        this->dt_c += state_add(28);
        this->ddt += state_add(29);

        this->anc += state_add.block<3,1>(30,0);
	};

    Matrix<double, DIM_STATE_GNSS+9, 1> operator-(const StatesGroupwithGNSS2& b)
	{
        Matrix<double, DIM_STATE_GNSS+9, 1> a;
        M3D rotd(b.rot_end.transpose() * this->rot_end);
        a.block<3,1>(0,0)  = Log(rotd);
        a.block<3,1>(3,0)  = this->pos_end - b.pos_end;
        a.block<3,1>(6,0)  = this->vel_end - b.vel_end;
        a.block<3,1>(9,0)  = this->omg     - b.omg;
        a.block<3,1>(12,0) = this->acc     - b.acc;
        a.block<3,1>(15,0) = this->gravity - b.gravity;
        a.block<3,1>(18,0) = this->bias_g  - b.bias_g;
        a.block<3,1>(21,0) = this->bias_a  - b.bias_a;
        a(24) = this->yaw_enu_local - b.yaw_enu_local;
        a(25) = this->dt_g - b.dt_g;
        a(26) = this->dt_r - b.dt_r;
        a(27) = this->dt_e - b.dt_e;
        a(28) = this->dt_c - b.dt_c;
        a(29) = this->ddt - b.ddt;

        a.block<3,1>(30,0) = this->anc - b.anc;
		return a;
	};

    void propagate_state(const double & dt)
    {
        this->rot_end = this->rot_end * Exp(this->omg, dt);
        V3D acc_s(this->rot_end * this->acc + this->gravity); 
        this->pos_end += this->vel_end * dt + 0.5 * acc_s * dt * dt;
        this->vel_end += acc_s * dt;
        this->dt_g += this->ddt * dt;
        this->dt_r += this->ddt * dt;
        this->dt_e += this->ddt * dt;
        this->dt_c += this->ddt * dt;
    }

    void propagate_cov(const double & dt) // check if it is right
    {
        M3D && acc_avr_skew  = skew_sym_mat(this->acc);
        M3D && F_vel_wrt_rot = this->rot_end * acc_avr_skew;
        M3D && F_exp = Exp(this->omg, - dt); //A
        M3D && F_vwr = - F_vel_wrt_rot * dt; //B
        M3D && F_rwa = this->rot_end * dt; //C
        M3D Eyedt = Eye3d * dt;
        
        MD(33,33) F_x_base;
        F_x_base.setIdentity();
        F_x_base.block<3,3>(0,0)  = F_exp;
        F_x_base.block<3,3>(0,9)  = Eyedt;
        F_x_base.block<3,3>(3,6)  = Eyedt;
        F_x_base.block<3,3>(6,0)  = F_vwr; //- F_vel_wrt_rot * dt; // B
        F_x_base.block<3,3>(6,12) = F_rwa; //this->rot_end * dt; // C
        F_x_base.block<3,3>(6,15) = Eyedt;
        F_x_base(25,29) = dt;
        F_x_base(26,29) = dt;
        F_x_base(27,29) = dt;
        F_x_base(28,29) = dt;

        double dt_sq = dt * dt;
        
        this->cov = F_x_base * this->cov * F_x_base.transpose();
        this->cov.block<6,6>(18,18).diagonal() += this->Q.block<6,6>(15,15).diagonal() * dt_sq;
        // this->cov.block<6,6>(24,24).diagonal() += this->Q.block<6,6>(6,6).diagonal() * dt_sq;
        this->cov.block<9,9>(6,6).diagonal()   += this->Q.block<9,9>(6,6).diagonal() * dt_sq;
        this->cov.block<6,6>(0,0).diagonal()   += this->Q.block<6,6>(0,0).diagonal() * dt_sq;
        this->cov.block<3,3>(30,30).diagonal()   += this->Q.block<3,3>(22,22).diagonal() * dt_sq;
        this->cov(25,25) += Q(21,21) * dt_sq;
        this->cov(26,26) += Q(21,21) * dt_sq;
        this->cov(27,27) += Q(21,21) * dt_sq;
        this->cov(28,28) += Q(21,21) * dt_sq;
        this->cov(29,29) += Q(21,21) * dt_sq;
        this->cov(24,24) += Q(25,25) * dt_sq;
        // this->cov.block<6,6>(9,9).diagonal()   += this->Q.block<6,6>(0,0).diagonal() * dt_sq;
    }

    StatesGroupwithGNSS2& copy_state_wo_gnss(const StatesGroup2& b)
	{
        this->rot_end = b.rot_end;
		this->pos_end = b.pos_end;
        this->vel_end = b.vel_end;
        this->acc = b.acc;
        this->omg = b.omg;
        this->bias_g  = b.bias_g;
        this->bias_a  = b.bias_a;
        this->gravity = b.gravity;
        this->cov.block<21,21>(0,0) = b.cov; // wrong
        return *this;
	};

    void set_cov_init(VD(DIM_STATE_GNSS+9) diag_init)
    {
        this->cov.setIdentity();
        this->cov.diagonal() = diag_init;
    }

    void set_Q(VD(26) Q_init_diag)
    {
        this->Q.setIdentity();
        this->Q.diagonal() = Q_init_diag;
    }

    void set_Q_ddt(double Q_ddt)
    {
        this->Q(21, 21) = Q_ddt;
    }

	M3D rot_end;      // the estimated attitude (rotation matrix) at the end lidar point
    V3D pos_end;      // the estimated position at the end lidar point (world frame)
    V3D vel_end;      // the estimated velocity at the end lidar point (world frame)
    V3D omg;          // the estimated angular velocity
    V3D acc;          // the estimated acceleration
    V3D gravity;      // the estimated gravity acceleration
    V3D bias_g;       // gyroscope bias
    V3D bias_a;       // accelerator bias
    double yaw_enu_local;
    double dt_g;
    double dt_r;
    double dt_e;
    double dt_c;
    double ddt;
    V3D anc;
    MD(DIM_STATE_GNSS+9, DIM_STATE_GNSS+9)  cov;     // states covariance
    MD(26, 26)  Q;    // Process noise
};


struct StatesGroupwithGNSS
{
    StatesGroupwithGNSS() {
		this->rot_end = M3D::Identity();
		this->pos_end = Zero3d;
        this->vel_end = Zero3d;
        this->bias_g  = Zero3d;
        this->bias_a  = Zero3d;
        this->gravity = Zero3d;
        this->yaw_enu_local = 0.0;
        this->dt_g = 0.0;
        this->dt_r = 0.0;
        this->dt_e = 0.0;
        this->dt_c = 0.0;
        this->ddt = 0.0;
        this->anc = Zero3d;
        this->cov     = MD(DIM_STATE_GNSS+3,DIM_STATE_GNSS+3)::Identity() * INIT_COV;
        this->cov.block<9,9>(9,9) = MD(9,9)::Identity() * 0.00001;
	};

    StatesGroupwithGNSS(const StatesGroupwithGNSS& b) {
		this->rot_end = b.rot_end;
		this->pos_end = b.pos_end;
        this->vel_end = b.vel_end;
        this->bias_g  = b.bias_g;
        this->bias_a  = b.bias_a;
        this->gravity = b.gravity;
        this->cov     = b.cov;

        this->yaw_enu_local = b.yaw_enu_local;
        this->dt_g = b.dt_g;
        this->dt_r = b.dt_r;
        this->dt_e = b.dt_e;
        this->dt_c = b.dt_c;
        this->ddt = b.ddt;
        
        this->anc = b.anc;
	};

    StatesGroupwithGNSS& operator=(const StatesGroupwithGNSS& b)
	{
        this->rot_end = b.rot_end;
		this->pos_end = b.pos_end;
        this->vel_end = b.vel_end;
        this->bias_g  = b.bias_g;
        this->bias_a  = b.bias_a;
        this->gravity = b.gravity;
        this->cov     = b.cov;

        this->yaw_enu_local = b.yaw_enu_local;
        this->dt_g = b.dt_g;
        this->dt_r = b.dt_r;
        this->dt_e = b.dt_e;
        this->dt_c = b.dt_c;
        this->ddt = b.ddt;

        this->anc = b.anc;
        return *this;
	};

    StatesGroupwithGNSS& copy_from(const StatesGroupwithGNSS2& b)
	{
        this->rot_end = b.rot_end;
		this->pos_end = b.pos_end;
        this->vel_end = b.vel_end;
        this->bias_g  = b.bias_g;
        this->bias_a  = b.bias_a;
        this->gravity = b.gravity;

        this->yaw_enu_local = b.yaw_enu_local;
        this->dt_g = b.dt_g;
        this->dt_r = b.dt_r;
        this->dt_e = b.dt_e;
        this->dt_c = b.dt_c;
        this->ddt = b.ddt;

        this->anc = b.anc;
        // this->cov     = b.cov.block<DIM_STATE,DIM_STATE>(0,0);
        // this->Q       = b.Q.block<12, 12>(0, 0);
        return *this;
	};

    StatesGroupwithGNSS& copy_state_wo_gnss(const StatesGroup& b)
	{
        this->rot_end = b.rot_end;
		this->pos_end = b.pos_end;
        this->vel_end = b.vel_end;
        this->bias_g  = b.bias_g;
        this->bias_a  = b.bias_a;
        this->gravity = b.gravity;
        this->cov.block<15,15>(0,0) = b.cov;
        return *this;
	};

    StatesGroupwithGNSS operator+(const Matrix<double, DIM_STATE_GNSS+3, 1> &state_add)
	{
        StatesGroupwithGNSS a;
		a.rot_end = this->rot_end * Exp(state_add(0,0), state_add(1,0), state_add(2,0));
		a.pos_end = this->pos_end + state_add.block<3,1>(3,0);
        a.vel_end = this->vel_end + state_add.block<3,1>(6,0);
        a.bias_g  = this->bias_g  + state_add.block<3,1>(9,0);
        a.bias_a  = this->bias_a  + state_add.block<3,1>(12,0);
        a.gravity = this->gravity + state_add.block<3,1>(15,0);
        
        a.yaw_enu_local = this->yaw_enu_local + state_add(18);
        a.dt_g = this->dt_g + state_add(19);
        a.dt_r = this->dt_r + state_add(20);
        a.dt_e = this->dt_e + state_add(21);
        a.dt_c = this->dt_c + state_add(22);
        a.ddt = this->ddt + state_add(23);

        a.anc = this->anc + state_add.block<3,1>(24,0);
        a.cov     = this->cov;
		return a;
	};

    StatesGroupwithGNSS& operator+=(const Matrix<double, DIM_STATE_GNSS+3, 1> &state_add)
	{
        this->rot_end = this->rot_end * Exp(state_add(0,0), state_add(1,0), state_add(2,0));
		this->pos_end += state_add.block<3,1>(3,0);
        this->vel_end += state_add.block<3,1>(6,0);
        this->bias_g  += state_add.block<3,1>(9,0);
        this->bias_a  += state_add.block<3,1>(12,0);
        this->gravity += state_add.block<3,1>(15,0);
        this->yaw_enu_local += state_add(18);
        this->dt_g += state_add(19);
        this->dt_r += state_add(20);
        this->dt_e += state_add(21);
        this->dt_c += state_add(22);
        this->ddt += state_add(23);
        this->anc += state_add.block<3,1>(24,0);
		return *this;
	};

    Matrix<double, DIM_STATE_GNSS+3, 1> operator-(const StatesGroupwithGNSS& b)
	{
        Matrix<double, DIM_STATE_GNSS+3, 1> a;
        M3D rotd(b.rot_end.transpose() * this->rot_end);
        a.block<3,1>(0,0)  = Log(rotd);
        a.block<3,1>(3,0)  = this->pos_end - b.pos_end;
        a.block<3,1>(6,0)  = this->vel_end - b.vel_end;
        a.block<3,1>(9,0)  = this->bias_g  - b.bias_g;
        a.block<3,1>(12,0) = this->bias_a  - b.bias_a;
        a.block<3,1>(15,0) = this->gravity - b.gravity;
        a(18) = this->yaw_enu_local - b.yaw_enu_local;
        a(19) = this->dt_g - b.dt_g;
        a(20) = this->dt_r - b.dt_r;
        a(21) = this->dt_e - b.dt_e;
        a(22) = this->dt_c - b.dt_c;
        a(23) = this->ddt - b.ddt;
        a.block<3,1>(24,0) = this->anc - b.anc;
		return a;
	};

    void resetpose()
    {
        this->rot_end = M3D::Identity();
		this->pos_end = Zero3d;
        this->vel_end = Zero3d;
    }

    void set_cov_init(VD(DIM_STATE_GNSS+3) diag_init)
    {
        this->cov.setIdentity();
        this->cov.diagonal() = diag_init;
    }

    void set_Q(VD(17) Q_init_diag)
    {
        this->Q.setIdentity();
        this->Q.diagonal() = Q_init_diag;
    }

    void set_Q_ddt(double Q_ddt)
    {
        this->Q(12,12) = Q_ddt;
    }

	M3D rot_end;      // the estimated attitude (rotation matrix) at the end lidar point
    V3D pos_end;      // the estimated position at the end lidar point (world frame)
    V3D vel_end;      // the estimated velocity at the end lidar point (world frame)
    V3D bias_g;       // gyroscope bias
    V3D bias_a;       // accelerator bias
    V3D gravity;      // the estimated gravity acceleration
    double yaw_enu_local;
    double dt_g;
    double dt_r;
    double dt_e;
    double dt_c;
    double ddt;
    V3D anc;
    Matrix<double, DIM_STATE_GNSS+3, DIM_STATE_GNSS+3>  cov;     // states covariance
    MD(17, 17)  Q;
};

#endif