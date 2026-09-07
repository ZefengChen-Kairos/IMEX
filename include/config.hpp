#ifndef CONFIG_HPP
#define CONFIG_HPP

// Compile-time switches, physical parameters and global solver state.
// Globals are C++17 inline variables: one definition shared by all TUs.

#include <iostream>
#include <vector>
#include <cmath>
#include <float.h>
#include <fstream>
#include <ctime>
#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <functional>
#include <iomanip>
#include <random>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>
#include <Eigen/Core>
#include <filesystem>
#include <sstream>


#define NUM_VARS 8
enum Formula{EnergyFormula, PressureFormula};
enum MagField{StaggerCT, UnStaggerCT, PotentialA, NoDivFree};
enum Test{test0, test1, test2, test3, test4, test5, test6, test7,test8,test9,test10, shocktest, convergencetest,test2d0, test1y, testOrszag_Tang, testblastwave, testblastwave1D,testRotor, testFieldLoop, testKHInstability, testdoubleshear, testRK_Riemann, testPlasmaEquilibrium, testAdvectedPlasmaEquilibrium, testMHD_vortex, testblastwave2, testShockCloud, testMHDJet,RKEOS1,RKEOS2,RKEOS3,RKEOS4};
enum Reconstruction{NoRecon, Minmod, VanLeer, ENO3, WENO5};
enum shock_threshold_Left{one, val_left};
enum shock_threshold_Right{zero, val_right};
enum shock_detector{Yes_add, No_First, No_Second};
enum SolvePressure{no_correction, once_correction, all_correction};
enum dtCalculation{Fullydt, Alfvendt, IMEXdt};
enum x_boundary_condition{x_Transmissive, x_Reflective, x_Periodic, x_Dirichlet};
enum y_boundary_condition{y_Transmissive, y_Reflective, y_Periodic, y_Dirichlet, y_JetInflow};
enum Add_Pvisc{MaterialSpeedp, AlfvenSpeedp, Noaddvisp};
enum Energy_Update{ConservationofU, DefinitionofU};
enum Flux_Diss{Local_Diss, Global_Diss};
enum ReconState{PrimState, ConsState};

inline dtCalculation caltimestep = IMEXdt;
inline Formula formula = PressureFormula;
inline MagField magfield = NoDivFree;
inline Add_Pvisc add_pvisc = MaterialSpeedp;
inline Flux_Diss flux_diss = Global_Diss;
inline bool GLF = true;
inline bool add_bvisc = true;    // 磁 A 子系统数值粘性:论文一直需要;系数=面平均快磁声速(对齐 CASTOR add_bvisc=.TRUE.)
inline Energy_Update U_update = DefinitionofU;
inline Test test = testMHD_vortex;
const int nxCells = 10;
const int nyCells = 10;
const int nzCells = 1;
const int ghost_cells = 3;
const int RKs = 2;
inline double CFL = 0.8;
const int full_step = 1;     // 对齐 CASTOR:特征值 magnetosonic 只在 step 1,之后材料速度(硬编码 timestep==1)
inline double A0 = 1e-3;
inline Reconstruction recon = Minmod;
inline ReconState reconstate = ConsState;
inline int SpaceOrder;
inline x_boundary_condition x_bc;
inline y_boundary_condition y_bc;
inline SolvePressure PressureSolver = all_correction;

inline bool add_RHSpvisc = false;
inline bool add_RHSBvisc = false;
inline bool Energy_Modify = false;

inline double GMREStol = 1e-14;
inline double GMRESiter = 1e32;
inline int GMRESrest = 30;
inline bool add_restart = false;

inline shock_detector add_shock_detector = No_Second;
inline double shock_kappa1 = 0.001;
inline double shock_threshold = 0.9;
inline shock_threshold_Left threshold_Left = val_left;
inline shock_threshold_Right threshold_Right = zero;

inline int nxPoints = nxCells + 1, nyPoints = nyCells + 1, nzPoints = nzCells + 1;
inline double t, time_stop, dt, adt, adt2, fully_dt, imex_dt;
inline int timestep = 0;
inline double x_0 = -0.5, x_1 = 0.5, y_0 = -1, y_1 = 1,z_0 = -1, z_1 = 1, B_max_x = 0, B_max_y = 0, U_max_x = 0,U_max_y = 0,a_max_x = 0,a_max_y = 0, a_max = 0;
inline std::string Path;
// Root directory for all output. Relative to the working directory by default;
// override at run time with the MHD_OUTPUT_ROOT environment variable or argv[1].
inline std::string output_root = "output";

inline double dx = (x_1 - x_0)/(nxPoints-1);
inline double dy = (y_1 - y_0)/(nyPoints-1);
inline double dz = (z_1 - z_0)/(nzPoints-1);
inline double ch_x,ch_y, cp_x, cp_y, cA_x, cA_y, cA, vx_max, vy_max, ch, v_max;


// ===================== 压力保正:Casulli 嵌套 Newton(piecewise-linear V) =====================
// 思路:把 pressure subsystem 的对角项 V(p)=p/(gamma-1) 在 p<p_floor 处改成大斜率 Mbig 的直线,
//       于是 V 全局非线性 -> 用 Newton 外层 + GMRES 内层求解;大斜率把 p 弹回 p_floor 之上。
inline bool   floor_on   = true;          // 总开关:false 则完全退回你原来的线性单次 GMRES
inline double p_floor    = 1e-7;          // 压力地板
inline double Mbig_fac   = 1e6;           // 大斜率倍率(相对 1/(gamma-1));GMRES 病态就调小
inline double theta = 0.9;   // 正性目标：保证 p >= theta * p_floor (Cor. 1)
inline double Mbig  = 0.0;   // 本次 solve 实际用的陡坡斜率——每次 solve 重算
// 必须和 p_freeze 同作用域，AX 要读它
inline int    Pnewton_max = 20;           // Newton 最大迭代
inline double Pnewton_tol = 1e-8;        // Newton 残差容差
inline Eigen::VectorXd p_freeze;          // 当前 Newton 迭代值,供 Jacobian(=AX)读取斜率
inline int    g_max_floored = 0;          // 全程单次 solve 内 floored 格点数的最大值(§4 统计用)
inline double g_min_p_ever  = 1e300;      // 全程最小压力(校验 p >= theta*p_floor)


// ============================================================================================
// =====================  EOS 抽象层(Ideal Gas / Redlich-Kwong)  =============================
// 论文核心:EOS 的全部非线性只通过对角项 V(p)=rho*e(rho,p) 及其斜率 dV/dp 进入隐式系统,
// 耦合矩阵 T(焓 Laplacian)与 EOS 无关。因此移植只需把 V(p)、声速、con<->pri 改成 EOS-aware。
// 理想气体路径保持逐字节不变(Vbase_slope==1/(gamma-1), 声速==sqrt(gamma p/rho))。
enum EOSType { IDEAL_GAS, REDLICH_KWONG };
inline EOSType eos_type = IDEAL_GAS;                 // 顶层开关;RK 测试在 main 里自动置 REDLICH_KWONG

// ---- Ideal Gas 参数----
inline double gamma_plasma = 5.0/3.0;
// ---- Redlich-Kwong 参数----
//double m_alfa = 0.5, m_b = 0.5, m_R = 0.4, m_cv = 1.0;
inline double m_alfa = 0.05, m_b = 0.1, m_R = 0.4, m_cv = 1.0;

// ---- RK 热/卡 EOS 的 Newton 残差(由 1D 代码移植)----
//   热 EOS:  p = m_R T/(v-b) - alfa/sqrt(T)/v/(v+b),   v = 1/rho
//   卡 EOS:  e = m_cv T + (3/2) alfa/b * ln(v/(v+b)) / sqrt(T)
inline double RK_fNR_T (double T,double v,double p){ return m_R*T/(v-m_b) - m_alfa/sqrt(T)/v/(v+m_b) - p; }
inline double RK_dfNR_T(double T,double v,double /*p*/){ return m_R/(v-m_b) + 0.5*m_alfa/T/sqrt(T)/v/(v+m_b); }
inline double RK_fNR_p (double T,double v,double e){ double U=log(v/(v+m_b)); return m_cv*T + 1.5*m_alfa/m_b*U/sqrt(T) - e; }
inline double RK_dfNR_p(double T,double v,double /*e*/){ double U=log(v/(v+m_b)); return m_cv - 0.75*m_alfa/m_b*U/T/sqrt(T); }


// 由 (rho,p) 解温度 T
inline double REDKWONG_computeTfromRhoAndP(double rho, double p){
    double v = 1.0/rho;
    double T_old, T_new = (1.0/rho - m_b)*(p + m_alfa*rho*rho)/m_R;
    const double eps = 1e-14;
    do{ T_old = T_new; T_new = T_old - RK_fNR_T(T_old,v,p)/RK_dfNR_T(T_old,v,p);
    }while(fabs(T_new - T_old)/fabs(T_old) > eps);
    return T_new;
}
// e(rho,p)
inline double REDKWONG_e_from_rho_p(double rho, double p){
    double v = 1.0/rho;
    double T = REDKWONG_computeTfromRhoAndP(rho,p);
    double U = log(v/(v+m_b));
    return m_cv*T + 1.5*m_alfa/m_b*U/sqrt(T);
}
// p(rho,e):先由卡 EOS 解 T,再代入热 EOS
inline double REDKWONG_p_from_rho_e(double rho, double e){
    double v = 1.0/rho;
    double T_old, T_new = (e + rho*m_alfa)/m_cv;
    const double eps = 1e-14;
    do{ T_old = T_new; T_new = T_old - RK_fNR_p(T_old,v,e)/RK_dfNR_p(T_old,v,e);
    }while(fabs(T_new - T_old)/fabs(T_old) > eps);
    double a = m_alfa/sqrt(T_new);
    return m_R*T_new/(v-m_b) - a/v/(v+m_b);
}
// RK 声速 a = sqrt(dp/drho|_s)(由 1D REDKWONG_computeSoundSpeedFromEoS 移植)
inline double REDKWONG_soundspeed(double rho, double p){
    double v = 1.0/rho;
    double T = REDKWONG_computeTfromRhoAndP(rho,p);
    double U = log(v/(v+m_b));
    double p_T = m_R/(v-m_b) + m_alfa/(2.0*T*sqrt(T)*v*(v+m_b));
    double e_T = m_cv - 0.75*m_alfa/m_b*U/(T*sqrt(T));
    double p_v = -m_R*T/pow((v-m_b),2) + m_alfa*(2.0*v+m_b)/(sqrt(T)*v*v*pow((v+m_b),2));
    double e_v = 1.5*m_alfa/m_b*(1.0/v - 1.0/(v+m_b))/sqrt(T);
    double a2 = p/(rho*rho)*(p_T/e_T) + v*v*(e_v - e_T*(p_v/p_T));
    return sqrt(a2);
}

// ---- 统一 EOS 接口:其余代码全部经由这些函数,EOS 切换只改 eos_type ----
inline double eos_p_from_rho_e(double rho, double e){
    return (eos_type==IDEAL_GAS) ? (gamma_plasma - 1.0)*rho*e
                                 : REDKWONG_p_from_rho_e(rho, e);
}
inline double eos_e_from_rho_p(double rho, double p){
    return (eos_type==IDEAL_GAS) ? p/((gamma_plasma - 1.0)*rho)
                                 : REDKWONG_e_from_rho_p(rho, p);
}
inline double eos_sound_a(double rho, double p){
    return (eos_type==IDEAL_GAS) ? sqrt(gamma_plasma*p/rho)
                                 : REDKWONG_soundspeed(rho, p);
}
inline double eos_sound_a2(double rho, double p){
    double a = eos_sound_a(rho, p); return a*a;
}
// 快磁声速 cf:用 a^2 通用表达,Bn 为法向磁场分量
inline double eos_fast_speed(double rho, double p, double Bn, double B2){
    double a2 = eos_sound_a2(rho, p);
    double t  = a2 + B2/rho;
    double disc = t*t - 4.0*a2*Bn*Bn/rho;     // 4 gamma p Bn^2/rho^2 = 4 a^2 Bn^2/rho
    if(disc < 0) disc = 0;
    return sqrt(0.5*(t + sqrt(disc)));
}

// ===================== 压力保正:Casulli 嵌套 Newton(piecewise V) =====================
// V(p) = rho*e(rho,p) 为隐式压力系统的对角项。理想气体 V(p)=p/(gamma-1) 不依赖 rho;
// RK 的 V(p)=rho*e(rho,p) 依赖 rho,故 V_mod/V_slope 改为 rho-aware。
// 保正:在物理斜率之下设一个大斜率 M(p<p_floor 切换),把 p 弹回 p_floor 之上(论文 (58)-(59))。

// 物理 V(p)=rho*e(rho,p)
inline double Vbase(double p, double rho){
    return (eos_type==IDEAL_GAS) ? p/(gamma_plasma - 1.0)
                                 : rho*REDKWONG_e_from_rho_p(rho, p);
}
// 物理 dV/dp = rho*de/dp:理想气体解析,RK 用中心差分(对齐 1D, small=1e-8)
inline double Vbase_slope(double p, double rho){
    if(eos_type==IDEAL_GAS) return 1.0/(gamma_plasma - 1.0);
    const double d = 1e-8;
    double ep = REDKWONG_e_from_rho_p(rho, p + d);
    double em = REDKWONG_e_from_rho_p(rho, p - d);
    return rho*(ep - em)/(2.0*d);
}
//inline double V_mod(double p, double rho){                  // 带保正地板的 V(p)
//  if(!floor_on) return Vbase(p, rho);
//  double M = Mbig_fac/(gamma_plasma - 1.0);               // 大斜率(与理想气体一致的量纲)
//  return (p >= p_floor) ? Vbase(p, rho)
//                        : Vbase(p_floor, rho) + M*(p - p_floor);
//}
//
//inline double V_slope(double p, double rho){                // dV/dp(Newton Jacobian 对角)
//  if(!floor_on) return Vbase_slope(p, rho);
//  double M = Mbig_fac/(gamma_plasma - 1.0);
//  return (p >= p_floor) ? Vbase_slope(p, rho) : M;
//}

inline double V_mod (double p, double rho) {
	if(!floor_on) return Vbase(p, rho);
	return (p >= p_floor) ? Vbase(p, rho)
	                      : Vbase(p_floor, rho) + Mbig * (p - p_floor);
}
inline double V_slope(double p, double rho) {
	if(!floor_on) return Vbase_slope(p, rho);
	return (p >= p_floor) ? Vbase_slope(p, rho) : Mbig;
}
// AX 和残差里读斜率统一走 V_slope(p_freeze[i], rho[i]);
// Mbig 由 Pressure_NestedNewton 每次求解前按 Cor.1 重算。

// AX 算子模式开关:false=线性载体算子 L(对角=1/(gamma-1));true=Newton Jacobian(对角=V_slope)
inline bool g_AX_jacobian = false;

// --- Dirichlet ambient state (blast wave); E/A computed with gamma_plasma at runtime ---
inline double rho_amb = 1.0, p_amb = 0.1, Bx_amb = 70.0/sqrt(4.0*M_PI);
inline bool g_bc_prim = false;   // true while bound_cond_u is applied to a PRIMITIVE grid (slot 4 = p, not E)
inline double g_fixed_dt = 0.0;   // 0 = 用真实 IMEX dt(生产);>0 仅为对 CASTOR 强制固定 dt
inline double E_amb_val(){ return p_amb/(gamma_plasma - 1.0) + 0.5*Bx_amb*Bx_amb; } // ambient total energy
// ============================================================================================

// ===================== MHD Jet (Example 4.8, Mach 800/2000/10000) =====================
// 底边 |x|<0.05 为 inflow strip(Dirichlet jet state),strip 外及其余三边为 outflow(透射)。
// AX 隐式侧:outflow 由 M_P/M_B 的 clamp 天然给出齐次 Neumann,不需额外处理;
// inflow strip 采用 ghost-cell Dirichlet:算子中 ghost 未知量取 0(齐次),已知值 lift 到 RHS。
inline double jet_B0  = sqrt(20000.0);   // sqrt(200)/sqrt(2000)/sqrt(20000)  <-> beta = 1e-2/1e-3/1e-4
inline double jet_v   = 10000.0;         // 800 / 2000 / 10000
inline double jet_p   = 1.0;
inline double jet_rho = 1.4;           // = gamma(在 test 配置里随 gamma_plasma 赋值)
inline bool jet_inlet(int i){                 // cell-center 判定: |x_c| < 0.05
	if(y_bc != y_JetInflow) return false;
	double xc = x_0 + (i - ghost_cells + 0.5) * dx;
	return fabs(xc) < 0.05;
}
inline bool jet_inlet_face(int i){            // x 面 (i+1/2) 判定,用于 staggered Bx_p0
	if(y_bc != y_JetInflow) return false;
	double xf = x_0 + (i - ghost_cells + 1.0) * dx;
	return fabs(xf) < 0.05 - 1e-12;
}
inline double jet_E_val(){                    // jet 守恒总能(含动能!v=800 时 KE~4.5e5)
	return jet_p/(gamma_plasma - 1.0) + 0.5*jet_rho*jet_v*jet_v + 0.5*jet_B0*jet_B0;
}
// ============================================================================================

#endif // CONFIG_HPP
