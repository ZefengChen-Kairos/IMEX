#ifndef DECLARATIONS_HPP
#define DECLARATIONS_HPP

// Forward declarations shared by all translation units.

#include "types.hpp"

using namespace std;
double getAlfvenSpeed(double *pri);
void copystate(double *state_1, double *state_0);
void copy(Grid &state_1, Grid &state_0);
void f_convective(double *f, double *con);
void f_pressure(double *f, double *con);
void f_magnetic(double *f, double *con);
void g_convective(double *g, double *con);
void g_pressure(double *g, double *con);
void g_magnetic(double *g, double *con);
void f_full_con(double *f, double *con);
void g_full_con(double *g, double *con);
void ini_con2pri(double *pri, double *con);
void ini_pri2con(double *con, double *pri);
void upri2ucon(Grid &u_con, Grid &u_pri);
void ucon2upri(Grid &u_pri,Grid &u_con);
void half_time(double *halfL, double *halfR, double *L, double *R);
double shock_filter(double val);
void MagneticFlux_x(double *Q, double *flux);
void MagneticFlux_y(double *Q, double *flux);
void f_KineticEnergyFlux(double *Q, double *flux);
void g_KineticEnergyFlux(double *Q, double *flux);
void MagneticFlux_Update_Energy(Grid &u_con, Grid &u_test);
void MagneticFlux_Pressure_x(double *Q, double *flux);
void MagneticFlux_Pressure_y(double *Q, double *flux);
void MagneticFlux_Update_Pressure(Grid &u_con, Grid &u_test, Grid &Qe_pri);
void RusanovK_x(double *QL, double *QR, double *flux);
void RusanovK_y(double *QL, double *QR, double *flux);
void RusanovK_Update(Grid &u_con, Grid &Qs, Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR);
void RusanovK_Pressure_x(double *QL, double *QR, double *flux);
void RusanovK_Pressure_y(double *QL, double *QR, double *flux);
void RusanovK_Update_Pressure(Grid &u_con, Grid &Qs, Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR);
double ini_calch_x(double *pri);
double ini_calch_y(double *pri);
void RusanovC_x(double *QL, double *QR, double *flux);
void RusanovC_y(double *QL, double *QR, double *flux);
void RusanovC_Update(Grid &u_con, Grid &Qs, Grid &w_xL, Grid &w_xR, Grid &w_yL, Grid &w_yR,double h);
void TVDReconstruction(Grid &u_con, Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR);
double calecon(double *con);
double calepri(double *pri);
double calcf(double* pri);
vector<double> ConvergenceTestFunction(double x);
vector<double> ConvergenceTestExact(double x);
vector<double> GaussLegendreNumericalIntegration(const double& xc, const double& dx, function<vector<double>(double)> func, int nGaussPoints);
void u_initialization(Grid &u, Array &Bx_p0, Array &By_0p);
void bound_cond_u(Grid &u);
void bound_cond_A(Grid &A);
double caldt(Grid &u_pri);
void update(Grid &u_con,Grid &flux, char dir);
void calch(Grid &u_pri);
void calca(Grid &u_pri);
double ini_calch(double *u_pri);
void calcp(Grid &u_pri);
void store(Grid &u_con, Array& vxn,Array& vyn, Array& Bxn, Array& Byn);
void bound_cond_array(Array &u);
double calenthalpy(double *con);
double calkinetic(double *con);
double calmagenergy(double *con);
void MagneticEigenvalues(Grid &u_pri);
void EnergyEigenvalues(Grid &u_pri);
void Staggered_CT(Grid& u_con, Array &Bx_p0, Array &By_0p, double h);
void A_initial(Grid &A);
inline double minmod(double a, double b);
void Unstaggered_CT(Grid& u_con, Grid& A,Grid& A_s, double h);
int M_P(int i, int j);
int M_B(int i, int j, int k);
void Magnetic_Subsystem_AX(Eigen::VectorXd& Ax, Eigen::VectorXd& x, Grid& u_con, Grid& Qe,Grid& Qs,double h);
void Magnetic_Subsystem_potential_AX(Eigen::VectorXd& Ax_out, Eigen::VectorXd& x,Grid& u_con, Grid& Qe, Grid& Qs, double h);
void Pressure_Subsystem_AX(Eigen::VectorXd& Ax, Eigen::VectorXd& x, Grid& u_con, Grid& Qe,Grid& Qs,double h);
void Pressure_NestedNewton(Eigen::VectorXd& p, Eigen::VectorXd& b, Grid& u_con, Grid& Qe, Grid& Qs, double h, const int dim);
void Energy_Subsystem_AX(Eigen::VectorXd& Ax, Eigen::VectorXd& x, Grid& u_con, Grid& Qe,Grid& Qs,double h);
void GMRES_Solver2(Eigen::VectorXd&x,Eigen::VectorXd& b, std::function<void(Eigen::VectorXd&, Eigen::VectorXd&, Grid&, Grid&, Grid&,double)> matrix_vector_product, Grid& u_con, Grid& Qe, Grid& Qs, double h, const int dim);
void GMRES_Solver(Eigen::VectorXd &x, Eigen::VectorXd &b, std::function<void(Eigen::VectorXd&, Eigen::VectorXd&, Grid&, Grid&, Grid&, double)> Ax_func,Grid &u_con, Grid &Qe, Grid &Qs, double h,const int dim, std::function<void(Eigen::VectorXd&, const Eigen::VectorXd&)> precond = {});
double htilde(double hL, double hR, double mL, double mR);
void PreMagSubsystem_Energy(Grid &u_con, Grid &Qs, Grid &Qe,Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR, Array &Bx_p0, Array &By_0p, double h, const int rks);
void MagSubsystem_UpdateA(Grid &u_con, Grid &Qs, Grid &Qe,Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR, Grid &A,Grid &A_s,double h, const int rks);
void PreSubsystem_Energy(Grid &u_con, Grid &Qs, Grid &Qe,Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR, double h, const int rks);
void PreMagSubsystem_Energy_potential(Grid &u_con, Grid &Qs, Grid &Qe,Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR, Array &Bx_p0, Array &By_0p, Grid &A,double h, const int rks);
void PreSubsystem_Pressure(Grid &u_con, Grid &Qs, Grid &Qe,Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR, double h, const int rks);
void PreMagSubsystem_Pressure(Grid &u_con, Grid &Qs, Grid &Qe,Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR, Array &Bx_p0, Array &By_0p, Array &A,double h, const int rks);
void PreMagSubsystem_Pressure_potential(Grid &u_con, Grid &Qs, Grid &Qe,Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR, Array &Bx_p0, Array &By_0p, Grid &A,double h, const int rks);
void MagSubsystem_UpdateB(Grid &u_con, Grid &Qs, Grid &Qe,Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR, Grid &A,double h, const int rks);
double cal_divergence_error(Grid &u_con, Array &Bx_0p, Array &By_p0);
void IntermediatePlotIntermediatePlot(Grid& u_con, Grid &u_pri,  const int iteration, const int plotINTERVAL, const std::vector<double>& IntermediatePlotTimes, bool isFinalTime);
void printVTKHeader(std::ofstream& outFile, int& nstep);
void printVTKVariable(std::ofstream& outFile, Grid& u_con, Grid& u_pri, int var, std::string variableName);
void printAllVariablesToVTK(Grid& u_con,Grid& u_pri, std::string testname, int nstep, bool isFinalTime);
void Staggered_CT_Consistent(Grid& u_con, Array &Bx_p0, Array &By_0p,Array &Bx_p0_Qs, Array &By_0p_Qs, double h);
void CollectMagFluxes(Grid &u_con, Grid &Qe, Grid &Qs, double h);
void bound_cond_staggered(Array &u, int axis);




// --- declarations for definitions that had none in the original single file ---
void   ShockDetector(Array &shock_eta, Grid &u_con);
void   Energy_flux_x(double *QL, double *QR, double *flux);
void   Energy_flux_y(double *QL, double *QR, double *flux);
void   Energy_flux_Update(Grid &u_con, Grid &Qs, Grid &w_xL, Grid &w_xR, Grid &w_yL, Grid &w_yR, double h);
void   Pressure_flux_x(double *QL, double *QR, double *flux);
void   Pressure_flux_y(double *QL, double *QR, double *flux);
void   Pressure_flux_Update(Grid &u_con, Grid &Qs, Grid &w_xL, Grid &w_xR, Grid &w_yL, Grid &w_yR, double h);
double caltotalenthalpy(double *con);
void   Staggered_CT(Grid& u_con, Array &Bx_p0, Array &By_0p, Array &Bx_p0_Qs, Array &By_0p_Qs, double h);
void   CT_CornerAssemble(Array &Ez_pp, Array &Ez_p0, Array &Ez_0p, Array &Ez_tilde);
void   CT_EMF_Stats(const char* name, Array &Ez_pp);
void   Diagnose_CT_EMF(Grid& u_con);
void   Unstaggered_CT_Consistent(Grid& u_con, Grid& A, Grid& A_s, double h);

#endif // DECLARATIONS_HPP
