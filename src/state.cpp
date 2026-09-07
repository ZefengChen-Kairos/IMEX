// Primitive/conservative conversions, physical fluxes, wave speeds.

#include "declarations.hpp"

double getAlfvenSpeed(double *pri){
	double rho = pri[0];
	double cA = pri[5]  / sqrt(rho);
	return cA;
}

void copystate(double *state_1, double *state_0){
	for(int i = 0; i < NUM_VARS; i++){
		state_1[i] = state_0[i];
	}
}

void copy(Grid &state_1, Grid &state_0){
	for(int i = 0; i < nxCells + 2 * ghost_cells; i++){
		for(int j = 0; j < nyCells + 2 * ghost_cells; j++){
			for(int k = 0; k < NUM_VARS; k++){
				state_1(i,j,k) = state_0(i,j,k);
			}
		}
	}
}

void f_convective(double *f, double *con){
	// calculate f
	double rho = con[0];
	double vx = con[1] / con[0];
	double vy = con[2] / con[0];
	double vz = con[3] / con[0];
	double pri[8];
	ini_con2pri(pri, con);
	double k = calkinetic(con);
	f[0] = rho * vx;
	f[1] = rho * vx * vx;
	f[2] = rho * vy * vx;
	f[3] = rho * vz * vx;
	f[4] = k * vx;
	f[5] = 0;
	f[6] = 0;
	f[7] = 0;
}

void f_pressure(double *f, double *con){
	// calculate f
	double vx = con[1] / con[0];
	double U = con[4];
	double pri[8];
	ini_con2pri(pri, con);
	double p = pri[4];
	f[0] = 0;
	f[1] = p;
	f[2] = 0;
	f[3] = 0;
	f[4] = (U + p) * vx;
	f[5] = 0;
	f[6] = 0;
	f[7] = 0;
}

void f_magnetic(double *f, double *con){
	// calculate f
	double vx = con[1] / con[0];
	double vy = con[2] / con[0];
	double vz = con[3] / con[0];
	double Bx = con[5];
	double By = con[6];
	double Bz = con[7];
	double pri[8];
	ini_con2pri(pri, con);
	double B_square = Bx * Bx + By * By + Bz * Bz;
	double vB = vx * Bx + vy * By + vz * Bz;
	f[0] = 0;
	f[1] = 0.5 * B_square - Bx * Bx;
	f[2] = -By * Bx;
	f[3] = -Bz * Bx;
	f[4] = B_square * vx - vB * Bx;
	f[5] = 0;
	f[6] = By * vx - vy * Bx;
	f[7] = Bz * vx - vz * Bx;
}


void g_convective(double *g, double *con){
	// calculate f
	//	// calculate g
	//	// based on conservative value roh rohv E
	double rho = con[0];
	double vx = con[1] / con[0];
	double vy = con[2] / con[0];
	double vz = con[3] / con[0];
	double pri[NUM_VARS];
	ini_con2pri(pri, con);
	
	double k = calkinetic(con);
	g[0] = rho * vy;
	g[1] = rho * vx * vy;
	g[2] = rho * vy * vy;
	g[3] = rho * vz * vy;
	g[4] = k * vy;
	g[5] = 0;
	g[6] = 0;
	g[7] = 0;
}

void g_pressure(double *g, double *con){
	// calculate f
	//	// calculate g
	//	// based on conservative value roh rohv E
	double vy = con[2] / con[0];
	double U = con[4];
	double pri[NUM_VARS];
	ini_con2pri(pri, con);
	double p = pri[4];
	
	g[0] = 0;
	g[1] = 0;
	g[2] = p;
	g[3] = 0;
	g[4] = (U + p) * vy;
	g[5] = 0;
	g[6] = 0;
	g[7] = 0;
}

void g_magnetic(double *g, double *con){
	// calculate f
	double vx = con[1] / con[0];
	double vy = con[2] / con[0];
	double vz = con[3] / con[0];
	double Bx = con[5];
	double By = con[6];
	double Bz = con[7];
	double pri[NUM_VARS];
	ini_con2pri(pri, con);
	double B_square = Bx * Bx + By * By + Bz * Bz;
	double vB = vx * Bx + vy * By + vz * Bz;
	
	
	g[0] = 0;
	g[1] = -Bx * By;
	g[2] = 0.5 * B_square - By * By;
	g[3] = -Bz * By;
	g[4] = B_square * vy - vB * By;
	g[5] = Bx * vy - vx * By;
	g[6] = 0;
	g[7] = Bz * vy - vz * By;
}

void f_full_con(double *f, double *con){
	double rho = con[0];
	double vx = con[1] / con[0];
	double vy = con[2] / con[0];
	double vz = con[3] / con[0];
	double Bx = con[5];
	double By = con[6];
	double Bz = con[7];
	double pri[NUM_VARS];
	ini_con2pri(pri, con);
	double p = pri[4];
	double B_square = Bx * Bx + By * By + Bz * Bz;
	double vB = vx * Bx + vy * By + vz * Bz;
	double kinetic_energy = 0.5 * rho * (vx * vx + vy * vy + vz * vz);
	double e = calecon(con);
	double enthalpy = e + p/rho;
	
	f[0] = rho * vx;
	f[1] = rho * vx * vx + p + 0.5 * B_square - Bx * Bx;
	f[2] = rho * vy * vx - By * Bx;
	f[3] = rho * vz * vx - Bz * Bx;
	f[4] = (kinetic_energy + B_square) * vx - vB * Bx + enthalpy * rho * vx;
	f[5] = 0;
	f[6] = By * vx - vy * Bx;
	f[7] = Bz * vx - vz * Bx;
}

void g_full_con(double *g, double *con){
	// calculate g
	// based on conservative value roh rohv E
	double rho = con[0];
	double vx = con[1] / con[0];
	double vy = con[2] / con[0];
	double vz = con[3] / con[0];
	double U = con[4];
	double Bx = con[5];
	double By = con[6];
	double Bz = con[7];
	double pri[NUM_VARS];
	ini_con2pri(pri, con);
	double p = pri[4];
	double B_square = Bx * Bx + By * By + Bz * Bz;
	double vB = vx * Bx + vy * By + vz * Bz;
	g[0] = rho * vy;
	g[1] = rho * vx * vy - Bx * By;
	g[2] = rho * vy * vy + p + 0.5 * B_square - By * By;
	g[3] = rho * vz * vy - Bz * By;
	g[4] = (U + p + 0.5 * B_square) * vy - vB * By;
	g[5] = Bx * vy - vx * By;
	g[6] = 0;
	g[7] = Bz * vy - vz * By;
}


void ini_con2pri(double *pri, double *con){
	// convert conservative (rho rhov E) to primitive (rho v p)
	// rho v p
	double e = calecon(con);
	pri[0] = con[0];
	pri[1] = con[1] / con[0];
	pri[2] = con[2] / con[0];
	pri[3] = con[3] / con[0];
	pri[4] = eos_p_from_rho_e(con[0], e);
	pri[5] = con[5];
	pri[6] = con[6];
	pri[7] = con[7];
}

void ini_pri2con(double *con, double *pri){
	// convert primitive (rho v p) to conservative  (rho rhov E)
	double e = calepri(pri);
	// rho rhov E
	double vx = pri[1];
	double vy = pri[2];
	double vz = pri[3];
	double Bx = pri[5];
	double By = pri[6];
	double Bz = pri[7];
	double B_square = Bx * Bx + By * By + Bz * Bz;
	double v_square = vx * vx + vy * vy + vz * vz;
	con[0] = pri[0];
	con[1] = pri[0] * pri[1];
	con[2] = pri[0] * pri[2];
	con[3] = pri[0] * pri[3];
	con[4] = pri[0] * e + 0.5 * pri[0] * v_square + 0.5 * B_square;
	con[5] = pri[5];
	con[6] = pri[6];
	con[7] = pri[7];
}

void upri2ucon(Grid &u_con, Grid &u_pri){
	double ini_u_con[NUM_VARS], ini_u_pri[NUM_VARS];
	for(int i = 0; i < nxCells + 2 * ghost_cells; i++){
		for(int j = 0; j < nyCells + 2 * ghost_cells; j++){
			for(int k = 0; k < NUM_VARS; k++){
				ini_u_pri[k] = u_pri(i,j,k);
			}
			ini_pri2con(ini_u_con, ini_u_pri);
			for(int k = 0; k < NUM_VARS; k++){
				u_con(i, j, k) = ini_u_con[k];
			}
		}
	}
}

void ucon2upri(Grid &u_pri,Grid &u_con){
	double ini_u_con[NUM_VARS], ini_u_pri[NUM_VARS];
	for(int i = 0; i < nxCells + 2 * ghost_cells; i++){
		for(int j = 0; j < nyCells + 2 * ghost_cells; j++){
			for(int k = 0; k < NUM_VARS; k++){
				ini_u_con[k] = u_con(i, j, k);
			}
			ini_con2pri(ini_u_pri, ini_u_con);
			for(int k = 0; k < NUM_VARS; k++){
				u_pri(i, j, k) = ini_u_pri[k];
			}
		}
	}
}

void MagneticFlux_x(double *Q, double *flux){
	double B2 = Q[5] * Q[5] + Q[6] * Q[6] + Q[7] * Q[7];
	double vB   = ( Q[1]*Q[5] + Q[2]*Q[6] + Q[3]*Q[7] ) / Q[0];
	flux[0] = 0;
	flux[1] = (2 - gamma_plasma) * B2 * 0.5 - Q[5] * Q[5];
	flux[2] = -Q[5] * Q[6];
	flux[3] = -Q[5] * Q[7];
	flux[4] = Q[1] / Q[0] * B2 * 0.5 - Q[5] * vB;
	flux[5] = 0;
	flux[6] = 0;
	flux[7] = 0;
}

void MagneticFlux_y(double *Q, double *flux){
	double B2 = Q[5] * Q[5] + Q[6] * Q[6] + Q[7] * Q[7];
	double vB   = ( Q[1]*Q[5] + Q[2]*Q[6] + Q[3]*Q[7] ) / Q[0];
	flux[0] = 0;
	flux[1] = -Q[5] * Q[6];
	flux[2] = (2 - gamma_plasma) * B2 * 0.5 - Q[6] * Q[6];
	flux[3] = -Q[6] * Q[7];
	flux[4] = Q[2] / Q[0] * B2 * 0.5 - Q[6] * vB;
	flux[5] = 0;
	flux[6] = 0;
	flux[7] = 0;
}

void MagneticFlux_Pressure_x(double *Q, double *flux){
	double B2 = Q[5] * Q[5] + Q[6] * Q[6] + Q[7] * Q[7];
	double vB = ( Q[1]*Q[5] + Q[2]*Q[6] + Q[3]*Q[7] ) / Q[0];
	flux[0] = 0;
	flux[1] = B2 * 0.5 - Q[5] * Q[5];
	flux[2] = -Q[5] * Q[6];
	flux[3] = -Q[5] * Q[7];
	flux[4] = Q[1] / Q[0] * B2 - Q[5] * vB;
	flux[5] = 0;
	flux[6] = 0;
	flux[7] = 0;
}

void MagneticFlux_Pressure_y(double *Q, double *flux){
	double B2 = Q[5] * Q[5] + Q[6] * Q[6] + Q[7] * Q[7];
	double vB   = ( Q[1]*Q[5] + Q[2]*Q[6] + Q[3]*Q[7] ) / Q[0];
	flux[0] = 0;
	flux[1] = -Q[5] * Q[6];
	flux[2] = B2 * 0.5 - Q[6] * Q[6];
	flux[3] = -Q[6] * Q[7];
	flux[4] = Q[2] / Q[0] * B2 - Q[6] * vB;
	flux[5] = 0;
	flux[6] = 0;
	flux[7] = 0;
}


void f_KineticEnergyFlux(double *Q, double *flux){
	flux[0] = 0;
	flux[1] = -(gamma_plasma - 1) * 0.5 * (Q[1] * Q[1] + Q[2] * Q[2] + Q[3] * Q[3]) / Q[0];
	flux[2] = 0;
	flux[3] = 0;
	flux[4] = 0;
	flux[5] = 0;
	flux[6] = 0;
	flux[7] = 0;
}

void g_KineticEnergyFlux(double *Q, double *flux){
	flux[0] = 0;
	flux[1] = 0;
	flux[2] = -(gamma_plasma - 1) * 0.5 * (Q[1] * Q[1] + Q[2] * Q[2] + Q[3] * Q[3]) / Q[0];
	flux[3] = 0;
	flux[4] = 0;
	flux[5] = 0;
	flux[6] = 0;
	flux[7] = 0;
}

double ini_calch_x(double *pri){
	double rho = pri[0];
	double p = pri[4];
	double Bx = pri[5];
	double By = pri[6];
	double Bz = pri[7];
	double B2 = Bx * Bx + By * By + Bz * Bz;
	
	double cf_x = eos_fast_speed(rho, p, Bx, B2);
	
	return fabs(pri[1]) + cf_x;
}

double ini_calch_y(double *pri){
	double rho = pri[0];
	double p = pri[4];
	double Bx = pri[5];
	double By = pri[6];
	double Bz = pri[7];
	double B2 = Bx * Bx + By * By + Bz * Bz;
	
	
	double cf_y = eos_fast_speed(rho, p, By, B2);
	
	return fabs(pri[2]) + cf_y;
}
