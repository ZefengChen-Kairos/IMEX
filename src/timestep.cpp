// Time-step control, conservative update, enthalpy/eigenvalue helpers.

#include "declarations.hpp"

double caldt(Grid &u_pri){
	vx_max = 0, vy_max = 0, v_max = 0;
	double dt;
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			if(fabs(u_pri(i,j,1)) > fabs(vx_max)){
				vx_max = fabs(u_pri(i,j,1));
			}
			if(fabs(u_pri(i,j,2)) > fabs(vy_max)){
				vy_max = fabs(u_pri(i,j,2));
			}
			if(sqrt(u_pri(i,j,1) * u_pri(i,j,1) + u_pri(i,j,2) * u_pri(i,j,2) + u_pri(i,j,3) * u_pri(i,j,3)) > fabs(v_max)){
				v_max = sqrt(u_pri(i,j,1) * u_pri(i,j,1) + u_pri(i,j,2) * u_pri(i,j,2) + u_pri(i,j,3) * u_pri(i,j,3));
			}
		}
	}
	
	//timestep == 1
	switch(caltimestep){
		case Fullydt:
			a_max_x = ch_x;
			a_max_y = ch_y;
			a_max = ch;
			break;
		case Alfvendt:
			if(timestep < full_step + 1){
				a_max_x = ch_x;
				a_max_y = ch_y;
				a_max = ch;
				//std::cout << "ch: "<< ch << endl;
			}else{
				a_max_x = cA_x;
				a_max_y = cA_y;
				a_max = cA;
			}
			break;
		case IMEXdt:
			if(timestep < full_step + 1){
				a_max_x = ch_x;
				a_max_y = ch_y;
				a_max = ch;
				//std::cout << "ch: "<< ch << endl;
			}else{
				a_max_x = vx_max;
				a_max_y = vy_max;
				a_max = v_max;
				//dt = CFL * 1 / (vx_max / dx);
			}
			break;
	}
	dt = CFL * 1 / (2 * a_max / fmin(dx,dy));
	fully_dt = CFL * 1 / (2 * ch / fmin(dx,dy));
	imex_dt = CFL * 1 / (2 * a_max / fmin(dx,dy));
	static double prev_dt = 0.0;
	if(timestep > 1 && prev_dt > 0.0) dt = fmin(dt, 2.0 * prev_dt);  // CASTOR enlarge_dt = 2.0
	prev_dt = dt;
	if(timestep > full_step){
		a_max_x = vx_max;
		a_max_y = vy_max;
		a_max = v_max;
	}
	return dt;
}

void update(Grid &u_con,Grid &flux, char dir){
	if(dir == 'x'){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
				for(int k = 0; k < NUM_VARS; k++){
					u_con(i, j, k) = u_con(i, j, k) - (adt/dx) * (flux(i, j, k) - flux(i - 1, j, k));
				}
			}
		}
	}else if(dir == 'y'){
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
				for(int k = 0; k < NUM_VARS; k++){
					u_con(i, j, k) = u_con(i, j, k) - (adt/dy) * (flux(i, j, k) - flux(i, j - 1, k));
				}
			}
		}
	}
}

void calch(Grid &u_pri){
	double ch_grid_x = 1e-8, ch_grid_y = 1e-8, ch_grid = 1e-8;
	for(int i = 0; i < nxCells + 2 * ghost_cells; i++){
		for(int j = 0; j < nyCells + 2 * ghost_cells; j++){
			double rho = u_pri(i,j,0);
			double p = u_pri(i,j,4);
			double Bx = u_pri(i,j,5);
			double By = u_pri(i,j,6);
			double Bz = u_pri(i,j,7);
			double B2 = Bx * Bx + By * By + Bz * Bz;
			double B = sqrt(B2);
			double v = sqrt(u_pri(i,j,1) * u_pri(i,j,1) + u_pri(i,j,2) * u_pri(i,j,2) + u_pri(i,j,3) * u_pri(i,j,3));
			double cf_x = eos_fast_speed(rho, p, Bx, B2);
			double cf_y = eos_fast_speed(rho, p, By, B2);
			
			double cf = eos_fast_speed(rho, p, B, B2);
			
			if(ch_grid_x < u_pri(i,j,1) + cf_x){
				ch_grid_x = u_pri(i,j,1) + cf_x;
			}
			if(ch_grid_y < u_pri(i,j,2) + cf_y){
				ch_grid_y = u_pri(i,j,2) + cf_y;
			}
			if(ch_grid < v + cf){
				ch_grid = v + cf;
			}
		}
	}
	ch_x = ch_grid_x;
	ch_y = ch_grid_y;
	ch = ch_grid;
}

void calca(Grid &u_pri){
	double ca_grid_x = 1e-8;
	double ca_grid_y = 1e-8;
	double ca_grid = 1e-8;
	double pri[NUM_VARS];
	for(int i = 0; i < nxCells + 2 * ghost_cells; i++){
		for(int j = 0; j < nyCells + 2 * ghost_cells; j++){
			
			for(int k = 0; k < NUM_VARS; k++){
				pri[k] = u_pri(i,j,k);
			}
			double B = sqrt(pri[5] * pri[5] + pri[6] * pri[6] + pri[7] * pri[7]);
			if(ca_grid_x < pri[5]  / sqrt(pri[0])){
				ca_grid_x = pri[5]  / sqrt(pri[0]);
			}
			if(ca_grid_y < pri[6]  / sqrt(pri[0])){
				ca_grid_y = pri[6]  / sqrt(pri[0]);
			}
			if(ca_grid < B / sqrt(pri[0])){
				ca_grid = B  / sqrt(pri[0]);
			}
		}
	}
	cA_x = ca_grid_x;
	cA_y = ca_grid_y;
	cA = ca_grid;
}

void calcp(Grid &/*u_pri*/){
	cp_x = sqrt(0.18 * ch_x);
	cp_y = sqrt(0.18 * ch_y);
}

void store(Grid &u_con, Array& vxn,Array& vyn, Array& Bxn, Array& Byn){
	for(int i = 0; i < nxCells + 2 * ghost_cells; i++){
		for(int j = 0; j < nyCells + 2 * ghost_cells; j++){
			vxn(i,j) = u_con(i,j,1) / u_con(i,j,0);
			vyn(i,j) = u_con(i,j,2) / u_con(i,j,0);
			Bxn(i,j) = u_con(i,j,5);
			Byn(i,j) = u_con(i,j,6);
		}
	}
}

void bound_cond_array(Array &u){
	switch (x_bc) {
		case x_Transmissive:case x_Dirichlet:
			for(int j = 0; j < nyCells + ghost_cells * 2; j++){
				for(int k = 0; k < ghost_cells; k++){
					u(k,j) = u(ghost_cells,j);
					u(nxCells + ghost_cells + k,j) = u(nxCells + ghost_cells - 1,j);
				}
			}
			break;
		case x_Reflective:
			for(int j = 0; j < nyCells + ghost_cells * 2; j++){
				for(int ghost = 0; ghost < ghost_cells; ghost++){
					u(ghost,j) = -u(2 * ghost_cells - 1 - ghost,j);
					u(nxCells + ghost_cells + ghost,j) = -u(nxCells + ghost_cells - 1 - ghost,j);
				}
			}
			break;
		case x_Periodic:
			for(int j = 0; j < nyCells + ghost_cells * 2; j++){
				for(int ghost = 0; ghost < ghost_cells; ghost++){
					u(ghost,j) = u(nxCells + ghost,j);
					u(nxCells + ghost_cells + ghost,j) = u(ghost_cells + ghost,j);
				}
			}
			break;
	}
	
	switch(y_bc){
		case y_Transmissive: case y_Dirichlet: case y_JetInflow:
			// jet: 通用数组(shock_eta/Ez 等)透射即可;需要 strip Dirichlet 的量
			// (rhov*_star / final_pressure)在调用点单独覆盖
			for(int i = 0; i < nxCells + ghost_cells * 2; i++){
				for(int ghost = 0; ghost < ghost_cells; ghost++){
					u(i,ghost) = u(i,ghost_cells);
					u(i,nyCells + ghost_cells + ghost) = u(i,nyCells + ghost_cells - 1);
				}
			}
			break;
		case y_Reflective:
			for(int i = 0; i < nxCells + ghost_cells * 2; i++){
				for(int ghost = 0; ghost < ghost_cells; ghost++){
					u(i,ghost) = -u(i,2 * ghost_cells - 1 - ghost);
					u(i,nyCells + ghost_cells + ghost) = -u(i,nyCells + ghost_cells - 1 - ghost);
				}
			}
			break;
		case y_Periodic:
			for(int i = 0; i < nxCells + ghost_cells * 2; i++){
				for(int ghost = 0; ghost < ghost_cells; ghost++){
					u(i,ghost) = u(i,nyCells + ghost);
					u(i,nyCells + ghost_cells + ghost) = u(i,ghost_cells + ghost);
				}
			}
			break;
	}
}

double calenthalpy(double *con){
	double e = calecon(con);
	double pri[NUM_VARS];
	ini_con2pri(pri, con);
	double h = e + pri[4] / pri[0];
	return h;
}

double caltotalenthalpy(double *con){
	double pri[NUM_VARS];
	ini_con2pri(pri, con);
	double h = (con[4] + pri[4]) / pri[0];
	return h;
}

double calkinetic(double *con){
	double rho = con[0];
	double vx = con[1] / con[0];
	double vy = con[2] / con[0];
	double vz = con[3] / con[0];
	double kinetic_energy = 0.5 * rho*(vx * vx + vy * vy + vz * vz);
	//double kinetic_energy = 0.5 * (con[1] * con[1] + con[2] * con[2] + con[3] * con[3]) / (con[0]);
	return kinetic_energy;
}

double calmagenergy(double *con){
	double mag_energy = 0.5 * (con[5] * con[5] + con[6] * con[6] + con[7] * con[7]);
	return mag_energy;
}

void MagneticEigenvalues(Grid &u_pri){
	double bmax_x = 0;
	double bmax_y = 0;
	double b_x,b_y;
	for(int i = 0 ; i < nxCells + 2 * ghost_cells; i++){
		for(int j = 0; j < nyCells + 2 * ghost_cells; j++){
			b_x = 0.5 * (fabs(u_pri(i,j,1)) + sqrt(u_pri(i,j,1) * u_pri(i,j,1) + 4 * (u_pri(i,j,5) * u_pri(i,j,5) + u_pri(i,j,6) * u_pri(i,j,6) + u_pri(i,j,7) * u_pri(i,j,7)) / u_pri(i,j,0)  ));
			
			b_y = 0.5 * (fabs(u_pri(i,j,2)) + sqrt(u_pri(i,j,2) * u_pri(i,j,2) + 4 * (u_pri(i,j,5) * u_pri(i,j,5) + u_pri(i,j,6) * u_pri(i,j,6) + u_pri(i,j,7) * u_pri(i,j,7)) / u_pri(i,j,0)  ));
			if(fabs(bmax_x) < fabs(b_x)){
				bmax_x = b_x;
			}
			if(fabs(bmax_y) < fabs(b_y)){
				bmax_y = b_y;
			}
		}
	}
	B_max_x = bmax_x;
	B_max_y = bmax_y;
}

void EnergyEigenvalues(Grid &u_pri){
	double umax_x = 0, umax_y = 0;
	double u_x, u_y;
	for(int i = 0 ; i < nxCells + 2 * ghost_cells; i++){
		for(int j = 0; j < nyCells + 2 * ghost_cells; j++){
			double c = eos_sound_a(u_pri(i,j,0), u_pri(i,j,4));
			double m = 0.5 * (u_pri(i,j,5) * u_pri(i,j,5) + u_pri(i,j,6) * u_pri(i,j,6) + u_pri(i,j,7) * u_pri(i,j,7));
			double k = 0.5 * (u_pri(i,j,1) * u_pri(i,j,1) + u_pri(i,j,2) * u_pri(i,j,2) + u_pri(i,j,3) * u_pri(i,j,3));
			u_x = 0.5 * (u_pri(i,j,1) + sqrt(u_pri(i,j,1) * u_pri(i,j,1) + 4 * (c * c - (gamma_plasma - 1) * (m + k + u_pri(i,j,1) * u_pri(i,j,1)) )));
			u_y = 0.5 * (u_pri(i,j,2) + sqrt(u_pri(i,j,2) * u_pri(i,j,2) + 4 * (c * c - (gamma_plasma - 1) * (m + k + u_pri(i,j,2) * u_pri(i,j,2)) )));
			if(umax_x < u_x){
				umax_x = u_x;
			}
			if(umax_y < u_y){
				umax_y = u_y;
			}
		}
	}
	U_max_x = umax_x;
	U_max_y = umax_y;
}
