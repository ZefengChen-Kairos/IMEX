// Shock detector and the explicit numerical fluxes (energy, pressure, Rusanov).

#include "declarations.hpp"

void ShockDetector(Array &shock_eta, Grid &u_con) {
	// 转换守恒变量到原始变量 (rho, u, v, w, p, Bx, By, Bz)
	Grid u_pri(nxCells, nyCells);
	ucon2upri(u_pri, u_con);
	
	// 中间数组:先存 chi_tilde (公式 4.40),再做邻居扩展 (公式 4.41)
	Array shock_eta_tilde(nxCells, nyCells);
	
	// ---- 第一步:计算 chi_tilde (公式 4.40) ----
	for (int i = ghost_cells; i < nxCells + ghost_cells; i++) {
		for (int j = ghost_cells; j < nyCells + ghost_cells; j++) {
			
			// 3x3 邻域内的快磁声速 c_f = sqrt(a^2 + c_A^2),取最小值
			//   a   = sqrt(gamma p / rho)        (声速)
			//   c_A = |B| / sqrt(rho)            (Alfven 速)
			//   u_pri: 0=rho, 4=p, 5=Bx, 6=By, 7=Bz
			double cf_min_nbr = 1.0e300;
			for (int di = -1; di <= 1; di++) {
				for (int dj = -1; dj <= 1; dj++) {
					double rho_nbr = u_pri(i + di, j + dj, 0);
					double p_nbr   = u_pri(i + di, j + dj, 4);
					double Bx_nbr  = u_pri(i + di, j + dj, 5);
					double By_nbr  = u_pri(i + di, j + dj, 6);
					double Bz_nbr  = u_pri(i + di, j + dj, 7);
					
					double a2  = eos_sound_a2(rho_nbr, p_nbr);
					double cA2 = (Bx_nbr * Bx_nbr + By_nbr * By_nbr + Bz_nbr * Bz_nbr)
					/ rho_nbr;
					double cf  = std::sqrt(a2 + cA2);
					
					if (cf < cf_min_nbr) cf_min_nbr = cf;
				}
			}
			
			// 速度散度 div.u = du/dx + dv/dy (中心差分; w/vz 不进 2D 散度)
			//   u_pri(.,.,1)=u, u_pri(.,.,2)=v
			double div_u = (u_pri(i + 1, j, 1) - u_pri(i - 1, j, 1)) / (2.0 * dx)
			+ (u_pri(i, j + 1, 2) - u_pri(i, j - 1, 2)) / (2.0 * dy);
			
			// 公式 (4.40):
			//   chi_tilde = min[1, max(0, -(h*div.u + k1*cf_min) / (k1*cf_min))]
			double h           = std::min(dx, dy);
			double numerator   = h * div_u + shock_kappa1 * cf_min_nbr;
			double denominator = shock_kappa1 * cf_min_nbr;
			
			shock_eta_tilde(i, j) = std::min(1.0,
				std::max(0.0, -numerator / denominator));
		}
	}
	
	// 边界条件 (让邻域扩展可以访问 ghost cell)
	bound_cond_array(shock_eta_tilde);
	
	// ---- 第二步:邻居扩展 (公式 4.41) ----
	//   chi_{i,j} = max over 3x3 of chi_tilde
	for (int i = ghost_cells; i < nxCells + ghost_cells; i++) {
		for (int j = ghost_cells; j < nyCells + ghost_cells; j++) {
			double chi_max = 0.0;
			for (int di = -1; di <= 1; di++) {
				for (int dj = -1; dj <= 1; dj++) {
					if (shock_eta_tilde(i + di, j + dj) > chi_max) {
						chi_max = shock_eta_tilde(i + di, j + dj);
					}
				}
			}
			shock_eta(i, j) = chi_max;
//			shock_eta(i, j) = 1.0;
//			shock_eta(i, j) = 0.0;
		}
	}
	bound_cond_array(shock_eta);
}

void Energy_flux_x(double *QL, double *QR, double *flux){
	double fCL[NUM_VARS], fCR[NUM_VARS];
	
	f_convective(fCL, QL);
	f_convective(fCR, QR);
	fCL[4] = 0;
	fCR[4] = 0;
	double fKL[NUM_VARS], fKR[NUM_VARS];
	
	f_KineticEnergyFlux(QL, fKL);
	f_KineticEnergyFlux(QR, fKR);
	
	double fBL[NUM_VARS], fBR[NUM_VARS];
	MagneticFlux_x(QL, fBL);
	MagneticFlux_x(QR, fBR);
	
	double fL[NUM_VARS], fR[NUM_VARS];
	for(int k = 0; k < NUM_VARS; k++){
		fL[k] = fCL[k] + fBL[k] + fKL[k];
		fR[k] = fCR[k] + fBR[k] + fKR[k];
	}
	
	
	
	double priL[NUM_VARS], priR[NUM_VARS];
	ini_con2pri(priL, QL);
	ini_con2pri(priR, QR);
	
	double s_max = 0;
	
	if(GLF){
		s_max = a_max;
	}else{
		// 通量耗散始终用材料速度(对齐 CASTOR RusanovMHD_B);磁声速 startup 只用于 dt(caldt),不进通量
		s_max = fmax(fabs(priL[1]), fabs(priR[1]));
	}
	
	for(int k = 0; k < NUM_VARS; k++){
		flux[k] = 0.5 * (fL[k] + fR[k]) - 0.5 * s_max * (QR[k] - QL[k]);
	}
}

void Energy_flux_y(double *QL, double *QR, double *flux){
	double fCL[NUM_VARS], fCR[NUM_VARS];
	
	g_convective(fCL, QL);
	g_convective(fCR, QR);
	fCL[4] = 0;
	fCR[4] = 0;
	double fKL[NUM_VARS], fKR[NUM_VARS];
	
	g_KineticEnergyFlux(QL, fKL);
	g_KineticEnergyFlux(QR, fKR);
	
	double fBL[NUM_VARS], fBR[NUM_VARS];
	MagneticFlux_y(QL, fBL);
	MagneticFlux_y(QR, fBR);
	
	double fL[NUM_VARS], fR[NUM_VARS];
	for(int k = 0; k < NUM_VARS; k++){
		fL[k] = fCL[k] + fBL[k] + fKL[k];
		fR[k] = fCR[k] + fBR[k] + fKR[k];
	}
	
	
	double priL[NUM_VARS], priR[NUM_VARS];
	ini_con2pri(priL, QL);
	ini_con2pri(priR, QR);
	
	double s_max = 0;
	
	
	
	if(GLF){
		s_max = a_max;
	}else{
		// 通量耗散始终用材料速度(对齐 CASTOR RusanovMHD_B);磁声速 startup 只用于 dt(caldt),不进通量
		s_max = fmax(fabs(priL[2]), fabs(priR[2]));
	}
	
	for(int k = 0; k < NUM_VARS; k++){
		flux[k] = 0.5 * (fL[k] + fR[k]) - 0.5 * s_max * (QR[k] - QL[k]);
	}
}

void Energy_flux_Update(Grid &u_con, Grid &Qs, Grid &w_xL, Grid &w_xR, Grid &w_yL, Grid &w_yR,double h){
	for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			double Q_xL[NUM_VARS], Q_xR[NUM_VARS], fp[NUM_VARS], fm[NUM_VARS];
			double Q_yL[NUM_VARS], Q_yR[NUM_VARS], gp[NUM_VARS], gm[NUM_VARS];
			for(int k = 0; k < NUM_VARS; k++){
				Q_xL[k] = w_xR(i,j,k);
				Q_xR[k] = w_xL(i + 1,j,k);
			}
			Energy_flux_x(Q_xL, Q_xR, fp);
			
			
			for(int k = 0; k < NUM_VARS; k++){
				Q_yL[k] = w_yR(i,j,k);
				Q_yR[k] = w_yL(i,j + 1,k);
				
			}
			Energy_flux_y(Q_yL, Q_yR, gp);
			
			
			for(int k = 0; k < NUM_VARS; k++){
				Q_xL[k] = w_xR(i - 1,j,k);
				Q_xR[k] = w_xL(i,j,k);
			}
			Energy_flux_x(Q_xL, Q_xR, fm);
			
			for(int k = 0; k < NUM_VARS; k++){
				Q_yL[k] = w_yR(i,j - 1,k);
				Q_yR[k] = w_yL(i,j,k);
			}
			Energy_flux_y(Q_yL, Q_yR, gm);
			
			for(int k = 0; k < 5; k++){
				u_con(i,j,k) = Qs(i,j,k) - h / dx * (fp[k] - fm[k]) - h / dy * (gp[k] - gm[k]);
			}
		}
	}
}

void Pressure_flux_x(double *QL, double *QR, double *flux){
	double fCL[NUM_VARS], fCR[NUM_VARS];
	f_convective(fCL, QL);
	f_convective(fCR, QR);
	
	double fBL[NUM_VARS], fBR[NUM_VARS];
	MagneticFlux_Pressure_x(QL, fBL);
	MagneticFlux_Pressure_x(QR, fBR);
	
	double fL[NUM_VARS], fR[NUM_VARS];
	for(int k = 0; k < NUM_VARS; k++){
		fL[k] = fCL[k] + fBL[k];
		fR[k] = fCR[k] + fBR[k];
	}
	
	double priL[NUM_VARS], priR[NUM_VARS];
	ini_con2pri(priL, QL);
	ini_con2pri(priR, QR);
	
	double s_max = 0;
	
	if(GLF){
		s_max = a_max;
	}else{
		// 通量耗散始终用材料速度(对齐 CASTOR RusanovMHD_B);磁声速 startup 只用于 dt(caldt),不进通量
		s_max = fmax(fabs(priL[1]), fabs(priR[1]));
	}
	
	
	for(int k = 0; k < NUM_VARS; k++){
		if(k != 4){
			flux[k] = 0.5 * (fL[k] + fR[k]) - 0.5 * s_max * (QR[k] - QL[k]);
		}else{
			double KL = 0.5 * (QL[1] * QL[1] + QL[2] * QL[2] + QL[3] * QL[3]) / QL[0];
			double KR = 0.5 * (QR[1] * QR[1] + QR[2] * QR[2] + QR[3] * QR[3]) / QR[0];
			flux[k] = 0.5 * (fL[k] + fR[k]) - 0.5 * s_max * (KR - KL);
			flux[k] = 0.5 * (fL[k] + fR[k]) - 0.5 * s_max * (QR[k] - QL[k]);
		}
		
//		flux[k] = 0.5 * (fL[k] + fR[k]) - 0.5 * s_max * (QR[k] - QL[k]);
	}
}

void Pressure_flux_y(double *QL, double *QR, double *flux){
	double fCL[NUM_VARS], fCR[NUM_VARS];
	
	g_convective(fCL, QL);
	g_convective(fCR, QR);
	
	double fBL[NUM_VARS], fBR[NUM_VARS];
	MagneticFlux_Pressure_y(QL, fBL);
	MagneticFlux_Pressure_y(QR, fBR);
	
	double fL[NUM_VARS], fR[NUM_VARS];
	for(int k = 0; k < NUM_VARS; k++){
		fL[k] = fCL[k] + fBL[k];
		fR[k] = fCR[k] + fBR[k];
	}
	
	
	double priL[NUM_VARS], priR[NUM_VARS];
	ini_con2pri(priL, QL);
	ini_con2pri(priR, QR);
	
	double s_max = 0;
	
	if(GLF){
		s_max = a_max;
	}else{
		// 通量耗散始终用材料速度(对齐 CASTOR RusanovMHD_B);磁声速 startup 只用于 dt(caldt),不进通量
		s_max = fmax(fabs(priL[2]), fabs(priR[2]));
	}
	
	for(int k = 0; k < NUM_VARS; k++){
		if(k != 4){
			flux[k] = 0.5 * (fL[k] + fR[k]) - 0.5 * s_max * (QR[k] - QL[k]);
		}else{
			double KL = 0.5 * (QL[1] * QL[1] + QL[2] * QL[2] + QL[3] * QL[3]) / QL[0];
			double KR = 0.5 * (QR[1] * QR[1] + QR[2] * QR[2] + QR[3] * QR[3]) / QR[0];
			flux[k] = 0.5 * (fL[k] + fR[k]) - 0.5 * s_max * (KR - KL);
			flux[k] = 0.5 * (fL[k] + fR[k]) - 0.5 * s_max * (QR[k] - QL[k]);
		}
		
//		flux[k] = 0.5 * (fL[k] + fR[k]) - 0.5 * s_max * (QR[k] - QL[k]);
	}
}


void Pressure_flux_Update(Grid &u_con, Grid &Qs, Grid &w_xL, Grid &w_xR, Grid &w_yL, Grid &w_yR,double h){
	for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			double Q_xL[NUM_VARS], Q_xR[NUM_VARS], fp[NUM_VARS], fm[NUM_VARS];
			double Q_yL[NUM_VARS], Q_yR[NUM_VARS], gp[NUM_VARS], gm[NUM_VARS];
			for(int k = 0; k < NUM_VARS; k++){
				Q_xL[k] = w_xR(i,j,k);
				Q_xR[k] = w_xL(i + 1,j,k);
			}
			Pressure_flux_x(Q_xL, Q_xR, fp);
			
			
			for(int k = 0; k < NUM_VARS; k++){
				Q_yL[k] = w_yR(i,j,k);
				Q_yR[k] = w_yL(i,j + 1,k);
				
			}
			Pressure_flux_y(Q_yL, Q_yR, gp);
			
			for(int k = 0; k < NUM_VARS; k++){
				Q_xL[k] = w_xR(i - 1,j,k);
				Q_xR[k] = w_xL(i,j,k);
			}
			Pressure_flux_x(Q_xL, Q_xR, fm);
			
			for(int k = 0; k < NUM_VARS; k++){
				Q_yL[k] = w_yR(i,j - 1,k);
				Q_yR[k] = w_yL(i,j,k);
			}
			Pressure_flux_y(Q_yL, Q_yR, gm);
			
			for(int k = 0; k < 5; k++){
				u_con(i,j,k) = Qs(i,j,k) - h / dx * (fp[k] - fm[k]) - h / dy * (gp[k] - gm[k]);
			}
		}
	}
}


void RusanovC_x(double *QL, double *QR, double *flux){
	double fL[NUM_VARS], fR[NUM_VARS];
	
	f_convective(fL, QL);
	f_convective(fR, QR);
	
	switch(formula){
		case EnergyFormula:
			fL[4] = 0;
			fR[4] = 0;
			break;
		case PressureFormula:
			break;
	}
	
	double priL[NUM_VARS], priR[NUM_VARS];
	ini_con2pri(priL, QL);
	ini_con2pri(priR, QR);
	
	double s_max = 0;
	
	if(GLF){
		s_max = a_max;
	}else{
		// 通量耗散始终用材料速度(对齐 CASTOR RusanovMHD_B);磁声速 startup 只用于 dt(caldt),不进通量
		s_max = fmax(fabs(priL[1]), fabs(priR[1]));
	}
	
	for(int k = 0; k < NUM_VARS; k++){
		if(k != 4){
			flux[k] = 0.5 * (fL[k] + fR[k]) - 0.5 * s_max * (QR[k] - QL[k]);
		}else{
			double KL = 0.5 * (QL[1] * QL[1] + QL[2] * QL[2] + QL[3] * QL[3]) / QL[0];
			double KR = 0.5 * (QR[1] * QR[1] + QR[2] * QR[2] + QR[3] * QR[3]) / QR[0];
			flux[k] = 0.5 * (fL[k] + fR[k]) - 0.5 * s_max * (KR - KL);
			flux[k] = 0.5 * (fL[k] + fR[k]) - 0.5 * s_max * (QR[k] - QL[k]);
		}
	}
}

void RusanovC_y(double *QL, double *QR, double *flux){
	double gL[NUM_VARS], gR[NUM_VARS];
	g_convective(gL, QL);
	g_convective(gR, QR);
	
	switch(formula){
		case EnergyFormula:
			gL[4] = 0;
			gR[4] = 0;
			break;
		case PressureFormula:
			break;
	}
	
	double priL[NUM_VARS], priR[NUM_VARS];
	ini_con2pri(priL, QL);
	ini_con2pri(priR, QR);
	
	double s_max = 0;
	
	if(GLF){
		s_max = a_max;
	}else{
		// 通量耗散始终用材料速度(对齐 CASTOR RusanovMHD_B);磁声速 startup 只用于 dt(caldt),不进通量
		s_max = fmax(fabs(priL[2]), fabs(priR[2]));
	}
	
	for(int k = 0; k < NUM_VARS; k++){
		if(k != 4){
			flux[k] = 0.5 * (gL[k] + gR[k]) - 0.5 * s_max * (QR[k] - QL[k]);
		}else{
			double KL = 0.5 * (QL[1] * QL[1] + QL[2] * QL[2] + QL[3] * QL[3]) / QL[0];
			double KR = 0.5 * (QR[1] * QR[1] + QR[2] * QR[2] + QR[3] * QR[3]) / QR[0];
			flux[k] = 0.5 * (gL[k] + gR[k]) - 0.5 * s_max * (KR - KL);
			flux[k] = 0.5 * (gL[k] + gR[k]) - 0.5 * s_max * (QR[k] - QL[k]);
		}
	}
}


void RusanovC_Update(Grid &u_con, Grid &Qs, Grid &w_xL, Grid &w_xR, Grid &w_yL, Grid &w_yR,double h){
	for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			double Q_xL[NUM_VARS], Q_xR[NUM_VARS], fp[NUM_VARS], fm[NUM_VARS];
			double Q_yL[NUM_VARS], Q_yR[NUM_VARS], gp[NUM_VARS], gm[NUM_VARS];
			for(int k = 0; k < NUM_VARS; k++){
				Q_xL[k] = w_xR(i,j,k);
				Q_xR[k] = w_xL(i + 1,j,k);
				
			}
			RusanovC_x(Q_xL, Q_xR, fp);
			
			for(int k = 0; k < NUM_VARS; k++){
				Q_yL[k] = w_yR(i,j,k);
				Q_yR[k] = w_yL(i,j + 1,k);
				
			}
			RusanovC_y(Q_yL, Q_yR, gp);
			
			
			for(int k = 0; k < NUM_VARS; k++){
				Q_xL[k] = w_xR(i - 1,j,k);
				Q_xR[k] = w_xL(i,j,k);
			}
			RusanovC_x(Q_xL, Q_xR, fm);
			
			for(int k = 0; k < NUM_VARS; k++){
				Q_yL[k] = w_yR(i,j - 1,k);
				Q_yR[k] = w_yL(i,j,k);
			}
			RusanovC_y(Q_yL, Q_yR, gm);
			
			switch(formula){
				case EnergyFormula:
					for(int k = 0; k < 5; k++){
						u_con(i,j,k) = Qs(i,j,k) - h / dx * (fp[k] - fm[k]) - h / dy * (gp[k] - gm[k]);
					}
					break;
				case PressureFormula:
					for(int k = 0; k < 5; k++){
						u_con(i,j,k) = Qs(i,j,k) - h / dx * (fp[k] - fm[k]) - h / dy * (gp[k] - gm[k]);
					}
					break;
			}
		}
	}
}
