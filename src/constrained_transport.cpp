// Staggered and unstaggered CT, vector-potential initialisation, EMF diagnostics.

#include "declarations.hpp"

// CT-consistent BC for a staggered NORMAL-face field.
//   axis=0 -> Bx_p0 (face at i+1/2, normal = x);  axis=1 -> By_0p (face at j+1/2, normal = y).
// 法向(normal):边界法向面 index ghost_cells-1 / N+ghost-1 是 CT 的真实 DOF,必须保留;
//               只把外侧 ghost 面由边界面零阶外推(否则破坏边界 cell 的 div-free)。
// 切向(tangential):普通 transmissive 复制。Periodic 方向:周期复制。
// 通用 bound_cond_array(cell-centered transmissive)会把边界法向面也冲掉 -> 边界 divB 指数增长。
void bound_cond_staggered(Array &u, int axis){
	// ---------- x-direction ----------
	if(x_bc == x_Periodic){
		for(int j = 0; j < nyCells + ghost_cells*2; j++)
			for(int g = 0; g < ghost_cells; g++){
				u(g,j) = u(nxCells+g,j);
				u(nxCells+ghost_cells+g,j) = u(ghost_cells+g,j);
			}
	}else if(axis==0){ // x is NORMAL (Bx_p0): preserve boundary face ghost_cells-1
		for(int j = 0; j < nyCells + ghost_cells*2; j++){
			for(int k = 0; k < ghost_cells-1; k++) u(k,j) = u(ghost_cells-1,j);
			for(int k = 0; k < ghost_cells;   k++) u(nxCells+ghost_cells+k,j) = u(nxCells+ghost_cells-1,j);
		}
	}else{ // x is TANGENTIAL (By_0p): plain transmissive
		for(int j = 0; j < nyCells + ghost_cells*2; j++)
			for(int k = 0; k < ghost_cells; k++){
				u(k,j) = u(ghost_cells,j);
				u(nxCells+ghost_cells+k,j) = u(nxCells+ghost_cells-1,j);
			}
	}
	// ---------- y-direction ----------
	if(y_bc == y_Periodic){
		for(int i = 0; i < nxCells + ghost_cells*2; i++)
			for(int g = 0; g < ghost_cells; g++){
				u(i,g) = u(i,nyCells+g);
				u(i,nyCells+ghost_cells+g) = u(i,ghost_cells+g);
			}
	}else if(axis==1){ // y is NORMAL (By_0p): preserve boundary face ghost_cells-1
		for(int i = 0; i < nxCells + ghost_cells*2; i++){
			for(int k = 0; k < ghost_cells-1; k++) u(i,k) = u(i,ghost_cells-1);
			for(int k = 0; k < ghost_cells;   k++) u(i,nyCells+ghost_cells+k) = u(i,nyCells+ghost_cells-1);
		}
	}else{ // y is TANGENTIAL (Bx_p0): plain transmissive
		for(int i = 0; i < nxCells + ghost_cells*2; i++)
			for(int k = 0; k < ghost_cells; k++){
				u(i,k) = u(i,ghost_cells);
				u(i,nyCells+ghost_cells+k) = u(i,nyCells+ghost_cells-1);
			}
	}
	// ---------- MHD jet: 底边 inflow strip 覆盖 ----------
	// 注意: 不碰边界法向面 j=ghost_cells-1 本身(由 CT 演化;strip 内 ghost 态给对后
	// Ez = vx*By - vy*Bx = 0,该面自动保持 B0,硬钉反而向第一排内点注入 divB)。
	if(y_bc == y_JetInflow){
		if(axis == 1){          // By_0p: 边界面以下的 ghost 面 = jet_B0
			for(int i = 0; i < nxCells + ghost_cells*2; i++)
				if(jet_inlet(i))
					for(int k = 0; k < ghost_cells - 1; k++) u(i,k) = jet_B0;
		}else{                  // Bx_p0(底边切向): ghost 行 = 0(jet Bx=0),按面心 x 判 strip
			for(int i = 0; i < nxCells + ghost_cells*2 - 1; i++)
				if(jet_inlet_face(i))
					for(int k = 0; k < ghost_cells; k++) u(i,k) = 0.0;
		}
	}
}

void Staggered_CT(Grid& u_con, Array &Bx_p0, Array &By_0p, Array &Bx_p0_Qs, Array &By_0p_Qs,double h){
	Array Ez_0p(nxCells,nyCells), Ez_p0(nxCells,nyCells), Ez_bar_pp(nxCells,nyCells),Ez_pp(nxCells,nyCells),Ez_tilde(nxCells,nyCells);
	double Q_L[NUM_VARS], Q_R[NUM_VARS], fL[NUM_VARS], fR[NUM_VARS], gL[NUM_VARS], gR[NUM_VARS];
	Grid w_xL(nxCells,nyCells), w_xR(nxCells,nyCells),w_yL(nxCells,nyCells), w_yR(nxCells,nyCells);
	
	TVDReconstruction(u_con,w_xL,w_xR,w_yL,w_yR);
	
	for(int i = 0; i < nxCells + 2 * ghost_cells - 1; i++){
		for(int j = 0; j < nyCells + 2 * ghost_cells - 1; j++){
			for(int k = 0; k < NUM_VARS; k++){
				Q_L[k] = w_xR(i,j,k);
				Q_R[k] = w_xL(i + 1,j,k);
			}
			
			f_full_con(fL,Q_L);
			f_full_con(fR,Q_R);
			
			Ez_p0(i,j) = -(0.5 * (fL[6] + fR[6]) - 0.5 * fmax(fabs(Q_L[1]/Q_L[0]),fabs(Q_R[1]/Q_R[0])) * (Q_R[6] - Q_L[6]));
			Ez_p0(i,j) = -(0.5 * (fL[6] + fR[6]));

			for(int k = 0; k < NUM_VARS; k++){
				Q_L[k] = w_yR(i,j,k);
				Q_R[k] = w_yL(i,j + 1,k);
			}
			
			g_full_con(gL, Q_L);
			g_full_con(gR, Q_R);
		
			Ez_0p(i,j) = (0.5 * (gL[5] + gR[5]) - 0.5 * fmax(fabs(Q_L[2]/Q_L[0]),fabs(Q_R[2]/Q_R[0])) * (Q_R[5] - Q_L[5]));
			Ez_0p(i,j) = (0.5 * (gL[5] + gR[5]));
		}
	}
	
	bound_cond_array(Ez_p0);
	bound_cond_array(Ez_0p);
	
	for(int j = 0; j < nyCells + 2 * ghost_cells - 2; j++){
		for(int i = 0; i < nxCells + 2 * ghost_cells - 2; i++){
			Ez_bar_pp(i,j) = 0.25 * (Ez_p0(i,j) + Ez_p0(i,j + 1) + Ez_0p(i,j) + Ez_0p(i + 1,j));
		}
	}
	
	for(int j = 0; j < nyCells + 2 * ghost_cells; j++){
		for(int i = 0; i < nxCells + 2 * ghost_cells; i++){
			Ez_tilde(i,j) = (u_con(i,j,2) / u_con(i,j,0) * u_con(i,j,5) - u_con(i,j,1) / u_con(i,j,0) * u_con(i,j,6));
		}
	}
	bound_cond_array(Ez_bar_pp);
	
	
	for(int j = 1; j < nyCells + 2 * ghost_cells - 2; j++){
		for(int i = 1; i < nxCells + 2 * ghost_cells - 2; i++){
			double term_1 = 0.5 * ((Ez_p0(i,j) - Ez_tilde(i,j)) / (dx / 2.0) + (Ez_p0(i,j + 1) - Ez_tilde(i,j + 1)) / (dx / 2.0));
			double term_2 = 0.5 * ((Ez_tilde(i + 1,j) - Ez_p0(i,j)) / (dx / 2.0) + (Ez_tilde(i + 1,j + 1) - Ez_p0(i,j + 1)) / (dx / 2.0));
			double term_3 = 0.5 * ((Ez_0p(i,j) - Ez_tilde(i,j)) / (dy / 2.0) + (Ez_0p(i + 1,j) - Ez_tilde(i + 1,j)) / (dy / 2.0));
			double term_4 = 0.5 * ((Ez_tilde(i,j + 1) - Ez_0p(i,j)) / (dy / 2.0) + (Ez_tilde(i + 1,j + 1) - Ez_0p(i + 1,j)) / (dy / 2.0));
			
			Ez_pp(i,j) = Ez_bar_pp(i,j) + dx / 8.0 * (term_1 - term_2) + dy / 8.0 * (term_3 - term_4);
		}
	}
	
	//bound_cond_array(Ez_pp);
	for(int j = ghost_cells - 1; j < nyCells + ghost_cells; j++){
		for(int i = ghost_cells - 1; i < nxCells + ghost_cells; i++){
			Bx_p0(i,j) = Bx_p0_Qs(i,j) - h / dy * (Ez_pp(i,j) - Ez_pp(i,j - 1));
			By_0p(i,j) = By_0p_Qs(i,j) + h / dx * (Ez_pp(i,j) - Ez_pp(i - 1,j));
		}
	}
	bound_cond_staggered(Bx_p0, 0);   // x-normal: 保留边界法向面(CT div-free)
	bound_cond_staggered(By_0p, 1);   // y-normal
	for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			u_con(i,j,5) = 0.5 * (Bx_p0(i - 1,j) + Bx_p0(i,j));
			u_con(i,j,6) = 0.5 * (By_0p(i,j - 1) + By_0p(i,j));
		}
	}
	bound_cond_u(u_con);
}

// ===== CT EMF 静态诊断:在精确初值上把角点 Ez_pp 分解为三个来源 =====
// 角点组装是 (Ez_p0, Ez_0p, Ez_tilde) 的线性算子,所以可以分三次调用逐项归因
void CT_CornerAssemble(Array &Ez_pp, Array &Ez_p0, Array &Ez_0p, Array &Ez_tilde){
	Array Ez_bar(nxCells, nyCells);
	for(int i = 0; i < nxCells + 2*ghost_cells; i++)
		for(int j = 0; j < nyCells + 2*ghost_cells; j++) Ez_bar(i,j) = 0.0;
	for(int j = 0; j < nyCells + 2*ghost_cells - 2; j++)
		for(int i = 0; i < nxCells + 2*ghost_cells - 2; i++)
		Ez_bar(i,j) = 0.25*(Ez_p0(i,j)+Ez_p0(i,j+1)+Ez_0p(i,j)+Ez_0p(i+1,j));
	bound_cond_array(Ez_bar);
	for(int j = 1; j < nyCells + 2*ghost_cells - 2; j++){
		for(int i = 1; i < nxCells + 2*ghost_cells - 2; i++){
			double t1 = 0.5*((Ez_p0(i,j)-Ez_tilde(i,j))/(dx/2.0)+(Ez_p0(i,j+1)-Ez_tilde(i,j+1))/(dx/2.0));
			double t2 = 0.5*((Ez_tilde(i+1,j)-Ez_p0(i,j))/(dx/2.0)+(Ez_tilde(i+1,j+1)-Ez_p0(i,j+1))/(dx/2.0));
			double t3 = 0.5*((Ez_0p(i,j)-Ez_tilde(i,j))/(dy/2.0)+(Ez_0p(i+1,j)-Ez_tilde(i+1,j))/(dy/2.0));
			double t4 = 0.5*((Ez_tilde(i,j+1)-Ez_0p(i,j))/(dy/2.0)+(Ez_tilde(i+1,j+1)-Ez_0p(i+1,j))/(dy/2.0));
			Ez_pp(i,j) = Ez_bar(i,j) + dx/8.0*(t1-t2) + dy/8.0*(t3-t4);
		}
	}
}

void CT_EMF_Stats(const char* name, Array &Ez_pp){
	double l2E=0, mxE=0, l2B=0, mxB=0;
	for(int j = ghost_cells; j < nyCells + ghost_cells; j++)
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			double e = Ez_pp(i,j);
			l2E += dx*dy*e*e; if(fabs(e) > mxE) mxE = fabs(e);
			double dBx = -(Ez_pp(i,j) - Ez_pp(i,j-1)) / dy;   // 面 (i+1/2, j) 的注入率
			double dBy =  (Ez_pp(i,j) - Ez_pp(i-1,j)) / dx;   // 面 (i, j+1/2)
			double m2 = dBx*dBx + dBy*dBy;
			l2B += dx*dy*m2; if(sqrt(m2) > mxB) mxB = sqrt(m2);
		}
	printf("[CT-DIAG] %-8s  Ez: L2=%.3e max=%.3e | dB/dt: L2=%.3e max=%.3e\n",
		name, sqrt(l2E), mxE, sqrt(l2B), mxB);
}

void Diagnose_CT_EMF(Grid& u_con){
	Array Cp0(nxCells,nyCells), C0p(nxCells,nyCells);      // 中心平均部分
	Array Rp0(nxCells,nyCells), R0p(nxCells,nyCells);      // Rusanov 耗散部分
	Array Ez_tilde(nxCells,nyCells), Zero(nxCells,nyCells);
	Array EppC(nxCells,nyCells), EppR(nxCells,nyCells), EppT(nxCells,nyCells), EppTot(nxCells,nyCells);
	for(int i = 0; i < nxCells + 2*ghost_cells; i++)
		for(int j = 0; j < nyCells + 2*ghost_cells; j++){
			Cp0(i,j)=C0p(i,j)=Rp0(i,j)=R0p(i,j)=Zero(i,j)=Ez_tilde(i,j)=0.0;
			EppC(i,j)=EppR(i,j)=EppT(i,j)=EppTot(i,j)=0.0;
		}
	double Q_L[NUM_VARS],Q_R[NUM_VARS],fL[NUM_VARS],fR[NUM_VARS],gL[NUM_VARS],gR[NUM_VARS];
	Grid w_xL(nxCells,nyCells), w_xR(nxCells,nyCells), w_yL(nxCells,nyCells), w_yR(nxCells,nyCells);
	TVDReconstruction(u_con, w_xL, w_xR, w_yL, w_yR);
	for(int i = 0; i < nxCells + 2*ghost_cells - 1; i++){
		for(int j = 0; j < nyCells + 2*ghost_cells - 1; j++){
			for(int k = 0; k < NUM_VARS; k++){ Q_L[k]=w_xR(i,j,k); Q_R[k]=w_xL(i+1,j,k); }
			f_full_con(fL,Q_L); f_full_con(fR,Q_R);
			Cp0(i,j) = -(0.5*(fL[6]+fR[6]));
			Rp0(i,j) =  0.5*fmax(fabs(Q_L[1]/Q_L[0]),fabs(Q_R[1]/Q_R[0]))*(Q_R[6]-Q_L[6]);
			for(int k = 0; k < NUM_VARS; k++){ Q_L[k]=w_yR(i,j,k); Q_R[k]=w_yL(i,j+1,k); }
			g_full_con(gL,Q_L); g_full_con(gR,Q_R);
			C0p(i,j) =  0.5*(gL[5]+gR[5]);
			R0p(i,j) = -0.5*fmax(fabs(Q_L[2]/Q_L[0]),fabs(Q_R[2]/Q_R[0]))*(Q_R[5]-Q_L[5]);
		}
	}
	bound_cond_array(Cp0); bound_cond_array(C0p);
	bound_cond_array(Rp0); bound_cond_array(R0p);
	for(int j = 0; j < nyCells + 2*ghost_cells; j++)
		for(int i = 0; i < nxCells + 2*ghost_cells; i++)
		Ez_tilde(i,j) = u_con(i,j,2)/u_con(i,j,0)*u_con(i,j,5)
	- u_con(i,j,1)/u_con(i,j,0)*u_con(i,j,6);
	CT_CornerAssemble(EppC, Cp0, C0p, Zero);       // 中心平均的贡献
	CT_CornerAssemble(EppR, Rp0, R0p, Zero);       // Rusanov 的贡献
	CT_CornerAssemble(EppT, Zero, Zero, Ez_tilde); // GS05 修正里 Ez_tilde 的贡献
	for(int j = 1; j < nyCells + 2*ghost_cells - 2; j++)
		for(int i = 1; i < nxCells + 2*ghost_cells - 2; i++)
		EppTot(i,j) = EppC(i,j) + EppR(i,j) + EppT(i,j);
	CT_EMF_Stats("central", EppC);
	CT_EMF_Stats("rusanov", EppR);
	CT_EMF_Stats("tilde",   EppT);
	CT_EMF_Stats("TOTAL",   EppTot);
}

void A_initial(Grid &A){
	// 注意: i 和 j 都要覆盖全部网格(包括ghost cells)
	for (int i = 0; i < nxCells + 2 * ghost_cells; i++) {
		for (int j = 0; j < nyCells + 2 * ghost_cells; j++) {
			A(i,j,0) = 0;
			A(i,j,1) = 0;
			double Ax_val = 0.0, Ay_val = 0.0, Az_val = 0.0;
			double x = x_0 + (i - ghost_cells + 0.5) * dx;
			double y = y_0 + (j - ghost_cells + 0.5) * dy;
			switch(test) {
				case test0: 
				case shocktest:
				case convergencetest: 
				case test2d0:
				case testdoubleshear:
					break;
				case testMHDJet:
					// B = (0, B0, 0)  ->  Az = Bx*y - By*x = -B0 * x (线性, bound_cond_A 外推精确)
					Az_val = -jet_B0 * x;
					break;
				// ============================================================
				// 1D x方向 Riemann Problems
				// Az = Bx*y - By*(x - x_d),  Ay = Bz*(x - x_d),  Ax = 0
				// Bx 在两侧相同 (div B = 0 的1D约束)
				// ============================================================
				
				case test1: {  // x_d = 0, Bx=0.75
					double Bx0 = 0.75;
					double By_here = (x < 0) ? 1.0 : -1.0;
					Az_val = Bx0 * y - By_here * (x - 0.0);
					Ay_val = 0.0;  // Bz = 0
					break;
				}
					
				case test2: {  // x_d = -0.1, Bx = 2/sqrt(4pi) ≈ 0.564190
					double xd = -0.1;
					double Bx0 = 0.564190;
					if (x < xd) {
						Az_val = Bx0 * y - 1.015541 * (x - xd);
						Ay_val = 0.564190 * (x - xd);
					} else {
						Az_val = Bx0 * y - 1.135262 * (x - xd);
						Ay_val = 0.564923 * (x - xd);
					}
					break;
				}
					
				case test3: {  // x_d = -0.1, Bx = 1.1
					double xd = -0.1;
					double Bx0 = 1.1;
					if (x < xd) {
						Az_val = Bx0 * y - 1.0 * (x - xd);
						Ay_val = 0.0 * (x - xd);
					} else {
						Az_val = Bx0 * y - 0.785887 * (x - xd);
						Ay_val = 0.618370 * (x - xd);
					}
					break;
				}
					
				case test4: {  // x_d = 0, Bx = 1.3
					double Bx0 = 1.3;
					double By_here = (x < 0) ? 1.0 : -1.0;
					Az_val = Bx0 * y - By_here * (x - 0.0);
					Ay_val = 0.0;  // Bz = 0
					break;
				}
					
				case test5: {  // x_d = 0
					double s4p = sqrt(4.0 * M_PI);
					double Bx0 = 0.05 / s4p;
					if (x < 0) {
						Az_val = Bx0 * y - (-2.0 / s4p) * x;
						Ay_val = (-1.0 / s4p) * x;
					} else {
						Az_val = Bx0 * y - (2.0 / s4p) * x;
						Ay_val = (1.0 / s4p) * x;
					}
					break;
				}
					
				case test6: {  // x_d = 0, B相同两侧 → A连续
					double s4p = sqrt(4.0 * M_PI);
					double Bx0 = 4.0 / s4p;
					double By0 = 4.0 / s4p;
					double Bz0 = 1.0 / s4p;
					Az_val = Bx0 * y - By0 * x;
					Ay_val = Bz0 * x;
					break;
				}
					
				case test7: {  // x_d = 0
					double s4p = sqrt(4.0 * M_PI);
					double Bx0 = 1.0 / s4p;
					double Bz0 = 1.0 / s4p;  // 两侧相同
					if (x < 0) {
						Az_val = Bx0 * y - (-1.0 / s4p) * x;  // By_L = -1/s4p
					} else {
						Az_val = Bx0 * y - (1.0 / s4p) * x;   // By_R = +1/s4p
					}
					Ay_val = Bz0 * x;
					break;
				}
					
					// ============================================================
					// 1D y方向 Riemann Problem
					// Az = -By*x + Bx*(y - y_d),  Ax = -Bz*(y - y_d),  Ay = 0
					// By 在两侧相同 (div B = 0 的1D y约束)
					// ============================================================
					
				case test1y: {  // y_d = 0, By = 0.75
					double By0 = 0.75;
					double Bx_here = (y < 0) ? 1.0 : -1.0;
					Az_val = -By0 * x + Bx_here * (y - 0.0);
					Ay_val = 0.0;
					break;
				}
				case testOrszag_Tang:
					Az_val = cos(2.0*M_PI*y)/(2.0*M_PI) + cos(4.0*M_PI*x)/(4.0*M_PI);
					break;
				
				case testblastwave:
					Az_val = (1.0/sqrt(4.0*M_PI) * 70.0) * y;
					break;
				case testShockCloud:
					if(x < 0.6){
						Az_val = -2.1826182 * x;
						Ay_val = -2.1826182 * x;
					}else{
						Az_val = -1.30957092 - 0.56418958 * (x-0.6);
						Ay_val = -1.30957092 + 0.56418958 * (x-0.6);
					}
					break;
				case testblastwave2:
					Az_val = (1.0/sqrt(4.0*M_PI) * 1000.0) * y;
					break;
				case testRotor:
					Az_val = (2.5/sqrt(4.0*M_PI)) * y;
					break;
				
				case testFieldLoop: {
					double C = A0 / sqrt(4.0*M_PI);
					double R = 0.3;
					double r = sqrt(x*x + y*y);
					Az_val = (r < R) ? C*(R-r) : 0.0;
					break;
				}
					
				case testKHInstability:
					Az_val = 0.1 * A0 * y;
					break;
				
				case testRK_Riemann: {  // x-Riemann,By 跨 x=0 跳变(3.6->4)/s4p;Bx,Bz 两侧相同
					double s4p = sqrt(4.0 * M_PI);
					double Bx0 = 2.0 / s4p;
					double Bz0 = 2.0 / s4p;
					if (x <= 0) Az_val = Bx0 * y - (3.6 / s4p) * x;
					else        Az_val = Bx0 * y - (4.0  / s4p) * x;
					Ay_val = Bz0 * x;
					break;
				}
				case testPlasmaEquilibrium:
				case testAdvectedPlasmaEquilibrium:
					Az_val = 0.0;     // 面内 Bx=By=0 -> Az=0;Bz 作为独立分量演化,不经 Az
					break;
				case testMHD_vortex: {  // 换成你实际的 enum 名字
					double r2 = x*x + y*y;
					Az_val = (1.0 / (2.0*M_PI)) * exp(0.5 * (1.0 - r2));
					// Bz = 0，所以 Ax = Ay = 0，不用改
					break;
				}
				default:
					break;
				
			}
			
			A(i, j, 0) = Ax_val;
			A(i, j, 1) = Ay_val;
			A(i, j, 2) = Az_val;
		}
	}
	bound_cond_A(A);
}

inline double minmod(double a, double b) {
	if (a * b <= 0.0) return 0.0;
	return (fabs(a) < fabs(b)) ? a : b;
}

// ============================================================================
//  Unstaggered CT — Rossmanith (MPACT) 式高分辨率迎风势平流
//  直接整体替换原 Unstaggered_CT 函数, 签名不变。
//
//  变化: 不再经由 Godunov 角点 EMF (重构+Rusanov) 更新 A, 而是按论文 3.6.2 /
//  Rossmanith 的 Step III 把 A 当作被 u^{n+1} 平流的标量, 用 LeVeque 波传播型
//  高分辨率迎风格式离散  ∂A/∂t + u ∂A/∂x + v ∂A/∂y = 0:
//    一阶迎风 fluctuation  +  带限制器的二阶 Lax–Wendroff 修正
//  耗散作用在 A 的限制器上 (光滑区 O(Δx²)), 不再作用在 B 的胞心跳变上。
//
//  依赖: 你已有的 inline minmod(a,b), bound_cond_A, bound_cond_u。
//  注意: 内部对 A_s 调用一次 bound_cond_A 以刷新 ghost (差分需要 ±2 邻居,
//        ghost_cells=3 足够)。若 A_s 的 ghost 已在别处维护, 该调用无副作用。
//  稳定性: 二阶修正含 (1 − |u|h/Δ) 因子, 要求方向 Courant |u|h/Δ ≤ 1;
//        材料 CFL 取步长时自动满足 (CFL=0.45 < 1)。若你在某测试用更大的
//        等效 Courant, 修正因子变负仍稳定但精度降; 可在此 clamp 到 [0,1]。
// ============================================================================

inline double vanleer_lim(double a, double b){
	if(a*b <= 0.0) return 0.0;
	return 2.0*a*b/(a + b);
}

void Unstaggered_CT(Grid& u_con, Grid& A, Grid& A_s, double h) {
	
	bound_cond_A(A_s);   // 差分模板需要 A_s 的 ghost 值
	
	// ===== Step III: 高分辨率迎风平流  A^{n+1} = A_s - h * (u·∇A) =====
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			
			// ---------- x 方向 ----------
			// 面速度 (u^{n+1}, 来自动量): 左面 i-1/2, 右面 i+1/2
			double uL = 0.5*(u_con(i-1,j,1)/u_con(i-1,j,0) + u_con(i,j,1)/u_con(i,j,0));
			double uR = 0.5*(u_con(i,j,1)/u_con(i,j,0)   + u_con(i+1,j,1)/u_con(i+1,j,0));
			
			// 面上的跳变 (用 A_s, 与阶段基一致)
			double dAL  = A_s(i,j,2)   - A_s(i-1,j,2);      // i-1/2
			double dAR  = A_s(i+1,j,2) - A_s(i,j,2);        // i+1/2
			double dALL = A_s(i-1,j,2) - A_s(i-2,j,2);      // i-3/2
			double dARR = A_s(i+2,j,2) - A_s(i+1,j,2);      // i+3/2
			
			// 一阶迎风 fluctuation: u⁺·ΔA_{i-1/2} + u⁻·ΔA_{i+1/2}
			double fluct_x = fmax(uL, 0.0)*dAL + fmin(uR, 0.0)*dAR;
			
			// 二阶限制修正 (LeVeque): F̃_{f} = ½|u_f|(1-|u_f|h/dx)·minmod(ΔA_f, ΔA_upwind)
			double corrL = 0.5*fabs(uL)*(1.0 - fabs(uL)*h/dx)
			* vanleer_lim(dAL, (uL >= 0.0 ? dALL : dAR));
			double corrR = 0.5*fabs(uR)*(1.0 - fabs(uR)*h/dx)
			* vanleer_lim(dAR, (uR >= 0.0 ? dAL : dARR));
			
			double dAdt_x = fluct_x/dx + (corrR - corrL)/dx;
			
			// ---------- y 方向 ----------
			double vL = 0.5*(u_con(i,j-1,2)/u_con(i,j-1,0) + u_con(i,j,2)/u_con(i,j,0));
			double vR = 0.5*(u_con(i,j,2)/u_con(i,j,0)   + u_con(i,j+1,2)/u_con(i,j+1,0));
			
			double dAB  = A_s(i,j,2)   - A_s(i,j-1,2);      // j-1/2
			double dAT  = A_s(i,j+1,2) - A_s(i,j,2);        // j+1/2
			double dABB = A_s(i,j-1,2) - A_s(i,j-2,2);
			double dATT = A_s(i,j+2,2) - A_s(i,j+1,2);
			
			double fluct_y = fmax(vL, 0.0)*dAB + fmin(vR, 0.0)*dAT;
			
			double corrB = 0.5*fabs(vL)*(1.0 - fabs(vL)*h/dy)
			* vanleer_lim(dAB, (vL >= 0.0 ? dABB : dAT));
			double corrT = 0.5*fabs(vR)*(1.0 - fabs(vR)*h/dy)
			* vanleer_lim(dAT, (vR >= 0.0 ? dAB : dATT));
			
			double dAdt_y = fluct_y/dy + (corrT - corrB)/dy;
			
			// ---------- 更新 ----------
			A(i,j,2) = A_s(i,j,2) - h * (dAdt_x + dAdt_y);
		}
	}
	bound_cond_A(A);
	
	// ===== Step IV: B = curl A (与原实现相同) =====
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			u_con(i,j,5) =  (A(i,j+1,2) - A(i,j-1,2)) / (2.0*dy);
			u_con(i,j,6) = -(A(i+1,j,2) - A(i-1,j,2)) / (2.0*dx);
		}
	}
	bound_cond_u(u_con);
}

void Unstaggered_CT_Consistent(Grid& u_con, Grid& A, Grid& A_s, double h) {
	Array Ez_p0(nxCells,nyCells), Ez_0p(nxCells,nyCells);
	Array Ez_tilde(nxCells,nyCells), Ez_pp(nxCells,nyCells);
	
	// 面 Ez = 一致化通量的对偶 (与 Staggered_CT_Consistent 完全相同)
	for(int i = 0; i < nxCells + 2*ghost_cells; i++){
		for(int j = 0; j < nyCells + 2*ghost_cells; j++){
			Ez_p0(i,j) = -(*gFluxBy_p0)(i,j);
			Ez_0p(i,j) = +(*gFluxBx_0p)(i,j);
			Ez_pp(i,j) = 0.0;
		}
	}
	
	// 胞心参考 Ez (GS05 修正用), u_con 此刻 = (u^{n+1}, B̃^{n+1})
	for(int j = 0; j < nyCells + 2*ghost_cells; j++)
		for(int i = 0; i < nxCells + 2*ghost_cells; i++)
		Ez_tilde(i,j) = u_con(i,j,2)/u_con(i,j,0)*u_con(i,j,5)
	- u_con(i,j,1)/u_con(i,j,0)*u_con(i,j,6);
	
	// GS05 角点装配 (复用)
	CT_CornerAssemble(Ez_pp, Ez_p0, Ez_0p, Ez_tilde);
	
	// 胞心 Ez = 四角平均, 更新 A  (结构与你原版 Step 2 相同)
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			double Ez_cell = 0.25*(Ez_pp(i,j) + Ez_pp(i-1,j)
				+ Ez_pp(i,j-1) + Ez_pp(i-1,j-1));
			A(i,j,2) = A_s(i,j,2) - h * Ez_cell;
		}
	}
	bound_cond_A(A);
	
	// B = curl A  (与你原版 Step 3 相同)
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			u_con(i,j,5) =  (A(i,j+1,2) - A(i,j-1,2)) / (2.0*dy);
			u_con(i,j,6) = -(A(i+1,j,2) - A(i-1,j,2)) / (2.0*dx);
		}
	}
	bound_cond_u(u_con);
}
