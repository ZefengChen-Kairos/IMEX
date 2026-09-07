// Test configuration and the IMEX time-marching driver.

#include "declarations.hpp"

int main(int argc, char *argv[]){
	if(const char* env = std::getenv("MHD_OUTPUT_ROOT")) output_root = env;
	else if(argc > 1) output_root = argv[1];

	switch (test) {
		case test0:
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.2;
			gamma_plasma = 1.4;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case test1: case test1y:
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.1;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case test2:
			x_0 = -0.5, x_1 = 0.5, y_0 = -1, y_1 = 1,z_0 = -1, z_1 = 1;
			time_stop = 0.2;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case test3:
			//x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.15;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case test4:
			//x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.16;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case test5:
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.04;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case test6:
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.03;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case test7:
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.25;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case test8:
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.16;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case test9:
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.08;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case test10:
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.16;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case RKEOS1:
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.1;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case RKEOS2:
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.2;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case RKEOS3:
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.008;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case RKEOS4:
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.1;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		
		case shocktest:
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.8;
			gamma_plasma = 1.4;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case convergencetest:
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 1.0;
			x_bc = x_Periodic;
			y_bc = y_Periodic;
			break;
		case test2d0:
			x_0 = 0, x_1 = 1, y_0 = 0, y_1 = 1,z_0 = -1, z_1 = 1;
			time_stop = 0.25;
			gamma_plasma = 1.4;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case testOrszag_Tang:
			x_0 = 0, x_1 = 1, y_0 = 0, y_1 = 1,z_0 = -1, z_1 = 1;
			time_stop = 1.0;
			x_bc = x_Periodic;
			y_bc = y_Periodic;
			gamma_plasma = 1.4;
			break;
		case testblastwave:
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.01;
			gamma_plasma = 1.4;
			x_bc = x_Transmissive;   // match CASTOR: rho/v/p Dirichlet(ambient) + B transmissive
			y_bc = y_Transmissive;
			break;
		case testShockCloud:
			x_0 = 0, x_1 = 1, y_0 = 0, y_1 = 1,z_0 = -1, z_1 = 1;
			time_stop =  0.06;
			gamma_plasma = 5.0/3.0;
			x_bc = x_Transmissive;   // match CASTOR: rho/v/p Dirichlet(ambient) + B transmissive
			y_bc = y_Transmissive;
			break;
		case testblastwave2:
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.001;
			gamma_plasma = 1.4;
			x_bc = x_Transmissive;   // match CASTOR: rho/v/p Dirichlet(ambient) + B transmissive
			y_bc = y_Transmissive;
			break;
		case testblastwave1D:
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.01;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case testRotor:
			x_0 = -0.5, x_1 = 0.5, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 0.25;
			gamma_plasma = 1.4;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case testFieldLoop:
			x_0 = -1, x_1 = 1, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = 1;
			gamma_plasma = 5.0/3.0;
			x_bc = x_Periodic;
			y_bc = y_Periodic;
			break;
		case testKHInstability:
			x_0 = 0, x_1 = 2, y_0 = -0.5, y_1 = 0.5,z_0 = -1, z_1 = 1;
			time_stop = (4.8 / A0 ) / 6.0;
			gamma_plasma = 1.4;
			x_bc = x_Periodic;
			y_bc = y_Periodic;
			break;
		case testdoubleshear:
			x_0 = 0, x_1 = 1, y_0 = 0, y_1 = 1,z_0 = -1, z_1 = 1;
			time_stop = 1.8;
			gamma_plasma = 1.4;
			x_bc = x_Periodic;
			y_bc = y_Periodic;
			break;
		case testRK_Riemann:
			x_0 = -0.5, x_1 = 0.5, y_0 = -1, y_1 = 1, z_0 = -1, z_1 = 1;
			time_stop = 0.1;
			x_bc = x_Transmissive;
			y_bc = y_Transmissive;
			break;
		case testPlasmaEquilibrium:
			x_0 = -3.0, x_1 = 3.0, y_0 = -3.0, y_1 = 3.0, z_0 = -1, z_1 = 1;
			time_stop = 2.0;
			gamma_plasma = 5.0/3.0;
			x_bc = x_Transmissive;     // 静态平衡,零流;透射即可
			y_bc = y_Transmissive;
			break;
		case testAdvectedPlasmaEquilibrium:
			x_0 = -3.0, x_1 = 3.0, y_0 = -3.0, y_1 = 3.0, z_0 = -1, z_1 = 1;
			time_stop = 400.0;         // 平流绕域一圈
			gamma_plasma = 5.0/3.0;
			x_bc = x_Periodic;         // 平流版用周期(让环绕回来)
			y_bc = y_Periodic;
			break;
		case testMHD_vortex:
			x_0 = -10.0, x_1 = 10.0, y_0 = -10.0, y_1 = 10.0, z_0 = -1, z_1 = 1;
			time_stop = 1000.0;         // 平流绕域一圈
			gamma_plasma = 5.0/3.0;
			x_bc = x_Periodic;         // 平流版用周期(让环绕回来)
			y_bc = y_Periodic;
			break;
		case testMHDJet:
			// Example 4.8: Mach 800 MHD jet, [-0.5,0.5]x[0,1.5], 400x600
			// (nxCells/nyCells 是编译期常量,记得顶部改 400/600)
			x_0 = -0.5, x_1 = 0.5, y_0 = 0.0, y_1 = 1.5, z_0 = -1, z_1 = 1;
			time_stop = 0.00015;                 // Mach 800 标准终止时间
			gamma_plasma = 1.4;
			jet_rho = gamma_plasma;            // 入流密度 = gamma
			jet_p   = 1.0;
			x_bc = x_Transmissive;             // 左右: outflow
			y_bc = y_JetInflow;                // 底: |x|<0.05 inflow, 其余 outflow;顶: outflow
			break;
	}
	
	dx = (x_1 - x_0)/(nxPoints-1);
	dy = (y_1 - y_0)/(nyPoints-1);
	dz = (z_1 - z_0)/(nzPoints-1);

	t = 0, dt = 0;
	Grid u_pri(nxCells, nyCells), u_con(nxCells, nyCells), u_con_n(nxCells, nyCells), Qe(nxCells, nyCells),Qe_pri(nxCells, nyCells), Qs(nxCells, nyCells), w_xL(nxCells, nyCells), w_xR(nxCells, nyCells),w_yL(nxCells, nyCells), w_yR(nxCells, nyCells),A(nxCells,nyCells);
	
	// 新增A的IMEX相关变量
	Grid A_n(nxCells, nyCells);           // 时间步开头的A
	Tensor HsA(RKs, nxCells, nyCells), HsB(RKs, nxCells, nyCells);    // A的时间导数历史，需要3分量
	Array Bx_p0(nxCells,nyCells), By_0p(nxCells,nyCells),Bx_p0_n(nxCells,nyCells), By_0p_n(nxCells,nyCells), Az(nxCells, nyCells), Bx(nxCells,nyCells),By(nxCells,nyCells), Bz(nxCells,nyCells);
	gFluxBy_p0 = new Array(nxCells, nyCells);
	gFluxBx_0p = new Array(nxCells, nyCells);
	Array Bx_p0_Qe(nxCells,nyCells), By_0p_Qe(nxCells,nyCells),Bx_p0_Qs(nxCells,nyCells), By_0p_Qs(nxCells,nyCells);
	u_initialization(u_pri, Bx_p0, By_0p);
	bound_cond_staggered(Bx_p0, 0);   // x-normal: 保留边界法向面(CT div-free)
	bound_cond_staggered(By_0p, 1); 
	upri2ucon(u_con, u_pri);
	A_initial(A);
	
	
	std::vector<double> bRK(RKs);
	Array Ae(RKs, RKs), Ai(RKs, RKs);
	Tensor Hs(RKs, nxCells, nyCells);
	Array shock_eta(nxCells, nyCells);
	std::vector<std::pair<int, double>> div_error;
	std::vector<std::pair<int, double>> B_inf;
	std::vector<std::pair<int, double>> dt_comp;
	const int plotINTERVAL = 10000000;
	
	std::vector<double> IntermediatePlotTimes = {
		0.25 * time_stop,
		1.0/3.0 * time_stop,
		0.5 * time_stop,
		2.0/3.0 * time_stop,
		0.75 * time_stop,
		1.00 * time_stop,
	};
	
	bool isFinalTime = false;
	
	IntermediatePlotIntermediatePlot(u_con, u_pri,  timestep, plotINTERVAL,  IntermediatePlotTimes, isFinalTime);

	if(RKs == 1){
		Ae(0,0) = 0.0;
		Ai(0,0) = 1.0;
		bRK[0]  = 1.0;
	}else if(RKs == 2){
		double bb = 1.0 - 1.0 / sqrt(2.0);
		double cc = 1.0 / (2.0 * bb);
		Ae(0,0) = 0.0;
		Ae(0,1) = 0.0;
		
		Ae(1,0) = cc;
		Ae(1,1) = 0.0;
		
		Ai(0,0) = bb;
		Ai(0,1) = 0.0;
		
		Ai(1,0) = 1 - bb;
		Ai(1,1) = bb;
		
		bRK[0] = 1.0 - bb;
		bRK[1] = bb;
	}else if(RKs == 3){
		double r12 = 0.5, r13=1.0/3.0, r14=1.0/4.0;
		Ae(0,0) = 0.0;
		Ae(0,1) = 0.0;
		Ae(0,2) = 0.0;
		
		Ae(1,0) = r12;
		Ae(1,1) = 0.0;
		Ae(1,2) = 0.0;
		
		Ae(2,0) = r12;
		Ae(2,1) = r12;
		Ae(2,2) = 0.0;
		
		Ai(0,0) = r14;
		Ai(0,1) = 0.0;
		Ai(0,2) = 0.0;
		
		Ai(1,0) = 0.0;
		Ai(1,1) = r14;
		Ai(1,2) = 0.0;
		
		Ai(2,0) = r13;
		Ai(2,1) = r13;
		Ai(2,2) = r13;
		
		bRK[0] = r13;
		bRK[1] = r13;
		bRK[2] = r13;
	}
	div_error.emplace_back(timestep, cal_divergence_error(u_con, Bx_p0, By_0p));
	double B_inf_m = 0;
	double divB = 0;
	for(int i = 1; i < nxCells + 2*ghost_cells - 1; i++){
		for(int j = 1; j < nyCells + 2*ghost_cells - 1; j++){
			// fix: 原循环 i=0..n+2g-1 访问 i±1 越界(-O2 下随机段错);收缩一圈
			divB = (u_con(i+1,j,5) - u_con(i-1,j,5)) / (2.0*dx) + (u_con(i,j+1,6) - u_con(i,j-1,6)) / (2.0*dy);
			if(fabs(divB) > B_inf_m){
				B_inf_m = fabs(divB);
			}
		}
	}
	B_inf.emplace_back(timestep, B_inf_m);
	
	do{
		timestep++;
		g_bc_prim = true;  bound_cond_u(u_pri);  g_bc_prim = false;  // primitive ghosts (slot 4 = p)
		calch(u_pri);
		calcp(u_pri);
		calca(u_pri);
		
		dt = fmin(caldt(u_pri), time_stop - t);
		if(g_fixed_dt > 0.0) dt = fmin(g_fixed_dt, time_stop - t);   // 对比用:强制固定 dt(=CASTOR 的 FIXDT)
		if(fabs(dt - (time_stop - t)) < 1e-8){
			isFinalTime = true;
		}
		
		Grid A_n_snapshot(nxCells, nyCells);
		copy(A_n_snapshot, A);
		
		
		MagneticEigenvalues(u_pri);
		EnergyEigenvalues(u_pri);
		t = t + dt;
		//std::cout << B_max << " " << U_max << endl;
		
		upri2ucon(u_con, u_pri);
		
		copy(u_con_n,u_con);
		
		for(int i = 0; i < nxCells + 2*ghost_cells; i++){
			for(int j = 0; j < nyCells + 2*ghost_cells; j++){
				Bx_p0_n(i,j) = Bx_p0(i,j); 
				By_0p_n(i,j) = By_0p(i,j);
			}
		}
		
		for(int rks = 0; rks < RKs; rks++){
			copy(Qe,u_con_n);
			copy(Qs,u_con_n);
			
			Grid A_s(nxCells, nyCells), A_e(nxCells, nyCells);
			copy(A_s, A_n_snapshot);
			copy(A_e, A_n_snapshot);
			
			for(int i = 0; i < nxCells + 2*ghost_cells; i++){
				for(int j = 0; j < nyCells + 2*ghost_cells; j++){
					Bx_p0_Qe(i,j) = Bx_p0_n(i,j); 
					By_0p_Qe(i,j) = By_0p_n(i,j);
					Bx_p0_Qs(i,j) = Bx_p0_n(i,j);
					By_0p_Qs(i,j) = By_0p_n(i,j);
				}
			}
						
			ucon2upri(Qe_pri, Qe);
			
			for(int r = 0; r < rks; r++){
				for(int i = 0; i < nxCells + 2 * ghost_cells; i++){
					for(int j = 0; j < nyCells + 2 * ghost_cells; j++){
						for(int k = 0; k < NUM_VARS; k++){
							Qe(i,j,k) = Qe(i,j,k) + Ae(rks,r) * dt * Hs(r,i,j,k);
							Qs(i,j,k) = Qs(i,j,k) + Ai(rks,r) * dt * Hs(r,i,j,k);
						}
					}
				}
			}
			for(int r = 0; r < rks; r++){
				for(int i = 0; i < nxCells + 2*ghost_cells; i++){
					for(int j = 0; j < nyCells + 2*ghost_cells; j++){
						for(int k = 0; k < 3; k++){   // A 只有 3 分量
							A_e(i,j,k) += Ae(rks,r) * dt * HsA(r,i,j,k);
							A_s(i,j,k) += Ai(rks,r) * dt * HsA(r,i,j,k);
						}
						Bx_p0_Qe(i,j) += Ae(rks,r) * dt * HsB(r,i,j,0); 
						By_0p_Qe(i,j) += Ae(rks,r) * dt * HsB(r,i,j,1);
						Bx_p0_Qs(i,j) += Ai(rks,r) * dt * HsB(r,i,j,0);
						By_0p_Qs(i,j) += Ai(rks,r) * dt * HsB(r,i,j,1);
					}
				}
			}
			
						
			copy(A, A_s);
			
			adt =  Ai(rks,rks) * dt;
			
			
			TVDReconstruction(Qe,w_xL,w_xR,w_yL,w_yR);
			
			// Convective Subsystem F^c
			RusanovC_Update(u_con, Qs, w_xL,w_xR,w_yL,w_yR, adt);
			bound_cond_u(u_con);

			// Magnetic Subsystem F^B
			if(nxCells == 1 || nyCells == 1){
				MagSubsystem_UpdateB(u_con, Qs, Qe, w_xL, w_xR, w_yL, w_yR, A, adt,rks);
			}else{
				switch(magfield){
					case StaggerCT: case UnStaggerCT: case NoDivFree:
						MagSubsystem_UpdateB(u_con, Qs, Qe, w_xL, w_xR, w_yL, w_yR, A, adt, rks);
						break;
					case PotentialA:
						MagSubsystem_UpdateA(u_con, Qs, Qe, w_xL, w_xR, w_yL, w_yR, A,A_s, adt, rks);
						break;
				}
			}
			bound_cond_u(u_con);
			
			// Store Magnetic Field
			for(int i = 0; i < 2 * ghost_cells + nxCells; i++){
				for(int j = 0; j < 2 * ghost_cells + nyCells; j++){
					Bx(i,j) = u_con(i,j,5);
					By(i,j) = u_con(i,j,6);
					Bz(i,j) = u_con(i,j,7);
				}
			}
			
			Grid u_test(nxCells, nyCells);
			for(int i = 0; i < 2 * ghost_cells + nxCells; i++){
				for(int j = 0; j < 2 * ghost_cells + nyCells; j++){
					w_xL(i,j,5) = u_con(i,j,5);
					w_xL(i,j,6) = u_con(i,j,6);
					w_xL(i,j,7) = u_con(i,j,7);
					
					w_xR(i,j,5) = u_con(i,j,5);
					w_xR(i,j,6) = u_con(i,j,6);
					w_xR(i,j,7) = u_con(i,j,7);
					
					w_yL(i,j,5) = u_con(i,j,5);
					w_yL(i,j,6) = u_con(i,j,6);
					w_yL(i,j,7) = u_con(i,j,7);
					
					w_yR(i,j,5) = u_con(i,j,5);
					w_yR(i,j,6) = u_con(i,j,6);
					w_yR(i,j,7) = u_con(i,j,7);
					
				}
			}
			
			switch(formula){
				case EnergyFormula:
					Energy_flux_Update(u_con, Qs, w_xL, w_xR, w_yL, w_yR, adt);
					break;
				case PressureFormula:
					Pressure_flux_Update(u_con, Qs, w_xL, w_xR, w_yL, w_yR, adt);
					break;
			}
			bound_cond_u(u_con); 
			
			// Pressure Subsystem F^p + Momentum Update
			switch(formula){
				case EnergyFormula:
					PreSubsystem_Energy(u_con, Qs, Qe, w_xL, w_xR, w_yL, w_yR, adt,rks);
					break;
				case PressureFormula:
					PreSubsystem_Pressure(u_con, Qs, Qe, w_xL, w_xR, w_yL, w_yR, adt,rks);
					break;
			}
			
			
			
					
			
			bound_cond_u(u_con); 
			
			
			// divergence-free
			if(nxCells != 1 && nyCells != 1){
				switch(magfield){
					case StaggerCT:
						//Staggered_CT(u_con, Bx_p0, By_0p, Bx_p0_Qs,By_0p_Qs,adt);
						Staggered_CT_Consistent(u_con, Bx_p0, By_0p, Bx_p0_Qs,By_0p_Qs,adt);
						break;
					case UnStaggerCT:
//						Unstaggered_CT(u_con, A, A_s,adt);
						Unstaggered_CT_Consistent(u_con, A, A_s,adt);
						break;
					case PotentialA: case NoDivFree:
						break;
				}
			}
			
//			for(int i = 0; i < 2 * ghost_cells + nxCells; i++){
//				for(int j = 0; j < 2 * ghost_cells + nyCells; j++){
//					//					u_con(i,j,0) = Qe(i,j,0); 
//					//					u_con(i,j,1) = Qe(i,j,1); 
//					//					u_con(i,j,2) = Qe(i,j,2); 
//					//					u_con(i,j,3) = Qe(i,j,3); 
//					//					u_con(i,j,4) = Qe(i,j,4); 
//					u_con(i,j,5) = Qe(i,j,5); 
//					u_con(i,j,6) = Qe(i,j,6); 
//					u_con(i,j,7) = Qe(i,j,7); 
//				}
//			}
			
			bound_cond_u(u_con);
			
			// Modify Energy
			if(Energy_Modify){
				for(int i = ghost_cells; i < ghost_cells + nxCells; i++){
					for(int j = ghost_cells; j < ghost_cells + nyCells; j++){
						u_con(i,j,4) = u_con(i,j,4) + ((u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6) + u_con(i,j,7) * u_con(i,j,7)) - (Bx(i,j) * Bx(i,j) + By(i,j) * By(i,j) + Bz(i,j) * Bz(i,j))) / 2.0;
					}
				}
			}
			
			bound_cond_u(u_con);
			ucon2upri(u_pri, u_con);
			
			for(int i = 0; i < nxCells + 2 * ghost_cells; i++){
				for(int j = 0; j < nyCells + 2 * ghost_cells; j++){
					for(int k = 0; k < NUM_VARS; k++){
						Hs(rks,i,j,k) = (u_con(i,j,k) - Qs(i,j,k)) / adt;
					}
				}
			}
			
			for(int i = 0; i < nxCells + 2*ghost_cells; i++){
				for(int j = 0; j < nyCells + 2*ghost_cells; j++){
					for(int k = 0; k < 3; k++){
						HsA(rks,i,j,k) = (A(i,j,k) - A_s(i,j,k)) / adt;
					}
					HsB(rks,i,j,0) = (Bx_p0(i,j) - Bx_p0_Qs(i,j)) / adt;
					HsB(rks,i,j,1) = (By_0p(i,j) - By_0p_Qs(i,j)) / adt;
				}
			}
			
				
			bound_cond_u(u_con);
		}
				
		ucon2upri(u_pri, u_con);
		
		std::cout<<"[timestep, dt, time] = "  << timestep << " "<<dt << " " << t << ", Divergence Error: "<< cal_divergence_error(u_con, Bx_p0, By_0p)<< endl;
		{   // ===== 跨代码对比诊断(归一化无关量) =====
			double emin=1e30, mass=0, intE=0, mmax=0; int im=-1, jm=-1;
			for(int i=ghost_cells;i<nxCells+ghost_cells;i++) for(int j=ghost_cells;j<nyCells+ghost_cells;j++){
				double rho=u_con(i,j,0);
				double k=0.5*(u_con(i,j,1)*u_con(i,j,1)+u_con(i,j,2)*u_con(i,j,2)+u_con(i,j,3)*u_con(i,j,3))/rho;
				double m=0.5*(u_con(i,j,5)*u_con(i,j,5)+u_con(i,j,6)*u_con(i,j,6)+u_con(i,j,7)*u_con(i,j,7));
				double e=u_con(i,j,4)-k-m;
				mass+=rho; intE+=e; if(m>mmax)mmax=m; if(e<emin){emin=e;im=i;jm=j;}
			}
			double xe=x_0+(im-ghost_cells+0.5)*dx, ye=y_0+(jm-ghost_cells+0.5)*dy;
			std::cout<<"DIAG ts="<<timestep<<" t="<<t<<" emin="<<emin<<" at("<<xe<<","<<ye<<")"
			<<" mass="<<mass*dx*dy<<" intE="<<intE*dx*dy<<" magE_max="<<mmax<<std::endl;
		}
		IntermediatePlotIntermediatePlot(u_con, u_pri,  timestep, plotINTERVAL,  IntermediatePlotTimes, isFinalTime);
		div_error.emplace_back(timestep, cal_divergence_error(u_con, Bx_p0, By_0p));
		
		
		double B_inf_m = 0;
		double divB = 0;
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
				divB = (u_con(i+1,j,5) - u_con(i-1,j,5)) / (2.0*dx) + (u_con(i,j+1,6) - u_con(i,j-1,6)) / (2.0*dy);
				if(fabs(divB) > B_inf_m){
					B_inf_m = fabs(divB);
				}
			}
		}
		B_inf.emplace_back(timestep, B_inf_m);
		dt_comp.emplace_back(timestep, imex_dt/fully_dt);
		
		Diagnose_CT_EMF(u_con);   // t=0 精确初值上的 CT 残余分解
	}while(t < time_stop);
	IntermediatePlotIntermediatePlot(u_con, u_pri,  timestep, plotINTERVAL,  IntermediatePlotTimes, isFinalTime);
	
	ofstream outFilerho(Path + "/RK_rho.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFilerho << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << u_con(i,j,0) << endl;
		}
	}
	outFilerho.close();
	
	ofstream outFilevx(Path + "/RK_vx.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFilevx << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << u_pri(i,j,1) << endl;
		}
	}
	outFilevx.close();
	
	ofstream outFilevy(Path + "/RK_vy.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFilevy << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << u_pri(i,j,2) << endl;
		}
	}
	outFilevy.close();
	
	ofstream outFileE(Path + "/RK_E.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFileE << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << u_con(i,j,4) << endl;
		}
	}
	outFileE.close();
	
	
	ofstream outFilep(Path + "/RK_p.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFilep << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << (u_pri(i,j,4)) << endl;
		}
	}
	outFilep.close();
	
	ofstream outFileBx(Path + "/RK_Bx.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFileBx << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << (u_pri(i,j,5)) << endl;
		}
	}
	outFileBx.close();

	ofstream outFileBy(Path + "/RK_By.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFileBy << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << (u_pri(i,j,6)) << endl;
		}
	}
	outFileBy.close();
	
	ofstream outFileBz(Path + "/RK_Bz.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFileBz << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << (u_pri(i,j,7)) << endl;
		}
	}
	outFileBz.close();
	
	
	ofstream outFilev(Path + "/RK_v_mag.dat");
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			//			outFile << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6))/u_con(i,j,7)<< endl;
			
			outFilev << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << sqrt(u_pri(i,j,1) * u_pri(i,j,1) + u_pri(i,j,2) * u_pri(i,j,2)) << endl;
		}
	}
	outFilev.close();
	
	
	ofstream outFiledt_comp(Path + "/RK_dt_comp.dat");
	for(const auto& p : dt_comp)
		outFiledt_comp << p.first << " " << p.second << "\n";
	outFiledt_comp.close();
	
	ofstream outFilediv(Path + "/RK_div_error.dat");
	for(const auto& p : div_error)
		outFilediv << p.first << " " << p.second << "\n";
	outFilediv.close();
	
	ofstream outFileB_inf(Path + "/RK_B_inf.dat");
	for(const auto& p : B_inf)
		outFileB_inf << p.first << " " << p.second << "\n";
	outFileB_inf.close();
	
	if(eos_type == REDLICH_KWONG){
		Array T(nxCells,nyCells);
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
				T(i,j) = REDKWONG_computeTfromRhoAndP(u_pri(i,j,0), u_pri(i,j,4));
			}
		}
		bound_cond_array(T);
		ofstream outFileT(Path + "/RK_T.dat");
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
				outFileT << x_0 + dx * (i - ghost_cells + 0.5) << " " << y_0 + dy * (j - ghost_cells + 0.5) << " " << T(i,j) << endl;
			}
		}
		outFileT.close();
	}
}
