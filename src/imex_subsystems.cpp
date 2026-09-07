// IMEX subsystem assembly and updates; divergence-error diagnostic.

#include "declarations.hpp"

void MagSubsystem_UpdateA(Grid &u_con, Grid &Qs, Grid &Qe,Grid &/*w_xL*/, Grid &/*w_xR*/,Grid &/*w_yL*/, Grid &/*w_yR*/, Grid &A,Grid &A_s,double h, const int /*rks*/){
	Array shock_eta(nxCells,nyCells);
	Grid Qe_pri(nxCells, nyCells);
	Grid Qs_pri(nxCells, nyCells);
	ucon2upri(Qe_pri, Qe);
	ucon2upri(Qs_pri, Qs);
	Eigen::VectorXd b_U_vect(nxCells*nyCells), b_vect(3 * nxCells * nyCells), A_vect(3 * nxCells * nyCells), U_vect(nxCells * nyCells);
	
	
	[[maybe_unused]] typedef Eigen::SparseMatrix<double> SpMat;
	[[maybe_unused]] typedef Eigen::Triplet<double> T;
	
	
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			
			// B^n at cell (i,j)
			double Bxn = Qe(i,j,5);
			double Byn = Qe(i,j,6);
			double Bzn = Qe(i,j,7);
			
			// ρ^{n+1}
			double rho = u_con(i,j,0);
			
			// S = v^{**} - h/ρ * ∇p^n  (逐点, 用中心差分梯度G)
			double Sx = u_con(i,j,1) / rho
			- h / rho * (Qe_pri(i+1,j,4) - Qe_pri(i-1,j,4)) / (2.0 * dx);
			double Sy = u_con(i,j,2) / rho
			- h / rho * (Qe_pri(i,j+1,4) - Qe_pri(i,j-1,4)) / (2.0 * dy);
			double Sz = u_con(i,j,3) / rho;
			
			// A^{**} = A^n - h * (B^n × S)
			double Ax_star = A_s(i,j,0) - h * (Byn * Sz - Bzn * Sy);
			double Ay_star = A_s(i,j,1) - h * (Bzn * Sx - Bxn * Sz);
			double Az_star = A_s(i,j,2) - h * (Bxn * Sy - Byn * Sx);
			
			// 存入RHS向量
			b_vect((i - ghost_cells) + nxCells * (j - ghost_cells))                         = Ax_star;
			b_vect((i - ghost_cells) + nxCells * (j - ghost_cells) + nxCells * nyCells)     = Ay_star;
			b_vect((i - ghost_cells) + nxCells * (j - ghost_cells) + 2 * nxCells * nyCells) = Az_star;
			
			// 初始猜测 = A^n
			A_vect((i - ghost_cells) + nxCells * (j - ghost_cells))                         = A(i,j,0);
			A_vect((i - ghost_cells) + nxCells * (j - ghost_cells) + nxCells * nyCells)     = A(i,j,1);
			A_vect((i - ghost_cells) + nxCells * (j - ghost_cells) + 2 * nxCells * nyCells) = A(i,j,2);
		}
	}

	GMRES_Solver(A_vect, b_vect, Magnetic_Subsystem_potential_AX, u_con, Qe, Qs, h, 3*nxCells*nyCells);
	
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			A(i,j,0) = A_vect((i - ghost_cells) + nxCells * (j - ghost_cells));
			A(i,j,1) = A_vect((i - ghost_cells) + nxCells * (j - ghost_cells) + nxCells * nyCells);
			A(i,j,2) = A_vect((i - ghost_cells) + nxCells * (j - ghost_cells) + 2 * nxCells * nyCells);
		}
	}
	
	//	bound_cond_u(A, 3);
	bound_cond_A(A);
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			u_con(i,j,5) = 1 / (2.0 * dy) * (A(i,j + 1,2) - A(i,j - 1,2));
			u_con(i,j,6) = -1 / (2.0 * dx) * (A(i + 1,j,2) - A(i - 1,j,2));
			u_con(i,j,7) = 1 / (2.0 * dx) * (A(i + 1,j,1) - A(i - 1,j,1)) - 1 / (2.0 * dy) * (A(i,j + 1,0) - A(i,j - 1,0));
			
		}
	}
	bound_cond_u(u_con);
}

void PreSubsystem_Pressure(Grid &u_con, Grid &Qs, Grid &Qe,Grid &/*w_xL*/, Grid &/*w_xR*/,Grid &/*w_yL*/, Grid &/*w_yR*/, double h, const int rks){
	Array shock_eta(nxCells,nyCells);
	Grid Qe_pri(nxCells, nyCells);
	Grid Qs_pri(nxCells, nyCells);
	ucon2upri(Qe_pri, Qe);
	ucon2upri(Qs_pri, Qs);
	Eigen::VectorXd b_p_vect(nxCells*nyCells), b_Bt_vect(3 * nxCells * nyCells), B_vect(3 * nxCells * nyCells), p_vect(nxCells * nyCells);
	
	
	[[maybe_unused]] typedef Eigen::SparseMatrix<double> SpMat;
	[[maybe_unused]] typedef Eigen::Triplet<double> T;
	
	
	
	bound_cond_u(u_con);
	
	Array rhovx_star(nxCells, nyCells);
	Array rhovy_star(nxCells, nyCells);
	Array rhovz_star(nxCells, nyCells);
	
	for(int i = ghost_cells; i < ghost_cells + nxCells; i++){
		for(int j = ghost_cells; j < ghost_cells + nyCells; j++){
			rhovx_star(i,j) = u_con(i,j,1);
			rhovy_star(i,j) = u_con(i,j,2);
			rhovz_star(i,j) = u_con(i,j,3);
			
		}
	}
	bound_cond_array(rhovx_star);
	bound_cond_array(rhovy_star);
	bound_cond_array(rhovz_star);
	// jet inflow strip: 隐式 RHS 的 div(h_hat * rho v) 必须看到入流动量 rho_jet*v_jet,
	// 否则透射 copy 会把喷口处的隐式焓通量散度当成零梯度
	if(y_bc == y_JetInflow){
		for(int i = 0; i < nxCells + ghost_cells*2; i++){
			if(!jet_inlet(i)) continue;
			for(int g = 0; g < ghost_cells; g++){
				rhovx_star(i,g) = 0.0;
				rhovy_star(i,g) = jet_rho * jet_v;
				rhovz_star(i,g) = 0.0;
			}
		}
	}
	
	
//	std::cout << "rk = " <<rks + 1 <<", First Internal: "<<endl;
	vector<T> triplets_p;
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			double B2_00 = pow(u_con(i,j,5), 2.0) + pow(u_con(i,j,6), 2.0) + pow(u_con(i,j,7), 2.0);
			//			double B2_00 = pow(Qe(i,j,5), 2.0) + pow(Qe(i,j,6), 2.0) + pow(Qe(i,j,7), 2.0);
			double con00_n[NUM_VARS];
			double conm0_n[NUM_VARS];
			double conp0_n[NUM_VARS];
			double con0m_n[NUM_VARS];
			double con0p_n[NUM_VARS];
			
			
			for(int k = 0; k < NUM_VARS; k++){
				con00_n[k] = Qe(i,j,k);
				conm0_n[k] = Qe(i - 1,j,k);
				conp0_n[k] = Qe(i + 1,j,k);
				con0m_n[k] = Qe(i,j - 1,k);
				con0p_n[k] = Qe(i,j + 1,k);
			}
			
			
			double h00 = calenthalpy(con00_n) * Qe(i,j,0) / u_con(i,j,0);
			double hp0 = calenthalpy(conp0_n) * Qe(i + 1,j,0) / u_con(i + 1,j,0);
			double hm0 = calenthalpy(conm0_n) * Qe(i - 1,j,0) / u_con(i - 1,j,0);
			double h0p = calenthalpy(con0p_n) * Qe(i,j + 1,0) / u_con(i,j + 1,0);
			double h0m = calenthalpy(con0m_n) * Qe(i,j - 1,0) / u_con(i,j - 1,0);
			
			double p_star = u_con(i,j,4) - 0.5 * B2_00;
			double RHS;

			RHS = p_star - 0.5 * (rhovx_star(i,j) * rhovx_star(i,j) + rhovy_star(i,j) * rhovy_star(i,j) + rhovz_star(i,j) * rhovz_star(i,j)) / u_con(i,j,0) - h * (1 / dx * (0.5 * (hp0 + h00) * 0.5 * (rhovx_star(i + 1,j) + rhovx_star(i,j)) - 0.5 * (hm0 + h00) * 0.5 * (rhovx_star(i - 1,j) + rhovx_star(i,j))) + (1 / dy * (0.5 * (h0p + h00) * 0.5 * (rhovy_star(i,j + 1) + rhovy_star(i,j)) - 0.5 * (h0m + h00) * 0.5 * (rhovy_star(i,j - 1) + rhovy_star(i,j)))));
			
//			if(internal < 1e-6){
//				std::cout << internal << " ";
//			}

			b_p_vect(M_P(i, j)) = RHS;

			// 强 Dirichlet:边界 cell 的 RHS = p_amb(配算子恒等行 Ax=pii -> 解出 p=p_amb)
			if( ((x_bc==x_Dirichlet)&&(i==ghost_cells||i==nxCells+ghost_cells-1)) || ((y_bc==y_Dirichlet)&&(j==ghost_cells||j==nyCells+ghost_cells-1)) ) b_p_vect(M_P(i, j)) = p_amb;

			// jet inflow strip: ghost-cell Dirichlet lifting(配 Pressure_Subsystem_AX 的 p_0m=0)
			// 系数 = 算子中 p_{j-1} 的系数取负: +h^2/dy^2 * 0.5*(h00+h0m) * jet_p
			if(j == ghost_cells && jet_inlet(i)) b_p_vect(M_P(i, j)) += h*h/(dy*dy) * 0.5*(h00 + h0m) * jet_p;

			p_vect(M_P(i, j)) = Qe_pri(i,j,4);
		}
	}

//	std::cout << endl;


	Pressure_NestedNewton(p_vect, b_p_vect, u_con, Qe, Qs, h, nxCells*nyCells);
	
//	std::cout << "rk = " <<rks + 1<< ", First Solution: "<<std::endl;
//	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
//		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
//			if(p_vect(M_P(i,j))  < 1e-6){
//				std::cout << p_vect(M_P(i,j))  << " ";
//			}
//		}
//	}
//	std::cout << endl;
	
	for(int i = ghost_cells; i < ghost_cells + nxCells; i++){
		for(int j = ghost_cells; j < ghost_cells + nyCells; j++){
			double fL, fR,	flux_xm, flux_xp;
			fL = p_vect(M_P(i - 1,j));
			fR = p_vect(M_P(i,j));

			flux_xm = 0.5 * (fL + fR);

			fL = p_vect(M_P(i,j));
			fR = p_vect(M_P(i + 1,j));

			flux_xp = 0.5 * (fL + fR);

			u_con(i,j,1) = rhovx_star(i,j) - h / dx * (flux_xp - flux_xm);

			double gL, gR, flux_ym, flux_yp;
			// jet inflow strip: 底面压力通量用 Dirichlet ghost p=jet_p(M_P clamp 会退化成零梯度)
			gL = (j == ghost_cells && jet_inlet(i)) ? jet_p : p_vect(M_P(i,j - 1));
			gR = p_vect(M_P(i,j));

			flux_ym = 0.5 * (gL + gR);

			gL = p_vect(M_P(i,j));
			gR = p_vect(M_P(i,j + 1));

			flux_yp = 0.5 * (gL + gR);

			u_con(i,j,2) = rhovy_star(i,j) - h / dy * (flux_yp - flux_ym);
		}
	}
	
	
	bound_cond_u(u_con);
	
	
	
//	std::cout << "rk = " << rks + 1<< ", Second Internal: "<<endl;
	if(PressureSolver == all_correction || (PressureSolver == once_correction && rks == RKs - 1)){
		//std::cout <<"Second Pressure Solved "<< "rk: " << rks + 1 << endl;
		vector<T> triplets_p_second;
		for(int i = ghost_cells; i < ghost_cells + nxCells; i++){
			for(int j = ghost_cells; j < ghost_cells + nyCells; j++){
				double B2_00 = pow(u_con(i,j,5), 2.0) + pow(u_con(i,j,6), 2.0) + pow(u_con(i,j,7), 2.0);
				
				double con00_n[NUM_VARS];
				double conm0_n[NUM_VARS];
				double conp0_n[NUM_VARS];
				double con0m_n[NUM_VARS];
				double con0p_n[NUM_VARS];
				
				
				for(int k = 0; k < NUM_VARS; k++){
					con00_n[k] = Qe(i,j,k);
					conm0_n[k] = Qe(i - 1,j,k);
					conp0_n[k] = Qe(i + 1,j,k);
					con0m_n[k] = Qe(i,j - 1,k);
					con0p_n[k] = Qe(i,j + 1,k);
				}
				
				
				double h00 = calenthalpy(con00_n) * Qe(i,j,0) / u_con(i,j,0);
				double hp0 = calenthalpy(conp0_n) * Qe(i + 1,j,0) / u_con(i + 1,j,0);
				double hm0 = calenthalpy(conm0_n) * Qe(i - 1,j,0) / u_con(i - 1,j,0);
				double h0p = calenthalpy(con0p_n) * Qe(i,j + 1,0) / u_con(i,j + 1,0);
				double h0m = calenthalpy(con0m_n) * Qe(i,j - 1,0) / u_con(i,j - 1,0);
				
				double p_star = u_con(i,j,4) - 0.5 * B2_00;
				
				
				double RHS;
				
				RHS = p_star - 0.5 * (u_con(i,j,1) * u_con(i,j,1) + u_con(i,j,2) * u_con(i,j,2) + u_con(i,j,3) * u_con(i,j,3)) / u_con(i,j,0) - h * (1 / dx * (0.5 * (hp0 + h00) * 0.5 * (rhovx_star(i + 1,j) + rhovx_star(i,j)) - 0.5 * (hm0 + h00) * 0.5 * (rhovx_star(i - 1,j) + rhovx_star(i,j))) + (1 / dy * (0.5 * (h0p + h00) * 0.5 * (rhovy_star(i,j + 1) + rhovy_star(i,j)) - 0.5 * (h0m + h00) * 0.5 * (rhovy_star(i,j - 1) + rhovy_star(i,j)))));
				
//				if(internal < 1e-8){
//					std::cout << internal << " ";
//				}
				b_p_vect(M_P(i, j)) = RHS;
				// 强 Dirichlet:边界 cell RHS = p_amb(配算子恒等行)
				if( ((x_bc==x_Dirichlet)&&(i==ghost_cells||i==nxCells+ghost_cells-1)) || ((y_bc==y_Dirichlet)&&(j==ghost_cells||j==nyCells+ghost_cells-1)) ) b_p_vect(M_P(i, j)) = p_amb;

				// jet inflow strip: ghost-cell Dirichlet lifting(同第一次求解)
				if(j == ghost_cells && jet_inlet(i)) b_p_vect(M_P(i, j)) += h*h/(dy*dy) * 0.5*(h00 + h0m) * jet_p;

			}
		}
//		std::cout<<endl;

		Pressure_NestedNewton(p_vect, b_p_vect, u_con, Qe, Qs, h, nxCells*nyCells);
		
//		std::cout << "p Solution: " << std::endl;
//		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
//			for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
//				if(p_vect(M_P(i,j))  < 0){
//					std::cout << p_vect(M_P(i,j))  << " ";
//				}
//			}
//		}
//		std::cout << endl;
	}
//	std::cout << "rk = " <<rks + 1<< ", Internal Energy: "<< endl;
	switch(U_update){
		case DefinitionofU:
			for(int i = ghost_cells; i < ghost_cells + nxCells; i++){
				for(int j = ghost_cells; j < ghost_cells + nyCells; j++){
					u_con(i,j,4) = Vbase(p_vect(M_P(i,j)), u_con(i,j,0)) + 0.5 * (u_con(i,j,1) * u_con(i,j,1) + u_con(i,j,2) * u_con(i,j,2) + u_con(i,j,3) * u_con(i,j,3)) / u_con(i,j,0) + 0.5 * (u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6) + u_con(i,j,7) * u_con(i,j,7));
					
					
//					if(internal_energy < 0){
//						std::cout << "E_total=" << u_con(i,j,4) << " ";
//						std::cout << "KE=" << (0.5 * (u_con(i,j,1)*u_con(i,j,1) + u_con(i,j,2)*u_con(i,j,2) + u_con(i,j,3)*u_con(i,j,3))/u_con(i,j,0)) << " ";
//						std::cout << "ME=" << (0.5 * (u_con(i,j,5)*u_con(i,j,5) + u_con(i,j,6)*u_con(i,j,6) + u_con(i,j,7)*u_con(i,j,7))) << " ";
//						std::cout << "IE=" << internal_energy;
//					}
				}
			}
//			std::cout<<endl;
			break;
		case ConservationofU:
			Array final_pressure(nxCells, nyCells);
			for(int i = ghost_cells; i < ghost_cells + nxCells; i++){
				for(int j = ghost_cells; j < ghost_cells + nyCells; j++){
					final_pressure(i,j) = p_vect(M_P(i,j));
				}
			}
			bound_cond_array(final_pressure);
			// jet inflow strip: 底 ghost 压力 = jet_p(能量通量 h_0m 用)
			if(y_bc == y_JetInflow){
				for(int ii = 0; ii < nxCells + ghost_cells*2; ii++){
					if(!jet_inlet(ii)) continue;
					for(int gg = 0; gg < ghost_cells; gg++) final_pressure(ii,gg) = jet_p;
				}
			}
			// filled in below but never read -- see README, 'Known issues'
			[[maybe_unused]] double priL[NUM_VARS], priR[NUM_VARS];
			for(int i = ghost_cells; i < ghost_cells + nxCells; i++){
				for(int j = ghost_cells; j < ghost_cells + nyCells; j++){
					double fL, fR,	flux_xm, flux_xp;
					double gL, gR, flux_ym, flux_yp;
					
					double e_00 = eos_e_from_rho_p(u_con(i,j,0),     final_pressure(i,j));
					double e_p0 = eos_e_from_rho_p(u_con(i + 1,j,0), final_pressure(i + 1,j));
					double e_m0 = eos_e_from_rho_p(u_con(i - 1,j,0), final_pressure(i - 1,j));
					double e_0p = eos_e_from_rho_p(u_con(i,j + 1,0), final_pressure(i,j + 1));
					double e_0m = eos_e_from_rho_p(u_con(i,j - 1,0), final_pressure(i,j - 1));
					
					double h_00 = e_00 + final_pressure(i,j) / u_con(i,j,0);
					double h_p0 = e_p0 + final_pressure(i + 1,j) / u_con(i + 1,j,0);
					double h_m0 = e_m0 + final_pressure(i - 1,j) / u_con(i - 1,j,0);
					double h_0p = e_0p + final_pressure(i,j + 1) / u_con(i,j + 1,0);
					double h_0m = e_0m + final_pressure(i,j - 1) / u_con(i,j - 1,0);
					
					
					priL[0] = u_con(i - 1,j,0); priL[4] = final_pressure(i - 1,j);
					priL[1] = u_con(i - 1,j,1) / u_con(i - 1,j,0); priL[2] = u_con(i - 1,j,2) / u_con(i - 1,j,0); priL[3] = u_con(i - 1,j,3) / u_con(i - 1,j,0);
					priL[5] = u_con(i - 1,j,5); priL[6] = u_con(i - 1,j,6); priL[7] = u_con(i - 1,j,7);
					
					priR[0] = u_con(i,j,0); priR[4] = final_pressure(i,j);
					priR[1] = u_con(i,j,1) / u_con(i,j,0); priR[2] = u_con(i,j,2) / u_con(i,j,0); priR[3] = u_con(i,j,3) / u_con(i,j,0);
					priR[5] = u_con(i,j,5); priR[6] = u_con(i,j,6); priR[7] = u_con(i,j,7);
					fL = h_m0 * u_con(i-1,j,1);
					fR = h_00 * u_con(i,j,1);
					
					flux_xm = 0.5 * (fL + fR);
					
					priL[0] = u_con(i,j,0); priL[4] = final_pressure(i,j);
					priL[1] = u_con(i,j,1) / u_con(i,j,0); priL[2] = u_con(i,j,2) / u_con(i,j,0); priL[3] = u_con(i,j,3) / u_con(i,j,0);
					priL[5] = u_con(i,j,5); priL[6] = u_con(i,j,6); priL[7] = u_con(i,j,7);
					
					priR[0] = u_con(i + 1,j,0); priR[4] = final_pressure(i + 1,j);
					priR[1] = u_con(i + 1,j,1) / u_con(i + 1,j,0); priR[2] = u_con(i + 1,j,2) / u_con(i + 1,j,0); priR[3] = u_con(i + 1,j,3) / u_con(i + 1,j,0);
					priR[5] = u_con(i + 1,j,5); priR[6] = u_con(i + 1,j,6); priR[7] = u_con(i + 1,j,7);
					fL = h_00 * u_con(i,j,1);
					fR = h_p0 * u_con(i+1,j,1);
					
					flux_xp = 0.5 * (fL + fR);
					
					priL[0] = u_con(i,j - 1,0); priL[4] = final_pressure(i,j - 1);
					priL[1] = u_con(i,j - 1,1) / u_con(i,j - 1,0); priL[2] = u_con(i,j - 1,2) / u_con(i,j - 1,0); priL[3] = u_con(i,j - 1,3) / u_con(i,j - 1,0);
					priL[5] = u_con(i,j - 1,5); priL[6] = u_con(i,j - 1,6); priL[7] = u_con(i,j - 1,7);
					
					priR[0] = u_con(i,j,0); priR[4] = final_pressure(i,j);
					priR[1] = u_con(i,j,1) / u_con(i,j,0); priR[2] = u_con(i,j,2) / u_con(i,j,0); priR[3] = u_con(i,j,3) / u_con(i,j,0);
					priR[5] = u_con(i,j,5); priR[6] = u_con(i,j,6); priR[7] = u_con(i,j,7);
					
					gL = h_0m * u_con(i,j - 1,2);
					gR = h_00 * u_con(i,j,2);
					
					flux_ym = 0.5 * (gL + gR);
					
					priL[0] = u_con(i,j,0); priL[4] = final_pressure(i,j);
					priL[1] = u_con(i,j,1) / u_con(i,j,0); priL[2] = u_con(i,j,2) / u_con(i,j,0); priL[3] = u_con(i,j,3) / u_con(i,j,0);
					priL[5] = u_con(i,j,5); priL[6] = u_con(i,j,6); priL[7] = u_con(i,j,7);
					
					priR[0] = u_con(i,j + 1,0); priR[4] = final_pressure(i,j + 1);
					priR[1] = u_con(i,j + 1,1) / u_con(i,j + 1,0); priR[2] = u_con(i,j + 1,2) / u_con(i,j + 1,0); priR[3] = u_con(i,j + 1,3) / u_con(i,j + 1,0);
					priR[5] = u_con(i,j + 1,5); priR[6] = u_con(i,j + 1,6); priR[7] = u_con(i,j + 1,7);
					
					gL = h_00 * u_con(i,j,2);
					gR = h_0p * u_con(i,j + 1,2);
					
					flux_yp = 0.5 * (gL + gR);
					
					u_con(i,j,4) = u_con(i,j,4) - h / (dx) * (flux_xp - flux_xm) - h / (dy) * (flux_yp - flux_ym);
					
					
					double internal_energy = u_con(i,j,4) - (0.5 * (u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6) + u_con(i,j,7) * u_con(i,j,7))) - (0.5 * (u_con(i,j,1) * u_con(i,j,1) + u_con(i,j,2) * u_con(i,j,2) + u_con(i,j,3) * u_con(i,j,3)) / u_con(i,j,0) );
					
					if(internal_energy < 0){
						std::cout << "E_total=" << u_con(i,j,4) << " ";
						std::cout << "KE=" << (0.5 * (u_con(i,j,1)*u_con(i,j,1) + u_con(i,j,2)*u_con(i,j,2) + u_con(i,j,3)*u_con(i,j,3))/u_con(i,j,0)) << " ";
						std::cout << "ME=" << (0.5 * (u_con(i,j,5)*u_con(i,j,5) + u_con(i,j,6)*u_con(i,j,6) + u_con(i,j,7)*u_con(i,j,7))) << " ";
						std::cout << "IE=" << internal_energy;
					}
				}
			}
			std::cout<<endl;
			break;
	}
	
	bound_cond_u(u_con);
	
	
	b_p_vect.resize(0);
	b_Bt_vect.resize(0);
	B_vect.resize(0);
	p_vect.resize(0);
	bound_cond_u(u_con);
}

void PreSubsystem_Energy(Grid &u_con, Grid &Qs, Grid &Qe,Grid &/*w_xL*/, Grid &/*w_xR*/,Grid &/*w_yL*/, Grid &/*w_yR*/, double h, const int rks){
	Grid Qe_pri(nxCells, nyCells);
	ucon2upri(Qe_pri, Qe);
	double U_star;
	Eigen::VectorXd b_U_vect(nxCells*nyCells), U_vect(nxCells * nyCells);
	[[maybe_unused]] typedef Eigen::SparseMatrix<double> SpMat;
	[[maybe_unused]] typedef Eigen::Triplet<double> T;
	std::cout << "rk = " <<rks + 1 <<", RHS: "<<endl;
	vector<T> triplets_U;
	SpMat M_U(nxCells * nyCells, nxCells * nyCells);
	
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			double con00_n[NUM_VARS];
			double conm0_n[NUM_VARS];
			double conp0_n[NUM_VARS];
			double con0m_n[NUM_VARS];
			double con0p_n[NUM_VARS];
			
			for(int k = 0; k < NUM_VARS; k++){
				con00_n[k] = Qe(i,j,k);
				conm0_n[k] = Qe(i - 1,j,k);
				conp0_n[k] = Qe(i + 1,j,k);
				con0m_n[k] = Qe(i,j - 1,k);
				con0p_n[k] = Qe(i,j + 1,k);
			}
			
			double h_00 = caltotalenthalpy(con00_n) *  con00_n[0]/ u_con(i,j,0);
			double h_p0 = caltotalenthalpy(conp0_n) *  conp0_n[0]/ u_con(i+1,j,0);
			double h_m0 = caltotalenthalpy(conm0_n) *  conm0_n[0]/ u_con(i-1,j,0);
			double h_0p = caltotalenthalpy(con0p_n) *  con0p_n[0]/ u_con(i,j+1,0);
			double h_0m = caltotalenthalpy(con0m_n) *  con0m_n[0]/ u_con(i,j-1,0);
			
			U_star = u_con(i,j,4) - 0.5 * h / dx * (h_p0 * u_con(i + 1,j,1) - h_m0 * u_con(i - 1,j,1)) - 0.5 * h / dy * (h_0p * u_con(i,j + 1,2) - h_0m * u_con(i,j - 1,2));


			b_U_vect(M_P(i, j)) = U_star * dx * dy;

			// Dirichlet lifting: subtract operator applied to (interior=0, ghost=E_bc)
			if(x_bc == x_Dirichlet || y_bc == y_Dirichlet || y_bc == y_JetInflow){
				double hhp_x = 0.5*(h_00+h_p0), hhm_x = 0.5*(h_00+h_m0);
				double hhp_y = 0.5*(h_00+h_0p), hhm_y = 0.5*(h_00+h_0m);
				// ghost-cell Dirichlet 1×(配 getE ghost=0);钉内能=p_amb -> E = p_amb/(g-1) + 1/2 B_ghost^2
				// (用透射 B,使内能恒=ambient,B 偏离 ambient 时也一致)
				// jet strip: E_ghost = jet_E_val()(含动能 0.5*rho*v^2,不可省!)
				auto getEc = [&](int ii, int jj) -> double {
					bool xg = (ii < ghost_cells) || (ii >= nxCells + ghost_cells);
					bool yg = (jj < ghost_cells) || (jj >= nyCells + ghost_cells);
					if(jj < ghost_cells && jet_inlet(ii)) return jet_E_val();   // 与 getE 的 strip 分支严格配对
					bool g = (x_bc == x_Dirichlet && xg) || (y_bc == y_Dirichlet && yg);
					if(!g) return 0.0;
					return p_amb/(gamma_plasma-1.0)
					     + 0.5*(u_con(ii,jj,5)*u_con(ii,jj,5)+u_con(ii,jj,6)*u_con(ii,jj,6)+u_con(ii,jj,7)*u_con(ii,jj,7));
				};
				double L_const = - (gamma_plasma - 1) * h * h * (
					1.0/(dx*dx) * ( hhp_x * getEc(i+1,j) + hhm_x * getEc(i-1,j) )
					+ 1.0/(dy*dy) * ( hhp_y * getEc(i,j+1) + hhm_y * getEc(i,j-1) )
				);
				L_const *= dx * dy;
				b_U_vect(M_P(i, j)) -= L_const;
			}
			U_vect(M_P(i,j)) = u_con(i,j,4); // 初始猜测
			
			if(U_star < 1e-6){
				std::cout << U_star << " ";
			}
		}
		
	}
	std::cout << endl;
	// Diagonal (Jacobi) preconditioner for the energy solve (matches CASTOR precopDiag_SIMHD_energy)
	auto precE = [&](Eigen::VectorXd& out, const Eigen::VectorXd& in){
		out.resize(in.size());
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
				double c00[NUM_VARS], cp0[NUM_VARS], cm0[NUM_VARS], c0p[NUM_VARS], c0m[NUM_VARS];
				for(int k = 0; k < NUM_VARS; k++){
					c00[k]=Qe(i,j,k); cp0[k]=Qe(i+1,j,k); cm0[k]=Qe(i-1,j,k); c0p[k]=Qe(i,j+1,k); c0m[k]=Qe(i,j-1,k);
				}
				double h00 = caltotalenthalpy(c00)*c00[0]/u_con(i,j,0);
				double hp0 = caltotalenthalpy(cp0)*cp0[0]/u_con(i+1,j,0);
				double hm0 = caltotalenthalpy(cm0)*cm0[0]/u_con(i-1,j,0);
				double h0p = caltotalenthalpy(c0p)*c0p[0]/u_con(i,j+1,0);
				double h0m = caltotalenthalpy(c0m)*c0m[0]/u_con(i,j-1,0);
				double hhp_x=0.5*(h00+hp0), hhm_x=0.5*(h00+hm0);
				double hhp_y=0.5*(h00+h0p), hhm_y=0.5*(h00+h0m);
				double diag = 1.0 + (gamma_plasma-1.0)*h*h*( (hhp_x+hhm_x)/(dx*dx) + (hhp_y+hhm_y)/(dy*dy) );
				out(M_P(i,j)) = in(M_P(i,j)) / (diag * dx * dy);
			}
		}
	};
	GMRES_Solver(U_vect, b_U_vect, Energy_Subsystem_AX, u_con, Qe, Qs, h, nxCells*nyCells, precE);

	for(int i = ghost_cells; i < ghost_cells + nxCells; i++){
		for(int j = ghost_cells; j < ghost_cells + nyCells; j++){
			u_con(i,j,4) = U_vect(M_P(i, j));
		}
	}
	bound_cond_u(u_con);
	for(int i = ghost_cells; i < ghost_cells + nxCells; i++){
		for(int j = ghost_cells; j < ghost_cells + nyCells; j++){
			u_con(i,j,1) -= h * (1/(2*dx)) * (gamma_plasma-1) * (u_con(i+1,j,4) - u_con(i-1,j,4));
			u_con(i,j,2) -= h * (1/(2*dy)) * (gamma_plasma-1) * (u_con(i,j+1,4) - u_con(i,j-1,4));
		}
	}
	
	bound_cond_u(u_con);
	
	std::cout << "rk = " <<rks + 1 <<", Internal Energy: "<<endl;
	
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			double internal = u_con(i,j,4) - 0.5*(u_con(i,j,1)*u_con(i,j,1) + u_con(i,j,2)*u_con(i,j,2) + u_con(i,j,3)*u_con(i,j,3))/u_con(i,j,0) -  0.5*(u_con(i,j,5)*u_con(i,j,5) + u_con(i,j,6)*u_con(i,j,6) + u_con(i,j,7)*u_con(i,j,7));
			if(internal < 1e-6){
				std::cout << internal << " ";
			}
		}
	}
	std::cout << endl;
}


void MagSubsystem_UpdateB(Grid &u_con, Grid &Qs, Grid &Qe,Grid &/*w_xL*/, Grid &/*w_xR*/,Grid &/*w_yL*/, Grid &/*w_yR*/, Grid &/*A*/,double h, const int /*rks*/){
	double Bx_star, By_star, Bz_star;
	Grid Qe_pri(nxCells, nyCells);
	Grid Qs_pri(nxCells, nyCells);
	ucon2upri(Qe_pri, Qe);
	ucon2upri(Qs_pri, Qs);
	Eigen::VectorXd b_U_vect(nxCells*nyCells), b_Bt_vect(3 * nxCells * nyCells), B_vect(3 * nxCells * nyCells), U_vect(nxCells * nyCells);
	
	
	[[maybe_unused]] typedef Eigen::SparseMatrix<double> SpMat;
	[[maybe_unused]] typedef Eigen::Triplet<double> T;
	
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			double rho_p0 = 0.5 * (u_con(i,j,0) + u_con(i + 1,j,0));
			double rho_m0 = 0.5 * (u_con(i,j,0) + u_con(i - 1,j,0));
			double rho_0p = 0.5 * (u_con(i,j,0) + u_con(i,j + 1,0));
			double rho_0m = 0.5 * (u_con(i,j,0) + u_con(i,j - 1,0));
			
			double rhovxstar_p0 = 0.5 * (u_con(i,j,1) + u_con(i + 1,j,1));
			double rhovxstar_m0 = 0.5 * (u_con(i,j,1) + u_con(i - 1,j,1));
			double rhovxstar_0p = 0.5 * (u_con(i,j,1) + u_con(i,j + 1,1));
			double rhovxstar_0m = 0.5 * (u_con(i,j,1) + u_con(i,j - 1,1));
			
			double rhovystar_p0 = 0.5 * (u_con(i,j,2) + u_con(i + 1,j,2));
			double rhovystar_m0 = 0.5 * (u_con(i,j,2) + u_con(i - 1,j,2));
			double rhovystar_0p = 0.5 * (u_con(i,j,2) + u_con(i,j + 1,2));
			double rhovystar_0m = 0.5 * (u_con(i,j,2) + u_con(i,j - 1,2));
			
			double rhovzstar_p0 = 0.5 * (u_con(i,j,3) + u_con(i + 1,j,3));
			double rhovzstar_m0 = 0.5 * (u_con(i,j,3) + u_con(i - 1,j,3));
			double rhovzstar_0p = 0.5 * (u_con(i,j,3) + u_con(i,j + 1,3));
			double rhovzstar_0m = 0.5 * (u_con(i,j,3) + u_con(i,j - 1,3));
			
			
			double p_pp = 0.25 * (Qe_pri(i,j,4) + Qe_pri(i + 1,j,4) + Qe_pri(i,j + 1,4) + Qe_pri(i + 1,j + 1,4));
			double p_pm = 0.25 * (Qe_pri(i,j,4) + Qe_pri(i + 1,j,4) + Qe_pri(i,j - 1,4) + Qe_pri(i + 1,j - 1,4));
			double p_mp = 0.25 * (Qe_pri(i,j,4) + Qe_pri(i - 1,j,4) + Qe_pri(i,j + 1,4) + Qe_pri(i - 1,j + 1,4));
			double p_mm = 0.25 * (Qe_pri(i,j,4) + Qe_pri(i - 1,j,4) + Qe_pri(i,j - 1,4) + Qe_pri(i - 1,j - 1,4));
			
			double Bxn_p0 = 0.5 * (Qe(i,j,5) + Qe(i + 1,j,5));
			double Bxn_m0 = 0.5 * (Qe(i,j,5) + Qe(i - 1,j,5));
			double Bxn_0p = 0.5 * (Qe(i,j,5) + Qe(i,j + 1,5));
			double Bxn_0m = 0.5 * (Qe(i,j,5) + Qe(i,j - 1,5));
			
			double Byn_p0 = 0.5 * (Qe(i,j,6) + Qe(i + 1,j,6));
			double Byn_m0 = 0.5 * (Qe(i,j,6) + Qe(i - 1,j,6));
			double Byn_0p = 0.5 * (Qe(i,j,6) + Qe(i,j + 1,6));
			double Byn_0m = 0.5 * (Qe(i,j,6) + Qe(i,j - 1,6));
			
			double Bzn_p0 = 0.5 * (Qe(i,j,7) + Qe(i + 1,j,7));
			double Bzn_m0 = 0.5 * (Qe(i,j,7) + Qe(i - 1,j,7));
			double Bzn_0p = 0.5 * (Qe(i,j,7) + Qe(i,j + 1,7));
			double Bzn_0m = 0.5 * (Qe(i,j,7) + Qe(i,j - 1,7));
			
			
			
			Bx_star = Qs(i,j,5) - (h / dy) * ((Bxn_0p / rho_0p * (rhovystar_0p - (h / dy) * (Qe_pri(i,j + 1,4) - Qe_pri(i,j,4))) - Byn_0p / rho_0p * (rhovxstar_0p - h / dx * (p_pp - p_mp)  )) - (Bxn_0m / rho_0m * (rhovystar_0m - (h / dy) * (Qe_pri(i,j,4) - Qe_pri(i,j - 1,4))) - Byn_0m / rho_0m * (rhovxstar_0m - h / dx * (p_pm - p_mm)  )) ); 
			
			By_star = Qs(i,j,6) - (h / dx) * ((Byn_p0 / rho_p0 * (rhovxstar_p0 - (h / dx) * (Qe_pri(i + 1,j,4) - Qe_pri(i,j,4))) - Bxn_p0 / rho_p0 * (rhovystar_p0 - h / dy * (p_pp - p_pm) ) ) - (Byn_m0/rho_m0 * (rhovxstar_m0 - (h / dx) * (Qe_pri(i,j,4) - Qe_pri(i - 1,j,4))) - Bxn_m0 / rho_m0 * (rhovystar_m0 - h / dy * (p_mp - p_mm)  ) ));
			
			double Bz_star_1 = ((Bzn_p0 / rho_p0 * (rhovxstar_p0 - h / dx * (Qe_pri(i + 1,j,4) - Qe_pri(i,j,4))) - Bxn_p0 / rho_p0 * (rhovzstar_p0) ) - (Bzn_m0 / rho_m0 * (rhovxstar_m0 - h / dx * (Qe_pri(i,j,4) - Qe_pri(i - 1,j,4))) - Bxn_m0 / rho_m0 * (rhovzstar_m0)) );
			
			double Bz_star_2 = ((Bzn_0p / rho_0p * (rhovystar_0p - h / dy * (Qe_pri(i,j + 1,4) - Qe_pri(i,j,4))) - Byn_0p / rho_0p * (rhovzstar_0p) ) - (Bzn_0m / rho_0m * (rhovystar_0m - h / dy * (Qe_pri(i,j,4) - Qe_pri(i,j - 1,4))) - Byn_0m / rho_0m * (rhovzstar_0m)) );
			
			Bz_star = Qs(i,j,7) - (h / dx) * Bz_star_1 - (h / dy) * Bz_star_2;
			
			
			//Laplacian Operator
			b_Bt_vect((i - ghost_cells) + nxCells * (j - ghost_cells)) =  (Bx_star);
			b_Bt_vect((i - ghost_cells) + nxCells * (j - ghost_cells) + nxCells * nyCells) = (By_star);
			b_Bt_vect((i - ghost_cells) + nxCells * (j - ghost_cells) + 2 * nxCells * nyCells) = (Bz_star);
			
			B_vect((i - ghost_cells) + nxCells * (j - ghost_cells)) = Qe(i,j,5);
			B_vect((i - ghost_cells) + nxCells * (j - ghost_cells) + nxCells * nyCells) = Qe(i,j,6);
			B_vect((i - ghost_cells) + nxCells * (j - ghost_cells) + 2 * nxCells * nyCells) = Qe(i,j,7);
			
		}
	}
	GMRES_Solver(B_vect, b_Bt_vect, Magnetic_Subsystem_AX, u_con, Qe, Qs, h, 3*nxCells*nyCells);
	
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			u_con(i,j,5) = B_vect((i - ghost_cells) + nxCells * (j - ghost_cells));
			u_con(i,j,6) = B_vect((i - ghost_cells) + nxCells * (j - ghost_cells) + nxCells * nyCells);
			u_con(i,j,7) = B_vect((i - ghost_cells) + nxCells * (j - ghost_cells) + 2 * nxCells * nyCells);
		}
	}
	bound_cond_u(u_con);
	
	CollectMagFluxes(u_con, Qe, Qs, h);
}

double cal_divergence_error(Grid &u_con, Array &Bx_0p, Array &By_p0){
	// 每种 div-free 方案测它自己保证为零的散度定义,且测当前真实 B(不是旧的交错场):
	//   StaggerCT                        -> 交错面场 (Bx_p0,By_0p) 的有限体积散度
	//   PotentialA / UnStaggerCT / NoDivFree -> cell-centered B = curl A 的散度,用 CASTOR
	//        的 ComputeCellGradientA 同款梯度算子:内部中心差分,域边界线性外推=单侧差分。
	// 注意:边界**不能**用中心差分去碰 bound_cond_u 的透射 ghost B(≠curl A),否则
	//       边界 div(curl)≠0 且随时间增长(透射 ghost 与 curl A 不一致)。CASTOR 用同一个
	//       外推梯度算子构造 B 和测 divB,故边界也保持一致。
	// (旧版固定用 Bx_p0/By_0p,对后两种方案是没更新的初始值 -> 假的 0)
	double error = 0;
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			double d;
			if(magfield == StaggerCT){
				d = (Bx_0p(i,j) - Bx_0p(i-1,j)) / dx + (By_p0(i,j) - By_p0(i,j-1)) / dy;
			}else{
				// CASTOR ComputeCellGradientA: 内部中心差分;域边界线性外推 -> 单侧差分
				double dBxdx, dBydy;
				if(i == ghost_cells)                    dBxdx = (u_con(i+1,j,5) - u_con(i,j,5)) / dx;
				else if(i == nxCells + ghost_cells - 1) dBxdx = (u_con(i,j,5) - u_con(i-1,j,5)) / dx;
				else                                    dBxdx = (u_con(i+1,j,5) - u_con(i-1,j,5)) / (2.0*dx);
				if(j == ghost_cells)                    dBydy = (u_con(i,j+1,6) - u_con(i,j,6)) / dy;
				else if(j == nyCells + ghost_cells - 1) dBydy = (u_con(i,j,6) - u_con(i,j-1,6)) / dy;
				else                                    dBydy = (u_con(i,j+1,6) - u_con(i,j-1,6)) / (2.0*dy);
				d = dBxdx + dBydy;
			}
			error += dx * dy * d * d;
		}
	}
	return sqrt(error);
}
