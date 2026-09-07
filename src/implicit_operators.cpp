// Matrix-free operators for the implicit subsystems, nested Newton, GMRES.

#include "declarations.hpp"

int M_P(int i, int j){
	int coord = 0;
	int I = i, J = j;   // the switches below cover every enumerator and overwrite these
	
	if(i < ghost_cells){
		switch(x_bc){
			case x_Transmissive: case x_Dirichlet:
				I = ghost_cells;
				break;
			case x_Reflective:
				I = 2 * ghost_cells - 1 - i;
				break;
			case x_Periodic:
				I = nxCells + i;
				break;
		}
	}else if(i > nxCells + ghost_cells - 1){
		switch(x_bc){
			case x_Transmissive: case x_Dirichlet:
				I = nxCells + ghost_cells - 1;
				break;
			case x_Reflective:
				I = 2 * (nxCells + ghost_cells) - i - 1;
				break;
			case x_Periodic:
				I = i - nxCells;
				break;
		}
	}else{
		I = i;
	}
	
	if(j < ghost_cells){
		switch(y_bc){
			case y_Transmissive: case y_Dirichlet: case y_JetInflow:
				J = ghost_cells;
				break;
			case y_Reflective:
				J = 2 * ghost_cells - 1 - j;
				break;
			case y_Periodic:
				J = nyCells + j;
				break;
		}
	}else if(j > nyCells + ghost_cells - 1){
		switch(y_bc){
			case y_Transmissive: case y_Dirichlet: case y_JetInflow:
				J = nyCells + ghost_cells - 1;
				break;
			case y_Reflective:
				J = 2 * (nyCells + ghost_cells) - j - 1;
				break;
			case y_Periodic:
				J = j - nyCells;
				break;
		}
	}else{
		J = j;
	}
	coord = (I - ghost_cells) + nxCells * (J - ghost_cells);
	return coord;
}

int M_B(int i, int j, int k){
	int coord = 0;
	int I = i, J = j;   // the switches below cover every enumerator and overwrite these
	
	if(i < ghost_cells){
		switch(x_bc){
			case x_Transmissive: case x_Dirichlet:
				I = ghost_cells;
				break;
			case x_Reflective:
				I = 2 * ghost_cells - 1 - i;
				break;
			case x_Periodic:
				I = nxCells + i;
				break;
		}
	}else if(i > nxCells + ghost_cells - 1){
		switch(x_bc){
			case x_Transmissive: case x_Dirichlet:
				I = nxCells + ghost_cells - 1;
				break;
			case x_Reflective:
				I = 2 * (nxCells + ghost_cells) - i - 1;
				break;
			case x_Periodic:
				I = i - nxCells;
				break;
		}
	}else{
		I = i;
	}
	
	if(j < ghost_cells){
		switch(y_bc){
			case y_Transmissive: case y_Dirichlet: case y_JetInflow:
				J = ghost_cells;
				break;
			case y_Reflective:
				J = 2 * ghost_cells - 1 - j;
				break;
			case y_Periodic:
				J = nyCells + j;
				break;
		}
	}else if(j > nyCells + ghost_cells - 1){
		switch(y_bc){
			case y_Transmissive: case y_Dirichlet: case y_JetInflow:
				J = nyCells + ghost_cells - 1;
				break;
			case y_Reflective:
				J = 2 * (nyCells + ghost_cells) - j - 1;
				break;
			case y_Periodic:
				J = j - nyCells;
				break;
		}
	}else{
		J = j;
	}
	
	if(k == 5){
		coord = (I - ghost_cells) + nxCells * (J - ghost_cells);
	}else if(k == 6){
		coord = (I - ghost_cells) + nxCells * (J - ghost_cells) + nxCells * nyCells;
	}else if(k == 7){
		coord = (I - ghost_cells) + nxCells * (J - ghost_cells) + 2 * nxCells * nyCells;
	}else{
		std::cout << "Error" << endl;
	}
	
	return coord;
}

void Magnetic_Subsystem_AX(Eigen::VectorXd& Ax, Eigen::VectorXd& x, Grid& u_con, Grid& Qe,Grid& /*Qs*/,double h){
	int reflct_xL, reflct_yL, reflct_xR, reflct_yR;
	Grid Qe_pri(nxCells, nyCells);
	ucon2upri(Qe_pri, Qe);
	
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			reflct_xL = 1, reflct_yL = 1,reflct_xR = 1,reflct_yR = 1;
			int xL, xR, yL, yR;
			switch(x_bc){
				case x_Transmissive: case x_Dirichlet:
					xL = max(ghost_cells,i-1);
					xR = min(nxCells + ghost_cells - 1,i+1);
					break;
				case x_Reflective:
					xL = max(ghost_cells,i-1);
					xR = min(nxCells + ghost_cells - 1,i+1);
					
					if(xL == ghost_cells){
						reflct_xL = -1;
					}
					
					if(xR == nxCells + ghost_cells - 1){
						reflct_xR = -1;
					}
					break;
				case x_Periodic:
					if(i == ghost_cells){
						xL = nxCells + ghost_cells - 1;
					}else{
						xL = i - 1;
					}
					
					if(i == nxCells + ghost_cells - 1){
						xR = ghost_cells;
					}else{
						xR = i + 1;
					}
					break;
			}
			
			switch(y_bc){
				case y_Transmissive: case y_Dirichlet: case y_JetInflow:
					yL = max(ghost_cells,j-1);
					yR = min(nyCells + ghost_cells - 1,j+1);
					break;
				case y_Reflective:
					yL = max(ghost_cells,j-1);
					yR = min(nyCells + ghost_cells - 1,j+1);
					
					if(yL == ghost_cells){
						reflct_yL = -1;
					}
					
					if(yR == nyCells + ghost_cells - 1){
						reflct_yR = -1;
					}
					break;
				case y_Periodic:
					if(j == ghost_cells){
						yL = nyCells + ghost_cells - 1;
					}else{
						yL = j - 1;
					}
					
					if(j == nyCells + ghost_cells - 1){
						yR = ghost_cells;
					}else{
						yR = j + 1;
					}
					break;
			}
			
			double rho_p0 = 0.5 * (u_con(i,j,0) + u_con(i + 1,j,0));
			double rho_m0 = 0.5 * (u_con(i,j,0) + u_con(i - 1,j,0));
			double rho_0p = 0.5 * (u_con(i,j,0) + u_con(i,j + 1,0));
			double rho_0m = 0.5 * (u_con(i,j,0) + u_con(i,j - 1,0));
			
			double Bxn_p0 = 0.5 * (Qe(i,j,5) + Qe(i + 1,j,5));
			double Bxn_m0 = 0.5 * (Qe(i,j,5) + Qe(i - 1,j,5));
			double Bxn_0p = 0.5 * (Qe(i,j,5) + Qe(i,j + 1,5));
			double Bxn_0m = 0.5 * (Qe(i,j,5) + Qe(i,j - 1,5));
			double Bxn_pp = 0.25 * (Qe(i,j,5) + Qe(i + 1,j,5) + Qe(i,j + 1,5) + Qe(i + 1,j + 1,5));
			double Bxn_pm = 0.25 * (Qe(i,j,5) + Qe(i + 1,j,5) + Qe(i,j - 1,5) + Qe(i + 1,j - 1,5));
			double Bxn_mp = 0.25 * (Qe(i,j,5) + Qe(i - 1,j,5) + Qe(i,j + 1,5) + Qe(i - 1,j + 1,5));
			double Bxn_mm = 0.25 * (Qe(i,j,5) + Qe(i - 1,j,5) + Qe(i,j - 1,5) + Qe(i - 1,j - 1,5));
			
			double Byn_p0 = 0.5 * (Qe(i,j,6) + Qe(i + 1,j,6));
			double Byn_m0 = 0.5 * (Qe(i,j,6) + Qe(i - 1,j,6));
			double Byn_0p = 0.5 * (Qe(i,j,6) + Qe(i,j + 1,6));
			double Byn_0m = 0.5 * (Qe(i,j,6) + Qe(i,j - 1,6));
			double Byn_pp = 0.25 * (Qe(i,j,6) + Qe(i + 1,j,6) + Qe(i,j + 1,6) + Qe(i + 1,j + 1,6));
			double Byn_pm = 0.25 * (Qe(i,j,6) + Qe(i + 1,j,6) + Qe(i,j - 1,6) + Qe(i + 1,j - 1,6));
			double Byn_mp = 0.25 * (Qe(i,j,6) + Qe(i - 1,j,6) + Qe(i,j + 1,6) + Qe(i - 1,j + 1,6));
			double Byn_mm = 0.25 * (Qe(i,j,6) + Qe(i - 1,j,6) + Qe(i,j - 1,6) + Qe(i - 1,j - 1,6));
			
			double Bzn_p0 = 0.5 * (Qe(i,j,7) + Qe(i + 1,j,7));
			double Bzn_m0 = 0.5 * (Qe(i,j,7) + Qe(i - 1,j,7));
			double Bzn_0p = 0.5 * (Qe(i,j,7) + Qe(i,j + 1,7));
			double Bzn_0m = 0.5 * (Qe(i,j,7) + Qe(i,j - 1,7));
			double Bzn_pp = 0.25 * (Qe(i,j,7) + Qe(i + 1,j,7) + Qe(i,j + 1,7) + Qe(i + 1,j + 1,7));
			double Bzn_pm = 0.25 * (Qe(i,j,7) + Qe(i + 1,j,7) + Qe(i,j - 1,7) + Qe(i + 1,j - 1,7));
			double Bzn_mp = 0.25 * (Qe(i,j,7) + Qe(i - 1,j,7) + Qe(i,j + 1,7) + Qe(i - 1,j + 1,7));
			double Bzn_mm = 0.25 * (Qe(i,j,7) + Qe(i - 1,j,7) + Qe(i,j - 1,7) + Qe(i - 1,j - 1,7));
			
			double Bx_pp = 0.25 * (x(M_B(i,j,5)) + x(M_B(i+1,j,5)) * reflct_xR + x(M_B(i,j+1,5)) * reflct_yR + x(M_B(i+1,j+1,5)) * reflct_xR * reflct_yR);
			double Bx_pm = 0.25 * (x(M_B(i,j,5)) + x(M_B(i+1,j,5)) * reflct_xR + x(M_B(i,j-1,5)) * reflct_yL + x(M_B(i+1,j-1,5)) * reflct_xR * reflct_yL);
			double Bx_mp = 0.25 * (x(M_B(i,j,5)) + x(M_B(i-1,j,5)) * reflct_xL + x(M_B(i,j+1,5)) * reflct_yR + x(M_B(i-1,j+1,5)) * reflct_xL * reflct_yR);
			double Bx_mm = 0.25 * (x(M_B(i,j,5)) + x(M_B(i-1,j,5)) * reflct_xL + x(M_B(i,j-1,5)) * reflct_yL + x(M_B(i-1,j-1,5)) * reflct_xL * reflct_yL);
			
			double By_pp = 0.25 * (x(M_B(i,j,6)) + x(M_B(i+1,j,6)) * reflct_xR + x(M_B(i,j+1,6)) * reflct_yR + x(M_B(i+1,j+1,6)) * reflct_xR * reflct_yR);
			double By_pm = 0.25 * (x(M_B(i,j,6)) + x(M_B(i+1,j,6)) * reflct_xR + x(M_B(i,j-1,6)) * reflct_yL + x(M_B(i+1,j-1,6)) * reflct_xR * reflct_yL);
			double By_mp = 0.25 * (x(M_B(i,j,6)) + x(M_B(i-1,j,6)) * reflct_xL + x(M_B(i,j+1,6)) * reflct_yR + x(M_B(i-1,j+1,6)) * reflct_xL * reflct_yR);
			double By_mm = 0.25 * (x(M_B(i,j,6)) + x(M_B(i-1,j,6)) * reflct_xL + x(M_B(i,j-1,6)) * reflct_yL + x(M_B(i-1,j-1,6)) * reflct_xL * reflct_yL);
			
			double Bz_pp = 0.25 * (x(M_B(i,j,7)) + x(M_B(i+1,j,7)) * reflct_xR + x(M_B(i,j+1,7)) * reflct_yR + x(M_B(i+1,j+1,7)) * reflct_xR * reflct_yR);
			double Bz_pm = 0.25 * (x(M_B(i,j,7)) + x(M_B(i+1,j,7)) * reflct_xR + x(M_B(i,j-1,7)) * reflct_yL + x(M_B(i+1,j-1,7)) * reflct_xR * reflct_yL);
			double Bz_mp = 0.25 * (x(M_B(i,j,7)) + x(M_B(i-1,j,7)) * reflct_xL + x(M_B(i,j+1,7)) * reflct_yR + x(M_B(i-1,j+1,7)) * reflct_xL * reflct_yR);
			double Bz_mm = 0.25 * (x(M_B(i,j,7)) + x(M_B(i-1,j,7)) * reflct_xL + x(M_B(i,j-1,7)) * reflct_yL + x(M_B(i-1,j-1,7)) * reflct_xL * reflct_yL);
			
			//Laplacian Operator
//			double nu_x = B_max_x * dx;
//			double nu_y = B_max_y * dy;
			
			
			double Bx_1 = Bxn_0p / rho_0p * (1 / dx * (Bxn_pp * By_pp - Bxn_mp * By_mp) - 1 / dy * ((0.5 * (Qe(i, j + 1, 5) * x(M_B(i,j + 1,5)) * reflct_yR + Qe(i, j + 1, 6) * x(M_B(i,j + 1,6)) * reflct_yR + Qe(i, j + 1, 7) * x(M_B(i,j + 1,7)) * reflct_yR ) - Qe(i, j + 1, 6) * x(M_B(i,j + 1,6)) * reflct_yR) - (0.5 * (Qe(i, j, 5) * x(M_B(i,j,5)) + Qe(i, j, 6) * x(M_B(i,j,6)) + Qe(i, j, 7) * x(M_B(i,j,7))) - Qe(i, j, 6) * x(M_B(i,j,6))) ) );
			
			double Bx_2 = Byn_0p / rho_0p * (1 / dx * ((0.5 * (Bxn_pp * Bx_pp + Byn_pp * By_pp + Bzn_pp * Bz_pp) - Bxn_pp * Bx_pp) - (0.5 * (Bxn_mp * Bx_mp + Byn_mp * By_mp + Bzn_mp * Bz_mp) - Bxn_mp * Bx_mp )) - 1 / dy * (Qe(i,j + 1,6) * x(M_B(i,j + 1,5)) * reflct_yR - Qe(i,j,6) * x(M_B(i,j,5)))  );
			
			double Bx_3 = Bxn_0m / rho_0m * (1 / dx * (Bxn_pm * By_pm - Bxn_mm * By_mm) - 1 / dy * ((0.5 *(Qe(i, j, 5) * x(M_B(i,j,5)) + Qe(i, j, 6) * x(M_B(i,j,6)) + Qe(i, j, 7) * x(M_B(i,j,7))) - Qe(i, j, 6) * x(M_B(i, j, 6))) - (0.5 * (Qe(i, j - 1, 5) * x(M_B(i,j - 1,5)) * reflct_yL + Qe(i, j - 1, 6) * x(M_B(i,j - 1,6))  * reflct_yL + Qe(i, j - 1, 7) * x(M_B(i,j - 1,7))  * reflct_yL ) - Qe(i, j - 1, 6) * x(M_B(i,j - 1,6))  * reflct_yL) ) );
			
			double Bx_4 =  Byn_0m / rho_0m * (1 / dx * ((0.5 * (Bxn_pm * Bx_pm + Byn_pm * By_pm + Bzn_pm * Bz_pm) - Bxn_pm * Bx_pm) - (0.5 * (Bxn_mm * Bx_mm + Byn_mm * By_mm + Bzn_mm * Bz_mm) - Bxn_mm * Bx_mm)) - 1 / dy * (Qe(i,j,6) * x(M_B(i,j,5)) - Qe(i,j - 1,6) * x(M_B(i,j - 1,5)) * reflct_yL)  );
			
			
			Ax(M_B(i,j,5)) = x(M_B(i,j,5)) + (h * h/ dy) * ((Bx_1 + Bx_2) - (Bx_3 + Bx_4));

			double By_1 = Byn_p0 / rho_p0 * (-1 / dx * ((0.5 * ((Qe(i + 1, j, 5) * x(M_B(i + 1,j,5)) * reflct_xR + Qe(i + 1, j, 6) * x(M_B(i+1,j,6)) * reflct_xR + Qe(i + 1, j, 7) * x(M_B(i + 1, j, 7)) * reflct_xR )) - Qe(i + 1, j, 5) * x(M_B(i + 1, j, 5)) * reflct_xR ) - (0.5 * ((Qe(i, j, 5) * x(M_B(i, j, 5)) + Qe(i, j, 6) * x(M_B(i, j, 6)) + Qe(i, j, 7) * x(M_B(i, j, 7)))) - Qe(i, j, 5) * x(M_B(i, j, 5)) )) + 1 / dy * (Bxn_pp * By_pp - Bxn_pm * By_pm));
			
			double By_2 = Bxn_p0 / rho_p0 * (-1 / dx * (Qe(i + 1, j, 5) * x(M_B(i + 1, j, 6)) * reflct_xR - Qe(i, j, 5) * x(M_B(i, j, 6))) + 1 / dy * ((0.5 * (Bxn_pp * Bx_pp + Byn_pp * By_pp + Bzn_pp * Bz_pp) - Byn_pp * By_pp) - (0.5 * (Bxn_pm * Bx_pm + Byn_pm * By_pm + Bzn_pm * Bz_pm ) - Byn_pm * By_pm )));
			
			double By_3 = Byn_m0 / rho_m0 * (-1 / dx * ((0.5 * (Qe(i, j, 5) * x(M_B(i,j,5)) + Qe(i, j, 6) * x(M_B(i,j,6)) + Qe(i, j, 7) * x(M_B(i, j, 7)) ) - Qe(i, j, 5) * x(M_B(i, j, 5))) - (0.5 * ((Qe(i - 1, j, 5) * x(M_B(i - 1, j, 5)) * reflct_xL + Qe(i - 1, j, 6) * x(M_B(i - 1, j, 6)) * reflct_xL + Qe(i - 1, j, 7) * x(M_B(i - 1, j, 7)) * reflct_xL)) - Qe(i - 1, j, 5) * x(M_B(i - 1, j, 5)) * reflct_xL )) + 1 / dy * (Bxn_mp * By_mp - Bxn_mm * By_mm));
			
			double By_4 = Bxn_m0 / rho_m0 * (-1 / dx * (Qe(i, j, 5) * x(M_B(i, j, 6)) - Qe(i - 1, j, 5) * x(M_B(i - 1, j, 6)) * reflct_xL) + 1 / dy * ((0.5 * (Bxn_mp * Bx_mp + Byn_mp * By_mp + Bzn_mp * Bz_mp) - Byn_mp * By_mp)  - (0.5 * (Bxn_mm * Bx_mm + Byn_mm * By_mm + Bzn_mm * Bz_mm ) - Byn_mm * By_mm )));
			
			Ax(M_B(i,j,6)) = x(M_B(i,j,6)) + (h * h / dx) * ((By_1 + By_2) - (By_3 + By_4));
			

			double Bz_1 = Bzn_p0 / rho_p0 * (-1/dx * ((0.5 * (Qe(i + 1,j,5) * x(M_B(i + 1,j,5)) * reflct_xR + Qe(i + 1,j,6) * x(M_B(i + 1,j,6))* reflct_xR + Qe(i + 1,j,7) * x(M_B(i + 1,j,7))* reflct_xR) - Qe(i + 1,j,5) * x(M_B(i + 1,j,5))* reflct_xR ) - (0.5 * (Qe(i,j,5) * x(M_B(i,j,5)) + Qe(i,j,6) * x(M_B(i,j,6)) + Qe(i,j,7) * x(M_B(i,j,7))) - Qe(i,j,5) * x(M_B(i,j,5)))) + 1 / dy * (Byn_pp * Bx_pp - Byn_pm * Bx_pm ) ); 
			
			double Bz_2 =  Bxn_p0 / rho_p0 * (-1 / dx * (Qe(i + 1,j,5) * x(M_B(i + 1, j, 7)) * reflct_xR - Qe(i,j,5) * x(M_B(i, j, 7))) - 1 / dy * (Byn_pp * Bz_pp - Byn_pm * Bz_pm));
			
			double Bz_3 = Bzn_m0 / rho_m0 * (-1/dx * ((0.5 * (Qe(i,j,5) * x(M_B(i,j,5)) + Qe(i,j,6) * x(M_B(i,j,6)) + Qe(i,j,7) * x(M_B(i,j,7))) - Qe(i,j,5) * x(M_B(i,j,5))) - (0.5 * (Qe(i - 1,j,5) * x(M_B(i - 1,j,5)) * reflct_xL + Qe(i - 1,j,6) * x(M_B(i - 1,j,6)) * reflct_xL + Qe(i - 1,j,7) * x(M_B(i - 1,j,7))* reflct_xL) - Qe(i-1,j,5) * x(M_B(i-1,j,5)) * reflct_xL )) + 1 / dy * (Byn_mp * Bx_mp - Byn_mm * Bx_mm ) ); 
			
			double Bz_4 =  Bxn_m0 / rho_m0 * (-1 / dx * (Qe(i,j,5) * x(M_B(i, j, 7)) - Qe(i - 1,j,5) * x(M_B(i - 1, j, 7)) * reflct_xL ) - 1 / dy * (Byn_mp * Bz_mp - Byn_mm * Bz_mm));
			
			double Bz_5 = Bzn_0p / rho_0p * (1 / dx * (Bxn_pp * By_pp - Bxn_mp * By_mp) - 1 / dy * ((0.5 * (Qe(i,j+1,5) * x(M_B(i,j+1,5)) * reflct_yR + Qe(i,j+1,6) * x(M_B(i,j+1,6))* reflct_yR + Qe(i,j+1,7) * x(M_B(i,j+1,7)) * reflct_yR ) - Qe(i,j+1,6) * x(M_B(i,j+1,6))* reflct_yR)  - (0.5 * (Qe(i,j,5) * x(M_B(i,j,5)) + Qe(i,j,6) * x(M_B(i,j,6)) + Qe(i,j,7) * x(M_B(i,j,7)) ) - Qe(i,j,6) * x(M_B(i,j,6)))) );
			
			double Bz_6 = Byn_0p / rho_0p * (-1 / dx * (Bxn_pp * Bz_pp - Bxn_mp * Bz_mp) - 1 / dy * (Qe(i,j + 1,6) * x(M_B(i, j + 1, 7)) * reflct_yR - Qe(i,j,6) * x(M_B(i, j, 7))) );
			
			double Bz_7 = Bzn_0m / rho_0m * (1 / dx * (Bxn_pm * By_pm - Bxn_mm * By_mm) - 1 / dy * ((0.5 * (Qe(i,j,5) * x(M_B(i,j,5)) + Qe(i,j,6) * x(M_B(i,j,6)) + Qe(i,j,7) * x(M_B(i,j,7))) - Qe(i,j,6) * x(M_B(i,j,6))) - (0.5 * (Qe(i,j - 1,5) * x(M_B(i,j - 1,5)) * reflct_yL + Qe(i,j - 1,6) * x(M_B(i,j - 1,6)) * reflct_yL + Qe(i,j - 1,7) * x(M_B(i,j - 1,7)) * reflct_yL ) - Qe(i,j - 1,6) * x(M_B(i,j - 1,6))* reflct_yL)) );
			
			double Bz_8 = Byn_0m / rho_0m * (-1 / dx * (Bxn_pm * Bz_pm - Bxn_mm * Bz_mm) - 1 / dy * (Qe(i,j,6) * x(M_B(i, j, 7)) - Qe(i,j - 1,6) * x(M_B(i, j - 1, 7))* reflct_yL) );
			

			
			Ax(M_B(i,j,7)) = x(M_B(i,j,7)) + h *h  / dx * ((Bz_1 + Bz_2) - (Bz_3 + Bz_4)) + h * h/ dy * ((Bz_5 + Bz_6) - (Bz_7 + Bz_8));
			
			double lap_Bx, lap_By, lap_Bz;
			double nu_x = 0.0, nu_y = 0.0;   // switch below covers every enumerator
			
			if(GLF){
				nu_x = B_max_x * dx;
				nu_y = B_max_y * dy;
				lap_Bx = nu_x/(dx * dx) * (x(M_B(i-1,j,5)) * reflct_xL - 2.0*x(M_B(i,j,5)) + x(M_B(i+1,j,5)) * reflct_xR) + nu_y/(dy * dy) * (x(M_B(i,j - 1,5)) * reflct_yL - 2.0*x(M_B(i,j,5)) + x(M_B(i,j + 1,5)) * reflct_yR);
				Ax(M_B(i,j,5)) = Ax(M_B(i,j,5)) - h * lap_Bx;
				
				lap_By = nu_x/(dx * dx) * (x(M_B(i - 1,j,6)) * reflct_xL - 2.0*x(M_B(i,j,6)) + x(M_B(i + 1,j,6)) * reflct_xR) + nu_y/(dy * dy) * (x(M_B(i,j - 1,6)) * reflct_yL - 2.0*x(M_B(i,j,6)) + x(M_B(i,j + 1,6)) * reflct_yR);
				
				Ax(M_B(i,j,6)) = Ax(M_B(i,j,6)) - h * lap_By;
				
				lap_Bz = nu_x/(dx * dx) * (x(M_B(i - 1,j,7)) * reflct_xL - 2.0*x(M_B(i,j,7)) + x(M_B(i + 1,j,7)) * reflct_xR) + nu_y/(dy * dy) * (x(M_B(i,j - 1,7)) * reflct_yL - 2.0*x(M_B(i,j,7)) + x(M_B(i,j + 1,7)) * reflct_yR);
				Ax(M_B(i,j,7)) = Ax(M_B(i,j,7)) - h * lap_Bz;
			}else{
				double b_x = 0.5 * (fabs(Qe_pri(i,j,1)) + sqrt(Qe_pri(i,j,1) * Qe_pri(i,j,1) + 4 * (Qe_pri(i,j,5) * Qe_pri(i,j,5) + Qe_pri(i,j,6) * Qe_pri(i,j,6) + Qe_pri(i,j,7) * Qe_pri(i,j,7)) / Qe_pri(i,j,0)  ));
				
				double b_y = 0.5 * (fabs(Qe_pri(i,j,2)) + sqrt(Qe_pri(i,j,2) * Qe_pri(i,j,2) + 4 * (Qe_pri(i,j,5) * Qe_pri(i,j,5) + Qe_pri(i,j,6) * Qe_pri(i,j,6) + Qe_pri(i,j,7) * Qe_pri(i,j,7)) / Qe_pri(i,j,0)  ));
				
				double b_x_m = 0.5 * (fabs(Qe_pri(i - 1,j,1)) + sqrt(Qe_pri(i - 1,j,1) * Qe_pri(i - 1,j,1) + 4 * (Qe_pri(i - 1,j,5) * Qe_pri(i - 1,j,5) + Qe_pri(i - 1,j,6) * Qe_pri(i - 1,j,6) + Qe_pri(i - 1,j,7) * Qe_pri(i - 1,j,7)) / Qe_pri(i - 1,j,0)  ));
				
				double b_y_m = 0.5 * (fabs(Qe_pri(i,j - 1,2)) + sqrt(Qe_pri(i,j - 1,2) * Qe_pri(i,j - 1,2) + 4 * (Qe_pri(i,j - 1,5) * Qe_pri(i,j - 1,5) + Qe_pri(i,j - 1,6) * Qe_pri(i,j - 1,6) + Qe_pri(i,j - 1,7) * Qe_pri(i,j - 1,7)) / Qe_pri(i,j - 1,0)  ));
				
				double b_x_p = 0.5 * (fabs(Qe_pri(i + 1,j,1)) + sqrt(Qe_pri(i + 1,j,1) * Qe_pri(i + 1,j,1) + 4 * (Qe_pri(i + 1,j,5) * Qe_pri(i + 1,j,5) + Qe_pri(i + 1,j,6) * Qe_pri(i + 1,j,6) + Qe_pri(i + 1,j,7) * Qe_pri(i + 1,j,7)) / Qe_pri(i + 1,j,0)  ));
				
				double b_y_p = 0.5 * (fabs(Qe_pri(i,j + 1,2)) + sqrt(Qe_pri(i,j + 1,2) * Qe_pri(i,j + 1,2) + 4 * (Qe_pri(i,j + 1,5) * Qe_pri(i,j + 1,5) + Qe_pri(i,j + 1,6) * Qe_pri(i,j + 1,6) + Qe_pri(i,j + 1,7) * Qe_pri(i,j + 1,7)) / Qe_pri(i,j + 1,0)  ));
				
				double nu_x_p = fmax(fabs(b_x), fabs(b_x_p)) * dx;
				double nu_x_m = fmax(fabs(b_x), fabs(b_x_m)) * dx;
				double nu_y_p = fmax(fabs(b_y), fabs(b_y_p)) * dy;
				double nu_y_m = fmax(fabs(b_y), fabs(b_y_m)) * dy;
				
				lap_Bx = 1 /(dx * dx) * ( nu_x_p * (x(M_B(i+1,j,5)) * reflct_xR - x(M_B(i,j,5))) - nu_x_m * (x(M_B(i,j,5)) - x(M_B(i-1,j,5))* reflct_xL)) + 1 /(dy * dy) * ( nu_y_p * (x(M_B(i,j + 1,5)) * reflct_yR -  x(M_B(i,j,5)))  - nu_y_m * (x(M_B(i,j,5)) - x(M_B(i,j - 1,5)) * reflct_yL));
				Ax(M_B(i,j,5)) = Ax(M_B(i,j,5)) - h * lap_Bx;
				
				lap_By = 1 /(dx * dx) * ( nu_x_p * (x(M_B(i+1,j,6)) * reflct_xR - x(M_B(i,j,6))) - nu_x_m * (x(M_B(i,j,6)) - x(M_B(i-1,j,6))* reflct_xL)) + 1 /(dy * dy) * ( nu_y_p * (x(M_B(i,j + 1,6)) * reflct_yR -  x(M_B(i,j,6)))  - nu_y_m * (x(M_B(i,j,6)) - x(M_B(i,j - 1,6)) * reflct_yL));
				
				Ax(M_B(i,j,6)) = Ax(M_B(i,j,6)) - h * lap_By;
				
				lap_Bz = 1 /(dx * dx) * ( nu_x_p * (x(M_B(i+1,j,7)) * reflct_xR - x(M_B(i,j,7))) - nu_x_m * (x(M_B(i,j,7)) - x(M_B(i-1,j,7))* reflct_xL)) + 1 /(dy * dy) * ( nu_y_p * (x(M_B(i,j + 1,7)) * reflct_yR -  x(M_B(i,j,7)))  - nu_y_m * (x(M_B(i,j,7)) - x(M_B(i,j - 1,7)) * reflct_yL));
				Ax(M_B(i,j,7)) = Ax(M_B(i,j,7)) - h * lap_Bz;
			}
		}
	}
}

void Magnetic_Subsystem_potential_AX(
	Eigen::VectorXd& Ax_out, Eigen::VectorXd& x,
	Grid& u_con, Grid& Qe, Grid& Qs, double h
) {
	Grid Qe_pri(nxCells, nyCells);
	ucon2upri(Qe_pri, Qs);   // A-viscosity 特征值用 Qs(隐式累加器),对齐 CASTOR avisc=MagneticEigenvalues(Qs)

	// GLF: A-viscosity 系数全局化 —— 对齐 CASTOR (SI_MHD.f90 avisc(iDim,:)=MAXVAL(avisc(iDim,1:nElem)))
	// 每级从 Qs(=Qe_pri)在 real cell 上取方向全局 max,而非主循环里步首 u_pri 的 B_max_x/y。
	double gmax_bx = 0.0, gmax_by = 0.0;
	if(GLF){
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++)
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			double Bsq = Qe_pri(i,j,5)*Qe_pri(i,j,5)+Qe_pri(i,j,6)*Qe_pri(i,j,6)+Qe_pri(i,j,7)*Qe_pri(i,j,7);
			double bxv = 0.5*(fabs(Qe_pri(i,j,1))+sqrt(Qe_pri(i,j,1)*Qe_pri(i,j,1)+4.0*Bsq/Qe_pri(i,j,0)));
			double byv = 0.5*(fabs(Qe_pri(i,j,2))+sqrt(Qe_pri(i,j,2)*Qe_pri(i,j,2)+4.0*Bsq/Qe_pri(i,j,0)));
			if(bxv > gmax_bx) gmax_bx = bxv;
			if(byv > gmax_by) gmax_by = byv;
		}
	}

	// ★ 把GMRES向量复制到Grid, 做线性外推BC
	Grid AG(nxCells, nyCells);
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++)
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			AG(i,j,0) = x((i-ghost_cells) + nxCells*(j-ghost_cells));
			AG(i,j,1) = x((i-ghost_cells) + nxCells*(j-ghost_cells) + nxCells*nyCells);
			AG(i,j,2) = x((i-ghost_cells) + nxCells*(j-ghost_cells) + 2*nxCells*nyCells);
		}
	bound_cond_A(AG);
	
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			double Bxn_00 = Qe(i,j,5),   Byn_00 = Qe(i,j,6),   Bzn_00 = Qe(i,j,7);
			double Bxn_P0 = Qe(i+1,j,5), Byn_P0 = Qe(i+1,j,6), Bzn_P0 = Qe(i+1,j,7);
			double Bxn_M0 = Qe(i-1,j,5), Byn_M0 = Qe(i-1,j,6), Bzn_M0 = Qe(i-1,j,7);
			double Bxn_0P = Qe(i,j+1,5), Byn_0P = Qe(i,j+1,6), Bzn_0P = Qe(i,j+1,7);
			double Bxn_0M = Qe(i,j-1,5), Byn_0M = Qe(i,j-1,6), Bzn_0M = Qe(i,j-1,7);
			double rho_00 = u_con(i,j,0);
			
			// B^{n+1} = curl(A^{n+1}) — 全部从AG读
			
			double bx_P0 = (AG(i+1,j+1,2) - AG(i+1,j-1,2)) / (2.0*dy);
			double by_P0 = -(AG(i+2,j,2) - AG(i,j,2)) / (2.0*dx);
			double bz_P0 = (AG(i+2,j,1) - AG(i,j,1)) / (2.0*dx)
			- (AG(i+1,j+1,0) - AG(i+1,j-1,0)) / (2.0*dy);
			
			double bx_M0 = (AG(i-1,j+1,2) - AG(i-1,j-1,2)) / (2.0*dy);
			double by_M0 = -(AG(i,j,2) - AG(i-2,j,2)) / (2.0*dx);
			double bz_M0 = (AG(i,j,1) - AG(i-2,j,1)) / (2.0*dx)
			- (AG(i-1,j+1,0) - AG(i-1,j-1,0)) / (2.0*dy);
			
			double bx_0P = (AG(i,j+2,2) - AG(i,j,2)) / (2.0*dy);
			double by_0P = -(AG(i+1,j+1,2) - AG(i-1,j+1,2)) / (2.0*dx);
			double bz_0P = (AG(i+1,j+1,1) - AG(i-1,j+1,1)) / (2.0*dx)
			- (AG(i,j+2,0) - AG(i,j,0)) / (2.0*dy);
			
			double bx_0M = (AG(i,j,2) - AG(i,j-2,2)) / (2.0*dy);
			double by_0M = -(AG(i+1,j-1,2) - AG(i-1,j-1,2)) / (2.0*dx);
			double bz_0M = (AG(i+1,j-1,1) - AG(i-1,j-1,1)) / (2.0*dx)
			- (AG(i,j,0) - AG(i,j-2,0)) / (2.0*dy);
			
			// ==============================================================
			// 张量 T^B(α,β):  α 是动量分量(用 B_new=b), β 是通量方向(用 B_n)
			// 对应 Fortran MagneticFluxSI_momentum:
			//   FB(α,β) = B28Pe·δ_{αβ} - B_n^β · B_new^α     (μ0=1)
			//   其中 B28Pe = 0.5 · (B_n · B_new)
			// ==============================================================
			
			// —— x 方向通量 (β=x)，用在 P0, M0 位置 ——
			// 需要：T(2,1)=Txx, T(3,1)=Tyx, T(4,1)=Tzx
			double bdot_P0 = Bxn_P0*bx_P0 + Byn_P0*by_P0 + Bzn_P0*bz_P0;
			double Txx_P0 = bdot_P0/2.0 - Bxn_P0*bx_P0;   // 对角：α=β=x, 对称不变
			double Tyx_P0 = -Bxn_P0*by_P0;                 // α=y, β=x:  -B_n^x · B_new^y
			double Tzx_P0 = -Bxn_P0*bz_P0;                 // α=z, β=x:  -B_n^x · B_new^z
			
			double bdot_M0 = Bxn_M0*bx_M0 + Byn_M0*by_M0 + Bzn_M0*bz_M0;
			double Txx_M0 = bdot_M0/2.0 - Bxn_M0*bx_M0;
			double Tyx_M0 = -Bxn_M0*by_M0;
			double Tzx_M0 = -Bxn_M0*bz_M0;
			
			// —— y 方向通量 (β=y)，用在 0P, 0M 位置 ——
			// 需要：T(2,2)=Txy, T(3,2)=Tyy, T(4,2)=Tzy
			double bdot_0P = Bxn_0P*bx_0P + Byn_0P*by_0P + Bzn_0P*bz_0P;
			double Txy_0P = -Byn_0P*bx_0P;                 // α=x, β=y:  -B_n^y · B_new^x
			double Tyy_0P = bdot_0P/2.0 - Byn_0P*by_0P;   // 对角：α=β=y, 对称不变
			double Tzy_0P = -Byn_0P*bz_0P;                 // α=z, β=y:  -B_n^y · B_new^z
			
			double bdot_0M = Bxn_0M*bx_0M + Byn_0M*by_0M + Bzn_0M*bz_0M;
			double Txy_0M = -Byn_0M*bx_0M;
			double Tyy_0M = bdot_0M/2.0 - Byn_0M*by_0M;
			double Tzy_0M = -Byn_0M*bz_0M;
			
			// 散度: divT_α = ∂_x T(α,x) + ∂_y T(α,y)
			double divT_x = (Txx_P0-Txx_M0)/(2.0*dx) + (Txy_0P-Txy_0M)/(2.0*dy);
			double divT_y = (Tyx_P0-Tyx_M0)/(2.0*dx) + (Tyy_0P-Tyy_0M)/(2.0*dy);
			double divT_z = (Tzx_P0-Tzx_M0)/(2.0*dx) + (Tzy_0P-Tzy_0M)/(2.0*dy);
			
			double Fx = h/rho_00 * divT_x;
			double Fy = h/rho_00 * divT_y;
			double Fz = h/rho_00 * divT_z;
			
			double cross_x = Byn_00*Fz - Bzn_00*Fy;
			double cross_y = Bzn_00*Fx - Bxn_00*Fz;
			double cross_z = Bxn_00*Fy - Byn_00*Fx;
			
			// Laplacian — 也从AG读
			double lap_Ax = 0.0, lap_Ay = 0.0, lap_Az = 0.0;
			if(GLF){
				double nu_x = gmax_bx * dx;   // 每级 Qs real-cell 全局 max(对齐 CASTOR GLF avisc)
				double nu_y = gmax_by * dy;
				lap_Ax = nu_x/(dx*dx)*(AG(i-1,j,0)-2.0*AG(i,j,0)+AG(i+1,j,0))
				+ nu_y/(dy*dy)*(AG(i,j-1,0)-2.0*AG(i,j,0)+AG(i,j+1,0));
				lap_Ay = nu_x/(dx*dx)*(AG(i-1,j,1)-2.0*AG(i,j,1)+AG(i+1,j,1))
				+ nu_y/(dy*dy)*(AG(i,j-1,1)-2.0*AG(i,j,1)+AG(i,j+1,1));
				lap_Az = nu_x/(dx*dx)*(AG(i-1,j,2)-2.0*AG(i,j,2)+AG(i+1,j,2))
				+ nu_y/(dy*dy)*(AG(i,j-1,2)-2.0*AG(i,j,2)+AG(i,j+1,2));
			}else{
				double Bsq = Qe_pri(i,j,5)*Qe_pri(i,j,5)+Qe_pri(i,j,6)*Qe_pri(i,j,6)+Qe_pri(i,j,7)*Qe_pri(i,j,7);
				double b_x   = 0.5*(fabs(Qe_pri(i,j,1))+sqrt(Qe_pri(i,j,1)*Qe_pri(i,j,1)+4.0*Bsq/Qe_pri(i,j,0)));
				double b_x_p = 0.5*(fabs(Qe_pri(i+1,j,1))+sqrt(Qe_pri(i+1,j,1)*Qe_pri(i+1,j,1)+4.0*(Qe_pri(i+1,j,5)*Qe_pri(i+1,j,5)+Qe_pri(i+1,j,6)*Qe_pri(i+1,j,6)+Qe_pri(i+1,j,7)*Qe_pri(i+1,j,7))/Qe_pri(i+1,j,0)));
				double b_x_m = 0.5*(fabs(Qe_pri(i-1,j,1))+sqrt(Qe_pri(i-1,j,1)*Qe_pri(i-1,j,1)+4.0*(Qe_pri(i-1,j,5)*Qe_pri(i-1,j,5)+Qe_pri(i-1,j,6)*Qe_pri(i-1,j,6)+Qe_pri(i-1,j,7)*Qe_pri(i-1,j,7))/Qe_pri(i-1,j,0)));
				double b_y   = 0.5*(fabs(Qe_pri(i,j,2))+sqrt(Qe_pri(i,j,2)*Qe_pri(i,j,2)+4.0*Bsq/Qe_pri(i,j,0)));
				double b_y_p = 0.5*(fabs(Qe_pri(i,j+1,2))+sqrt(Qe_pri(i,j+1,2)*Qe_pri(i,j+1,2)+4.0*(Qe_pri(i,j+1,5)*Qe_pri(i,j+1,5)+Qe_pri(i,j+1,6)*Qe_pri(i,j+1,6)+Qe_pri(i,j+1,7)*Qe_pri(i,j+1,7))/Qe_pri(i,j+1,0)));
				double b_y_m = 0.5*(fabs(Qe_pri(i,j-1,2))+sqrt(Qe_pri(i,j-1,2)*Qe_pri(i,j-1,2)+4.0*(Qe_pri(i,j-1,5)*Qe_pri(i,j-1,5)+Qe_pri(i,j-1,6)*Qe_pri(i,j-1,6)+Qe_pri(i,j-1,7)*Qe_pri(i,j-1,7))/Qe_pri(i,j-1,0)));
				// 面平均(对齐 CASTOR ComputeCellLaplace: nup=0.5(avisc_i+avisc_{i+1}));原来用 max 不对
				double nu_xp = 0.5*(fabs(b_x)+fabs(b_x_p))*dx;
				double nu_xm = 0.5*(fabs(b_x)+fabs(b_x_m))*dx;
				double nu_yp = 0.5*(fabs(b_y)+fabs(b_y_p))*dy;
				double nu_ym = 0.5*(fabs(b_y)+fabs(b_y_m))*dy;
				for(int k = 0; k < 3; k++){
					double lk = 1.0/(dx*dx)*(nu_xp*(AG(i+1,j,k)-AG(i,j,k))-nu_xm*(AG(i,j,k)-AG(i-1,j,k)))
					+ 1.0/(dy*dy)*(nu_yp*(AG(i,j+1,k)-AG(i,j,k))-nu_ym*(AG(i,j,k)-AG(i,j-1,k)));
					if(k == 0)      lap_Ax = lk;
					else if(k == 1) lap_Ay = lk;
					else            lap_Az = lk;
				}
			}
			// 输出 — 直接用平坦索引, 不用M_B
			// A-solve numerical viscosity: OFF by default to match CASTOR (add_bvisc=.FALSE.).
			// 在低 β 强爆炸前沿,这个粘性会把 A 解搅坏 -> B 偏高 -> e=E-B^2/2 负压。
			if(!add_bvisc){ lap_Ax = 0.0; lap_Ay = 0.0; lap_Az = 0.0; }
			int idx = (i-ghost_cells) + nxCells*(j-ghost_cells);
			Ax_out(idx)                         = AG(i,j,0) - h*cross_x - h*lap_Ax;
			Ax_out(idx + nxCells*nyCells)        = AG(i,j,1) - h*cross_y - h*lap_Ay;
			Ax_out(idx + 2*nxCells*nyCells)      = AG(i,j,2) - h*cross_z - h*lap_Az;
		}
	}
}

void Pressure_Subsystem_AX(Eigen::VectorXd& Ax, Eigen::VectorXd& x, Grid& u_con, Grid& Qe,Grid& /*Qs*/,double h){
	Grid Qe_pri(nxCells, nyCells);
	ucon2upri(Qe_pri, Qe);
	Array shock_eta(nxCells,nyCells);
	ShockDetector(shock_eta, Qe);
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
			
			
			double h00 = calenthalpy(con00_n) * Qe(i,j,0) / u_con(i,j,0);
			double hp0 = calenthalpy(conp0_n) * Qe(i + 1,j,0) / u_con(i + 1,j,0);
			double hm0 = calenthalpy(conm0_n) * Qe(i - 1,j,0) / u_con(i - 1,j,0);
			double h0p = calenthalpy(con0p_n) * Qe(i,j + 1,0) / u_con(i,j + 1,0);
			double h0m = calenthalpy(con0m_n) * Qe(i,j - 1,0) / u_con(i,j - 1,0);

			// Dirichlet: homogeneous part (ghost = -mirror) — 照抄 epsilon_BC 的 get_x_unknown
			// Dirichlet: -镜像(对齐参照 2DEuler_epsilon_BC get_x_unknown);lifting 用 1×(get_x_const)
			auto getp = [&](int ii, int jj) -> double {
				int II = ii, JJ = jj; double s = 1.0;
				if(x_bc == x_Dirichlet){
					if(ii < ghost_cells){ II = 2*ghost_cells - 1 - ii; s = -s; }
					else if(ii >= nxCells + ghost_cells){ II = 2*(nxCells+ghost_cells) - 1 - ii; s = -s; }
				}
				if(y_bc == y_Dirichlet){
					if(jj < ghost_cells){ JJ = 2*ghost_cells - 1 - jj; s = -s; }
					else if(jj >= nyCells + ghost_cells){ JJ = 2*(nyCells+ghost_cells) - 1 - jj; s = -s; }
				}
				return s * x(M_P(II, JJ));
			};
			double pii = x(M_P(i,j));

			// 强 Dirichlet:边界 cell 在 LHS 钉死 p=p_amb(恒等行),配 RHS b_p_vect=p_amb
			if( ((x_bc==x_Dirichlet)&&(i==ghost_cells||i==nxCells+ghost_cells-1)) || ((y_bc==y_Dirichlet)&&(j==ghost_cells||j==nyCells+ghost_cells-1)) ){ Ax(M_P(i,j)) = pii; continue; }
			// eq (52): Lagrange degree-2 three-point stencil for d/dx(h dp/dx)
			// V(p) 对角项:floor_on 时用冻结斜率 V_slope(p_freeze)(此函数即作 Newton Jacobian);
			// 否则退回原线性斜率 1/(gamma-1)。耦合项(焓 Laplacian + 黏性)完全不变。
			double Vdiag = g_AX_jacobian ? V_slope( p_freeze(M_P(i,j)), u_con(i,j,0) ) : (1.0/(gamma_plasma - 1.0));
//			Ax(M_P(i,j)) = Vdiag * pii - h * h * (
//				1 / (dx * dx) * ( (0.75*hm0 + 0.25*hp0)*getp(i-1,j) - (hm0 + hp0)*pii + (0.25*hm0 + 0.75*hp0)*getp(i+1,j) )
//			  + 1 / (dy * dy) * ( (0.75*h0m + 0.25*h0p)*getp(i,j-1) - (h0m + h0p)*pii + (0.25*h0m + 0.75*h0p)*getp(i,j+1) ) );
			
			// jet inflow strip: ghost-cell Dirichlet 齐次部分(ghost 未知量=0);
			// 已知值 jet_p 的 lifting 在 PreSubsystem_Pressure 的 RHS 组装处(与本系数严格配对)
			double p_0m = (j == ghost_cells && jet_inlet(i)) ? 0.0 : x(M_P(i,j-1));
			Ax(M_P(i,j)) = Vdiag * pii - h * h * (1 / (dx * dx) * (0.5 * (hp0 + h00) * (x(M_P(i + 1, j)) - x(M_P(i,j))) - 0.5 * (h00 + hm0) * (x(M_P(i,j)) - x(M_P(i-1,j)))) + 1 / (dy * dy) * (0.5 * (h0p + h00) * (x(M_P(i,j+1)) - x(M_P(i,j))) - 0.5 * (h0m + h00) * (x(M_P(i,j)) - p_0m)));
			
			//Laplacian Operator
			double nu_x = 0.0, nu_y = 0.0;   // switch below covers every enumerator
			double nu_x_p = 0.0, nu_y_p = 0.0;   // switch below covers every enumerator
			double nu_x_m = 0.0, nu_y_m = 0.0;
			
			double lap_p;
			if(GLF){
				switch(add_pvisc){
					case MaterialSpeedp:
						nu_x = vx_max * dx;
						nu_y = vy_max * dy;
						break;
					case AlfvenSpeedp:
						nu_x = cA_x * dx;
						nu_y = cA_y * dy;
						break;
					case Noaddvisp:
						nu_x = 0;
						nu_y = 0;
						break;
				}
				lap_p = nu_x/(dx * dx) * ( 0.5 * (shock_eta(i + 1,j) + shock_eta(i,j)) * (getp(i+1,j) - pii) - 0.5 * (shock_eta(i,j) + shock_eta(i - 1,j)) * (pii - getp(i - 1,j))) + nu_y/(dy * dy) * ( 0.5 * (shock_eta(i,j + 1) + shock_eta(i,j)) * (getp(i,j + 1) - pii) - 0.5 * (shock_eta(i,j) + shock_eta(i,j - 1)) * (pii - getp(i,j-1)));
			}else{
				switch(add_pvisc){
					case MaterialSpeedp:
						nu_x_p = fmax(fabs(Qe_pri(i,j,1)), fabs(Qe_pri(i + 1,j,1))) * dx;
						nu_x_m = fmax(fabs(Qe_pri(i,j,1)), fabs(Qe_pri(i - 1,j,1))) * dx;
						nu_y_p = fmax(fabs(Qe_pri(i,j,2)), fabs(Qe_pri(i,j + 1,2))) * dy;
						nu_y_m = fmax(fabs(Qe_pri(i,j,2)), fabs(Qe_pri(i,j - 1,2))) * dy;
						break;
					case AlfvenSpeedp:
						nu_x_p = fmax(fabs(Qe_pri(i,j,5)  / sqrt(Qe_pri(i,j,0))), fabs(Qe_pri(i + 1,j,5)  / sqrt(Qe_pri(i + 1,j,0)))) * dx;
						nu_x_m = fmax(fabs(Qe_pri(i,j,5)  / sqrt(Qe_pri(i,j,0))), fabs(Qe_pri(i - 1,j,5)  / sqrt(Qe_pri(i - 1,j,0)))) * dx;
						nu_y_p = fmax(fabs(Qe_pri(i,j,6)  / sqrt(Qe_pri(i,j,0))), fabs(Qe_pri(i,j + 1,6)  / sqrt(Qe_pri(i,j + 1,0)))) * dy;
						nu_y_m = fmax(fabs(Qe_pri(i,j,6)  / sqrt(Qe_pri(i,j,0))), fabs(Qe_pri(i,j - 1,6)  / sqrt(Qe_pri(i,j - 1,0)))) * dy;
						break;
					case Noaddvisp:
						nu_x_p = 0;
						nu_x_m = 0;
						nu_y_p = 0;
						nu_y_m = 0;
						break;
				}
				lap_p =  1/(dx * dx) * ( 0.5 * (shock_eta(i + 1,j) + shock_eta(i,j)) * nu_x_p *(getp(i+1,j) - pii) - 0.5 * (shock_eta(i,j) + shock_eta(i - 1,j)) * nu_x_m *(pii - getp(i - 1,j))) + 1/(dy * dy) * ( 0.5 * (shock_eta(i,j + 1) + shock_eta(i,j)) * nu_y_p * (getp(i,j + 1) - pii) - 0.5 * (shock_eta(i,j) + shock_eta(i,j - 1)) * nu_y_m * (pii - getp(i,j-1)));
			}
			// 黏性项(仅激波处激活的耗散)统一按局部 dV/dp = d(rho*e)/dp 缩放:
			// 理想气体 Vbase_slope == 1/(gamma-1)(逐字节不变),RK == rho*de/dp(冻结于 p_freeze,保持算子线性)。
			double visc_scale = Vbase_slope(Qe_pri(i,j,4), u_con(i,j,0));
			Ax(M_P(i,j)) = Ax(M_P(i,j)) - visc_scale * h * lap_p;
//			Ax(M_P(i,j)) = Ax(M_P(i,j)) - h * lap_p;
		}
	}
}

// ===================== 嵌套 Newton:保正压力求解 =====================
// 解非线性系统  g(p) = V_mod(p) + T p - b = 0,  其中 T = 焓 Laplacian + 黏性(线性,Qe 冻结)。
// 复用 Pressure_Subsystem_AX,由 g_AX_jacobian 选择算子形态:
//   g_AX_jacobian=false 时返回线性载体算子 L p = (1/(gamma-1)) p + T p;
//   g_AX_jacobian=true  时返回 Newton Jacobian  J x = V_slope(p_freeze,rho) x + T x。
// 残差 g = L(p) - (1/(gamma-1)) p + V_mod(p,rho) - b。载体常数 1/(gamma-1) 精确抵消,
// 提取出 EOS 无关的耦合 T;故对理想气体与 RK 同一套机制(论文 Eq.(52)-(53))。
void Pressure_NestedNewton(Eigen::VectorXd& p, Eigen::VectorXd& b,
                           Grid& u_con, Grid& Qe, Grid& Qs, double h, const int dim){
	const double a = 1.0/(gamma_plasma - 1.0);          // 载体常数(精确抵消)
	Eigen::VectorXd g(dim), dp(dim), neg_g(dim), Lp(dim);
	// Jacobi 预条件(对齐 precE 的做法): 对角 = V 斜率 + h^2 焓 Laplacian 对角。
	// 焓部分与模式无关,预计算一次;V 斜率按 g_AX_jacobian/p_freeze 在应用时取。
	Eigen::VectorXd Tdiag(dim);
	{
		Grid Qe_pri_pc(nxCells, nyCells);
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
				double c00[NUM_VARS], cp0[NUM_VARS], cm0[NUM_VARS], c0p[NUM_VARS], c0m[NUM_VARS];
				for(int k = 0; k < NUM_VARS; k++){
					c00[k]=Qe(i,j,k); cp0[k]=Qe(i+1,j,k); cm0[k]=Qe(i-1,j,k); c0p[k]=Qe(i,j+1,k); c0m[k]=Qe(i,j-1,k);
				}
				double h00 = calenthalpy(c00)*Qe(i,j,0)/u_con(i,j,0);
				double hp0 = calenthalpy(cp0)*Qe(i+1,j,0)/u_con(i+1,j,0);
				double hm0 = calenthalpy(cm0)*Qe(i-1,j,0)/u_con(i-1,j,0);
				double h0p = calenthalpy(c0p)*Qe(i,j+1,0)/u_con(i,j+1,0);
				double h0m = calenthalpy(c0m)*Qe(i,j-1,0)/u_con(i,j-1,0);
				Tdiag(M_P(i,j)) = h*h*( (0.5*(hp0+h00) + 0.5*(hm0+h00))/(dx*dx)
				                      + (0.5*(h0p+h00) + 0.5*(h0m+h00))/(dy*dy) );
			}
		}
	}
	auto precP = [&](Eigen::VectorXd& out, const Eigen::VectorXd& in){
		out.resize(in.size());
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
				int idx = M_P(i,j);
				bool bnd = ((x_bc==x_Dirichlet)&&(i==ghost_cells||i==nxCells+ghost_cells-1))
				         || ((y_bc==y_Dirichlet)&&(j==ghost_cells||j==nyCells+ghost_cells-1));
				if(bnd){ out(idx) = in(idx); continue; }   // 恒等行
				double Vd = g_AX_jacobian ? V_slope(p_freeze(idx), u_con(i,j,0)) : a;
				out(idx) = in(idx) / (Vd + Tdiag(idx));
			}
		}
	};
	// 仅在 "理想气体 + 关保正" 时退回单次线性 GMRES;RK 是非线性 EOS,必须走 Newton。
	if(!floor_on && eos_type==IDEAL_GAS){
		g_AX_jacobian = false;           // 线性算子 L(对角=1/(gamma-1))
		GMRES_Solver(p, b, Pressure_Subsystem_AX, u_con, Qe, Qs, h, dim, precP);
		return;
	}
	// ---- Cor.1:按本次 b 自适应选陡坡斜率 Mbig(先验保证 p >= theta*p_floor)----
	if(floor_on){
		double bmin = b.minCoeff();
		double Vpe;
		if(eos_type == IDEAL_GAS){
			Vpe = p_floor/(gamma_plasma - 1.0);
		}else{                                  // RK:V 逐格,取上界 max_i rho_i*e(rho_i,p_floor)
			Vpe = p_floor/(gamma_plasma - 1.0); // EOS 反解失败时的兜底
			for(int i = ghost_cells; i < nxCells + ghost_cells; i++)
				for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
					double v = Vbase(p_floor, u_con(i,j,0));
					if(std::isfinite(v) && v > Vpe) Vpe = v;
				}
		}
		double M_cor  = std::max(0.0, Vpe - bmin) / ((1.0 - theta)*p_floor);
		double M_base = Mbig_fac/(gamma_plasma - 1.0);   // 旧固定值作下限(光滑算例逐位不变)
		Mbig = std::max(M_base, M_cor);
	}
	p_freeze = p;                        // Newton 初值 = 当前压力
	int it = 0;
	double res0 = -1.0;                  // 首个残差,做相对判据基准
	for(it = 0; it < Pnewton_max; it++){
		// ---- 残差 g(p_freeze):用线性载体算子 L ----
		g_AX_jacobian = false;
		Pressure_Subsystem_AX(Lp, p_freeze, u_con, Qe, Qs, h);
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
				int idx = M_P(i,j);
				bool bnd = ((x_bc==x_Dirichlet)&&(i==ghost_cells||i==nxCells+ghost_cells-1))
				         || ((y_bc==y_Dirichlet)&&(j==ghost_cells||j==nyCells+ghost_cells-1));
				if(bnd) g(idx) = Lp(idx) - b(idx);                          // 恒等行: p - p_amb
				else    g(idx) = Lp(idx) - a*p_freeze(idx) + V_mod(p_freeze(idx), u_con(i,j,0)) - b(idx);
			}
		}
		double res = g.norm();
		if(it == 0) res0 = res;
		// 相对判据 + 绝对地板:jet 等大尺度问题(p~1e5, KE~4.5e5)下,
		// 纯绝对 1e-10 等价于要求相对 ~1e-21,机器精度不可达 -> 永远 hit max iters
		if(res < Pnewton_tol * res0 + Pnewton_tol){
//			std::cout << "  [P-Newton] converged it=" << it << " res=" << res << std::endl;
			break;
		}
		// ---- Jacobian 解:J dp = -g(g_AX_jacobian=true,AX 即 Jacobian) ----
		neg_g = -g;
		dp.setZero();
		g_AX_jacobian = true;
		GMRES_Solver(dp, neg_g, Pressure_Subsystem_AX, u_con, Qe, Qs, h, dim, precP);
		p_freeze += dp;
//		std::cout << "  [P-Newton] it=" << it << " res=" << res << " |dp|=" << dp.norm() << std::endl;
	}
	g_AX_jacobian = false;               // 复位,避免影响后续(理想气体)线性求解
	if(it == Pnewton_max)
		std::cout << "  [P-Newton] WARNING: hit max iters (" << Pnewton_max << ")" << std::endl;
	p = p_freeze;
	// ---- 正性诊断:floored 计数 + Cor.1 运行时校验(理论成立则永不触发)----
	if(floor_on){
		int n_floored = 0; double pmin = 1e300;
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++)
			for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
				double pv = p(M_P(i,j));
				if(pv < p_floor) ++n_floored;
				if(pv < pmin)    pmin = pv;
			}
		if(n_floored > g_max_floored) g_max_floored = n_floored;
		if(pmin < g_min_p_ever)       g_min_p_ever  = pmin;
		if(pmin < theta*p_floor - 1e-14)
			std::cout << "  [POSITIVITY VIOLATED] pmin=" << pmin
			          << " < theta*p_floor=" << theta*p_floor << std::endl;
	}
}
// ====================================================================

void Energy_Subsystem_AX(Eigen::VectorXd& Ax, Eigen::VectorXd& x, Grid& u_con, Grid& Qe,Grid& /*Qs*/,double h){
	Grid Qe_pri(nxCells, nyCells);
	ucon2upri(Qe_pri, Qe);
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
			
			double hhp_x = 0.5 * (h_00 + h_p0);
			double hhm_x = 0.5 * (h_00 + h_m0);
			double hhp_y = 0.5 * (h_00 + h_0p);
			double hhm_y = 0.5 * (h_00 + h_0m);

			// Dirichlet: homogeneous part (ghost = -mirror) — 照抄 epsilon_BC 的 get_x_unknown
			auto getE = [&](int ii, int jj) -> double {   // CASTOR ghost-cell Dirichlet: 齐次 ghost=0(同 getp,避免 -镜像×2 角点漏)
				bool xg = (ii < ghost_cells) || (ii >= nxCells + ghost_cells);
				bool yg = (jj < ghost_cells) || (jj >= nyCells + ghost_cells);
				if(x_bc == x_Dirichlet && xg) return 0.0;
				if(y_bc == y_Dirichlet && yg) return 0.0;
				if(jj < ghost_cells && jet_inlet(ii)) return 0.0;   // jet strip 底 ghost: Dirichlet 齐次;lift 在 getEc
				return x(M_P(ii, jj));
			};
			double Ei = x(M_P(i,j));
			Ax(M_P(i,j)) = Ei
			- (gamma_plasma - 1) * h * h * (
				1.0/(dx*dx) * ( hhp_x * (getE(i+1,j) - Ei) - hhm_x * (Ei - getE(i-1,j)) )
				+ 1.0/(dy*dy) * ( hhp_y * (getE(i,j+1) - Ei) - hhm_y * (Ei - getE(i,j-1)) )
			);

			Ax(M_P(i,j)) = Ax(M_P(i,j)) * dx * dy;
		}
	}
}


void GMRES_Solver(
	Eigen::VectorXd &x,
	Eigen::VectorXd &b,
	std::function<void(Eigen::VectorXd&, Eigen::VectorXd&, Grid&, Grid&, Grid&, double)> Ax_func,
	Grid &u_con, Grid &Qe, Grid &Qs, double h,
	const int dim,
	std::function<void(Eigen::VectorXd&, const Eigen::VectorXd&)> precond)
{
	const int JMAX = 25;
	const int outmax = dim;
	// left preconditioner: applyM(out,in) = M^{-1} in  (identity if no precond)
	auto applyM = [&](Eigen::VectorXd& out, const Eigen::VectorXd& in){
		if(precond) precond(out, in); else out = in;
	};

	std::vector<Eigen::VectorXd> v(JMAX + 1, Eigen::VectorXd(dim));
	std::vector<Eigen::VectorXd> Av(JMAX + 1, Eigen::VectorXd(dim));
	Eigen::VectorXd w(dim), Ax(dim), r0(dim), tmp(dim);
	
	std::vector<std::vector<double>> H(JMAX + 2, std::vector<double>(JMAX + 2, 0.0));
	std::vector<double> alpha(JMAX + 2, 0.0);
	std::vector<double> gam(JMAX + 2, 0.0);
	std::vector<double> s(JMAX + 2, 0.0);
	std::vector<double> c(JMAX + 2, 0.0);
	
	double beta;
	int j = 0;
	double tol_eff = GMREStol;   // 相对判据: 首个 outer 的初始残差 x GMREStol + 绝对地板
	
	for (int outer = 0; outer < outmax; outer++) {
		
		Ax.setZero();
		Ax_func(Ax, x, u_con, Qe, Qs, h);
		tmp = b - Ax;
		applyM(r0, tmp);            // r0 = M^{-1}(b - A x)

		gam[0] = r0.norm();
		if(outer == 0) tol_eff = GMREStol * gam[0] + GMREStol;

		if (gam[0] <= tol_eff) {
//			std::cout << "GMRES converged: outer=" << outer
//			<< " residual=" << gam[0] << std::endl;
			return;
		}

		v[0] = r0 / gam[0];
		Ax.setZero();
		Ax_func(Ax, v[0], u_con, Qe, Qs, h);
		applyM(Av[0], Ax);          // Av[0] = M^{-1} A v[0]
		
		// Reset H
		for (int ii = 0; ii < JMAX + 2; ii++)
			std::fill(H[ii].begin(), H[ii].end(), 0.0);
		
		for (j = 0; j < JMAX; j++) {
			
			// Arnoldi: modified Gram-Schmidt
			for (int i = 0; i <= j; i++) {
				H[i][j] = v[i].dot(Av[j]);
			}
			w = Av[j];
			for (int i = 0; i <= j; i++) {
				w -= H[i][j] * v[i];
			}
			H[j + 1][j] = w.norm();
			
			// Apply previous Givens rotations
			for (int i = 0; i < j; i++) {
				double tmp1 = c[i + 1] * H[i][j] + s[i + 1] * H[i + 1][j];
				double tmp2 = s[i + 1] * H[i][j] - c[i + 1] * H[i + 1][j];
				H[i][j]     = tmp1;
				H[i + 1][j] = tmp2;
			}
			
			// New Givens rotation
			beta = std::sqrt(H[j][j] * H[j][j] + H[j + 1][j] * H[j + 1][j]);
			s[j + 1] = H[j + 1][j] / beta;
			c[j + 1] = H[j][j] / beta;
			H[j][j]  = beta;
			gam[j + 1] = s[j + 1] * gam[j];
			gam[j]     = c[j + 1] * gam[j];
			
			if (gam[j + 1] > tol_eff) {
				v[j + 1] = w / H[j + 1][j];
				Ax.setZero();
				Ax_func(Ax, v[j + 1], u_con, Qe, Qs, h);
				applyM(Av[j + 1], Ax);   // Av[j+1] = M^{-1} A v[j+1]
			} else {
				// Converged: back-substitute and update x
				for (int i = j; i >= 0; i--) {
					alpha[i] = gam[i];
					for (int k = i + 1; k <= j; k++) {
						alpha[i] -= H[i][k] * alpha[k];
					}
					alpha[i] /= H[i][i];
				}
				for (int i = 0; i <= j; i++) {
					x += alpha[i] * v[i];
				}
//				std::cout << "GMRES converged: outer=" << outer
//				<< " inner=" << j + 1
//				<< " residual=" << gam[j + 1] << std::endl;
				return;
			}
		}
		
		// Not converged: back-substitute and restart
		std::fill(alpha.begin(), alpha.end(), 0.0);
		for (int i = JMAX - 1; i >= 0; i--) {
			alpha[i] = gam[i];
			for (int k = i + 1; k < JMAX; k++) {
				alpha[i] -= H[i][k] * alpha[k];
			}
			alpha[i] /= H[i][i];
		}
		for (int i = 0; i < JMAX; i++) {
			x += alpha[i] * v[i];
		}
		
		if (outer % 10 == 0 || outer < 3) {
//			std::cout << "GMRES restart " << outer
//			<< " residual=" << gam[JMAX] << std::endl;
		}
	}
	
	std::cout << "WARNING: GMRES reached max outer iterations, residual=" << gam[0] << std::endl;
}
