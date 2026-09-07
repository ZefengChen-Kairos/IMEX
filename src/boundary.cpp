// Boundary conditions for the cell-centred state and the vector potential.

#include "declarations.hpp"

void bound_cond_u(Grid &u){
	switch (x_bc) {
		case x_Transmissive:
			for(int j = 0; j < nyCells + ghost_cells * 2; j++){
				for(int i = 0; i < ghost_cells; i++){
					for(int k = 0; k < NUM_VARS; k++){
						u(i, j, k) = u(ghost_cells, j, k);
						u(nxCells + ghost_cells + i, j, k) = u(nxCells + ghost_cells - 1, j, k);
					}
				}
			}
			break;
		case x_Reflective:
			for(int j = 0; j < nyCells + ghost_cells * 2; j++){
				for(int ghost = 0; ghost < ghost_cells; ghost++){
					u(ghost,j,0) = u(2 * ghost_cells - 1 - ghost,j,0);
					u(nxCells + ghost_cells + ghost,j,0) = u(nxCells + ghost_cells - 1 - ghost,j,0);
					
					u(ghost,j,1) = -u(2 * ghost_cells - 1 - ghost,j,1);
					u(nxCells + ghost_cells + ghost,j,1) = -u(nxCells + ghost_cells - 1 - ghost,j,1);
					
					u(ghost,j,2) = -u(2 * ghost_cells - 1 - ghost,j,2);
					u(nxCells + ghost_cells + ghost,j,2) = -u(nxCells + ghost_cells - 1 - ghost,j,2);
					
					u(ghost,j,3) = -u(2 * ghost_cells - 1 - ghost,j,3);
					u(nxCells + ghost_cells + ghost,j,3) = -u(nxCells + ghost_cells - 1 - ghost,j,3);
					
					u(ghost,j,4) = u(2 * ghost_cells - 1 - ghost,j,4);
					u(nxCells + ghost_cells + ghost,j,4) = u(nxCells + ghost_cells - 1 - ghost,j,4);
					
					u(ghost,j,5) = -u(2 * ghost_cells - 1 - ghost,j,5);
					u(nxCells + ghost_cells + ghost,j,5) = -u(nxCells + ghost_cells - 1 - ghost,j,5);
					
					u(ghost,j,6) = -u(2 * ghost_cells - 1 - ghost,j,6);
					u(nxCells + ghost_cells + ghost,j,6) = -u(nxCells + ghost_cells - 1 - ghost,j,6);
					
					u(ghost,j,7) = -u(2 * ghost_cells - 1 - ghost,j,7);
					u(nxCells + ghost_cells + ghost,j,7) = -u(nxCells + ghost_cells - 1 - ghost,j,7);
				}
			}
			break;
		case x_Periodic:
			for(int j = 0; j < nyCells + ghost_cells * 2; j++){
				for(int ghost = 0; ghost < ghost_cells; ghost++){
					for(int k = 0; k < NUM_VARS; k++){
						u(ghost,j,k) = u(nxCells + ghost,j,k);
						u(nxCells + ghost_cells + ghost,j,k) = u(ghost_cells + ghost,j,k);
					}
				}
			}
			break;
		case x_Dirichlet:{
			for(int j = 0; j < nyCells + ghost_cells * 2; j++){
				for(int g = 0; g < ghost_cells; g++){
					int L = g, R = nxCells + ghost_cells + g;
					u(L,j,0)=rho_amb; u(L,j,1)=0; u(L,j,2)=0; u(L,j,3)=0;   // rho,v: Dirichlet ambient (v=0)
					u(R,j,0)=rho_amb; u(R,j,1)=0; u(R,j,2)=0; u(R,j,3)=0;
					for(int k = 5; k < NUM_VARS; k++){      // B: transmissive (copy interior) -> match CASTOR
						u(L,j,k) = u(ghost_cells, j, k);
						u(R,j,k) = u(nxCells + ghost_cells - 1, j, k);
					}
					// slot4: 钉压强 p_amb;守恒模式 E = p_amb/(g-1)+1/2 B_ghost^2(E 跟透射 B 浮动,内能恒=ambient)
					if(g_bc_prim){ u(L,j,4)=p_amb; u(R,j,4)=p_amb; }
					else{
						u(L,j,4)=p_amb/(gamma_plasma-1.0)+0.5*(u(L,j,5)*u(L,j,5)+u(L,j,6)*u(L,j,6)+u(L,j,7)*u(L,j,7));
						u(R,j,4)=p_amb/(gamma_plasma-1.0)+0.5*(u(R,j,5)*u(R,j,5)+u(R,j,6)*u(R,j,6)+u(R,j,7)*u(R,j,7));
					}
				}
			}
			break;
		}
	}

	switch(y_bc){
		case y_Transmissive:
			for(int i = 0; i < nxCells + ghost_cells * 2; i++){
				for(int j = 0; j < ghost_cells; j++){
					for(int k = 0; k < NUM_VARS; k++){
						u(i, j, k) = u(i, ghost_cells, k);
						u(i,nyCells + ghost_cells + j,k) = u(i, nyCells + ghost_cells - 1, k);
					}
				}
			}
			break;
		case y_Reflective:
			for(int i = 0; i < nxCells + ghost_cells * 2; i++){
				for(int ghost = 0; ghost < ghost_cells; ghost++){
					u(i,ghost,0) = u(i,2 * ghost_cells - 1 - ghost,0);
					u(i,nyCells + ghost_cells + ghost,0) = u(i,nyCells + ghost_cells - 1 - ghost,0);
					
					u(i,ghost,1) = -u(i,2 * ghost_cells - 1 - ghost,1);
					u(i,nyCells + ghost_cells + ghost,1) = -u(i,nyCells + ghost_cells - 1 - ghost,1);
					
					u(i,ghost,2) = -u(i,2 * ghost_cells - 1 - ghost,2);
					u(i,nyCells + ghost_cells + ghost,2) = -u(i,nyCells + ghost_cells - 1 - ghost,2);
					
					u(i,ghost,3) = -u(i,2 * ghost_cells - 1 - ghost,3);
					u(i,nyCells + ghost_cells + ghost,3) = -u(i,nyCells + ghost_cells - 1 - ghost,3);
					
					u(i,ghost,4) = u(i,2 * ghost_cells - 1 - ghost,4);
					u(i,nyCells + ghost_cells + ghost,4) = u(i,nyCells + ghost_cells - 1 - ghost,4);
					
					u(i,ghost,5) = -u(i,2 * ghost_cells - 1 - ghost,5);
					u(i,nyCells + ghost_cells + ghost,5) = -u(i,nyCells + ghost_cells - 1 - ghost,5);
					
					u(i,ghost,6) = -u(i,2 * ghost_cells - 1 - ghost,6);
					u(i,nyCells + ghost_cells + ghost,6) = -u(i,nyCells + ghost_cells - 1 - ghost,6);
					
					u(i,ghost,7) = -u(i,2 * ghost_cells - 1 - ghost,7);
					u(i,nyCells + ghost_cells + ghost,7) = -u(i,nyCells + ghost_cells - 1 - ghost,7);
				}
			}
			break;
		case y_Periodic:
			for(int i = 0; i < nxCells + ghost_cells * 2; i++){
				for(int ghost = 0; ghost < ghost_cells; ghost++){
					for(int k = 0; k < NUM_VARS; k++){
						u(i,ghost,k) = u(i,nyCells + ghost,k);
						u(i,nyCells + ghost_cells + ghost,k) = u(i,ghost_cells + ghost,k);
					}
				}
			}
			break;
		case y_Dirichlet:{
			for(int i = 0; i < nxCells + ghost_cells * 2; i++){
				for(int g = 0; g < ghost_cells; g++){
					int B = g, T = nyCells + ghost_cells + g;
					u(i,B,0)=rho_amb; u(i,B,1)=0; u(i,B,2)=0; u(i,B,3)=0;
					u(i,T,0)=rho_amb; u(i,T,1)=0; u(i,T,2)=0; u(i,T,3)=0;
					for(int k = 5; k < NUM_VARS; k++){      // B: transmissive (copy interior) -> match CASTOR
						u(i,B,k) = u(i, ghost_cells, k);
						u(i,T,k) = u(i, nyCells + ghost_cells - 1, k);
					}
					// slot4: 钉压强 p_amb;守恒模式 E = p_amb/(g-1)+1/2 B_ghost^2
					if(g_bc_prim){ u(i,B,4)=p_amb; u(i,T,4)=p_amb; }
					else{
						u(i,B,4)=p_amb/(gamma_plasma-1.0)+0.5*(u(i,B,5)*u(i,B,5)+u(i,B,6)*u(i,B,6)+u(i,B,7)*u(i,B,7));
						u(i,T,4)=p_amb/(gamma_plasma-1.0)+0.5*(u(i,T,5)*u(i,T,5)+u(i,T,6)*u(i,T,6)+u(i,T,7)*u(i,T,7));
					}
				}
			}
			break;
		}
		case y_JetInflow:{
			// 底边: |x|<0.05 -> Dirichlet jet 态 (rho,v,B,p)=(jet_rho, 0,jet_v,0, 0,jet_B0,0, jet_p)
			//       strip 外 -> outflow(透射)。顶边: outflow(透射)。
			for(int i = 0; i < nxCells + ghost_cells * 2; i++){
				for(int g = 0; g < ghost_cells; g++){
					int T = nyCells + ghost_cells + g;
					for(int k = 0; k < NUM_VARS; k++)               // 顶: 透射
						u(i,T,k) = u(i, nyCells + ghost_cells - 1, k);
					if(jet_inlet(i)){                                // 底 strip: jet Dirichlet
						u(i,g,0) = jet_rho;
						u(i,g,5) = 0.0;  u(i,g,6) = jet_B0;  u(i,g,7) = 0.0;
						if(g_bc_prim){                               // 原始量: slot1-3=v, slot4=p
							u(i,g,1) = 0.0;  u(i,g,2) = jet_v;  u(i,g,3) = 0.0;
							u(i,g,4) = jet_p;
						}else{                                       // 守恒量: slot1-3=rho v, slot4=E
							u(i,g,1) = 0.0;  u(i,g,2) = jet_rho * jet_v;  u(i,g,3) = 0.0;
							u(i,g,4) = jet_E_val();
						}
					}else{                                           // 底 strip 外: 透射
						for(int k = 0; k < NUM_VARS; k++)
							u(i,g,k) = u(i, ghost_cells, k);
					}
				}
			}
			break;
		}
	}
}

void bound_cond_A(Grid &A) {
	// ===== x方向 =====
	switch (x_bc) {
		case x_Transmissive: case x_Dirichlet: {
			// 线性外推 (论文 Eq.63); 对线性环境 A_z=Bx*y 精确，等价于 Dirichlet
			for(int j = 0; j < nyCells + ghost_cells * 2; j++){
				for(int k = 0; k < 3; k++){
					// 左边界: 斜率从最内两个cell取
					double slope_L = A(ghost_cells + 1, j, k) - A(ghost_cells, j, k);
					for(int g = 0; g < ghost_cells; g++){
						A(ghost_cells - 1 - g, j, k) = A(ghost_cells, j, k) - (g + 1) * slope_L;
					}
					// 右边界
					int iR = nxCells + ghost_cells - 1;
					double slope_R = A(iR, j, k) - A(iR - 1, j, k);
					for(int g = 0; g < ghost_cells; g++){
						A(iR + 1 + g, j, k) = A(iR, j, k) + (g + 1) * slope_R;
					}
				}
			}
			break;
		}
		case x_Reflective: {
			// A 在反射边界的行为需要特别处理
			// 对于 Dirichlet 类型的边界, 也用线性外推
			for(int j = 0; j < nyCells + ghost_cells * 2; j++){
				for(int k = 0; k < 3; k++){
					double slope_L = A(ghost_cells + 1, j, k) - A(ghost_cells, j, k);
					for(int g = 0; g < ghost_cells; g++){
						A(ghost_cells - 1 - g, j, k) = A(ghost_cells, j, k) - (g + 1) * slope_L;
					}
					int iR = nxCells + ghost_cells - 1;
					double slope_R = A(iR, j, k) - A(iR - 1, j, k);
					for(int g = 0; g < ghost_cells; g++){
						A(iR + 1 + g, j, k) = A(iR, j, k) + (g + 1) * slope_R;
					}
				}
			}
			break;
		}
		case x_Periodic: {
			// 周期BC: 正常复制
			for(int j = 0; j < nyCells + ghost_cells * 2; j++){
				for(int ghost = 0; ghost < ghost_cells; ghost++){
					for(int k = 0; k < 3; k++){
						A(ghost, j, k) = A(nxCells + ghost, j, k);
						A(nxCells + ghost_cells + ghost, j, k) = A(ghost_cells + ghost, j, k);
					}
				}
			}
			break;
		}
	}
	// ===== y方向 =====
	switch (y_bc) {
		case y_Transmissive: case y_Dirichlet: case y_JetInflow: {
			// jet: Az=-B0*x 沿 y 为常数,线性外推精确;strip 内 Ez=0 -> A 在入口不动
			for(int i = 0; i < nxCells + ghost_cells * 2; i++){
				for(int k = 0; k < 3; k++){
					// 下边界
					double slope_L = A(i, ghost_cells + 1, k) - A(i, ghost_cells, k);
					for(int g = 0; g < ghost_cells; g++){
						A(i, ghost_cells - 1 - g, k) = A(i, ghost_cells, k) - (g + 1) * slope_L;
					}
					// 上边界
					int jR = nyCells + ghost_cells - 1;
					double slope_R = A(i, jR, k) - A(i, jR - 1, k);
					for(int g = 0; g < ghost_cells; g++){
						A(i, jR + 1 + g, k) = A(i, jR, k) + (g + 1) * slope_R;
					}
				}
			}
			break;
		}
		case y_Reflective: {
			for(int i = 0; i < nxCells + ghost_cells * 2; i++){
				for(int k = 0; k < 3; k++){
					double slope_L = A(i, ghost_cells + 1, k) - A(i, ghost_cells, k);
					for(int g = 0; g < ghost_cells; g++){
						A(i, ghost_cells - 1 - g, k) = A(i, ghost_cells, k) - (g + 1) * slope_L;
					}
					int jR = nyCells + ghost_cells - 1;
					double slope_R = A(i, jR, k) - A(i, jR - 1, k);
					for(int g = 0; g < ghost_cells; g++){
						A(i, jR + 1 + g, k) = A(i, jR, k) + (g + 1) * slope_R;
					}
				}
			}
			break;
		}
		case y_Periodic: {
			for(int i = 0; i < nxCells + ghost_cells * 2; i++){
				for(int ghost = 0; ghost < ghost_cells; ghost++){
					for(int k = 0; k < 3; k++){
						A(i, ghost, k) = A(i, nyCells + ghost, k);
						A(i, nyCells + ghost_cells + ghost, k) = A(i, ghost_cells + ghost, k);
					}
				}
			}
			break;
		}
	}
	// ---------- MHD jet: 底边 inflow strip -> A 偶反射 ----------
	// 偶反射使入口面 dAz/dy = 0 精确成立 => Bx|face = 0(严格 Dirichlet);
	// By = -dAz/dx 由入流维持在 B0。线性外推(上面)等价于对 Bx 透射,strip 内不够强。
	// 该规则同时作用于 Magnetic_Subsystem_potential_AX 里的 AG(算子内闭合)与求解后的
	// B=curl(A) 重构,两处自动一致;且反射是关于 A 的线性映射,算子仍为合法线性算子。
	if(y_bc == y_JetInflow){
		for(int i = 0; i < nxCells + ghost_cells * 2; i++){
			if(!jet_inlet(i)) continue;
			for(int k = 0; k < 3; k++)
				for(int g = 0; g < ghost_cells; g++)
					A(i, ghost_cells - 1 - g, k) = A(i, ghost_cells + g, k);
		}
	}
}
