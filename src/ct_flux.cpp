// EMF-consistent staggered CT built from the implicit magnetic fluxes.

#include "declarations.hpp"

void CollectMagFluxes(Grid &u_con, Grid &Qe, Grid &Qs, double h){
	Grid Qe_pri(nxCells, nyCells);
	ucon2upri(Qe_pri, Qe);
	Array &Fx = *gFluxBy_p0;
	Array &Fy = *gFluxBx_0p;
	
	for(int i = 0; i < nxCells + 2*ghost_cells; i++)
		for(int j = 0; j < nyCells + 2*ghost_cells; j++){ Fx(i,j) = 0.0; Fy(i,j) = 0.0; }
	
	for(int i = 1; i < nxCells + 2*ghost_cells - 1; i++){
		for(int j = 1; j < nyCells + 2*ghost_cells - 1; j++){
			
			// ======== 面/角平均: 与 B* 装配 (4928-4999) 和 AX (3676-3746) 同式 ========
			double rho_p0 = 0.5*(u_con(i,j,0) + u_con(i+1,j,0));
			double rho_0p = 0.5*(u_con(i,j,0) + u_con(i,j+1,0));
			
			double rvx_p0 = 0.5*(u_con(i,j,1) + u_con(i+1,j,1));
			double rvy_p0 = 0.5*(u_con(i,j,2) + u_con(i+1,j,2));
			double rvx_0p = 0.5*(u_con(i,j,1) + u_con(i,j+1,1));
			double rvy_0p = 0.5*(u_con(i,j,2) + u_con(i,j+1,2));
			
			double p_pp = 0.25*(Qe_pri(i,j,4)+Qe_pri(i+1,j,4)+Qe_pri(i,j+1,4)+Qe_pri(i+1,j+1,4));
			double p_pm = 0.25*(Qe_pri(i,j,4)+Qe_pri(i+1,j,4)+Qe_pri(i,j-1,4)+Qe_pri(i+1,j-1,4));
			double p_mp = 0.25*(Qe_pri(i,j,4)+Qe_pri(i-1,j,4)+Qe_pri(i,j+1,4)+Qe_pri(i-1,j+1,4));
			
			// B^n (Qe) 的面/角平均
			double Bxn_p0 = 0.5*(Qe(i,j,5) + Qe(i+1,j,5));
			double Byn_p0 = 0.5*(Qe(i,j,6) + Qe(i+1,j,6));
			double Bxn_0p = 0.5*(Qe(i,j,5) + Qe(i,j+1,5));
			double Byn_0p = 0.5*(Qe(i,j,6) + Qe(i,j+1,6));
			
			double Bxn_pp = 0.25*(Qe(i,j,5)+Qe(i+1,j,5)+Qe(i,j+1,5)+Qe(i+1,j+1,5));
			double Bxn_pm = 0.25*(Qe(i,j,5)+Qe(i+1,j,5)+Qe(i,j-1,5)+Qe(i+1,j-1,5));
			double Bxn_mp = 0.25*(Qe(i,j,5)+Qe(i-1,j,5)+Qe(i,j+1,5)+Qe(i-1,j+1,5));
			double Byn_pp = 0.25*(Qe(i,j,6)+Qe(i+1,j,6)+Qe(i,j+1,6)+Qe(i+1,j+1,6));
			double Byn_pm = 0.25*(Qe(i,j,6)+Qe(i+1,j,6)+Qe(i,j-1,6)+Qe(i+1,j-1,6));
			double Byn_mp = 0.25*(Qe(i,j,6)+Qe(i-1,j,6)+Qe(i,j+1,6)+Qe(i-1,j+1,6));
			double Bzn_pp = 0.25*(Qe(i,j,7)+Qe(i+1,j,7)+Qe(i,j+1,7)+Qe(i+1,j+1,7));
			double Bzn_pm = 0.25*(Qe(i,j,7)+Qe(i+1,j,7)+Qe(i,j-1,7)+Qe(i+1,j-1,7));
			double Bzn_mp = 0.25*(Qe(i,j,7)+Qe(i-1,j,7)+Qe(i,j+1,7)+Qe(i-1,j+1,7));
			
			// B̃^{n+1} (u_con 解) 的角平均 —— AX 里 x(M_B(...)) 的对应物
			double Bx1_pp = 0.25*(u_con(i,j,5)+u_con(i+1,j,5)+u_con(i,j+1,5)+u_con(i+1,j+1,5));
			double Bx1_pm = 0.25*(u_con(i,j,5)+u_con(i+1,j,5)+u_con(i,j-1,5)+u_con(i+1,j-1,5));
			double Bx1_mp = 0.25*(u_con(i,j,5)+u_con(i-1,j,5)+u_con(i,j+1,5)+u_con(i-1,j+1,5));
			double By1_pp = 0.25*(u_con(i,j,6)+u_con(i+1,j,6)+u_con(i,j+1,6)+u_con(i+1,j+1,6));
			double By1_pm = 0.25*(u_con(i,j,6)+u_con(i+1,j,6)+u_con(i,j-1,6)+u_con(i+1,j-1,6));
			double By1_mp = 0.25*(u_con(i,j,6)+u_con(i-1,j,6)+u_con(i,j+1,6)+u_con(i-1,j+1,6));
			double Bz1_pp = 0.25*(u_con(i,j,7)+u_con(i+1,j,7)+u_con(i,j+1,7)+u_con(i+1,j+1,7));
			double Bz1_pm = 0.25*(u_con(i,j,7)+u_con(i+1,j,7)+u_con(i,j-1,7)+u_con(i+1,j-1,7));
			double Bz1_mp = 0.25*(u_con(i,j,7)+u_con(i-1,j,7)+u_con(i,j+1,7)+u_con(i-1,j+1,7));
			
			// ======== x 面 (i+1/2, j): By 的通量  →  Ez_p0 = -Fx ========
			// F^n 部分: By_star (5005 行) 第一个大括号, 原样
			double Fn_x = Byn_p0/rho_p0 * (rvx_p0 - (h/dx)*(Qe_pri(i+1,j,4) - Qe_pri(i,j,4)))
			- Bxn_p0/rho_p0 * (rvy_p0 - (h/dy)*(p_pp - p_pm));
			
			// 隐式部分: AX 的 By_1 + By_2 (3764-3766 行), x(M_B(a,b,k)) → u_con(a,b,k)
			double By_1 = Byn_p0/rho_p0 * (
				-1.0/dx * ( (0.5*( Qe(i+1,j,5)*u_con(i+1,j,5)
					+ Qe(i+1,j,6)*u_con(i+1,j,6)
					+ Qe(i+1,j,7)*u_con(i+1,j,7)) - Qe(i+1,j,5)*u_con(i+1,j,5))
					- (0.5*( Qe(i,j,5)*u_con(i,j,5)
						+ Qe(i,j,6)*u_con(i,j,6)
						+ Qe(i,j,7)*u_con(i,j,7)) - Qe(i,j,5)*u_con(i,j,5)) )
				+ 1.0/dy * (Bxn_pp*By1_pp - Bxn_pm*By1_pm) );
			
			double By_2 = Bxn_p0/rho_p0 * (
				-1.0/dx * (Qe(i+1,j,5)*u_con(i+1,j,6) - Qe(i,j,5)*u_con(i,j,6))
				+ 1.0/dy * ( (0.5*(Bxn_pp*Bx1_pp + Byn_pp*By1_pp + Bzn_pp*Bz1_pp) - Byn_pp*By1_pp)
					- (0.5*(Bxn_pm*Bx1_pm + Byn_pm*By1_pm + Bzn_pm*Bz1_pm) - Byn_pm*By1_pm) ) );
			
			// 隐式 LLF 扩散通量 (可表示为 EMF 的分量): -ν_x·∂By/∂x |_{i+1/2}
			double nu_x;
			if(GLF){
				nu_x = B_max_x * dx;
			}else{
				double b_i  = 0.5*(fabs(Qe_pri(i,j,1))   + sqrt(Qe_pri(i,j,1)*Qe_pri(i,j,1)
					+ 4.0*(Qe_pri(i,j,5)*Qe_pri(i,j,5)+Qe_pri(i,j,6)*Qe_pri(i,j,6)+Qe_pri(i,j,7)*Qe_pri(i,j,7))/Qe_pri(i,j,0)));
				double b_ip = 0.5*(fabs(Qe_pri(i+1,j,1)) + sqrt(Qe_pri(i+1,j,1)*Qe_pri(i+1,j,1)
					+ 4.0*(Qe_pri(i+1,j,5)*Qe_pri(i+1,j,5)+Qe_pri(i+1,j,6)*Qe_pri(i+1,j,6)+Qe_pri(i+1,j,7)*Qe_pri(i+1,j,7))/Qe_pri(i+1,j,0)));
				nu_x = fmax(fabs(b_i), fabs(b_ip)) * dx;
			}
			double Fdiff_x = - nu_x * (u_con(i+1,j,6) - u_con(i,j,6)) / dx;
			
//			Fx(i,j) = Fn_x + h*(By_1 + By_2);
			Fx(i,j) = Fn_x + h*(By_1 + By_2) + Fdiff_x;
			
			// ======== y 面 (i, j+1/2): Bx 的通量  →  Ez_0p = +Fy ========
			// F^n 部分: Bx_star (5003 行) 的 (…)_0p 组, 原样
			double Fn_y = Bxn_0p/rho_0p * (rvy_0p - (h/dy)*(Qe_pri(i,j+1,4) - Qe_pri(i,j,4)))
			- Byn_0p/rho_0p * (rvx_0p - (h/dx)*(p_pp - p_mp));
			
			// 隐式部分: AX 的 Bx_1 + Bx_2 (3753-3755 行)
			double Bx_1 = Bxn_0p/rho_0p * (
				1.0/dx * (Bxn_pp*By1_pp - Bxn_mp*By1_mp)
				- 1.0/dy * ( (0.5*( Qe(i,j+1,5)*u_con(i,j+1,5)
					+ Qe(i,j+1,6)*u_con(i,j+1,6)
					+ Qe(i,j+1,7)*u_con(i,j+1,7)) - Qe(i,j+1,6)*u_con(i,j+1,6))
					- (0.5*( Qe(i,j,5)*u_con(i,j,5)
						+ Qe(i,j,6)*u_con(i,j,6)
						+ Qe(i,j,7)*u_con(i,j,7)) - Qe(i,j,6)*u_con(i,j,6)) ) );
			
			double Bx_2 = Byn_0p/rho_0p * (
				1.0/dx * ( (0.5*(Bxn_pp*Bx1_pp + Byn_pp*By1_pp + Bzn_pp*Bz1_pp) - Bxn_pp*Bx1_pp)
					- (0.5*(Bxn_mp*Bx1_mp + Byn_mp*By1_mp + Bzn_mp*Bz1_mp) - Bxn_mp*Bx1_mp) )
				- 1.0/dy * (Qe(i,j+1,6)*u_con(i,j+1,5) - Qe(i,j,6)*u_con(i,j,5)) );
			
			double nu_y;
			if(GLF){
				nu_y = B_max_y * dy;
			}else{
				double b_j  = 0.5*(fabs(Qe_pri(i,j,2))   + sqrt(Qe_pri(i,j,2)*Qe_pri(i,j,2)
					+ 4.0*(Qe_pri(i,j,5)*Qe_pri(i,j,5)+Qe_pri(i,j,6)*Qe_pri(i,j,6)+Qe_pri(i,j,7)*Qe_pri(i,j,7))/Qe_pri(i,j,0)));
				double b_jp = 0.5*(fabs(Qe_pri(i,j+1,2)) + sqrt(Qe_pri(i,j+1,2)*Qe_pri(i,j+1,2)
					+ 4.0*(Qe_pri(i,j+1,5)*Qe_pri(i,j+1,5)+Qe_pri(i,j+1,6)*Qe_pri(i,j+1,6)+Qe_pri(i,j+1,7)*Qe_pri(i,j+1,7))/Qe_pri(i,j+1,0)));
				nu_y = fmax(fabs(b_j), fabs(b_jp)) * dy;
			}
			double Fdiff_y = - nu_y * (u_con(i,j+1,5) - u_con(i,j,5)) / dy;
			
//			Fy(i,j) = Fn_y + h*(Bx_1 + Bx_2);
			Fy(i,j) = Fn_y + h*(Bx_1 + Bx_2) + Fdiff_y;
		}
	}
	bound_cond_array(Fx);
	bound_cond_array(Fy);
	
	// ---- 自检: 通量差分 vs 胞心解的实际增量 ----
	// 残差 = Laplacian 中不可表示为 EMF 的分量 (Bx 的 x 向 / By 的 y 向扩散)
	// 预期 ~O(ν·ΔB) 的小量且随网格收敛; 若为 O(1), 说明通量转录有误, 必须先修。
	if(g_ct_selfcheck){
		double res_max = 0.0;
		for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
			for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
				double By_pred = Qs(i,j,6) - (h/dx)*(Fx(i,j) - Fx(i-1,j));
				double Bx_pred = Qs(i,j,5) - (h/dy)*(Fy(i,j) - Fy(i,j-1));
				double r = fmax(fabs(By_pred - u_con(i,j,6)), fabs(Bx_pred - u_con(i,j,5)));
				if(r > res_max) res_max = r;
			}
		}
		printf("[CT-CONS] flux-vs-cell residual max = %.3e (expected ~O(nu*dB), NOT machine zero)\n", res_max);
	}
}

// 替换 Staggered_CT (安装步骤[5])。角点装配复用诊断里的 CT_CornerAssemble。
void Staggered_CT_Consistent(Grid& u_con, Array &Bx_p0, Array &By_0p,
	Array &Bx_p0_Qs, Array &By_0p_Qs, double h){
		Array Ez_p0(nxCells,nyCells), Ez_0p(nxCells,nyCells);
		Array Ez_tilde(nxCells,nyCells), Ez_pp(nxCells,nyCells);
		
		// Riccardo eq (43): EMF = 通量对偶
		for(int i = 0; i < nxCells + 2*ghost_cells; i++){
			for(int j = 0; j < nyCells + 2*ghost_cells; j++){
				Ez_p0(i,j) = -(*gFluxBy_p0)(i,j);
				Ez_0p(i,j) = +(*gFluxBx_0p)(i,j);
				Ez_pp(i,j) = 0.0;
			}
		}
		
		// eq (45): 胞心 Ez 用 (u^{n+1}, B̃^{n+1}) —— CT 调用时 u_con 恰是这个状态
		for(int j = 0; j < nyCells + 2*ghost_cells; j++)
			for(int i = 0; i < nxCells + 2*ghost_cells; i++)
			Ez_tilde(i,j) = u_con(i,j,2)/u_con(i,j,0)*u_con(i,j,5)
		- u_con(i,j,1)/u_con(i,j,0)*u_con(i,j,6);
		
		// eq (41)/(42)/(44): GS05 角点装配 —— 复用诊断代码里的 CT_CornerAssemble
		CT_CornerAssemble(Ez_pp, Ez_p0, Ez_0p, Ez_tilde);
		
		// eq (39): 面场演化 (与原 Staggered_CT 相同)
		for(int j = ghost_cells - 1; j < nyCells + ghost_cells; j++){
			for(int i = ghost_cells - 1; i < nxCells + ghost_cells; i++){
				Bx_p0(i,j) = Bx_p0_Qs(i,j) - h / dy * (Ez_pp(i,j) - Ez_pp(i,j - 1));
				By_0p(i,j) = By_0p_Qs(i,j) + h / dx * (Ez_pp(i,j) - Ez_pp(i - 1,j));
			}
		}
		bound_cond_staggered(Bx_p0, 0);
		bound_cond_staggered(By_0p, 1);
		
		// eq (46): 胞心 = 面平均, 覆盖预测子 B̃^{n+1}
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
				u_con(i,j,5) = 0.5 * (Bx_p0(i - 1,j) + Bx_p0(i,j));
				u_con(i,j,6) = 0.5 * (By_0p(i,j - 1) + By_0p(i,j));
			}
		}
		bound_cond_u(u_con);
	}
