// TVD reconstruction, limiters and Gauss-Legendre quadrature helpers.

#include "declarations.hpp"

void TVDReconstruction(Grid &u_con, Grid &w_xL, Grid &w_xR,Grid &w_yL, Grid &w_yR){
	double c = 0.0;   // limited slope; every reconstruction branch below assigns it
	double q0, q1, q2,q0L, q1L, q2L,q0R, q1R, q2R, wr0, wr1, wr2, wr0_bar, wr1_bar, wr2_bar, beta0, beta1, beta2, d0, d1, d2;
	double diff1, diff2, diff3, diff4;
	Grid u_pri(nxCells,nyCells);
	ucon2upri(u_pri, u_con);
	Grid u(nxCells,nyCells);
	
	for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
		for(int i = ghost_cells - 1; i < nxCells + ghost_cells + 1; i++){
			for(int k = 0; k < NUM_VARS; k++){
				switch (reconstate) {
					case ConsState:
						u(i,j,k) = u_con(i,j,k);
						break;
					case PrimState:
						u(i,j,k) = u_pri(i,j,k);
						break;
				}
			}
		}
	}
	bound_cond_u(u);
		
		
	for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
		for(int i = ghost_cells - 1; i < nxCells + ghost_cells + 1; i++){
			for(int k = 0; k < NUM_VARS; k++){
				double slope_xR = u(i + 1,j,k) - u(i,j,k);
				double slope_xL = u(i,j,k) - u(i - 1,j,k);
				
				double slope_yR = u(i,j + 1,k) - u(i,j,k);
				double slope_yL = u(i,j,k) - u(i,j - 1,k);
				switch (recon) {
					case NoRecon:
						SpaceOrder = 1;
						w_xL(i,j,k) = u(i,j,k);
						w_xR(i,j,k) = u(i,j,k);
						w_yL(i,j,k) = u(i,j,k);
						w_yR(i,j,k) = u(i,j,k);
						break;
					case Minmod:
						SpaceOrder = 2;
						if(slope_xL * slope_xR <= 0){
							c = 0;
						}else if(fabs(slope_xR) <= fabs(slope_xL)){
							c = slope_xR;
						}else if(fabs(slope_xL) < fabs(slope_xR)){
							c = slope_xL;
						}
						w_xL(i,j,k) = u(i,j,k) - 0.5 * c;
						w_xR(i,j,k) = u(i,j,k) + 0.5 * c;
						
						if(test == testMHD_vortex){
							c = 0.5 * (slope_xL + slope_xR);
						}
						
						if(slope_yL * slope_yR <= 0){
							c = 0;
						}else if(fabs(slope_yR) <= fabs(slope_yL)){
							c = slope_yR;
						}else if(fabs(slope_yL) < fabs(slope_yR)){
							c = slope_yL;
						}
						
						if(test == testMHD_vortex){
							c = 0.5 * (slope_yL + slope_yR);
						}
						w_yL(i,j,k) = u(i,j,k) - 0.5 * c;
						w_yR(i,j,k) = u(i,j,k) + 0.5 * c;
						break;
					case VanLeer:
						SpaceOrder = 2;
						if(fabs(slope_xR) < 1e-16){
							c = 0;
						}else{
							double r = slope_xL/slope_xR;
							if(r < 0){
								c = 0;
							}else{
								c = fmin(2/(1 + r), (2 * r)/(1 + r));
							}
						}
						w_xL(i,j,k) = u(i,j,k) - 0.5 * c * slope_xL;
						w_xR(i,j,k) = u(i,j,k) + 0.5 * c * slope_xR;
						
						if(fabs(slope_yR) < 1e-16){
							c = 0;
						}else{
							double r = slope_yL/slope_yR;
							if(r < 0){
								c = 0;
							}else{
								c = fmin(2/(1 + r), (2 * r)/(1 + r));
							}
						}
						w_yL(i,j,k) = u(i,j,k) - 0.5 * c * slope_yL;
						w_yR(i,j,k) = u(i,j,k) + 0.5 * c * slope_yR;
						break;
					case ENO3:
						SpaceOrder = 3;
						q0L = -1.0/6.0 * u(i - 2, j, k) + 5.0/6.0 * u(i - 1, j, k) + 1.0/3.0 * u(i, j, k);
						q1L = 1.0/3.0 * u(i - 1, j, k) + 5.0/6.0 * u(i, j, k) - 1.0/6.0 * u(i + 1, j, k);
						q2L = 11.0/6.0 * u(i, j, k) - 7.0/6.0 * u(i + 1, j, k) + 1.0/3.0 * u(i + 2, j, k);
						
						
						q0R = 1.0/3.0 * u(i - 2, j, k) - 7.0/6.0 * u(i - 1, j, k) + 11.0/6.0 * u(i, j, k);
						q1R = -1.0/6.0 * u(i - 1, j, k) + 5.0/6.0 * u(i, j, k) + 1.0/3.0 * u(i + 1, j, k);
						q2R = 1.0/3.0 * u(i, j, k) + 5.0/6.0 * u(i + 1, j, k) - 1.0/6.0 * u(i + 2, j, k);
						
						d0 = 1.0/10.0;
						d1 = 6.0/10.0;
						d2 = 3.0/10.0;
						
						
						diff1 = fabs(u(i, j, k) - u(i - 1, j, k));
						diff2 = fabs(u(i + 1, j, k) - u(i, j, k));
						
						if (diff1 >= diff2){
							diff3 = fabs(abs(u(i + 2, j, k) - u(i + 1, j, k)) - fabs(u(i + 1, j, k) - u(i, j, k)));
							diff4 = fabs(abs(u(i + 1, j, k) - u(i, j, k)) - fabs(u(i, j, k) - u(i - 1, j, k)));
							if (diff3 >= diff4){
								w_xL(i, j, k) = q1L;
								w_xR(i, j, k) = q1R;
							}else{
								w_xL(i, j, k) = q2L;
								w_xR(i, j, k) = q2R;
							}
						}else{
							diff3 = fabs(abs(u(i, j, k) - u(i - 1, j, k)) - fabs(u(i - 1, j, k) - u(i - 2, j, k)));
							diff4 = fabs(abs(u(i + 1, j, k) - u(i, j, k)) - fabs(u(i, j, k) - u(i - 1, j, k)));
							if (diff4 >= diff3){
								w_xL(i, j, k) = q0L;
								w_xR(i, j, k) = q0R;
							}else{
								w_xL(i, j, k) = q1L;
								w_xR(i, j, k) = q1R;
							}
						}
						
						q0L = -1.0/6.0 * u(i, j - 2, k) + 5.0/6.0 * u(i, j - 1, k) + 1.0/3.0 * u(i, j, k);
						q1L = 1.0/3.0 * u(i, j - 1, k) + 5.0/6.0 * u(i, j, k) - 1.0/6.0 * u(i, j + 1, k);
						q2L = 11.0/6.0 * u(i, j, k) - 7.0/6.0 * u(i, j + 1, k) + 1.0/3.0 * u(i, j + 2, k);
						
						
						q0R = 1.0/3.0 * u(i, j - 2, k) - 7.0/6.0 * u(i, j - 1, k) + 11.0/6.0 * u(i, j, k);
						q1R = -1.0/6.0 * u(i, j - 1, k) + 5.0/6.0 * u(i, j, k) + 1.0/3.0 * u(i, j + 1, k);
						q2R = 1.0/3.0 * u(i, j, k) + 5.0/6.0 * u(i, j + 1, k) - 1.0/6.0 * u(i, j + 2, k);
						
						d0 = 1.0/10.0;
						d1 = 6.0/10.0;
						d2 = 3.0/10.0;
						
						
						diff1 = fabs(u(i, j, k) - u(i, j - 1, k));
						diff2 = fabs(u(i, j + 1, k) - u(i, j, k));
						
						if (diff1 >= diff2){
							diff3 = fabs(abs(u(i, j + 2, k) - u(i, j + 1, k)) - fabs(u(i, j + 1, k) - u(i, j, k)));
							diff4 = fabs(abs(u(i, j + 1, k) - u(i, j, k)) - fabs(u(i, j, k) - u(i, j - 1, k)));
							if (diff3 >= diff4){
								w_yL(i, j, k) = q1L;
								w_yR(i, j, k) = q1R;
							}else{
								w_yL(i, j, k) = q2L;
								w_yR(i, j, k) = q2R;
							}
						}else{
							diff3 = fabs(abs(u(i, j, k) - u(i, j - 1, k)) - fabs(u(i, j - 1, k) - u(i, j - 2, k)));
							diff4 = fabs(abs(u(i, j + 1, k) - u(i, j, k)) - fabs(u(i, j, k) - u(i, j - 1, k)));
							if (diff4 >= diff3){
								w_yL(i, j, k) = q0L;
								w_yR(i, j, k) = q0R;
							}else{
								w_yL(i, j, k) = q1L;
								w_yR(i, j, k) = q1R;
							}
						}
						break;
					case WENO5:
						SpaceOrder = 4;
						q0 = -1.0/6.0 * u(i - 2, j, k) + 5.0/6.0 * u(i - 1, j, k) + 1.0/3.0 * u(i, j, k);
						q1 = 1.0/3.0 * u(i - 1, j, k) + 5.0/6.0 * u(i, j, k) - 1.0/6.0 * u(i + 1, j, k);
						q2 = 11.0/6.0 * u(i, j, k) - 7.0/6.0 * u(i + 1, j, k) + 1.0/3.0 * u(i + 2, j, k);
						
						beta0 = 13.0/12.0 * pow((u(i - 2, j, k) - 2 * u(i - 1, j, k) + u(i, j, k)),2.0) + 1.0/4.0 * pow((u(i - 2, j, k) - 4 * u(i - 1, j, k) + 3 * u(i, j, k)),2.0);
						beta1 = 13.0/12.0 * pow((u(i - 1, j, k) - 2 * u(i, j, k) + u(i + 1, j, k)),2.0) + 1.0/4.0 * pow((u(i - 1, j, k) - u(i + 1, j, k)),2.0);
						beta2 = 13.0/12.0 * pow((u(i, j, k) - 2 * u(i + 1, j, k) + u(i + 2, j, k)),2.0) + 1.0/4.0 * pow((3 * u(i, j, k) - 4 * u(i + 1, j, k) + u(i + 2, j, k)),2.0);
						
						d0 = 3.0/10.0;
						d1 = 6.0/10.0;
						d2 = 1.0/10.0;
						
						wr0_bar = d0 / pow((1e-6 + beta0),2.0);
						wr1_bar = d1 / pow((1e-6 + beta1),2.0);
						wr2_bar = d2 / pow((1e-6 + beta2),2.0);
						
						wr0 = wr0_bar / (wr0_bar + wr1_bar + wr2_bar);
						wr1 = wr1_bar / (wr0_bar + wr1_bar + wr2_bar);
						wr2 = wr2_bar / (wr0_bar + wr1_bar + wr2_bar);
						
						w_xL(i, j, k) = wr0 * q0 + wr1 * q1 + wr2 * q2;
						
						q0 = 1.0/3.0 * u(i - 2, j, k) - 7.0/6.0 * u(i - 1, j, k) + 11.0/6.0 * u(i, j, k);
						q1 = -1.0/6.0 * u(i - 1, j, k) + 5.0/6.0 * u(i, j, k) + 1.0/3.0 * u(i + 1, j, k);
						q2 = 1.0/3.0 * u(i, j, k) + 5.0/6.0 * u(i + 1, j, k) - 1.0/6.0 * u(i + 2, j, k);
						
						d0 = 1.0/10.0;
						d1 = 6.0/10.0;
						d2 = 3.0/10.0;
						
						wr0_bar = d0 / pow((1e-6 + beta0),2.0);
						wr1_bar = d1 / pow((1e-6 + beta1),2.0);
						wr2_bar = d2 / pow((1e-6 + beta2),2.0);
						
						wr0 = wr0_bar / (wr0_bar + wr1_bar + wr2_bar);
						wr1 = wr1_bar / (wr0_bar + wr1_bar + wr2_bar);
						wr2 = wr2_bar / (wr0_bar + wr1_bar + wr2_bar);
						
						w_xR(i, j, k) = wr0 * q0 + wr1 * q1 + wr2 * q2;
						
						
						q0 = -1.0/6.0 * u(i, j - 2, k) + 5.0/6.0 * u(i, j - 1, k) + 1.0/3.0 * u(i, j, k);
						q1 = 1.0/3.0 * u(i, j - 1, k) + 5.0/6.0 * u(i, j, k) - 1.0/6.0 * u(i, j + 1, k);
						q2 = 11.0/6.0 * u(i, j, k) - 7.0/6.0 * u(i, j + 1, k) + 1.0/3.0 * u(i, j + 2, k);
						
						beta0 = 13.0/12.0 * pow((u(i, j - 2, k) - 2 * u(i, j - 1, k) + u(i, j, k)),2.0) + 1.0/4.0 * pow((u(i, j - 2, k) - 4 * u(i, j - 1, k) + 3 * u(i, j, k)),2.0);
						beta1 = 13.0/12.0 * pow((u(i, j - 1, k) - 2 * u(i, j, k) + u(i, j + 1, k)),2.0) + 1.0/4.0 * pow((u(i, j - 1, k) - u(i, j + 1, k)),2.0);
						beta2 = 13.0/12.0 * pow((u(i, j, k) - 2 * u(i, j + 1, k) + u(i, j + 2, k)),2.0) + 1.0/4.0 * pow((3 * u(i, j, k) - 4 * u(i, j + 1, k) + u(i, j + 2, k)),2.0);
						
						d0 = 3.0/10.0;
						d1 = 6.0/10.0;
						d2 = 1.0/10.0;
						
						wr0_bar = d0 / pow((1e-6 + beta0),2.0);
						wr1_bar = d1 / pow((1e-6 + beta1),2.0);
						wr2_bar = d2 / pow((1e-6 + beta2),2.0);
						
						wr0 = wr0_bar / (wr0_bar + wr1_bar + wr2_bar);
						wr1 = wr1_bar / (wr0_bar + wr1_bar + wr2_bar);
						wr2 = wr2_bar / (wr0_bar + wr1_bar + wr2_bar);
						
						w_yL(i, j, k) = wr0 * q0 + wr1 * q1 + wr2 * q2;
						
						q0 = 1.0/3.0 * u(i, j - 2, k) - 7.0/6.0 * u(i, j - 1, k) + 11.0/6.0 * u(i, j, k);
						q1 = -1.0/6.0 * u(i, j - 1, k) + 5.0/6.0 * u(i, j, k) + 1.0/3.0 * u(i, j + 1, k);
						q2 = 1.0/3.0 * u(i, j, k) + 5.0/6.0 * u(i, j + 1, k) - 1.0/6.0 * u(i, j + 2, k);
						
						d0 = 1.0/10.0;
						d1 = 6.0/10.0;
						d2 = 3.0/10.0;
						
						wr0_bar = d0 / pow((1e-6 + beta0),2.0);
						wr1_bar = d1 / pow((1e-6 + beta1),2.0);
						wr2_bar = d2 / pow((1e-6 + beta2),2.0);
						
						wr0 = wr0_bar / (wr0_bar + wr1_bar + wr2_bar);
						wr1 = wr1_bar / (wr0_bar + wr1_bar + wr2_bar);
						wr2 = wr2_bar / (wr0_bar + wr1_bar + wr2_bar);
						
						w_yR(i, j, k) = wr0 * q0 + wr1 * q1 + wr2 * q2;
						break;
				}
			}
		}
	}
	
	switch (reconstate) {
		case ConsState:
			break;
		case PrimState:
			upri2ucon(w_xL, w_xL);
			upri2ucon(w_xR, w_xR);
			upri2ucon(w_yL, w_yL);
			upri2ucon(w_yR, w_yR);
			break;
	}
	
	bound_cond_u(w_xL);
	bound_cond_u(w_xR);
	bound_cond_u(w_yL);
	bound_cond_u(w_yR);
}



double calecon(double *con){
	double rho = con[0];
	double vx = con[1] / con[0];
	double vy = con[2] / con[0];
	double vz = con[3] / con[0];
	double U = con[4];
	double Bx = con[5];
	double By = con[6];
	double Bz = con[7];
	double B_square = Bx * Bx + By * By + Bz * Bz;
	double v_square = vx * vx + vy * vy + vz * vz;
	double e = (U - 0.5 * rho * v_square - 0.5 * B_square)/rho;
	return e;
}

double calepri(double *pri){
	double e = eos_e_from_rho_p(pri[0], pri[4]);
	return e;
}

double calcf(double* pri){
	double Bx = pri[5];
	double By = pri[6];
	double Bz = pri[7];
	double cs = eos_sound_a(pri[0], pri[4]);
	double ca_x = fabs(Bx * Bx + By * By + Bz * Bz) / sqrt(pri[0]);
	double ca_y = fabs(Bx * Bx + By * By + Bz * Bz) / sqrt(pri[0]);
	double ca_z = fabs(Bx * Bx + By * By + Bz * Bz) / sqrt(pri[0]);
	double cf_x = sqrt(0.5 * ((cs * cs + ca_x * ca_x) + sqrt((cs * cs + ca_x * ca_x) * (cs * cs + ca_x * ca_x) - ((4 * (cs * cs * Bx * Bx))/pri[0]))));
	double cf_y = sqrt(0.5 * ((cs * cs + ca_y * ca_y) + sqrt((cs * cs + ca_y * ca_y) * (cs * cs + ca_y * ca_y) - ((4 * (cs * cs * By * By))/pri[0]))));
	double cf_z = sqrt(0.5 * ((cs * cs + ca_z * ca_z) + sqrt((cs * cs + ca_z * ca_z) * (cs * cs + ca_z * ca_z) - ((4 * (cs * cs * Bz * Bz))/pri[0]))));
	double cf = max({cf_x, cf_y, cf_z});
	return cf;
}

vector<double> ConvergenceTestFunction(double x){
	vector<double> f(8);
	f[0] = 1;
	f[1] = 0;
	f[2] = -0.2 * cos(2 * M_PI * x);
	f[3] = 0.2 * sin(2 * M_PI * x);
	f[4] = 1;
	f[5] = 1;
	f[6] = 0.2 * cos(2 * M_PI * x);
	f[7] = -0.2 * sin(2 * M_PI * x);
	return f;
}

vector<double> ConvergenceTestExact(double x){
	vector<double> f(NUM_VARS);
	f[0] = 1;
	f[1] = 0;
	f[2] = -0.2 * cos(2 * M_PI * x);
	f[3] = 0.2 * sin(2 * M_PI * x);
	f[4] = 1;
	f[5] = 1;
	f[6] = 0.2 * cos(2 * M_PI * x);
	f[7] = -0.2 * sin(2 * M_PI * x);
	return f;
}

vector<double> GaussLegendreNumericalIntegration(const double& xc, const double& dx, function<vector<double>(double)> func, int nGaussPoints){
	vector<double> Weights(nGaussPoints);
	vector<double> Location(nGaussPoints);
	
	// Select Gauss-Legendre nodes & weights
	switch(nGaussPoints){
		case 1:
			Location = {0.0};
			Weights  = {1.0};
			break;
		case 2:
			Location = {-dx/2.0/sqrt(3.0), dx/2.0/sqrt(3.0)};
			Weights  = {0.5, 0.5};
			break;
		case 3:
			Location = {-dx/2.0*sqrt(3.0/5.0), 0.0, dx/2.0*sqrt(3.0/5.0)};
			Weights  = {5.0/18.0, 8.0/18.0, 5.0/18.0};
			break;
		case 4:
			Location = {-dx/2.0*0.8611363115940526,
				-dx/2.0*0.3399810435848563,
				+dx/2.0*0.3399810435848563,
				+dx/2.0*0.8611363115940526};
			
			Weights  = {0.5*0.3478548451374538,
				0.5*0.6521451548625461,
				0.5*0.6521451548625461,
				0.5*0.3478548451374538};
			break;
		case 5:
			Location = {-dx/2.0*0.9061798459386640,
				-dx/2.0*0.5384693101056831,
				0.0,
				+dx/2.0*0.5384693101056831,
				+dx/2.0*0.9061798459386640};
			
			Weights  = {0.5*0.2369268850561891,
				0.5*0.4786286704993665,
				0.5*0.5688888888888889,
				0.5*0.4786286704993665,
				0.5*0.2369268850561891};
			break;
		case 6:
			Location = {-dx/2.0*0.9324695142031521,
				-dx/2.0*0.6612093864662645,
				-dx/2.0*0.2386191860831969,
				+dx/2.0*0.2386191860831969,
				+dx/2.0*0.6612093864662645,
				+dx/2.0*0.9324695142031521};
			
			Weights  = {0.5*0.1713244923791704,
				0.5*0.3607615730481386,
				0.5*0.4679139345726910,
				0.5*0.4679139345726910,
				0.5*0.3607615730481386,
				0.5*0.1713244923791704};
			break;
		default:
			std::cerr << "Number of Gauss-Legendre nodes not implemented\n";
			std::exit(EXIT_FAILURE);
	}
	vector<double> integral;
	integral.resize(func(xc).size(), 0.0);
	
	for(int i = 0; i < nGaussPoints; i++){
		double xGauss = xc + Location[i];
		double w      = Weights[i];
		vector<double> f = func(xGauss);
		for(size_t k = 0; k < f.size(); k++){
			integral[k] += w * f[k];
		}
	}
	
	return integral;
}
