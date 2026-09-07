// Initial data for every test case in the Test enum.

#include "declarations.hpp"

void u_initialization(Grid &u, Array &Bx_p0, Array &By_0p){
	vector<double> result;
	// 分界半径：你需要指定;
	// 平滑厚度：你需要指定;
	double eta;
	vector<double> Weights(2);
	vector<double> Location(2);
	double Mx;
	for(int i = ghost_cells; i < nxCells + ghost_cells; i++){
		for(int j = ghost_cells; j < nyCells + ghost_cells; j++){
			double x = x_0 + (i - ghost_cells + 0.5) * dx;
			double y = y_0 + (j - ghost_cells + 0.5) * dy;
			switch (test) {
				case test0:
					if(x < 0){
						u(i,j,0) = 1;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 1;
						u(i, j, 5) = 0;
						u(i, j, 6) = 0;
						u(i, j, 7) = 0;
					}else{
						u(i,j,0) = 0.125;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 0.1;
						u(i, j, 5) = 0;
						u(i, j, 6) = 0;
						u(i, j, 7) = 0;
					}
					break;
				case test1:
					if(x < 0){
						u(i,j,0) = 1;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 1;
						u(i, j, 5) = 0.75;
						u(i, j, 6) = 1;
						u(i, j, 7) = 0;
					}else{
						u(i,j,0) = 0.125;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 0.1;
						u(i, j, 5) = 0.75;
						u(i, j, 6) = -1;
						u(i, j, 7) = 0;
					}
					break;
				case test1y:
					if(y < 0){
						u(i,j,0) = 1;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 1;
						u(i, j, 5) = 1;
						u(i, j, 6) = 0.75;
						u(i, j, 7) = 0;
					}else{
						u(i,j,0) = 0.125;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 0.1;
						u(i, j, 5) = -1;
						u(i, j, 6) = 0.75;
						u(i, j, 7) = 0;
					}
					break;
				case test2:
					if(x < -0.1){
						u(i,j,0) = 1.08;
						u(i,j,1) = 1.2;
						u(i,j,2) = 0.01;
						u(i,j,3) = 0.5;
						u(i,j,4) = 0.95;
						u(i, j, 5) = 0.564190;
						u(i, j, 6) = 1.015541;
						u(i, j, 7) = 0.564190;
					}else{
						u(i,j,0) = 0.9891;
						u(i,j,1) = -0.0131;
						u(i,j,2) = 0.0269;
						u(i,j,3) = 0.010037;
						u(i,j,4) =  0.97159;
						u(i, j, 5) = 0.564190;
						u(i, j, 6) = 1.135262;
						u(i, j, 7) = 0.564923;
					}
					break;
				case test3:
					if(x < -0.1){
						u(i,j,0) = 1.7;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 1.7;
						u(i, j, 5) = 1.1;
						u(i, j, 6) = 1;
						u(i, j, 7) = 0;
					}else{
						u(i,j,0) = 0.2;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = -1.496891;
						u(i,j,4) =  0.2;
						u(i, j, 5) = 1.1;
						u(i, j, 6) = 0.785887;
						u(i, j, 7) = 0.618370;
					}
					break;
				case test4:
					if(x < 0){
						u(i,j,0) = 1;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 1;
						u(i, j, 5) = 1.3;
						u(i, j, 6) = 1;
						u(i, j, 7) = 0;
					}else{
						u(i,j,0) = 0.4;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) =  0.4;
						u(i, j, 5) = 1.3;
						u(i, j, 6) = -1;
						u(i, j, 7) = 0;
					}
					break;
				case test5:
					if(x < 0){
						u(i,j,0) = 0.15;
						u(i,j,1) = 21.55;
						u(i,j,2) = 1.0;
						u(i,j,3) = 1.0;
						u(i,j,4) = 0.28;
						u(i, j, 5) = 0.05 / sqrt(4 * M_PI);
						u(i, j, 6) = -2.0 / sqrt(4 * M_PI);
						u(i, j, 7) = -1.0 / sqrt(4 * M_PI);
					}else{
						u(i,j,0) = 0.1;
						u(i,j,1) = -26.45;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) =  0.1;
						u(i, j, 5) = 0.05 / sqrt(4 * M_PI);
						u(i, j, 6) = 2.0 / sqrt(4 * M_PI);
						u(i, j, 7) = 1.0 / sqrt(4 * M_PI);
					}
					break;
				case test6:
					if(x < 0){
						u(i,j,0) = 1.0;
						u(i,j,1) = 36.87;
						u(i,j,2) = -0.115;
						u(i,j,3) = -0.0386;
						u(i,j,4) = 1.0;
						u(i, j, 5) = 4.0 / sqrt(4 * M_PI);
						u(i, j, 6) = 4.0 / sqrt(4 * M_PI);
						u(i, j, 7) = 1.0 / sqrt(4 * M_PI);
					}else{
						u(i,j,0) = 1;
						u(i,j,1) = -36.87;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 1;
						u(i, j, 5) = 4 / sqrt(4 * M_PI);
						u(i, j, 6) = 4 / sqrt(4 * M_PI);
						u(i, j, 7) = 1 / sqrt(4 * M_PI);
					}
					break;
				case test7:
					if(x < 0){
						u(i,j,0) = 1.0 / (4 * M_PI);
						u(i,j,1) = -1;
						u(i,j,2) = 1;
						u(i,j,3) = -1;
						u(i,j,4) = 1.0;
						u(i, j, 5) = 1.0 / sqrt(4 * M_PI);
						u(i, j, 6) = -1.0 / sqrt(4 * M_PI);
						u(i, j, 7) = 1.0 / sqrt(4 * M_PI);
					}else{
						u(i,j,0) = 1 / (4 * M_PI);
						u(i,j,1) = -1;
						u(i,j,2) = -1;
						u(i,j,3) = -1;
						u(i,j,4) = 1;
						u(i, j, 5) = 1 / sqrt(4 * M_PI);
						u(i, j, 6) = 1 / sqrt(4 * M_PI);
						u(i, j, 7) = 1 / sqrt(4 * M_PI);
					}
					break;
				case test8:
					if(x < 0){
						u(i,j,0) = 1.0;
						u(i,j,1) = -2.0;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 0.45;
						u(i, j, 5) = 0;
						u(i, j, 6) = 0.5;
						u(i, j, 7) = 0;
					}else{
						u(i,j,0) = 1.0;
						u(i,j,1) = 2.0;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 0.45;
						u(i, j, 5) = 0;
						u(i, j, 6) = 0.5;
						u(i, j, 7) = 0;
					}
					{
						// 抹宽约 4 个 cell(δ 越大越平滑);
						double vx  = 1 * sqrt(dx*dx);     // -2 -> +2 平滑过渡,去掉 x=0 速度间断
						u(i,j,0) = 1.0;
						u(i,j,1) = vx;
						u(i,j,2) = 0; u(i,j,3) = 0;
						u(i,j,4) = 0.45;
						u(i,j,5) = 0; u(i,j,6) = 0.5; u(i,j,7) = 0;
						break;
					}
					break;
				case test9:
					if(x < 0){
						u(i,j,0) = 1.0;
						u(i,j,1) = 10;
						u(i,j,2) = 1;
						u(i,j,3) = -1;
						u(i,j,4) = 20.0;
						u(i, j, 5) = 5.0 / sqrt(4 * M_PI);
						u(i, j, 6) = 5.0 / sqrt(4 * M_PI);
						u(i, j, 7) = 0;
					}else{
						u(i,j,0) = 1.0;
						u(i,j,1) = -10;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 1.0;
						u(i, j, 5) = 5 / sqrt(4 * M_PI);
						u(i, j, 6) = 5 / sqrt(4 * M_PI);
						u(i, j, 7) = 0;
					}
					break;
				case test10:
					if(x < 0){
						u(i,j,0) = 1.0;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 1.0;
						u(i, j, 5) = 0.7;
						u(i, j, 6) = 0;
						u(i, j, 7) = 0;
					}else{
						u(i,j,0) = 0.3;
						u(i,j,1) = 0.0;
						u(i,j,2) = 0.0;
						u(i,j,3) = 1.0;
						u(i,j,4) = 0.2;
						u(i, j, 5) = 0.7;
						u(i, j, 6) = 1.0;
						u(i, j, 7) = 0;
					}
					break;
				case RKEOS1:
					if(x < 0){
						u(i,j,0) = 1.0;
						u(i,j,1) = 1.0;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 2.0;
						u(i, j, 5) = 0.0;
						u(i, j, 6) = 0;
						u(i, j, 7) = 0;
					}else{
						u(i,j,0) = 1.0;
						u(i,j,1) = -1.0;
						u(i,j,2) = 0.0;
						u(i,j,3) = 0.0;
						u(i,j,4) = 1.0;
						u(i, j, 5) = 0.0;
						u(i, j, 6) = 0.0;
						u(i, j, 7) = 0.0;
					}
					break;
				case RKEOS2:
					if(x < 0){
						u(i,j,0) = 1.0;
						u(i,j,1) = 0.0;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 1.0;
						u(i, j, 5) = 0.0;
						u(i, j, 6) = 0;
						u(i, j, 7) = 0;
					}else{
						u(i,j,0) = 0.125;
						u(i,j,1) = 0.0;
						u(i,j,2) = 0.0;
						u(i,j,3) = 0.0;
						u(i,j,4) = 0.1;
						u(i, j, 5) = 0.0;
						u(i, j, 6) = 0.0;
						u(i, j, 7) = 0.0;
					}
					break;
				case RKEOS3:
					if(x < 0.1){
						u(i,j,0) = 1.0;
						u(i,j,1) = 0.0;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 1000;
						u(i, j, 5) = 0.0;
						u(i, j, 6) = 0;
						u(i, j, 7) = 0;
					}else{
						u(i,j,0) = 1.0;
						u(i,j,1) = 0.0;
						u(i,j,2) = 0.0;
						u(i,j,3) = 0.0;
						u(i,j,4) = 0.01;
						u(i, j, 5) = 0.0;
						u(i, j, 6) = 0.0;
						u(i, j, 7) = 0.0;
					}
					break;
				case RKEOS4:
					if(x < 0){
						u(i,j,0) = 1.0;
						u(i,j,1) = 0.0;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 2.0;
						u(i, j, 5) = 0.0;
						u(i, j, 6) = 0;
						u(i, j, 7) = 0;
					}else{
						u(i,j,0) = 1.5;
						u(i,j,1) = 0.0;
						u(i,j,2) = 0.0;
						u(i,j,3) = 0.0;
						u(i,j,4) = 1.0;
						u(i, j, 5) = 0.0;
						u(i, j, 6) = 0.0;
						u(i, j, 7) = 0.0;
					}
					break;
				case shocktest:
					if(x < 0){
						u(i,j,0) = 1;
						u(i,j,1) = 2;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 0.1;
						u(i, j, 5) = 0;
						u(i, j, 6) = 0;
						u(i, j, 7) = 0;
					}else{
						u(i,j,0) = 1;
						u(i,j,1) = -2;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 0.1;
						u(i, j, 5) = 0;
						u(i, j, 6) = 0;
						u(i, j, 7) = 0;
					}
					break;
				case convergencetest:
					//					u(i,j,0) = 1;
					//					u(i,j,1) = 0;
					//					u(i,j,2) = -0.2 * cos(2 * M_PI * x);
					//					u(i,j,3) = 0.2 * sin(2 * M_PI * x);
					//					u(i,j,4) = 1;
					//					u(i, j, 5) = 1;
					//					u(i, j, 6) = 0.2 * cos(2 * M_PI * x);
					//					u(i, j, 7) = -0.2 * sin(2 * M_PI * x);
					
					result = GaussLegendreNumericalIntegration(x, dx, ConvergenceTestFunction, 2 );
					u(i,j,0) = result[0];
					u(i,j,1) = result[1];
					u(i,j,2) = result[2];
					u(i,j,3) = result[3];
					u(i,j,4) = result[4];
					u(i, j, 5) = result[5];
					u(i, j, 6) = result[6];
					u(i, j, 7) = result[7];
					break;
				case test2d0:
					if(sqrt(x * x + y * y) <= 0.4){
						u(i,j,0) = 1;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 1;
						u(i, j, 5) = 0;
						u(i, j, 6) = 0;
						u(i, j, 7) = 0;
					}else{
						u(i,j,0) = 0.125;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 0.1;
						u(i, j, 5) = 0;
						u(i, j, 6) = 0;
						u(i, j, 7) = 0;
					}
					break;
				case testOrszag_Tang:
					u(i, j, 0) = gamma_plasma * gamma_plasma;
					u(i, j, 1) = -sin(2 * M_PI * y);
					u(i, j, 2) = sin(2 * M_PI * x);
					u(i, j, 3) = 0;
					u(i, j, 4) = gamma_plasma;
					u(i, j, 5) = -sin(2 * M_PI * y);
					u(i, j, 6) = sin(4 * M_PI * x);
					u(i, j, 7) = 0;
					break;
				case testblastwave: {
					double rad = sqrt(x*x + y*y);
					double rD = 0.1;
					double prsL = 1000.0, prsR = 0.1;
					double eps = 1 * sqrt(dx*dx + dy*dy);   // smoothing width
					
					u(i,j,0) = 1.0;
					u(i,j,1) = 0; u(i,j,2) = 0; u(i,j,3) = 0;
					u(i,j,4) = 0.5*(prsR+prsL) + 0.5*(prsR-prsL) * erf((rad - rD)/eps);
					u(i,j,5) = 70.0 / sqrt(4.0 * M_PI);
					u(i,j,6) = 0.0;
					u(i,j,7) = 0.0;
					break;
				}
				case testblastwave2: {
					double rad = sqrt(x*x + y*y);
					double rD = 0.1;
					double prsL = 10000.0, prsR = 0.1;
					double eps = 1 * sqrt(dx*dx + dy*dy);   // smoothing width
					
					u(i,j,0) = 1.0;
					u(i,j,1) = 0; u(i,j,2) = 0; u(i,j,3) = 0;
					u(i,j,4) = 0.5*(prsR+prsL) + 0.5*(prsR-prsL) * erf((rad - rD)/eps);
					u(i,j,5) = 1000.0 / sqrt(4.0 * M_PI);
					u(i,j,6) = 0.0;
					u(i,j,7) = 0.0;
					break;
				}
				case testShockCloud:{
					if(x < 0.6){
						u(i,j,0) = 3.86859;
						u(i,j,1) = 0; u(i,j,2) = 0; u(i,j,3) = 0;
						u(i,j,4) = 167.345;
						u(i,j,5) = 0;
						u(i,j,6) = 2.1826182;
						u(i,j,7) = -2.1826182;
					}else{
						u(i,j,0) = 1;
						u(i,j,1) = -11.2536; u(i,j,2) = 0; u(i,j,3) = 0;
						u(i,j,4) = 1;
						u(i,j,5) = 0;
						u(i,j,6) = 0.56418958;
						u(i,j,7) = 0.56418958;
						if(sqrt((x - 0.8) * (x - 0.8) + (y - 0.5) * (y - 0.5)) < 0.15){
							u(i,j,0) = 10;
						}
					}
					break;
				}
				case testblastwave1D:					
					if(fabs(x) <= 0.1){
						u(i,j,0) = 1;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 1000;
						u(i, j, 5) = 1 / sqrt(4 * M_PI) * 70;
						u(i, j, 6) = 0;
						u(i, j, 7) = 0;
					}else{
						u(i,j,0) = 1;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 0.1;
						u(i, j, 5) = 1 / sqrt(4 * M_PI) * 70;
						u(i, j, 6) = 0;
						u(i, j, 7) = 0;
					}
//					
					break;
				case testRK_Riemann:
					// Redlich-Kwong EOS MHD shock tube(对应 1D 代码 test8 / 论文 4.3.1),x 方向 Riemann
					if(x <= 0){
						u(i,j,0) = 1;
						u(i,j,1) = 1;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 2;
						u(i,j,5) = 0;
						u(i,j,6) = 0;
						u(i,j,7) = 0;
					}else{
						u(i,j,0) = 1.0;
						u(i,j,1) = -1.0;
						u(i,j,2) = 0.0;
						u(i,j,3) = 0.0;
						u(i,j,4) = 1.0;
						u(i,j,5) = 0;
						u(i,j,6) = 0;
						u(i,j,7) = 0;
					}
					break;
				case testRotor:
					if(sqrt((x) * (x) + (y) * (y)) <= 0.1){
						u(i,j,0) = 10;
						u(i,j,1) = -10 * y;
						u(i,j,2) = 10 * x;
						u(i,j,3) = 0;
						u(i,j,4) = 1;
						u(i, j, 5) = 2.5 / sqrt(4 * M_PI);
						u(i, j, 6) = 0;
						u(i, j, 7) = 0;
					}else{
						u(i,j,0) = 1;
						u(i,j,1) = 0;
						u(i,j,2) = 0;
						u(i,j,3) = 0;
						u(i,j,4) = 1;
						u(i, j, 5) = 2.5 / sqrt(4 * M_PI);
						u(i, j, 6) = 0;
						u(i, j, 7) = 0;
					}
					
					break;
				case testFieldLoop:
					if(sqrt(x * x + y * y) < 0.3){
						u(i,j,0) = 1;
						u(i,j,1) = 2;
						u(i,j,2) = 1;
						u(i,j,3) = 0;
						u(i,j,4) = 1e5;
						u(i, j, 5) = (-A0 * y / sqrt(x * x + y * y) ) / sqrt(4 * M_PI);
						u(i, j, 6) = (A0 * x / sqrt(x * x + y * y) ) / sqrt(4 * M_PI);
						u(i, j, 7) = 0;
					}else{
						u(i,j,0) = 1;
						u(i,j,1) = 2;
						u(i,j,2) = 1;
						u(i,j,3) = 0;
						u(i,j,4) = 1e5;
						u(i, j, 5) = 0;
						u(i, j, 6) = 0;
						u(i, j, 7) = 0;
					}
					break;
				case testKHInstability:
					Mx = A0;
					if (y >= -9.0/32.0 && y < -7.0/32.0) {
						eta = 0.5 * (1.0 + sin(16.0 * M_PI * (y + 0.25)));
					}
					else if (y >= -7.0/32.0 && y < 7.0/32.0) {
						eta = 1.0;
					}
					else if (y >= 7.0/32.0 && y < 9.0/32.0) {
						eta = 0.5 * (1.0 - sin(16.0 * M_PI * (y - 0.25)));
					}
					else {
						eta = 0.0;
					}
					u(i,j,0) = 1.4;                         // rho = gamma = 1.4
					u(i,j,1) = Mx * (1.0 - 2.0 * eta);        // u = Mx(1 - 2 eta(y))
					u(i,j,2) = 0.1 * Mx * sin(2.0 * M_PI * x);// v = 0.1 Mx sin(2 pi x)
					u(i,j,3) = 0.0;                           // w = 0
					u(i,j,4) = 1.0;                           // p = 1
					u(i,j,5) = 0.1 * Mx;                      // Bx = 0.1 Mx (uniform)
					u(i,j,6) = 0.0;                           // By = 0
					u(i,j,7) = 0.0;   
					break;
				case testdoubleshear:
					u(i,j,0) = 1;
					if(y < 0.5){
						u(i,j,1) = tanh(30 * (y - 0.25));
					}else{
						u(i,j,1) = tanh(30 * (0.75 - y));
					}
					u(i,j,2) = 0.05 * sin(2 * M_PI * x);
					u(i,j,4) = 1e5 / gamma_plasma;
					break;
				case testPlasmaEquilibrium:
					{
						double B0 = 8.0, R = 0.9;
						double r  = sqrt(x*x + y*y);
						double Bz = (r <= R) ? B0 * (r*r)/(R*R) : B0;     // 纵向场 Bz(r)
						double p  = (r <= R) ? 0.5*(B0*B0 - Bz*Bz) + 0.1  // 力平衡压力
						: 0.1;
						
						u(i,j,0) = 1.0;          // rho
						u(i,j,1) = 0;         // u  (vx)
						u(i,j,2) = 0;         // v  (vy)
						u(i,j,3) = 0.0;          // w
						u(i,j,4) = p;            // p
						u(i,j,5) = 0.0;          // Bx
						u(i,j,6) = 0.0;          // By
						u(i,j,7) = Bz;           // Bz  ← 场全在 z 方向
					}
					break;
				case testAdvectedPlasmaEquilibrium:
					{
						double B0 = 8.0, R = 0.9;
						double r  = sqrt(x*x + y*y);
						double Bz = (r <= R) ? B0 * (r*r)/(R*R) : B0;     // 纵向场 Bz(r)
						double p  = (r <= R) ? 0.5*(B0*B0 - Bz*Bz) + 0.1  // 力平衡压力
						: 0.1;
						
						u(i,j,0) = 1.0;          // rho
						u(i,j,1) = 0.04;         // u  (vx)
						u(i,j,2) = 0.04;         // v  (vy)
						u(i,j,3) = 0.0;          // w
						u(i,j,4) = p;            // p
						u(i,j,5) = 0.0;          // Bx
						u(i,j,6) = 0.0;          // By
						u(i,j,7) = Bz;           // Bz  ← 场全在 z 方向
					}
					break;
				case testMHD_vortex:{
					const double eps      = 1.0;              // 论文里定义的 ε
					const double mu_tilde = sqrt(4.0 * M_PI);  // 论文里定义的 μ̃
					
					double r2 = x*x + y*y;
					
					double expHalf = exp(0.5 * (1.0 - r2));  // e^{0.5(1-r^2)}，用于 u, v, Bx, By
					double expFull = exp(1.0 - r2);          // e^{(1-r^2)}，用于压强磁场修正项
					
					double rho0 = 1.0;
					double u0   = -eps / (2.0*M_PI) * expHalf * y;
					double v0   =  eps / (2.0*M_PI) * expHalf * x;
					double Bx0  = -1.0 / (2.0*M_PI) * expHalf * y;
					double By0  =  1.0 / (2.0*M_PI) * expHalf * x;
					
					double p0 = 1.0
					+ (1.0/(8.0*M_PI)) * (mu_tilde/(2.0*M_PI)) * (mu_tilde/(2.0*M_PI)) * (1.0 - r2) * expFull
					- 0.5 * (eps/(2.0*M_PI)) * (eps/(2.0*M_PI)) * expFull;
					
					u(i,j,0) = rho0;
					u(i,j,1) = u0;
					u(i,j,2) = v0;
					u(i,j,3) = 0.0;
					u(i,j,4) = p0;
					u(i,j,5) = Bx0;
					u(i,j,6) = By0;
					u(i,j,7) = 0.0;
					break;
				}

				case testMHDJet:{
					// ambient: (rho, v, B, p) = (0.1*gamma, 0,0,0, 0,B0,0, 1)
					u(i,j,0) = 0.1 * gamma_plasma;
					u(i,j,1) = 0.0;
					u(i,j,2) = 0.0;
					u(i,j,3) = 0.0;
					u(i,j,4) = 1.0;
					u(i,j,5) = 0.0;
					u(i,j,6) = jet_B0;
					u(i,j,7) = 0.0;
					break;
				}

			}
		}
	}
	// 此处 u 为原始量 (main 里传入 u_pri):jet 的 Dirichlet ghost 必须按原始量填 (slot4=p)
	if(y_bc == y_JetInflow){ g_bc_prim = true; bound_cond_u(u); g_bc_prim = false; }
	else bound_cond_u(u);
	
//	for(int i = 0; i < nxCells + 2 * ghost_cells - 1; i++){
//		for(int j = 0; j < nyCells + 2 * ghost_cells - 1; j++){
//			Bx_p0(i,j) = 0.5 * (u(i,j,5) + u(i + 1,j,5));
//			By_0p(i,j) = 0.5 * (u(i,j,6) + u(i,j + 1,6));
//		}
//	}
	// ↓↓↓ 从这里开始替换(原 2478–2483 行的 0.5*(...) 平均循环)↓↓↓
	if(test == testMHD_vortex){
		auto Az_corner = [](int i, int j) {
			double x = x_0 + (i - ghost_cells + 1.0) * dx;
			double y = y_0 + (j - ghost_cells + 1.0) * dy;
			double r2 = x*x + y*y;
			return (1.0 / (2.0*M_PI)) * exp(0.5 * (1.0 - r2));
		};
		for(int i = 0; i < nxCells + 2*ghost_cells; i++)
			for(int j = 1; j < nyCells + 2*ghost_cells; j++)
			Bx_p0(i,j) =  (Az_corner(i,j) - Az_corner(i,j-1)) / dy;
		for(int i = 1; i < nxCells + 2*ghost_cells; i++)
			for(int j = 0; j < nyCells + 2*ghost_cells; j++)
			By_0p(i,j) = -(Az_corner(i,j) - Az_corner(i-1,j)) / dx;
		// 胞心 B 改为面平均,与 CT 口径一致
		for(int i = 1; i < nxCells + 2*ghost_cells; i++)
			for(int j = 1; j < nyCells + 2*ghost_cells; j++){
				u(i,j,5) = 0.5 * (Bx_p0(i-1,j) + Bx_p0(i,j));
				u(i,j,6) = 0.5 * (By_0p(i,j-1) + By_0p(i,j));
			}
		bound_cond_u(u);
	}else{
		// 其他测试暂时保留原来的胞心平均
		for(int i = 0; i < nxCells + 2 * ghost_cells - 1; i++)
			for(int j = 0; j < nyCells + 2 * ghost_cells - 1; j++){
				Bx_p0(i,j) = 0.5 * (u(i,j,5) + u(i + 1,j,5));
				By_0p(i,j) = 0.5 * (u(i,j,6) + u(i,j + 1,6));
			}
	}
	// ↑↑↑ 替换到这里 ↑↑↑
	
	

	bound_cond_staggered(Bx_p0, 0);   // x-normal: 保留边界法向面(CT div-free)
	bound_cond_staggered(By_0p, 1); 
}
