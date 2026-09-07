// Intermediate plotting and VTK output.

#include "declarations.hpp"

void IntermediatePlotIntermediatePlot(Grid& u_con, Grid &u_pri,  const int iteration, const int plotINTERVAL, const std::vector<double>& IntermediatePlotTimes, bool isFinalTime){
	bool isPlotTime = false;
	std::string testname;
	switch(test){
		case test0: 
			testname = "test0";
			break;
		case test1: 
			testname = "test1";
			break;
		case test2: 
			testname = "test2";
			break;
		case test3: 
			testname = "test3";
			break;
		case test4: 
			testname = "test4";
			break;
		case test5:
			testname = "test5";
			break;
		case test6:
			testname = "test6";
			break;
		case test7:
			testname = "test7";
			break;
		case test8:
			testname = "test8";
			break;
		case test9:
			testname = "test9";
			break;
		case test10:
			testname = "test10";
			break;
		case RKEOS1:
			testname = "RKEOS1";
			break;
		case RKEOS2:
			testname = "RKEOS2";
			break;
		case RKEOS3:
			testname = "RKEOS3";
			break;
		case RKEOS4:
			testname = "RKEOS4";
			break;
		case shocktest:
			testname = "shocktest";
			break;
		case convergencetest:
			testname = "convergencetest";
			break;
		case test2d0:
			testname = "test2d0";
			break;
		case test1y:
			testname = "test1y";
			break;
		case testOrszag_Tang:
			testname = "Orszag_Tang";
			break;
		case testblastwave:
			testname = "blast_wave";
			break;
		case testblastwave2:
			testname = "blast_wave_2";
			break;
		case testMHDJet:
			testname = std::string("MHD_jet") + " JetB0 = "+std::to_string(jet_B0*jet_B0) + " JetM0 = " + std::to_string(jet_v);
			break;
		case testShockCloud:
			testname = "testShockCloud";
			break;
		case testblastwave1D:
			testname = "blast_wave1D";
			break;
		case testRotor:
			testname = "Rotor";
			break;
		case testFieldLoop:
			testname = "FieldLoop" + std::to_string(A0) ;
			break;
		case testKHInstability:
			testname = "KHInstability" + std::to_string(A0);
			break;
		case testdoubleshear:
			testname = "doubleshear";
			break;
		case testRK_Riemann:
			testname = "RK_Riemann";
			break;
		case testPlasmaEquilibrium:
			testname = "PlasmaEquilibrium";
			break;
		case testAdvectedPlasmaEquilibrium:
			testname = "AdvectedPlasmaEquilibrium";
			break;
		case testMHD_vortex:
			testname = "MHD_vortex";
			break;
	}
	
	for(int k=0; k<(int)IntermediatePlotTimes.size(); k++){
		double Tdiff = std::abs(IntermediatePlotTimes[k]-t);
		if(Tdiff <= dt/2.0){isPlotTime = true; break;}
	}
	if(iteration == 0) isPlotTime = true;
	if(iteration % plotINTERVAL == 0) isPlotTime = true;
	if(isPlotTime || isFinalTime){
		std::cout<<"Intermediate plot"<< std::endl;
		printAllVariablesToVTK(u_con,u_pri, testname,iteration,isFinalTime);
	}
}

/* VISIT OUTPUT FILE CREATION */

/*Prints the header section for a VTK file. */


void printVTKHeader(std::ofstream& outFile, int& nstep){
	outFile << "# vtk DataFile Version 3.0" << std::endl;
	outFile << "vtk output" << std::endl;
	outFile << "ASCII" << std::endl;
	outFile << "DATASET RECTILINEAR_GRID" << std::endl;
	outFile << "FIELD FieldData 2" << std::endl;
	outFile << "CYCLE 1 1 int" << std::endl;
	outFile << nstep << std::endl;
	outFile << "TIME 1 1 double" << std::endl;
	outFile << t << std::endl;
	
	outFile << "DIMENSIONS "<< nxCells + 1 << " " << nyCells + 1 << " 1" << std::endl;
	
	outFile << "X_COORDINATES " << nxCells + 1 <<" double" << std::endl;
	
	for(int i=ghost_cells; i<nxCells + ghost_cells ; i++){
		outFile << x_0+(((double)i)-ghost_cells +0.5)*dx << " ";
	}
	
	outFile <<  std::endl;
	
	outFile << "Y_COORDINATES " << nyCells + 1 <<" double" << std::endl;
	
	for(int j=ghost_cells; j<nyCells+ghost_cells ; j++){
		outFile << y_0 + (((double)j)-ghost_cells + 0.5 )*dy << " ";
	}
	
	outFile << std::endl;
	
	outFile << "Z_COORDINATES 1 double" << std::endl;
	
	outFile << "0" << std::endl;
	
	outFile << "CELL_DATA " << nxCells*nyCells << std::endl;
	
	return;
}

/** Prints a scalar variable into a VTK file.
*/
void printVTKVariable(std::ofstream& outFile, Grid& u_con, Grid& u_pri, int var, std::string variableName){
	outFile << variableName << " 1 " << nxCells * nyCells << " double" << std::endl;
	
	int printWidth = 24;
	
	for(int j = 0; j < nyCells; j++){
		for(int i = 0; i < nxCells; i++){
			double varPrint;
			if(var < NUM_VARS){
				varPrint = u_con(i,j,var);
			}else if(var < NUM_VARS + 4){
				varPrint = u_pri(i,j,var - NUM_VARS + 1);
			}else{
				varPrint = sqrt(u_con(i,j,5) * u_con(i,j,5) + u_con(i,j,6) * u_con(i,j,6)  + u_con(i,j,7) * u_con(i,j,7) );
			}
			outFile << std::setw(printWidth) << std::setprecision(std::numeric_limits<double>::digits10 + 1) << varPrint << " ";
		}
		outFile << std::endl;
	}
	
	outFile << std::endl;
}

/**Prints as many variables as desired to a single VTK file.
*/
void printAllVariablesToVTK(Grid& u_con,Grid& u_pri, std::string testname, int nstep, bool /*isFinalTime*/){
	string formula_type;
	string div_type;
	string update_type;
	string dt_type;
	string eos;
	string positive_preserving;
	switch (recon) {
		case NoRecon:
			SpaceOrder = 1;
			break;
		case Minmod: case VanLeer:
			SpaceOrder = 2;
			break;
		case ENO3:
			SpaceOrder = 3;
			break;
		case WENO5:
			SpaceOrder = 4;
			break;
	}
	
	switch (formula) {
		case EnergyFormula:
			formula_type = "EnergyFormula";
			update_type = "";
			break;
		case PressureFormula:
			formula_type = "PressureFormula";
			switch(U_update){
				case ConservationofU:
					update_type = "conservation";
					break;
				case DefinitionofU:
					update_type = "definition";
					break;
			}
			break;
	}
	
	switch(magfield){
		case StaggerCT:
			div_type = "StaggerCT";
			break;
		case UnStaggerCT:
			div_type = "UnStaggerCT";
			break;
		case PotentialA:
			div_type = "Potential";
			break;
		case NoDivFree:
			div_type = "NoDivFree";
			break;
	}
	
	switch(caltimestep){
		case Fullydt:
			dt_type = "Fullydt";
			break;
		case Alfvendt:
			dt_type = "Alfvendt";
			break;
		case IMEXdt:
			dt_type = "IMEXdt" + std::to_string(full_step);
			break;
	}
	
	switch(eos_type){
		case IDEAL_GAS:
			eos = "Ideal Gas";
			break;
		case REDLICH_KWONG:
			eos = "Redlich Kwong";
			break;
	}
	
	if(floor_on){
		positive_preserving = "PreservingON";
	}else{
		positive_preserving = "PreservingOFF";
	}
	if(nxCells == 1 || nyCells == 1){
		div_type = "";
	}
	
	
	
	std::ostringstream ss;
	ss << std::fixed << std::setprecision(0) << 100*CFL;
	std::string s = ss.str();  
	
	Path = output_root + "/" + eos  +"/" +positive_preserving +"/"+ testname + "/" + formula_type + "/"+ div_type + "/"+ update_type + "/" + std::to_string(nxCells) + "x"+ std::to_string(nyCells) + "_"+ dt_type+"_RK" + std::to_string(RKs) + "_Order" + std::to_string(SpaceOrder) +"_CFL0" + s;
	
	filesystem::create_directories(Path);
	
	std::string fileString = Path + "/" + "plt_" + testname; // "./Plot/" + name;
	
	fileString.append("_");
	fileString.append(to_string(nstep).substr(0,6));
	//if(isFinalTime) fileString.append("_FINAL");
	fileString.append(".vtk");
	
	std::ofstream outFile(fileString.c_str());
		
	
	printVTKHeader(outFile,nstep);
	
	std::vector< std::tuple<int,std::string,int,int,int> > variablesToPrint;
	
	variablesToPrint.push_back(std::tuple<int,std::string,int,int,int>(0,std::string("rho"),0,0,0));
	variablesToPrint.push_back(std::tuple<int,std::string,int,int,int>(1,std::string("mom_X"),0,0,0));
	variablesToPrint.push_back(std::tuple<int,std::string,int,int,int>(2,std::string("mom_Y"),0,0,0));
	variablesToPrint.push_back(std::tuple<int,std::string,int,int,int>(3,std::string("mom_Z"),0,0,0));
	variablesToPrint.push_back(std::tuple<int,std::string,int,int,int>(4,"U",0,0,0));
	variablesToPrint.push_back(std::tuple<int,std::string,int,int,int>(5,std::string("B_X"),0,0,0));
	variablesToPrint.push_back(std::tuple<int,std::string,int,int,int>(6,std::string("B_Y"),0,0,0));
	variablesToPrint.push_back(std::tuple<int,std::string,int,int,int>(7,std::string("B_Z"),0,0,0));
	variablesToPrint.push_back(std::tuple<int,std::string,int,int,int>(8,"velX",0,0,0));
	variablesToPrint.push_back(std::tuple<int,std::string,int,int,int>(9,"velY",0,0,0));
	variablesToPrint.push_back(std::tuple<int,std::string,int,int,int>(10,"velZ",0,0,0));
	variablesToPrint.push_back(std::tuple<int,std::string,int,int,int>(11,"p",0,0,0));
	variablesToPrint.push_back(std::tuple<int,std::string,int,int,int>(12,"B",0,0,0));
	
	
	outFile << "FIELD FieldData " << variablesToPrint.size() << std::endl;
	
	for(unsigned int i = 0;i <variablesToPrint.size(); i++){
		printVTKVariable(outFile,u_con,u_pri,std::get<0>(variablesToPrint[i]),std::get<1>(variablesToPrint[i]));
	}
	
	
	outFile << "POINT_DATA " << (nxCells+1)*(nyCells+1) << std::endl;
	
	return;
}

// ============================================================================
//  EMF-consistent Staggered CT  (Riccardo JCP eq.43 的通量对偶思路)
//  适用: PressureFormula 分支
//
//  原理: 面 EMF 不再用 f_full_con + Rusanov 独立重构, 而是直接取磁子系统
//  胞心更新实际使用的面通量, 变号得到 (E_z = -F_[By] / +G_[Bx]):
//      通量 = F^n 部分 (B* 装配里的 Qe 项, 你代码 5003/5005 行)
//           + h * G^{n+1} 部分 (Magnetic_Subsystem_AX 作用在解 B̃^{n+1} 上, 3753-3766 行)
//           + 隐式 LLF 扩散通量中可写成 EMF 的两个分量 (3811-3833 行)
//  这样 CT 继承隐式求解的 Alfvén 尺度阻尼, 而不是环外另造显式 EMF。
//
//  已知限制(设计使然, 不是 bug):
//  - Laplacian 中 Bx 的 x 向扩散 / By 的 y 向扩散无法写成 EMF 形式
//    (正是胞心算子破坏 div-free 的部分), CT 丢弃之。因此自检残差不为
//    机器零, 其量级 = 这两个分量 ~ O(ν·ΔB), 打印监控即可。
//  - 边界直接用 ghost 值 (bound_cond_u 已填), 周期/透射下与 AX 的
//    reflct 逻辑严格等价; 反射边界的符号细节未逐一对齐, vortex 不受影响。
//
//  安装步骤 (共 4 处):
//  [1] 本文件 SECTION A 的全局声明 → 放到你全局参数区 (GLF 声明附近)
//  [2] SECTION B 的两个函数 → 放到 Staggered_CT 定义之后
//      (依赖你已有的 CT_CornerAssemble —— 诊断代码里那个, 保留它)
//  [3] main 里 Array 声明区之后加:
//        gFluxBy_p0 = new Array(nxCells, nyCells);
//        gFluxBx_0p = new Array(nxCells, nyCells);
//  [4] 在包含 GMRES_Solver(B_vect, b_Bt_vect, Magnetic_Subsystem_AX, ...)
//      的那个函数末尾 —— B 写回 u_con 的循环和 bound_cond_u(u_con) 之后 ——
//      加一行:
//        CollectMagFluxes(u_con, Qe, Qs, h);
//      (h 用该函数收到的同一个 h; 千万不能放在别处, 此刻 u_con 的动量
//       还是 RHS 装配用的 (ρu)*, B 槽已是 B̃^{n+1}, 两者都必须是这个状态)
//  [5] main 的 switch(magfield) 里:
//        case StaggerCT:
//            Staggered_CT_Consistent(u_con, Bx_p0, By_0p, Bx_p0_Qs, By_0p_Qs, adt);
//            break;
//      (原 Staggered_CT 保留不删, 便于 A/B 对比)
//  文件头声明区补:
//        void CollectMagFluxes(Grid &u_con, Grid &Qe, Grid &Qs, double h);
//        void Staggered_CT_Consistent(Grid& u_con, Array &Bx_p0, Array &By_0p,
//                                     Array &Bx_p0_Qs, Array &By_0p_Qs, double h);
// ============================================================================



// ---------------------------- SECTION B: 函数 ------------------------------

// 在磁子系统解出 B̃^{n+1} 并写回 u_con 之后立刻调用 (安装步骤[4])。
// 此刻: u_con(·,·,1..2) = (ρu)^* (与 RHS 装配一致), u_con(·,·,5..7) = B̃^{n+1}。
