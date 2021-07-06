#include "Elasticity2D.h"

//------Without PML---------//

Elasticity2D::Elasticity2D()
{
	//initializing all parameters 
	readData(); //Vp,Vs,rho, hx, hz, tau, ...

	//memory alloc
	lam = new realval[Nx*Nz];						// постоянная Ламе
	mu = new realval[Nx*Nz];						// модуль сдвига

	//искомые величины на текущем шаге по времени
	sigma_xx = new realval[Nx*Nz];					// нормальное напряжение упругого материала по x
	sigma_zz = new realval[Nx*Nz];					// нормальное напряжение упругого материала по z
	sigma_xz = new realval[(Nx - 1)*(Nz - 1)];		// касательное напряжение упругого материала
	v_x = new realval[(Nx - 1)*Nz];					// скорость твердх частиц по x
	v_z = new realval[Nx*(Nz - 1)];					// скорость твердых частиц по z

	//for viscoelasticity
	R_xx_old = new realval*[Nx*Nz];
	R_zz_old = new realval*[Nx*Nz];
	R_xz_old = new realval*[(Nx-1)*(Nz-1)];

	R_xx_new = new realval*[Nx*Nz];
	R_zz_new = new realval*[Nx*Nz];
	R_xz_new = new realval*[(Nx - 1)*(Nz - 1)];

	tau_11 = new realval[Nx*Nz];
	tau_13 = new realval[Nx*Nz];
	tau_33 = new realval[Nx*Nz];
	tau_55 = new realval[(Nx - 1)*(Nz - 1)];

	src = new realval[Nx*Nz];

	//vector issues
	/*lam.resize(Nx*Nz);
	mu.resize(Nx*Nz);
	sigma_xx.resize(Nx*Nz);
	sigma_zz.resize(Nx*Nz);
	sigma_xz.resize((Nx - 1)*(Nz - 1));
	v_x.resize((Nx - 1)*Nz);
	v_z.resize(Nx*(Nz - 1));

	src.resize(Nx*Nz);*/


	rx = tau / hx;
	rz = tau / hz;

	//for rendering (output)
	coordinate();

	// Initial data (zero)

	int i, j;
#pragma omp parallel num_threads(4)
	{
#pragma omp for private(i,j) schedule(guided)
		for (j = 0; j < Nz; j++)
		{
			for (i = 0; i < Nx; i++)
			{
				sigma_xx[j*Nx + i] = 0;
				sigma_zz[j*Nx + i] = 0;
				src[j*Nx + i] = 0;
				if (i < Nx - 1)
				{
					v_x[j*(Nx - 1) + i] = 0;
				}
				if (j < Nz - 1)
				{
					v_z[j*Nx + i] = 0;
				}
				if ((j < Nz - 1) && (i < Nx - 1))
				{
					sigma_xz[j*(Nx - 1) + i] = 0;
				}
			}
		}
	}

	// set source points
	for (i = 0; i < Nx; i++)
	{
		src[inpData.z_src*Nx + i] = 1;
	}



	//mu, lambda
#pragma omp parallel
	{
#pragma omp for private(i,j) schedule(guided)
		for (j = 0; j < Nz; j++)
		{
			for (i = 0; i < Nx; i++)
			{
				mu[j*Nx + i] = rho[j*Nx + i] * Vs[j*Nx + i] * Vs[j*Nx + i];
				lam[j*Nx + i] = rho[j*Nx + i] * Vp[j*Nx + i] * Vp[j*Nx + i] - 2 * mu[j*Nx + i];
			}
		}
	}

	cur_time = 0;
	cur_time_str = 0;

}

Elasticity2D::~Elasticity2D()
{
	delete[]lam;
	delete[]mu;
	delete[]sigma_xx;
	delete[]sigma_xz;
	delete[]sigma_zz;
	delete[]v_x;
	delete[]v_z;
	delete[]src;
	delete[]rho;
	delete[]Vp;
	delete[]Vs;

	delete[]R_xx_old;//memory variables
	delete[]R_zz_old;
	delete[]R_xz_old;

	delete[]R_xx_new;//memory variables
	delete[]R_zz_new;
	delete[] R_xz_new;
	
	delete[] tau_11;
	delete[] tau_13;
	delete[] tau_33;
	delete[] tau_55;
}

void Elasticity2D::readData()
{
	ifstream parameters;
	parameters.open("INPUT.txt", ios::in);

	if (parameters.is_open())
	{
		double tmp; //smt i don't need
		parameters >> inpData.v0;
		cout << "v0 = " << inpData.v0 << "\n";
		parameters.ignore(1000, '\n');
		parameters >> inpData.Time;
		cout << "Time = " << inpData.Time << "\n";
		parameters.ignore(1000, '\n');
		parameters >> inpData.srcplace;
		cout << "srcplace = " << inpData.srcplace << "\n";
		parameters.ignore(1000, '\n');
		parameters >> inpData.z_src;
		cout << "z_src = " << inpData.z_src << "\n";
		parameters.ignore(1000, '\n');
		parameters >> inpData.z_rec_1;
		cout << "z_rec_1 = " << inpData.z_rec_1 << "\n";
		parameters.ignore(1000, '\n');
		parameters >> inpData.z_rec_2;
		cout << "z_rec_2 = " << inpData.z_rec_2 << "\n";
		parameters.ignore(1000, '\n');
		parameters >> inpData.Nz_PML_l;
		cout << "Nz_PML_l = " << inpData.Nz_PML_l << "\n";
		parameters.ignore(1000, '\n');
		parameters >> inpData.Nz_PML_r;
		cout << "Nz_PML_r = " << inpData.Nz_PML_r << "\n";
		parameters.ignore(1000, '\n');
		parameters >> tmp;
		parameters.ignore(1000, '\n');
		parameters >> inpData.sigma_xx_rec;
		cout << "sigma_xx_rec = " << inpData.sigma_xx_rec << "\n";
		parameters.ignore(1000, '\n');
		parameters >> inpData.sigma_zz_rec;
		cout << "sigma_zz_rec = " << inpData.sigma_zz_rec << "\n";
		parameters.ignore(1000, '\n');
		parameters >> inpData.sigma_xz_rec;
		cout << "sigma_xz_rec = " << inpData.sigma_xz_rec << "\n";
		parameters.ignore(1000, '\n');
		parameters >> tmp;
		parameters.ignore(1000, '\n');
		parameters >> tmp;
		parameters.ignore(1000, '\n');
		parameters >> inpData.v_x_rec;
		cout << "v_x_rec = " << inpData.v_x_rec << "\n";
		parameters.ignore(1000, '\n');
		parameters >> inpData.v_z_rec;
		cout << "v_z_rec = " << inpData.v_z_rec << "\n";
		parameters.ignore(1000, '\n');
		parameters >> tmp;
		parameters.ignore(1000, '\n');
		parameters >> inpData.sigma_xx_snap;
		cout << "sigma_xx_snap = " << inpData.sigma_xx_snap << "\n";
		parameters.ignore(1000, '\n');
		parameters >> inpData.sigma_zz_snap;
		cout << "sigma_zz_snap = " << inpData.sigma_zz_snap << "\n";
		parameters.ignore(1000, '\n');
		parameters >> inpData.sigma_xz_snap;
		cout << "sigma_xz_snap = " << inpData.sigma_xz_snap << "\n";
		parameters.ignore(1000, '\n');
		parameters >> tmp;
		parameters.ignore(1000, '\n');
		parameters >> tmp;
		parameters.ignore(1000, '\n');
		parameters >> inpData.v_x_snap;
		cout << "v_x_snap = " << inpData.v_x_snap << "\n";
		parameters.ignore(1000, '\n');
		parameters >> inpData.v_z_snap;
		cout << "v_z_snap = " << inpData.v_z_snap << "\n";
		parameters.ignore(1000, '\n');
		parameters >> inpData.NSnaps;
		cout << "NSnaps = " << inpData.NSnaps << "\n";
		parameters.ignore(1000, '\n');
	}
	else
	{
		std::cout << "INPUT.txt failed." << std::endl;
	}
	parameters.close();

	inpData.Time = inpData.Time + 3 / inpData.v0; // add time for signal to form

	ifstream grid;
	grid.open("grid.bin", ios::binary | ios::in);
	if (grid.is_open())
	{
		grid.read((char*)&Nx, sizeof(int));
		grid.read((char*)&Nz, sizeof(int));
		grid.read((char*)&hx, sizeof(realval));
		grid.read((char*)&hz, sizeof(realval));
		grid.read((char*)&tau, sizeof(realval));
	}
	else
	{
		std::cout << "grid.bin failed." << std::endl;
	}
	grid.close();

	cout << "Nx="<< Nx << endl;
	cout << "Nz=" << Nz << endl;

	fflush(stdout);


	rho = new realval[Nx*Nz];
	Vp = new realval[Nx*Nz];
	Vs = new realval[Nx*Nz];

	/*rho.resize(Nx*Nz);
	Vp.resize(Nx*Nz);
	Vs.resize(Nx*Nz);*/

	ifstream Frho;
	Frho.open("rho.bin", ios::binary | ios::in);
	if (!Frho.is_open())
	{
		std::cout << "rho.bin failed." << std::endl;
	}

	ifstream FVp;
	FVp.open("Vp.bin", ios::binary | ios::in);
	if (!FVp.is_open())
	{
		std::cout << "Vp.bin failed." << std::endl;
	}

	ifstream FVs;
	FVs.open("Vs.bin", ios::binary | ios::in);
	if (!FVs.is_open())
	{
		std::cout << "Vs.bin failed." << std::endl;
	}

	for (i = 0; i < Nx*Nz; i++)
	{
		Frho.read((char*)&rho[i], sizeof hx); 
		FVp.read((char*)&Vp[i], sizeof hx);
		FVs.read((char*)&Vs[i], sizeof hx);    
	}
	Frho.close();
	FVp.close();
	FVs.close();

	//read the tau time of relacsation and tensor tau



}

void Elasticity2D:: filename(char* name, int i, char* F)
{
	char num[150];
	sprintf_s(num, "%d", i);
	int size_F = sizeof(F);
	int size_num = sizeof(i);
	strcpy_s(F, strlen(F) + strlen(name)+1, name);
	strcat_s(F, strlen(F) + strlen(num) + 1, num);
	const char* s = ".bin";
	strcat_s(F, strlen(F) + strlen(s) + 1, s);
}

realval Elasticity2D::average_x(realval A[], int i, int j, int Nx)
{
	return 0.5*(A[j*Nx + i] + A[j*Nx + i + 1]);
}



realval Elasticity2D::average_z(realval A[], int i, int j, int Nx)
{
	return 0.5*(A[j*Nx + i] + A[(j + 1)*Nx + i]);
}

realval Elasticity2D::average_xz(realval A[], int i, int j, int Nx)
{
	realval av;
	if ((A[j*Nx + i] < 10) || (A[j*Nx + i + 1] < 10) || (A[(j + 1)*Nx + i] < 10) || (A[(j + 1)*Nx + i + 1] < 10))
		return 0;
	else
	{
		av = 0.25*(1 / A[j*Nx + i] + 1 / A[j*Nx + i + 1] + 1 / A[(j + 1)*Nx + i] + 1 / A[(j + 1)*Nx + i + 1]);
		return 1 / av;
	}
}

realval Elasticity2D::average_x_vec(std::vector<realval>&A, int i, int j, int Nx)
{
	return 0.5*(A[j*Nx + i] + A[j*Nx + i + 1]);
}



realval Elasticity2D::average_z_vec(std::vector<realval>&A, int i, int j, int Nx)
{
	return 0.5*(A[j*Nx + i] + A[(j + 1)*Nx + i]);
}

realval Elasticity2D::average_xz_vec(std::vector<realval>&A, int i, int j, int Nx)
{
	realval av;
	if ((A[j*Nx + i] < 10) || (A[j*Nx + i + 1] < 10) || (A[(j + 1)*Nx + i] < 10) || (A[(j + 1)*Nx + i + 1] < 10))
		return 0;
	else
	{
		av = 0.25*(1 / A[j*Nx + i] + 1 / A[j*Nx + i + 1] + 1 / A[(j + 1)*Nx + i] + 1 / A[(j + 1)*Nx + i + 1]);
		return 1 / av;
	}
}


void Elasticity2D::coordinate()
{
	realval *x_int = new realval[Nx];						// координаты "целых" узлов по x
	realval *z_int = new realval[Nz];						// координаты "целых" узлов по z
	realval *x_med = new realval[Nx - 1];					// координаты "дробных" узлов по x
	realval *z_med = new realval[Nz - 1];

	// Определение координат узлов
	for (i = 0; i < Nx; i++)
	{
		x_int[i] = i * hx;
	}

	for (i = 0; i < Nx - 1; i++)
	{
		x_med[i] = 0.5*hx + i * hx;
	}

	for (i = 0; i < Nz; i++)
	{
		z_int[i] = i * hz;
	}

	for (i = 0; i < Nz - 1; i++)
	{
		z_med[i] = 0.5*hz + i * hz;
	}

	ofstream x_int_res;
	x_int_res.open("x_int.bin", ios::binary | ios::out);

	ofstream z_int_res;
	z_int_res.open("z_int.bin", ios::binary | ios::out);

	ofstream x_med_res;
	x_med_res.open("x_med.bin", ios::binary | ios::out);

	ofstream z_med_res;
	z_med_res.open("z_med.bin", ios::binary | ios::out);

	for (i = 0; i < Nx; i++)
	{
		x_int_res.write((char*)&x_int[i], sizeof rx);
	}

	for (i = 0; i < Nz; i++)
	{
		z_int_res.write((char*)&z_int[i], sizeof rx);
	}

	for (i = 0; i < Nx - 1; i++)
	{
		x_med_res.write((char*)&x_med[i], sizeof rx);
	}

	for (i = 0; i < Nz - 1; i++)
	{
		z_med_res.write((char*)&z_med[i], sizeof rx);
	}

	x_int_res.close();
	z_int_res.close();
	x_med_res.close();
	z_med_res.close();

}

void Elasticity2D::outputData()
{
	// запись данных в приемниках в файлы

	// создание снэпов
	if (inpData.NSnaps > 0)
	{
		DtPlot = inpData.Time / (float)inpData.NSnaps;
	}

	if (cur_time==0 || ((fabs(cur_time_str - counter * DtPlot) <= tau) && (counter <= inpData.NSnaps)))
	{
		ofstream v_x_res;
		ofstream v_z_res;
		ofstream q_x_res;
		ofstream q_z_res;
		ofstream sigma_xx_res;
		ofstream sigma_zz_res;
		ofstream sigma_xz_res;
		ofstream p_res;

		if (inpData.v_x_snap == 1)
		{
		   strcpy_s(path1, sizeof(path1), "v_x_"); // не для MPI
	
			filename(path1, counter, FullPath);
			v_x_res.open(FullPath, ios::binary | ios::out);                  // не для MPI
			v_x_res.write((char*)v_x, (sizeof(rx)*(Nx - 1)*Nz));
			//v_x_res.write((char*)&v_x, (sizeof(rx)*(Nx - 1)*Nz));
			v_x_res.close();
		}

		if (inpData.v_z_snap == 1)
		{
			strcpy_s(path1, sizeof(path1), "v_z_"); // не для MPI
		
			filename(path1, counter, FullPath);
			v_z_res.open(FullPath, ios::binary | ios::out);                  // не для MPI
			v_z_res.write((char*)v_z, (sizeof(rx)*(Nz - 1)*Nx));
			//v_z_res.write((char*)&v_z, (sizeof(rx)*(Nz - 1)*Nx));
			v_z_res.close();
		}

		if (inpData.sigma_xx_snap == 1)
		{
			strcpy_s(path1, sizeof(path1), "sigma_xx_"); // не для MPI
			
			filename(path1, counter, FullPath);
			sigma_xx_res.open(FullPath, ios::binary | ios::out);                  // не для MPI
			sigma_xx_res.write((char*)sigma_xx, (sizeof(rx)*Nx*Nz));
			//sigma_xx_res.write((char*)&sigma_xx, (sizeof(rx)*Nx*Nz));
			sigma_xx_res.close();
		}

		if (inpData.sigma_zz_snap == 1)
		{
			strcpy_s(path1, sizeof(path1), "sigma_zz_"); // не для MPI
			
			filename(path1, counter, FullPath);
			sigma_zz_res.open(FullPath, ios::binary | ios::out);                  // не для MPI
			sigma_zz_res.write((char*)sigma_zz, (sizeof(rx)*Nx*Nz));
			//sigma_zz_res.write((char*)&sigma_zz, (sizeof(rx)*Nx*Nz));
			sigma_zz_res.close();
		}

		if (inpData.sigma_xz_snap == 1)
		{
			strcpy_s(path1, sizeof(path1), "sigma_xz_"); // не для MPI
			
			filename(path1, counter, FullPath);
			sigma_xz_res.open(FullPath, ios::binary | ios::out);                  // не для MPI
			sigma_xz_res.write((char*)sigma_xz, (sizeof(rx)*(Nx - 1)*(Nz - 1)));
			//sigma_xz_res.write((char*)&sigma_xz, (sizeof(rx)*(Nx - 1)*(Nz - 1)));
			sigma_xz_res.close();
		}
		counter = counter + 1;
	}
}


void Elasticity2D::Elasticity()
{
	//time_scale output
	ofstream time_scale_vel;
	strcpy_s(path1, sizeof(path1), "time_scale_vel.bin"); // не для MPI
	time_scale_vel.open(path1, ios::binary | ios::out);

	ofstream time_scale_str;
	strcpy_s(path1, sizeof(path1), "time_scale_str.bin"); // не для MPI
	time_scale_str.open(path1, ios::binary | ios::out);

	time_scale_vel.write((char*)&cur_time, sizeof rx);
	time_scale_str.write((char*)&cur_time_str, sizeof rx);

	//zero step data output
	outputData();

	//_rec ouput 

	//first receiver
	ofstream v_x_rec_1;
	strcpy_s(path1, sizeof(path1), "v_x_rec_1.bin"); // не для MPI
	v_x_rec_1.open(path1, ios::binary | ios::out);

	ofstream v_z_rec_1;
	strcpy_s(path1, sizeof(path1), "v_z_rec_1.bin"); // не для MPI
	v_z_rec_1.open(path1, ios::binary | ios::out);

	ofstream sigma_xx_rec_1;
	strcpy_s(path1, sizeof(path1), "sigma_xx_rec_1.bin"); // не для MPI
	sigma_xx_rec_1.open(path1, ios::binary | ios::out);

	ofstream sigma_zz_rec_1;
	strcpy_s(path1, sizeof(path1), "sigma_zz_rec_1.bin"); // не для MPI
	sigma_zz_rec_1.open(path1, ios::binary | ios::out);

	ofstream sigma_xz_rec_1;
	strcpy_s(path1, sizeof(path1), "sigma_xz_rec_1.bin"); // не для MPI
	sigma_xz_rec_1.open(path1, ios::binary | ios::out);

	//second receiver
	ofstream v_x_rec_2;
	strcpy_s(path1, sizeof(path1), "v_x_rec_2.bin"); // не для MPI
	v_x_rec_2.open(path1, ios::binary | ios::out);

	ofstream v_z_rec_2;
	strcpy_s(path1, sizeof(path1), "v_z_rec_2.bin"); // не для MPI
	v_z_rec_2.open(path1, ios::binary | ios::out);

	ofstream sigma_xx_rec_2;
	strcpy_s(path1, sizeof(path1), "sigma_xx_rec_2.bin"); // не для MPI
	sigma_xx_rec_2.open(path1, ios::binary | ios::out);

	ofstream sigma_zz_rec_2;
	strcpy_s(path1, sizeof(path1), "sigma_zz_rec_2.bin"); // не для MPI
	sigma_zz_rec_2.open(path1, ios::binary | ios::out);

	ofstream sigma_xz_rec_2;
	strcpy_s(path1, sizeof(path1), "sigma_xz_rec_2.bin"); // не для MPI
	sigma_xz_rec_2.open(path1, ios::binary | ios::out);

	if (inpData.v_x_rec == 1)
	{
		for (i = 0; i < Nx - 1; i++)
		{
			v_x_rec_1.write((char*)&v_x[inpData.z_rec_1 * (Nx - 1) + i], sizeof rx);
			v_x_rec_2.write((char*)&v_x[inpData.z_rec_2 * (Nx - 1) + i], sizeof rx);
		}
	}

	if (inpData.v_z_rec == 1)
	{
		for (i = 0; i < Nx; i++)
		{
			v_z_rec_1.write((char*)&v_z[inpData.z_rec_1 * Nx + i], sizeof rx);
			v_z_rec_2.write((char*)&v_z[inpData.z_rec_2 * Nx + i], sizeof rx);
		}
	}

	if (inpData.sigma_xx_rec == 1)
	{
		for (i = 0; i < Nx; i++)
		{
			sigma_xx_rec_1.write((char*)&sigma_xx[inpData.z_rec_1 * Nx + i], sizeof rx);
			sigma_xx_rec_2.write((char*)&sigma_xx[inpData.z_rec_2 * Nx + i], sizeof rx);
		}
	}

	if (inpData.sigma_zz_rec == 1)
	{
		for (i = 0; i < Nx; i++)
		{
			sigma_zz_rec_1.write((char*)&sigma_zz[inpData.z_rec_1 * Nx + i], sizeof rx);
			sigma_zz_rec_2.write((char*)&sigma_zz[inpData.z_rec_2 * Nx + i], sizeof rx);
		}
	}

	if (inpData.sigma_xz_rec == 1)
	{
		for (i = 0; i < Nx - 1; i++)
		{
			sigma_xz_rec_1.write((char*)&sigma_xz[inpData.z_rec_1 * (Nx - 1) + i], sizeof rx);
			sigma_xz_rec_2.write((char*)&sigma_xz[inpData.z_rec_2 * (Nx - 1) + i], sizeof rx);
		}
	}


	realval v0 = inpData.v0;
	int i, j;
	realval f_src, c1, c2, c3, c4;

	while (cur_time < inpData.Time)
	{
		cur_time_str = cur_time + 0.5*tau;
		cur_time = cur_time + tau;

		time_scale_vel.write((char*)&cur_time, sizeof rx);
		time_scale_str.write((char*)&cur_time_str, sizeof rx);

		//memory variables
		realval** swap;
		swap = R_xx_new;
		R_xx_new = R_xx_old;
		R_xx_old = swap;

#pragma omp parallel
		{
#pragma omp for private(i,j) shedule(guided)
			for (j = 1; j < Nz - 1; j++)
			{
				for (i = 1; i < Nx - 1; i++)
				{
					(*R_xx_new)[j*Nx + i] = (1 / (2 * tau_sigma + tau))*(2 * tau_sigma* (*R_xx_old)[j*Nx + i] - tau * (*R_xx_old)[j*Nx + i] - tau_11[j*Nx + i] * 2 * rx*(v_x[j*(Nx - 1) + i] - v_x[j*(Nx - 1) + i - 1]) + tau_13[j*Nx + i] * 2 * rz* (v_z[j*Nx + i] - v_z[(j - 1)*Nx + i]));
					(*R_zz_new)[j*Nx + i] = (1 / (2 * tau_sigma + tau))*(2 * tau_sigma* (*R_zz_old)[j*Nx + i] - tau * (*R_zz_old)[j*Nx + i] - tau_13[j*Nx + i] * 2 * rx*(v_x[j*(Nx - 1) + i] - v_x[j*(Nx - 1) + i - 1]) + tau_33[j*Nx + i] * 2 * rz* (v_z[j*Nx + i] - v_z[(j - 1)*Nx + i]));

				}
			}
		}
		//periodic condition
#pragma omp parallel
		{
#pragma omp for private(j) schedule(guided)
			for (j = 1; j < Nz - 1; j++)
			{
				(*R_xx_new)[j*Nx] = (1 / (2 * tau_sigma + tau))*(2 * tau_sigma* (*R_xx_old)[j*Nx] - tau * (*R_xx_old)[j*Nx] - tau_11[j*Nx] * 2 * rx*(v_x[j*(Nx - 1)] - v_x[j*(Nx - 1) + (Nx - 2)]) + tau_13[j*Nx] * 2 * rz* (v_z[j*Nx] - v_z[(j - 1)*Nx]));
				(*R_zz_new)[j*Nx] = (1 / (2 * tau_sigma + tau))*(2 * tau_sigma* (*R_zz_old)[j*Nx] - tau * (*R_zz_old)[j*Nx] - tau_13[j*Nx] * 2 * rx*(v_x[j*(Nx - 1)] - v_x[j*(Nx - 1) + (Nx - 2)]) + tau_33[j*Nx] * 2 * rz* (v_z[j*Nx] - v_z[(j - 1)*Nx]));
				(*R_xx_new)[j*Nx + (Nx - 1)] = (*R_xx_new)[j*Nx];
				(*R_zz_new)[j*Nx + (Nx - 1)] = (*R_zz_new)[j*Nx];
			}
		}

#pragma omp parallel
		{
#pragma omp for private(i,j) schedule(guided)
			for (j = 0; j < Nz - 1; j++)
			{
				for (i = 0; i < Nx - 1; i++)
				{
					(*R_xz_new)[j*(Nx - 1) + i] = (1 / (2 * tau_sigma + tau))*(2 * tau_sigma*(*R_xz_old)[j*(Nx - 1) + i] - tau * (*R_xz_old)[j*(Nx - 1) + i] - tau_55[j*(Nx - 1) + i] * 2 * (rx*(v_z[j*Nx + i + 1] - v_z[j*Nx + i]) + rz * (v_x[(j + 1)*(Nx - 1) + i] - v_x[j*(Nx - 1) + i])));
				}
			}
		}

		// Without PML
						// нахождение напряжений и давления жидкости на шаге времени n+1/2
						#pragma omp parallel
						{
							#pragma omp for private(i,j,f_src) schedule(guided)
							for (j = 1; j < Nz - 1; j++)
							{
								for (i = 1; i < Nx - 1; i++)
								{
									// формула источника
									f_src = (1 - 2 * (PI*inpData.v0*(cur_time_str - 3 / v0))*(PI*v0*(cur_time_str - 3 / v0)))*exp(-((PI*v0*(cur_time_str - 3 / v0))*(PI*v0*(cur_time_str - 3 / v0))));
									sigma_xx[j*Nx + i] = sigma_xx[j*Nx + i] + rx * (lam[j*Nx + i] + 2 * mu[j*Nx + i])*(v_x[j*(Nx - 1) + i] - v_x[j*(Nx - 1) + i - 1]) + rz * lam[j*Nx + i] * (v_z[j*Nx + i] - v_z[(j - 1)*Nx + i]) + src[j*Nx + i] * f_src;
									sigma_zz[j*Nx + i] = sigma_zz[j*Nx + i] + rx * lam[j*Nx + i] * (v_x[j*(Nx - 1) + i] - v_x[j*(Nx - 1) + i - 1]) + rz * (lam[j*Nx + i] + 2 * mu[j*Nx + i]) * (v_z[j*Nx + i] - v_z[(j - 1)*Nx + i]) + src[j*Nx + i] * f_src;
									//viscoelasiticity
									sigma_xx[j*Nx + i] = sigma_xx[j*Nx + i] - (tau_sigma / tau)*((*R_xx_new)[j*Nx + i] - (*R_xx_old)[j*Nx + i]);
									sigma_zz[j*Nx + i] = sigma_zz[j*Nx + i] - (tau_sigma / tau)*((*R_zz_new)[j*Nx + i] - (*R_zz_old)[j*Nx + i]);

								}
							}
						}
						//периодические граничные условия
						#pragma omp parallel
						{
							#pragma omp for private(j,f_src) schedule(guided)
							for (j = 1; j < Nz - 1; j++)
							{
								// формула источника
								f_src = (1 - 2 * (PI*v0*(cur_time_str - 3 / v0))*(PI*v0*(cur_time_str - 3 / v0)))*exp(-((PI*v0*(cur_time_str - 3 / v0))*(PI*v0*(cur_time_str - 3 / v0))));
								sigma_xx[j*Nx] = sigma_xx[j*Nx] + rx * (lam[j*Nx] + 2 * mu[j*Nx])*(v_x[j*(Nx - 1)] - v_x[j*(Nx - 1) + (Nx - 2)]) + rz * lam[j*Nx] * (v_z[j*Nx] - v_z[(j - 1)*Nx]) + src[j*Nx] * f_src;
								sigma_zz[j*Nx] = sigma_zz[j*Nx] + rx * lam[j*Nx] * (v_x[j*(Nx - 1)] - v_x[j*(Nx - 1) + (Nx - 2)]) + rz * (lam[j*Nx] + 2 * mu[j*Nx]) * (v_z[j*Nx] - v_z[(j - 1)*Nx]) + src[j*Nx] * f_src;
								//viscoelasticity
								sigma_xx[j*Nx] = sigma_xx[j*Nx] - (tau_sigma / tau)*((*R_xx_new)[j*Nx] - (*R_xx_old)[j*Nx]);
								sigma_zz[j*Nx] = sigma_zz[j*Nx] - (tau_sigma / tau)*((*R_zz_new)[j*Nx] - (*R_zz_old)[j*Nx]);

								sigma_xx[j*Nx + (Nx - 1)] = sigma_xx[j*Nx];
								sigma_zz[j*Nx + (Nx - 1)] = sigma_zz[j*Nx];
							}
						}
		#pragma omp parallel
						{
		#pragma omp for private(i,j) schedule(guided)
							for (j = 0; j < Nz - 1; j++)
							{
								for (i = 0; i < Nx - 1; i++)
								{
									sigma_xz[j*(Nx - 1) + i] = sigma_xz[j*(Nx - 1) + i] + average_xz(mu, i, j, Nx)*(rx*(v_z[j*Nx + i + 1] - v_z[j*Nx + i]) + rz*(v_x[(j + 1)*(Nx - 1) + i] - v_x[j*(Nx - 1) + i]));
									//viscoelasticity
									sigma_xz[j*(Nx - 1) + i] = sigma_xz[j*(Nx - 1) + i] - (tau_sigma / tau)*((*R_xz_new)[j*(Nx - 1) + i] - (*R_xz_old)[j*(Nx - 1) + i]);
								}
							}
						}
						// нахождение скорости жидкости и твердых частиц на n+1 (горизонтальная компонента)
		#pragma omp parallel
						{
		#pragma omp for private(i,j,c1,c2,c3,c4) schedule(guided)
							for (j = 1; j < Nz - 1; j++)
							{
								for (i = 0; i < Nx - 1; i++)
								{
									c3 = 1 / average_x(rho, i, j, Nx);
									v_x[j*(Nx - 1) + i] = v_x[j*(Nx - 1) + i] + c3 * (rx*(sigma_xx[j*Nx + i + 1] - sigma_xx[j*Nx + i]) + rz * (sigma_xz[j*(Nx - 1) + i] - sigma_xz[(j - 1)*(Nx - 1) + i]));
								}
							}
						}
		
						// нахождение скорости жидкости и твердых частиц на n+1 (вертикальная компонента)
		#pragma omp parallel
						{
		#pragma omp for private(i,j,c1,c2,c3,c4) schedule(guided)
							for (j = 0; j < Nz - 1; j++)
							{
								for (i = 1; i < Nx - 1; i++)
								{
									c3 = 1 / average_z(rho, i, j, Nx);
									v_z[j*Nx + i] = v_z[j*Nx + i] + c3 * (rx*(sigma_xz[j*(Nx - 1) + i] - sigma_xz[j*(Nx - 1) + i - 1]) + rz * (sigma_zz[(j + 1)*Nx + i] - sigma_zz[j*Nx + i]));
								}
							}
						}
						//Периодические граничные условия
		#pragma omp parallel
						{
		#pragma omp for private(j,c1,c2,c3,c4) schedule(guided)
							for (j = 0; j < Nz - 1; j++)
							{
								c3 = 1 / average_z(rho, 0, j, Nx);
								v_z[j*Nx] = v_z[j*Nx] + c3 * (rx*(sigma_xz[j*(Nx - 1)] - sigma_xz[j*(Nx - 1) + (Nx - 2)]) + rz * (sigma_zz[(j + 1)*Nx] - sigma_zz[j*Nx]));
								v_z[j*Nx + (Nx - 1)] = v_z[j*Nx];

							}
						}


						outputData();

						if (inpData.v_x_rec == 1)
						{
							for (i = 0; i < Nx - 1; i++)
							{
								v_x_rec_1.write((char*)&v_x[inpData.z_rec_1 * (Nx - 1) + i], sizeof rx);
								v_x_rec_2.write((char*)&v_x[inpData.z_rec_2 * (Nx - 1) + i], sizeof rx);
							}
						}

						if (inpData.v_z_rec == 1)
						{
							for (i = 0; i < Nx; i++)
							{
								v_z_rec_1.write((char*)&v_z[inpData.z_rec_1 * Nx + i], sizeof rx);
								v_z_rec_2.write((char*)&v_z[inpData.z_rec_2 * Nx + i], sizeof rx);
							}
						}

						if (inpData.sigma_xx_rec == 1)
						{
							for (i = 0; i < Nx; i++)
							{
								sigma_xx_rec_1.write((char*)&sigma_xx[inpData.z_rec_1 * Nx + i], sizeof rx);
								sigma_xx_rec_2.write((char*)&sigma_xx[inpData.z_rec_2 * Nx + i], sizeof rx);
							}
						}

						if (inpData.sigma_zz_rec == 1)
						{
							for (i = 0; i < Nx; i++)
							{
								sigma_zz_rec_1.write((char*)&sigma_zz[inpData.z_rec_1 * Nx + i], sizeof rx);
								sigma_zz_rec_2.write((char*)&sigma_zz[inpData.z_rec_2 * Nx + i], sizeof rx);
							}
						}

						if (inpData.sigma_xz_rec == 1)
						{
							for (i = 0; i < Nx - 1; i++)
							{
								sigma_xz_rec_1.write((char*)&sigma_xz[inpData.z_rec_1 * (Nx - 1) + i], sizeof rx);
								sigma_xz_rec_2.write((char*)&sigma_xz[inpData.z_rec_2 * (Nx - 1) + i], sizeof rx);
							}
						}

	}

	time_scale_vel.close();
	time_scale_str.close();

}