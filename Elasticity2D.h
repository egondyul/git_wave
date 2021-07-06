#pragma once
#include <cmath>
#include <iostream> 
#include <fstream> 
#include <stdio.h>
#include <assert.h>
#include <omp.h>
#include <string.h>
#include <sstream>
#include <time.h>
#include<vector>

using namespace std;

typedef float realval;
typedef long myint; //for WINDOWS

class Elasticity2D
{
	struct inputData
	{
		realval v0; //frequency
		realval Time; // time of propagation
		int srcplace; //source placement: 0-fluid; 1-solid
		int z_src; // source z-node
		int z_rec_1; // receiver 1 z-node
		int z_rec_2; // receiver 2 z-node
		int Nz_PML_l; // left PML length
		int Nz_PML_r; // right PML length


		int sigma_xx_rec;  // normal stress in x-direction record lines data: 1-record; 0-do not record
		int sigma_zz_rec;  // normal stress in z-direction record lines data: 1-record; 0-do not record
		int sigma_xz_rec;  // shear stress record lines data: 1-record; 0-do not record
		int v_x_rec;       // x-component of solid velocity record lines data: 1-record; 0-do not record
		int v_z_rec;       // z-component of solid velocity record lines data: 1-record; 0-do not record
		int sigma_xx_snap; // normal stress in x-direction snapshots: 1-make; 0-do not make
		int sigma_zz_snap; // normal stress in z-direction snapshots: 1-make; 0-do not make
		int sigma_xz_snap; // shear stress snapshots: 1-make; 0-do not make
		int v_x_snap;      // x-component of solid velocity snapshots: 1-make; 0-do not make
		int v_z_snap;      // z-component of solid velocity snapshots: 1-make; 0-do not make
		int NSnaps;        // number of snapshots during propagation time (through equal time intervals)

	};

	void filename(char* name, int i, char* F);
	realval average_x(realval A[], int i, int j, int Nx);
	realval average_z(realval A[], int i, int j, int Nx);
	realval average_xz(realval A[], int i, int j, int Nx);

	realval average_x_vec(std::vector<realval> &A, int i, int j, int Nx);
	realval average_z_vec(std::vector<realval> &A, int i, int j, int Nx);
	realval average_xz_vec(std::vector<realval> &A, int i, int j, int Nx);
	void readData();
	//for outout files
	void outputData(); //sigma, v and time_scale
	void coordinate();

public:
	Elasticity2D();
	~Elasticity2D();
	void Elasticity();
private:
	inputData inpData;
	int Nx;
	int Nz;
	realval hx;
	realval hz;
	realval tau;
	char path1[150]="Data/";
	char FullPath[150];

private:
	//from input files
	realval *rho; //centres
	realval *Vp;
	realval *Vs;
	//unknown yet
	realval *lam;
	realval *mu;
	realval *sigma_xx;					// ���������� ���������� �������� ��������� �� x
	realval *sigma_zz ;					// ���������� ���������� �������� ��������� �� z
	realval *sigma_xz;		// ����������� ���������� �������� ���������
	realval *v_x;					// �������� ������ ������ �� x
	realval *v_z;

	//for viscoelasticity 
	realval** R_xx_old;//memory variables
	realval** R_zz_old;
	realval** R_xz_old;

	realval** R_xx_new;//memory variables
	realval** R_zz_new;
	realval** R_xz_new;
	//from unput files
	realval tau_sigma;
	realval* tau_11;
	realval* tau_13;
	realval* tau_33;
	realval* tau_55;


	//coordinate of source
	realval *src;

	//std::vector<realval> src;

	int i, j;
	int counter=1;
	realval cur_time, cur_time_str;			// ����� �� ������� ���� ��� ��������� � ���������� ��������������
	realval rx, rz;							// tau/hx, tau/hz
	realval DtPlot;                       // ��������������� ��������� � ������� ��������� � ��� �������
	const realval PI = 3.14159265f;
	int flag = 0;

};

