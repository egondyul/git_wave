function GenerateData(Session, FRACTYPE,Frequency,LayerWidth,maxVp,Dx,Nx,Time)
%FRACTYPE:
%1-Homogeneous layers, 2-parallel layers, 3-perpendiculat layers
folder=strcat(num2str(Session),'\',num2str(FRACTYPE),'\',num2str(Frequency));
eval(['!mkdir ' folder]);
folder=strcat(folder,'\');
folderData=strcat(folder,'Data');
eval(['!mkdir ' folderData]);
folderData=strcat(folder, 'Data\');
WaveLength = ceil(ceil(maxVp*100/Frequency)/100/Dx)     % Length of the signal
%% LayerSize for usual cases

LayerSize = 0*WaveLength;                            % width of heterogeneous layer

PML = 0;                                    % width of PML layer
% 
% Nz=4*WaveLength
% Size=Nx*Nz;
% Dz=Dx;
% Lz=Nz*Dz;
Nz = PML+6*WaveLength+LayerSize+PML
Size = Nx*Nz;
Dz = Dx;
LX = 3*Dx
LZ = (PML+6*WaveLength+LayerSize+PML)*Dz

Dx
%=============== Нахождение максимума фазовой скорости ====================
Dt = Dx*Dz/1300/(Dx+Dz)   % Gondyul: нужно сделать скриптик, который будет находить подходящий Dt
% Dt=5*10^(-9);
nu=100;%[kg/m s]=[Pa s]
rho = 2500;
% wmax=((32*16*pi^(4)*nu^2)/(rho^(2)*Dx^4))*(1+(1+(WaveLength^2*rho^(2)*Dx^(4))/(8^(2)*2^(4)*pi^(4)*nu^(4)))^(1/2))
% vmax=(1/((rho)^(1/2)))*(WaveLength^(2)+4*wmax^(2)*nu^(2))^(1/4)
% Dt=(1/vmax)%*((1/(Dx*Dx))+(1/(Dx*Dx)))^(1/2)
% 
% kMax=2/Dx
% lambda=1500^2*1000
% wMax2=32*pi^4*kMax^4*nu^2*(1+sqrt(1+lambda^2*rho^2/(8^2*pi^4*kMax^4*nu^4)))/rho^2
% VMax=(lambda^2+4*wMax2*nu^2)^(1/4)/sqrt(rho)
% Dt=sqrt(2/Dx^2)/VMax
%==========================================================================

%% ГЕНЕРАЦИЯ ТРЕЩИН - С УМЕНЬШЕНИЕМ КОНЦЕНРТАЦИИ В ПРАВОЙ ЧАСТИ УЧАСТКА
FracConc = 0.0625;                % концентрация трещин - отношение числа узлов трещин к общему кол-ву узлов области
Zstart = PML+WaveLength*4;   % номер узла по z, с которого начинаются трещины
Zend = PML+WaveLength*4+LayerSize;     % номер узла по z, на котором заканчиваются трещины
% Zstart=round(Nz/2);
% Zend=round(3*Nz/4);
%% DISPLAY GEOMETRIC PROPERTIES
PML
 Source = PML + WaveLength %important src_data
%Source=Zstart+2*WaveLength; %just for curious
Receiver1 = PML + 2*WaveLength
Receiver2 = PML + 5*WaveLength + LayerSize 
Zstart
Zend
DomainSize = Nz
%Time = LZ / 1000;           % Gondyul: вместо 1000 нужно что-то типа минимальной фазовой скорости в среде
LayerSize

% Source=round(Nz/4)-round(Nz/8);
% Receiver1=round(Nz/4)+round(Nz/8);
% Receiver2=Zend+round(Nz/8);
% Zstart
% Zend
% DomainSize = Nz
%Time = LZ / 1000;           % Gondyul: вместо 1000 нужно что-то типа минимальной фазовой скорости в среде
LayerSize

FINPUT = fopen([folder 'INPUT.txt'],'w');
FINPUTData=fopen([folderData 'INPUT.txt'],'w');
formatSpec = '%d               # v0                frequency\r\n';
fprintf(FINPUT, formatSpec, Frequency);
fprintf(FINPUTData, formatSpec, Frequency);
formatSpec = '%f              # Time              time\r\n';
fprintf(FINPUT, formatSpec, Time);
fprintf(FINPUTData, formatSpec, Time);
formatSpec = '%d                   # srcplace          source placement: 0-fluid, 1-solid\r\n';
fprintf(FINPUT, formatSpec, 1);
fprintf(FINPUTData, formatSpec, 1);
formatSpec = '%d                 # z_src             source z node\r\n';
fprintf(FINPUT, formatSpec, Source);
fprintf(FINPUTData, formatSpec, Source);
formatSpec = '%d                 # z_rec_1           receiver 1 z node\r\n';
fprintf(FINPUT, formatSpec, Receiver1);
fprintf(FINPUTData, formatSpec, Receiver1);
formatSpec = '%d                # z_rec_2           receiver 2 z node\r\n';
fprintf(FINPUT, formatSpec, Receiver2);
fprintf(FINPUTData, formatSpec, Receiver2);
formatSpec = '%d                 # Nz_PML_l          left PML length\r\n';
fprintf(FINPUT, formatSpec, PML);
formatSpec = '%d                 # Nz_PML_r          right PML length\r\n';
fprintf(FINPUT, formatSpec, PML);
formatSpec = '%d                   # p_rec\r\n';
fprintf(FINPUT, formatSpec, 1);
formatSpec = '%d                   # sigma_xx_rec\r\n';
fprintf(FINPUT, formatSpec, 1);
formatSpec = '%d                   # sigma_zz_rec\r\n';
fprintf(FINPUT, formatSpec, 1);
formatSpec = '%d                   # sigma_xz_rec\r\n';
fprintf(FINPUT, formatSpec, 1);
formatSpec = '%d                   # q_x_rec\r\n';
fprintf(FINPUT, formatSpec, 1);
formatSpec = '%d                   # q_z_rec\r\n';
fprintf(FINPUT, formatSpec, 1);
formatSpec = '%d                   # v_x_rec\r\n';
fprintf(FINPUT, formatSpec, 1);
formatSpec = '%d                   # v_z_rec\r\n';
fprintf(FINPUT, formatSpec, 1);
formatSpec = '%d                   # p_snap\r\n';
fprintf(FINPUT, formatSpec, 1);
formatSpec = '%d                   # sigma_xx_snap\r\n';
fprintf(FINPUT, formatSpec, 1);
formatSpec = '%d                   # sigma_zz_snap\r\n';
fprintf(FINPUT, formatSpec, 1);
formatSpec = '%d                   # sigma_xz_snap\r\n';
fprintf(FINPUT, formatSpec, 1);
formatSpec = '%d                   # q_x_snap\r\n';
fprintf(FINPUT, formatSpec, 1);
formatSpec = '%d                   # q_z_snap\r\n';
fprintf(FINPUT, formatSpec, 1);
formatSpec = '%d                   # v_x_snap\r\n';
fprintf(FINPUT, formatSpec, 1);
formatSpec = '%d                   # v_z_snap\r\n';
fprintf(FINPUT, formatSpec, 1);
formatSpec = '%d                   # NSnaps\r\n';
fprintf(FINPUT, formatSpec, 20);
formatSpec = '%d                   # MPIFlag            1-MPI 0-NO MPI\r\n';
fprintf(FINPUT, formatSpec, 0);
formatSpec = '%d                   # dt / tau [ T2 = ceil( T / # ) ]\r\n';
fprintf(FINPUT, formatSpec, 5);
fclose(FINPUT);


%for folderData with wavelength and Vp
FINPUTData=fopen([folderData 'INPUT.txt'],'w');
formatSpec = '%d               # v0                frequency\r\n';
fprintf(FINPUTData, formatSpec, Frequency);
formatSpec = '%f              # Time              time\r\n';
fprintf(FINPUTData, formatSpec, Time);
formatSpec = '%d                   # srcplace          source placement: 0-fluid, 1-solid\r\n';
fprintf(FINPUTData, formatSpec, 1);
formatSpec = '%d                 # z_src             source z node\r\n';
fprintf(FINPUTData, formatSpec, Source);
formatSpec = '%d                 # z_rec_1           receiver 1 z node\r\n';
fprintf(FINPUTData, formatSpec, Receiver1);
formatSpec = '%d                # z_rec_2           receiver 2 z node\r\n';
fprintf(FINPUTData, formatSpec, Receiver2);
formatSpec = '%d                # WaveLength           number of points per wavelength\r\n';
fprintf(FINPUTData, formatSpec, WaveLength);
fclose(FINPUTData);

NFracs = 0;       % число трещин
FracPointNumber = 0;              % честное число точек трещин
AAA(1:Nx*Nz) = 0;                 % А - индикатор трещин: если в узле А = 1, то там материал трещины, 0 - матрицы.


%% ЗАПОЛНЕНИЕ ЧАСТИ ОБЛАСТИ МАТЕРИАЛОМ ТРЕЩИН
if (FRACTYPE == 1)
    AAA((Zstart-1)*Nx+1:(Zend)*Nx)=1;
end
%% ЗАПОЛНЕНИЕ СЛОИСТЫМ МАТЕРИАЛОМ
% ПЕРПЕНДИКУЛЯРНЫЕ СЛОИ
if (FRACTYPE == 3)
    j = Zstart - 1;
    while (j <= Zend - LayerWidth)
       AAA(j*Nx + 1:(j+LayerWidth)*Nx) = 1; 
       j = j + LayerWidth*2;
    end
end

% ПАРАЛЛЕЛЬНЫЕ СЛОИ
if (FRACTYPE == 2)
    i =LayerWidth+1;
    while (i < Nx - LayerWidth)
        for l = 0:LayerWidth-1
           for j = Zstart - 1: Zend - 1
               AAA(j*Nx + i + l) = 1;
           end
        end
        i = i + 2*LayerWidth;
    end
    for j = Zstart: Zend
       AAA(j*Nx) = AAA(j);
    end
end

if(FRACTYPE==9)
    length=floor(Nz/50);
    delta=3*length;
    nl=floor(Nz/(length+delta))
    tmp=1; %delta+length
    j=0;
    while(j<=Nz-1)
         AAA(j*Nx + 1:(j+length)*Nx) = 1; 
       j = j + length*2;
    end
    for j=1:Nz-1
        for i=1:Nx-1
            if(i>Nx/3&&i<2*Nx/3)
                AAA(j*Nx+i)=1;
            end
                
        end
    end
end


%% ЗАПИСЬ ПАРАМЕТРОВ МОДЕЛИ

Fgrid = fopen([folder 'grid.bin'],'w');
fwrite(Fgrid,Nx,'int');
fwrite(Fgrid,Nz,'int');
fwrite(Fgrid,Dx,'double');
fwrite(Fgrid,Dz,'double');
fwrite(Fgrid,Dt,'double');
fclose(Fgrid);

FgridData = fopen([folderData 'grid.bin'],'w');
fwrite(FgridData,Nx,'int');
fwrite(FgridData,Nz,'int');
fwrite(FgridData,Dx,'double');
fwrite(FgridData,Dz,'double');
fwrite(FgridData,Dt,'double');
fclose(FgridData);

%% create pphysical properties files
MatPropExport_ViscoElast(AAA, Frequency, Nx, Nz, folder);

