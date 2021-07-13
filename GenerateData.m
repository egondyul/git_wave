function GenerateData(FRACTYPE,Frequency,LayerWidth,folder,path,name,maxVp,Dx,Nx,Time)
% FRACTYPE = 7; % 1-PARALLEL; 2-PERPENDICULAR; 3-NONINTERSECTING; 4-INTERSECTING; 5-PARALLEL LAYERS; 6-PERPENDICULAR LAYERS; 7-HOMOGENEOUS LAYER; 8-IMPORT LAYER
% easy domain properties
folderData=strcat(folder,'Data');
eval(['!mkdir ' folderData]);
folderData=strcat(folder, 'Data\');
WaveLength = ceil(ceil(maxVp*100/Frequency)/100/Dx)     % Length of the signal
%% LayerSize for usual cases
if (FRACTYPE < 8)
    LayerSize = 4*WaveLength;                            % width of heterogeneous layer
%     LayerSize = ceil(0.5*LayerSize + rand(1)*0.5*LayerSize);  % addition for good estimation
end
if (FRACTYPE == 8)
    LayerSize = min(4*WaveLength,2000);
%     LayerSize = ceil(0.5*LayerSize + rand(1)*0.5*LayerSize);  % addition for good estimation
end
PML = WaveLength;                                    % width of PML layer

Nz = PML+6*WaveLength+LayerSize+PML
Size = Nx*Nz;
Dz = Dx;
LX = 3*Dx
LZ = (PML+6*WaveLength+LayerSize+PML)*Dz

Dx
%=============== Нахождение максимума фазовой скорости ====================
Dt = Dx*Dz/1300/(Dx+Dz);   % Gondyul: нужно сделать скриптик, который будет находить подходящий Dt
% nu=10;%[kg/m s]=[Pa s]
% rho = 2540;
% wmax=((32*16*pi^(4)*nu^2)/(rho^(2)*Dx^4))*(1+(1+(WaveLength^2*rho^(2)*Dx^(4))/(8^(2)*2^(4)*pi^(4)*nu^(4)))^(1/2))
% vmax=(1/((rho)^(1/2)))*(WaveLength^(2)+4*wmax^(2)*nu^(2))^(1/4)
% Dt=(1/vmax)%*((1/(Dx*Dx))+(1/(Dx*Dx)))^(1/2)
%==========================================================================

%% ГЕНЕРАЦИЯ ТРЕЩИН - С УМЕНЬШЕНИЕМ КОНЦЕНРТАЦИИ В ПРАВОЙ ЧАСТИ УЧАСТКА
FracConc = 0.0625;                % концентрация трещин - отношение числа узлов трещин к общему кол-ву узлов области
Zstart = PML+WaveLength*4;   % номер узла по z, с которого начинаются трещины
Zend = PML+WaveLength*4+LayerSize;     % номер узла по z, на котором заканчиваются трещины

%% DISPLAY GEOMETRIC PROPERTIES
PML
Source = PML + WaveLength
Receiver1 = PML + 2*WaveLength
Receiver2 = PML + 5*WaveLength + LayerSize 
Zstart
Zend
DomainSize = Nz
Time = LZ / 1000;           % Gondyul: вместо 1000 нужно что-то типа минимальной фазовой скорости в среде
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
fclose(FINPUTData);

NFracs = 0;       % число трещин
FracPointNumber = 0;              % честное число точек трещин
AAA(1:Nx*Nz) = 0;                 % А - индикатор трещин: если в узле А = 1, то там материал трещины, 0 - матрицы.
%% СЧИТЫВАНИЕ С МОДЕЛИ С ГРУППАМИ ТРЕЩИН
if (FRACTYPE == 8)
    %% for bin
    frec = fopen([path name],'r');
    REC = fread(frec,inf,'int');
    fclose(frec);

    for j=1:length(REC)
        if (REC(j)>0)
            REC(j) = 1;
        end
    end
    AAA((Zstart-1)*Nx+1:(Zstart+LayerSize-1)*Nx) = REC(1:Nx*LayerSize)';
end
%% СЧИТЫВАНИЕ МОДЕЛИ И ЕЕ РАСТЯГИВАНИЕ ОТРАЖЕНИЕМ ИСХОДНОЙ
if (FRACTYPE == 9)
    frec = fopen([path name],'r');
    REC = fread(frec,inf,'int');
    fclose(frec);
    REC = reshape(REC,[500 length(REC)/500]);
    REC = REC(:,101:1900);
    REC(:,1801:3600) = REC(:,1800:-1:1);
    for j=1:29
        REC(:,3600*j+1:3600*(j+1)) = REC(:,3600:-1:1);
    end
%     REC = reshape(REC,[LayerSize*500 1]);
    REC = reshape(REC,[size(REC,2)*500 1]);
    
    for j=1:length(REC)
        if (REC(j)>0)
            REC(j) = 1;
        end
    end
    AAA((Zstart-1)*Nx+1:(Zstart+LayerSize-1)*Nx) = REC(1:Nx*LayerSize)';    
    
end
%% ЗАПОЛНЕНИЕ ЧАСТИ ОБЛАСТИ МАТЕРИАЛОМ ТРЕЩИН
if (FRACTYPE == 7)
    AAA((Zstart-1)*Nx+1:(Zend)*Nx)=1;
end
%% ЗАПОЛНЕНИЕ СЛОИСТЫМ МАТЕРИАЛОМ
% ПЕРПЕНДИКУЛЯРНЫЕ СЛОИ
if (FRACTYPE == 6)
    j = Zstart - 1;
    while (j <= Zend - LayerWidth)
       AAA(j*Nx + 1:(j+LayerWidth)*Nx) = 1; 
       j = j + LayerWidth*2;
    end
end

% ПАРАЛЛЕЛЬНЫЕ СЛОИ
if (FRACTYPE == 5)
    i = 1;
    while (i < Nx - LayerWidth)
        for l = 0:LayerWidth-1
           for j = Zstart - 1: Zend - 1
               AAA(j*Nx + i + l) = 1;
           end
        end
        i = i + 2*LayerWidth;
    end
    for j = Zstart: Zend
       AAA(j*Nx) = 1;
    end
end
%%
if (FRACTYPE < 5)
while (FracPointNumber < FracConc*Nx*(Zend - Zstart))  
%% ПЕРЕСЕКАЮЩИЕСЯ ТРЕЩИНЫ 
if (FRACTYPE == 4)
 FracNx = 2;   % число узлов в трещине по оси X
 FracNz = 15;   % число узлов в трещине по оси Z
 Xfrac = round(Nx*rand + 1);     
 Zfrac = floor(unidrnd(Zend - Zstart - FracNz)+Zstart); 
 XFracs(NFracs + 1) = Xfrac - 1;    % записываем координаты хвостов трещин в отдельные массивы
 ZFracs(NFracs + 1) = Zfrac - 1;    % отнимаем 1 чтобы перейти к индексам в C++
 Orient (NFracs + 1) = 1;
 
 for i=Xfrac:Xfrac+FracNx-1
     for j=Zfrac:Zfrac+FracNz-1
         
         if(i==1)
             AAA((j-1)*Nx+1) = 1;
             AAA((j-1)*Nx+Nx) = 1;
         end;
         
         if((i<Nx)&&(i>1))
             AAA((j-1)*Nx+i) = 1;
         end;
         
         %if(i == Nx)
         if (abs(i - Nx) < 1)                         
             AAA((j-1)*Nx+i) = 1;
             AAA((j-2)*Nx+i+1) = 1;
         end
         
         if(i>Nx)             
            AAA((j-2)*Nx+i+1) = 1;
         end            
     end
 end

 FracNx = 15;   % число узлов в трещине по оси X
 FracNz = 2;   % число узлов в трещине по оси Z
 Xfrac = round(rand*FracNx+Xfrac-(FracNx-FracNz));
 if(Xfrac < 1)
    Xfrac = Nx + Xfrac; 
 end
 Zfrac = round(rand*(FracNx-FracNz)+Zfrac);
 XFracs(NFracs + 2) = Xfrac - 1;    % записываем координаты хвостов трещин в отдельные массивы
 ZFracs(NFracs + 2) = Zfrac - 1;    % отнимаем 1 чтобы перейти к индексам в C++
 Orient (NFracs + 2) = 2;
 for i=Xfrac:Xfrac+FracNx-1
     for j=Zfrac:Zfrac+FracNz-1
         if(i==1)
             AAA((j-1)*Nx+1) = 1;
             AAA((j-1)*Nx+Nx) = 1;
         end;
         
         if((i<Nx)&&(i>1))
            
            AAA((j-1)*Nx+i) = 1;
         end;
         
         if (abs(i - Nx) < 1)
                     
            AAA((j-1)*Nx+i) = 1;
            AAA((j-2)*Nx+i+1) = 1;
         end
         
         if(i>Nx) 
            AAA((j-2)*Nx+i+1) = 1;
         end;  
     end
 end
NFracs = NFracs + 2;
end

%% ТРЕЩИНЫ ОДИНАКОВОЙ ОРИЕНТАЦИИ
% ПАРАЛЛЕЛЬНЫЕ ТРЕЩИНЫ
if (FRACTYPE == 1)
FracNx = 2;   % число узлов в трещине по оси X
 FracNz = 15;   % число узлов в трещине по оси Z
ConCheck = 1;
while (ConCheck > 0)
    ConCheck = 0;
 Xfrac = round(Nx*rand + 1);      
 Zfrac = round(((Zend-FracNz)-Zstart)*rand + Zstart);

  for i=Xfrac-1:Xfrac+FracNx
     for j=Zfrac-1:Zfrac+FracNz
         if ((j > 0)&&(j < Nz + 1))
         if(i==0)
             ConCheck = ConCheck + AAA((j-1)*Nx+Nx-1);
         end
         
         if(i==1)
             ConCheck = ConCheck + AAA((j-1)*Nx+1);
             ConCheck = ConCheck + AAA((j-1)*Nx+Nx);
         end
         
         if((i<Nx)&&(i>1))
             ConCheck = ConCheck + AAA((j-1)*Nx+i);
         end
         
         if (abs(i - Nx) < 1)                       
             ConCheck = ConCheck + AAA((j-1)*Nx+i);
             ConCheck = ConCheck + AAA((j-2)*Nx+i+1);
         end
         
         if(i>Nx)              
            ConCheck = ConCheck + AAA((j-2)*Nx+i+1);
         end 
         end
     end
  end
 if (ConCheck == 0)   
 for i=Xfrac:Xfrac+FracNx-1
     for j=Zfrac:Zfrac+FracNz-1
         
         if(i==1)
             AAA((j-1)*Nx+1) = 1;
             AAA((j-1)*Nx+Nx) = 1;
         end;
         
         if((i<Nx)&&(i>1))
             AAA((j-1)*Nx+i) = 1;
         end;
         
         if (abs(i - Nx) < 1)                   
             AAA((j-1)*Nx+i) = 1;
             AAA((j-2)*Nx+i+1) = 1;
         end
         
         if(i>Nx)             
            AAA((j-2)*Nx+i+1) = 1;
         end;            
     end
 end
 end
end
 XFracs(NFracs + 1) = Xfrac - 1;    % записываем координаты хвостов трещин в отдельные массивы
 ZFracs(NFracs + 1) = Zfrac - 1;    % отнимаем 1 чтобы перейти к индексам в C++
 Orient (NFracs + 1) = 1;           % ДЛЯ ГОРИЗОНТ
 NFracs = NFracs + 1;
 for i=Xfrac:Xfrac+FracNx-1
     for j=Zfrac:Zfrac+FracNz-1
         
         if(i==1)
             AAA((j-1)*Nx+1) = 1;
             AAA((j-1)*Nx+Nx) = 1;
         end;
         
         if((i<Nx)&&(i>1))
             AAA((j-1)*Nx+i) = 1;
         end;
         
         if (abs(i - Nx) < 1)                      
             AAA((j-1)*Nx+i) = 1;
             AAA((j-2)*Nx+i+1) = 1;
         end
         
         if(i>Nx)               
            AAA((j-2)*Nx+i+1) = 1;
         end;            
     end
 end
end
%  
%  % ДЛЯ РАВНОМЕРНОГО РАСПРЕДЕЛЕНИЯ ТРЕЩИН ДВУХ ОРИЕНТАЦИИ
%% ВЕРТИКАЛЬНЫЕ ТРЕЩИНЫ
if (FRACTYPE == 2)
 FracNx = 15;   % число узлов в трещине по оси X
 FracNz = 2;   % число узлов в трещине по оси Z
ConCheck = 1;
while (ConCheck > 0)
    ConCheck = 0;
 Xfrac = round(Nx*rand + 1);      
 Zfrac = round(((Zend-FracNz)-Zstart)*rand + Zstart);
 
  for i=Xfrac-1:Xfrac+FracNx
     for j=Zfrac-1:Zfrac+FracNz
         if ((j > 0)&&(j < Nz + 1))
         if(i==0)
             ConCheck = ConCheck + AAA((j-1)*Nx+Nx-1);
         end
         
         if(i==1)
             ConCheck = ConCheck + AAA((j-1)*Nx+1);
             ConCheck = ConCheck + AAA((j-1)*Nx+Nx);
         end
         
         if((i<Nx)&&(i>1))
             ConCheck = ConCheck + AAA((j-1)*Nx+i);
         end
         
         if (abs(i - Nx) < 1)                    
             ConCheck = ConCheck + AAA((j-1)*Nx+i);
             ConCheck = ConCheck + AAA((j-2)*Nx+i+1);
         end
         
         if(i>Nx)              
            ConCheck = ConCheck + AAA((j-2)*Nx+i+1);
         end 
         end
     end
  end
 if (ConCheck == 0)   
 for i=Xfrac:Xfrac+FracNx-1
     for j=Zfrac:Zfrac+FracNz-1
         
         if(i==1)
             AAA((j-1)*Nx+1) = 1;
             AAA((j-1)*Nx+Nx) = 1;
         end;
         
         if((i<Nx)&&(i>1))
             AAA((j-1)*Nx+i) = 1;
         end;
         
         if (abs(i - Nx) < 1)                  
             AAA((j-1)*Nx+i) = 1;
             AAA((j-2)*Nx+i+1) = 1;
         end
         
         if(i>Nx)              
            AAA((j-2)*Nx+i+1) = 1;
         end;            
     end
 end
 end
end
 XFracs(NFracs + 1) = Xfrac - 1;    % записываем координаты хвостов трещин в отдельные массивы
 ZFracs(NFracs + 1) = Zfrac - 1;    % отнимаем 1 чтобы перейти к индексам в C++
 Orient (NFracs + 1) = 2;           % ДЛЯ ВЕРТИК
 NFracs = NFracs + 1;
 for i=Xfrac:Xfrac+FracNx-1
     for j=Zfrac:Zfrac+FracNz-1
         
         if(i==1)
             AAA((j-1)*Nx+1) = 1;
             AAA((j-1)*Nx+Nx) = 1;
         end;
         
         if((i<Nx)&&(i>1))
             AAA((j-1)*Nx+i) = 1;
         end;
         
         if (abs(i - Nx) < 1)                      
             AAA((j-1)*Nx+i) = 1;
             AAA((j-2)*Nx+i+1) = 1;
         end
         
         if(i>Nx)             
            AAA((j-2)*Nx+i+1) = 1;
         end;            
     end
 end
end

%% НЕПЕРЕСЕКАЮЩИЕСЯ ТРЕЩИНЫ
if (FRACTYPE == 3)
FracNx = 2;   % число узлов в трещине по оси X
FracNz = 15;   % число узлов в трещине по оси Z
ConCheck = 1;
while (ConCheck > 0)
    ConCheck = 0;
 Xfrac = round(Nx*rand + 1);      
 Zfrac = round(((Zend-FracNz)-Zstart)*rand + Zstart);
 
  for i=Xfrac-1:Xfrac+FracNx
     for j=Zfrac-1:Zfrac+FracNz
         if ((j > 0)&&(j < Nz + 1))
         if(i==0)
             ConCheck = ConCheck + AAA((j-1)*Nx+Nx-1);
         end
         
         if(i==1)
             ConCheck = ConCheck + AAA((j-1)*Nx+1);
             ConCheck = ConCheck + AAA((j-1)*Nx+Nx);
         end
         
         if((i<Nx)&&(i>1))
             ConCheck = ConCheck + AAA((j-1)*Nx+i);
         end
         
         if (abs(i - Nx) < 1)                       
             ConCheck = ConCheck + AAA((j-1)*Nx+i);
             ConCheck = ConCheck + AAA((j-2)*Nx+i+1);
         end
         
         if(i>Nx)             
            ConCheck = ConCheck + AAA((j-2)*Nx+i+1);
         end 
         end
     end
  end
 if (ConCheck == 0)   
 for i=Xfrac:Xfrac+FracNx-1
     for j=Zfrac:Zfrac+FracNz-1
         
         if(i==1)
             AAA((j-1)*Nx+1) = 1;
             AAA((j-1)*Nx+Nx) = 1;
         end;
         
         if((i<Nx)&&(i>1))
             AAA((j-1)*Nx+i) = 1;
         end;
         
         if (abs(i - Nx) < 1)                        
             AAA((j-1)*Nx+i) = 1;
             AAA((j-2)*Nx+i+1) = 1;
         end
         
         if(i>Nx)              
            AAA((j-2)*Nx+i+1) = 1;
         end;            
     end
 end
 end
end
 XFracs(NFracs + 1) = Xfrac - 1;    % записываем координаты хвостов трещин в отдельные массивы
 ZFracs(NFracs + 1) = Zfrac - 1;    % отнимаем 1 чтобы перейти к индексам в C++
 Orient (NFracs + 1) = 1; 


 FracNx = 15;   % число узлов в трещине по оси X
 FracNz = 2;   % число узлов в трещине по оси Z
ConCheck = 1;
while (ConCheck > 0)
    ConCheck = 0;
 Xfrac = round(Nx*rand + 1);     
 Zfrac = round(((Zend-FracNz)-Zstart)*rand + Zstart);
 
  for i=Xfrac-1:Xfrac+FracNx
     for j=Zfrac-1:Zfrac+FracNz
         if ((j > 0)&&(j < Nz + 1))
         if(i==0)
             ConCheck = ConCheck + AAA((j-1)*Nx+Nx-1);
         end
         
         if(i==1)
             ConCheck = ConCheck + AAA((j-1)*Nx+1);
             ConCheck = ConCheck + AAA((j-1)*Nx+Nx);
         end
         
         if((i<Nx)&&(i>1))
             ConCheck = ConCheck + AAA((j-1)*Nx+i);
         end
         
         if (abs(i - Nx) < 1)                     
             ConCheck = ConCheck + AAA((j-1)*Nx+i);
             ConCheck = ConCheck + AAA((j-2)*Nx+i+1);
         end
         
         if(i>Nx)              
            ConCheck = ConCheck + AAA((j-2)*Nx+i+1);
         end 
         end
     end
  end
 if (ConCheck == 0)
 for i=Xfrac:Xfrac+FracNx-1
     for j=Zfrac:Zfrac+FracNz-1
         
         if(i==1)
             AAA((j-1)*Nx+1) = 1;
             AAA((j-1)*Nx+Nx) = 1;
         end;
         
         if((i<Nx)&&(i>1))
             AAA((j-1)*Nx+i) = 1;
         end;
         
         if (abs(i - Nx) < 1)                       
             AAA((j-1)*Nx+i) = 1;
             AAA((j-2)*Nx+i+1) = 1;
         end
         
         if(i>Nx)             
            AAA((j-2)*Nx+i+1) = 1;
         end;            
     end
 end
 end
end
 XFracs(NFracs + 2) = Xfrac - 1;    % записываем координаты хвостов трещин в отдельные массивы
 ZFracs(NFracs + 2) = Zfrac - 1;    % отнимаем 1 чтобы перейти к индексам в C++
 Orient (NFracs + 2) = 2;
 NFracs = NFracs + 2;
end
%% ПОДСЧЕТ КОНЦЕНТРАЦИИ
FracPointNumber = sum(AAA((Zstart-1)*Nx + 1:Zend*Nx));
end
end

%% ЗАПИСЬ ОРИЕНТАЦИЙ И КООРДИНАТ ТРЕЩИН
if (FRACTYPE < 5)
    %For Orientation
    OrientName = 'Orient.bin'; 
    FOrient = fopen(OrientName,'w'); % massive of 1s and 2s: 1 - horizontal fracture vertex with minimal x and z coordinates, 2 - vertical fracture vertex with minimal x and z coordinates, 0 - otherwise
    fwrite(FOrient,Orient,'int');
    fclose(FOrient); 

    % For fractures x-coordinate
    XFracsName = 'XFracs.bin'; 
    FXFracs = fopen(XFracsName,'w'); % massive of 1s and 2s: 1 - horizontal fracture vertex with minimal x and z coordinates, 2 - vertical fracture vertex with minimal x and z coordinates, 0 - otherwise
    fwrite(FXFracs,XFracs,'int');
    fclose(FXFracs);

    % For fractures z-coordinate
    ZFracsName = 'ZFracs.bin'; 
    FZFracs = fopen(ZFracsName,'w'); % massive of 1s and 2s: 1 - horizontal fracture vertex with minimal x and z coordinates, 2 - vertical fracture vertex with minimal x and z coordinates, 0 - otherwise
    fwrite(FZFracs,ZFracs,'int');
    fclose(FZFracs);

    NFracsName = 'NFracs.bin'; 
    FNFracs = fopen(NFracsName,'w'); % number of fractures
    fwrite(FNFracs,NFracs,'int');
    fclose(FNFracs); 
end
 

%% ЗАПИСЬ ПАРАМЕТРОВ МОДЕЛИ

Fgrid = fopen([folder 'grid.bin'],'w');
fwrite(Fgrid,Nx,'int');
fwrite(Fgrid,Nz,'int');
fwrite(Fgrid,Dx,'float');
fwrite(Fgrid,Dz,'float');
fwrite(Fgrid,Dt,'float');
fclose(Fgrid);

FgridData = fopen([folderData 'grid.bin'],'w');
fwrite(FgridData,Nx,'int');
fwrite(FgridData,Nz,'int');
fwrite(FgridData,Dx,'float');
fwrite(FgridData,Dz,'float');
fwrite(FgridData,Dt,'float');
fclose(FgridData);

%% create pphysical properties files
MatPropExport_ViscoElast(AAA, Frequency, Nx, Nz, folder);
