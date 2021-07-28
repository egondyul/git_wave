function MatPropExport_ViscoElast(AAA,Frequency, Nx,Nz,folder)
ax = (Nx-1)*Nz;
az = Nx*(Nz-1);
axz = (Nx-1)*(Nz-1);
nn = Nx*Nz
%% properties (1 - background, 2 - fractures)
 tau_sigma = 0;%10^(-2)/(4*pi*Frequency)
nu=0.1;
tau0=0;%nu*(4*pi*Frequency)/(10^(-2))
% % background
% rho(1) = 2450;
% Vp=6000;
% Vssolid=4000;
% c55(1) = Vssolid^2*rho(1)
% c11(1) = Vp^2*rho(1)-c55(1)
% c33(1) = c11(1);
% c55(1) = c11(1);
% c13(1) = c11(1);
% tau11(1) = 0;%10^(-10);
% tau33(1) = 0;%10^(-10);
% tau55(1) = 0;%10^(-10);
% tau13(1) = 0;%10^(-10);
% 
% %Vp=sqrt(c11(1)/rho(1));
% 
% Vpf=1500;
% %fractures
% rho(2) = 1000;
% c11(2) = Vpf^2*rho(2);
% c33(2) = Vpf^2*rho(2);
% c55(2) = 0;
% c13(2) = Vpf^2*rho(2);
% tau11(2) = tau0;
% tau33(2)=tau0;
% tau55(2) = tau0;
% tau13(2) = 0;

% rho(1) = 2500;
% c11(1) = 4*10^9;
% c33(1) = 4*10^9;
% c55(1) = 1*10^9;
% c13(1) = 2*10^9;
% tau11(1) = 0;%10^(-10);
% tau33(1) = 0;%10^(-10);
% tau55(1) = 0;%10^(-10);
% tau13(1) = 0;%10^(-10);
% 
% Vp=sqrt(c11(1)/rho(1));

Vpf=1500;
%fractures
rho(2) = 1000;
c11(2) = Vpf^2*rho(2);
c33(2) = Vpf^2*rho(2);
c55(2) = 0;
c13(2) =Vpf^2*rho(2);
tau11(2) = 2*tau0;
tau33(2)=2*tau0;
tau55(2) = tau0;
tau13(2) = 0;

rho(1) = 1000;
c11(1) = Vpf^2*rho(2);
c33(1) = Vpf^2*rho(2);
c55(1) = 0;
c13(1) = Vpf^2*rho(2);
tau11(1) = tau0;%10^(-10);
tau33(1) = tau0;%10^(-10);
tau55(1) = tau0;%10^(-10);
tau13(1) = 0;%10^(-10);

Vp=sqrt(c11(1)/rho(1));

%% save files
folderData=strcat(folder, 'Data\');
FINPUTData=fopen([folderData 'Vp.txt'],'w');
formatSpec = '%f               # Vp                estimation of velocity\r\n';
fprintf(FINPUTData, formatSpec, Vp);
fclose(FINPUTData);

Frho_f = fopen([folder 'c11.bin'],'w');
fwrite(Frho_f,(1-AAA).*ones(1,nn)*c11(1) + AAA.*ones(1,nn)*c11(2),'double');
fclose(Frho_f);

Frho_f = fopen([folder 'c33.bin'],'w');
fwrite(Frho_f,(1-AAA).*ones(1,nn)*c33(1) + AAA.*ones(1,nn)*c33(2),'double');
fclose(Frho_f);

Frho_f = fopen([folder 'c13.bin'],'w');
fwrite(Frho_f,(1-AAA).*ones(1,nn)*c13(1) + AAA.*ones(1,nn)*c13(2),'double');
fclose(Frho_f);

Frho_f = fopen([folder 'tau11.bin'],'w');
fwrite(Frho_f,(1-AAA).*ones(1,nn)*tau11(1) + AAA.*ones(1,nn)*tau11(2),'double');
fclose(Frho_f);

Frho_f = fopen([folder 'tau33.bin'],'w');
fwrite(Frho_f,(1-AAA).*ones(1,nn)*tau33(1) + AAA.*ones(1,nn)*tau33(2),'double');
fclose(Frho_f);

Frho_f = fopen([folder 'tau13.bin'],'w');
fwrite(Frho_f,(1-AAA).*ones(1,nn)*tau13(1) + AAA.*ones(1,nn)*tau13(2),'double');
fclose(Frho_f);

Frho_f = fopen([folder 'tau_sigma.bin'],'w');
fwrite(Frho_f,tau_sigma,'double');
fclose(Frho_f);

%% average properties
param = rho;
param_model = (1-AAA).*ones(1,nn)*param(1) + AAA.*ones(1,nn)*param(2);
param_model = reshape(param_model,Nx,Nz);
param_averX =  param_model(2:Nx,:) + param_model(1:Nx-1,:);
param_averX = param_averX * 0.5;
F = fopen([folder 'rho_x.bin'],'w');
fwrite(F,reshape(param_averX,[1,ax]),'double');
fclose(F);
param_averZ =  param_model(:,2:Nz) + param_model(:,1:Nz-1);
param_averZ = param_averZ * 0.5;
F = fopen([folder 'rho_z.bin'],'w');
fwrite(F,reshape(param_averZ,[1,az]),'double');
fclose(F);

figure;
  imagesc((1:Nz)*0.002, (1:Nx-1)*0.002,param_averX); 
  
param = c55;
param_model = (1-AAA).*ones(1,nn)*param(1) + AAA.*ones(1,nn)*param(2);
param_model = reshape(param_model,Nx,Nz);
param_aver=zeros(Nx-1,Nz-1);
for i=1:Nz-1
    for j=1:Nx-1
        tmp1=param_model(j,i);
        tmp2=param_model(j+1,i);
        tmp3=param_model(j,i+1);
        tmp4=param_model(j+1,i+1);
        if(tmp1~=0&&tmp2~=0&&tmp3~=0&&tmp4~=0)
            param_aver(j,i)=1/tmp1+1/tmp2+1/tmp3+1/tmp4;
        end
    end
end
param_averXz_c55=zeros(Nx-1,Nz-1);
for i=1:Nz-1
    for j=1:Nx-1
        if(param_aver(j,i)~=0)
            param_averXz_c55(j,i)=4/param_aver(j,i);
        end
    end
end
figure;
  imagesc((1:Nz-1)*0.002, (1:Nx-1)*0.002,param_averXz_c55); 

% param_averZ=zeros(Nx,Nz-1);
% for i=1:Nz-1
%    for j=1:Nx
%       tmp1=param_model(j,i); 
%       tmp2=param_model(j,i+1);
%       if(tmp1~=0&&tmp2~=0)
%          param_averZ(j,i)=1/tmp1+1/tmp2; 
%       end
%    end
% end
% param_averXz=zeros(Nx-1,Nz-1);
% for i=1:Nz-1
%    for j=1:Nx-1
%        tmp1=param_averZ(j,i);
%        tmp2=param_averZ(j+1,i);
%       param_averXz(j,i)= 1/tmp1+1/tmp2;
%    end
% end
% param_averXz=4./param_averXz;

F = fopen([folder 'c55_xz.bin'],'w');
axz = (Nx-1)*(Nz-1);
fwrite(F,reshape(param_averXz_c55,[1,axz]),'double');
fclose(F);

param = tau55;
param_model = (1-AAA).*ones(1,nn)*param(1) + AAA.*ones(1,nn)*param(2);
param_model = reshape(param_model,Nx,Nz);
param_aver=zeros(Nx-1,Nz-1);
for i=1:Nz-1
    for j=1:Nx-1
        tmp1=param_model(j,i);
        tmp2=param_model(j+1,i);
        tmp3=param_model(j,i+1);
        tmp4=param_model(j+1,i+1);
        if(tmp1~=0&&tmp2~=0&&tmp3~=0&&tmp4~=0)
            param_aver(j,i)=1/tmp1+1/tmp2+1/tmp3+1/tmp4;
        end
    end
end
param_averXz_tau55=zeros(Nx-1,Nz-1);
for i=1:Nz-1
    for j=1:Nx-1
         if(param_aver(j,i)~=0)
            param_averXz_tau55(j,i)=4/param_aver(j,i);
        end
    end
end

for j=2:Nz-2
    for i=2:Nx-2
         if(param_averXz_c55(i,j)==0)
            if(param_averXz_c55(i-1,j)~=0&&param_averXz_c55(i,j-1)~=0&&param_averXz_c55(i,j+1)==0)
                param_averXz_tau55(i,j)=param_averXz_tau55(i+1,j+1);
            elseif(param_averXz_c55(i-1,j)~=0&&param_averXz_c55(i,j-1)==0&&param_averXz_c55(i,j+1)==0)
                param_averXz_tau55(i,j)=param_averXz_tau55(i+1,j);
            elseif(param_averXz_c55(i-1,j)~=0&&param_averXz_c55(i,j+1)~=0&&param_averXz_c55(i,j-1)==0)
                param_averXz_tau55(i,j)=param_averXz_tau55(i+1,j-1);
            elseif(param_averXz_c55(i+1,j)~=0&&param_averXz_c55(i,j-1)~=0&&param_averXz_c55(i,j+1)==0)
                param_averXz_tau55(i,j)=param_averXz_tau55(i-1,j+1);
            elseif(param_averXz_c55(i+1,j)~=0&&param_averXz_c55(i,j+1)==0&&param_averXz_c55(i,j-1)==0)
                param_averXz_tau55(i,j)=param_averXz_tau55(i-1,j);
            elseif(param_averXz_c55(i+1,j)~=0&&param_averXz_c55(i,j+1)~=0&&param_averXz_c55(i,j-1)==0)
                param_averXz_tau55(i,j)=param_averXz_tau55(i-1,j-1);
            elseif(param_averXz_c55(i-1,j)==0&&param_averXz_c55(i+1,j)==0&&param_averXz_c55(i,j+1)~=0&&param_averXz_c55(i,j-1)==0)
                param_averXz_tau55(i,j)=param_averXz_tau55(i,j-1);
            elseif(param_averXz_c55(i-1,j)==0&&param_averXz_c55(i+1,j)==0&&param_averXz_c55(i,j+1)==0&&param_averXz_c55(i,j-1)~=0)
                param_averXz_tau55(i,j)=param_averXz_tau55(i,j+1);
            end
        end
    end
end
figure;
  imagesc((1:Nz-1)*0.002, (1:Nx-1)*0.002,param_averXz_tau55); 

% param_averZ=zeros(Nx,Nz-1);
% for i=1:Nz-1
%    for j=1:Nx
%       tmp1=param_model(j,i); 
%       tmp2=param_model(j,i+1);
%       if(tmp1~=0&&tmp2~=0)
%          param_averZ(j,i)=1/tmp1+1/tmp2; 
%       end
%    end
% end
% param_averXz=zeros(Nx-1,Nz-1);
% for i=1:Nz-1
%    for j=1:Nx-1
%        tmp1=param_averZ(j,i);
%        tmp2=param_averZ(j+1,i);
%       param_averXz(j,i)= 1/tmp1+1/tmp2;
%    end
% end
% param_averXz=4./param_averXz;

F = fopen([folder 'tau55_xz.bin'],'w');
axz = (Nx-1)*(Nz-1);
fwrite(F,reshape(param_averXz_tau55,[1,axz]),'double');
fclose(F);
