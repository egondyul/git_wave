function MatPropExport_ViscoElast(AAA,Frequency, Nx,Nz,folder)
ax = (Nx-1)*Nz;
az = Nx*(Nz-1);
axz = (Nx-1)*(Nz-1);
nn = Nx*Nz
%% properties (1 - background, 2 - fractures)
 tau_sigma = 0;%10^(-2)/(4*pi*Frequency);

% background
rho(1) = 2500;
c11(1) = 4*10^9;
c33(1) = 4*10^9;
c55(1) = 1*10^9;
c13(1) = 2*10^9;
tau11(1) = 0;%10^(-10);
tau33(1) = 0;%10^(-10);
tau55(1) = 0;%10^(-10);
tau13(1) = 0;%10^(-10);

%fractures
rho(2) = 2500;
c11(2) = 4*10^9;
c33(2) = 4*10^9;
c55(2) = 1*10^9;
c13(2) = 2*10^9;
tau11(2) = 0;%10^(-11);
tau33(2) = 0;%10^(-11);
tau55(2) = 0;%10^(-11);
tau13(2) = 0;%10^(-11);
%% save files
Frho_f = fopen([folder 'c11.bin'],'w');
fwrite(Frho_f,(1-AAA).*ones(1,nn)*c11(1) + AAA.*ones(1,nn)*c11(2),'float');
fclose(Frho_f);

Frho_f = fopen([folder 'c33.bin'],'w');
fwrite(Frho_f,(1-AAA).*ones(1,nn)*c33(1) + AAA.*ones(1,nn)*c33(2),'float');
fclose(Frho_f);

Frho_f = fopen([folder 'c13.bin'],'w');
fwrite(Frho_f,(1-AAA).*ones(1,nn)*c13(1) + AAA.*ones(1,nn)*c13(2),'float');
fclose(Frho_f);

Frho_f = fopen([folder 'tau11.bin'],'w');
fwrite(Frho_f,(1-AAA).*ones(1,nn)*tau11(1) + AAA.*ones(1,nn)*tau11(2),'float');
fclose(Frho_f);

Frho_f = fopen([folder 'tau33.bin'],'w');
fwrite(Frho_f,(1-AAA).*ones(1,nn)*tau33(1) + AAA.*ones(1,nn)*tau33(2),'float');
fclose(Frho_f);

Frho_f = fopen([folder 'tau13.bin'],'w');
fwrite(Frho_f,(1-AAA).*ones(1,nn)*tau13(1) + AAA.*ones(1,nn)*tau13(2),'float');
fclose(Frho_f);

Frho_f = fopen([folder 'tau_sigma.bin'],'w');
fwrite(Frho_f,tau_sigma,'float');
fclose(Frho_f);

%% average properties
param = rho;
param_model = (1-AAA).*ones(1,nn)*param(1) + AAA.*ones(1,nn)*param(2);
param_model = reshape(param_model,Nx,Nz);
param_averX =  param_model(2:Nx,:) + param_model(1:Nx-1,:);
param_averX = param_averX * 0.5;
F = fopen([folder 'rho_x.bin'],'w');
fwrite(F,reshape(param_averX,[1,ax]),'float');
fclose(F);
param_averZ =  param_model(:,2:Nz) + param_model(:,1:Nz-1);
param_averZ = param_averZ * 0.5;
F = fopen([folder 'rho_z.bin'],'w');
fwrite(F,reshape(param_averZ,[1,az]),'float');
fclose(F);

param = c55;
param_model = (1-AAA).*ones(1,nn)*param(1) + AAA.*ones(1,nn)*param(2);
param_model = reshape(param_model,Nx,Nz);
param_averZ=zeros(Nx,Nz-1);
for i=1:Nz-1
   for j=1:Nx
      tmp1=param_model(j,i); 
      tmp2=param_model(j,i+1);
      if(tmp1~=0&&tmp2~=0)
         param_averZ(j,i)=1/tmp1+1/tmp2; 
      end
   end
end
param_averXz=zeros(Nx-1,Nz-1);
for i=1:Nz-1
   for j=1:Nx-1
       tmp1=param_averZ(j,i);
       tmp2=param_averZ(j+1,i);
      param_averXz(j,i)= 1/tmp1+1/tmp2;
   end
end
param_averXz=4./param_averXz;

F = fopen([folder 'c55_xz.bin'],'w');
axz = (Nx-1)*(Nz-1);
fwrite(F,reshape(param_averXz,[1,axz]),'float');
fclose(F);

param = tau55;
param_model = (1-AAA).*ones(1,nn)*param(1) + AAA.*ones(1,nn)*param(2);
param_model = reshape(param_model,Nx,Nz);
param_averZ=zeros(Nx,Nz-1);
for i=1:Nz-1
   for j=1:Nx
      tmp1=param_model(j,i); 
      tmp2=param_model(j,i+1);
      if(tmp1~=0&&tmp2~=0)
         param_averZ(j,i)=1/tmp1+1/tmp2; 
      end
   end
end
param_averXz=zeros(Nx-1,Nz-1);
for i=1:Nz-1
   for j=1:Nx-1
       tmp1=param_averZ(j,i);
       tmp2=param_averZ(j+1,i);
      param_averXz(j,i)= 1/tmp1+1/tmp2;
   end
end
param_averXz=4./param_averXz;

F = fopen([folder 'tau55_xz.bin'],'w');
axz = (Nx-1)*(Nz-1);
fwrite(F,reshape(param_averXz,[1,axz]),'float');
fclose(F);

% param = c55;
% param_model = (1-AAA).*ones(1,nn)*param(1) + AAA.*ones(1,nn)*param(2);
% param_model = reshape(param_model,Nx,Nz);
% param_averX =  1./param_model(2:Nx,:) + 1./param_model(1:Nx-1,:);
% param_averXZ = 1./param_averX(:,2:Nz) + 1./param_averX(:,1:Nz-1);
% param_averXZ = 4./param_averXZ;
% F = fopen([folder 'c55_xz.bin'],'w');
% fwrite(F,reshape(param_averXZ,[1,axz]),'float');
% fclose(F);
% 
% param = tau55;
% param_model = (1-AAA).*ones(1,nn)*param(1) + AAA.*ones(1,nn)*param(2);
% param_model = reshape(param_model,Nx,Nz);
% param_averX =  1./param_model(2:Nx,:) + 1./param_model(1:Nx-1,:);
% param_averXZ = 1./param_averX(:,2:Nz) + 1./param_averX(:,1:Nz-1);
% param_averXZ = 4./param_averXZ;
% F = fopen([folder 'tau55_xz.bin'],'w');
% fwrite(F,reshape(param_averXZ,[1,axz]),'float');
% fclose(F);
