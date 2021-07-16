%-------------rEC------------%
Frequency = 10000 ;
Fgrid = fopen('grid.bin','r');
Nx = fread(Fgrid,1,'int');
Nz = fread(Fgrid,1,'int');
dx = fread(Fgrid,1,'float');
dz = fread(Fgrid,1,'float');
dt = fread(Fgrid,1,'float');
fclose(Fgrid);
    
FINPUT = fopen('INPUT.txt','r');
freq = fscanf(FINPUT,'%d');
trash = fgetl(FINPUT);
time = fscanf(FINPUT,'%f');
trash=fgetl(FINPUT);
src=fscanf(FINPUT,'%d');
trash=fgetl(FINPUT);
z_src=fscanf(FINPUT,'%d');
trash=fgetl(FINPUT);
z_rec_1=fscanf(FINPUT,'%d');
trash=fgetl(FINPUT);
z_rec_2=fscanf(FINPUT,'%d');
fclose(FINPUT);
snpTimes = (1:20)*time/20;

 field = 'v_z_rec_1';
 frec = fopen([field '.bin'],'r'); % enter name of the file
 REC = fread(frec,[Nx,inf],'float');
 fclose(frec);
 
  field2 = 'v_z_rec_2';
 frec2 = fopen([field2 '.bin'],'r'); % enter name of the file
 REC2 = fread(frec2,[Nx,inf],'float');
 fclose(frec2);
 
 t_array=0:dt:time+3/freq+dt;
 plot(REC(1,:),t_array);
 hold on
plot(REC2(1,:),t_array);
 
s=z_rec_2*dz-z_rec_1*dz;
 
% %--------------SNAPSHOTS-------------%
% 
    field = 'sigma_xx_20';
    name='\sigma_{xx}, i=20';
    frec = fopen([field '.bin'],'r'); % enter name of the file
    REC = fread(frec,[Nx,inf],'float');
    fclose(frec);
    h=figure;
    dz = dx;
    imagesc((1:Nz)*dz, (1:Nx)*dx,REC); % for snapshot data
    Lz=Nz*dz;
    %axis equal;   
    set(gca,'FontSize',14);
    colormap copper;
    colorbar;
    title(name);
    xlabel('z, m','FontSize',14);
    ylabel('x, m','FontSize',14);

    caxis([-0.25 0.25]) %for sigma
    %caxis([-3.8*10^(-8) 3.8*10^(-8)]) %for v_x and v_z

    saveas(h,field, 'png');
    %close(h);

%-----------PARAMETERS-------------%

% field='rho_x';
% frec=fopen([field '.bin'],'r');
% REC=fread(frec,[Nx-1,inf],'float');
% fclose(frec);
% h=figure;
% dz = dx;
% imagesc((1:Nz)*dz, (1:Nx)*dx,REC); % for snapshot data
% %axis equal;   
% set(gca,'FontSize',14);
% colormap copper;
% colorbar;
% title(field);
% xlabel('z, m','FontSize',14);
% ylabel('x, m','FontSize',14);
%  saveas(h,field, 'png');