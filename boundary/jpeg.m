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

flagRec=1;
flagSnapSigma=0;
SnapName='sigma_zz_';
Name='sigma_{zz}, i=';

flagSnapV=0;
SnapNameV='v_z_';
NameV='v_z,i=';

flagParameters=0;

if(flagRec==1)
    field = 'v_z_rec_1';
     frec = fopen([field '.bin'],'r'); % enter name of the file
     REC = fread(frec,[Nx,inf],'float');
     fclose(frec);

      field2 = 'v_z_rec_2';
     frec2 = fopen([field2 '.bin'],'r'); % enter name of the file
     REC2 = fread(frec2,[Nx,inf],'float');
     fclose(frec2);
 
     h=figure;
     t_array=0:dt:time+3/freq+dt;
     plot(REC(1,:),t_array);
     hold on
     plot(REC2(1,:),t_array);
     saveas(h,'trassa', 'png');
     close(h);
     
    s=z_rec_2*dz-z_rec_1*dz;
end

% % %--------------SNAPSHOTS-------------%
if(flagSnapSigma==1)
    for k=1:20
        field=strcat(SnapName,num2str(k));
        name=strcat(Name,num2str(k));
         frec = fopen([field '.bin'],'r'); % enter name of the file
         nx=Nx;
         nz=Nz;
         if(strcmp(SnapName,'sigma_xz_'))
             nx=Nx-1;
             nz=Nz-1;
         end
         
    REC = fread(frec,[nx,inf],'float');
    fclose(frec);
    h=figure;
    dz = dx;
    imagesc((1:nz)*dz, (1:nx)*dx,REC); % for snapshot data
    %axis equal;   
    set(gca,'FontSize',14);
    colormap copper;
    colorbar;
    title(name);
    xlabel('z, m','FontSize',14);
    ylabel('x, m','FontSize',14);
    if(strcmp(SnapName,'sigma_xx_')||strcmp(SnapName,'sigma_zz_')||strcmp(SnapName,'sigma_xz_'))
            caxis([-0.25 0.25]) %for sigma
    end

    if(strcmp(SnapName,'v_x_')||strcmp(SnapName,'v_z_'))
            caxis([-1.5*10^(-8) 1.5*10^(-8)]) %for v_x and v_z
    end
    saveas(h,field, 'png');
    close(h);
    end
end

if(flagSnapV==1)
    for k=1:20
            field=strcat(SnapNameV,num2str(k));
            name=strcat(NameV,num2str(k));
             frec = fopen([field '.bin'],'r'); % enter name of the file
             nx=Nx;
             nz=Nz;
         if(strcmp(SnapNameV,'v_x_'))
                 nx=Nx-1;
         end
         if(strcmp(SnapNameV,'v_z_'))
                 nz=Nz-1;
         end
          REC = fread(frec,[nx,inf],'float');
        fclose(frec);
        h=figure;
        dz = dx;
        imagesc((1:nz)*dz, (1:nx)*dx,REC); % for snapshot data
        %axis equal;   
        set(gca,'FontSize',14);
        colormap copper;
        colorbar;
        title(name);
        xlabel('z, m','FontSize',14);
        ylabel('x, m','FontSize',14);

        if(strcmp(SnapName,'v_x_')||strcmp(SnapName,'v_z_'))
                caxis([-1.5*10^(-8) 1.5*10^(-8)]) %for v_x and v_z
        end
        saveas(h,field, 'png');
        close(h);
    end
end
% %     field = 'sigma_zz_7';
% %     name='\sigma_{zz}, i=7';
%     field = 'v_z_20';
%     name='v_{z}, i=20';
%     frec = fopen([field '.bin'],'r'); % enter name of the file
%     REC = fread(frec,[Nx,inf],'float');
%     fclose(frec);
%     h=figure;
%     dz = dx;
%     imagesc((1:Nz-1)*dz, (1:Nx)*dx,REC); % for snapshot data
%     Lz=Nz*dz;
%     %axis equal;   
%     set(gca,'FontSize',14);
%     colormap copper;
%     colorbar;
%     title(name);
%     xlabel('z, m','FontSize',14);
%     ylabel('x, m','FontSize',14);
% 
%     %caxis([-0.25 0.25]) %for sigma
%     caxis([-1.5*10^(-8) 1.5*10^(-8)]) %for v_x and v_z
% 
%     saveas(h,field, 'png');
%     %close(h);

%-----------PARAMETERS-------------%

% field='tau55_xz';
% frec=fopen([field '.bin'],'r');
% REC=fread(frec,[Nx-1,inf],'float');
% fclose(frec);
% h=figure;
% dz = dx;
% imagesc((1:Nz-1)*dz, (1:Nx-1)*dx,REC); % for snapshot data
% %axis equal;   
% set(gca,'FontSize',14);
% colormap copper;
% colorbar;
% title(field);
% xlabel('z, m','FontSize',14);
% ylabel('x, m','FontSize',14);
%  saveas(h,field, 'png');