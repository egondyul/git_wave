function out=ResExport(freq, metcase, path)

out.freq = freq; 
path1=strcat(path,num2str(metcase), '/', num2str(freq),'/Data/');


%% read grid data
Fn = 'grid.bin';
F = fopen([path1 Fn],'r');
Nx = fread(F,1,'int');
trash = fread(F,1,'int');
trash = fread(F,1,'double');
dz = fread(F,1,'double');
dt = fread(F,1,'double');
fclose(F);

out.dz = dz;

%% read source and receivers coordinates
Fn = [path1 'INPUT.txt'];
F = fopen(Fn);
tline = fgetl(F);
tline = fgetl(F);
tline = fgetl(F);
out.src = fscanf(F, '%d');
tline = fgetl(F);
out.rec1 = fscanf(F, '%d');
tline = fgetl(F);
out.rec2 = fscanf(F, '%d');
tline = fgetl(F);
out.WaveLength = fscanf(F, '%d');
fclose(F);

PML = out.WaveLength;
LayerSize = out.rec2 - out.rec1 - 3*out.WaveLength;
out.src = out.src*dz;
out.rec1 = out.rec1*dz;
out.rec2 = out.rec2*dz;
out.interface = (PML + out.WaveLength*4)*dz;
out.LayerWidth = LayerSize*dz;


%% read the data (v_z)
p_rec_1 = fopen([path1 'v_z_rec_1.bin'],'r');     
PREC1 = fread(p_rec_1,[Nx,inf],'double');
fclose(p_rec_1);

p_rec_2 = fopen([path1 'v_z_rec_2.bin'],'r');
PREC2 = fread(p_rec_2,[Nx,inf],'double');
fclose(p_rec_2);

out.time = ((1:1:size(PREC1,2))*dt).';
    
%% average the traces
nn = Nx-1;
out.recnum = nn;
MP1 = mean(PREC1(1:nn,:));      
MP2 = mean(PREC2(1:nn,:));

Fn = [path1 'Vp.txt'];
F = fopen(Fn);
out.Vp = fscanf(F, '%f');
tline = fgetl(F);
fclose(F);

out.trace1=MP1.';
out.trace2=MP2.';
 