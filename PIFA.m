% EXAMPLE / generated file for openEMS from FreeCAD
%
% This is generated file
%
% FreeCAD to OpenEMS plugin by Lubomir Jagos
%

close all
clear
clc

%% setup the simulation
unit = 1e-3; % all length in mm
physical_constants;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                substrate.width
%  _______________________________________________    __ substrate.
% | A                        ifa.l                |\  __    thickness
% | |ifa.e         __________________________     | |
% | |             |    ___  _________________| w2 | |
% | |       ifa.h |   |   ||                      | |
% |_V_____________|___|___||______________________| |
% |                .w1   .wf\                     | |
% |                   |.fp|  \                    | |
% |                       |    feed point         | |
% |                       |                       | | substrate.length
% |<- substrate.width/2 ->|                       | |
% |                                               | |
% |_______________________________________________| |
%  \_______________________________________________\|
%
% Note: It's not checked whether your settings make sense, so check
%       graphical output carefully.
%
substrate.width  = 80;             % width of substrate
substrate.length = 80;             % length of substrate
substrate.thickness = 1.5;         % thickness of substrate
substrate.cells = 4;               % use 4 cells for meshing substrate

ifa.h  = 8;            % height of short circuit stub
ifa.l  = 22.5;         % length of radiating element
ifa.w1 = 4;            % width of short circuit stub
ifa.w2 = 2.5;          % width of radiating element
ifa.wf = 1;            % width of feed element
ifa.fp = 4;            % position of feed element relative to short
                       %  circuit stub
ifa.e  = 10;           % distance to edge

%setup feeding
feed.R = 50;     %feed resistance

%% switches & options...
postprocessing_only = 0;
draw_3d_pattern = 0; % this may take a while...
use_pml = 0;         % use pml boundaries instead of mur

currDir = strrep(pwd(), '\', '\\');
display(currDir);

%LuboJ, JUST TO SEE RESULT
openEMS_opts = '--no-simulation';

%% prepare simulation folder
Sim_Path = 'tmp';
Sim_CSX = 'IFA.xml';

if isfolder(Sim_Path)==1    %If there is a folder
  if length(dir(Sim_Path)) > 2  %if folder contain files then we delete them
    delete(strcat(Sim_Path,'\*'))
  endif
else
  mkdir(Sim_Path); % create empty simulation folder
endif

%% setup FDTD parameter & excitation function
max_timesteps = 60000;%6750;
min_decrement = 0.001; % equivalent to -50 dB
FDTD = InitFDTD( 'NrTS', max_timesteps);%, 'EndCriteria', min_decrement );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXCITATION Gauss_excit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f0 = 2.5e9;
fc = 1.0e9;
FDTD = SetGaussExcite(FDTD, f0, fc );
max_res = c0 / (f0 + fc) / 20;
BC = {"MUR","MUR","MUR","MUR","MUR","MUR"}; % boundary conditions
%BC = {"PML_8" "PML_8" "PML_8" "PML_8" "PML_8" "PML_8"};
FDTD = SetBoundaryCond( FDTD, BC );

CSX = InitCSX();

CSX = AddMetal( CSX, 'PEC' );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MESH variable init
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mesh.x = [];
mesh.y = [];
mesh.z = [];

%mesh.x = [mesh.x (-60:0.5:60)];%outside the F shape but still inside the dielectric
mesh.x = [mesh.x (-60:5:-10) + 0];%outside the F shape but still inside the dielectric
mesh.x = [mesh.x (-10:1:-1) + 0];%upper part of the F shape
mesh.x = [mesh.x (-1:0.5:3) + 0]; %Region of the feed: 1 mesh inside the feed
mesh.x = [mesh.x (3:1:16.5) + 0]; %Lower part of the F shape
mesh.x = [mesh.x (16.5:5:60) + 0]; %

%mesh.y = [mesh.y (-60:0.5:60)];%outside the F shape but still inside the dielectric
mesh.y = [mesh.y (-60:5:28) + 0];
mesh.y = [mesh.y (28:0.5:32) + 0];%Region of the excitation: One mesh inside the excitation
mesh.y = [mesh.y (32:1:60) + 0];

mesh.z = [mesh.z (-70:0.4:70) 1.5];% (1.3:0.1:1.8)];
CSX = DefineRectGrid(CSX, unit, mesh);

##SimBox = [substrate.width*2 substrate.length*2 150];
##%initialize the mesh with the "air-box" dimensions
##mesh.x = [-SimBox(1)/2 SimBox(1)/2];
##mesh.y = [-SimBox(2)/2 SimBox(2)/2];
##mesh.z = [-SimBox(3)/2 SimBox(3)/2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATERIAL - Substrate_PIFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = AddMaterial( CSX, 'Substrate_PIFA' );
CSX = SetMaterialProperty( CSX, 'Substrate_PIFA', 'Epsilon', 4.3, 'Kappa', 0.00058609);
CSX = ImportSTL(CSX, 'Substrate_PIFA',10000, [currDir '/PIFA_Substrate_gen_model.stl'],'Transform',{'Scale', 1});
% add extra cells to discretize the substrate thickness
##mesh.z = [linspace(0,substrate.thickness,substrate.cells+1) mesh.z];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATERIAL - Copper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = AddMetal( CSX, 'Copper' );
CSX = ImportSTL(CSX, 'Copper',9900, [currDir '/PIFA_Copper_gen_model.stl'],'Transform',{'Scale', 1});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATERIAL - GND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = AddMetal( CSX, 'GND' );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRID PRIORITIES GENERATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

##ifa_mesh = DetectEdges(CSX, [], 'SetProperty','Copper');
##mesh.x = [mesh.x SmoothMeshLines(ifa_mesh.x, 0.05)];
##mesh.y = [mesh.y SmoothMeshLines(ifa_mesh.y, 0.05)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LUMPED PART Gauss_excit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% apply the excitation & resist as a current source
tl = [0,substrate.length/2-ifa.e,substrate.thickness];   % translate
start = [0 0 0] + tl;
stop  = start + [ifa.wf 0.5 0];
[CSX port] = AddLumpedPort(CSX, 5 ,1 ,feed.R, start, stop, [0 1 0], true);

%% finalize the mesh
% generate a smooth mesh with max. cell size: lambda_min / 20
##mesh = DetectEdges(CSX, mesh);
##mesh = SmoothMesh(mesh, c0 / (f0+fc) / unit / 20);
##CSX = DefineRectGrid(CSX, unit, mesh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NF2FF PROBES GRIDLINES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );
CSXGeomPlot( [Sim_Path '/' Sim_CSX] );

if (postprocessing_only==0)
    %% run openEMS
    RunOpenEMS( Sim_Path, Sim_CSX);
end

%% postprocessing & do the plots
freq = linspace( f0-fc, f0+fc, 501 );
port = calcPort(port, Sim_Path, freq);

Zin = port.uf.tot ./ port.if.tot;
s11 = port.uf.ref ./ port.uf.inc;
P_in = real(0.5 * port.uf.tot .* conj( port.if.tot )); % antenna feed power

% plot feed point impedance
figure
plot( freq/1e6, real(Zin), 'k-', 'Linewidth', 2 );
hold on
grid on
plot( freq/1e6, imag(Zin), 'r--', 'Linewidth', 2 );
title( 'feed point impedance' );
xlabel( 'frequency f / MHz' );
ylabel( 'impedance Z_{in} / Ohm' );
legend( 'real', 'imag' );

% plot reflection coefficient S11
figure
plot( freq/1e6, 20*log10(abs(s11)), 'k-', 'Linewidth', 2 );
hold on
grid on
title( 'reflection coefficient S_{11}' );
xlabel( 'frequency f / MHz' );
ylabel( 'reflection coefficient |S_{11}|' );
legend( 'S11' );

drawnow

%% NFFF contour plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find resonance frequncy from s11
f_res_ind = find(s11==min(s11));
f_res = freq(f_res_ind);

%%
disp( 'calculating 3D far field pattern and dumping to vtk (use Paraview to visualize)...' );
thetaRange = (0:2:180);
phiRange = (0:2:360) - 180;
nf2ff = CalcNF2FF(nf2ff, Sim_Path, f_res, thetaRange*pi/180, phiRange*pi/180,'Verbose',1,'Outfile','3D_Pattern.h5');

figure
plotFF3D(nf2ff)
hold on
grid on

% display power and directivity
disp( ['radiated power: Prad = ' num2str(nf2ff.Prad) ' Watt']);
disp( ['directivity: Dmax = ' num2str(nf2ff.Dmax) ' (' num2str(10*log10(nf2ff.Dmax)) ' dBi)'] );
disp( ['efficiency: nu_rad = ' num2str(100*nf2ff.Prad./real(P_in(f_res_ind))) ' %']);

E_far_normalized = nf2ff.E_norm{1} / max(nf2ff.E_norm{1}(:)) * nf2ff.Dmax;
DumpFF2VTK([Sim_Path '/3D_Pattern.vtk'],E_far_normalized,thetaRange,phiRange,1e-3);

