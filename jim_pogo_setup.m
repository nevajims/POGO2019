
% POGO_model =   setup_POGO_model(reshaped_proc_data,mode_no )

load('rpd_8.mat')


whichFrame = 1;
whichSig = 1;



SAFE_MESH = reshaped_proc_data.mesh;
mesh_details = get_element_details(SAFE_MESH);

mode_no = 4 ;  %   mode 4 is longitudinal


dx = mesh_details.mean_edge_length;
dz = dx ;
nz = ceil(1/dx);    %make it a 1m length

POGO_model = convert_SAFE_POGO(SAFE_MESH , dz , nz , 1 );

% courant = 0.3; % mix of sizes so be conservative with stability   
CFL = 0.3;     %courant

% define material properties
% STEEL PROPERTIES - NOT USED
% E = 210e9;
% G = 80e9;
% nu = E/2/G-1;
% rho = 8000;

% COPPER PROPERTIES
E   =        117e9        ;       % youngs_modulus 
nu  =        0.35         ;       % poissons_ratio
rho =        8960         ;       % density
G   =        E/(2*(1+nu)) ;       % Shear Modulus  

%sound speeds
c0  = sqrt(E*(1-nu)/(rho*(1+nu)*(1-2*nu)))  ; % m/s
cSh = sqrt(E/(2*rho*(1+nu)))                ; % m/s   

POGO_model.dt = (dx/c0) * CFL ;
POGO_model.nt = round(6.6667e-4/POGO_model.dt) ;

% want 30e-6 s long simulation
% length of signal  6 cycles at 100,000 kHz     =   6.0000e-05  s
%  assume its for a 1 m long rod
%  time to go 2m     2/3000          6.667 -4
%  define frequency


chosen_freq = 100e3;
% number of cycles in toneburst
nCyc = 6;

%get number of elements
nEls = size(POGO_model.elNodes,2);


%define the material references (set them all to material 1)
POGO_model.matTypeRefs = ones(nEls,1);
%then define material 1 with the properties from before
POGO_model.matTypes{1}.paramsType = 0;
POGO_model.matTypes{1}.paramValues = [E, nu, rho];

%[m3d] = genPogoHannSignal(m3d, nCyc, freq, whichSig, whichFrame );

POGO_model = recreate_mode_shape(POGO_model , nCyc , chosen_freq , reshaped_proc_data , mode_no , whichFrame ,whichSig);


% how do you apply this signal to the 





%extrudePogo - extrude a 2D Pogo mesh into a 3D model
%   [ m3d ] = extrudePogo( m, nz, dz )
%
%m - 2d input model
%nz - number of desired points in the z (out of plane) direction
%dz - spacing in the z direction
%m3d - 3d output model
%
%z positions run from zero up to (nz-1)*dz


