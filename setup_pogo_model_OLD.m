
function POGO_model =   setup_pogo_model(reshaped_proc_data , mode_no, rod_length_m, chosen_freq, nCyc )

% create a notes txt file   with input parameters in and save in the
% directory
% make the save directory an input parameter




% mode_no = 4 ;  %   mode 4 is longitudinal in safe model
% load('rpd_8.mat')
% chosen_freq = 100e3   in Hz
% number of cycles in toneburst
% nCyc = 6                      ;
%   put the 3 material properties in here 2

whichFrame = 1;
whichSig = 1;   % 0 for force 1 for displ

SAFE_MESH = reshaped_proc_data.mesh;
mesh_details = get_element_details(SAFE_MESH);
dx = mesh_details.mean_edge_length;
dz = dx ;
nz = ceil(rod_length_m/dx);    %make it a 1m length

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
c0  = sqrt(E*(1-nu)/(rho*(1+nu)*(1-2*nu))) ; % m/s
disp(num2str(c0))

cSh = sqrt(E/(2*rho*(1+nu)))               ; % m/s   

POGO_model.dt = (dx/c0) * CFL                 ;


POGO_model.nt = round(6.6667e-4/POGO_model.dt);
%  time to go 2m     2/3000          6.667 -4
% want 30e-6 s long simulation
% length of signal  6 cycles at 100,000 kHz     =   6.0000e-05  s
%  assume its for a 1 m long rod
%  define frequency

%get number of elements
nEls = size(POGO_model.elNodes,2);

rPoint1 = 1;
rPoint2  = floor(size(POGO_model.nodePos,2)/3) ;
rPoint3  = floor(size(POGO_model.nodePos,2)/2) ;
rPoint4  = floor(2*size(POGO_model.nodePos,2)/3) ;
rPoint5  = size(POGO_model.nodePos,2);

POGO_model = recreate_mode_shape(POGO_model , nCyc , chosen_freq , reshaped_proc_data , mode_no , whichFrame ,whichSig);
%define the material references (set them all to material 1)
POGO_model.matTypeRefs = ones(nEls,1);
%then define material 1 with the properties from before
POGO_model.matTypes{1}.paramsType = 0;
POGO_model.matTypes{1}.paramValues = [E, nu, rho];
%[m3d] = genPogoHannSignal(m3d, nCyc, freq, whichSig, whichFrame );
POGO_model.nDims = 3 ;
POGO_model.elTypes = { struct('name','C3D6R', 'paramsType', 0) };
POGO_model.elTypeRefs = ones(nEls,1);
POGO_model.nDofPerNode = 3;
POGO_model.measFreq = 2;
POGO_model.measStart = 1;

% just get al thge nodes in the middle plane of nodes

[~,mid_index] =   min(abs(POGO_model.nodePos(3,:)-max(POGO_model.nodePos(3,:))/2));
mid_value = POGO_model.nodePos(3,mid_index);
mid_node_vals =  find(POGO_model.nodePos(3,:) == mid_value);
front_node_vals =  find(POGO_model.nodePos(3,:) == 0);
end_node_vals =  find(POGO_model.nodePos(3,:) == max(POGO_model.nodePos(3,:)));



measset_vals = [front_node_vals,mid_node_vals,end_node_vals];

POGO_model.measSets{1}.name = 'main';
POGO_model.measSets{1}.measNodes = sort([measset_vals,measset_vals,measset_vals]);  
%POGO_model.measSets{1}.measNodes = [rPoint1, rPoint1, rPoint1, rPoint2, rPoint2, rPoint2];
POGO_model.measSets{1}.measDof = repmat([1,2,3],[1,length(measset_vals)]);

POGO_model.fieldStoreIncs = 1:25:100*floor (POGO_model.nt/100);  

savePogoInp('\\nas\Public\POGO\CW5\model.pogo-inp',POGO_model);
save( '\\nas\Public\POGO\CW5\POGO_model',  'POGO_model')

%savePogoInp('model.pogo-inp', POGO_model);


end %function POGO_model =   setup_POGO_model(reshaped_proc_data,mode_no )






