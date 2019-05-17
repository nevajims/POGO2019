
% create a specific signal and then apply to a specific node and DOF (degree of Freedon)
% plot it out
 


load('rpd_8.mat')
m = convertJMesh2PogoMesh(reshaped_proc_data);
mesh_details = get_element_details(reshaped_proc_data.mesh);

nNodes1 = size(m.nodePos,2);

dx = mesh_details.mean_edge_length;
dz = dx ;
disp(['mean element size = ',num2str( dx*1E3),' mm.'])
CFL = 0.3;     %courant

% COPPER PROPERTIES
E                  =        117e9        ;       % youngs_modulus 
nu                 =        0.35         ;       % poissons_ratio
rho                =        8960         ;       % density
G                  =        E/(2*(1+nu)) ;       % Shear Modulus  
Length_to_travel_m =        2            ;       % distance the wave needs to travel during the simulation

%sound speeds
c0  = sqrt(E*(1-nu)/(rho*(1+nu)*(1-2*nu)))  ; % m/s
cSh = sqrt(E/(2*rho*(1+nu)))                ; % m/s   

disp(['c0 = ', num2str(c0),', cSh = ', num2str(cSh) , '.'])

model_.dt = (dx/c0) * CFL ;

%  want 30e-6 s long simulation
%  length of signal  6 cycles at 100,000 kHz     =   6.0000e-05  s
%  assume its for a 1 m long rod
%  time to go 2m     2/3000          6.667 -4
%  define frequency

total_time_of_simulation_s =  Length_to_travel_m/c0    ;
model_.nt = round(6.6667e-4/model_.dt)                 ;


whichSig     = 1       ;   %   These are the default values (will have to increment when multiple signals are create
whichFrame   = 1       ;   %   These are the default values  
nCyc         = 6       ;   %
freq         = 100E3   ;   %   100 kHz

Window_type  = 1    ;          % Hanning window
%model_.nt           ;           % this is the total number of points not just the amount in th cycle
%model_.dt           ;

[model_] = genPogoHannSignal(model_, nCyc, freq, whichSig, whichFrame );

aa  = model_.frames{whichFrame}.sigs{whichSig}.sig  ;   %   
nt_ = model_.frames{whichFrame}.ntSig                 ;   % these should not be unecessary
dt_ = model_.frames{whichFrame}.dtSig                 ;   % these should not be unecessary
 
disp(['size of signal = ' ,num2str(length(aa))])

figure (12)
plot( [0:dt_:dt_*(nt_-1)],aa,'.')



