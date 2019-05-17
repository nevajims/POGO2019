%  make the excitation a  fraction of the maximum value in the x / y and z
%  directions  
%--------------------------------------------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------------------------------------------
% compare the signal from the 'extract excitation with the signal from 


file_name        = 'rpd_16.mat'  ;        %    the safe solution structure
window_type      = 1             ;        %  
mode_no          = 4             ;        %    4   = longitudinal
node_plot_index  = 1             ;        %    which node to apply the excitation to
doplot           = 1             ;        %    plot for the extract excitation function 
x_y_z            = 1             ;        %    x_y_z = 1(x)              x_y_z = 2(y)            x_y_z = 3(z)
points_per_trace = 110;


load(file_name)

model_old = convertJMesh2PogoMesh(reshaped_proc_data);
mesh_details = get_element_details(reshaped_proc_data.mesh);

nNodes1 = size(model_old.nodePos,2);

sigType  = 1;     %  force (1 = displacemnt)


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
model_old.dt = (dx/c0) * CFL ;


%  want 30e-6 s long simulation
%  length of signal  6 cycles at 100,000 kHz     =   6.0000e-05  s
%  assume its for a 1 m long rod
%  time to go 2m     2/3000          6.667 -4
%  define frequency
total_time_of_simulation_s =  Length_to_travel_m/c0    ;
model_old.nt = round(6.6667e-4/model_old.dt)                 ;

whichSig         = 1                 ;   %   These are the default values (will have to increment when multiple signals are create
whichFrame       = 1                 ;   %   These are the default values  
chosen_freq      = 100E3             ;   %   hz
chosen_freq_kHz  = chosen_freq/1E3   ;   %   khz
N_cycles         = 6       ;   %

%  need the mode shape from the   /  x_y_z  / node number /   mode_no  / chosen_freq

if  chosen_freq > min (reshaped_proc_data.freq(:,mode_no)) && chosen_freq <  max (reshaped_proc_data.freq(:,mode_no)) 

freq_index = max(find((reshaped_proc_data.freq(:,mode_no)-chosen_freq)<=0)); 
disp(['freq index = ' , num2str(freq_index)])

switch(x_y_z)
    case(1)
mode_shape  = reshaped_proc_data.ms_x( node_plot_index ,   freq_index  ,   mode_no )  ;
    case(2)
mode_shape  = reshaped_proc_data.ms_y( node_plot_index ,   freq_index  ,   mode_no )  ;
    case(3)
mode_shape  = reshaped_proc_data.ms_z( node_plot_index ,   freq_index  ,   mode_no )  ;


end %switch(x_y_z)

max_abs_sig =  max([max(abs(reshaped_proc_data.ms_x(:,freq_index  ,   mode_no))), max(abs(reshaped_proc_data.ms_y(:,freq_index  ,   mode_no ))), max(abs(reshaped_proc_data.ms_z(:,freq_index  ,   mode_no )))]);


% reshaped_proc_data(node_number ) 

%--------------------------------------------------------------------------------------------------------------------------------------------
[model_new] = genPogoHannSignal_jim(model_old , N_cycles , chosen_freq , mode_shape , max_abs_sig , whichSig , whichFrame );

[Pogo_model_input]     =   extract_excitation(N_cycles , points_per_trace , file_name , window_type , mode_no , chosen_freq_kHz , node_plot_index, doplot );

new_model_def          = recreate_mode_shape(model_old, N_cycles , chosen_freq , reshaped_proc_data , mode_no , whichFrame, sigType);

% [Pogo_model_input]     =   extract_excitation(N_cycles , points_per_trace , file_name , window_type , mode_no , chosen_freq_kHz , node_plot_index, doplot );
%  Pogo_model_input
%--------------------------------------------------------------------------------------------------------------------------------------------

else
disp(['chosen freq is not in the range of frequencies inthe file  min = ',num2str(round(min (reshaped_proc_data.freq(:,mode_no)/1000)))  ,'kHz, max = ',num2str(round(max (reshaped_proc_data.freq(:,mode_no)/1000)) ),' kHz.'  ])
end %if  chosen_freq_kHz > min (reshaped_proc_data.freq(:,mode_no)) && chosen_freq_kHz <  max (reshaped_proc_data.freq(:,mode_no))

% plot the POGO version


aa  = new_model_def.frames{whichFrame}.sigs{whichSig}.sig    ;   %   
nt_ = new_model_def.frames{whichFrame}.ntSig                 ;   % these should not be unecessary
dt_ = new_model_def.frames{whichFrame}.dtSig                 ;   % these should not be unecessary
disp(['length of signal = ' ,num2str(length(aa))])
0

figure (2)
plot( [0:dt_ *1000000 : 1000000 * dt_*(nt_-1)] , aa , '.-')
xlabel('Time (s) ')
ylabel('Signal ')

title(['max signal = ' ,num2str(max(aa)),' min signal = ' ,num2str(min(aa))])















