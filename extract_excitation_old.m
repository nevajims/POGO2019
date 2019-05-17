% ****DONE  Load the RPD file
% ****DONE  Choose the mode
% ****DONE  how to go from a complex number to a movement
% ****DONE  go from 0 to 2*pi
% ****DONE  get the starting angle from the angle of the complex number
% ****DONE  create a hann window of the sme size 
% ****DONE  produce a file
% ****DONE  put the two in the same plot
% ****DONE


% Pogo_model_input
% put amps as a % of all the nodes
% create a mesh animator
% create an app from
% old



% N_cycles , del_t,  sim_t  

 function  [Pogo_model_input]     =   extract_excitation_old(N_cycles , points_per_trace , file_name , window_type , mode_no , chosen_freq_kHz , node_plot_index, doplot )

%function  [Pogo_model_input]     =   extract_excitation(N_cycles , del_t,  sim_t  , file_name , window_type , mode_no , excit_freq_kHz , node_plot_index, doplot )



%  window_type  = 1 (hann window)
%  window_type  = anything else (no window)

% [Pogo_model_input]     =   extract_excitation(6 , 4.6499e-08,  6.6667e-4  , 'rpd_3.mat' , 1 , 4 , 100 ,38, 1 )

% excit_freq_kHz = 100;
% del_t           =  4.6499e-08 ;
% sim_t           =  6.6667e-4;  
% Create the signal the add on the zero end to it
% del_t  =   the time step
% must be a minumum of two samples per cycle
% sim_t  =   the simulation time -  must be long enough for the entire
% excitation signal to fit in

% need the delta t and the and a 
% [Pogo_model_input]     =   extract_excitation(6 , 20 , 'rpd_16.mat' , 1 , 4 , 100000 ,1)

chosen_freq_kHz                                  = chosen_freq_kHz*1000                       ; 
%pulse_time                                  = 1/excit_freq * N_cycles                   ; 
%sample_freq                                 = 1/del_t                                   ;

%if pulse_time < sim_t
%if sample_freq > 2 * excit_freq 
   
%points_per_trace                            =  excit_freq / sample_freq                 

file_name_removed                            = file_name                                ; 
file_name_removed(find( file_name =='_'))    = ' '                                      ;

freq_range                                   = 0:2*pi/points_per_trace:N_cycles*2*pi    ;
hann_w                                       = (1-cos(freq_range/N_cycles))/2           ;
none_                                        = ones (size(freq_range))                  ;

switch(window_type)
    case(1)
        wind      = hann_w;
        wind_name = 'Hann';   

    case(2)
        wind       = none_;
        wind_name  = 'None';

    otherwise
        wind       = none_;
        wind_name  = 'None';

end %switch(window_type)


del_t = 1000000/chosen_freq_kHz ;   % us
load(file_name)

if  chosen_freq_kHz > min (reshaped_proc_data.freq(:,mode_no)) && chosen_freq_kHz <  max (reshaped_proc_data.freq(:,mode_no)) 
freq_index = max(find((reshaped_proc_data.freq(:,4)-chosen_freq_kHz)<=0)); 

time_ = del_t * freq_range/(2*pi);  %   in us

for node_index = 1:size(reshaped_proc_data.ms_x,1)

ms_x = reshaped_proc_data.ms_x( node_index ,   freq_index  ,   mode_no );
amp_x = abs(ms_x);
phase_x = angle(ms_x);
sig_x(:,node_index) = amp_x*sin(freq_range + phase_x).*wind;
 
ms_y = reshaped_proc_data.ms_y(node_index ,   freq_index  ,   mode_no );
amp_y = abs(ms_y);
phase_y = angle(ms_y);
sig_y(:,node_index) = amp_y*sin(freq_range + phase_y).*wind;

ms_z = reshaped_proc_data.ms_z( node_index ,   freq_index  ,   mode_no );
amp_z = abs(ms_z);
phase_z = angle(ms_z);
sig_z(:,node_index) = amp_z*sin(freq_range + phase_z).*wind;

end %for node_index = 1:size(reshaped_proc_data.ms_x,1)

Pogo_model_input.mesh                   = reshaped_proc_data.mesh          ;
Pogo_model_input.sig_x                  = sig_x                            ;
Pogo_model_input.sig_y                  = sig_y                            ;
Pogo_model_input.sig_z                  = sig_z                            ;
Pogo_model_input.time_                  = time_                            ; 


Pogo_model_input.stats.N_cycles         = N_cycles                         ;
Pogo_model_input.stats.window_type      = window_type                      ;
Pogo_model_input.stats.file_name        = file_name                        ;
Pogo_model_input.stats.points_per_trace = points_per_trace                 ;
Pogo_model_input.stats.number_of_nodes  = size(reshaped_proc_data.ms_x,1)  ;
Pogo_model_input.stats.wind_name        = wind_name                        ;
Pogo_model_input.stats.mode_no          = mode_no                          ;
Pogo_model_input.stats.excit_freq_MHz  = chosen_freq_kHz/1000                 ;

if doplot ==1

    
if node_plot_index <= size(reshaped_proc_data.ms_x,1)
  
figure
max_sig = max(max([sig_x(:,node_plot_index),sig_y(:,node_plot_index),sig_z(:,node_plot_index)]));
min_sig = min(min([sig_x(:,node_plot_index),sig_y(:,node_plot_index),sig_z(:,node_plot_index)]));


max_all_sig = max(max([sig_x,sig_y,sig_z]));


subplot(3,2,1)
plot(time_,100*sig_x(:,node_plot_index)/max_sig ,'.')
ylim([-100 100])
ylabel('ampl (%)')
set(gca,'xtick',[])
title(['X: A = ', num2str(round( 100*max(sig_x(:,node_plot_index))/max_sig)),' % (',num2str(round( 100*max(sig_x(:,node_plot_index))/max_all_sig)),'%), ,P = ',num2str(round(180*phase_x/pi)),' deg' ])

subplot(3,2,3)
plot(time_,100*sig_y(:,node_plot_index)/max_sig ,'.')
ylim([-100 100])
ylabel('ampl (%)')
set(gca,'xtick',[])
title(['Y: A = ', num2str(round( 100*max(sig_y(:,node_plot_index))/max_sig)),' %,(',num2str(round( 100*max(sig_y(:,node_plot_index))/max_all_sig)),'%), P = ',num2str(round(180*phase_y/pi)),' deg' ])

subplot(3,2,5)
plot(time_,100*sig_z(:,node_plot_index)/max_sig ,'.')
ylim([-100 100])
xlabel ('time (us)')
ylabel('ampl (%)')
title(['Z: A = ', num2str(round( 100*max(sig_z(:,node_plot_index))/max_sig)  ),' %,(',num2str(round( 100*max(sig_z(:,node_plot_index))/max_all_sig)),'%), P = ',num2str(round(180*phase_z/pi)),' deg' ])

suptitle(['(Mode = ',num2str(mode_no),', Freq = ' ,num2str(chosen_freq_kHz/1000),'kHz) Node = ',num2str(node_plot_index),' /',num2str(size(reshaped_proc_data.ms_x,1)),', win = ',wind_name   ,', (',file_name_removed ,').'])
disp(['Mode = ',num2str(mode_no),', Freq = ' ,num2str(chosen_freq_kHz/1000), ' MHz (',num2str(N_cycles),' cycles ',wind_name,').'])

subplot(3,2,[2,4,6])
fv.Vertices = Pogo_model_input.mesh.nd.pos;
fv.Faces = Pogo_model_input.mesh.el.nds;
patch(fv, 'FaceColor', 'w');
%patch(fv);
axis equal;
axis off;

hold on
plot (Pogo_model_input.mesh.nd.pos(node_plot_index,1),Pogo_model_input.mesh.nd.pos(node_plot_index,2),'o','markersize',20,'LineWidth',3, 'color','r')
else
disp(['chosen node to plot does not exist ( this mesh has ',num2str(size(reshaped_proc_data.ms_x,1)),' nodes).'])   
end % if node_plot_index <= size(reshaped_proc_data.ms_x,1)
end % if doplot ==1

else
disp(['chosen freq is not in the range of frequencies inthe file  min = ',num2str(round(min (reshaped_proc_data.freq(:,mode_no)/1000)))  ,'kHz, max = ',num2str(round(max (reshaped_proc_data.freq(:,mode_no)/1000)) ),' kHz.'  ])
Pogo_model_input = 'void';
end %if  chosen_freq_kHz > min (reshaped_proc_data.freq(:,mode_no)) && chosen_freq_kHz <  max (reshaped_proc_data.freq(:,mode_no))

%%else
%disp ('sample_freq too low')    
%end %if sample_freq > 2 * excit_freq 

%else
%disp('pulse_time > sim_t')
%end % if pulse_time > sim_t

end % function  [   ]     =   extract_excitation(  )




