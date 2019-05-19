function new_model_def = recreate_mode_shape(model_def, N_cycles ,chosen_freq,reshaped_proc_data,mode_no,whichFrame,sigType)

% Confirm that node - numbers dont change between the safe model and the pogo model face
% model_def         -  the defined model structure for the pogo analysis
% N_cycles          - numbero of cycles for the excitation pulse
% chosen_freq       -
% reshaped_proc_data
% mode_no
% whichFrame
% sigType  0 =forc
% inputs       --------------
% outputs      --------------
% look with x/y and Z
% loop with all the modes in


freq_index = max(find((reshaped_proc_data.freq(:,mode_no)-chosen_freq)<=0)); 
disp(['freq index = ' , num2str(freq_index)])
max_abs_sig =  max([max(abs(reshaped_proc_data.ms_x(:,freq_index , mode_no))), max(abs(reshaped_proc_data.ms_y(:,freq_index  ,   mode_no ))), max(abs(reshaped_proc_data.ms_z(:,freq_index  ,   mode_no )))]);
whichSig  = 0;


for index_1 = 1: size(reshaped_proc_data.ms_x,1) % increment for the number of nodes in the mesh
for index_2 = 1:3  %  increment for x_y nad z     
    whichSig = whichSig + 1;
    
switch(index_2)
    case(1)
mode_shape  = reshaped_proc_data.ms_x(index_1, freq_index , mode_no);
    case(2)
mode_shape  = reshaped_proc_data.ms_y(index_1, freq_index , mode_no);
    case(3)
mode_shape  = reshaped_proc_data.ms_z(index_1, freq_index , mode_no);
end %switch(x_y_z)

[model_def] = genPogoHannSignal_jim(model_def , N_cycles , chosen_freq , mode_shape , max_abs_sig , whichSig , whichFrame );


% say what nodes to apply it to
model_def.frames{whichFrame}.sigs{whichSig}.nodeSpec = index_1 ;

% and which degrees of freedom

switch(index_2)
   case(1)
model_def.frames {whichFrame}.sigs {whichSig}.dofSpec = [1];
   case(2)
model_def.frames {whichFrame}.sigs {whichSig}.dofSpec = [2];
   case(3)
model_def.frames {whichFrame}.sigs {whichSig}.dofSpec = [3];
end %switch(index_2)

% set the signal amplitudes in each direction
model_def.frames{whichFrame}.sigs {whichSig}.sigAmps = [1];  %   dont change the amplitude as it has already been set by the (normalised) abs value of the mode shape at that node
% set the type to 0, for force
model_def.frames{whichFrame}.sigs {whichSig}.sigType = sigType;   %  0 for force 1 for displacement

end %for index_2 = 1:3  %  increment for x_y nad z     
end  %for index_1 = 1: size(reshaped_proc_data.ms_x,1)


new_model_def = model_def;


%  need a LOOP FOR x_y_z
%  NEED A SECOND LOOP FOR THE NODES
%  A value increments for both of these







end %function new_model = recreate_mode_shape(   )
