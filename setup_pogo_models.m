function POGO_model =   setup_pogo_models(batch_types , freqs , modes) 


%batch_types
%(1) multiple models with    1 freq and mode
%(2) 1 model  with  multiple
%(3)


P_W_D = pwd  ;
cd('m:\SAFE_solutions')
cd(P_W_D)

current_folder_info =  dir('*')                                                                                                  ;
isdirs ={current_folder_info.isdir}                                                                                              ;
dirnames ={current_folder_info.name}                                                                                             ;
dummy = cell2mat(isdirs)                                                                                                         ;
dir_inds = find(dummy == 1)                                                                                                      ;
dummy2 = dirnames(dir_inds)                                                                                                      ;
all_dir_names = {dummy2{3:end}}                                                                                                  ;
choice = listdlg('PromptString' , 'Select mesh group to solve' , 'SelectionMode' , 'single'  , 'ListString' , all_dir_names)     ;
cd(all_dir_names{choice})


valid = 1 ; 
reason_not_valid= '' ;

switch (batch_types)

    case(1)
% multiple files but single others
if length(freqs) == 1 && length(modes) ==1

    
else
    
    
end %if length(freqs) == 1 && length(modes) ==1

    
    case(2)        

        
    case(3)

        
end %switch (batch_types)


%-----------------------------------------------------
% POGO_model =   setup_pogo_models(                        SAFE_SOLUTION ,   MODE_NO,  chosen_freq,  nCyc)
% SPECIFICALLY FOR RECREATING SAFE MODELS


% either choose

% multiple safe solutions
%
%
 
% materials types
% (1) -  steel
% (2) -  copper
%  for a specific SAFE solution ::
%  Use the SAFE mesh and extrude it
%  Use the modal solution to the SAFE model 
%  select a mode number /  frequency /  number of cycles  and create an  excitation on the end nodes
%  Input the other input parameters for the solutins::     
%  choose the format of the solutions blocks (needs to be limited)
% TO DO-------------------------------
% Setup 
% automatically create a NOTES file
% save to the NAS directory  
% put in a user defined name rather than  'model'-     use the SAFE model name plus the mode and frequency
% plus    
% POST PROCESSING --     this needs to be checked  measuresets
% TO DO-------------------------------

% make the save directory an input parameter
% MODE_NO = 4 ;  %   mode 4 is longitudinal in safe model
% load('rpd_8.mat')
% chosen_freq = 100e3   in Hz
% number of cycles in toneburst
% nCyc = 6                      ;
% put the 3 material properties in here 2

%-------------------------------------------------------------------------------------
% Variables -----------------------------------------------------
%-------------------------------------------------------------------------------------
MODEL_LENGTH  = 2 ;   % METRES  
nCyc          = 6 ;   % pulse length
%-------------------------------------------------------------------------------------

%-------------------------------------------------------------------------------------
% MESH -----------------------------------------------------
%-------------------------------------------------------------------------------------
SAFE_MESH = SAFE_SOLUTION.mesh;
mesh_details = get_element_details(SAFE_MESH);
dx = mesh_details.mean_edge_length;
dz = dx ;
nz = ceil(MODEL_LENGTH/dx);    %make it a 1m length
POGO_model = convert_SAFE_POGO(SAFE_MESH , dz , nz , 1 );
nEls = size(POGO_model.elNodes,2);
%-------------------------------------------------------------------------------------

%-------------------------------------------------------------------------------------
%  Input the material properties ---------------------------------     
%-------------------------------------------------------------------------------------
material_index = 2; %   2 is copper
[POGO_model, ~ , speed_cSh ]  = get_material_properties(material_index);
%------------------------------------------------------------------------------------

%-------------------------------------------------------------------------------------
%  Time length of the model ---------------------------------     
%-------------------------------------------------------------------------------------
CFL = 0.3;     % courant = 0.3; % mix of sizes so be conservative with stability   
time_there_and_back = MODEL_LENGTH / speed_cSh ;
POGO_model.dt = (dx/speed_cSh) * CFL;
POGO_model.nt = round(time_there_and_back/POGO_model.dt);
%------------------------------------------------------------------------------------

%-------------------------------------------------------------------------------------
%  Input the element properties
%-------------------------------------------------------------------------------------
POGO_model.nDims = 3 ;
POGO_model.elTypes = { struct('name','C3D6R', 'paramsType', 0) };
POGO_model.elTypeRefs = ones(nEls,1);
POGO_model.nDofPerNode = 3;
POGO_model.measFreq = 2;
POGO_model.measStart = 1;
%-------------------------------------------------------------------------------------

%-------------------------------------------------------------------------------------
% Input the excitation
%-------------------------------------------------------------------------------------
whichFrame = 1;
whichSig = 1;  
POGO_model = recreate_mode_shape(POGO_model , nCyc , chosen_freq , SAFE_SOLUTION , MODE_NO , whichFrame ,whichSig);
%-------------------------------------------------------------------------------------

%-------------------------------------------------------------------------------------
% Specific  Post processing data sets 
%-------------------------------------------------------------------------------------
[POGO_model.measSets,POGO_model.fieldStoreIncs] =   get_post_processing_sets(output_index);
%-------------------------------------------------------------------------------------

%  save_file_name  -   should incorparate the   frequency , mode ,

file_and_path_name = ['\\nas\POGO_solutions\', save_file_name ];

savePogoInp([file_and_path_name,'.pogo-inp' ],POGO_model);
save       (file_and_path_name ,'POGO_model')





%-------------------------------------------------------------------------------------
%savePogoInp('model.pogo-inp', POGO_model);
end %function POGO_model =   setup_POGO_model(SAFE_SOLUTION,MODE_NO )
%------------------------------------------------------------------------------------------------------------------------------------------------
%------------------------------------------------------------------------------------------------------------------------------------------------

function  POGO_mesh = convert_SAFE_POGO(SAFE_MESH , del_z , no_z , doplot )  

%  SAFE_MESH 
%  del_z  -  z spacing
%  no_z   - number of node slices
%  doplot -   plot 2d     face
%  to display the mesh-  ignor the elements onthe side and create square
%  ones -  create triangular ones on the front from the safe mesh and then merge the two together 
%  create a new stucture
%  get it working with the shaped meshes



if no_z > 1
SAFE_NODES      =   SAFE_MESH.nd.pos        ;
SAFE_EL_NODES   =   SAFE_MESH.el.nds        ; 
nodePos = [SAFE_NODES'  ;  zeros(1,length(SAFE_NODES))] ; 

for index = 1 : no_z-1
new_layer = [SAFE_NODES';  (del_z*index)*ones(1,length(SAFE_NODES))];
nodePos = [nodePos,new_layer];  
end %for index = 1 : no_z-1

POGO_mesh.nodePos =  nodePos;

%now to the element nodes
nodes_per_layer             = size(SAFE_NODES   ,1)     ; 
element_nodes_per_layer     = size(SAFE_EL_NODES,1)     ;
nodes_per_element           = size(SAFE_EL_NODES,2)*2   ;
elements_nodes =  zeros(nodes_per_layer*(no_z - 1) , nodes_per_element);
for index = 1: no_z-1
elements_nodes( (index-1)*element_nodes_per_layer+1:index*element_nodes_per_layer, :)    = [SAFE_EL_NODES + ((index-1)*nodes_per_layer)    , SAFE_EL_NODES + (index*nodes_per_layer)];
end %for index = 1: no_z-1
POGO_mesh.elNodes = elements_nodes';
SAFE_MESH = get_ordered_edge_nodes(SAFE_MESH,0); 
% POGO_display_element_nodes
%(1)  Face =  (safe nodes with a NAN on the end)
PDEN_1 = [SAFE_MESH.el.nds,NaN(size(SAFE_MESH.el.nds,1),1)];
%size = 394     4
%(2)  edge
% PDEN_2
display_el_no = 0;
for index = 1 : no_z-1
for index_2 = 1 : length (SAFE_MESH.ordered_edge_nodes)
display_el_no = display_el_no +1 ;

if index_2 ~= length (SAFE_MESH.ordered_edge_nodes)
PDEN_2(display_el_no,:) = [(SAFE_MESH.ordered_edge_nodes(index_2) +(length( SAFE_MESH.nd.pos)*(index-1)  )), (SAFE_MESH.ordered_edge_nodes(index_2+1) +(length( SAFE_MESH.nd.pos)*(index-1))),...
    (SAFE_MESH.ordered_edge_nodes(index_2+1) +(length( SAFE_MESH.nd.pos)*(index))), (SAFE_MESH.ordered_edge_nodes(index_2)+(length( SAFE_MESH.nd.pos)*(index)))];
else
PDEN_2(display_el_no,:) = [(SAFE_MESH.ordered_edge_nodes(index_2) +(length( SAFE_MESH.nd.pos)*(index-1)  )), (SAFE_MESH.ordered_edge_nodes(1) +(length( SAFE_MESH.nd.pos)*(index-1))),...
(SAFE_MESH.ordered_edge_nodes(1) +(length( SAFE_MESH.nd.pos)*(index))), (SAFE_MESH.ordered_edge_nodes(index_2)+(length( SAFE_MESH.nd.pos)*(index)))];
  
end %if index_2 ~= length (SAFE_MESH.ordered_edge_nodes)
end %for index_2 = 1 : length (SAFE_MESH.ordered_edge_nodes)
end %for index = 1 : no_z-1

PDEN_3 = [SAFE_MESH.el.nds+(length( SAFE_MESH.nd.pos)*(no_z-1)),NaN(size(SAFE_MESH.el.nds,1),1)];

POGO_mesh.POGO_display_element_nodes = [PDEN_1;PDEN_2;PDEN_3]; 


% need number of 


if doplot == 1
figure;
fv.Vertices =  POGO_mesh.nodePos'   ;
fv.Faces    =  POGO_mesh.POGO_display_element_nodes;
patch(fv, 'FaceColor', 'c');
%axis equal;
%axis off;
view([-135 35]);
%axis equal
end %if doplot == 1


else
disp('this is a single layer of nodes')
POGO_mesh = 'void'  ;
end


end %function  POGO_mesh = convert_SAFE_POGO(SAFE_MESH , del_z , no_z ,doplot)  

function mesh_ = get_ordered_edge_nodes(mesh_,do_plot) 

for index = 1 : length (mesh_.nd.pos)
lens(index) =  length(find(mesh_.el.nds == index));
end %for index = 1 : length (mesh_.nd.pos)

edge_node_nos = find (lens <= 4);

for index = 1:length(edge_node_nos)
[elements_temp,~] = find(mesh_.el.nds == edge_node_nos(index));

if length(elements_temp)==2
elements_temp = [elements_temp;NaN;NaN];
end %if length(elements_temp)==2

if length(elements_temp)==3
elements_temp = [elements_temp;NaN];
end %if length(elements_temp)==2

elements_on_node(index,:) = elements_temp'  ;

end %for index = 1:length(mesh_.el.nds)

[~,f_index] = min(mesh_.nd.pos(edge_node_nos))      ;
[~,node_node_index] = min(mesh_.nd.pos(edge_node_nos,1));
first_node = edge_node_nos(node_node_index)                 ;

% find the elements that have a line on the edge
unique_elements_on_edge = unique(elements_on_node(~isnan(elements_on_node)));

% now how to find the first element with lin on the edge
for index = 1 : length(unique_elements_on_edge)
no_nodes_per_element(index) =  length(find( elements_on_node == unique_elements_on_edge(index)));
end %for index = 1 : length(elements_on_edge)

elements_with_line_on_edge  = unique_elements_on_edge(find(no_nodes_per_element ==2))  ;
elements_on_first_node =  elements_on_node(first_node,find(~isnan(elements_on_node(first_node,:))==1));

if length(elements_on_first_node)==3

count_ = 0;    
for index = 1:3
if length(find(elements_with_line_on_edge == elements_on_first_node(index)))
count_ = count_+1;
line_elements_on_first_node(count_) = elements_on_first_node(index);
end %if length(find(elements_with_line_on_edge == elements_on_first_node(index)))
end %for index = 1:3
elseif length(elements_on_first_node) == 2
line_elements_on_first_node = elements_on_first_node;
else
disp('shouldnt get here')    
end

% now find the two connected nodes

[node_nos_1,~] = find(elements_on_node == line_elements_on_first_node(1));
[node_nos_2,~] = find(elements_on_node == line_elements_on_first_node(2));
node_list  = [node_nos_1;node_nos_2];
connected_nodes =  node_list(find(node_list~=first_node ));

[~,second_node_index] = max(mesh_.nd.pos(connected_nodes,2));

second_node = connected_nodes(second_node_index);

ordered_node_list = [first_node;second_node];  
% elements_with_line_on_edge


for index = 3 : length(edge_node_nos)
% find elements inprevios node
prev_element_possibilities =  elements_on_node(ordered_node_list(index-2),:) ;
prev_element_possibilities =  prev_element_possibilities(find(~isnan(prev_element_possibilities)));

%keyboard
next_element_possibilities =  elements_on_node(ordered_node_list(index-1),:) ;
next_element_possibilities =  next_element_possibilities(find(~isnan(next_element_possibilities)));


for index_2 = 1:length(next_element_possibilities)
% must be an edge line element and not a prev element
 
if  ~sum(prev_element_possibilities == next_element_possibilities(index_2))==1  && length (find(elements_with_line_on_edge == next_element_possibilities(index_2))) == 1
[new_node_pos,~] = find(elements_on_node== next_element_possibilities(index_2));    
ordered_node_list(index) =  new_node_pos(find(new_node_pos~=ordered_node_list(index-1)));

end % if  ~sum(prev_element_possibilities == next_element_possibilities(index_2))==1  && length (find(elements_with_line_on_edge == next_element_possibilities(index_2))) == 1

end %for index_2 = 1:length(next_element_possibilities)
end % for index = 3 : length(edge_node_nos)

mesh_.ordered_edge_nodes = edge_node_nos(ordered_node_list)                            ;

if do_plot ==1
% added in the fist node at the end so that the line goes all the way round

plot(  mesh_.nd.pos([mesh_.ordered_edge_nodes,mesh_.ordered_edge_nodes(1)],1)   ,mesh_.nd.pos([mesh_.ordered_edge_nodes,mesh_.ordered_edge_nodes(1)],2)  , 'x-')

hold on
plot(mesh_.nd.pos(first_node,1),mesh_.nd.pos(first_node,2)  ,'rx', 'MarkerSize',15)

axis equal
end % if do_plot ==1

end %function ordered_edge_node = get_mesh_.ordered_edge_nodes(mesh_,do_plot) 

function [POGO_model,  c0, cSh]     =    get_material_properties(material_index)

switch(material_index)

    case(1)  % Steel
               % STEEL PROPERTIES
E = 210e9;
nu = 0.32;
rho = 8050;   % kg/m3


    case(2)  % Copper
               % COPPER PROPERTIES
E   =        117e9        ;       % youngs_modulus 
nu  =        0.35         ;       % poissons_ratio
rho =        8960         ;       % density kg/m3

end %switch(material_index)

G   =        E/(2*(1+nu)) ; % shear modulus

%sound speeds
c0  = sqrt(E*(1-nu)/(rho*(1+nu)*(1-2*nu))) ; % m/s

% disp(num2str(c0))
cSh = sqrt(E/(2*rho*(1+nu)))               ; % m/s   

POGO_model.matTypeRefs = ones(nEls,1);     %then define material 1 with the properties from before
POGO_model.matTypes{1}.paramsType = 0;
POGO_model.matTypes{1}.paramValues = [E, nu, rho];



if show_vals  == 1

        
    
end %if show_vals  == 1


end

function new_model_def              = recreate_mode_shape(model_def, N_cycles ,chosen_freq,reshaped_proc_data,mode_no,whichFrame,sigType)

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

function [measSets,fieldStoreIncs]  =   get_post_processing_sets(output_index)

% rPoint1 = 1; rPoint2  = floor(size(POGO_model.nodePos,2)/3) ; rPoint3  = floor(size(POGO_model.nodePos,2)/2) ; rPoint4  = floor(2*size(POGO_model.nodePos,2)/3) ; rPoint5  = size(POGO_model.nodePos,2);
% POGO_model.measSets{1}.measNodes = [rPoint1, rPoint1, rPoint1, rPoint2, rPoint2, rPoint2];
% just get al thge nodes in the middle plane of nodes
switch (output_index)
    case(1)    
% this needs to be checked  below
 [~,mid_index] =   min(            abs(     POGO_model.nodePos(3,:)- max(POGO_model.nodePos(3,:))       /2))      ;
% this needs to be checked  below

mid_value = POGO_model.nodePos(3,mid_index);
mid_node_vals     =  find(POGO_model.nodePos(3,:) == mid_value);
front_node_vals   =  find(POGO_model.nodePos(3,:) == 0);
end_node_vals     =  find(POGO_model.nodePos(3,:) == max(POGO_model.nodePos(3,:)));
measset_vals      =  [front_node_vals,mid_node_vals,end_node_vals];

measSets{1}.name = 'main';
measSets{1}.measNodes = sort([measset_vals,measset_vals,measset_vals]);  
measSets{1}.measDof = repmat([1,2,3],[1,length(measset_vals)]);
fieldStoreIncs = 1:25:100*floor (POGO_model.nt/100);  

end %switch (output_index)
end %function [measSets,fieldStoreIncs] =   get_post_processing_sets(   )
