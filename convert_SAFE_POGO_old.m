function  POGO_mesh = convert_SAFE_POGO_old(SAFE_MESH , del_z , no_z , doplot )  


%  SAFE_MESH 
%  The latest version of this is now embedded into the function
%  setup_pogo_model
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
