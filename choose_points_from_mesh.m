

function  Choose_points_from_mesh(POGO_MESH)

% Choose from the pogo mesh  / or from the results data file
% Create a sequential program to start and then convert to a windows program 


Selection_Options = strsplit(num2str(sort(unique(POGO_MESH.nodePos(3,:)))));
[indx,tf] = listdlg('PromptString','Select a file:','SelectionMode','single','ListString', Selection_Options);

if tf == 1
temp_= sort(unique(POGO_MESH.nodePos(3,:)));
chosen_z_value =  temp_(indx);
disp(['Chosen Z value of slice = ' , num2str(chosen_z_value),' mm.'])
end %if tf == 1

[~,Node_Indices] = find(POGO_MESH.nodePos(3,:) == chosen_z_value);


element_stack = [];

for index = 1:length(Node_Indices)
[~,temp]= find(POGO_MESH.elNodes(1:3,:)==Node_Indices(index));
element_stack = [element_stack;temp];
end %for index = 1:length(Node_Indices)
slice_elements = sort(unique(element_stack));

fv.Vertices =  POGO_MESH.nodePos'                           ;
fv.Faces    =  POGO_MESH.elNodes(1:3,slice_elements)'       ;

patch(fv, 'FaceColor', 'c')                   ;
axis equal


[chosen_x , chosen_y] = ginput(1)             ;


% find the closest node to this point


keyboard


% Now get all the nodes in that slice and find the elements
% First select the slice you want

end % function  Choose_points_from_mesh(POGO_MESH)



% display the mesh
% select points from it
% select the slice
% give the node position and 
% the results
% plot the x y and x displacement for that node
% give the phase / amplitude and plot the fft for a particular region of
% the results
% auto   - first / second
% manual -  

%--------


