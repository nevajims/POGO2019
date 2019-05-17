% convert SAFE Pogo 3d mesh into SFE surface mesh
function  SAFE_mesh = convert_POGO_SAFE(POGO_MESH, doplot)  

% make it for any type of element


SAFE_mesh.nd.pos = POGO_MESH.nodePos (1:2,find(POGO_MESH.nodePos(3,:) == POGO_MESH.nodePos(3,1)))'                                                                            ;
SAFE_mesh.el.nds   = POGO_MESH.elNodes (1:size(POGO_MESH.elNodes,1)/2, find (POGO_MESH.elNodes(1,:)   <=   max(find(POGO_MESH.nodePos(3,:) == POGO_MESH.nodePos(3,1)))))'         ;


if doplot == 1
figure;
fv.Vertices =  SAFE_mesh.nd.pos   ;
fv.Faces    =  SAFE_mesh.el.nds   ;
patch(fv, 'FaceColor', 'c');
axis equal;
axis off;
end %if doplot == 1

end %function  SAFE_mesh = convert_POGO_SAFE(POGO_MESH, doplot)  

