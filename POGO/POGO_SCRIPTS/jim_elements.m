JIM_DATA_PATH = 'C:\Users\Remote\Desktop\JIM\meshes\MESH_RAIL_56_1.mat';

load(JIM_DATA_PATH);

fv.Vertices = mesh.nd.pos;
fv.Faces = mesh.el.nds;
patch(fv, 'FaceColor', 'c');
axis equal;
axis off;