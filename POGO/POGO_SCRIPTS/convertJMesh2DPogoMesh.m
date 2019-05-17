function n = convertJMesh2DPogoMesh(jmesh)

mesh = jmesh.mesh;

n.nodePos = mesh.nd.pos';
n.elNodes = mesh.el.nds';
%CPE3 = plane strain triangle
n.elTypes = { struct('name','CPE3', 'paramsType', 0) };
n.nDims = 2;
n.nDofPerNode = 2;
%material type 1 (user defined elsewhere)
n.elTypeRefs = ones(1,size(mesh.el.nds,1 ))';

%compare with
%POGO_PATH = 'C:\Users\Remote\Desktop\POGO_SCRIPTS\';
%m = loadVtuMesh(strcat(POGO_PATH, 'temp.vtu') );
%{
  struct with fields:

        nodePos: [2×75046 double]
        elNodes: [3×148262 double]
        elTypes: {[1×1 struct]}
          nDims: 2
    nDofPerNode: 2
     elTypeRefs: [148262×1 double]
%}

end