function [ m3 ] = combinePogoMesh( m1, m2, n1, n2 )
%combinePogoMesh - combine two Pogo meshes together
%   [ m3 ] = combinePogoMesh( m1, m2, n1, n2 );
%
%m1, m2 - input models to be combined 
%         m1 is the 'master' - keeps any sources etc., these will be
%         dropped from 2
%m3 - output model
%n1 - stitching nodes in 1
%n2 - corresponding nodes in 2 in same position
%
%Written by P. Huthwaite, Aug 2017
%Not for distribution

if m1.nDims ~= m2.nDims
    error('Dimensions of the two models do not match')
end

ln1 = length(n1);
ln2 = length(n2);
if ln1 ~= ln2
    error('length(n1) ~= length(n2) - the node arrays should be equal in size')
end

nNodesM1 = size(m1.nodePos,2);
nNodesM2 = size(m2.nodePos,2);

%derive a node mapping for m2
nodeMap2r = 1:nNodesM2;
nodeMap2r(n2) = [];
nodeMap2 = zeros(nNodesM2,1);
nodeMap2(nodeMap2r) = 1:(nNodesM2-ln2);
nodeMap2 = nodeMap2+nNodesM1;
nodeMap2(nodeMap2 == nNodesM1) = n1;

%nodeMap2

m3 = m1;
elNodesTemp = m2.elNodes;
elNodesTemp(elNodesTemp == 0) = 1; %just to avoid indexing issues
elNodesTemp = nodeMap2(elNodesTemp);
elNodesTemp(m2.elNodes == 0) = 0;

nPerEl1 = size(m1.elNodes,1);
nPerEl2 = size(m2.elNodes,1);

if nPerEl2 > nPerEl1
    %need to expand
    nExtra = nPerEl2 - nPerEl1;
    nEls1 = size(m1.elNodes,2);
    m3.elNodes = [m3.elNodes; zeros(nExtra,nEls1)];
end
m3.elNodes = [m3.elNodes elNodesTemp];

nEt1 = length(m1.elTypes);
nEt2 = length(m2.elTypes);

for etCnt = 1:nEt2
    m3.elTypes{nEt1+etCnt} = m2.elTypes{etCnt};
end
elTypesTemp = m2.elTypeRefs+nEt1;
m3.elTypeRefs = [m1.elTypeRefs; elTypesTemp];

nodePosTemp = m2.nodePos;
nodePosTemp(:,n2) = [];
m3.nodePos = [m3.nodePos nodePosTemp];

if isfield(m2,'frames') || isfield(m2,'measSets') || isfield(m2,'nodeFix')
    warning('Non geometry related features found in m2. Ignoring.')
end

end

