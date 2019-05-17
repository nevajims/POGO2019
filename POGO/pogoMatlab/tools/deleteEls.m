function [ mDel ] = deleteEls( m, elsDelete )
%deleteEls - delete elements from a Pogo model
%   [ mDel ] = deleteEls( m, elsDelete )
%
% m - original model
% elsDelete - list of elements to be deleted
% mDel - model with elements deleted
%Written by P. Huthwaite, 2017
%Not to be distributed.

    checkNodes = m.elNodes(:,elsDelete);
    checkNodes = unique(checkNodes(:));
    checkNodes(checkNodes == 0) = [];
    
    nElNodes = size(m.elNodes,2);
    elRemap = 1:nElNodes;
    elRemap = setdiff(elRemap,elsDelete);
    newElNodes = m.elNodes(:,elRemap);
    newElNodes = unique(newElNodes);
    newElNodes(newElNodes == 0) = [];
    
    %is there anything in checknodes which isn't in the final list?
    nodesDelete = setdiff(checkNodes,newElNodes);
    
    nNodesOrig = size(m.nodePos,2);
    nodesRemap = 1:nNodesOrig;
    nodesRemap = setdiff(nodesRemap,nodesDelete);
    
    nodesRemapInv = zeros(nNodesOrig,1);
    nodesRemapInv(nodesRemap) = 1:length(nodesRemap);
       
    mDel = m;
    mDel.nodePos = m.nodePos(:,nodesRemap);
    n = m.elNodes(:,elRemap);
    %will have some 0 values if have mixed elements
    blanks = (n == 0);
    n(blanks) = 1;
    mDel.elNodes = nodesRemapInv(n);
    mDel.elNodes(blanks) = 0;
    mDel.matTypeRefs = m.matTypeRefs(elRemap);
    mDel.elTypeRefs = m.elTypeRefs(elRemap);
    if isfield(m,'orientRefs')
        mDel.orientRefs = m.orientRefs(elRemap);
    end
    if isfield(m,'fixNodes')
        mDel.fixNodes = nodesRemapInv(m.fixNodes);
        mDel.fixDof(mDel.fixNodes == 0) = [];
        mDel.fixNodes(mDel.fixNodes == 0) = [];
    end
    for sCnt = 1:length(m.measSets)
        mDel.measSets{sCnt}.measNodes = nodesRemapInv(m.measSets{sCnt}.measNodes);
        mDel.measSets{sCnt}.measDof(mDel.measSets{sCnt}.measNodes == 0) = [];
        mDel.measSets{sCnt}.measNodes(mDel.measSets{sCnt}.measNodes == 0) = [];
    end
    if isfield(m,'dofGroups')
        for n = 1:length(m.dofGroups)
            mDel.dofGroups{n}.nodeSpec = nodesRemapInv(m.dofGroups{n}.nodeSpec);
            mDel.dofGroups{n}.dofSpec(mDel.dofGroups{n}.nodeSpec == 0) = [];
            mDel.dofGroups{n}.nodeSpec(mDel.dofGroups{n}.nodeSpec == 0) = [];
        end
    end
    
    for fCnt = 1:length(m.frames)
        for sCnt = 1:length(m.frames{fCnt}.sigs)
            if ~isfield(mDel.frames{fCnt}.sigs{sCnt},'isDofGroup') || mDel.frames{fCnt}.sigs{sCnt}.isDofGroup == 0
                n = nodesRemapInv(mDel.frames{fCnt}.sigs{sCnt}.nodeSpec);
            
                %delete references to the deleted nodes
                mDel.frames{fCnt}.sigs{sCnt}.dofSpec(n == 0) = [];
                mDel.frames{fCnt}.sigs{sCnt}.sigAmps(n == 0) = [];

                n(n == 0) = [];
                mDel.frames{fCnt}.sigs{sCnt}.nodeSpec = n;
            end
        end
    end       
    
    %remove any unused materials
    nMats = length(m.matTypes);
    keepMats = unique(mDel.matTypeRefs);
    nMatsNew = length(keepMats);
    
    newMat = zeros(nMats,1);    %maps old materials to new ones
    newMat(keepMats) = 1:length(keepMats);
    mDel.matTypeRefs = newMat(mDel.matTypeRefs);
    mDel = rmfield(mDel,'matTypes');
    for mCnt = 1:nMatsNew
        mDel.matTypes{mCnt} = m.matTypes{keepMats(mCnt)};
    end
end

