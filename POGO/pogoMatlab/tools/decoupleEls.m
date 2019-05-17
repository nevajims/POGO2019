function [ mDec ] = decoupleEls( m, dupeNodes, dir )
%decoupleEls - decouple elements in a Pogo model
%   [ mDec ] = decoupleEls( m, dupeNodes, dir )
%
%m - input model structure - must contain nodePos, elNodes
%mDec - output model with decoupled elements
%dupeNodes - nodes to duplicate in order to separate elements
%dir - direction vector. Elements in this direction from the decoupled 
%node will be allocated the new node.
%
%Written by Peter Huthwaite, Dec 2016

    dupeNodes = dupeNodes(:);
    mDec = m;
    nNodesOrig = size(m.nodePos,2);
    nNodesTot = nNodesOrig+length(dupeNodes);
    mDec.nodePos = [m.nodePos m.nodePos(:,dupeNodes)];
    
    nEls = size(m.elNodes,2);
    nPerEl = size(m.elNodes,1);
    
    [ex, ey, ez] = getElCents(m);
    if m.nDims == 3
        elCents = [ex(:) ey(:) ez(:)].';
    else
        elCents = [ex(:) ey(:)].';
    end
    
    for eCnt = 1:nEls
        for nCnt = 1:nPerEl
            nRep = find(m.elNodes(nCnt,eCnt) == dupeNodes);
            if ~isempty(nRep)
                nPos = m.nodePos(:,dupeNodes(nRep));
                ePos = elCents(:,eCnt);
                if (dotprod((nPos(:)-ePos(:)).',dir(:))>0)
                    mDec.elNodes(nCnt,eCnt) = nRep+nNodesOrig;
                end
            end
        end
    end
    
    %inter = intersect(mDec.elNodes(:),1:nNodesTot);
    diff = setdiff(1:nNodesTot,mDec.elNodes(:));
    if ~isempty(diff)
        diff
        error('Some hanging nodes')
    end
end

