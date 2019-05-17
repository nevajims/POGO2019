function [ m, pr, pTh, pz ] = genGridPipe( ri, nr, nc, nz, dr, dz )
%genGridPipe - generate a Pogo mesh for a pipe
%   [ m ] = genGridPipe( ri, nr, nc, nz, dr, dTh, dz )
%
%m - generated model struct in Pogo format
%ri - inner radius
%nr, nc, nz - number of nodes in radial, circumferential and axial directions
%dr, dz - spacing in each direction (circ calculated automatically)
%pr, pTh, pz - output node coordinates in r, theta and z directions
%
%8 noded brick elements are generated
%order of nodes and elements is radially, circumferentially, then axially
%radially, inside to out
%
%Written by Peter Huthwaite, Dec 2016
m.nDims = 3;
m.nDofPerNode = 3;

%r fastest
%theta 
%z slowest

%make nodes
dTh = 2*pi/nc;

pr = repmat( dr*( (0:nr-1) ).', [1,nc,nz])+ri;
pTh = repmat(dTh*( (1:nc).'-(nc+1)/2).', [nr,1,nz] );
pz = repmat(dz*permute( (1:nz).'-(nz+1)/2,[2,3,1]), [nr,nc,1] );

pr = reshape(pr,[],1);
pTh = reshape(pTh,[],1);
pz = reshape(pz,[],1);

px = pr.*cos(pTh);
py = pr.*sin(pTh);

m.nodePos = [px(:), py(:), pz(:)].';

nEls = (nr-1)*(nc)*(nz-1);

nNodesPerEl = 8;

elNodes = zeros(nNodesPerEl,nEls);

cubeNodes = zeros(8,1);
elCnt = 0;
for elCntZ = 1:(nz-1)
    for elCntC = 1:(nc)
        for elCntR = 1:(nr-1)
                        
            cubeNodes(1) = elCntR-1+(elCntC-1)*nr + (elCntZ-1)*nr*nc;
            cubeNodes(2) = elCntR-1+1+(elCntC-1)*nr + (elCntZ-1)*nr*nc;
            cubeNodes(3) = elCntR-1+1+(elCntC-1+1)*nr + (elCntZ-1)*nr*nc;
            cubeNodes(4) = elCntR-1+(elCntC-1+1)*nr + (elCntZ-1)*nr*nc;
            
            cubeNodes(5) = elCntR-1+(elCntC-1)*nr + (elCntZ-1+1)*nr*nc;
            cubeNodes(6) = elCntR-1+1+(elCntC-1)*nr + (elCntZ-1+1)*nr*nc;
            cubeNodes(7) = elCntR-1+1+(elCntC-1+1)*nr + (elCntZ-1+1)*nr*nc;
            cubeNodes(8) = elCntR-1+(elCntC-1+1)*nr + (elCntZ-1+1)*nr*nc;
            
            if elCntC == nc
                cubeNodes(3) = elCntR-1+1+(0)*nr + (elCntZ-1)*nr*nc;
                cubeNodes(4) = elCntR-1+(0)*nr + (elCntZ-1)*nr*nc;
                cubeNodes(7) = elCntR-1+1+(0)*nr + (elCntZ-1+1)*nr*nc;
                cubeNodes(8) = elCntR-1+(0)*nr + (elCntZ-1+1)*nr*nc;
            end

            %hex:
            elCnt = elCnt+1;
            elNodes(1:8,elCnt) = cubeNodes(:);
        end
    end
end

m.elNodes = elNodes+1;

m.elTypes = cell(1,1);
m.elTypes{1}.name = 'C3D8R';
m.elTypes{1}.paramsType = 0;

m.elTypeRefs = ones(nEls,1);

end

