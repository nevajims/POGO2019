function [ field, fileVer, header ] = loadPogoField( fileName )
%loadPogoField - load field data from Pogo FE
%
% [ field, fileVer, header ] = loadPogoField( fileName )
%
%fileName - the file name
%fileVer - the file version
%header - the header for the file
%
%field - a struct containing:
%nodeLocs - the locations of the nodes, dims fast, node numbers slow
%times - the times at which field values are recorded
%ux, uy, uz - displacements at nodes in x, y and z (if available) directions
%nodeNums - which node numbers in the original mesh the values correspond to
%
% Written by P. Huthwaite, September 2012
% Updated to include 3D support, 28/3/12 -- PH
% Updated to load to struct April 2014, PH
% Updated to v 1.03 (potentially restrict nodes at which field is saved)
%       with nodeNums field, May 2016, PH
% Minor update (record file precision), March 2018, PH


addExt = 0;
if verLessThan('matlab','9.1')
    if isempty(strfind(fileName,'.')) %#ok<STREMP>
        addExt = 1;
    end
else
    if ~contains(fileName,'.')
        addExt = 1;
    end
end
if addExt
    fileName = [fileName '.pogo-field'];
end
    fid = fopen(fileName,'rb');
    if (fid == -1) 
        error('File %s could not be opened.', fileName)
    end

    header = deblank(fread(fid, 20, '*char')');
    
    if strcmp(header, '%pogo-field1.0')
        fileVer = 1;
    elseif strcmp(header, '%pogo-field1.02')
        fileVer = 1.02;
    elseif strcmp(header, '%pogo-field1.03')
        fileVer = 1.03;        
    else    
		disp(header)
        error('File is wrong format.')
    end
        
    prec = fread(fid, 1, 'int32');
    if prec ~= 4 && prec ~= 8
        error('Specified precision %d incorrect. Should be 4 or 8 bytes.', prec)
    end
    field.prec = prec;
    if prec == 4
    	precStr = 'float32';
    else
        precStr = 'float64';
    end
    nDims = fread(fid, 1, 'int32');
    nDofPerNode = 2;
    if fileVer >= 1.02
        nDofPerNode = fread(fid, 1, 'int32');
        if nDofPerNode ~= 2 && nDofPerNode ~= 3
            error('Unsupported value for nDofPerNode: %d',nDofPerNode);
        end
    end
    
    %nDofPerNode
    
    nNodes = fread(fid, 1, 'int32');
    
    
    nodeLocs = fread(fid, [nDims, nNodes], precStr);
    if fileVer >= 1.03
        field.nodeNums = fread(fid, nNodes, 'int32')+1;
    else
        field.nodeNums = 1:nNodes;
    end
    nFieldStores = fread(fid, 1, 'int32');
    
    times = zeros(nFieldStores,1);
    ux = zeros(nNodes, nFieldStores);
    uy = zeros(nNodes, nFieldStores);
    if nDofPerNode == 3
        uz = zeros(nNodes, nFieldStores);
    end
    for cnt = 1:nFieldStores
        times(cnt) = fread(fid, 1, precStr);
        ux(:,cnt) = fread(fid, nNodes, precStr);
        uy(:,cnt) = fread(fid, nNodes, precStr);
        if nDofPerNode == 3
            uz(:,cnt) = fread(fid, nNodes, precStr);
        end     
    end
    
    fclose(fid);

    field.nodeLocs = nodeLocs;
    field.times = times;
    field.ux = ux;
    field.uy = uy;
    if nDofPerNode == 3
        field.uz = uz;
    end
end

