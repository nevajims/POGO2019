function [ model ] = loadPogoInp( fileName )
% function [ model ] = loadPogoInp( fileName )
% see savePogoInp for full documentation

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
    fileName = [fileName '.pogo-inp'];
end

fid = fopen(fileName,'rb','n','UTF-8');
if (fid == -1) 
    error('File could not be opened.')
    %return;
end

rawRead = fread(fid, 20, 'uint8').';
nullTerm = find(rawRead == 0,1);
if ~isempty(nullTerm)
    rawRead(nullTerm:end) = 0;
end
rawRead = char(rawRead);
rawRead = deblank(rawRead);
        
header = rawRead;
header = deblank(header);
fileVer = -1;
if strcmp(header, '%pogo-inp1.0') == 1
    fileVer = 1;
end
if strcmp(header, '%pogo-inp1.01') == 1
    fileVer = 1.01;
end
if strcmp(header, '%pogo-inp1.02') == 1
    fileVer = 1.02;
end
if strcmp(header, '%pogo-inp1.03') == 1
    fileVer = 1.03;
end
if strcmp(header, '%pogo-inp1.04') == 1
    fileVer = 1.04;
end
if strcmp(header, '%pogo-inp1.05') == 1
    fileVer = 1.05;
end
if strcmp(header, '%pogo-inp1.06') == 1
    fileVer = 1.06;
end
if strcmp(header, '%pogo-inp1.07') == 1
    fileVer = 1.07;
end
if strcmp(header, '%pogo-inp1.08') == 1
    fileVer = 1.08;
end
if strcmp(header, '%pogo-inp1.09') == 1
    fileVer = 1.09;
end
if strcmp(header, '%pogo-inp1.10') == 1
    fileVer = 1.10;
end
% if strcmp(header, '%pogo-inp1.11') == 1
%     fileVer = 1.11;
% end
if fileVer == -1
    %fileVer = 1;
    
    %disp(fileVer)
    %disp(header)
    error('File is wrong format. header: %s.', header)
end
%disp(fileVer)
model.fileVer = fileVer;

maxUint = uint32(uint64(2)^32-1);

if fileVer >= 1.10
    dofSaveString = 'uint64'; 
else
    dofSaveString = 'int32'; 
end

model.prec = fread(fid, 1, 'int32');
if model.prec == 8
    precString = 'float64';
else
    precString = 'float32';
end
%model.prec
model.nDims = fread(fid, 1, 'int32');

if fileVer >= 1.02
    model.nDofPerNode = fread(fid, 1, 'int32');
else
    %nDofPerNode = -1;
end

rawRead = fread(fid, 1024, 'uint8').';
nullTerm = find(rawRead == 0,1);
if ~isempty(nullTerm)
    rawRead(nullTerm:end) = 0;
end
rawRead = char(rawRead);
rawRead = deblank(rawRead);
        
notes = rawRead;
if ~isempty(notes)
    model.notes = notes;
end


rawRead = fread(fid, 80, 'uint8').';
nullTerm = find(rawRead == 0,1);
if ~isempty(nullTerm)
    rawRead(nullTerm:end) = 0;
end
rawRead = char(rawRead);
rawRead(~isstrprop(rawRead,'alphanum')) = '';
rawRead = deblank(rawRead);
model.runName = rawRead;

model.nt = fread(fid, 1, 'uint32');
model.dt = fread(fid, 1, precString);

nNodes = fread(fid, 1, 'uint32');
model.nodePos = fread(fid, [model.nDims,nNodes], precString);

nEls = fread(fid, 1, 'uint32');
nNodesPerEl = fread(fid, 1, 'int32');
model.elTypeRefs = fread(fid, nEls, 'int32')+1;
model.matTypeRefs = fread(fid, nEls, 'int32')+1;
model.orientRefs = fread(fid, nEls, 'int32')+1;

if fileVer >= 1.10
    elNodesRead = fread(fid, [nNodesPerEl, nEls], 'uint32');
    model.elNodes = elNodesRead+1;
    model.elNodes(elNodesRead == uint64(2)^32-1) = 0;
else
    model.elNodes = fread(fid, [nNodesPerEl, nEls], 'int32')+1;
end

if fileVer >= 1.03
    nPmlSets = fread(fid, 1, 'int32');
    if nPmlSets ~= 0
        nPmlSets
        error('nPmlSets > 0, but PML reading unsupported by this version.')
    end

    nPmlParams = fread(fid, 1, 'int32');
    if nPmlParams ~= 0
        error('nPmlParams > 0, but PML reading unsupported by this version.')
    end
end


if feof(fid)
    %error('Unexpected EOF found')
    disp('Unexpected EOF found')
    dbstack
    fclose(fid);
    return
end

nElTypes = fread(fid,1,'int32');
model.elTypes = cell(nElTypes,1);
for eCnt = 1:nElTypes
    %rawRead = fread(fid, 20, '*char').';
    rawRead = fread(fid, 20, 'uint8').';
    elType = rawRead;

    %processing in case it's a messed up string!
    nullTerm = find(elType == 0,1);
    if ~isempty(nullTerm)
        elType(nullTerm:end) = 0;
    end
    elType = char(elType);
    elType(~isstrprop(elType,'alphanum')) = '';
    elType = deblank(elType);
    
    model.elTypes{eCnt}.name = elType;
    model.elTypes{eCnt}.paramsType = fread(fid,1,'int32');
    nElParams = fread(fid,1,'int32');
    if nElParams > 0
        model.elTypes{eCnt}.paramValues = fread(fid,nElParams,precString);
    end
end


if feof(fid)
    %error('Unexpected EOF found')
    disp('Unexpected EOF found')
    dbstack
    fclose(fid);
    return
end

nMatTypes = fread(fid,1,'int32');
model.matTypes = cell(nMatTypes,1);
for eCnt = 1:nMatTypes
	if fileVer >= 1.09
		model.matTypes{eCnt}.parent = fread(fid,1,'int32')+1;
	end
    model.matTypes{eCnt}.paramsType = fread(fid,1,'int32');
    nMatParams = fread(fid,1,'int32');
    model.matTypes{eCnt}.paramValues = fread(fid,nMatParams,precString);
    %return 
end

nOr = fread(fid,1,'int32');
if nOr > 0
    model.or = cell(nOr,1);
    for eCnt = 1:nOr
        model.or{eCnt}.paramsType = fread(fid,1,'int32');
        nOrParams = fread(fid,1,'int32');
        model.or{eCnt}.paramValues = fread(fid,nOrParams,precString);
    end
end

if feof(fid)
    %error('Unexpected EOF found')
    disp('Unexpected EOF found')
    dbstack
    fclose(fid);
    return
end

if fileVer >= 1.10
    nFixDof = fread(fid,1,'uint64');
else
    nFixDof = fread(fid,1,'uint32');
end

if nFixDof > 0
    fixDof = fread(fid,nFixDof,dofSaveString);
    model.fixDof = mod(fixDof,4)+1;
    model.fixNodes = floor(fixDof/4)+1;
end

if fileVer >= 1.06
    nTieSets = fread(fid,1,'int32');
    if nTieSets > 0
        model.ties = cell(nTieSets,1);

        for tCnt = 1:nTieSets
            model.ties{tCnt}.tieTransform = fread(fid,[model.nDofPerNode, model.nDofPerNode],precString);
            nTies = fread(fid,1,'int32');
            model.ties{tCnt}.masterNodes = fread(fid,[nTies,1],'int32');
            model.ties{tCnt}.slaveNodes = fread(fid,[nTies,1],'int32');
        end
    end
end
        
%DOF groups:
if fileVer >= 1.08
    nDofGroups = fread(fid,1,'int32');
    for dCnt = 1:nDofGroups
        nDof = fread(fid,1,'int32');
        dofSpec = fread(fid,nDof,dofSaveString);
        model.dofGroups{dCnt}.dofSpec = mod(dofSpec,4)+1;
        model.dofGroups{dCnt}.nodeSpec = floor(dofSpec/4)+1;
        model.dofGroups{dCnt}.dofWeight = fread(fid,nDof,precString);
    end
else
    nDofGroups = 0;
end

%frames:
if fileVer >= 1.07
    nFrames = fread(fid,1,'int32');
else
    nFrames = 1;
end

model.frames = cell(nFrames,1);
for fCnt = 1:nFrames
    %sigs:
    nSigs = fread(fid,1,'int32');
    model.frames{fCnt}.ntSig = fread(fid,1,'int32');
    model.frames{fCnt}.dtSig = fread(fid,1,precString);

    model.frames{fCnt}.sigs = cell(nSigs,1);
    for sCnt = 1:nSigs
        nDofForSig = fread(fid,1,'int32');
        if fileVer >= 1.01
            model.frames{fCnt}.sigs{sCnt}.sigType = fread(fid,1,'int32');
        end
        if fileVer >= 1.08
            model.frames{fCnt}.sigs{sCnt}.isDofGroup = fread(fid,1,'int8');
        else
            model.frames{fCnt}.sigs{sCnt}.isDofGroup = 0;
        end
        dofSpec = fread(fid,nDofForSig,dofSaveString);
        if model.frames{fCnt}.sigs{sCnt}.isDofGroup == 0
            model.frames{fCnt}.sigs{sCnt}.dofSpec = mod(dofSpec,4)+1;
            model.frames{fCnt}.sigs{sCnt}.nodeSpec = floor(dofSpec/4)+1;
        else
            model.frames{fCnt}.sigs{sCnt}.dofGroup = dofSpec+1;
        end

        model.frames{fCnt}.sigs{sCnt}.sigAmps = fread(fid,nDofForSig,precString);
        model.frames{fCnt}.sigs{sCnt}.sig = fread(fid,model.frames{fCnt}.ntSig,precString);
    end

    if feof(fid)
        %error('Unexpected EOF found')
        disp('Unexpected EOF found')
        dbstack
        fclose(fid);
        return
    end
end

%meas hist:
if fileVer >= 1.05
    nMeasSets = fread(fid,1,'int32');
    model.measSets = cell(nMeasSets,1);
    
    model.measFreq = fread(fid,1,'int32');
    model.measStart = fread(fid,1,'int32')+1;
    
    for mSetCnt = 1:nMeasSets
        
        rawRead = fread(fid, 20, 'uint8').';
        nullTerm = find(rawRead == 0,1);
        if ~isempty(nullTerm)
            rawRead(nullTerm:end) = 0;
        end
        rawRead = char(rawRead);
        rawRead(~isstrprop(rawRead,'alphanum')) = '';
        rawRead = deblank(rawRead);
        
        name = rawRead;%fread(fid, 20, '*char')';
        model.measSets{mSetCnt}.name = deblank(name);
    
        nMeas = fread(fid,1,'int32');
        if fileVer >= 1.08
            model.measSets{mSetCnt}.isDofGroup = fread(fid,1,'int8');
        else
            model.measSets{mSetCnt}.isDofGroup = 0;
        end
        
        measDof = fread(fid,nMeas,dofSaveString);
        if model.measSets{mSetCnt}.isDofGroup == 0
            model.measSets{mSetCnt}.measDof = mod(measDof,4)+1;
            model.measSets{mSetCnt}.measNodes = floor(measDof/4)+1;
        else
            model.measSets{mSetCnt}.dofGroup = measDof + 1;
        end
    end
    
    
else
    nMeas = fread(fid,1,'int32');
    model.measFreq = fread(fid,1,'int32');
    if fileVer >= 1.04
        model.measStart = fread(fid,1,'int32')+1;
    else
        model.measStart = 1;
    end
    model.measSets = cell(1);
    measDof = fread(fid,nMeas,dofSaveString);
    model.measSets{1}.name = 'main';
    model.measSets{1}.measDof = mod(measDof,4)+1;
    model.measSets{1}.measNodes = floor(measDof/4)+1;
end

%meas field:
nFieldStores = fread(fid,1,'int32');
if nFieldStores > 0
    model.fieldStoreIncs = fread(fid,nFieldStores,'int32')+1;
end 

if fileVer >= 1.05
    nFieldNodes = fread(fid,1,'uint32');
    if nFieldNodes ~= maxUint
        model.fieldStoreNodes = fread(fid,nFieldNodes,'uint32')+1;
    end
end

fclose(fid);

end

