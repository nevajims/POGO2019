function [ mAbs ] = addAbsBound( m, xLims, yLims, zLims, nAbsVals, K, c, f0, doSrm )
%addAbsBound - add an absorbing boundary to a Pogo input file
%   [ mAbs ] = addAbsBound( m, xLims, yLims, zLims, nAbsVals, K, c, f0, doSrm )
%
%m - model to add absorbing boundary to
%mAbs - model with absorbing boundary added
%xLims, yLims, zLims - limit definitions in each dimension (zLims ignored
%in 2D). Set each to [] if unused (or omit). 
%nAbsVals - number of 'layers' or material values to use. Defaults to 60.
%K - Maximum damping factor. A good parameter will be estimated if
%undefined.
%c - Material velocity used in calculating K. Omit and it will be
%calculated from the first material found in the absorbing region.
%f0 - Frequency used in calculating K. Will be estimated if not provided.
%doSrm - 1 use the stiffness reduction method (more efficient) [default], 0
%           don't use SRM
%
%NB - generally if performance is poor this is because boundary is not big
%enough
%
%Limits definitions
%------------------
%If four values contained in each array, defines the start and end points 
%for the two 'ramps' in order from smallest dimensional value to largest. i.e.
%
%Damping coeff
%   ^ 
% K | ------                                              ------
%   |       --                                          --
%   |         ----                                  ----
%   |             -----                        -----
%   |                  ------            ------
% 0 |-----------------------------------------------------------
%   |       A               B            C                D
%
%xLims(1) corresponds to A, xLims(2) to B and so on. If just two are given,
%these are xLims(1) corresponds to C, xLims(2) corresponds to D (i.e. just 
%a ramp upwards). The same holds for other dimensions.
%
%Generally, need around 1.5 lambda layer width for SRM, 3 lambda for ALID.
%
%Written by P. Huthwaite, 2017
%Updated with SRM, Feb 2018, PH

if nargin < 9 || isempty(doSrm)
    doSrm = 1;
end

if nargin < 8 || isempty(f0)
    %take the first signal and get the dominant frequency
    df = 1/(m.frames{1}.ntSig*m.frames{1}.dtSig);
    spect = abs(fft(m.frames{1}.sigs{1}.sig));
    spect(round(m.frames{1}.ntSig/2):end) = 0;
    [~, ind] = max(spect);
    f0 = (ind-1)*df;
end
if nargin < 7 || isempty(c)
    c = 0; %predict later from materials
end

if nargin < 6 || isempty(K)
    K = 0; %predict later from other params
end
if nargin < 5 || isempty(nAbsVals)
    nAbsVals = 60;
end

if nargin < 4 
    zLims = [];
end
if nargin < 3
    yLims = [];
end


%get element centroids
[xm, ym, zm] = getElCents(m);


l = -1;%length
%upper and lower limits
if isempty(xLims) 
    dampParamXL = zeros(size(xm));
    dampParamXU = zeros(size(xm));
elseif length(xLims) == 2
    %define start and end points
    xl = xLims(1);
    xu = xLims(2);
    l = abs(xl - xu);
    dampParamXL = (xm-xl)/(xu-xl);
    dampParamXL(dampParamXL < 0) = 0;
    dampParamXL(dampParamXL > 1) = 1;

    dampParamXU = zeros(size(xm));
elseif length(xLims) == 4
    %define up, down, down, up
    xl = xLims(2);
    xu = xLims(1);
    l = abs(xl - xu);
    dampParamXL = (xm-xl)/(xu-xl);
    dampParamXL(dampParamXL < 0) = 0;
    dampParamXL(dampParamXL > 1) = 1;

    xl = xLims(3);
    xu = xLims(4);
    dampParamXU = (xm-xl)/(xu-xl);
    dampParamXU(dampParamXU < 0) = 0;
    dampParamXU(dampParamXU > 1) = 1;
else
    error('xLims poorly defined')
end

if isempty(yLims) 
    dampParamYL = zeros(size(ym));
    dampParamYU = zeros(size(ym));
elseif length(yLims) == 2
    %define start and end points
    yl = yLims(1);
    yu = yLims(2);
    if l < 0
        l = abs(yl - yu);
    end
    dampParamYL = (ym-yl)/(yu-yl);
    dampParamYL(dampParamYL < 0) = 0;
    dampParamYL(dampParamYL > 1) = 1;

    dampParamYU = zeros(size(ym));
elseif length(yLims) == 4
    %define up, down, down, up
    yl = yLims(2);
    yu = yLims(1);
    if l < 0
        l = abs(yl - yu);
    end
    dampParamYL = (ym-yl)/(yu-yl);
    dampParamYL(dampParamYL < 0) = 0;
    dampParamYL(dampParamYL > 1) = 1;

    yl = yLims(3);
    yu = yLims(4);
    dampParamYU = (ym-yl)/(yu-yl);
    dampParamYU(dampParamYU < 0) = 0;
    dampParamYU(dampParamYU > 1) = 1;
else
    error('yLims poorly defined')
end


if m.nDims == 3
    if isempty(zLims) 
        dampParamZL = zeros(size(zm));
        dampParamZU = zeros(size(zm));
    elseif length(zLims) == 2
        %define start and end points
        zl = zLims(1);
        zu = zLims(2);
        if l < 0
            l = abs(zl - zu);
        end
        dampParamZL = (zm-zl)/(zu-zl);
        dampParamZL(dampParamZL < 0) = 0;
        dampParamZL(dampParamZL > 1) = 1;

        dampParamZU = zeros(size(zm));
    elseif length(zLims) == 4
        zl = zLims(2);
        zu = zLims(1);
        if l < 0
            l = abs(zl - zu);
        end
        dampParamZL = (zm-zl)/(zu-zl);
        dampParamZL(dampParamZL < 0) = 0;
        dampParamZL(dampParamZL > 1) = 1;

        zl = zLims(3);
        zu = zLims(4);
        dampParamZU = (zm-zl)/(zu-zl);
        dampParamZU(dampParamZU < 0) = 0;
        dampParamZU(dampParamZU > 1) = 1;
    else
        error('zLims poorly defined')
    end
end

if m.nDims == 2
    dampParam = max([dampParamXL(:) dampParamXU(:) dampParamYL(:) dampParamYU(:)].');
else
    dampParam = max([dampParamXL(:) dampParamXU(:) dampParamYL(:) dampParamYU(:) dampParamZL(:) dampParamZU(:)].');
    %dampParam = dampParamXL.*dampParamXU.*dampParamYL.*dampParamYU.*dampParamZL.*dampParamZU;
end

changeEls = find(dampParam < 1 & dampParam > 0);

if isempty(changeEls)
    error('No elements found to change')
end

%check materials present in the area already
changeMats = unique(m.matTypeRefs(changeEls));
if isempty(changeMats)
    error('No materials found to change')
end

nChangeMats = length(changeMats);

%get params from each material
if K == 0
    
    K = zeros(nChangeMats,1);
    E = zeros(nChangeMats,1);
    for matType = 1:nChangeMats
        matNum = changeMats(matType);
        if c == 0 && m.matTypes{matNum}.paramsType ~= 0
            error('Need a velocity value input if not using isotropic materials')
        end
        E(matType) = m.matTypes{matNum}.paramValues(1);
        nu = m.matTypes{matNum}.paramValues(2);
        rho = m.matTypes{matNum}.paramValues(3);
        c0 = sqrt(E(matType)*(1-nu)/(rho*(1+nu)*(1-2*nu)));
        % cSh = sqrt(E/(2*rho*(1+nu)));
        c = c0;


        if doSrm
            K(matType) = 2*pi*f0;
        else
            K(matType) = getDamping(0.01,3,l,c,f0);
        end
        fprintf('K parameter for material %d set to %6.4g per second\n',changeMats(matType),K(matType))
        %end
    end
end


nNewMats = nChangeMats*nAbsVals;
nInitMats = length(m.matTypes);

%define new materials
P = 3;
%alphaVals = ((1:nAbsVals)/nAbsVals).^3;
Xparam = ((1:nAbsVals)/nAbsVals);
newMats = cell(1,nNewMats);

for mCnt = 1:nChangeMats
    mNum = changeMats(mCnt);
    pType = m.matTypes{mNum}.paramsType;
    %get param values, stripping off damping if present
    if pType == 0
        prevParams = m.matTypes{mNum}.paramValues(1:3);
        denNum = 3;
    elseif pType == 1
        prevParams = m.matTypes{mNum}.paramValues(1:10);
        denNum = 10;
        if doSrm
            error('Only isotropic materials supported for SRM')
        end
    else
        error('Unrecognised material type')
    end
    for aCnt = 1:nAbsVals
        newMatNum = (mCnt-1)*nAbsVals+aCnt;
        newMats{newMatNum}.paramsType = pType;

        newParams = prevParams;
        
        k0 = 2*pi*f0/c0;
        
        cVal = Xparam(aCnt).^P*K(mCnt);
        
        if doSrm
            alphaMax = -log(0.01)/(k0*l);
            newParams(1) = prevParams(1)*exp(-(alphaMax*Xparam(aCnt).^P)*k0*(Xparam(aCnt)*l));
        end
        
   
        newMats{newMatNum}.paramValues = [newParams(:).' cVal];
                
        newMats{newMatNum}.parent = mNum;
    end
end
dampParam = reshape(dampParam,size(m.matTypeRefs));

newMatTypeRefs = nInitMats+(m.matTypeRefs(changeEls(:))-1)*nAbsVals...
    +round(dampParam(changeEls(:))*(nAbsVals)+0.5);

mAbs = m;

mAbs.matTypes = [m.matTypes(:).' newMats];
mAbs.matTypeRefs(changeEls) = newMatTypeRefs;


end

