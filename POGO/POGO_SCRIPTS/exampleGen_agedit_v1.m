clear

MESH = 0;
PATH = 'POGO_SCRIPTS\';

%define element size
dx = 0.2e-3;
%lambda = 4mm so this is 20 els per wavelength at 1 MHz

%define nodes along the bottom of the free meshed domain
nodex = (0.1:dx:0.2).';

%put into a points array, along with another point at the top to make a triangle
points = [ [ nodex zeros(length(nodex),1)]
    0.15 0.06].';

%add points inside for a rectangle to be deleted
rectPoints = [ 0.14 0.01
    0.16 0.01
    0.16 0.025
    0.14 0.025].';
%combine:
points = [points rectPoints];

nPointsOuter = length(nodex)+1;

%define the segments to joint the points around the outside together
segsOut = [ (1:nPointsOuter).' [2:nPointsOuter 1].'].';
%this matrix will be:
%1,2
%2,3
%3,4
% ...
%n-1, n
%n, 1

%segments for the internal rectangle
segsRect = [(nPointsOuter+1:nPointsOuter+4).' [(nPointsOuter+2:nPointsOuter+4) nPointsOuter+1].'].';
%combined segments
segs = [segsOut segsRect];
    
%one hole which will be deleted
holes = [0.148 0.015].';

if MESH==1

%save into a poly format
savePoly( strcat(PATH, 'temp.poly'), points, segs, holes );

%IMPORTANT
%This might need to be changed depending on your system setup
%Essentially you want to run pogoMesh with the file just generated
%and maximum length set to dx*1.5
%system(sprintf('pogoMesh temp.poly -l %f',dx*1.5))

%MANUALLY
%sprintf('pogoMesh temp.poly -s %f',dx)

sprintf('pogoMesh temp.poly -l %f',dx*1.5);


%having done the meshing can delete the file
%delete('temp.poly')

else

%load in mesh
m = loadVtuMesh(strcat(PATH, 'temp.vtu') );
%can delete mesh file now loaded
%delete('temp.vtu')

%generate a 2D grid
%centre locations for grid 
% - note that we want to make the top of this line up with the nodes from
% the free mesh
cx = 0.15;
cy = -0.05;
%points in x and y
nx = 0.1/dx+1;
ny = 0.1/dx+1;
%actually make grid
m2 = genGrid2D(nx,ny,dx,dx,cx,cy);

%get the overlapping node numbers from each model; the top line from 2, the
%bottom line from 1
n2 = (1:nx) + (ny-1)*nx;
n1 = 1:nx;

if 0
    %can plot the two models if desired
    figure
    plot(m.nodePos(1,:),m.nodePos(2,:),'k.')
    hold on
    plot(m2.nodePos(1,:),m2.nodePos(2,:),'bx')
    plot(m.nodePos(1,n1),m.nodePos(2,n1),'go')
    plot(m2.nodePos(1,n2),m2.nodePos(2,n2),'r+')
    axis equal
end

m3 = combinePogoMesh(m,m2,n1,n2);

%get number of nodes in model 1, the free meshed model
nNodes1 = size(m.nodePos,2);

%source and receiver points - offset by nNodes1 from the start, then simple
%arithmetic in the grid to set two locations
sPoint = nNodes1+(nx+1)/2+(ny+1)/2*nx;
rPoint = nNodes1+(nx+1)/2+(ny-5)*nx;

%define material properties
E = 210e9;
G = 80e9;
nu = E/2/G-1;
rho = 8000;

courant = 0.3; % mix of sizes so be conservative with stability

%define frequency
freq = 1500e3;
%number of cycles in toneburst
nCyc = 3;
%sound speeds
c0 = sqrt(E*(1-nu)/(rho*(1+nu)*(1-2*nu)));
cSh = sqrt(E/(2*rho*(1+nu)));

%get wavelengths
lambda = c0/freq;
lambdaSh = cSh/freq;

%set dt and nt
m3.dt = dx/c0*courant; 
%want 30e-6 s long simulation
m3.nt = round(30e-6/m3.dt);

%get number of elements
nEls = size(m3.elNodes,2);
%define the material references (set them all to material 1)
m3.matTypeRefs = ones(nEls,1);
%then define material 1 with the properties from before
m3.matTypes{1}.paramsType = 0;
m3.matTypes{1}.paramValues = [E, nu, rho];

%put a Hann signal on sig 1, frame 1
m3 = genPogoHannSignal(m3,nCyc,freq,1,1);

%say what nodes to apply it to
m3.frames{1}.sigs{1}.nodeSpec = [sPoint sPoint];
%and which degrees of freedom
m3.frames{1}.sigs{1}.dofSpec = [1 2];

%set the signal amplitudes in each direction
m3.frames{1}.sigs{1}.sigAmps = [-0.3 0.6];
%set the type to 0, for force
m3.frames{1}.sigs{1}.sigType = 0;

%define the measurement frequency and start time
m3.measFreq = 2;
m3.measStart = 1;

%define the measurement sets
m3.measSets{1}.name = 'main';
m3.measSets{1}.measNodes = [rPoint,rPoint];
m3.measSets{1}.measDof = [1 2];
%receive in x and y at node rPoint

%calculate time increments to store the field at
gap = round(m3.nt/80); % do at around 80 time steps
if gap < 1
    gap = 1;
end
m3.fieldStoreIncs = 1:gap:m3.nt;

%delete some elements
[ex,ey] = getElCents(m3); %get element centroids 
%find which elements to delete based on coordinates
elsDel = find(ex < 0.15 & ex > 0.14 & ey > -0.09 & ey < -0.08);
%delete them
m3 = deleteEls(m3,elsDel);

%put an absorbing boundary on the bottom edge
xLims = [];
yLims = [-0.1 -0.1+0.012  100 102];
%^only use first two values here - the last two are set to very large so
%that they are outside the model. Here we're using 3 wavelength wide
%boundaries.

m3 = addAbsBound(m3,xLims,yLims,[],[],[],c0,freq);

%duplicate some nodes to decouple some elements - model a very thin crack
px = m3.nodePos(1,:);
py = m3.nodePos(2,:);
nodesDupe = find(px < 0.1801 & px > 0.1799 & py > -0.05 & py < -0.04);
n = nodesDupe;
m3 = decoupleEls(m3,nodesDupe,[1 0]); 
%direction set to x direction

%save the file - this is done in the slightly older v1.07 format
%savePogoInp('example.pogo-inp',m3,1.07);
savePogoInp(strcat(PATH, 'example.pogo-inp'),m3);

end


