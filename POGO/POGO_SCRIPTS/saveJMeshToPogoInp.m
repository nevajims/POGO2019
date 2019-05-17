clear

%DOF
%Damping

%POGO_PATH = 'C:\Users\Alva\Desktop\TTi\POGO\pogoMatlab\POGO_SCRIPTS\';
%m = loadVtuMesh(strcat(POGO_PATH, 'temp.vtu') );

PATH = 'POGO_SCRIPTS\';
OUTFILE = '3Djimrail2.pogo-inp';
JIM_DATA_INPATH = 'JIM\meshes\MESH_RAIL_56_1.mat';

jmesh = load(JIM_DATA_INPATH);

%define element size (was dx = 0.2e-3;)
%take mean from mesh
ed = get_element_details(jmesh.mesh);
dx = ed.mean_edge_length;

jm = convertJMesh2DPogoMesh( jmesh );

%get point on top
top_value = max(jm.nodePos(2,:));
top_point = find(jm.nodePos(2,:) == top_value)

figure
plot(jm.nodePos(1,:),jm.nodePos(2,:),'k.')
hold on
plot(jm.nodePos(1,top_point),jm.nodePos(2,top_point),'ro')

z_length = 1;                   %rail length metres
n_slices = floor(z_length/dx);  %duplicate in z direction n times

m = extrudePogo(jm, n_slices, dx);
%C3D6R (triangle extruded to trianglular prism)

%sit source on top point on slices 3
top_point_10_slice = (10 - 1) * size(jm.nodePos,2) + top_point;
%sit receiver on top of mid slice
mid_slice = floor((n_slices + 1)/ 2)
top_point_mid_slice = (mid_slice - 1) * size(jm.nodePos,2) + top_point;

%MODEL VARIABLES
 
%source and receiver points
sPoint = top_point_10_slice;
rPoint = top_point_mid_slice;

%crude display rail nodes. Source node red. Receiver green.
figure
rotate3d on
scatter3(m.nodePos(1,:),m.nodePos(2,:),m.nodePos(3,:), 'b.')
hold on
scatter3(m.nodePos(1,sPoint),m.nodePos(2,sPoint),m.nodePos(3,sPoint), 'ro', 'MarkerFaceColor', 'r')
scatter3(m.nodePos(1,rPoint),m.nodePos(2,rPoint),m.nodePos(3,rPoint), 'go', 'MarkerFaceColor', 'g')

%define material properties (material Isotropic)
E = 210e9;
G = 80e9;
nu = E/2/G-1;
rho = 8000;

%Courant-Friedrichs-Lewy (CFL) condition (#elements / timestep)
courant = 0.3; % mix of sizes so be conservative with stability

%define frequency
freq = 60;
%number of cycles in toneburst
nCyc = 3;
%sound speeds
c0 = sqrt(E*(1-nu)/(rho*(1+nu)*(1-2*nu)));
cSh = sqrt(E/(2*rho*(1+nu)));

%get wavelengths
lambda = c0/freq;
lambdaSh = cSh/freq;

%Add model variable to model

%set dt and nt
m.dt = dx/c0*courant; 
%want 30e-6 s long simulation (would be 5 secs in our model)
duration_secs = 5
m.nt = round(duration_secs/m.dt);

%get number of elements
nEls = size(m.elNodes,2);
%define the material references (set them all to material 1)
m.matTypeRefs = ones(nEls,1);
%then define material 1 with the properties from before
%Type 0 = Isotropic (Parameters: E, nu, rho, [alpha (damping) optional])
m.matTypes{1}.paramsType = 0;
m.matTypes{1}.paramValues = [E, nu, rho];

%SIGNAL

%put a Hann signal on sig 1, frame 1
m = genPogoHannSignal(m,nCyc,freq,1,1);

%say what nodes to apply it to
%Pogo manual p47
%* NB if saving to v1.08 or earlier: while nodes can appear multiple times in nodeSpec,
%they must be applied to different DOF. Only a single signal can be applied to each
%DOF in the model.
%* From v1.09, these restrictions are lifted.
m.frames{1}.sigs{1}.nodeSpec = [sPoint sPoint sPoint];

%and which degrees of freedom
%ignoring y to make demo faster
m.frames{1}.sigs{1}.dofSpec = [1 2 3];

%set the signal amplitudes in each direction (x,z)
m.frames{1}.sigs{1}.sigAmps = [-0.3 0 0.6];
%set the type to 0 = force (1 = displacement)
m.frames{1}.sigs{1}.sigType = 0;

%define the measurement frequency and start time
m.measFreq = 2;
m.measStart = 1;

%define the measurement sets
m.measSets{1}.name = 'main';
m.measSets{1}.measNodes = [rPoint,rPoint,rPoint];
m.measSets{1}.measDof = [1 2 3];
%receive in x,z at node rPoint for demo simplicity

%calculate time increments to store the field at
gap = round(m.nt/80); % do at around 80 time steps
if gap < 1
    gap = 1;
end
m.fieldStoreIncs = 1:gap:m.nt;

%{
%put an absorbing boundary on the bottom edge
xLims = [];
yLims = [-0.1 -0.1+0.012  100 102];
%^only use first two values here - the last two are set to very large so
%that they are outside the model. Here we're using 3 wavelength wide
%boundaries.

m3 = addAbsBound(m3,xLims,yLims,[],[],[],c0,freq);
%}

%save the file - this is done in the slightly older v1.07 format

savePogoInp(strcat(PATH, OUTFILE),m);

%}

%}