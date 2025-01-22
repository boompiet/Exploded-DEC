% -------------------------------------------------------------------------
% Conduction analysis using extended DEC formulation
% -------------------------------------------------------------------------
% Purpose:      Solve conduction problems on Voronoi complexes and their 
%               low order skeletons using extended formulation of discrete
%               exterior calculus
%
% Institution:  King Fahd University of Petroleum and Minerals
%
% Author:       Dr Pieter Boom
% Date:         2025/01/07
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% preliminaries
% -------------------------------------------------------------------------
clc; close all; clear all

% -------------------------------------------------------------------------
% create Voronoi complex structure
% -------------------------------------------------------------------------
disp('>> create Voronoi complex structure'); vcmplx = setup_complex(4);
vcmplx(1).name='vertex'; vcmplx(2).name='edge';
vcmplx(3).name='face';   vcmplx(4).name='polyhedron';

% -------------------------------------------------------------------------
% read Neper .tess file for Voronoi complex
% -------------------------------------------------------------------------
filename = 'n1000-id1.tess';
disp(['>> read Neper .tess file (',filename,') for Voronoi complex'])
[vcmplx] = read_tess(filename,vcmplx);

% -------------------------------------------------------------------------
% get Voronoi geometry
% -------------------------------------------------------------------------
disp('>> get Voronoi geometry'); vcmplx = vgeo(vcmplx);
% boundary operators
vBO_1 = build_BO(vcmplx,1); vBO_2 = build_BO(vcmplx,2);
vBO_3 = build_BO(vcmplx,3); sf3 = abs(sum(vBO_3)'); 
% surface geometry indicators
se2 = abs(vBO_2)'*sf3; se2 = se2>1; sn1 = abs(vBO_1)'*se2; sn1 = sn1>1; 
sv4 = abs(vBO_3)*sf3; sv4 = sv4>1;

% -------------------------------------------------------------------------
% build Voronoi connectivity
% -------------------------------------------------------------------------
disp('>> build Voronoi connectivity'); vcmplx = vcon(vcmplx);

% -------------------------------------------------------------------------
% get Delaunay geometry
% -------------------------------------------------------------------------
disp('>> build Delaunay dual geometry'); vcmplx = dgeo(vcmplx);

% -------------------------------------------------------------------------
% create Forman extended complex structure
% -------------------------------------------------------------------------
disp('>> create Forman complex structure'); fcmplx = setup_complex(4);
fcmplx(1).name='vertex'; fcmplx(2).name='edge';
fcmplx(3).name='face';   fcmplx(4).name='polyhedron';

% -------------------------------------------------------------------------
% setup Forman extended complex from Voronoi complex
% -------------------------------------------------------------------------
disp('>> setup Forman complex from Voronoi complex')
fcmplx = v2f(fcmplx,vcmplx);

% -------------------------------------------------------------------------
% build Forman extended complex connectivity
% -------------------------------------------------------------------------
disp('>> build Forman connectivity'); fcmplx = fcon(fcmplx);
% boundary operators and curvature
fBO_1 = build_BO(fcmplx,1); omega = boundary_curvature(fcmplx);

%% -------------------------------------------------------------------------
% graphite model for conductivity and interface area of extended complex
% -------------------------------------------------------------------------
% Set conductivities along Forman edges
% *** this is currently setup as a constant value for each edge type, but
%     could easily be modified to have variable values
k1 = (1.0*10^0)*ones(sum(vcmplx(1).num(2).val),1); % V node-edge cond
k2 = (1.0*10^0)*ones(sum(vcmplx(2).num(3).val),1); % V edge-face cond
k3 = (1.0*10^0)*ones(sum(vcmplx(3).num(4).val),1); % V face-volume cond
k  = spdiags([k1; k2; k3],0,length([k1; k2; k3]),length([k1; k2; k3]));

% Set geometry parameters
% *** this is currently setup as a constant value for each dimension, but 
%     could easily be modified to have variable values
%-- Voronoi node volume
cnt1 = 1; cnt2 = vcmplx(1).num(1).val;
node_diam       = 10^-1*ones(vcmplx(1).num(1).val,1); 
node_vol        = pi/6*node_diam.^3./omega(cnt1:cnt2);
%-- Voronoi edge cross-sectional area
cnt1 = cnt2+1; cnt2 = cnt2+vcmplx(2).num(2).val;
edge_diam       = 10^-2*ones(vcmplx(2).num(2).val,1); 
node_edge_area  = pi/4*edge_diam.^2./omega(cnt1:cnt2);
%-- Voronoi face thickness
cnt1 = cnt2+1; cnt2 = cnt2+vcmplx(3).num(3).val;
face_thic       = 10^-3*ones(vcmplx(3).num(3).val,1)./omega(cnt1:cnt2);

% -------------------------------------------------------------------------
% get extended DEC geometry and Hodge stars
% -------------------------------------------------------------------------
disp('>> build extended DEC geometric operators'); 
vcmplx = dgeo_mod(vcmplx,face_thic,node_edge_area,node_vol);
%vcmplx(1).dvol_mod=vcmplx(1).dvol; vcmplx(2).dvol_mod=vcmplx(2).dvol;
%vcmplx(3).dvol_mod=vcmplx(3).dvol; vcmplx(4).dvol_mod=vcmplx(4).dvol;
[star1,star2] = dvgeo(vcmplx);
% diagonal matrix of geometric Hodge stars
ibs1 = spdiags(1./star1,0,length(star1),length(star1));
bs2  = spdiags(star2,0,length(star2),length(star2));

% -------------------------------------------------------------------------
% get specimen max/min dimensions
% -------------------------------------------------------------------------
max_x = max(vcmplx(1).cc(:,1)); min_x = min(vcmplx(1).cc(:,1));
max_y = max(vcmplx(1).cc(:,2)); min_y = min(vcmplx(1).cc(:,2));
max_z = max(vcmplx(1).cc(:,3)); min_z = min(vcmplx(1).cc(:,3));

% -------------------------------------------------------------------------
% sort arrays for plotting and calculation purposes
% -------------------------------------------------------------------------
disp('>> sort arrays for plotting and calculation purposes')
small=10^-8; z  = fcmplx(1).cc(:,3); 

% Forman vertices
pindx = [1:vcmplx(1).num(1).val];
eindx = pindx(end) + [1:vcmplx(2).num(2).val];
findx = eindx(end) + [1:vcmplx(3).num(3).val];
vindx = findx(end) + [1:vcmplx(4).num(4).val];
[pz,pzi] = sort(z(pindx)); [ez,ezi] = sort(z(eindx));
[fz,fzi] = sort(z(findx)); [vz,vzi] = sort(z(vindx));
vindx = findx(end) + vzi; findx = eindx(end) + fzi;
eindx = pindx(end) + ezi; pindx = pzi;

% Forman edges
fcmplx(2).cc = 0.5*abs(fBO_1)*fcmplx(1).cc;
peindx = [1:sum(vcmplx(1).num(2).val)];
efindx = peindx(end) + [1:sum(vcmplx(2).num(3).val)];
fvindx = efindx(end) + [1:sum(vcmplx(3).num(4).val)];
z2 = fcmplx(2).cc(:,3); [pez,pezi] = sort(z2(peindx)); 
[efz,efzi] = sort(z2(efindx));[fvz,fvzi] = sort(z2(fvindx));
fvindx = efindx(end) + fvzi; efindx = peindx(end) + efzi; peindx = pezi;

% flux points
zsf_v = zeros(size(fcmplx(1).cc(:,3)));  
zsf_v(fcmplx(1).cc(:,3)<min_z+small)=1;  
zsf_e = abs(fBO_1)*zsf_v;  zsf_e(zsf_e>1)=0;
osf_v = zeros(size(fcmplx(1).cc(:,3)));  
osf_v(fcmplx(1).cc(:,3)>max_z-small)=1;  
osf_e = abs(fBO_1)*osf_v;  osf_e(osf_e>1)=0;

% -------------------------------------------------------------------------
% setup system
% -------------------------------------------------------------------------
disp('>> setup system and solve problem')
int_indx = find(z~=max_z & z~=min_z);
A = ibs1 * fBO_1' * bs2 * k * fBO_1; rhs = - A*(z==max_z);
t_old = zeros(fcmplx(1).num(1).val,1); t_old(z==max_z) = 1;

% -------------------------------------------------------------------------
% integrate forward in time
% -------------------------------------------------------------------------
n = 30;              % number of time steps
h = 0.05;            % Time step size 

figure;
AI = speye(size(A(int_indx,int_indx)));
AR = A(int_indx,int_indx);
rhsR = rhs(int_indx);
for t=h:h:(h*n)
    disp(['Time: ',num2str(t),' of max time ', num2str(h*n)])
    t_new = (AI + h*AR) \ ( t_old(int_indx) + h*rhsR);
    t_old(int_indx) = t_new;

    % solution 2D
    subplot(2,2,1); plot(z(pindx),t_old(pindx),'k*',z2,z2./max(z2),'r-');
    xlabel('z-coordinate'); ylabel('Values'); title('Vertices'); 
    legend('comp','exact','Location','northwest');
    subplot(2,2,2); plot(z(eindx),t_old(eindx),'k*',z2,z2./max(z2),'r-');
    xlabel('z-coordinate'); ylabel('Values'); title('Edges'); 
    legend('comp','exact','Location','northwest');
    subplot(2,2,3); plot(z(findx),t_old(findx),'k*',z2,z2./max(z2),'r-');
    xlabel('z-coordinate'); ylabel('Values'); title('Faces'); 
    legend('comp','exact','Location','northwest');
    subplot(2,2,4); plot(z(vindx),t_old(vindx),'k*',z2,z2./max(z2),'r-');
    xlabel('z-coordinate'); ylabel('Values'); title('Volumes'); 
    legend('comp','exact','Location','northwest');

    % compute fluxes (effective diffusivity)
    flux = abs(bs2 * k * fBO_1 * t_old);
    zflux = zsf_e.*flux; oflux = osf_e.*flux;
    
    disp([sum(zflux),sum(oflux),abs(sum(zflux)-sum(oflux))]./(max_z-min_z))

    pause(1)
end

% -------------------------------------------------------------------------
% extra plots
% -------------------------------------------------------------------------
% flux 2D
figure;
flux = abs(bs2 * k * fBO_1 * t_old);
subplot(3,1,1); plot(z2(peindx),flux(peindx),'k*-');
xlabel('z-coordinate'); ylabel('Fluxes'); title('Vertices-Edges');
subplot(3,1,2); plot(z2(efindx),flux(efindx),'k*-');
xlabel('z-coordinate'); ylabel('Fluxes'); title('Edges-Faces');
subplot(3,1,3); plot(z2(fvindx),flux(fvindx),'k*-');
xlabel('z-coordinate'); ylabel('Fluxes'); title('Faces-Volumes');

% solution 3D
figure;
i = 1; j = vcmplx(1).num(1).val; subplot(2,2,1); set(gca,'FontSize',16); axis([min_x max_x min_y max_y min_z max_z]);
plot3c2(fcmplx(1).cc(i:j,1),fcmplx(1).cc(i:j,2),fcmplx(1).cc(i:j,3),t_old(i:j),'o',' ',0,1); title('Vertices');
i = j+1; j = i-1+vcmplx(2).num(2).val; subplot(2,2,2); set(gca,'FontSize',16); axis([min_x max_x min_y max_y min_z max_z]);
plot3c2(fcmplx(1).cc(i:j,1),fcmplx(1).cc(i:j,2),fcmplx(1).cc(i:j,3),t_old(i:j),'o',' ',0,1); title('Edges');
i = j+1; j = i-1+vcmplx(3).num(3).val; subplot(2,2,3); set(gca,'FontSize',16); axis([min_x max_x min_y max_y min_z max_z]);
plot3c2(fcmplx(1).cc(i:j,1),fcmplx(1).cc(i:j,2),fcmplx(1).cc(i:j,3),t_old(i:j),'o',' ',0,1); title('Faces');
i = j+1; j = i-1+vcmplx(4).num(4).val; subplot(2,2,4); set(gca,'FontSize',16); axis([min_x max_x min_y max_y min_z max_z]);
plot3c2(fcmplx(1).cc(i:j,1),fcmplx(1).cc(i:j,2),fcmplx(1).cc(i:j,3),t_old(i:j),'o',' ',0,1); title('Volumes');





