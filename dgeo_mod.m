function cmplx = dgeo_mod(cmplx,face_thic,node_edge_area,node_vol)
% -------------------------------------------------------------------------
% dgeo_mod.m
% -------------------------------------------------------------------------
% Purpose:      Compute Modified Delaunay geometry
%
% Pre:          Voronoi mesh has been computed
%
% Institution:  University of Manchester
% Group:        Mechanics and Physics of Solids 
%
% Author:       Dr Pieter Boom
% Date:         2021/12/14
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% polyhedra - Delaunay vertices
% - cc (done previous)
% - dvol
% -------------------------------------------------------------------------
cmplx(4).dvol = ones(cmplx(4).num(4).val,1);

% -------------------------------------------------------------------------
% faces - Delaunay edges
% - dvol
% -------------------------------------------------------------------------
cmplx(3).dvol_mod = zeros(cmplx(3).num(3).val,1);

for i=1:cmplx(3).num(3).val
    % Delaunay edge length is the norm of the difference between end points
    pnt = cmplx(4).cc(cmplx(3).bndop(4).indx(i,1),:);
    pnt = cmplx(3).bndop(4).sgn(i,1)*pnt;
    if (cmplx(3).num(4).val(i)==2)
        pnt = pnt + cmplx(3).bndop(4).sgn(i,2)*...
            cmplx(4).cc(cmplx(3).bndop(4).indx(i,2),:);
    else; pnt = pnt - cmplx(3).bndop(4).sgn(i,1)*cmplx(3).cc(i,:);
    end

    cmplx(3).dvol_mod(i) = min(face_thic(i),0.8*norm(pnt));
end

% -------------------------------------------------------------------------
% edges - Delaunay faces
% - dvol
% -------------------------------------------------------------------------
cmplx(2).dvol_mod = zeros(cmplx(2).num(2).val,1);

for i=1:cmplx(2).num(2).val
    % area is the sum of triangles partitions (face,edge,vertex)
    % circumcenters
    if cmplx(2).dvol(i)~=0
        for j=1:cmplx(2).num(3).val(i)
            k = cmplx(2).bndop(3).indx(i,j);
            d = norm(cmplx(2).cc(i,:) - cmplx(3).cc(k,:));
            sgn = sign(dot(cmplx(2).cc(i,:) - cmplx(3).cc(k,:),...
                cmplx(2).cc(i,:) - cmplx(3).bc(k,:)));

            cmplx(2).dvol_mod(i) = cmplx(2).dvol_mod(i) + ...
                2/pi*(min(abs(cmplx(2).dvol(i)),node_edge_area(i))./...
                abs(cmplx(2).dvol(i)))^(1/2)*d*sgn*cmplx(3).dvol_mod(k);
        end
    end
end

% -------------------------------------------------------------------------
% vertices - Delaunay tetra/hexahedra
% - dvol
% -------------------------------------------------------------------------
cmplx(1).dvol_mod = zeros(cmplx(1).num(1).val,1);

for i=1:cmplx(1).num(1).val
    % volume is the sum of hexahedron partitions (tet,face,edge,vertex)
    % circumcenters
    if cmplx(1).dvol(i) ~= 0
        for j=1:cmplx(1).num(2).val(i)
            k = cmplx(1).bndop(2).indx(i,j);
            d = norm(cmplx(1).cc(i,:) - cmplx(2).cc(k,:));
            sgn = sign(dot(cmplx(1).cc(i,:) - cmplx(2).cc(k,:),...
                cmplx(1).cc(i,:) - cmplx(2).bc(k,:)));

            cmplx(1).dvol_mod(i) = cmplx(1).dvol_mod(i) + ...
                2/pi*(min(abs(cmplx(1).dvol(i)),node_vol(i))./...
                abs(cmplx(1).dvol(i)))^(1/3)*d*sgn*cmplx(2).dvol_mod(k);
        end
    else
        disp('----')
    end
end



