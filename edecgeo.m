function [star1,star2] = edecgeo...
    (cmplx,node_vol,edge_vol,face_vol,vol_vol,...
    node_edge_area,edge_face_area,face_vol_area)
% -------------------------------------------------------------------------
% edecgeo.m
% -------------------------------------------------------------------------
% Purpose:      Compute the extended Delaunay/Voronoi Hodge star
%
% Institution:  King Fahd University of Petroleum and Minerals
%
% Author:       Dr Pieter Boom
% Date:         2025/01/07
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% star 1
% -------------------------------------------------------------------------
star1 = [node_vol; edge_vol; face_vol; vol_vol];

% -------------------------------------------------------------------------
% star 2
% -------------------------------------------------------------------------
% - edge type 1 - Vornoi node -> edge
cnt = 0;
for i=1:cmplx(2).num(2).val
    for j=1:cmplx(2).num(1).val(i)
        cnt = cnt+1;
        k = cmplx(2).bndop(1).indx(i,j);
        sgn = sign(dot(cmplx(2).bndop(1).sgn(i,j)*...
            (cmplx(1).cc(k,:)-cmplx(2).cc(i,:)),cmplx(2).dir(i,:)));
        
        % Voronoi edge area divided by node-edge Forman edge length
        star2(cnt,1) = sgn*node_edge_area(i)./norm(...
            cmplx(1).cc(k,:)-cmplx(2).cc(i,:));
    end
end

% - edge type 2 - Vornoi edge -> face
for i=1:cmplx(3).num(3).val
    for j=1:cmplx(3).num(2).val(i)
        cnt = cnt+1;
        k = cmplx(3).bndop(2).indx(i,j);
        sgn = sign(dot((cmplx(2).cc(k,:)-cmplx(3).cc(i,:)),...
            (cmplx(2).cc(k,:)-cmplx(3).bc(i,:)) ));
        
        % Voronoi (face thickness x edge length) divided by edge-face 
        % Forman edge length
        star2(cnt,1) = sgn*edge_face_area(i,k)./norm(...
            cmplx(2).cc(k,:)-cmplx(3).cc(i,:));
    end
end

% - edge type 3 - Vornoi face -> volume
for i=1:cmplx(4).num(4).val
    for j=1:cmplx(4).num(3).val(i)
        cnt = cnt+1;
        k = cmplx(4).bndop(3).indx(i,j);
        sgn = sign(dot((cmplx(3).cc(k,:)-cmplx(4).cc(i,:)),...
            (cmplx(3).cc(k,:)-cmplx(4).bc(i,:)) ));
        
        % Voronoi face area divided by face-vol Forman edge length
        star2(cnt,1) = sgn*face_vol_area(k)./norm(...
            cmplx(3).cc(k,:)-cmplx(4).cc(i,:));
    end
end





