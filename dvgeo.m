function [star1,star2] = dvgeo(vcmplx)
% -------------------------------------------------------------------------
% dvgeo.m
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
star1 = [...
    vcmplx(1).dvol_mod;                        % Delaunay cell
    vcmplx(2).dvol_mod.*vcmplx(2).vvol./3;     % Delaunay face / Voronoi edge
    vcmplx(3).dvol_mod.*vcmplx(3).vvol./3;     % Delaunay edge / Voronoi face
    vcmplx(4).vvol];                           % Voronoi cell

% -------------------------------------------------------------------------
% star 2
% -------------------------------------------------------------------------
% - edge type 1 - Vornoi node -> edge
cnt = 0;
for i=1:vcmplx(2).num(2).val
    for j=1:vcmplx(2).num(1).val(i)
        cnt = cnt+1;
        k = vcmplx(2).bndop(1).indx(i,j);
        sgn = sign(dot(vcmplx(2).bndop(1).sgn(i,j)*...
            (vcmplx(1).cc(k,:)-vcmplx(2).cc(i,:)),vcmplx(2).dir(i,:)));
        
        % *** Uncommenting this code will improve stability of time
        %     integration, but will eliminate the interaction between
        %     non-well-centered cells in the extended complex
        % if sgn<0
        %     star2(cnt,1) = 0;
        % else
            star2(cnt,1) = sgn*vcmplx(2).dvol_mod(i)...      % Delaunary face area
                ./norm(vcmplx(1).cc(k,:)-vcmplx(2).cc(i,:)); % Vornoi node -> edge length
        % end
    end
end

% - edge type 2 - Vornoi edge -> face
for i=1:vcmplx(3).num(3).val
    for j=1:vcmplx(3).num(2).val(i)
        cnt = cnt+1;
        k = vcmplx(3).bndop(2).indx(i,j);
        sgn = sign(dot((vcmplx(2).cc(k,:)-vcmplx(3).cc(i,:)),...
            (vcmplx(2).cc(k,:)-vcmplx(3).bc(i,:)) ))./2;
        
        % *** Uncommenting this code will improve stability of time
        %     integration, but will eliminate the interaction between
        %     non-well-centered cells in the extended complex
        % if sgn<0
        %     star2(cnt,1) = 0;
        % else
            star2(cnt,1) = sgn*vcmplx(3).dvol_mod(i)*vcmplx(2).vvol(k)... % Mixed Delaunay edge / Voronoi edge area
                ./norm(vcmplx(2).cc(k,:)-vcmplx(3).cc(i,:));              % Vornoi edge -> face length
        % end
    end
end

% - edge type 3 - Vornoi face -> volume
for i=1:vcmplx(4).num(4).val
    for j=1:vcmplx(4).num(3).val(i)
        cnt = cnt+1;
        k = vcmplx(4).bndop(3).indx(i,j);
        sgn = sign(dot((vcmplx(3).cc(k,:)-vcmplx(4).cc(i,:)),...
            (vcmplx(3).cc(k,:)-vcmplx(4).bc(i,:)) ));
        
        % *** Uncommenting this code will improve stability of time
        %     integration, but will eliminate the interaction between
        %     non-well-centered cells in the extended complex
        % if sgn<0
        %     star2(cnt,1) = 0;
        % else
            star2(cnt,1) = sgn*vcmplx(3).vvol(k)...          % Voronoi face area
                ./norm(vcmplx(3).cc(k,:)-vcmplx(4).cc(i,:)); % Vornoi face -> volume length
        % end
    end
end

% remove exact zero values
star1(star1==0) = 10^-16;


