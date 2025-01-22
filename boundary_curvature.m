function omega = boundary_curvature(cmplx)
% -------------------------------------------------------------------------
% boundary_curvature.m
% -------------------------------------------------------------------------
% Purpose:      Compute the extended Delaunay/Voronoi Hodge star
%
% Institution:  King Fahd University of Petroleum and Minerals
%
% Author:       Dr Pieter Boom
% Date:         2025/01/07
% -------------------------------------------------------------------------

omega = zeros(cmplx(1).num(1).val,1);
for i=1:cmplx(1).num(1).val
    % boundary correction based on a curvature
    for j=1:cmplx(1).num(4).val(i)
        k = cmplx(1).bndop(4).indx(i,j); cnt = 1;
        v(1,1:3) = - cmplx(1).bc(i,:);
        v(2,1:3) = - cmplx(1).bc(i,:);
        v(3,1:3) = - cmplx(1).bc(i,:);
        for m=1:cmplx(4).num(2).val(k)
            n = cmplx(4).bndop(2).indx(k,m);
            if cmplx(2).bndop(1).indx(n,1)==i
                v(cnt,1:3) = v(cnt,1:3) + ...
                    cmplx(1).bc(cmplx(2).bndop(1).indx(n,2),:);
                cnt = cnt + 1;
            elseif cmplx(2).bndop(1).indx(n,2)==i
                v(cnt,1:3) = v(cnt,1:3) + ...
                    cmplx(1).bc(cmplx(2).bndop(1).indx(n,1),:);
                cnt = cnt + 1; end; end
        if cnt ~= 4; disp(cnt); pause; end
        a = v(1,1:3); b = v(2,1:3); c = v(3,1:3);
        tmp = 2*atan(abs(dot(a,cross(b,c))) / ...
            (norm(a)*norm(b)*norm(c) + ...
            dot(a,b)*norm(c) + dot(b,c)*norm(a) + dot(a,c)*norm(b)));
        if tmp<0; tmp = tmp + 2*pi; end
        omega(i) = omega(i) + tmp;
    end 
    % compute fvol
    omega(i) = 4*pi()/omega(i);
end