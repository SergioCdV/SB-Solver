
function [d] = collision_constraint(r1, A1, v1, r2, A2, v2, iterations)
    % Compute vertex (convex hull) of both objects in the inertial space
    for i = 1:size(v1,2)
        v1(:,i) = A1 * v1(:,i);
    end

    for i = 1:size(v2,2)
        v2(:,i) = A2 * v2(:,i);
    end

    v1 = v1 + r1; 
    v2 = v2 + r2;

    % Compute the minimum distance through Gilbert-Johnson-Keerthi algorithm and the maximum distance through Expanding Polytope Algorithm
    [dp, dm, ~] = GJK(v1.', v2.', iterations);

    % Compute the Minkowski distance 
    d = norm(dm)-norm(dp);
end

%% Auxiliary functions 
function [a, b, flag] = GJK(shape1, shape2, iterations)
% GJK Gilbert-Johnson-Keerthi Collision detection implementation.
% Returns whether two convex shapes are are penetrating or not
% (true/false). Only works for CONVEX shapes.
%
% Inputs:
%   shape1: 
%   must have fields for XData,YData,ZData, which are the x,y,z 
%   coordinates of the vertices. Can be the same as what comes out of a 
%   PATCH object. It isn't required that the points form faces like patch
%   data. This algorithm will assume the convex hull of the x,y,z points
%   given.
%
%   shape2: 
%   Other shape to test collision against. Same info as shape1.
%   
%   iterations: 
%   The algorithm tries to construct a tetrahedron encompassing
%   the origin. This proves the objects have collided. If we fail within a
%   certain number of iterations, we give up and say the objects are not
%   penetrating. Low iterations means a higher chance of false-NEGATIVES
%   but faster computation. As the objects penetrate more, it takes fewer
%   iterations anyway, so low iterations is not a huge disadvantage.
%   
% Outputs:
%   flag:
%   true - objects collided
%   false - objects not collided
%
%  
%   This video helped me a lot when making this: https://mollyrocket.com/849
%   Not my video, but very useful.
%   
%   Matthew Sheen, 2016
%
    % Point 1 and 2 selection (line segment)
    v = [0.8 0.5 1];
    [a, b] = pickLine(v, shape2, shape1);

    % Point 3 selection (triangle)
    [a, b, c, ~] = pickTriangle(a, b, shape2, shape1, iterations);

    % Point 4 selection (tetrahedron)
    [a, b, c, d, flag] = pickTetrahedron(a, b, c, shape2, shape1, iterations);

    % EPA 
    if (flag)
        s = [a; b; c; d];                            % Initial simplex
        b = EPA(shape1, shape2, s, iterations);      % Maximum distance between the convex hulls
    else
        b = [];
    end
end

% Construct the first line of the simplex
function [a,b] = pickLine(v, shape1, shape2)
    b = support(shape2, shape1, v);
    a = support(shape2, shape1, -b);
end

function [a, b, c, flag] = pickTriangle(a ,b, shape1, shape2, IterationAllowed)
    % First try:
    ab = b-a;
    ao = -a;
    v = cross(cross(ab,ao),ab);  % v is perpendicular to ab pointing in the general direction of the origin.
    c = b;
    b = a;
    a = support(shape2, shape1, v);
    
    iter = 1;
    GoOn = true;
    while (GoOn && iter < IterationAllowed)
        % Time to check if we got it
        ab = b-a;
        ao = -a;
        ac = c-a;
        
        % Normal to face of triangle
        abc = cross(ab,ac);
        
        % Perpendicular to AB going away from triangle
        abp = cross(ab,abc);
        % Perpendicular to AC going away from triangle
        acp = cross(abc,ac);
        
        % First, make sure our triangle "contains" the origin in a 2d projection sense.
        % Is origin above (outside) AB?   
        if dot(abp,ao) > 0
            c = b;              % Throw away the furthest point and grab a new one in the right direction
            b = a;
            v = abp;            % Cross(Cross(ab,ao),ab);
        % Is origin above (outside) AC?
        elseif dot(acp, ao) > 0
            b = a;
            v = acp;            % Cross(Cross(ac,ao),ac);
        else
            % GoOn = false;
        end
        a = support(shape2,shape1,v);
        iter = iter + 1;
    end

    flag = ~GoOn;
end

function [a,b,c,d,flag] = pickTetrahedron(a,b,c,shape1,shape2,IterationAllowed)
    % Now, if we're here, we have a successful 2D simplex, and we need to check if the origin is inside a successful 3D simplex. So, is the origin above or below the triangle?
    flag = 0;
    ab = b-a;
    ac = c-a;

    % Normal to face of triangle
    abc = cross(ab,ac);
    ao = -a;

    if dot(abc, ao) > 0 % Above
        d = c;
        c = b;
        b = a;
        
        v = abc;
        a_new = support(shape2,shape1,v); % Tetrahedron new point

        if ~ismember(a_new', [a;b;c;d]')
            d = c;
            c = b;
            b = a;    
            a = a_new;
        end
        
    else % below
        d = b;
        b = a;
        v = -abc;
        a_new = support(shape2,shape1,v); % Tetrahedron new point

        if ~ismember(a_new', [a;b;c;d]')
            d = b;
            b = a;    
            a = a_new;
        end
    end

    iter = 1; 
    GoOn = true;
    while (GoOn && iter < IterationAllowed)
        % Check the tetrahedron
        ab = b-a;
        ao = -a;
        ac = c-a;
        ad = d-a;
        
        %We KNOW that the origin is not under the base of the tetrahedron based on
        %the way we picked a. So we need to check faces ABC, ABD, and ACD.
        
        % Normal to face of triangle
        abc = cross(ab,ac);
        
        % Above triangle ABC
        if dot(abc, ao) > 0 
            % No need to change anything, we'll just iterate again with this face as default
        else
            acd = cross(ac,ad);         % Normal to face of triangle
            
            % Above triangle ACD
            if dot(acd, ao) > 0 
                % Make this the new base triangle
                b = c;
                c = d;
                ab = ac;
                ac = ad;            
                abc = acd;     

            elseif dot(acd, ao) < 0
                adb = cross(ad,ab);     % Normal to face of triangle
                
                if dot(adb, ao) > 0     % Above triangle ADB
                    % Make this the new base triangle.
                    c = b;
                    b = d;              
                    ac = ab;
                    ab = ad;
                    abc = adb; 
                else
                    GoOn = false; 
                end
            else
                % The nearest point is one of the three points contained in the simplex
            end
        end
        
        % Try again, above
        if dot(abc, ao) > 0    
            v = abc;
            a_new = support(shape2,shape1,v);   % Tetrahedron new point

            if ~ismember(a_new', [a;b;c;d]')
                d = c;
                c = b;
                b = a;    
                a = a_new;
            end

        % Below
        else 
            v = -abc;
            a_new = support(shape2,shape1,v);   % Tetrahedron new point
            if ~ismember(a_new', [a;b;c;d]')
                d = b;
                b = a;    
                a = a_new;
            end
        end

        iter = iter + 1;
    end
end

function point = getFarthestInDir(shape, v)
    % Find the furthest point in a given direction for a shape
    XData = shape(:,1); % Making it more compatible with previous MATLAB releases.
    YData = shape(:,2);
    ZData = shape(:,3);
    dotted = XData*v(1) + YData*v(2) + ZData*v(3);
    [maxInCol, rowIdxSet] = max(dotted);
    [~, colIdx] = max(maxInCol);
    rowIdx = rowIdxSet(colIdx);
    point = [XData(rowIdx,colIdx), YData(rowIdx,colIdx), ZData(rowIdx,colIdx)];
end

% Support function to get the Minkowski difference
function [point] = support(shape1, shape2, v)
    point1 = getFarthestInDir(shape1, v);
    point2 = getFarthestInDir(shape2, -v);
    point = point1 - point2;
end

% Expanding Polytope Algorithm
function [dm] = EPA(shape1, shape2, s, maxIter)
    % Main iterations
    tol = 1e-4;                                 % Numerical stopping tolerance
    iter = 1;                                   % Initial iteration
    GoOn = true;                                % Convergence boolean
    dm = [];                                    % Maximum distance between the convex sets

    while (GoOn && iter < maxIter)
        [ds, e, index] = findClosestEdge(s);    % Find the closest edge
        p = support(shape1, shape2, e);         % Find the support function for the convex set
        d = dot(p,e);                           % Compute the distance

        % Convergence analysis
        if (d - ds < tol)
            dm = d * e;
            GoOn = false;
        else
            % Augment the simplex and the iterations
            s = [s(1:index,:); p; s(index + 1:end,:)];   
            iter = iter + 1;
        end
    end
end

function [dmin, nmin, index] = findClosestEdge(s)
    dmin = 1e20;
    nmin = [];
    index = [];
    for i = 1:size(s,1)
        for j = i + 1:size(s,1)
            a = s(i,:); 
            b = s(j,:);
            e = b-a;
            n = cross(e, cross(a,e));
            n = n / norm(n);
            d = dot(n,a);
            if (d < dmin)
                dmin = d;
                nmin = n; 
                index = j;
            end
        end
    end
end

function [reconstruct] = reconstruct_simplex(simplex, simplexFaces, extendedPoint)
    removalFaces = [];
    for i = 1:size(simplexFaces,1)
        face = simplexFaces(i,:);
        ab = face(1)-face(2);
        ac = face(1)-face(3); 
        n = cross(ab, ac);
        n = n/norm(n);
        a0 = [face(1) 0 0];
        if (dot(a0, n) > 0)
            n = -n;
        end

        if (dot(extendedPoint, n) > 0)
            removalFaces = [removalFaces; face];
        end
    end

    edges = []; 
    for i = 1:size(removalFaces,1)
        
    end
end