
function [q] = Matrix2Quat( Q )
    q = zeros(4,1); 
    q(4,1) = 0.5 * sqrt( 1 + Q(1,1) + Q(2,2) + Q(3,3) );

    if (q(4,1) ~= 0)
        Z = [Q(2,3) - Q(3,2); Q(3,1) - Q(1,3); Q(1,2) - Q(2,1)];
        q(1:3,1) = Z / (4 * q(4,1));
    else
        den = sqrt(  Q(1,2)^2 * Q(1,3)^2 + Q(1,2)^2 * Q(2,3)^2 + Q(1,3)^2 * Q(2,3)^2 );

        if (den ~= 0)
            q(1,1) = Q(1,3) * Q(1,2) / den;
            q(2,1) = Q(1,2) * Q(2,3) / den;
            q(3,1) = Q(1,3) * Q(2,3) / den;
        else
            a = 1;
        end
          
    end
end