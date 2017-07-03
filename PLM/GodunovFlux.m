% -------------------------------------------------%
% MATLAB code for Godunov-type fluxes              %
% Written by Dongwook Lee                          %
% -------------------------------------------------%

function [Flux]=GodunovFlux(s,uL,uR,PDEtype,q2)

    %Flux = 0;
    if (PDEtype == 1)
        % Upwind Godunov flux for constant advection
        if (s > 0)
            Flux = s*uL;
        else
            Flux = s*uR;
        end
    elseif (PDEtype == 2)
        % This is upwind Godunov flux for Burgers Eqn.
        if (uL >= uR)
            % Shock case
            if (s >= 0)
                Flux = 0.5*uL^2;
            elseif (s < 0)
                Flux = 0.5*uR^2;
            end
        else
            % Rarefaction case
            if (uR < 0.)
                Flux = 0.5*uR^2;
            elseif (uL > 0.)
                Flux = 0.5*uL^2;
            elseif ((uL < 0) & (uR > 0))
                Flux = 0;
            else
                Flux = 0;
            end
        end
    elseif (PDEtype == 3)
        % This is upwind Godunov flux for the Price model
        % Note that the flux function is concave!
        if (uL <= uR)
            % Shock case
            if (s >= 0)
                Flux = uL*(1+q2-2*q2*uL);
            elseif (s < 0)
                Flux = uR*(1+q2-2*q2*uR);
            end
        else %(uL>uR) 
            % Rarefaction case
            if (uL < 0.)
                Flux = uL*(1+q2-2*q2*uL);
            elseif (uR > 0.)
                Flux = uR*(1+q2-2*q2*uR);
            elseif ((uR < 0) & (uL > 0))
                Flux = 0;
            else
                Flux = 0;
            end
       end        
    end
end