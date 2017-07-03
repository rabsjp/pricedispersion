% -------------------------------------------------%
% MATLAB code for 3 types of boundary conditions   %
%     (1) outflow,                                 %
%     (2) periodic,                                %
%     (3) fixed                                    %
% Written by Dongwook Lee and JP Rabanal           %
% AMSC, UCSC                                       %
% -------------------------------------------------%

function [u]=applyBC(u,ibeg,iend,ngc,BCtype)

    if BCtype == 1
        %outflow BC
        u(1:ngc          )=u(ibeg);
        u(iend+1:iend+ngc)=u(iend);

    elseif BCtype == 2
        %periodic BC
        u(ibeg-ngc  )=u(iend-1);
        u(ibeg-ngc+1)=u(iend);
        u(iend+1    )=u(ibeg);
        u(iend+2    )=u(ibeg+1);
        
    elseif BCtype == 3
        %fixed BC
        u(1:ngc          )=0;
        u(iend+1:iend+ngc)=1;        
    end

end