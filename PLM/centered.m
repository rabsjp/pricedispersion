% -------------------------------------------------%
% MATLAB code for centered-non TVD slope limiter   %
%                                                  %
% Written by Dongwook Lee and JP Rabanal           %
% -------------------------------------------------%

function delta = centered(a,b)
    delta = .5*(a+b);
end