% -------------------------------------------------%
% MATLAB code for minmod slope limiter             %
%                                                  %
% Written by Prof. Dongwook Lee                    %
% AMSC, UCSC                                       %
% -------------------------------------------------%

function delta = minmod(a,b)
    delta = 0.5*(sign(a)+sign(b))*min(abs(a),abs(b));
end