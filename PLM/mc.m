% -------------------------------------------------%
% MATLAB code for mc slope limiter                 %
%                                                  %
% Written by Prof. Dongwook Lee                    %
% AMSC, UCSC                                       %
% -------------------------------------------------%

function delta = mc(a,b)
    delta = (sign(a)+sign(b))*min(min(abs(a),abs(b)),.25*abs(a+b));
end