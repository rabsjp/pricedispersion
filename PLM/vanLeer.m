% -------------------------------------------------%
% MATLAB code for van Leer's slope limiter         %
%                                                  %
% Written by Prof. Dongwook Lee                    %
% AMSC, UCSC                                       %
% -------------------------------------------------%

function delta = vanLeer(a,b)
    if (a*b <= 0.)
        delta = 0.;
    else
        delta = 2.*a*b/(a+b);
    end
end