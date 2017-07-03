function [dt] = ComputeDt(maxU,dx,cfl);
   
    if (maxU == 0.) 
        disp('Error in dt calculation: max lambda = 0')
    else
        dt=cfl*dx/abs(maxU);
    end
end