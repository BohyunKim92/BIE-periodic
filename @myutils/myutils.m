classdef myutils <handle
    %MYUTILS Class for defining static methods used for solving layerpotential problem

    methods(Static)
        %library with special functions
         v1 = JacobiTheta1(z,tau);
         val = dlogtheta(z,tau);
         d = dS2(tau);
         Geval = G(z,tau);
         [Gx,Gy] = gradG(z,tau);
            
         %other helper functions
         [inddiag] = dind(A);
         [tval] = ztot(z,c);
         [value] = evalf(fs,target,indxs);
         indxs = coarse_to_fine_indxs(fineN,coarseN)
         [ns] = copy_segment(s)
    end
end