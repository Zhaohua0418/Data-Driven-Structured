classdef simulator
    properties
        sys,A,B,C,D,G,nx,nu,ny,nd,sampler;
    end
    
    methods
        function obj = simulator(sys)
            obj.sys = sys;
            obj.A = sys.A;
            obj.B = sys.B;
            obj.C = sys.C;
            obj.D = sys.D;
            obj.G = sys.G;
            obj.nx = sys.nx;
            obj.nu = sys.nu;
            obj.ny = sys.ny;
            obj.nd = sys.nd;
            obj.sampler = struct('u', @() 2*rand(sys.nu, 1)-1,...
                                 'w', @() ball_sample(1, sys.nx)');
        end
        
        function out = run(obj, epsilon, T)
            x0 = zeros(obj.nx,1);
            xcurr = x0;
            X = [x0 zeros(obj.nx,T)];
            U = zeros(obj.nu, T);
            W = zeros(obj.nx, T);
            for i = 1:T
                %inputs
                ucurr = obj.sampler.u();
                
                %propagation 
                %this is where the noise enters (process ?)
                wcurr = obj.sampler.w()*epsilon;
                xnext = obj.A*xcurr + obj.B*ucurr + wcurr; 
               
                
                %storage
                X(:, i+1) = xnext;
                U(:, i) = ucurr;
                W(:, i) = wcurr;
                xcurr = xnext;
            end
            out = struct;
            out.Xp = X(:,2:T+1);
            out.Xm = X(:,1:T);
            out.Um = U;
            out.T = T;
            out.epsilon = epsilon;
            out.sys = obj.sys;
        end
    end
end

