classdef model_based_H2
    properties
        A,B,C,D,G,pattern,nx,nu,ny,nd;
    end
    
    methods
        function obj = model_based_H2(sys, pattern)
            obj.A = sys.A;
            obj.B = sys.B;
            obj.C = sys.C;
            obj.D = sys.D;
            obj.G = sys.G;
            obj.pattern = pattern;
            obj.nx = sys.nx;
            obj.nu = sys.nu;
            obj.ny = sys.ny;
            obj.nd = sys.nd;
        end
        
        function out = run(obj)
            lambda = 1;
            [prev_K, prev_X] = obj.find_optimal();
            while 1
                cvx_begin sdp
                variables K(obj.nu,obj.nx) X(obj.nx,obj.nx) Y(obj.nx,obj.nx) W(obj.ny,obj.ny) Z_plus(obj.nx,obj.nx)
                minimize trace(W) + lambda * trace(Z_plus)
                subject to
                         [W                     obj.C + obj.D*K;...
                          (obj.C + obj.D*K)'    Y] >= 0;

                         [X                  obj.A+obj.B*K          obj.G;...
                          (obj.A+obj.B*K)'   Y                      zeros(obj.nx,obj.nd);...
                          obj.G'             zeros(obj.nd,obj.nx)   eye(obj.nd, obj.nd)] >= 0;

                         Y <= inv(prev_X) - inv(prev_X)*(X-prev_X)*inv(prev_X) + Z_plus;
                            
                         Z_plus >= 0;

                         K.*(1-obj.pattern) == 0;
                cvx_end
                h2 = sqrt(trace(W));
                if norm(X-prev_X,"fro") < 0.01 && norm(X-inv(Y),"fro") < 0.01
                    out = struct;
                    out.K = K;
                    out.X = X;
                    out.Y = Y;
                    out.W = W;
                    out.Z_plus = Z_plus;
                    out.h2 = h2;
                    break
                end
                lambda = lambda * 2;
                prev_X = X;
            end
        end

        function [K,X] = find_optimal(obj)
            cvx_begin sdp
                    variables X(obj.nx,obj.nx) Y(obj.nu,obj.nx) W(obj.ny,obj.ny)
                    minimize trace(W)
                    subject to

                            [W                      obj.C*X + obj.D*Y;...
                             X*obj.C' + Y'*obj.D'    X] >= 0;

                            [X                      obj.A*X + obj.B*Y       obj.G;...
                             X*obj.A' + Y'*obj.B'   X                       zeros(obj.nx,obj.nd);...
                             obj.G'                 zeros(obj.nd,obj.nx)    eye(obj.nd, obj.nd)] >= 0;    
            cvx_end
            h2 = sqrt(trace(W));
            K = Y*inv(X);
        end
    end
end

