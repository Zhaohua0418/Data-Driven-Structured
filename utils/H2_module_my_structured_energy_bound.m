classdef H2_module_my_structured_energy_bound
    properties
        A,B,C,D,G,nx,nu,ny,nd,pattern
    end
    
    methods
        function obj = H2_module_my_structured_energy_bound(sys, pattern)
            obj.A = sys.A;
            obj.B = sys.B;
            obj.C = sys.C;
            obj.D = sys.D;
            obj.G = sys.G;
            obj.nx = sys.nx;
            obj.nu = sys.nu;
            obj.ny = sys.ny;
            obj.nd = sys.nd;
            obj.pattern = pattern;
        end
        
        function out = model_based_H2(obj)
            lambda = 1;
            % lambda = 1e7;
            prev_X = eye(obj.nx);
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

        function out = data_driven_H2(obj,info)
            % mod_info = info_transform(info);
            lambda = 1;
            prev_X = eye(obj.nx);
            while 1
                cvx_begin sdp
                    variables X(obj.nx,obj.nx) Y(obj.nx,obj.nx) K(obj.nu,obj.nx) W(obj.ny,obj.ny) Z_plus(obj.nx,obj.nx) a b
                    minimize trace(W) + lambda * trace(Z_plus)
                    subject to
                            Phi = blkdiag(info.T*info.epsilon^2*eye(obj.nx), -eye(info.T));
                            data_matrix = [eye(obj.nx)             info.X(:,2:info.T+1);...
                                   zeros(obj.nx,obj.nx)   -info.X(:,1:info.T);...
                                   zeros(obj.nu,obj.nx)   -info.U(:,1:info.T)];
                            SI = data_matrix*Phi*data_matrix';


                            [X-obj.G*obj.G'-b*eye(obj.nx)   zeros(obj.nx,obj.nx)    zeros(obj.nx,obj.nu)    zeros(obj.nx,obj.nx);...
                             zeros(obj.nx,obj.nx)'          zeros(obj.nx,obj.nx)    zeros(obj.nx,obj.nu)    eye(obj.nx,obj.nx);...
                             zeros(obj.nx,obj.nu)'          zeros(obj.nx,obj.nu)'   zeros(obj.nu,obj.nu)    K;...
                             zeros(obj.nx,obj.nx)'          eye(obj.nx,obj.nx)      K'                      Y] - a*blkdiag(SI,zeros(obj.nx,obj.nx))>=0;
            
                            [W          obj.C+obj.D*K;...
                             (obj.C+obj.D*K)'   Y]>=0;
            
                            Y <= inv(prev_X) - inv(prev_X)*(X-prev_X)*inv(prev_X) + Z_plus;
            
                            a>=0; b>0;
                            X>=0; W>=0; Z_plus >= 0;
                            
                            K.*(1-obj.pattern) == 0;
                cvx_end
                h2 = sqrt(trace(W));
                if norm(X-prev_X,"fro") < 0.01 && norm(X-inv(Y),"fro") < 0.01
                     out = struct;
                     out.X = X;
                     out.Y = Y;
                     out.K = K;
                     out.W = W;
                     out.Z_plus = Z_plus;
                     out.a = a;
                     out.b = b;
                     out.h2 = h2;
                     break
                end
                lambda = lambda * 2;
                prev_X = X;
            end
        end
    end
end

