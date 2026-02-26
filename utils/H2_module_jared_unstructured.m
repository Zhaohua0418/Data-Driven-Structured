classdef H2_module_jared_unstructured
    properties
        A,B,C,D,G,nx,nu,ny,nd
    end
    
    methods
        function obj = H2_module_jared_unstructured(sys)
            obj.A = sys.A;
            obj.B = sys.B;
            obj.C = sys.C;
            obj.D = sys.D;
            obj.G = sys.G;
            obj.nx = sys.nx;
            obj.nu = sys.nu;
            obj.ny = sys.ny;
            obj.nd = sys.nd;
        end
        
        function out = model_based_H2(obj)
           cvx_begin sdp
                    variables P(obj.nx,obj.nx) L(obj.nu,obj.nx) Q(obj.ny,obj.ny)
                    minimize trace(Q)
                    subject to

                            [Q                      obj.C*P + obj.D*L;...
                             P*obj.C' + L'*obj.D'    P] >= 0;

                            [P                      obj.A*P + obj.B*L       obj.G;...
                             P*obj.A' + L'*obj.B'   P                       zeros(obj.nx,obj.nd);...
                             obj.G'                 zeros(obj.nd,obj.nx)    eye(obj.nd, obj.nd)] >= 0;    
                            P>=0; Q>=0; 
            cvx_end
            h2 = sqrt(trace(Q));
            K = L*inv(P);
            out = struct;
            out.K = K;
            out.P = P;
            out.L = L;
            out.Q = Q;
            out.h2 = h2;
        end

        function out = data_driven_H2(obj,info)
            % mod_info = info_transform(info);
            Phi = blkdiag(info.T*info.epsilon^2*eye(obj.nx), -eye(info.T));
            data_matrix = [eye(obj.nx)             info.X(:,2:info.T+1);...
                   zeros(obj.nx,obj.nx)   -info.X(:,1:info.T);...
                   zeros(obj.nu,obj.nx)   -info.U(:,1:info.T)];
            SI = data_matrix*Phi*data_matrix';

            cvx_begin sdp
                    variables P(obj.nx,obj.nx) L(obj.nu,obj.nx) Q(obj.ny,obj.ny) a b
                    minimize trace(Q)
                    subject to
                            [P-obj.G*obj.G'-b*eye(obj.nx)   zeros(obj.nx,obj.nx)    zeros(obj.nx,obj.nu)    zeros(obj.nx,obj.nx);...
                             zeros(obj.nx,obj.nx)'          zeros(obj.nx,obj.nx)    zeros(obj.nx,obj.nu)    P;...
                             zeros(obj.nx,obj.nu)'          zeros(obj.nx,obj.nu)'   zeros(obj.nu,obj.nu)    L;...
                             zeros(obj.nx,obj.nx)'          P                       L'                      P] - a*blkdiag(SI,zeros(obj.nx,obj.nx))>=0;
            
                            [Q          obj.C*P+obj.D*L;...
                             (obj.C*P+obj.D*L)'   P]>=0;

                            a>=0; b>0;
                            P>=0; Q>=0; 
             cvx_end
             h2 = sqrt(trace(Q));
             K = L*inv(P);
             out = struct;
             out.P = P;
             out.L = L;
             out.K = K;
             out.Q = Q;
             out.a = a;
             out.b = b;
             out.h2 = h2;
        end
    end
end

