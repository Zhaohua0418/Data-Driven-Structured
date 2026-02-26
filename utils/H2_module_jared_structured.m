classdef H2_module_jared_structured
    properties
        A,B,C,D,G,nx,nu,ny,nd,pattern
    end
    
    methods
        function obj = H2_module_jared_structured(sys, pattern)
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
            cvx_begin sdp
                    variables P(obj.nx,obj.nx) L(obj.nu,obj.nx) Q(obj.ny,obj.ny) R(obj.nx,obj.nx)
                    minimize trace(Q)
                    subject to

                            [Q                      obj.C*R + obj.D*L;...
                             (obj.C*R + obj.D*L)'    R+R'-P] > 0;

                            [P-obj.G*obj.G'                      obj.A*R + obj.B*L ;...
                             (obj.A*R + obj.B*L)'                R+R'-P                 ] > 0;    
                            P>=0; Q>=0; 
                            L.*(1-obj.pattern) == 0;
                            Upsilon = sparse_pattern(obj.pattern);
                            R.*(1-Upsilon) == 0;
                            
            cvx_end
            h2 = sqrt(trace(Q));
            K = L*inv(R);
            out = struct;
            out.K = K;
            out.P = P;
            out.L = L;
            out.Q = Q;
            out.R = R;
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
            variables P(obj.nx,obj.nx) L(obj.nu,obj.nx) Q(obj.ny,obj.ny) R(obj.nx,obj.nx) a b
            minimize trace(Q)
            subject to
                            [P-obj.G*obj.G'-b*eye(obj.nx)   zeros(obj.nx,obj.nx)    zeros(obj.nx,obj.nu)    zeros(obj.nx,obj.nx);...
                             zeros(obj.nx,obj.nx)'          zeros(obj.nx,obj.nx)    zeros(obj.nx,obj.nu)    R;...
                             zeros(obj.nx,obj.nu)'          zeros(obj.nx,obj.nu)'   zeros(obj.nu,obj.nu)    L;...
                             zeros(obj.nx,obj.nx)'          R'                      L'                      R+R'-P] - a*blkdiag(SI,zeros(obj.nx,obj.nx))>=0;
            
                            [Q                      obj.C*R + obj.D*L;...
                             (obj.C*R + obj.D*L)'   R+R'-P]>=0;
           
                            a>=0; b>0;
                            P>=0; Q>=0; 
                            L.*(1-obj.pattern) == 0;
                            Upsilon = sparse_pattern(obj.pattern);
                            R.*(1-Upsilon) == 0;
            cvx_end
            h2 = sqrt(trace(Q));
            K = L*inv(R);
            out = struct;
            out.K = K;
            out.P = P;
            out.L = L;
            out.Q = Q;
            out.R = R;
            out.a = a;
            out.b = b;
            out.h2 = h2;
        end
    end
end

