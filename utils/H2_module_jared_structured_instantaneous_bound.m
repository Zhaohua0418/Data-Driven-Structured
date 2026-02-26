classdef H2_module_jared_structured_instantaneous_bound
    properties
        A,B,C,D,G,nx,nu,ny,nd,pattern
    end
    
    methods
        function obj = H2_module_jared_structured_instantaneous_bound(sys, pattern)
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

            cvx_begin sdp
            variables P(obj.nx,obj.nx) L(obj.nu,obj.nx) Q(obj.ny,obj.ny) R(obj.nx,obj.nx) tao(1,info.T) b
            minimize trace(Q)
            subject to
                            sum_over = 0;
                            for t = 1:info.T
                                left = [eye(obj.nx)            info.X(:,t+1);...
                                      zeros(obj.nx)           -info.X(:,t);...
                                      zeros(obj.nu,obj.nx)    -info.U(:,t);...
                                      zeros(obj.nx)            zeros(obj.nx,1)];
                                mid = blkdiag(info.epsilon^2*eye(obj.nx),-1);
                                sum_over = sum_over + tao(1,t) * left * mid * left';
                            end

                            [P-obj.G*obj.G'-b*eye(obj.nx)   zeros(obj.nx,obj.nx)    zeros(obj.nx,obj.nu)    zeros(obj.nx,obj.nx);...
                             zeros(obj.nx,obj.nx)'          zeros(obj.nx,obj.nx)    zeros(obj.nx,obj.nu)    R;...
                             zeros(obj.nx,obj.nu)'          zeros(obj.nx,obj.nu)'   zeros(obj.nu,obj.nu)    L;...
                             zeros(obj.nx,obj.nx)'          R'                      L'                      R+R'-P] - sum_over>=0;
            
                            [Q                      obj.C*R + obj.D*L;...
                             (obj.C*R + obj.D*L)'   R+R'-P]>=0;
           
                            tao>=0; b>0;
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
            out.tao = tao;
            out.b = b;
            out.h2 = h2;
        end
    end
end

