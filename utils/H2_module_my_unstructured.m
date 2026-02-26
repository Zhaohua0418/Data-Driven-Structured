classdef H2_module_my_unstructured
    properties
        A,B,C,D,G,nx,nu,ny,nd
    end
    
    methods
        function obj = H2_module_my_unstructured(sys)
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
            cvx_begin sdp
                    variables P(obj.nx,obj.nx) L(obj.nu,obj.nx) Q(obj.ny,obj.ny) tao(1,info.T) 
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



                            [P-obj.G*obj.G'   zeros(obj.nx,obj.nx)    zeros(obj.nx,obj.nu)    zeros(obj.nx,obj.nx);...
                             zeros(obj.nx,obj.nx)'          zeros(obj.nx,obj.nx)    zeros(obj.nx,obj.nu)    P;...
                             zeros(obj.nx,obj.nu)'          zeros(obj.nx,obj.nu)'   zeros(obj.nu,obj.nu)    L;...
                             zeros(obj.nx,obj.nx)'          P                       L'                      P] - sum_over>=0;

                            [Q          obj.C*P+obj.D*L;...
                             (obj.C*P+obj.D*L)'   P]>=0;

                            tao >= 0;
                            P>=0; Q>=0; 
             cvx_end
             h2 = sqrt(trace(Q));
             K = L*inv(P);
             out = struct;
             out.P = P;
             out.L = L;
             out.K = K;
             out.Q = Q;
             out.tao = tao;
             out.h2 = h2;
        end
    end
end

