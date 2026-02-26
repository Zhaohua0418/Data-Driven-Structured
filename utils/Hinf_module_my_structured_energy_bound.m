classdef Hinf_module_my_structured_energy_bound
    properties
        A,B,G,C,D,H,nx,nu,ny,nd,pattern
    end
    
    methods
        function obj = Hinf_module_my_structured_energy_bound(sys,pattern)
            obj.A = sys.A;
            obj.B = sys.B;
            obj.G = sys.G;
            obj.C = sys.C;
            obj.D = sys.D;
            obj.H = sys.H;
            obj.nx = sys.nx;
            obj.nu = sys.nu;
            obj.ny = sys.ny;
            obj.nd = sys.nd;
            obj.pattern = pattern;
        end
        
        function out = model_based_Hinf(obj)
            lambda = 1;
            prev_P = eye(obj.nx);
            while 1
                cvx_begin sdp
                    variables P(obj.nx,obj.nx) K(obj.nu,obj.nx) Y(obj.nx,obj.nx) hinf_square Z(obj.nx,obj.nx)
                    minimize hinf_square+lambda*trace(Z)
                    subject to
                            [P                      obj.A+obj.B*K           obj.G                   zeros(obj.nx,obj.ny);...
                            (obj.A+obj.B*K)'        Y                       zeros(obj.nx,obj.nd)    (obj.C+obj.D*K)';...
                            obj.G'                  zeros(obj.nd,obj.nx)    eye(obj.nd)             obj.H';...
                            zeros(obj.ny,obj.nx)    obj.C+obj.D*K           obj.H                   hinf_square*eye(obj.ny)] >0;
                            Y <= inv(prev_P) - inv(prev_P)*(P-prev_P)*inv(prev_P) + Z;
                            K.*(1-obj.pattern) == 0;
                            Z >= 0;
                cvx_end
                hinf = sqrt(hinf_square);
                if norm(P-prev_P,"fro") < 0.001 ...
                        % && norm(P-inv(Q),"fro") < 0.01 
                    out = struct;
                    out.K = K;
                    out.P = P;
                    out.Y = Y;
                    out.Z = Z;
                    out.hinf = hinf;
                    break
                end
                prev_P = P;
                if lambda < 1e8
                    lambda = lambda * 2;
                end    
            end
        end

        function out = data_driven_Hinf(obj,info)
            lambda = 1;
            prev_P = eye(obj.nx);
            while 1
                cvx_begin sdp
                    variables K(obj.nu,obj.nx) P(obj.nx,obj.nx) L(obj.nu,obj.nx) R(obj.nu+obj.nx,obj.nu+obj.nx) Y(obj.nx,obj.nx) hinf_square Z(obj.nx,obj.nx) a b
                    minimize hinf_square+lambda*trace(Z)
                    subject to
                            SI = gen_si(info);

                            upper_left = [P-obj.G*obj.G'-b*eye(obj.nx)  zeros(obj.nx,obj.nu+obj.nx);...
                                          zeros(obj.nu+obj.nx,obj.nx)       -R];
                            upper_right = [obj.G*obj.H';...
                                           R*[obj.C';...
                                              obj.D']];
                            bottom_right = [hinf_square*eye(obj.ny) - [obj.C obj.D]*R*[obj.C';...
                                                                                       obj.D'] - obj.H*obj.H'];
                            [upper_left     upper_right;...
                             upper_right'   bottom_right] - a*blkdiag(SI,zeros(obj.ny))>=0;
    
                            upper_right_2 = [eye(obj.nx);...
                                             K];
                            [R                  upper_right_2;...
                             upper_right_2'     Y] >= 0;
                            Y <= inv(prev_P) - inv(prev_P)*(P-prev_P)*inv(prev_P) + Z;
                            K.*(1-obj.pattern) == 0;
                            a>=0; b>0;
                            P>=0; R>=0; Z >= 0;
                cvx_end
                hinf = sqrt(hinf_square);
                if norm(P)>1e5
                    out = NaN;
                    break
                end
                if norm(P-prev_P,"fro") < 0.01 ...
                     % &&   norm(P-inv(Q),"fro") < 0.01
                    out = struct;
                    out.K = K;
                    out.P = P;
                    out.L = L;
                    out.R = R;
                    out.Y = Y;
                    out.hinf = hinf;
                    out.Z = Z;
                    out.a = a;
                    out.b = b;
                    break
                end
                prev_P = P;
                if lambda < 1e8
                    lambda = lambda * 2;
                end    
            end
        end
    end
end

