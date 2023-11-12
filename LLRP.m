% butterfly reproduce from paper
%
classdef LLRP < handle
    properties
        % constants
        hbar = 1.05e-34;
        e = 1.6e-19;
        a0 = 2.46e-10;
        sigma0 = [[1, 0]; [0, 1]];
        sigma1 = [[0, 1]; [1, 0]];
        sigma2 = [[0, -1i]; [1i, 0]];
        sigma3 = [[1, 0]; [0, -1]];
        sigma_plus;
        sigma_minus;
           
        % input parameters
        thetad;
        w;
        eta;
        vF;
        N_ll;

        % calculated variables
        thetar;
        thetar1;
        thetar2;
        KD;
        ktheta;
        w0;
        w1;
        zeta;
        T1;
        T2;
        T3;
        q1;
        q2;
        q3;

        % magnetic field related variable
        p;
        q;
        B;
        omegac;
        lB;
        Delta;
        H_size;
        layer_size;
        subband_size;
        atom_size;
        laguerre_matrix;
        f_matrix1;
        f_matrix2;
        f_matrix3;

        % after fix p, q, k1, k2
        H;
        H_info;

    end

    methods
        function obj = LLRP(nvargs)
            arguments
                nvargs.thetad = 2.5;
                nvargs.w = 110e-3 * 1.6e-19;
                nvargs.eta = 0.9;
                nvargs.vF = 1e6;
                nvargs.N_ll = 4;
            end
            obj.init(theta=nvargs.thetad, w=nvargs.w, eta=nvargs.eta,...
                vF=nvargs.vF,...
                N_ll=nvargs.N_ll);
        end

        function init(obj, nvargs)
            arguments
                obj;
                nvargs.thetad;
                nvargs.w;
                nvargs.eta;
                nvargs.vF;
                nvargs.N_ll;
            end
            if isfield(nvargs, 'thetad')
                obj.thetad = nvargs.thetad;
            end
            if isfield(nvargs, 'w')
                obj.w = nvargs.w;
            end
            if isfield(nvargs, 'eta')
                obj.eta = nvargs.eta;
            end
            if isfield(nvargs, 'vF')
                obj.vF = nvargs.vF;
            end
            if isfield(nvargs, 'N_ll')
                obj.N_ll = nvargs.N_ll;
            end

            obj.thetar = obj.thetad * pi / 180;
            obj.thetar1 = -obj.thetar/2;
            obj.thetar2 = obj.thetar/2;
            obj.KD = 4 * pi / (3 * obj.a0);
            obj.ktheta = 2 * obj.KD * sin(obj.thetar/2);
            obj.w0 = obj.eta * obj.w;
            obj.w1 = obj.w;
            obj.sigma_plus = (obj.sigma1 + 1i * obj.sigma2) / 2;
            obj.sigma_minus = (obj.sigma1 - 1i * obj.sigma2) / 2;
            obj.zeta = exp(1i * 2 * pi / 3);
            obj.T1 = obj.w0 * obj.sigma0 + obj.w1 * obj.sigma1;
            obj.T2 = obj.w0 * obj.zeta' * obj.sigma0 + obj.w1 * obj.sigma_plus +...
                obj.w1 * obj.zeta * obj.sigma_minus;
            obj.T3 = obj.w0 * obj.zeta * obj.sigma0 + obj.w1 * obj.sigma_plus +...
                obj.w1 * obj.zeta' * obj.sigma_minus;
            obj.q1 = obj.ktheta * [0, -1];
            obj.q2 = obj.ktheta * [sqrt(3), 1] / 2;
            obj.q3 = obj.ktheta * [-sqrt(3), 1] / 2;
        end

        function precompute_init(obj, p, q)
            obj.p = p;
            obj.q = q;
            obj.B = sqrt(3) / 4 * obj.ktheta^2 * obj.hbar / (obj.e * 2 * pi) * q / p;
            obj.omegac = sqrt(2 * obj.hbar * obj.e * obj.vF^2 * obj.B);
            obj.lB = sqrt(obj.hbar / (obj.e * obj.B));
            obj.Delta = sqrt(3) * obj.ktheta * obj.lB^2 / 2;
            obj.H_size = 2 * (2 * obj.N_ll + 1) * obj.q;
            obj.layer_size = (2 * obj.N_ll + 1) * obj.q;
            obj.subband_size = (2 * obj.N_ll + 1);
            obj.atom_size = 2;
            obj.H = zeros(obj.H_size);
            obj.generate_f_coeff();
        end

        function construct_H(obj, k1, k2, info_mode)
            arguments   
                obj;
                k1;
                k2;
                info_mode = false;
            end
            if info_mode
                obj.H_info = string(obj.H);
            end

            % iterate all matrix element (row, col)
            % now the level is layer -> q subbands -> AB -> landau n
            for row_d = 1:1:obj.H_size
                for col_d = 1:1:obj.H_size
                    row = row_d-1;
                    col = col_d-1;
                    % level label 0 or 1
                    l = floor(row/obj.layer_size);
                    lp = floor(col/obj.layer_size);
                    % subband label 0, 1, ... , q-1
                    j = floor(mod(row, obj.layer_size)/obj.subband_size); 
                    jp = floor(mod(col, obj.layer_size)/obj.subband_size);
                    % atom label 0, 1
                    alpha = -1;
                    alphap = -1;
                    n = -1;
                    np = -1;
                    if mod(mod(row, obj.layer_size), obj.subband_size) <= obj.N_ll
                        alpha = 0;
                        n = mod(mod(row, obj.layer_size), obj.subband_size);
                    else
                        alpha = 1;
                        n = mod(mod(row, obj.layer_size), obj.subband_size) -...
                            (obj.N_ll + 1);
                    end
                    if mod(mod(col, obj.layer_size), obj.subband_size) <= obj.N_ll
                        alphap = 0;
                        np = mod(mod(col, obj.layer_size), obj.subband_size);
                    else
                        alphap = 1;
                        np = mod(mod(col, obj.layer_size), obj.subband_size) -...
                            (obj.N_ll + 1);
                    end
                    
                    % TODO
                    % should consider delete the info mode later
                    % if info_mode
                    %     obj.H_info(row_d, col_d) = sprintf("<%d%d%d%d, %d%d%d%d>",...
                    %         l, j, alpha, n, lp, jp, alphap, np);
                    % end
                    
                    % give values
                    % intralayer
                    if l == lp && j == jp
                        if l == 0
                            cur_theta = obj.thetar1;
                        else
                            cur_theta = obj.thetar2;
                        end
                        if n == np + 1 && alpha == 0 && alphap == 1
                            obj.H(row_d, col_d) = -obj.omegac *...
                                exp(-1i*cur_theta)*...
                                sqrt(n); 
                            % if info_mode
                            %     obj.H_info(row_d, col_d) = obj.H_info(row_d, col_d)+...
                            %         sprintf('sq(%d)', n);
                            % end
                        end
                        if n == np - 1 && alpha == 1 && alphap == 0
                            obj.H(row_d, col_d) = -obj.omegac *...
                                exp(1i*cur_theta)*...
                                sqrt(np); 
                            % if info_mode
                            %     obj.H_info(row_d, col_d) = obj.H_info(row_d, col_d)+...
                            %         sprintf('sq(%d)', np);
                            % end
                        end
                    end
                    % interlayer
                    if l ~= lp && l == 1
                        cur_ele = 0.0;
                        cur_info = '';
                        if j == mod(jp, obj.q)
                            cur_ele = cur_ele + obj.T1(alpha+1, alphap+1) *...
                            obj.f_matrix1(n+1, np+1) *...
                            exp(-1i * k1 * obj.ktheta * obj.lB^2) *...
                            exp(-4 * pi * 1i * obj.p / obj.q * j);
                            % if info_mode
                            %     cur_info = append(cur_info, 't1');
                            % end
                        end
                        if j == mod(jp + 1, obj.q)
                            cur_ele = cur_ele + obj.T2(alpha+1, alphap+1) *...
                            obj.f_matrix2(n+1, np+1) *...
                            exp(1i * k2 * obj.Delta) *...
                            exp(1i / 2 * k1 * obj.ktheta * obj.lB^2) *...
                            exp(1i * pi * obj.p / obj.q * (2 * jp - 1));
                            % if info_mode
                            %     cur_info = append(cur_info, 't2');
                            % end
                        end
                        if j == mod(jp - 1, obj.q)
                            cur_ele = cur_ele + obj.T3(alpha+1, alphap+1) *...
                            obj.f_matrix3(n+1, np+1) *...
                            exp(-1i * k2 * obj.Delta) *...
                            exp(1i / 2 * k1 * obj.ktheta * obj.lB^2) *...
                            exp(1i * pi * obj.p / obj.q * (2 * jp + 1));
                            % if info_mode
                            %     cur_info = append(cur_info, 't3');
                            % end
                        end
                        
                        obj.H(row_d, col_d) = cur_ele;
                        % if info_mode
                        %     obj.H_info(row_d, col_d) = append(obj.H_info(row_d, col_d), cur_info);
                        % end
                    end
                end
            end
            % for the top right block of H, use T' directly
            obj.H(1:obj.layer_size, (obj.layer_size+1):obj.H_size) =...
                obj.H((obj.layer_size+1):obj.H_size, 1:obj.layer_size)';
            % if info_mode
            %     obj.H_info(1:obj.layer_size, (obj.layer_size+1):obj.H_size) =...
            %         obj.H_info((obj.layer_size+1):obj.H_size, 1:obj.layer_size)';
            % end
        end

        function generate_f_coeff(obj)
            z1 = (obj.q1(1) + 1i * obj.q1(2)) / sqrt(2) * obj.lB;
            z2 = (obj.q2(1) + 1i * obj.q2(2)) / sqrt(2) * obj.lB;
            z3 = (obj.q3(1) + 1i * obj.q3(2)) / sqrt(2) * obj.lB;
            % all three z have the same norm
            z_square = norm(z1)^2;
            for n = 1:(obj.N_ll + 1)
                for m = 1:(obj.N_ll + 1)
                    obj.laguerre_matrix(n, m) = laguerreL(n-1, m-1, z_square);
                end
            end
            for n = 1:(obj.N_ll+1)
                for m = 1:n
                    common_factor = sqrt(factorial(m-1) / factorial(n-1)) * ...
                        exp(-z_square/2);
                    obj.f_matrix1(n, m) = common_factor * (-z1')^(n-m) *...
                        obj.laguerre_matrix(m, n-m+1);
                    obj.f_matrix2(n, m) = common_factor * (-z2')^(n-m) *...
                        obj.laguerre_matrix(m, n-m+1);
                    obj.f_matrix3(n, m) = common_factor * (-z3')^(n-m) *...
                        obj.laguerre_matrix(m, n-m+1);
                end
                for m = (n+1):(obj.N_ll+1)
                    common_factor = sqrt(factorial(n-1) / factorial(m-1)) * ...
                        exp(-z_square/2);
                    obj.f_matrix1(n, m) = common_factor * (z1')^(m-n) *...
                        obj.laguerre_matrix(n, m-n+1);
                    obj.f_matrix2(n, m) = common_factor * (z2')^(m-n) *...
                        obj.laguerre_matrix(n, m-n+1);
                    obj.f_matrix3(n, m) = common_factor * (z3')^(m-n) *...
                        obj.laguerre_matrix(n, m-n+1);
                end
            end
        end
    end
end