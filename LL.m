classdef LL < handle
    properties
        % constants
        hbar = 1.05e-34;
        e = 1.6e-19;

        % lattice parameters
        a = 1.4e-10;
        w;
        vF = 0.866e6;

        % tuning parameters
        % angle in degree
        thetad;
        % angle in radius
        thetar;
        thetar1;
        thetar2;
        % magnetic field 
        B;
        % cutoff landau level numbers
        N_ll;
        % p and q
        p;
        q;
        % derived variables
        lB;
        K_vector;
        Kx;
        K;
        KD;
        Delta;
        % magnetic zone range
        k1_max;
        k2_max;
        % T matrix, T(1) = T_1
        Ts;
        % cached function object
        tilde_T;
        % Hamiltonian in magnetic field
        H_size;
        layer_size;
        subband_size;
        atom_size;
        H;
        % a string matrix to show the information of H element
        % only works when info mode == true
        H_info;
    end
    
    methods
        function obj = LL(nvargs)
            arguments
                nvargs.thetad = 2.5;
                nvargs.w = 110e-3 * 1.6e-19;
                nvargs.vF = 0.866e6;
                nvargs.N_ll = 4;
            end
            obj.init(theta=nvargs.thetad, w=nvargs.w, vF=nvargs.vF,...
                N_ll=nvargs.N_ll);
        end

        function init(obj, nvargs)
            % after thetad is determined from iniput, 
            % set up all variables depending on thetad.
            arguments  
                obj;
                nvargs.thetad;
                nvargs.w;
                nvargs.vF;
                nvargs.N_ll;
            end
            if isfield(nvargs, 'thetad')
                obj.thetad = nvargs.thetad;
            end
            if isfield(nvargs, 'w')
                obj.w = nvargs.w;
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
            obj.K_vector = 2*pi/obj.a * [1/3; sqrt(3)/9];
            obj.Kx = obj.K_vector(1);
            obj.K = norm(obj.K_vector);
            obj.KD = 2*obj.K*sin(obj.thetar/2);
            % original T matrix
            w_phase = exp(1i*2*pi/3);
            T1 = [[1, 1]; [1, 1]];
            T2 = [[1, w_phase]; [w_phase', 1]];
            T3 = [[1, w_phase']; [w_phase, 1]];
            obj.Ts = zeros(2,2,3);
            obj.Ts(:,:,1) = T1;
            obj.Ts(:,:,2) = T2;
            obj.Ts(:,:,3) = T3;
            obj.tilde_T = zeros(2,2,3);


        end
        
        function precompute_init(obj, p, q)
            % TODO
            % do we really need the coprime?
            % if not coprime, just fold the bands
            % if gcd(p, q) ~= 1
            %     pq = gcd(p, q);
            %     p = p / pq;
            %     q = q / pq;
            % end
            obj.p = p;
            obj.q = q;
            % determine B first, given p and q
            obj.B = 3*sqrt(3)/4 * obj.KD^2 * obj.hbar/obj.e * obj.q / ...
                (obj.p * 2*pi);
            obj.lB = sqrt(obj.hbar/(obj.e*obj.B));
            obj.Delta = sqrt(3)/2 * obj.lB^2 * obj.KD;
            % magnetic zone range
            obj.k1_max = obj.Delta / (obj.lB^2);
            obj.k2_max = 2*pi / (obj.q*obj.Delta);
            % for a known magnetic field, we can fix the cached tilde_T
            obj.tilde_T = memoize(@obj.t_coeff);
           
            obj.H_size = 2 * 2 * obj.N_ll * obj.q;
            obj.layer_size = 2 * obj.N_ll * obj.q;
            obj.subband_size = 2 * obj.N_ll;
            obj.atom_size = 2;
            obj.H = zeros(obj.H_size);
        end

        % plot B dependent landau fan
        function plot_band_fix_k1k2(obj, k1, k2)
            % for now, fix q = 1
            q = 1;
            p_range = 2:1:100;
            B_list = zeros(1, length(p_range));
            % need to do a precompute to update the q inside the object
            E_list = zeros(4*obj.N_ll*q, length(p_range));
            for ind = 1:1:length(p_range)
                p = p_range(ind);
                obj.precompute_init(p, q);
                obj.construct_H(k1, k2, false);
                B_list(ind) = obj.B;
                E_list(:, ind) = eig(obj.H);
            end
            fig1 = figure();
            plot(axes(fig1), B_list, E_list);
            fig2 = figure();
            plot(axes(fig2), sqrt(B_list), E_list);
        end

        function plot_band_fix_pq(obj, p, q)
            obj.precompute_init(p, q);
            % for now, fix a k1 point
            k1 = 0 * obj.k1_max/4;
            % give a k2 range
            num_k2_point = 100;
            k2_step = obj.k2_max / num_k2_point;
            k2_range = 0:k2_step:obj.k2_max;

            E_list = zeros(obj.H_size, length(k2_range));
            for ind = 1:1:length(k2_range)
                k2 = k2_range(ind);
                obj.construct_H(k1, k2, false);
                E_list(:, ind) = eig(obj.H);
            end
            fig1 = figure();
            ax1 = axes(fig1);
            plot(ax1, k2_range, E_list);
            title(ax1, sprintf('k1 = %d, B = %d', k1, obj.B));
        end
        
        function plot_butterfly(obj)
            
        end
        % general test of all functions
        function test(obj)
            k1 = 0;
            k2 = 0;
            % test if H is correctly set up
            obj.precompute_init(10, 1);
            obj.construct_H(k1, k2, true);
            % test ll fan plot
            obj.plot_band_fix_k1k2(0, 0);
            % test E versus k2 plot
            obj.plot_band_fix_pq(10, 1);
        end
    
        % plot T coeffs
        function tilde_T = t_coeff(obj, n, np)
            % np > n 
            % use cached function to store results of (n, np)
            % need to test if any efficiency is achieved
            tilde_T = zeros(2,2,3);
            if np >= n
                common_scale = sqrt(2^np * factorial(n) / (2^n * factorial(np)))*...
                    obj.w;
                tilde_T(:,:,1) = common_scale * obj.Ts(:,:,1) *...
                    ((sqrt(3)/4 - 1i/4)*obj.lB*obj.KD)^(np -n) *...
                    exp(-1/4*(obj.lB*obj.KD)^2) *...
                    exp(1i*sqrt(3)/2*obj.lB^2*obj.KD*obj.Kx) *...
                    laguerreL(n, np - n, 1/2*obj.lB^2*obj.KD^2);
                tilde_T(:,:,2) = common_scale * obj.Ts(:,:,2) *...
                    ((1i/2)*obj.lB*obj.KD)^(np -n) *...
                    exp(-1/4*(obj.lB*obj.KD)^2) *...
                    laguerreL(n, np - n, 1/2*obj.lB^2*obj.KD^2);
                tilde_T(:,:,3) = common_scale * obj.Ts(:,:,3) *...
                    ((-sqrt(3)/4 - 1i/4)*obj.lB*obj.KD)^(np -n) *...
                    exp(-1/4*(obj.lB*obj.KD)^2) *...
                    exp(-1i*sqrt(3)/2*obj.lB^2*obj.KD*obj.Kx) *...
                    laguerreL(n, np - n, 1/2*obj.lB^2*obj.KD^2);
            else
                common_scale = sqrt(2^n * factorial(np) / (2^np * factorial(n)))*...
                    obj.w;
                tilde_T(:,:,1) = common_scale * obj.Ts(:,:,1) *...
                    ((-sqrt(3)/4 - 1i/4)*obj.lB*obj.KD)^(n -np) *...
                    exp(-1/4*(obj.lB*obj.KD)^2) *...
                    exp(1i*sqrt(3)/2*obj.lB^2*obj.KD*obj.Kx) *...
                    laguerreL(np, n - np, 1/2*obj.lB^2*obj.KD^2);
                tilde_T(:,:,2) = common_scale * obj.Ts(:,:,2) *...
                    ((1i/2)*obj.lB*obj.KD)^(n -np) *...
                    exp(-1/4*(obj.lB*obj.KD)^2) *...
                    laguerreL(np, n - np, 1/2*obj.lB^2*obj.KD^2);
                tilde_T(:,:,3) = common_scale * obj.Ts(:,:,3) *...
                    ((sqrt(3)/4 - 1i/4)*obj.lB*obj.KD)^(n -np) *...
                    exp(-1/4*(obj.lB*obj.KD)^2) *...
                    exp(-1i*sqrt(3)/2*obj.lB^2*obj.KD*obj.Kx) *...
                    laguerreL(np, n - np, 1/2*obj.lB^2*obj.KD^2);
            end
            
        end
    
        % construct H given k1, k2, and rewrite obj.H
        function construct_H(obj, k1, k2, info_mode)
            arguments
                obj;
                k1;
                k2;
                info_mode = false;
            end
            % basis l = 1, 2, j = 0, n = 0, 1, 2,..., N
            if info_mode
                obj.H_info = string(obj.H);
            end
            % iterate all matrix element (row, col)
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
                    % landau label 0, 1, ..., N - 1 
                    n = floor(mod(mod(row, obj.layer_size), obj.subband_size)/...
                        obj.atom_size);
                    np = floor(mod(mod(col, obj.layer_size), obj.subband_size)/...
                        obj.atom_size);
                    % atom label 0, 1
                    alpha = mod(mod(mod(row, obj.layer_size), obj.subband_size),...
                        obj.atom_size);
                    alphap = mod(mod(mod(col, obj.layer_size), obj.subband_size),...
                        obj.atom_size);
        
                    % TODO
                    % should consider delete the info mode later
                    if info_mode
                        obj.H_info(row_d, col_d) = sprintf("<%d%d%d%d, %d%d%d%d>",...
                            l, j, n, alpha, lp, jp, np, alphap);
                    end
                    
                    % give values
                    % intralayer
                    if l == lp && j == jp
                        if l == 0
                            cur_theta = obj.thetar1;
                        else
                            cur_theta = obj.thetar2;
                        end
                        if n == np + 1 && alpha == 0 && alphap == 1
                            obj.H(row_d, col_d) = sqrt(2)*obj.hbar*obj.vF/obj.lB*...
                                1i*exp(1i*cur_theta)*...
                                sqrt(n); 
                            if info_mode
                                obj.H_info(row_d, col_d) = obj.H_info(row_d, col_d)+...
                                    sprintf('sq(%d)', n);
                            end
                        end
                        if n == np - 1 && alpha == 1 && alphap == 0
                            obj.H(row_d, col_d) = sqrt(2)*obj.hbar*obj.vF/obj.lB*...
                                (-1i)*exp(-1i*cur_theta)*...
                                sqrt(n+1); 
                            if info_mode
                                obj.H_info(row_d, col_d) = obj.H_info(row_d, col_d)+...
                                    sprintf('sq(%d)', n+1);
                            end
                        end
                    end
                    % interlayer
                    if l ~= lp
                        cur_ele = 0;
                        cur_tilde_T = obj.tilde_T(n, np);
                        cur_info = '';
                        if j == mod(jp + 1, obj.q)
                            cur_ele = cur_ele + cur_tilde_T(alpha+1, alpha+1, 1)*...
                            exp(-1i*k2*obj.Delta);
                            if info_mode
                                cur_info = append(cur_info, 't1');
                            end
                        end
                        if j == jp
                            cur_ele = cur_ele + cur_tilde_T(alpha+1, alphap+1, 2)*...
                            exp(1i*3/2*obj.KD*(k1*obj.lB^2+j*obj.Delta));
                            if info_mode
                                cur_info = append(cur_info, 't2');
                            end
                        end
                        if j == mod(jp - 1, obj.q)
                            cur_ele = cur_ele + cur_tilde_T(alpha+1, alphap+1, 3)*...
                            exp(1i*k2*obj.Delta);
                            if info_mode
                                cur_info = append(cur_info, 't3');
                            end
                        end
                        if l == 0
                            obj.H(row_d, col_d) = cur_ele;
                            if info_mode
                                obj.H_info(row_d, col_d) = append(obj.H_info(row_d, col_d), cur_info);
                            end
                        end
                    end
                end
            end
            % for the bottom left block of H, use T' directly
            obj.H((obj.layer_size+1):obj.H_size, 1:obj.layer_size) =...
                obj.H(1:obj.layer_size, (obj.layer_size+1):obj.H_size)';
            if info_mode
                obj.H_info((obj.layer_size+1):obj.H_size, 1:obj.layer_size) =...
                    obj.H_info(1:obj.layer_size, (obj.layer_size+1):obj.H_size)';
            end
        end
    end
end