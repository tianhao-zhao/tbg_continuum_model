classdef CM < handle
    % a continuum model handle stores all vector values 
    % provides h_inter calculation function
    
    properties
        % almost fixed constants
        tol = 1e-8;

        % parameters
        theta;
        w;
        cutoff;
        a;
        vF;
        hbar;

        % calculated values
        rtheta;
        rm1;
        rm2;
        K;
        Kp;
        b1;
        b2;
        K1;
        K2;
        theta_K1;
        theta_K2;
        % G_j, j=1,2,3
        Gjs;
        KD;
        % g_j, j=1,2,3
        gjs;
        tjs;
        w_phase;
        Tjs;
        g_network1;
        g_network2;
        m1;
        m2;
        H;
    end
    
    methods
        function obj = CM(nvargs)
            arguments
                nvargs.theta = 1.1;
                nvargs.w = 1.6e-19*110e-3;
                nvargs.cutoff = 8;
                nvargs.a = 1.4e-10;
                nvargs.vF = 0.866e6;
                nvargs.hbar = 1.05e-34;
            end
            obj.init(theta=nvargs.theta, w=nvargs.w, cutoff=nvargs.cutoff,...
                a=nvargs.a, vF=nvargs.vF, hbar=nvargs.hbar);
        end
        
        function init(obj, nvargs)
            arguments
                obj;
                nvargs.theta;
                nvargs.w;
                nvargs.cutoff;
                nvargs.a;
                nvargs.vF;
                nvargs.hbar;
            end
            if isfield(nvargs, 'theta')
                obj.theta = nvargs.theta;
            end
            if isfield(nvargs, 'w')
                obj.w = nvargs.w;
            end
            if isfield(nvargs, 'cutoff')
                obj.cutoff = nvargs.cutoff;
            end
            if isfield(nvargs, 'a')
                obj.a = nvargs.a;
            end
            if isfield(nvargs, 'vF')
                obj.vF = nvargs.vF;
            end
            if isfield(nvargs, 'hbar')
                obj.hbar = nvargs.hbar;
            end
            % assume all property are defined
            % use property to calculate everything ready
            obj.rtheta = obj.theta / 180 * pi;
            obj.rm1 = obj.rotation(-obj.rtheta/2);
            obj.rm2 = obj.rotation(obj.rtheta/2);
            obj.K = 2*pi/obj.a * [1/3; sqrt(3)/9];
            obj.Kp = 2*pi/obj.a * [1/3; -sqrt(3)/9];
            obj.b1 = 2*pi/obj.a * [1/3; -sqrt(3)/3];
            obj.b2 = 2*pi/obj.a * [1/3; sqrt(3)/3];
            obj.K1 = obj.rm1 * obj.K;
            obj.K2 = obj.rm2 * obj.K;
            obj.theta_K1 = atan(obj.K1(2)/obj.K1(1));
            obj.theta_K2 = atan(obj.K2(2)/obj.K2(1));
            G1 = [0; 0];
            G2 = obj.b2;
            G3 = obj.b1 + obj.b2;
            obj.Gjs = [G1, G2, G3];
            obj.KD = obj.K2 - obj.K1;
            g1 = (obj.rm1 - obj.rm2) * G1;
            g2 = (obj.rm1 - obj.rm2) * G2;
            g3 = (obj.rm1 - obj.rm2) * G3;
            obj.gjs = [g1, g2, g3];
            obj.tjs = obj.KD + obj.gjs;
            obj.w_phase = exp(1i*2*pi/3);
            T1 = [[1, 1]; [1, 1]];
            T2 = [[1, obj.w_phase]; [obj.w_phase', 1]];
            T3 = [[1, obj.w_phase']; [obj.w_phase, 1]];
            obj.Tjs(:, :, 1) = T1;
            obj.Tjs(:, :, 2) = T2;
            obj.Tjs(:, :, 3) = T3;
            
            obj.build_g_network;
            obj.H = zeros(2*(obj.m1 + obj.m2));
            obj.build_h_inter;
        end

        function build_g_network(obj)
            max_range = ceil(obj.cutoff * 2);
            [x_range_m, y_range_m] = meshgrid(-max_range:max_range);
            x_range = x_range_m(:);
            y_range = y_range_m(:);
            g_candidates1 = [obj.gjs(:,2), obj.gjs(:,3)] * [x_range, y_range]';
            g_candidates2 = g_candidates1 - [obj.KD];

            % TODO
            % only consider distance now
            % add option to make all vector 6 fold symmetrical
            obj.g_network1 = g_candidates1(:, vecnorm(g_candidates1, 2, 1) <= ...
                obj.cutoff * norm(obj.KD));
            obj.g_network2 = g_candidates2(:, vecnorm(g_candidates2, 2, 1) <= ...
                obj.cutoff * norm(obj.KD));
            obj.m1 = size(obj.g_network1, 2);
            obj.m2 = size(obj.g_network2, 2);
        end

        function build_h_inter(obj, dimension)
            % TODO
            % change code to be compatible with dimen > 1
            arguments
                obj;
                dimension = 1;
            end
            T_dim = size(obj.Tjs, 1);
            H_inter = zeros(obj.m1*T_dim, obj.m2*T_dim);
            for i = 1:obj.m1
                for j = 1:obj.m2
                    for t_ind = 1:size(obj.tjs, 2)
                        if obj.compare_vectors(obj.g_network1(:,i) - obj.g_network2(:, j),...
                                obj.tjs(:, t_ind), obj.tol)
                            H_inter((1+(i-1)*T_dim):i*T_dim,...
                            (1+(j-1)*T_dim):j*T_dim) = obj.w * obj.Tjs(:, :, t_ind);
                        end
                    end
                end
            end
            obj.H(1:(2*obj.m1), (1+2*obj.m1):(2*(obj.m1+obj.m2))) = H_inter;
            obj.H((1+2*obj.m1):(2*(obj.m1+obj.m2)), 1:(2*obj.m1)) = H_inter';
        end
        
        function build_h_intra1(obj, q, dimension)
            arguments
                obj;
                q;
                dimension = 1;
            end
            q = q(:);
            for i = 1:obj.m1
                cur_q = q + obj.g_network1(:, i);
                obj.H(2*i-1, 2*i) = obj.hbar * obj.vF *...
                    (cur_q(1) - 1i * cur_q(2)) * exp(1i * obj.theta_K1);
                obj.H(2*i, 2*i-1) = obj.hbar * obj.vF *...
                    (cur_q(1) + 1i * cur_q(2)) * exp(-1i * obj.theta_K1);
            end
        end

        function build_h_intra2(obj, q, dimension)
            arguments
                obj;
                q;
                dimension = 1;
            end
            q = q(:);
            for i = 1:obj.m2
                cur_q = q + obj.g_network2(:, i);
                obj.H(2*i-1+2*obj.m1, 2*i+2*obj.m1) = obj.hbar * obj.vF *...
                    (cur_q(1) - 1i * cur_q(2)) * exp(1i * obj.theta_K2);
                obj.H(2*i+2*obj.m1, 2*i-1+2*obj.m1) = obj.hbar * obj.vF *...
                    (cur_q(1) + 1i * cur_q(2)) * exp(-1i * obj.theta_K2);
            end
        end

        function E = calculate_e(obj, q)
            % rotate q, won't rotate if size(q) is 2x2
            if size(q, 1) > 2 
                q = q';
            elseif size(q, 1) < 2
                q = q';
            end
            E = zeros(2*(obj.m1+obj.m2), size(q, 2));
            for i = 1:size(q, 2)
                obj.build_h_intra1(q(:, i));
                obj.build_h_intra2(q(:, i));
                e = eig(obj.H);
                E(:, i) = e;
            end
        end

        function [fig, ax] = plot_kgammamk(obj)
            fig = figure();
            ax = axes(fig);
            q_k = [0; 0];
            q_gamma = norm(obj.KD)*[1/2; sqrt(3)/2];
            q_m = -1/6 * (obj.gjs(:, 2) + obj.gjs(:, 3));
            q_kp = -1/3 * (obj.gjs(:, 2) + obj.gjs(:,3));
            q_num_of_steps = 100;
            q_steps = 1/q_num_of_steps:1/q_num_of_steps:1;
            q = [q_k, q_k+(q_gamma-q_k)*q_steps, q_gamma+(q_m-q_gamma)*q_steps,...
                q_m+(q_kp-q_m)*q_steps];
            q_len = [0, norm(q_gamma-q_k)*q_steps, norm(q_gamma-q_k)+norm(q_m-q_gamma)*q_steps,...
                norm(q_gamma-q_k)+norm(q_m-q_gamma)+norm(q_kp-q_m)*q_steps];
            E = obj.calculate_e(q);
            E_len = size(E, 1);
            for i = (E_len/2-4):(E_len/2+5)
                plot(ax, q_len, E(i, :));
                hold(ax, 'on');
            end
        end

        function plot_3d(obj)
            [xx, yy] = meshgrid(norm(obj.KD)*0.1*(-15:15));
            sx = size(xx);
            xx = xx(:);
            yy = yy(:);
            q = [xx, yy];
            E = obj.calculate_e(q);
            fig = figure();
            ax = axes(fig);
            E_len = size(E, 1);
            for i = (E_len/2-1):(E_len/2+2)
                surf(ax, reshape(xx, sx), reshape(yy, sx), reshape(E(i, :), sx));
                hold(ax, 'on');
            end
        end
    
        function info(obj)
            disp(obj.theta);
            disp(obj.w / (1.6e-16) * 1000);
        end

        function run_test(obj)
            obj.build_g_network;
            fig1 = figure();
            ax1 = axes(fig1);
            scatter(ax1, obj.g_network1(1, :), obj.g_network1(2, :));
            hold (ax1, 'on');
            scatter(ax1, obj.g_network2(1, :), obj.g_network2(2, :));
            axis (ax1, 'equal');
            % disp(obj.H)
            %obj.plot_3d;
            obj.plot_kgammamk;
        end
    end
    methods(Static)
        function rm = rotation(theta)
            rm = [[cos(theta), -sin(theta)]; [sin(theta), cos(theta)]];
        end

        function res = compare_vectors(v1, v2, tol)
            if ~all(size(v1) == size(v2))
                res = false;
                return;
            else
                dv = v1 - v2;
                if sumsqr(dv) > tol
                    res = false;
                    return;
                end
            end
            res = true;
        end
    end
end

