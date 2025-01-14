classdef CM < handle
    % a continuum model handle stores all vector values 
    % provides h_inter calculation function
    
    properties
        % almost fixed constants
        tol = 1e-8;
        electron = 1.6e-19;

        % parameters
        theta;
        w;
        % corrugation and strain correction
        eta;
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
        m;
        H;

        % magnetic field 
        % cutoff number of hermite states |n>
        N_B;
        % a and a^dagger operators
        a_B;
        ad_B;
        H_B;
    end
    
    methods
        function obj = CM(nvargs)
            arguments
                nvargs.theta = 1.1;
                nvargs.w = 1.6e-19*110e-3;
                nvargs.eta = 1;
                nvargs.cutoff = 8;
                nvargs.a = 1.4e-10;
                nvargs.vF = 0.866e6;
                nvargs.hbar = 1.05e-34;
                nvargs.N_B = 40;
            end
            obj.init(theta=nvargs.theta, w=nvargs.w, eta = nvargs.eta, cutoff=nvargs.cutoff,...
                a=nvargs.a, vF=nvargs.vF, hbar=nvargs.hbar, N_B=nvargs.N_B);
        end
        
        function init(obj, nvargs)
            arguments
                obj;
                nvargs.theta;
                nvargs.w;
                nvargs.eta;
                nvargs.cutoff;
                nvargs.a;
                nvargs.vF;
                nvargs.hbar;
                nvargs.N_B;
            end
            if isfield(nvargs, 'theta')
                obj.theta = nvargs.theta;
            end
            if isfield(nvargs, 'w')
                obj.w = nvargs.w;
            end
            if isfield(nvargs, 'eta')
                obj.eta = nvargs.eta;
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
            if isfield(nvargs, 'N_B')
                obj.N_B = nvargs.N_B;
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
            T1 = [[obj.eta*1, 1]; [1, obj.eta*1]];
            T2 = [[obj.eta*1, obj.w_phase]; [obj.w_phase', obj.eta*1]];
            T3 = [[obj.eta*1, obj.w_phase']; [obj.w_phase, obj.eta*1]];
            obj.Tjs(:, :, 1) = T1;
            obj.Tjs(:, :, 2) = T2;
            obj.Tjs(:, :, 3) = T3;
            
            obj.build_g_network;
            obj.H = zeros(2*(obj.m1 + obj.m2));
            obj.build_h_inter;

            % % prepare for magnetic calculation
            % obj.a_B = diag(sqrt(1:obj.N_B-1), -1);
            % obj.ad_B = diag(sqrt(1:obj.N_B-1), 1);
            % obj.H_B = zeros(2*obj.m*obj.N_B);
            % obj.build_h_b_inter();
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
            obj.m = obj.m1 + obj.m2;
        end

        function build_h_inter(obj)
            arguments
                obj;
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
            obj.H(1:(2*obj.m1), (1+2*obj.m1):(2*(obj.m))) = H_inter;
            obj.H((1+2*obj.m1):(2*(obj.m)), 1:(2*obj.m1)) = H_inter';
        end
        
        function build_h_intra1(obj, q)
            arguments
                obj;
                q;
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

        function build_h_intra2(obj, q)
            arguments
                obj;
                q;
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
            E = zeros(2*(obj.m), size(q, 2));
            for i = 1:size(q, 2)
                obj.build_h_intra1(q(:, i));
                obj.build_h_intra2(q(:, i));
                e = eig(obj.H);
                E(:, i) = e;
            end
        end
        
        % TODO
        % find a way to utilize conjugate and quickly build matrix

        % see LL.m for nonzero field calculation
        % % functions for calculating landau levels under magnetic field
        % % B unit: T
        % % for layer 1, it is a and ad
        % % for layer 2, it is a - kd and ad - kd
        % % for now, rotation is not included e^(itheta_K)
        % function build_h_b_intra1(obj, B)
        %     lc = sqrt(obj.hbar/(obj.electron*B));
        %     for i = 1:obj.m1
        %         obj.H_B((1:obj.N_B)+2*(i-1)*obj.N_B, (obj.N_B+1:2*obj.N_B)+2*(i-1)*obj.N_B) =...
        %             obj.hbar*obj.vF*(sqrt(2)/lc * obj.a_B);
        %         obj.H_B((obj.N_B+1:2*obj.N_B)+2*(i-1)*obj.N_B, (1:obj.N_B)+2*(i-1)*obj.N_B) =...
        %             obj.hbar*obj.vF*(sqrt(2)/lc * obj.ad_B);
        %     end
        % end
        % 
        % function build_h_b_intra2(obj, B)
        %     lc = sqrt(obj.hbar/(obj.electron*B));
        %     for i = 1:obj.m2
        %         obj.H_B((1:obj.N_B)+2*(i-1+obj.m1)*obj.N_B, (obj.N_B+1:2*obj.N_B)+2*(i-1+obj.m1)*obj.N_B) =...
        %             obj.hbar*obj.vF*(sqrt(2)/lc * obj.a_B -...
        %             (obj.KD(1) - 1i*obj.KD(2))*eye(obj.N_B));
        %         obj.H_B((obj.N_B+1:2*obj.N_B)+2*(i-1+obj.m1)*obj.N_B, (1:obj.N_B)+2*(i-1+obj.m1)*obj.N_B) =...
        %             obj.hbar*obj.vF*(sqrt(2)/lc * obj.ad_B -...
        %             (obj.KD(1) + 1i*obj.KD(2))*eye(obj.N_B));
        %     end
        % end
        % 
        % % this is not dependent on B
        % % should only be called once
        % % same code with build_h_inter
        % function build_h_b_inter(obj)
        %     T_dim = size(obj.Tjs, 1);
        %     H_inter_B = zeros(obj.m1*T_dim*obj.N_B, obj.m2*T_dim*obj.N_B);
        %     for i = 1:obj.m1
        %         for j = 1:obj.m2
        %             for t_ind = 1:size(obj.tjs, 2)
        %                 if obj.compare_vectors(obj.g_network1(:,i) - obj.g_network2(:, j),...
        %                         obj.tjs(:, t_ind), obj.tol)
        %                     H_inter_B((1+(i-1)*T_dim*obj.N_B):i*T_dim*obj.N_B,...
        %                     (1+(j-1)*T_dim*obj.N_B):j*T_dim*obj.N_B) =...
        %                     obj.w * [obj.Tjs(1, 1, t_ind) * eye(obj.N_B),...
        %                     obj.Tjs(1, 2, t_ind) * eye(obj.N_B);...
        %                     obj.Tjs(2, 1, t_ind) * eye(obj.N_B),...
        %                     obj.Tjs(2, 2, t_ind) * eye(obj.N_B)];
        %                 end
        %             end
        %         end
        %     end
        %     obj.H_B(1:(2*obj.m1*obj.N_B), (1+2*obj.m1*obj.N_B):(2*obj.m*obj.N_B)) = H_inter_B;
        %     obj.H_B((1+2*obj.m1*obj.N_B):(2*obj.m*obj.N_B), 1:(2*obj.m1*obj.N_B)) = H_inter_B';
        % end
        % 
        % function E = calculate_e_b(obj, Bs)
        %     Bs = Bs(:);
        %     E = zeros(size(obj.H_B, 1), numel(Bs));
        %     for i = 1:numel(Bs)
        %         obj.build_h_b_intra1(Bs(i));
        %         obj.build_h_b_intra2(Bs(i));
        %         E(:, i) = eig(obj.H_B);
        %     end
        % end

        % functions for plotting
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
            num_of_bands_shown = floor(min(5, E_len/2));
            for i = (E_len/2-num_of_bands_shown+1):(E_len/2+num_of_bands_shown)
                plot(ax, q_len, E(i, :)/obj.electron, Color='k');
                hold(ax, 'on');
            end
            line_ticks = q_len([1, q_num_of_steps*[1:3]]);
            xticks(ax, line_ticks);
            xticklabels(ax, {'K_M', '\Gamma', 'M_M', 'K''_M'});
            xlim([0, q_len(length(q_len))]);
            xline(line_ticks(2:3), Color=[0.5, 0.5, 0.5]);
            plot(ax, q_len(1:q_num_of_steps), q_len(1:q_num_of_steps)*...
                obj.hbar*obj.vF/obj.electron, '--', Color='r');
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
            num_of_bands_shown = floor(min(2, E_len/2));
            for i = (E_len/2-num_of_bands_shown + 1):(E_len/2+num_of_bands_shown)
                surf(ax, reshape(xx, sx), reshape(yy, sx), reshape(E(i, :), sx));
                hold(ax, 'on');
            end
        end
    
        function info(obj)
            disp('\theta =');
            disp(obj.theta);
            disp('coupling w in eV =');
            disp(obj.w / (1.6e-16) * 1000);
            disp('vF in m/s =');
            disp(obj.vF)
            disp('corrugation correction ratio =')
            disp(obj.eta);
            % disp(obj.N_B);
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
            obj.plot_3d;
            obj.plot_kgammamk;

            % tic;
            % Bs = 1:20;
            % E = obj.calculate_e_b(Bs);
            % fig_b = figure();
            % ax_b = axes(fig_b);
            % for i = 1:size(E, 1)
            %     plot(ax_b, Bs, E(i, :)/obj.electron);
            %     hold(ax_b, 'on');
            % end
            % toc;
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

