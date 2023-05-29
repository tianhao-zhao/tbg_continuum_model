classdef TBG < handle
    %TBG Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % computational arguments
        % g network cutoff multiplier number
        network_cutoff_number = 8;
        % double comparing tolerance
        tol = 1e-2;

        % physics constants
        % distance between nearest A and B atoms
        a = 1;
        % single layer graphene vf
        vf = 1;
        hbar = 1;
        
        % first order fourier component 
        w = 0.05;
        % hopping matrix
        w_phase = exp(1i*2*pi/3);
        T1;
        T2;
        T3;
        Ts;

        % real space lattice parameters
        % lattice unit vectors
        a1;
        a2;
        % displacement of nearest neighbours
        a_delta1;
        a_delta2;
        a_delta3;

        % reciprocal space unit vectors
        b1;
        b2;
        % K and K' (K prime) vectors
        K;
        Kp;

        % commensurate index and angle
        cm;
        cn;
        ctheta;
        % rotation matrix
        % note that layer2 is rotated counterclockwise ctheta
        rm1;
        rm2;
        
        % rotated unit vectors
        % a1l1: a1 vector in layer1
        a1l1;
        a2l1;
        a1l2;
        a2l2;
        % a_delta1l1: displacement delta1 in layer1
        a_delta1l1;
        a_delta2l1;
        a_delta3l1;
        
        a_delta1l2;
        a_delta2l2;
        a_delta3l2;

        % super lattice unit vectors
        A1;
        A2;
        % super lattice reciprocal unit vectors
        B1;
        B2;

        % K and Kp in layer1 and 2
        K1;
        K1p;
        K2;
        K2p;

        % unit vectors
        KD;
        % Moire space unit vectors, they are just B1 and B2
        g1;
        g2;
        % Moire space Ks and Ks' (Ksp)
        Ks;
        Ksp;
        
        % displacement of q(k) to its three nearest neighbours
        g_deltas;
        g_network;

        % number of all g's in cutoff range
        M;

        % hamiltonian
        H;
        H_inter;
        H_intra1;
        H_intra2;

        % K as origin, to plot K-gamma-M-K
        K_origin;
        M_point;
        gamma;
        k_line;
        k_line_length;
        % how granular you want to plot
        % TODO add this to input params
        gran = (0:100) ./ 100;

        % energy along K-gamma-M-K line
        all_E_line;


    end
    
    methods
        function obj = TBG(nvargs)
            arguments
                nvargs.cm;
                nvargs.cn;
                nvargs.w (1,1) double = 0.0;
                nvargs.cutoff (1,1) int8 = 8;
                nvargs.display (1,1) logical = 0;
            end
            tic
            % handle inputs
            if isfield(nvargs, 'cm') && isfield(nvargs, 'cn')
                obj.cm = nvargs.cm;
                obj.cn = nvargs.cn;
            elseif ~isfield(nvargs, 'cm') && ~isfield(nvargs, 'cn')
                fprintf('Caution, using default (m, n) = (%d, %d)\n',...
                    obj.cm, obj.cn)
            else
                ME = MException('init:missingInput',...
                    'One of the commensurate index (cm or cn) is missing.');
                throw(ME);
            end
            % init
            init(obj, cm = nvargs.cm, cn = nvargs.cn, w = nvargs.w,...
                cutoff = nvargs.cutoff, display = nvargs.display);
            display_info(obj);  
            toc
        end
        
        function init(obj, nvargs)
            arguments
                obj;
                nvargs.cm;
                nvargs.cn;
                nvargs.w;
                nvargs.cutoff (1,1) double = 8;
                nvargs.display (1,1) logical = 0;
            end
            % handle inputs
            if isfield(nvargs, 'cm') && isfield(nvargs, 'cn')
                obj.cm = nvargs.cm;
                obj.cn = nvargs.cn;
            elseif ~isfield(nvargs, 'cm') && ~isfield(nvargs, 'cn')
                fprintf('Caution, using default or last set (m, n) = (%d, %d)\n',...
                    obj.cm, obj.cn)
            else
                ME = MException('init:missingInput',...
                    'One of the commensurate index (cm or cn) is missing.');
                throw(ME);
            end
            if isfield(nvargs, 'w')
                obj.w = nvargs.w;
            end
            if isfield(nvargs, 'cutoff')
                obj.network_cutoff_number = nvargs.cutoff;
            end
            
            % hopping matrix
            obj.T1 = obj.w * [[1, 1]; [1, 1]];
            obj.T2 = obj.w * [[1, obj.w_phase']; [obj.w_phase, 1]];
            obj.T3 = obj.w * [[1, obj.w_phase]; [obj.w_phase', 1]];
            obj.Ts = cat(3, obj.T1, obj.T2, obj.T3);
            
            % lattice unit vectors
            obj.a1 = obj.a * [3/2; -sqrt(3)/2];
            obj.a2 = obj.a * [3/2; sqrt(3)/2];
            % displacement of nearest neighbours
            obj.a_delta1 = 1/3*(obj.a1 + obj.a2);
            obj.a_delta2 = obj.a_delta1 - obj.a1;
            obj.a_delta3 = obj.a_delta1 - obj.a2;
            % reciprocal space unit vectors
            obj.b1 = 2 * pi / obj.a * [1/3; -sqrt(3)/3];
            obj.b2 = 2 * pi / obj.a * [1/3; sqrt(3)/3];
            % K and K' (K prime) vectors
            obj.K = 2 * pi / obj.a * [1/3; sqrt(3)/9];
            obj.Kp = 2 * pi / obj.a * [1/3; -sqrt(3)/9];

            % commensurate index and angle
            obj.ctheta = acos((obj.cn^2 + 4*obj.cn*obj.cm + obj.cm^2)/...
                (2 * (obj.cn^2 + obj.cn*obj.cm + obj.cm^2)));
            % rotation matrix
            % note that layer2 is rotated counterclockwise ctheta
            obj.rm1 = [[cos(-obj.ctheta/2), -sin(-obj.ctheta/2)];...
                [sin(-obj.ctheta/2), cos(-obj.ctheta/2)]];
            obj.rm2 = [[cos(obj.ctheta/2), -sin(obj.ctheta/2)];...
                [sin(obj.ctheta/2), cos(obj.ctheta/2)]];

            % rotated unit vectors
            % a1l1: a1 vector in layer1
            obj.a1l1 = obj.rm1 * obj.a1;
            obj.a2l1 = obj.rm1 * obj.a2;
            obj.a1l2 = obj.rm2 * obj.a1;
            obj.a2l2 = obj.rm2 * obj.a2;
            % a_delta1l1: displacement delta1 in layer1
            obj.a_delta1l1 = obj.rm1 * obj.a_delta1;
            obj.a_delta2l1 = obj.rm1 * obj.a_delta2;
            obj.a_delta3l1 = obj.rm1 * obj.a_delta3;
            
            obj.a_delta1l2 = obj.rm2 * obj.a_delta1;
            obj.a_delta2l2 = obj.rm2 * obj.a_delta2;
            obj.a_delta3l2 = obj.rm2 * obj.a_delta3;

            % super lattice unit vectors
            obj.A1 = obj.cn*obj.a1l1 + obj.cm*obj.a2l1;
            obj.A2 = -obj.cm*obj.a1l1 + (obj.cm+obj.cn)*obj.a2l1;
            % super lattice reciprocal unit vectors
            obj.B1 = (obj.cm+obj.cn)/(obj.cm^2+obj.cn^2+obj.cm*obj.cn) *...
                obj.rm1*obj.b1 + obj.cm/(obj.cm^2+obj.cn^2+obj.cm*obj.cn) *...
                obj.rm1*obj.b2;
            obj.B2 = -obj.cm/(obj.cm^2+obj.cn^2+obj.cm*obj.cn) *...
                obj.rm1*obj.b1 + obj.cn/(obj.cm^2+obj.cn^2+obj.cm*obj.cn) *...
                obj.rm1*obj.b2;

            % K and Kp in layer1 and 2
            obj.K1 = obj.rm1 * obj.K;
            obj.K1p = obj.rm1 * obj.Kp;
            obj.K2 = obj.rm2 * obj.K;
            obj.K2p = obj.rm2 * obj.Kp;
            
            % construct g vector network space
            % unit vectors
            obj.KD = obj.K1 - obj.K2;
            % Moire space unit vectors, they are just B1 and B2
            obj.g1 = [0; 2 * obj.KD(2)];
            obj.g2 = [3 * obj.KD(1); obj.KD(2)];
            % Moire space Ks and Ks' (Ksp)
            obj.Ks = [obj.KD(1); -obj.KD(2)];
            obj.Ksp = obj.KD;

            % displacement of q(k) to its three nearest neighbours
            obj.g_deltas = [obj.KD, obj.KD - obj.g1, obj.KD - obj.g2];
            % find all g vectors in cutoff range
            search_limit = ceil(3 * obj.network_cutoff_number);
            [search_range_x, search_range_y] = meshgrid(-search_limit:search_limit);
            obj.g_network = [obj.g1, obj.g2] *...
                [reshape(search_range_x, 1, []); reshape(search_range_y, 1, [])];
            obj.g_network = obj.g_network(:, vecnorm(obj.g_network + obj.KD, 2, 1)<=...
                obj.network_cutoff_number * norm(obj.KD));
            
            obj.M = size(obj.g_network, 2);
            
            % init H
            obj.H_inter = zeros(2*obj.M);
            
            % fix H inter
            calculate_hinter(obj, nvargs.display);
            % fix H
            obj.H = [zeros(2*obj.M), obj.H_inter;...
                obj.H_inter', zeros(2*obj.M)];

            % prepare K gamma M K line
            obj.K_origin = [0; 0];
            obj.M_point = (obj.Ks + obj.Ksp) / 2;
            % when origin is K point, gamma is Ksp
            obj.gamma = obj.Ksp;
            obj.k_line = [(obj.gamma - obj.K_origin) * obj.gran + obj.K_origin, ...
                (obj.M_point - obj.gamma) * obj.gran + obj.gamma, ...
                (obj.K_origin - obj.M_point) * obj.gran + obj.M_point];
            obj.k_line_length = [vecnorm((obj.gamma - obj.K_origin) * obj.gran, 2), ...
                vecnorm((obj.M_point - obj.gamma) * obj.gran, 2) + norm(obj.gamma - obj.K_origin), ...
                vecnorm((obj.K_origin - obj.M_point) * obj.gran, 2) + norm(obj.gamma - obj.K_origin) +...
                norm(obj.M_point - obj.gamma)];

            % prepare 3d mesh
            % TODO improve to hexagon shape

            obj.display_info();

        end

        % display parameters used to init the tbg instance
        function display_info(obj)
            fprintf('cm = %d, cn = %d\n', obj.cm, obj.cn)
            fprintf('cutoff number = %d\n', obj.network_cutoff_number)
            fprintf('w = %d\n', obj.w)
        end

        function display_all(obj)
        end

        function calculate_hinter(obj, display)
            arguments
                obj
                display (1,1) logical = 0;
            end
            if display
                fig_network = figure();
                ax_network = axes(fig_network);
                axis(ax_network, 'equal');
            end
            % TODO improve searching algo
            % now it is O(M^2) and with for loop
            for ii = 1:obj.M
                for jj = 1:obj.M
                    % ii in layer1, jj in layer2
                    dist = (obj.g_network(:, jj) + obj.KD) - obj.g_network(:, ii);
                    for dd = 1:3
                        if norm(obj.g_deltas(:, dd) - dist) <= obj.tol
                            % not sure if conjugate should be reversed 
                            obj.H_inter(ii*2-1:ii*2, jj*2-1:jj*2) = obj.Ts(:, :, dd)';
                            % for test if all connections are found
                            if display   
                                cur_line = [obj.g_network(:, ii), ...
                                    obj.g_network(:, jj) + obj.KD];
                                line(ax_network, cur_line(1, :), cur_line(2, :));
                            end
                        end
                    end
                end
            end
            
        end

        function calculate_h(obj, k)
            % TODO change the loop to vectorization
            
            k = k(:);
            for ii = 1:obj.M
                cur_k_l1 = k + obj.g_network(:, ii);
                cur_k_l2 = cur_k_l1 + obj.KD;
                obj.H(2*ii - 1, 2*ii) = obj.hbar * obj.vf * ...
                    (cur_k_l1(1) - 1i * cur_k_l1(2)) * exp(-1i * obj.ctheta/2);
                obj.H(2*ii, 2*ii - 1) = obj.hbar * obj.vf * ...
                    (cur_k_l1(1) + 1i * cur_k_l1(2)) * exp(1i * obj.ctheta/2);
                obj.H(2*ii - 1 + 2*obj.M, 2*ii + 2*obj.M) = obj.hbar * obj.vf * ...
                    (cur_k_l2(1) - 1i * cur_k_l2(2)) * exp(1i * obj.ctheta/2);
                obj.H(2*ii + 2*obj.M, 2*ii - 1 + 2*obj.M) = obj.hbar * obj.vf * ...
                    (cur_k_l2(1) + 1i * cur_k_l2(2)) * exp(-1i * obj.ctheta/2);
            end
        end

        function plot_line(obj)
            tic
            obj.all_E_line = zeros(4*obj.M, size(obj.k_line, 2));
            for ind = 1:size(obj.k_line, 2)
                k = obj.k_line(:, ind);
                obj.calculate_h(k);
                obj.all_E_line(:, ind) = eig(obj.H);
            end
            fig_line = figure();
            ax_line = axes(fig_line);
            for ii = 1:4*obj.M
                plot(ax_line, obj.k_line_length, obj.all_E_line(ii, :));
                hold on;
            end 
            hold off;
            fig_line_zoomin = figure();
            ax_line_zoomin = axes(fig_line_zoomin);
            for ii = 2*obj.M-5:2*obj.M+6
                plot(ax_line_zoomin, obj.k_line_length, obj.all_E_line(ii, :));
                hold on;
            end
            hold off;
            gran_size = size(obj.gran, 2);
            line_ticks = obj.k_line_length([1, gran_size * (1:3)]);
            for ax = [ax_line, ax_line_zoomin]
                xticks(ax, line_ticks);
                xticklabels(ax, {'K', '\Gamma', 'M', 'K'});
            end
            vf_e = (obj.all_E_line(2*obj.M+1, 5) - obj.all_E_line(2*obj.M+1, 1)) / ...
                (obj.k_line_length(5) - obj.k_line_length(1));
            title(ax_line_zoomin, sprintf('effective vf / vf_0 %d', vf_e));
            toc
        end

        function plot_3d(obj)
        end

    end
end

