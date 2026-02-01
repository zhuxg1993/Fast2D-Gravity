function [Fxz,Gxz] = cal_Fxz_mothed_4(mod, gobs, nodex, nodez)
% Initialize the output matrix
Fxz = zeros(size(gobs, 1), size(mod, 1) * size(mod, 2));

% Physical constants
G = 6.67430e-11;
scale_factor = 1e12;

% Threshold settings
eps_val = 1e-10;         % Singularity check threshold
dz_eps = 1e-10;          % Horizontal edge threshold

% Get grid dimensions
[NZ, NX] = size(mod);
NOTEX = NX + 1;
NOTEZ = NZ + 1;
num_points = size(gobs,1);  % Number of observation points

% Preallocate horizontal and vertical edge storage arrays
H_edge = zeros(NOTEZ, NX);  % Horizontal edge contributions (NOTEZ x NX)
V_edge = zeros(NZ, NOTEX);  % Vertical edge contributions (NZ x NOTEX)

% Preallocate node-based quantities
node_atan = zeros(NOTEZ, NOTEX); % Store atan2 values at each node
node_ln = zeros(NOTEZ, NOTEX);   % Store logarithmic values at each node
dis_x = zeros(NOTEZ, NOTEX);     % Store x-distance at each node
dis_z = zeros(NOTEZ, NOTEX);     % Store z-distance at each node

% ====== Record start time ======
start_time = tic;
% ===============================

for jj = 1:num_points
    
    xx = gobs(jj, 1);
    zz = gobs(jj, 2);
    
    % Precompute atan2 and logarithmic values at each node
    for j = 1:NOTEZ
        for i = 1:NOTEX
            dis_x(j, i) = nodex(j, i) - xx;
            dis_z(j, i) = nodez(j, i) - zz;
            
            % Singularity check
            if abs(dis_x(j, i)) < eps_val && abs(dis_z(j, i)) < eps_val
                node_atan(j, i) = 0.0;
                node_ln(j, i) = 0.0;
                continue;
            end
            
            % Compute atan2 value (note MATLAB atan2 argument order: y, x)
            node_atan(j, i) = atan2(dis_z(j, i), dis_x(j, i));
            
            % Compute logarithmic value
            node_ln(j, i) = log(dis_x(j, i)^2 + dis_z(j, i)^2);
        end
    end
    
    % Compute all top and bottom edge contributions (j: 1..NOTEZ, i: 1..NX)
    for j = 1:NOTEZ
        for i = 1:NX
            x1 = dis_x(j, i);
            z1 = dis_z(j, i);
            x2 = dis_x(j, i+1);
            z2 = dis_z(j, i+1);
            
            % Singularity check - skip this edge if an endpoint is close to the observation point
            if (abs(x1) < eps_val && abs(z1) < eps_val) || ...
                    (abs(x2) < eps_val && abs(z2) < eps_val)
                H_edge(j, i) = 0.0;
                continue;
            end
            
            % Horizontal edge check
            dz_edge = z2 - z1;
            is_level = abs(dz_edge) < dz_eps;
            
            if is_level
                % Horizontal edge contribution
                H_edge(j, i) = 0;
            else
                % General case formula
                c1 = node_atan(j, i);
                c2 = node_atan(j, i+1);
                
                % Adjust the angle to the range 0¨C2¦Ð
                if c1 < 0, c1 = c1 + 2*pi; end
                if c2 < 0, c2 = c2 + 2*pi; end
                
                % Handle angle difference wrapping
                delta_c = c2 - c1;
                if delta_c > pi
                    delta_c = delta_c - 2*pi;
                elseif delta_c < -pi
                    delta_c = delta_c + 2*pi;
                end
                
                % Compute the difference of logarithmic terms
                c3 = node_ln(j, i+1) - node_ln(j, i);
                
                % Compute edge vector
                dx = x2 - x1;
                dz = z2 - z1;
                denominator = dx^2 + dz^2;
                
                % Avoid division by zero
                if denominator < 1e-12
                    H_edge(j, i) = 0.0;
                    continue;
                end
                
                % Compute numerator term
                numerator = dz;
                
                % Compute the contribution of the current edge
                factor = dx * delta_c - 0.5 * dz * c3;
                H_edge(j, i) = 2.0 * G * (numerator / denominator) * factor * scale_factor;
            end
        end
    end
    
    % Compute all vertical edge contributions (j: 1..NZ, i: 1..NOTEX)
    for j = 1:NZ
        for i = 1:NOTEX
            x1 = dis_x(j, i);
            z1 = dis_z(j, i);
            x2 = dis_x(j+1, i);
            z2 = dis_z(j+1, i);
            
            % Singularity check
            if (abs(x1) < eps_val && abs(z1) < eps_val) || ...
                    (abs(x2) < eps_val && abs(z2) < eps_val)
                V_edge(j, i) = 0.0;
                continue;
            end
            
            % Vertical edge formula
            c3 = node_ln(j+1, i) - node_ln(j, i);
            V_edge(j, i) = -G * c3 * scale_factor;
        end
    end
    
    % Traverse all cells to assemble kernel values (column-major order)
    n = 1; % Grid index
    for i = 1:NX        % Column loop (X direction)
        for j = 1:NZ    % Row loop (Z direction)
            % Combine the contributions of the four edges:
            %   Top edge:    H_edge(j, i)
            %   Right edge:  V_edge(j, i+1)
            %   Bottom edge: H_edge(j+1, i) (negative)
            %   Left edge:   V_edge(j, i)   (negative)
            
            f = H_edge(j, i) ...        % Top edge (positive contribution)
                + V_edge(j, i+1) ...    % Right edge (positive contribution)
                - H_edge(j+1, i) ...    % Bottom edge (negative contribution)
                - V_edge(j, i);         % Left edge (negative contribution)
            
            Fxz(jj, n) = f;
            n = n + 1;
        end
    end
end

% F*m = g
Gxz = Fxz * reshape(mod,[],1);

% ====== Output results ======
computation_time = toc(start_time);
fprintf('Computation of gravity kernel matrix Fxz completed, time is: %.6f seconds.\n', computation_time);

end
