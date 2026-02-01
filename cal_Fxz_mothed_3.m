function [Fxz,Gxz] = cal_Fxz_mothed_3(mod, gobs, nodex, nodez)
% Initialize the output matrix
Fxz = zeros(size(gobs, 1), size(mod, 1) * size(mod, 2));

% Physical constants
G = 6.67430e-11;
scale_factor = 1e12;

% Singularity handling thresholds
eps_val = 1e-10;     % Singularity check threshold
dz_eps = 1e-10;      % Horizontal edge threshold

% Get grid dimensions
[grid_rows, grid_cols] = size(mod);
notez = grid_rows + 1;  % Number of horizontal edge rows (NZ+1)
notex = grid_cols + 1;  % Number of vertical edge columns (NX+1)

% Main iteration loop
tic;
for k = 1:size(gobs, 1) % Loop over observation points
    
    x_obs = gobs(k, 1);
    z_obs = gobs(k, 2);
    
    % Preallocate edge contribution matrices
    H_edge = zeros(notez, grid_cols);  % Horizontal edge contributions (NOTEZ ¡Á NX)
    V_edge = zeros(grid_rows, notex);  % Vertical edge contributions (NZ ¡Á NOTEX)
    
    % 1. Compute all horizontal edge contributions (top and bottom edges)
    for j = 1:notez  % j: 1..NZ+1
        for i = 1:grid_cols  % i: 1..NX
            % Get the two endpoints of the current edge
            x1 = nodex(j, i) - x_obs;
            z1 = nodez(j, i) - z_obs;
            x2 = nodex(j, i+1) - x_obs;
            z2 = nodez(j, i+1) - z_obs;
            
            % Singularity check - skip this edge if an endpoint is close to the observation point
            if (abs(x1) < eps_val && abs(z1) < eps_val) || ...
               (abs(x2) < eps_val && abs(z2) < eps_val)
                H_edge(j, i) = 0;
                continue;
            end
            
            % Horizontal edge check (dz ¡Ö 0)
            dz_edge = z2 - z1;
            if abs(dz_edge) < dz_eps
                % Horizontal edge formula
                H_edge(j, i) = 0;
            else
                % General case formula (consistent with the original code)
                c1 = atan2(z1, x1);
                c2 = atan2(z2, x2);
                if c1 < 0, c1 = c1 + 2*pi; end
                if c2 < 0, c2 = c2 + 2*pi; end
                
                % Handle angle difference wrapping
                delta_c = c2 - c1;
                if delta_c > pi
                    delta_c = delta_c - 2*pi;
                elseif delta_c < -pi
                    delta_c = delta_c + 2*pi;
                end
                
                % Compute the logarithmic term
                c3 = log((x2^2 + z2^2) / (x1^2 + z1^2));
                
                % Compute edge vector
                dx = x2 - x1;
                dz = z2 - z1;
                denominator = dx^2 + dz^2;
                
                % Avoid division by zero
                if denominator < 1e-12
                    H_edge(j, i) = 0;
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
    
    % 2. Compute all vertical edge contributions (left and right edges)
    for j = 1:grid_rows  % j: 1..NZ
        for i = 1:notex  % i: 1..NX+1
            % Get the coordinates of the two endpoints of the edge
            x1 = nodex(j, i) - x_obs;
            z1 = nodez(j, i) - z_obs;
            x2 = nodex(j+1, i) - x_obs;
            z2 = nodez(j+1, i) - z_obs;
            
            % Singularity check - skip edges containing the observation point
            if (abs(x1) < eps_val && abs(z1) < eps_val) || ...
               (abs(x2) < eps_val && abs(z2) < eps_val)
                V_edge(j, i) = 0;
                continue;
            end
            
            % Vertical edge formula (using x1 as the reference point)
            c3 = log((x2^2 + z2^2) / (x1^2 + z1^2));
            V_edge(j, i) = -G * c3 * scale_factor;
        end
    end
    
    % 3. Combine contributions of all grid cells (with corrected indexing)
    n = 1; % Grid index
    for i = 1:grid_cols % Grid column loop (1..NX)
        for j = 1:grid_rows  % Grid row loop (1..NZ)
            % Combine the contributions of the four edges:
            %   Top edge:    H_edge(j, i)
            %   Right edge:  V_edge(j, i+1)
            %   Bottom edge: H_edge(j+1, i) (negative)
            %   Left edge:   V_edge(j, i)   (negative)
            
            f = H_edge(j, i) ...      % Top edge
                + V_edge(j, i+1) ...  % Right edge
                - H_edge(j+1, i) ...  % Bottom edge (negative contribution)
                - V_edge(j, i);       % Left edge (negative contribution)
            
            Fxz(k, n) = f;
            n = n + 1;
        end
    end
end

% F*m = g
Gxz = Fxz * reshape(mod,[],1);

% Compute elapsed time
computation_time = toc;
fprintf('Computation of gravity kernel matrix Fxz completed, time is: %.6f seconds.\n', computation_time);
end
