%%%%%%%%%%%%%%%%% constants %%%%%%%%%%%%%%%%%%%%% 
total_channels = 340; 
height_BS = 20; % meters 
height_MS = 1.5; % meters 
MS_sensitivity = -95; % dBm 
traffic_intensity_per_user = 0.025; % Erlang 
path_loss_exponent = 4; 
f = 900; % MHz 
% Set the figure's position and size  
x = 250;     y = 200;   width = 850;  height = 350;  
%%%%%%%%%%%%%%% input parameters %%%%%%%%%%%%%%%% 
GOS = input('Enter Grade of Service (GOS): '); 
city_area = input('Enter City Area (in square kilometers): '); 
user_density = input('Enter User Density (users per square kilometer): '); 
SIR_min_dB = input('Enter Minimum Signal-to-Interference Ratio (SIR_min in dB): '); 
sectorization_method_i0 = input('Enter Sectorization Method (1 for 60° 
sectorization, 2 for 120° sectorization, 6 for omni-directional): '); 
%%%%%%%%%%%%%%%%%%%%%% design parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Calculate Cluster Size (N)    
[Cluster_size, reuse_ratio] = ClusterSize(SIR_min_dB, sectorization_method_i0, 
path_loss_exponent); 
% Calculate Traffic intensity per cell, and traffic intensity per sector 
if sectorization_method_i0 == 1 
no_sectors = 6;   
elseif sectorization_method_i0 == 2 
no_sectors = 3; 
else     
no_sectors = 1; 
end 
[A_cell, A_sector] = TrafficIntensity(GOS, no_sectors, total_channels, Cluster_size); 
% Calculate Number of cells & Cell radius 
total_num_cells = NumCells(A_cell, traffic_intensity_per_user, user_density, 
city_area); 
R = CellRadius(city_area, total_num_cells); 
% calculate Base station transmitted power 
CH = 0.8 + (1.1 * log10(f) - 0.7) * height_MS - 1.56 * log10(f); 
PL = 69.55 + 26.16 * log10(f) - 13.82 * log10(height_BS) - CH + (44.9 - 6.55 * 
log10(height_BS)) * log10(R); 
BSpower_dBm = MS_sensitivity + PL; 
disp('Design Parameters:'); 
disp(['1) Cluster Size: ' num2str(Cluster_size)]); 
disp(['2) Number of Cells: ' num2str(total_num_cells)]); 
disp(['3) Cell Radius: ' num2str(R) ' Km']); 
disp(['4) Traffic Intensity per Cell: ' num2str(A_cell) ' Erlang']); 
disp(['5) Traffic Intensity per Sector: ' num2str(A_sector) ' Erlang']); 
disp(['6) Base Station Transmitted Power: ' num2str(BSpower_dBm) ' dBm']);

% plot for the MS received power in dBm versus the receiver distance from the BS 
d = linspace(0, R, 100); 
CH = 0.8 + (1.1 * log10(f) - 0.7) * height_MS - 1.56 * log10(f); 
PL2 = 69.55 + 26.16 * log10(f) - 13.82 * log10(height_BS) - CH + (44.9 - 6.55 * 
log10(height_BS)) * log10(d); 
RX_power_dBm = BSpower_dBm - PL2; 
figure(1); 
plot(d, RX_power_dBm , 'LineWidth', 2); 
xlabel('Distance from Base Station (km)'); 
ylabel('MS Received Power (dBm)'); 
title('MS Received Power vs. Distance from Base Station'); 
grid on; 
% Set the figure's position and size 
set(gcf, 'Position', [x, y, width, height], 'Name', 'MS Received Power vs. Distance 
from Base Station');

%%%%%%%%%%%%%%%%%% functions %%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% Calculate Cluster Size (N) %%%%%%%%%%%%%%%%%%%%%
function [N, reuse_ratio] = ClusterSize(SIR_min_dB, i0, n)
    SIRmin_ratio = 10^(SIR_min_dB/10); 
    reuse_ratio = (SIRmin_ratio * i0)^(1/n)+1; % SNR=(D/R)^n/io
    N = ceil((reuse_ratio^2)/3);
    
    % Discrete values of cluster size
    max_index = ceil(sqrt(N));      % N = i^2 if k = 0  
    m = [];
    for i = 0:max_index
        for k = i:max_index
            m(end + 1) = i^2 + k^2 + i*k;
        end
    end
    z = sort(unique(m));
    zz = z;
    for j = 1:length(N)
        N(j) = min(zz(zz >= N(j)));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%% Calculate Traffic Intensity %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function [traffic_intensity_per_cell, traffic_intensity_per_sector] = TrafficIntensity(GOS, no_sectors, no_channels, N)
    K = floor(no_channels / N); % channels per cell
    C = floor(K / no_sectors);

  % Initial guess for fzero
    for A1 = 1:1000 % for loop to get close value to traffic intensity to use   it in fzero
       
        Pr = ( A1 ^ C / factorial ( C )) / sum ( A1 .^([0: C ]) ./ cumprod([0 ,(0: C -1)]+1) ) ;

        if GOS <= Pr
            break
       end
    end
    % Define the Erlang function
    Erlang = @(A) (A^C / factorial(C)) / sum(A.^((0:C)) ./ cumprod([1, (1:C)]));

    % Solve for traffic intensity using fzero
    traffic_intensity_per_sector = fzero(@(A) Erlang(A) - GOS, A1);

    % Calculate traffic intensity per cell
    traffic_intensity_per_cell = traffic_intensity_per_sector * no_sectors;
end

%%%%%%%%%%%%%%%%%%%%% num of cells %%%%%%%%%%%%%%%%%%%%
function [num_cells] = NumCells(traffic_intensity_per_cell, traffic_intensity_per_user, user_density, city_area)
    U = traffic_intensity_per_cell / traffic_intensity_per_user;  % number of users 
    cell_area = U / user_density;
    num_cells = ceil(city_area / cell_area);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%  Cell Radius (R) %%%%%%%%%%%%%%%%%%%%%
function [R] = CellRadius(city_area, num_cells)
    R = sqrt(city_area / (3/2 * sqrt(3) * num_cells)); % kilometers 
end
