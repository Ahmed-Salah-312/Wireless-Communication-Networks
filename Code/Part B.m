%%%%%%%%%%%%%%%%%%%%%%%%%PART_B%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%[1]Plot the cluster size versus SIR with range from 1dB : 30dB 
i0=[6,2,1];  
s=0; %counter for figures  
P=[]; %array of cluster size 
% Define the range of SIR_min values 
SIR_min_range_dB = 1:0.01:30; 
for i=1:length(i0)   
for j = 1:length(SIR_min_range_dB) 
[N , ~] = ClusterSize(SIR_min_range_dB(j), i0(i) , path_loss_exponent); 
P=[P,N]; %array of cluster size  
end 
%Plot cluster sizes for each no_sectors value 
figure(2+s); 
plot(SIR_min_range_dB, P, 'LineWidth', 2); 
xlabel('Minimum Signal-to-Interference Ratio (SIR_{min} dB)'); 
ylabel('Cluster Size'); 
title('Cluster Size vs. Minimum Signal-to-Interference Ratio for Different Number of Sectors'); 
grid on; 
hold on; 
% Set the figure's position and size 
set(gcf, 'Position', [x, y, width, height], 'Name', 'Cluster Size vs SIR'); 
P=[]; 
%figure(2+s); 
end  
legend('Omni-directional','120°sectorization','60° sectorization', 'Location', 'best');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%[2]for SIR_min = 19 dB & user density= 1400 
%[3]for SIR_min = 14 dB & user density= 1400 
X=[]; %number of cells 
Y=[]; % traffic intensity per cell  
GOS_Range=0.01:0.01:0.3; 
i=1;  
k=0; %counter for figures  
for SIR_min_dB=[14,19] 
for no_sectors=[1,3,6] 
[N, reuse_ratio] = ClusterSize(SIR_min_dB, i0(i) ,path_loss_exponent); 
for j=1:length(GOS_Range) 
[traffic_intensity_per_cell, traffic_intensity_per_sector] = TrafficIntensity(GOS_Range(j), no_sectors, total_channels, N); 
[num_cells] = NumCells(traffic_intensity_per_cell, 
traffic_intensity_per_user, user_density, city_area); 
X=[X,ceil(num_cells)]; 
Y=[Y,traffic_intensity_per_cell]; 
end 
figure(3+k); 
plot(GOS_Range,X, 'LineWidth', 2); 
title("Number of Cells vs GOS At SIR="+SIR_min_dB+"dB"); 
ylabel("Number Of Cells");  
xlabel("Grade Of Service(GOS)"); 
grid on; 
hold on; 
% Set the figure's position and size 
set(gcf, 'Position', [x, y, width, height],'Name',sprintf('Number of Cells vs GOS At %.d dB', SIR_min_dB)); 
figure(4+k); 
plot(GOS_Range,Y, 'LineWidth', 2); 
title("Traffic Intensity Per Cell vs GOS At SIR="+SIR_min_dB+"dB"); 
ylabel("Traffic Intensity Per Cell(Erlangs)"); 
xlabel("Grade Of Service(GOS)"); 
grid on; 
hold on; 
% Set the figure's position and size 
set(gcf, 'Position', [x, y, width, height], 'Name', sprintf('Traffic Intensity Per Cell vs GOS At %.d dB', SIR_min_dB)); 
X=[]; 
Y=[]; 
i=i+1; 
end 
figure(3+k); 
legend('Omni-directional','120°sectorization','60° sectorization', 'Location', 'best'); 
figure(4+k); 
legend('Omni-directional','120°sectorization','60° sectorization', 'Location', 'best'); 
k=k+2; 
i=1; 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  
%[4]for SIR_min = 19 dB & GOS= 2%  
%[5]for SIR_min = 14 dB &  GOS= 2%  
  
User_Density_Range=100:100:2000;   
M=[]; %number of cells 
L=[]; %  Cell_Radius  
m=1; % 
n=0; %counter for figures  
for SIR_min_dB=[14,19]  
    for no_sectors=[1,3,6] 
       [N, reuse_ratio] = ClusterSize(SIR_min_dB, i0(m) ,path_loss_exponent); 
         
       for j=1:length(User_Density_Range) 
           [traffic_intensity_per_cell, ~] = TrafficIntensity(GOS, no_sectors, total_channels, N); 
            [num_cells] = NumCells(traffic_intensity_per_cell, traffic_intensity_per_user, User_Density_Range(j), city_area); 
            num_cells=ceil(num_cells); 
            [R] = CellRadius(city_area, num_cells); 
            M=[M,ceil(num_cells)]; 
            L=[L,R];    
        end 
  
    figure(7+n); 
    plot(User_Density_Range,M, 'LineWidth', 2); 
    title("Number of Cells vs User Density At SIR="+SIR_min_dB+"dB"); 
    ylabel("Number Of Cells");  
    xlabel("User Density(users/km^2)"); 
    grid on; 
    hold on; 
     
    % Set the figure's position and size 
    set(gcf, 'Position', [x, y, width, height],'Name',sprintf('Number of Cells vs User_Density At %.d dB', SIR_min_dB)); 
  
    figure(8+n); 
    plot(User_Density_Range,L, 'LineWidth', 2); 
    title("Cell Radiusvs vs User Density At SIR="+SIR_min_dB+"dB"); 
    ylabel("Cell Radius(km)"); 
    xlabel("User Density(users/km^2)"); 
    grid on; 
    hold on; 
    % Set the figure's position and size 
     set(gcf, 'Position', [x, y, width, height], 'Name', sprintf('Cell Radius vs User Density At %.d dB', SIR_min_dB)); 
  
    M=[]; 
    L=[]; 
    m=m+1; 
     
    end 
     
    figure(7+n); 
    legend('Omni-directional','120°sectorization','60° sectorization', 'Location', 'best');
     figure(8+n); 
    legend('Omni-directional','120°sectorization','60° sectorization', 'Location', 'best'); 
    n=n+2; 
    m=1; 
end 
