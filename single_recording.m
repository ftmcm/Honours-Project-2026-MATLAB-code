% Calcium imaging analysis for a single recording session.
%
% This script processes Suite2p output from a single recording to:
%   1. Separate concatenated movies by voltage
%   2. Apply neuropil and photobleaching correction
%   3. Detect stimulus-evoked responses
%   4. Calculate response probability (AP proportion) per cell per voltage
%   5. Separate cells by distance from stimulation electrode
%   6. Generate input-output curves at single-cell and population level

%% Loading data
clear all
close all
load('suite2p/plane0/F.mat')
load('suite2p/plane0/Fneu.mat')
load('suite2p/plane0/iscell.mat')
load('suite2p/plane0/stat.mat')

%% User defined parameters
l_movie = 2999 
start_stim = 500 % frame at which stimulus starts
stim_interval = 500 % interval frames between stimuli
stimuli = [start_stim:stim_interval:l_movie-13];
n_movies = size(f_neuropil_true,2)/l_movie;
%% Filter to true cells only
true_cells = find(iscell(:,1)==1);
n_cells=size(true_cells,1);
f_true = F(true_cells,:);
f_neuropil_true = Fneu(true_cells,:);
%% Separate concatenated movies into individual movies
n_movies = round(n_movies);

for i = 1:n_movies
    start_frame = l_movie*(i-1)+1; %chop the movie
    final_frame = start_frame + l_movie-1
    win_this_movie=start_frame:final_frame;
    f_true_movie(:,:,i)=f_true(:,win_this_movie);    
    f_neuropil_true_movie(:,:,i) = f_neuropil_true(:,win_this_movie); 
end
%% Neuropil extraction
Fc_true = f_true_movie - 0.9*f_neuropil_true_movie;
dfof=Fc_true;
%% Extract voltages from file name 
% assumes that filename starts with voltage value e.g. 10V_recording.tif
folder = pwd;
files_tif = dir(fullfile(folder,'*.tif'));
files_tiff = dir(fullfile(folder,'*.tiff'));
abc = cat (1, files_tif, files_tiff);

fullname = {abc.name}';
voltage = sort(cellfun(@(s) sscanf(s, '%f'), fullname));
n_movies_vs_volt = [(1:n_movies).' voltage(:)]; 

this_folder = cd;
date_rec = regexp(this_folder,'(?<=/)\d+','match','once');
filename = sprintf('%s_voltage.mat', date_rec);
save(filename,"n_movies_vs_volt");

%% Z-score data
zdfof1 = zscore (dfof,0,2);
zdfof = smoothdata(zdfof1,"movmedian",20);

this_folder = cd;
date_rec = regexp(this_folder,'(?<=/)\d+','match','once');
filename = sprintf('%s_fluorescence_data',date_rec);
save(filename,"zdfof");

%% Plot activity for all cells
figure;

for i = 1:n_movies
    movie = dfof(:,:,i);
    movie_z = zscore(movie,0,2);
    for cell = 1:size (movie_z,1)
        subplot(1,9,i),plot(movie_z(cell,:)+10*cell); 
        hold on;
    end
    xline(stimuli);
    title(['Stimulus ',sprintf('(%d',voltage(i)),' V)'])
    xlabel('Frames (20Hz)')
    ylabel(['Cell (' sprintf('%d', length(true_cells)),'total)'])
end 

% Use this plot to visually inspect responses and determine which movie
% (voltage) shows the smallest clearly identifiable response.
% Set movie_min_resp in next section accordingly before running detection.

%% Detect response for stimulus
movie_min_resp = 1;

n_stimuli = length(stimuli);
AP_detect = NaN(length(true_cells),n_stimuli,n_movies);
col_stim = AP_detect;

% finds frame of max fluorescence 
for cell = 1: length(true_cells)
    for i = 1:n_stimuli
        for j = 1:n_movies
            response_win = stimuli(i)+1 : stimuli(i)+100;
            [X,col] = max(zdfof(cell, response_win, j));
            col_stim(cell, i, j) = response_win(1) + col - 1;
            AP_detect(cell,i,j) = 0;
        end
    end
end

% calculating deltaF_all; bsl_std_all
deltaF_all = NaN(length(true_cells),n_stimuli,n_movies);
bsl_std_all = NaN(length(true_cells),n_stimuli,n_movies);

for cell = 1:length(true_cells)
    for i = 1:n_stimuli
        for j = 1:n_movies
            col = col_stim(cell,i,j);
            start_F_max_win = max(col - 2, 1);
            end_F_max_win = min(col + 2, size(zdfof,2));
            F_max_win = start_F_max_win:end_F_max_win;
            F_max_mean = mean(zdfof(cell, F_max_win, j));
            F_bsl_win = stimuli(i)-80:stimuli(i)-5;
            F_bsl = mean(zdfof(cell,F_bsl_win,j));
            F_bsl_std = std(zdfof(cell,F_bsl_win,j));
            deltaF_all(cell, i, j) = F_max_mean - F_bsl;
            bsl_std_all(cell, i, j) = F_bsl_std;
        end
    end
end


% finding std and delta thresh
this_deltaF = deltaF_all(:,:,movie_min_resp);
median_deltaF = median(this_deltaF,'all');
[~, idx] = min(abs(this_deltaF - median_deltaF), [], 'all');
[row, stim_idx] = ind2sub(size(this_deltaF), idx);
F_bsl_win = stimuli(stim_idx)-80 : stimuli(stim_idx)-5;
F_bsl_std = std(zdfof(row, F_bsl_win, movie_min_resp));
std_thresh_og = median_deltaF / F_bsl_std;
delta_threshold = mean(median(deltaF_all(:,:,:), "all"));

%% Amplitude filter
deltaF_min_resp = mean(deltaF_all(:,:,movie_min_resp),"all");
amp_thresh = deltaF_min_resp/delta_threshold;

for cell = 1: length(true_cells)
    for i = 1:n_stimuli
        for j = 1:n_movies
                col = col_stim(cell,i,j);
                start_F_max_win = col - 3;
                end_F_max_win = col + 3;
                F_max_win = start_F_max_win:end_F_max_win;
                F_max_mean = mean (zdfof (cell, F_max_win, j));
                F_bsl_win = stimuli(i)-80:stimuli(i)-5;
                F_bsl = mean (zdfof (cell,F_bsl_win,j));
                if F_max_mean - F_bsl >= amp_thresh*delta_threshold
                    AP_detect(cell,i,j) = 1;
                end
            end 
        end
end


%% STD filter
% compares to bsl to remove noise
std_min_resp = mean(bsl_std_all(:,:,movie_min_resp),"all")
std_thresh =std_min_resp/std_thresh_og


for cell = 1: length(true_cells)
    for i = 1:n_stimuli
        for j = 1:n_movies
            if AP_detect (cell,i,j) == 1
                col = col_stim(cell,i,j);
                start_F_max_win = col - 2;
                end_F_max_win = col + 2;
                F_max_win = start_F_max_win:end_F_max_win;
                F_bsl_win = stimuli(i)-90:stimuli(i)-5;
                F_bsl = median(zdfof(cell,F_bsl_win,j));
                F_bsl_std = std(zdfof(cell,F_bsl_win,j));
                F_max_mean = mean (zdfof (cell, F_max_win, j));
                if F_max_mean - F_bsl < (std_thresh*F_bsl_std)
                    AP_detect(cell, i, j) = 0;
                end
            end
        end
    end
end
%% Calculation of AP proportion
AP_proportion = NaN(n_cells,n_movies);
for cell = 1:n_cells
    for i = 1:n_movies
       num_AP = sum(AP_detect(cell, :, i));
       AP_proportion(cell, i) = num_AP / length(stimuli);
    end
end

this_folder = cd;
date_rec = regexp(this_folder,'(?<=/)\d+','match','once');
filename = sprintf('%s_AP_proportion.mat', date_rec);

 save(filename,"AP_proportion");

%% Separating cells according to distance to electrode
% User has to determine these parameters
x_coord_electrode = max (x_y_coord_cell(:,2));
furthest_x_coord = min (x_y_coord_cell(:,2));

true_stat = stat(true_cells);
x_y_coord_cell = NaN(length(true_cells),3);
for cell = 1:length(true_cells)
    x_coord_cell = true_stat{1,cell}.med(1,2);
    y_coord_cell = true_stat{1,cell}.med(1,1);
    x_y_coord_cell (cell,:,:)= [cell,x_coord_cell,y_coord_cell];
end

size_field = round((x_coord_electrode-furthest_x_coord)/3);

far_cells = find (x_y_coord_cell(:,2)<=size_field+furthest_x_coord);
far_electrode = x_y_coord_cell (far_cells, :,:);

close_cells = find (x_y_coord_cell(:,2)>(2*size_field)+furthest_x_coord);
close_electrode = x_y_coord_cell (close_cells, :,:);

middle_cells = x_y_coord_cell(x_y_coord_cell(:,2) > size_field+furthest_x_coord & x_y_coord_cell(:,2) <= ((2*size_field)+furthest_x_coord));
middle_electrode = x_y_coord_cell (middle_cells, :,:);

this_folder = cd;
date_rec = regexp(this_folder,'(?<=/)\d+','match','once');

filename_close = sprintf('%s_close_electrode.mat', date_rec);
filename_middle = sprintf('%s_middle_electrode.mat', date_rec);
filename_far = sprintf('%s_far_electrode.mat', date_rec);

save(filename_close,'close_electrode');
save(filename_middle,'middle_electrode');
save(filename_far,'far_electrode');
         

%% Plot activity graph for specific cell
cell = 20;

for i = 1:n_movies
    movie = dfof(:,:,i);
    movie_z = zscore(movie,0,2);
        subplot(3,5,i),plot(zdfof(cell,:,i));
        %subplot(2,4,i),plot((smoothdata(movie_z(cell,:),"movmedian",20)));
        hold on;

    xline(stimuli);
    xlabel('Frames (20Hz)')
    ylabel(['Cell (' sprintf('%d', cell),')'])
    title(['Activity of cell at ' ,sprintf('%d', voltage(i)), ' V']) 
end 
%% Plot input-output graph for individual cell
cell = 20;
x = n_movies_vs_volt(:,2);
y = AP_proportion(cell,:);
figure;
plot_matrix = [x';y]';
plot(plot_matrix(:,1), plot_matrix(:,2),'color', [91 156 207]/256 ,'LineWidth',4);
    ylim ([-0.005 1.05]);
    tle = title(sprintf('Input-output curve for cell %d', cell(1,1)), 'color','white');
    xlabel('Voltage (V)')
    ylabel('Output (proportion of calcium transients)')
    ax = gca;
    box off
    ax.FontSize = 20;
    tle.FontSize = 30
    ax.XColor = 'white';
    ax.YColor = 'white';
    ax.Color = [34  34, 34]/256;
    set(gcf,'color',[34  34, 34]/256);
    ax.LineWidth = 4;
%% Plot input-output graph based on their distance
x = n_movies_vs_volt(:,2);
mean_y_close = mean(AP_proportion(close_cells,:));
mean_y_middle = mean(AP_proportion(middle_cells,:));
mean_y_far = mean(AP_proportion(far_cells,:));
blue = [91 156 207]/256;
red = [207 101 89]/256;
yellow = [207 179 89]/256
% gradient = [0.45  0.8 1.22];

figure;
hold on;
plot (x,mean_y_close,'color', blue,'LineWidth',4);
plot (x,mean_y_middle, 'color', red,'LineWidth',4);
plot (x,mean_y_far, 'color', yellow,'LineWidth',4);
% xline(voltage);
xlabel('Voltage (V)')
ylabel('Output (proportion of calcium transients)')
ax = gca;
box off
ylim ([-0.005 1.05]);
ax.FontSize = 20;
ax.XColor = 'white';
ax.YColor = 'white';
ax.Color = [34  34, 34]/256;
set(gcf,'color',[34  34, 34]/256);
ax.LineWidth = 4;
tle = title ('Input-output curve based on their distance from electrode', 'color','white')
    tle.FontSize = 30;
lgd = legend('Cells close to electrode','Cells at middle distance from electrode', 'Cells far from electrode');
    lgd.FontSize = 14;...
    lgd.LineWidth = 2;
    lgd.Color = 'black';
    lgd.EdgeColor = 'black'
axis square

%% Plot ideal response graph
cell = 15, j = 3

figure;
plot(zdfof(cell,420:999, j),'color', [91 156 207]/256 ,'LineWidth',2);
    xlabel('Frame')
    ylabel('Z-scored fluorescence')
    ax = gca;
    axis square
    box off
    ax.FontSize = 20;
    tle = title ('Example response to stimulus', 'Color','white');
        tle.FontSize = 30
    ax.XColor = 'white';
    ax.YColor = 'white';
    set(gcf,'color',[34  34, 34]/256);
    ax.Color = [34  34, 34]/256;
    ax.LineWidth = 4;

%% Plot activity for all voltages for cell
cell = 15;
for i = 1:n_movies
    movie = dfof(:,:,i);
    movie_z = zscore(movie,0,2);
        % subplot(1,n_movies,i),plot(zdfof(cell,:,i));
        subplot(n_movies/2,2,i),plot(zdfof(cell,:,i),'color', [91 156 207]/256 ,'LineWidth',2);
        hold on;
    xline(stimuli);
    xlabel('Frames (20Hz)')
    ylabel(['Cell (' sprintf('%d', cell),')'])
    ax = gca;
    ax.FontSize = 16;
    box off
    ax.XColor = 'white';
    ax.YColor = 'white';
    ax.Color = [34  34, 34]/256;
    set(gcf,'color',[34  34, 34]/256);
    ax.LineWidth = 4;
    tle = title(['Activity of cell at ' ,sprintf('%d', voltage(i)), ' V'],'Color','white')
        tle.FontSize = 25;
end  


%% Plot average response per stimulus

av_resp_zdfof = NaN(length(true_cells), size(zdfof,2));

for cell = 1: length(true_cells)
    for f= 1:size(zdfof,2)
        for j = 1:n_movies
            mean_f_value = mean(zdfof(cell,f,j));
            av_resp_zdfof(cell,f,j) = mean_f_value;
        end
    end
end


cell = 15; j = 3; 

figure;
hold on; 
plot(av_resp_zdfof(cell,400:999, j),'color', [91 156 207]/256 ,'LineWidth',4);
tle = title (['Average response per stimulus at ',sprintf('%d', voltage(j)), ' V'],'Color','white');
xlabel('Frame')
ylabel('Z-scored fluorescence')
ax = gca;
box off
ax.FontSize = 25;
ax.XColor = 'white';
ax.YColor = 'white';
ax.Color = [34  34, 34]/256;
set(gcf,'color',[34  34, 34]/256);
ax.LineWidth = 4;
tle.FontSize = 35;

%% Plot input-output graph for all cells
x = n_movies_vs_volt(:,2);
blue = [91 156 207]/256 
gradient = [0.2:1.03/length(true_cells):1.23];
gradient = gradient (1:length(true_cells));

figure;
hold on;
for cell = 1:length(true_cells)
    y = AP_proportion(cell,:);
    plot(x, y, 'LineWidth', 4);
    tle = title ('Input-output curve for all cells', 'color','white')
    xlabel('Voltage (V)')
    ylabel('Output (proportion of calcium transients)')
    ax = gca;
    box off
    axis square
    ylim ([-0.005 1.05]);
    ax.FontSize = 20;
    tle.FontSize = 30;
    ax.XColor = 'white';
    ax.YColor = 'white';
    ax.Color = [34  34, 34]/256;
    set(gcf,'color',[34  34, 34]/256);
    ax.LineWidth = 4;
end

%% Input-mean output graph
blue = [91 156 207]/256 

x = n_movies_vs_volt(:,2);
mean_y = mean(AP_proportion(:,:));
plot_matrix2 = [x';mean_y]';
figure;
hold on; 
plot (plot_matrix2(:,1),plot_matrix2(:,2),'color', [91 156 207]/256,'LineWidth', 4);
    tle = title ('Mean input-output curve', 'color', 'white')
    xlabel('Voltage (V)');
    ylabel('Mean output (proportion of calcium transients)');
    ax = gca;
    box off
    axis square
    ylim ([-0.005 1.05]);
    ax.FontSize = 20;
    tle.FontSize = 30;
    ax.XColor = 'white';
    ax.YColor = 'white';
    ax.Color = [34  34, 34]/256;
    set(gcf,'color',[34  34, 34]/256);
    ax.LineWidth = 4;