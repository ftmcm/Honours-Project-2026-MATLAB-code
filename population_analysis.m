%population_analysis
%% Loading data
clear all
close all
files = dir('*.mat');
allData = cell(size(files));

for ii=1:length(files)
  aFileData = load(files(ii).name);
  varName = fieldnames(aFileData);
  allData{ii} = aFileData.(varName{1});
end

clear ii; clear aFileData; clear varName;

n_files_per_recording = 9; 
n_recordings = length(files)/n_files_per_recording;

recordings = cell(n_recordings, n_files_per_recording);

for r = 1:n_recordings
    file_idx = (r-1)*n_files_per_recording + 1 : r*n_files_per_recording;
    recordings(r,:) = allData(file_idx);
end

% Order of files in folder, obtained from single_recording_analysis, check
% variable files to ensure correct order
AP_proportion = recordings(:,1);
close_electrode = recordings(:,2);
far_electrode = recordings(:,3);
fluorescence = recordings(:,4);
iscell = recordings(:,5);
middle_electrode = recordings(:,6);
F = recordings(:,7);
Fneu = recordings(:,8);
voltage = recordings(:,9);


voltage_value = cell(1, n_recordings);

for i = 1:n_recordings
    voltage_value{i} = voltage{i}(:,2);
end

voltage = voltage_value'
voltage_value = voltage_value'

n_cells = cellfun('size',AP_proportion,1);
n_voltages_array = cellfun('size',voltage,1); %number voltages per recording

% find number of unique voltage values
volt_vector = [];
start = 1;

for i = 1:n_recordings
    ending = start + length(voltage{i}(:)) - 1;
    volt_vector(start:ending,1) = voltage{i}(:);
    start = ending + 1;
end

volt_vector = unique(volt_vector)
n_volt = size(volt_vector,1) % number volts
max_n_volt = max(n_voltages_array)% max number volts in a recording

voltage_all = NaN(max_n_volt,n_recordings);

for r = 1:n_recordings
    for row = 1:max_n_volt
        if size(voltage_value{r,1},1) >= row
            voltage_all(row,r) = voltage_value{r,1}(row,1);
        end
    end
end

this_folder = pwd;
[~, folder_name] = fileparts(this_folder);

 %% Sorting data into usable arrays 

list_volt_values = unique(voltage_all);
sorted_voltage_all = NaN(n_volt,n_recordings);



for r = 1:n_recordings
    for row = 1:max_n_volt
        for volt = 1:n_volt
            if voltage_all(row,r) == list_volt_values(volt)
                sorted_voltage_all(volt,r) = list_volt_values(volt);
            end
        end
    end
end


sorted_AP_all = NaN(max(n_cells), n_volt, n_recordings);

for r = 1:n_recordings
    for cell = 1:n_cells(r)
        for col = 1:max_n_volt
            for volt = 1:n_volt
                if voltage_all(col,r) == list_volt_values(volt)
                    sorted_AP_all(cell,volt,r) = AP_proportion{r,1}(cell,col);
                end
            end
        end
    end
end

max_cells  = max(cellfun('size', fluorescence, 1));
max_frames = max(cellfun('size', fluorescence, 2));
max_movies = max(cellfun('size', fluorescence, 3));
n_frames = cellfun('size', fluorescence, 2);
n_movies = cellfun('size', fluorescence, 3);

sorted_fluorescence_all = NaN(max_cells, max_frames, n_volt, n_recordings);

for r = 1:n_recordings
    for cell = 1:n_cells(r)
        for frame = 1:n_frames(r)
            for movie = 1:n_movies(r)
                for volt = 1:n_volt
                    if voltage_all(movie,r) == list_volt_values(volt)
                        sorted_fluorescence_all(cell,frame, volt, r) = fluorescence{r,1}(cell,frame,movie);
                    end
                end
            end 
        end
    end
end


sorted_cells_all = NaN(max(n_cells),n_recordings);
for r = 1:n_recordings
    for cell = 1:n_cells(r)
        sorted_cells_all(cell,r) = cell;
    end
end
%% Sort iscell, F and Fneu, keeping only true cells (iscell(:,1) == 1)

F_all = F;
Fneu_all = Fneu;
iscell_all = iscell;
max_frames = max(cellfun('size', F_all, 2));

% sorted arrays
conc_sorted_F_all      = NaN(size(sorted_fluorescence_all,1), max_frames,n_recordings);
conc_sorted_Fneu_all   = NaN(size(sorted_fluorescence_all,1), max_frames,n_recordings);


for r = 1:n_recordings
    true_cells = iscell_all {r,1}(:,1) == 1;
     F_true = F_all{r,1}(true_cells,:);
     Fneu_true = Fneu_all{r,1}(true_cells,:);
     n_true = size(F_true,1);
     n_f = size(F_true,2);
     conc_sorted_F_all (1:n_true,1:n_f,r) = F_true;
     conc_sorted_Fneu_all (1:n_true,1:n_f,r) = Fneu_true;
end 

size_frames_recording = cellfun('size', fluorescence, 2);
conc_n_frames = cellfun('size', F_all, 2);
n_frames = cellfun('size', fluorescence, 2);

sorted_F_all   = NaN(size(sorted_fluorescence_all));
sorted_Fneu_all = NaN(size(sorted_fluorescence_all));

for r = 1:n_recordings
    n_frames_current = n_frames(r);
    for c = 1:n_cells(r)
        for movie = 1:n_movies(r)
            start_idx = (movie-1)*n_frames_current + 1;
            for volt = 1:n_volt
                if voltage_all(movie,r) == list_volt_values(volt)
                    for f = 1:n_frames_current
                        frame = start_idx + f - 1;
                        sorted_F_all(c,f,volt,r)    = conc_sorted_F_all(c, frame, r);
                        sorted_Fneu_all(c,f,volt,r) = conc_sorted_Fneu_all(c, frame, r);
                    end
                end
            end
        end
    end
end

%% Mean neuropil trace as bleaching estimate
stim_info = [ % [start frame where first stimuli was applied, frame interval between stimuli]
    200   300;
    300   300;
    500   500;
    500   500;
    500   500;
    500   500;
    500   500;
];


mean_neu    = mean(sorted_Fneu_all, 1, "omitnan");
sq_mean_neu = squeeze(mean_neu);   % [frames x volt x recordings]

% Build stim_info_matrix
n_frames  = cellfun('size', fluorescence, 2);

stim_info_matrix = NaN(n_recordings, 5); 

for r = 1:n_recordings
    stim_start = stim_info(r,1):stim_info(r,2):n_frames(r)-1;
    stim_info_matrix(r,1:length(stim_start)) = stim_start;
end
% 
% Compute bleaching-corrected fluorescence
% Formula: (neu_trace / mean(neu_trace[1:10])) * mean(cell_trace[1:10])
n_cells    = cellfun('size', AP_proportion, 1);
correct_Fc = NaN(size(sorted_fluorescence_all));

for r = 1:n_recordings
    for v = 1:n_volt
        neu_trace = sq_mean_neu(:, v, r);
        neu_early_mean = mean(neu_trace(1:10), 'omitnan');
        if isnan(neu_early_mean) || neu_early_mean == 0
            continue
        end
        neu_normalised = neu_trace / neu_early_mean;
        for c = 1:n_cells(r)
            cell_trace      = squeeze(sorted_fluorescence_all(c, :, v, r))';
            cell_early_mean = mean(cell_trace(1:10), 'omitnan');
            if isnan(cell_early_mean)
                continue
            end
            bleaching_estimate = neu_normalised * cell_early_mean;
            correct_Fc(c, :, v, r) = cell_trace - 0.9*bleaching_estimate;
        end
    end
end
%% AP detection
min_resp_volt = [5   6    20    10    25    20    20] % User defined from inspecting fluorescence traces                                                      % across all movies
n_stimuli = [6 6 5 5 5 5 5 ] % Used defined; number of stimuli on each movie
max_stimuli = 6;


for r = 1:n_recordings
    %new_sorted_AP_all(1:n_cells(r), :, r) = AP_prop_r;
    Fc_r      = squeeze(correct_Fc(1:n_cells(r), :, :, r));
    n_cells_r = n_cells(r);
    stimuli_r = stim_info_matrix(r, 1:n_stimuli(r));
    stimuli_r = stimuli_r(~isnan(stimuli_r));  
    n_stim_r  = n_stimuli(r)

    valid_volt = find(squeeze(any(any(~isnan(correct_Fc(1:n_cells(r), :, :, r)), 1), 2))');
    movie_min_resp = find(sorted_voltage_all(:,r) == min_resp_volt(r), 1);

    AP_detect   = NaN(n_cells_r, max_stimuli, n_volt);
    col_stim    = NaN(n_cells_r, max_stimuli, n_volt);
    deltaF_all  = NaN(n_cells_r, max_stimuli, n_volt);
    bsl_std_all = NaN(n_cells_r, max_stimuli, n_volt);

    % Find peak frame and baseline stats
    for c = 1:n_cells_r
        for i = 1:n_stim_r
            for v = valid_volt
                response_win = stimuli_r(i)+1 : min(stimuli_r(i)+80, size(Fc_r,2));
                this_trace   = Fc_r(c, response_win, v);
                [~, col] = max(this_trace);

                col_stim(c,i,v)  = response_win(1) + col - 1;
                AP_detect(c,i,v) = 0;

                col_now    = col_stim(c,i,v);
                F_max_win  = max(col_now-2,1) : min(col_now+2, size(Fc_r,2));
                F_bsl_win  = max(stimuli_r(i)-80,1) : max(stimuli_r(i)-10,1);

                F_max_mean = mean(Fc_r(c, F_max_win, v), 'omitnan');
                F_bsl      = mean(Fc_r(c, F_bsl_win, v), 'omitnan');
                F_bsl_std  = std(Fc_r(c, F_bsl_win, v), 0, 'omitnan');

                deltaF_all(c,i,v)  = F_max_mean - F_bsl;
                bsl_std_all(c,i,v) = F_bsl_std;
            end
        end
    end

    
    % thresholds
    this_deltaF     = deltaF_all(:,1:n_stim_r,movie_min_resp);
    median_deltaF   = median(this_deltaF, 'all', 'omitnan');
    [~, idx]        = min(abs(this_deltaF - median_deltaF), [], 'all');
    [row, stim_idx] = ind2sub(size(this_deltaF), idx);

    F_bsl_win_th    = max(stimuli_r(stim_idx)-80,1) : max(stimuli_r(stim_idx)-5,1);
    F_bsl_std_th    = std(Fc_r(row, F_bsl_win_th, movie_min_resp), 0, 'omitnan');
    std_thresh_og   = median_deltaF / F_bsl_std_th;

    delta_threshold = median(deltaF_all(:,1:n_stim_r,movie_min_resp), 'all', 'omitnan');
    deltaF_min_resp = mean(deltaF_all(:,1:n_stim_r,movie_min_resp), 'all', 'omitnan');
    amp_thresh      = deltaF_min_resp / delta_threshold;


    fprintf('Recording %d | amp_thresh: %.4f | delta_threshold: %.4f | amp_thresh*delta: %.4f | std_thresh: %.4f\n', ...
            r, amp_thresh, delta_threshold, amp_thresh*delta_threshold, std_thresh);

    % Amplitude filter
    for c = 1:n_cells_r
        for i = 1:n_stim_r
            for v = valid_volt
                col_now = col_stim(c,i,v);
                if isnan(col_now), continue, end

                F_max_win  = max(col_now-3,1) : min(col_now+3, size(Fc_r,2));
                F_bsl_win  = max(stimuli_r(i)-80,1) : max(stimuli_r(i)-5,1);

                F_max_mean = mean(Fc_r(c, F_max_win, v), 'omitnan');
                F_bsl      = mean(Fc_r(c, F_bsl_win, v), 'omitnan');

                if F_max_mean - F_bsl >= amp_thresh * delta_threshold
                    AP_detect(c,i,v) = 1;
                end
            end
        end
    end

    % STD filter
    for c = 1:n_cells_r
        for i = 1:n_stim_r
            for v = valid_volt
                if AP_detect(c,i,v) == 1
                    col_now = col_stim(c,i,v);
                    if isnan(col_now), continue, end

                    F_max_win  = max(col_now-2,1) : min(col_now+2, size(Fc_r,2));
                    F_bsl_win  = max(stimuli_r(i)-90,1) : max(stimuli_r(i)-5,1);

                    F_bsl      = median(Fc_r(c, F_bsl_win, v), 'omitnan');
                    F_bsl_std  = std(Fc_r(c, F_bsl_win, v), 0, 'omitnan');
                    F_max_mean = mean(Fc_r(c, F_max_win, v), 'omitnan');

                    if F_max_mean - F_bsl < std_thresh * F_bsl_std
                        AP_detect(c,i,v) = 0;
                    end
                end
            end
        end
    end

    % AP proportion [cells x volt]
    AP_prop_r = NaN(n_cells_r, n_volt);
    for c = 1:n_cells_r
        for v = valid_volt
            AP_prop_r(c,v) = sum(AP_detect(c,1:n_stim_r,v) == 1) / n_stim_r;
        end
    end
    new_sorted_AP_all(1:n_cells(r), :, r) = AP_prop_r;
end
% Important to compare code output with visually detectable responses. If
% they differ, own thresholds can be defined by user by adding the
% following onto the loop:
    % if r == X
    %   amp_thresh = Y
    %   std_thresh = Z
    % end 
 
%% Plots individual cell's activity and input-output curve
recording = 1;
cell = 10;
row = find(~isnan(sorted_voltage_all(:,recording)));
grey = [150 150 150]/256;
figure;
set(gcf, 'Position', [100 100 287 681.9]);
set(gcf, 'color', 'white');
hold on;
offset = 3;   % try 20, 50, or even 100
k = 1;
for i = row'
trace = squeeze(correct_Fc(cell,:,i,recording));
plot(smoothdata(trace,"movmean",30) + offset*k, 'color', grey, 'LineWidth', 3);
k = k + 1;
end
box off
ax = gca;
ax.FontSize = 20;
ax.XColor = 'white';
ax.YColor = 'white';
ax.Color = 'white';
ax.LineWidth = 2;

blue  = [0 70 173]/255;     
red   = [200 30 30]/255;    
green = [0 140 60]/255; 

valid = ~isnan(sorted_voltage_all(:,recording));

x = sorted_voltage_all(valid, recording);
y = squeeze(new_sorted_AP_all(cell, valid, recording));

figure;
hold on;

plot(x, y, ...
    'Color', red, ...
    'LineWidth', 4); 

% Axes styling
xlabel('Voltage (V)');
ylabel('Probability of calcium transients');

ylim([0 1.05]);
xlim([0 max(x)+5]);

ax = gca;
set(gcf,'color','white');
ax.Color = 'white';
ax.FontSize = 20;
ax.LineWidth = 2;
ax.XColor = 'black';
ax.YColor = 'black';

box off  

%% Plot all individual cells + population mean
% Build plot_matrix for all cells using new_sorted_AP_all

rows = size(sorted_voltage_all, 1);
total_cells = sum(n_cells);

plot_matrix = NaN(rows, 2, total_cells);

global_cell = 0;

for r = 1:n_recordings
    for c = 1:n_cells(r)
        global_cell = global_cell + 1;

        % x values = voltage
        plot_matrix(:,1,global_cell) = sorted_voltage_all(:,r);

        % y values = response for this cell
        plot_matrix(:,2,global_cell) = new_sorted_AP_all(c,:,r)';
    end
end


grey = [150 150 150]/256;

figure;
hold on;
set(gcf, 'Position', [100 100 933 600]);
set(gcf, 'color', 'white');

cols = lines(total_cells);   % different colour for each cell

for c = 1:total_cells
    x_cell = plot_matrix(:,1,c);
    y_cell = plot_matrix(:,2,c);

    valid = ~isnan(x_cell) & ~isnan(y_cell);

    lw = 0.5 + 1.5*rand;   % random line thickness
    fade = 0.35;           % makes colours a bit softer
    basecol = cols(c,:);
    col = basecol + (1 - basecol)*fade;

    plot(x_cell(valid), y_cell(valid), ...
        'Color', col, ...
        'LineWidth', lw);
end

mean_x = mean(plot_matrix(:,1,:), 3, 'omitnan');
mean_y = mean(plot_matrix(:,2,:), 3, 'omitnan');

plot(mean_x, mean_y, 'Color', 0.1*grey, 'LineWidth', 4);

xlabel('Voltage (V)');
ylabel('Probability of calcium transients');
ylim([-0.005 1.19]);

ax = gca;
box off
ax.FontSize = 20;
ax.XColor = 'black';
ax.YColor = 'black';
ax.Color = 'white';
ax.LineWidth = 2;

h1 = plot(nan, nan, 'Color', [0.75 0.75 0.75], 'LineWidth', 1.5);
h2 = plot(nan, nan, 'Color', 0.1*grey, 'LineWidth', 5);

lgd = legend([h1 h2], {'Individual cells', 'Population mean'}, ...
    'Location', 'northeast');

lgd.FontSize = 14;
lgd.TextColor = 'black';
lgd.LineWidth = 1.5;
lgd.Color = 'white';
lgd.EdgeColor = 'white';
plot (x,y)
%% Separates cells into distance-based groups
n_cells_close = cellfun('size', close_electrode, 1);
total_cells_close = sum(n_cells_close);
sorted_close_electrode = NaN(max(n_cells_close), 3, n_recordings);

for r = 1:n_recordings
    for row = 1:n_cells_close(r)
        for col = 1:3
            sorted_close_electrode(row,col,r) = close_electrode{r,1}(row,col);
        end
    end
end

n_cells_middle = cellfun('size', middle_electrode, 1);
total_cells_middle = sum(n_cells_middle);
sorted_middle_electrode = NaN(max(n_cells_middle), 3, n_recordings);

for r = 1:n_recordings
    for row = 1:n_cells_middle(r)
        for col = 1:3
            sorted_middle_electrode(row,col,r) = middle_electrode{r,1}(row,col);
        end
    end
end

n_cells_far = cellfun('size', far_electrode, 1);
total_cells_far = sum(n_cells_far);
sorted_far_electrode = NaN(max(n_cells_far), 3, n_recordings);

for r = 1:n_recordings
    for row = 1:n_cells_far(r)
        for col = 1:3
            sorted_far_electrode(row,col,r) = far_electrode{r,1}(row,col);
        end
    end
end
%% Calculates mean AP proportion for each distance group
mean_x = mean(plot_matrix(:,1,:), 3, "omitnan");

y_close_all = [];

for r = 1:n_recordings
    close_cells = sorted_close_electrode(:,1,r);
    close_cells = close_cells(~isnan(close_cells));
    y_close_all = [y_close_all; new_sorted_AP_all(close_cells,:,r)];
end

mean_y_close = mean(y_close_all, 1, "omitnan");

y_middle_all = []
for r = 1:n_recordings
    middle_cells = sorted_middle_electrode(:,1,r);
    middle_cells = middle_cells(~isnan(middle_cells));
    y_middle_all = [y_middle_all; new_sorted_AP_all(middle_cells,:,r)];
end
mean_y_middle = mean(y_middle_all, 1, "omitnan");


y_far_all = []
for r = 1:n_recordings
    far_cells = sorted_far_electrode(:,1,r);
    far_cells = far_cells(~isnan(far_cells));
    y_far_all = [y_far_all; new_sorted_AP_all(far_cells,:,r)];
end
mean_y_far = mean(y_far_all, 1, "omitnan");

%% Mean input output graph based on distance
blue  = [0 70 173]/255;     
red   = [200 30 30]/255;    
green = [0 140 60]/255; 

figure;
hold on;
plot (mean_x,mean_y_close,'color', blue,'LineWidth',4);
plot (mean_x,mean_y_middle, 'color', red,'LineWidth',4);
plot (mean_x,mean_y_far, 'color', green,'LineWidth',4);
xlabel('Voltage (V)')
ylabel('Probability of calcium transients')
ax = gca;
box off
 ylim ([-0.005 1.18]);
    ax.FontSize = 20;
    ax.XColor = 'black';
    ax.YColor = 'black';
    ax.Color = 'white';
    set(gcf,'color', 'white');
    ax.LineWidth = 2;
tle = title ('Input-output curve based on their distance from electrode', 'color','white')
    tle.FontSize = 30;
lgd = legend('Cells close to electrode','Cells at middle distance from electrode', 'Cells far from electrode');
    lgd.Location = "northwest"
    lgd.FontSize = 14;
    lgd.TextColor = 'black';
    lgd.LineWidth = 2;
    lgd.Color = 'white';
    lgd.EdgeColor = 'white'

    fig_file_name = sprintf('%s_distance_input_output', folder_name);
    savefig(fig_file_name) 


figure;
hold on;
plot (mean_x,mean_y_close,'-o', 'MarkerSize', 7, 'MarkerFaceColor',blue,'color', blue,'LineWidth',4);
plot (mean_x,mean_y_middle, '-o', 'MarkerSize', 7, 'MarkerFaceColor',red,'color', red,'LineWidth',4);
plot (mean_x,mean_y_far, '-o', 'MarkerSize', 7, 'MarkerFaceColor',green,'color', green,'LineWidth',4);
    

xlabel('Voltage (V)')
ylabel('Probability of calcium transients')
ax = gca;
box off
 ylim ([-0.005 1.18]);
    ax.FontSize = 20;
    ax.XColor = 'black';
    ax.YColor = 'black';
    ax.Color = 'white';
    set(gcf,'color', 'white');
    ax.LineWidth = 2;
tle = title ('Input-output curve based on their distance from electrode', 'color','white')
    tle.FontSize = 30;
lgd = legend('Cells close to electrode','Cells at middle distance from electrode', 'Cells far from electrode');
    lgd.Location = "northwest"
    lgd.FontSize = 14;
    lgd.TextColor = 'black';
    lgd.LineWidth = 2;
    lgd.Color = 'white';
    lgd.EdgeColor = 'white'

fig_file_name = sprintf('%s_distance_input_output_dots', folder_name);
    savefig(fig_file_name) 
%% Input output contrasting brain regions
blue  = [0 70 173]/255;     
red   = [200 30 30]/255;    
green = [0 140 60]/255; 

figure;
hold on; 
plot(LH_x, LH_y,'color', blue,'LineWidth',4)
plot(hor_x, hor_y,'color', red,'LineWidth',4)
plot(pons_mean_x, pons_mean_y,'color', green,'LineWidth',4)
    xlim ([0 40])
    xlabel('Voltage (V)');
    ylabel('Probability of calcium transients');
    ylim([-0.005 1.19]);
    ax = gca;
    set(gcf,'color', 'white');
    box off
    ax.FontSize = 20;
    ax.XColor = 'black';
    ax.YColor = 'black';
    ax.Color = 'white';
    ax.LineWidth = 2;

lgd = legend('Coronal',"Horizontal");
    lgd.Location = "northwest"
    lgd.FontSize = 14;
    lgd.TextColor = 'black';
    lgd.LineWidth = 2;
    lgd.Color = 'white';
    lgd.EdgeColor = 'white'
fig_file_name = '0V_40V_coronal_vs_horizontal';
savefig(fig_file_name) 

   
%% Calculating sigmoid formula for each cell
%
b_close   = NaN(size(y_close_all,1),1);
c_close   = NaN(size(y_close_all,1),1);
b_middle  = NaN(size(y_middle_all,1),1);
c_middle  = NaN(size(y_middle_all,1),1);
b_far     = NaN(size(y_far_all,1),1);
c_far     = NaN(size(y_far_all,1),1);



s = fittype('1/(1+exp(-b*(x-c)))', 'Coefficients',{'b','c'}, 'Independent','x', 'Dependent','y');
x = mean(sorted_voltage_all, 2, "omitnan")';
opts = fitoptions(s);
opts.StartPoint = [0.1, 1]; 
opts.Lower = [0, min(x)];
opts.Upper = [Inf, max(x)];


for cell = 1:size(y_close_all,1)
    y = y_close_all(cell,:);
    valid = ~isnan(x) & ~isnan(y);
    xv = x(valid)';
    yv = y(valid)';
    sfit =  fit(xv, yv, s, opts);
    b_close (cell) = sfit.b;
    c_close (cell) = sfit.c;
end

for cell = 1:size(y_middle_all,1)
    y = y_middle_all(cell,:);
    valid = ~isnan(x) & ~isnan(y);
    xv = x(valid)';
    yv = y(valid)';
    sfit =  fit(xv, yv, s, opts);
    b_middle (cell) = sfit.b;
    c_middle (cell) = sfit.c;
end

for cell = 1:size(y_far_all,1)
    y = y_far_all(cell,:);
    valid = ~isnan(x) & ~isnan(y);
    xv = x(valid)';
    yv = y(valid)';
    sfit =  fit(xv, yv, s, opts);
    b_far (cell) = sfit.b;
    c_far (cell) = sfit.c;
end

mean_y = mean(plot_matrix(:,2,:), 3, "omitnan");
mean_x = mean(plot_matrix(:,1,:), 3, "omitnan");

x = mean_x(:);
y = mean_y(:);

valid = ~isnan(x) & ~isnan(y);
xv = x(valid);
yv = y(valid);

sfit = fit(xv, yv, s, opts);
b_coronal = sfit.b;
c_coronal = sfit.c;
%% Statistical test and bar chart
[h,p,ci]=ttest2(c_close, c_far,"Tail","left")
[h,p,ci]=ttest2(c_close, c_far)

% bar chart for c-values
mean_c_close = mean(c_close);
mean_c_far = mean(c_far);

total_close = sum(n_cells_close)
total_far = sum(n_cells_far)

std_close = std(c_close);
std_far   = std(c_far);

sem_close = std_close / sqrt(total_close);
sem_far   = std_far / sqrt(total_far);

y = [mean_c_close mean_c_far]
figure
set(gcf, 'color', 'w')  

b = bar(y);
y_max = max(y + [sem_close sem_far]); 
b.FaceColor = [0.7 0.7 0.7];   
b.EdgeColor = 'k';             
b.LineWidth = 1.5;

hold on

errorbar([1 2], y, [sem_close sem_far], 'k', 'LineStyle', 'none', 'LineWidth', 1.5)  

ax = gca;
ax.Color = 'w';         
ax.LineWidth = 1.5;
ax.FontSize = 20;
ax.XColor = 'black';
ax.YColor = 'black';

xticks([1 2])
xticklabels({'Close','Far'})

ylabel('Voltage at 50% response', 'color', 'black')
xlabel('Distance group', 'color', 'black')

box off

