clear; close all;

addpath('C:\Users\smellor\Documents\GitHub\optitrack');

addpath('C:\Users\smellor\Documents\GitHub\spm-856cd60354d5a789b73e67de9e1faeb4163cf8fa\');
spm('defaults', 'eeg');

addpath('C:\Users\smellor\Documents\GitHub\BrewerMap')
colormap123 = colormap(flipud(brewermap(64,'RdBu')));
addpath('C:\Users\smellor\Documents\GitHub\icp');

addpath('C:\Users\smellor\Documents\GitHub\linspecer');

addpath('C:\Users\smellor\Documents\GitHub\MEGsurfer');

%% Format meta data to do analysis

delay = 10; % ms - neuro-1 delay between truth and recording

cd('E:\Data\Neuro1\Auditory\anonymised_for_sharing');

subIDs = {'sub-001', 'sub-002', 'sub-003'};

meta_data = table('Size', [0,8], 'VariableTypes', {'string', 'string', 'string', 'string', 'string', 'string', 'string', 'string'},...
    'VariableNames', {'sub', 'raw_data_name', 'raw_data_loc', 'analysed_data_loc', 'results_save_loc', 'OptiTrig', 'AudioTrig', 'stim_data_fname'});

% Format meta data table
for sub = 1:length(subIDs)

    % Set audio and optitrack trigger names
    if strcmp(subIDs{sub}, 'sub-001')
        optitrig = 'AI8';
        audiotrig = 'AI16';
    elseif strcmp(subIDs{sub}, 'sub-002')
        optitrig = 'AI16';
        audiotrig = 'AI8';
    elseif strcmp(subIDs{sub}, 'sub-003')
        optitrig = 'T3';
        audiotrig = 'A16';
    else
        error('Please set trigger channel names for participant %s', subIDs{sub})
    end

    % Search folders for files
    fpathRoot2 = fullfile(cd, 'rawData');

    % Find data files
    fpathRoot3 = fullfile(fpathRoot2, subIDs{sub}, 'meg');

    % Find raw data name
    listing = dir(fpathRoot3);
    listing = extractfield(listing, 'name');

    % Just choose lvm files
    lvmfiles = listing(endsWith(listing, '.lvm'));
            
    for bb = 1:length(lvmfiles)
        % Find corresponding stim data file
        fpath_stim = fullfile(fpathRoot2, subIDs{sub}, 'stim');
        stimfile = strrep(lvmfiles{bb}, '_meg.lvm', '_stim.csv');

        % Fill in table
        meta_data = [meta_data; {subIDs{sub}, lvmfiles{bb}, fpathRoot3, ...
            fullfile(cd, 'analysedData', subIDs{sub}), ...
            fullfile(cd, 'results', subIDs{sub}), optitrig, audiotrig, fullfile(fpath_stim, stimfile)}];
    end
end


clearvars -except meta_data colormap123 delay


%% Load data

for recording = 1:size(meta_data,1)

    if strcmp(meta_data{recording, "sub"}, 'sub-003')
        rad_ax = 'Z';
    else
        rad_ax = 'Y';
    end

    rawDataPath = char(meta_data{recording, "raw_data_loc"});
    analysedDataPath = char(meta_data{recording, "analysed_data_loc"});
    resultsPath = char(meta_data{recording, "results_save_loc"});

    if ~exist(analysedDataPath, "dir")
        mkdir(analysedDataPath);
    end
    if ~exist(resultsPath, "dir")
        mkdir(resultsPath);
    end
    
    cd(rawDataPath)
    
    fname = char(extractBefore(meta_data{recording, "raw_data_name"}, '_meg.lvm'));
    
    if isfile(fullfile(analysedDataPath, [fname, '_meg.mat']))
        D = spm_eeg_load(fullfile(analysedDataPath, [fname, '_meg.mat']));
    else
        S = [];
        S.positions = sprintf('%s_positions.tsv', fname);
        S.data = fullfile(cd, sprintf('%s_meg.lvm', fname));
        S.path = analysedDataPath;
        S.sMRI = fullfile(char(extractBefore(meta_data{recording, "raw_data_loc"},'meg')),...
            'anat', sprintf('%s.nii', meta_data{recording, 'sub'}));
        ctx_fname = fullfile(char(extractBefore(meta_data{recording, "raw_data_loc"},'meg')),...
            'anat', sprintf('%s.L.midthickness.4k_fs_LR.surf.gii', meta_data{recording, 'sub'}));
        combine_surfaces({ctx_fname, strrep(ctx_fname, '.L.', '.R.')}, ...
            fullfile(char(extractBefore(meta_data{recording, "raw_data_loc"},'meg')),...
            'anat', sprintf('%s_midthickness_cortex.gii', meta_data{recording, 'sub'})));
        S.cortex = fullfile(char(extractBefore(meta_data{recording, "raw_data_loc"},'meg')),...
            'anat', sprintf('%s_midthickness_cortex.gii', meta_data{recording, 'sub'}));
        S.iskull = fullfile(char(extractBefore(meta_data{recording, "raw_data_loc"},'meg')),...
            'anat', 'inner_skull.surf.gii');
        D = spm_opm_create(S);
    end
    
    cd(analysedDataPath)

    %% Preprocess data

    % PSD
    S = [];
    S.D = D;
    S.channels = D.chanlabels(indchantype(D, 'MEGMAG', 'GOOD'));
    S.plot = 1;
    S.triallength = 10e3;
    [po, freq] = spm_opm_psd(S);

    % Badchannels
    if strcmp(meta_data{recording, "sub"}, 'sub-001')
        D = badchannels(D, selectchannels(D, 'regexp_(^(27|8|31|59|6|5|19|41|24|30|17)-.*)'), 1);
    elseif strcmp(meta_data{recording, "sub"}, 'sub-002')
        D = badchannels(D, selectchannels(D, 'regexp_(^(2|19|28|14|31|59|1|4)-.*)'), 1);
        D = badchannels(D, selectchannels(D, 'regexp_(^5-.*-Y$)'), 1);
    elseif strcmp(meta_data{recording, "sub"}, 'sub-003')
        D = badchannels(D, selectchannels(D, 'regexp_(^(28|48|41)-.*)'), 1);
    else
        warning('Badchannels not set for participant %s', meta_data{recording, "sub"});
    end
    save(D);
    
    % Sync with optitrack
    if isfile(['opti_data_', D.fname])
        load(['opti_data_', D.fname]);
        D = spm_eeg_load(['t_', D.fname]);
    else
        cfg = [];
        fname = D.fname;
        cfg.filename = [rawDataPath, '\', extractBefore(fname, '_meg.'), '_optitrack.csv'];
        opti_data = readRigidBody(cfg);
        [opti_data, D] = syncOptitrackAndOPMdata(opti_data,D,'TriggerChannelName',meta_data{recording, "OptiTrig"});
        save(['opti_data_', D.fname], 'opti_data');
    end

    % Plot opm recordings with position and rotation
    figure;
    t = tiledlayout(3,1);
    nexttile; 
    plot(D.time, 1e-6*D(indchantype(D, 'MEGMAG', 'GOOD'),:,1), 'LineWidth', 3);
    C = linspecer(length(indchantype(D, 'MEGMAG', 'GOOD')));
    set(gca, 'ColorOrder', C);
    ylabel('B (nT)');
    ylim([-15 15]);
    xticklabels({});
    xlim([70, max(D.time)])
    set(gca, 'FontSize', 16);
    grid on; box on;
        
    nexttile;
    if strcmp(opti_data.cfg.LengthUnits, 'Meters')
        plot(D.time, opti_data.Scannercast.RigidBody{:,7:9}-opti_data.Scannercast.RigidBody{1,7:9}, 'LineWidth', 3);
    else
        plot(D.time, 1e-3*(opti_data.Scannercast.RigidBody{:,7:9}-opti_data.Scannercast.RigidBody{1,7:9}), 'LineWidth', 3);
    end
    ylabel({'Displace-'; 'ment (m)'});
    xlim([70, max(D.time)])
    xticklabels({});
    set(gca, 'FontSize', 24);
    ylim([-2 2]);
    grid on; box on;
    legend({'Left-Right', 'Up-Down', 'Door-Screen'}, 'location', 'eastoutside');
   
    nexttile;
    plot(D.time, 180*quat2eul(opti_data.Scannercast.RigidBody{:,[6,3:5]}, 'XYZ')/pi, 'LineWidth', 3);
    xlabel('Time (s)');
    ylabel({'Rotation';'(deg)'});
    xlim([70, max(D.time)])
    set(gca, 'FontSize', 24);
    ylim([-360 360]);
    grid on; box on;
    legend({'Pitch', 'Yaw', 'Roll'}, 'location', 'eastoutside');

    t.TileSpacing = 'compact';
    set(gcf, 'Position', [680   344   1172   652]);
    print(fullfile(meta_data{recording, "results_save_loc"}, ...
        sprintf('%s_all_time_series', extractBefore(meta_data{recording, "raw_data_name"}, '_meg.lvm'))),'-dpng','-r300');
    
    % Filter
    if isfile(['fff', D.fname])
        D = spm_eeg_load(['fff', D.fname]);
    else
        S = [];
        S.band = 'low';
        S.freq = 40;
        S.D = D;
        D = spm_eeg_filter(S);
            
        S = [];
        S.band = 'high';
        S.freq = 2;
        S.D = D;
        D = spm_eeg_filter(S);
    
        S = [];
        S.band = 'stop';
        S.freq = [49 51];
        S.D = D;
        S.order = 5;
        D = spm_eeg_filter(S);
    end
    
    % Do both HFC and AMM (separately) to test difference later
    DD = cell(1,5);
    DD{1} = D;

    % HFC
    if isfile(['h', D.fname])
        DD{2} = spm_eeg_load(['h', D.fname]);
    else
        S = [];
        S.D = D;
        DD{2} = spm_opm_hfc(S);
    end

    if isfile(['h2', D.fname])
        DD{3} = spm_eeg_load(['h2', D.fname]);
    else
        S = [];
        S.D = D;
        S.L = 2;
        S.prefix = 'h2';
        DD{3} = spm_opm_hfc(S);
    end

    % AMM without temporal extension
    if isfile(['m2', D.fname])
        DD{4} = spm_eeg_load(['m2', D.fname]);
    else
        S = [];
        S.D = D;
        S.corrLim = 1;
        S.reducerank = 0;
        S.prefix = 'm2';
        DD{4} = spm_opm_amm(S);
    end

    % AMM
    if isfile(['m', D.fname])
        DD{5} = spm_eeg_load(['m', D.fname]);
    else
        S = [];
        S.D = D;
        S.corrLim = 0.95;
        S.reducerank = 0;
        DD{5} = spm_opm_amm(S);
    end

    % Plot PSD
    figure; hold on; grid on; box on;
    co = colororder(gca);
    line_style = {'-', '--', ':', '-.'};
    ii = 0;
    pl = [];
    for proc_step = [1,2,4,5]
        ii = ii+1;
        S = [];
        S.D = DD{proc_step};
        S.channels = D.chanlabels(indchantype(D, 'MEGMAG', 'GOOD'));
        S.triallength = 10e3;
        [po, freq] = spm_opm_psd(S);

        mp = median(po,2);
        sem = 1.2533*std(po,[],2)./sqrt(size(po,2));

        fill([freq'; flipud(freq')], [min(mp-sem,[],2); flipud(max(mp+sem,[],2))], co(ii,:),...
            'linestyle', 'none', 'FaceAlpha', 0.4)
        pl(end+1) = plot(freq, mp, 'LineWidth', 2, 'LineStyle', line_style{ii}, 'color', co(ii,:));

    end
    set(gca,'yscale','log');
    set(gca, 'FontSize', 22);
    xlim([2 40]);
    ylim([10 1e4]);
    legend(pl, {'No spatial filter', 'HFC', 'AMM spatial', 'AMM with temporal'});
    xlabel('Frequency (Hz)');
    ylabel('PSD ($$fT\sqrt[-1]{Hz}$$)','interpreter','latex');

    save_name = sprintf('PSD_after_temporal_filtering_%s', ...
        extractBefore(meta_data{recording, "raw_data_name"}, '_meg.lvm'));
    print(fullfile(meta_data{recording, "results_save_loc"}, save_name),'-dpng','-r300');
    
    % Plot shielding factors
    figure; hold on; grid on; box on;
    ii = 1;
    pl = [];
    for proc_step = [2,4,5]
        ii = ii + 1;

        S = [];
        S.D1 = DD{1};
        S.D2 = DD{proc_step};
        S.channels = D.chanlabels(indchantype(D, 'MEGMAG', 'GOOD'));
        S.plot = 0;
        S.triallength = 10e3;
        S.dB = 1;
        [shield, freq] = spm_opm_rpsd(S);

        mp = median(shield,2);
        sem = 1.2533*std(shield,[],2)./sqrt(size(shield,2));

        fill([freq'; flipud(freq')], [min(mp-sem,[],2); flipud(max(mp+sem,[],2))], co(ii,:),...
            'linestyle', 'none', 'FaceAlpha', 0.4)
        pl(end+1) = plot(freq, mp, 'LineWidth', 2, 'color', co(ii,:), 'LineStyle', line_style{ii});
    end

    set(gca, 'FontSize', 22);
    xlim([2 40]);
    ylim([0 30]);
    legend(pl, {'HFC', 'AMM spatial', 'AMM with temporal'});
    xlabel('Frequency (Hz)');
    ylabel('Shielding Factor (dB)');

    save_name = sprintf('Shielding_factor_after_temporal_filtering_%s', extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'));
    print(fullfile(meta_data{recording, "results_save_loc"}, save_name),'-dpng','-r300');
    
    
    % Epoch
    for pp = 1:length(DD)

        D = DD{pp};

        if isfile(['e_', D.fname])
            D = spm_eeg_load(['e_', D.fname]);
        else
            S = [];
            S.D = D;
            S.timewin = [-200 500];
            S.condLabels = {'tone'};
            S.triggerChannels = {char(meta_data{recording, "AudioTrig"})};
            D = spm_opm_epoch_trigger(S);
        end
    
        % Set all epochs after 500 to bad
        goodTrials = indtrial(D, 'tone', 'GOOD');
        D = badtrials(D, goodTrials(501:end), 1);
        save(D);

        DD{pp} = D;
    end

    close all;
end

%% Find spatiotemporal clusters at sensor-level

clearvars -except meta_data colormap123 delay
rng(76);

for recording = 1:size(meta_data,1)
    fprintf('Recording no: %.f of %.f\n', recording, size(meta_data,1));
    cd(meta_data{recording, "analysed_data_loc"});
    start_string = {'e_ffft_', 'e_hffft_', 'e_h2ffft_', 'e_m2ffft_', 'e_mffft_'};
    for pp = 1:length(start_string)
        DD = spm_eeg_load(char(fullfile(meta_data{recording, "analysed_data_loc"}, ...
            strcat(start_string{pp}, extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'), '.mat'))));

        % Convert to fieldtrip
        data = ftraw(DD);

        % Just select good MEG channels
        cfg = [];
        cfg.channel = DD.chanlabels(indchantype(DD, 'MEGMAG', 'GOOD'));
        cfg.trials = indtrial(DD, 'tone', 'GOOD');
        data = ft_selectdata(cfg, data);

        % Find neighbours in topoplot layout
        cfg = [];
        cfg.method = 'distance';
        cfg.neighbourdist = 50;
        neighbours = ft_prepare_neighbours(cfg, data);

        % Create a zeros dataset for comparison
        data_zeros = data;
        data_zeros.trial = repmat({zeros(size(data.trial{1}))}, 1, size(data.trial,2));

        % Timelock
        cfg = [];
        cfg.keeptrials = 'yes';
        tl_data = ft_timelockanalysis(cfg, data);
        tl_zeros = ft_timelockanalysis(cfg, data_zeros);

        % Permutation test
        cfg = [];
        cfg.method = 'montecarlo';
        cfg.statistic = 'indepsamplesT';
        cfg.correctm = 'cluster';
        cfg.clusteralpha = 0.05;
        cfg.clusterstatistic = 'maxsum';
        cfg.minnbchan = 0; 
        cfg.neighbours = neighbours;
        cfg.tail = 0;
        cfg.alpha = 0.025;
        cfg.numrandomization = 500;
        cfg.latency = [50 150]*1e-3 + delay*1e-3; % Set time window to test over

        n_zeros  = size(tl_zeros.trial, 1);
        n_toi = size(tl_data.trial, 1);

        cfg.design = [ones(1,n_zeros), ones(1,n_toi)*2];
        cfg.ivar = 1;
        cfg.channel = DD.chanlabels(indchantype(DD, 'MEGMAG', 'GOOD'));
        [stat] = ft_timelockstatistics(cfg, tl_data, tl_zeros);

        % Average data
        cfg    = [];
        avg_dat  = ft_timelockanalysis(cfg, data);
        
        % Vector of all p-values associated with the clusters from ft_timelockstatistics.
        pos_cluster_pvals = [stat.posclusters(:).prob];

        % Which clusters are interesting to visualize
        pos_clust = find(pos_cluster_pvals < 0.025);
        pos = ismember(stat.posclusterslabelmat, pos_clust);

        % Negative clusters
        neg_cluster_pvals = [stat.negclusters(:).prob];
        neg_clust = find(neg_cluster_pvals < 0.025);
        neg = ismember(stat.negclusterslabelmat, neg_clust);

        % Plot topographies at every 10 ms
        timestep = 0.01; % timestep between time windows for each subplot (in seconds)
        sampling_rate = data.fsample; % Data has a temporal resolution of 300 Hz
        
        % Just select data between 50 and 150 ms
        tinds = find((stat.time*1e3 >= 50 + delay).*(stat.time*1e3 <= 150 + delay));
        sample_count  = length(tinds);
        j = 50e-3+delay*1e-3:timestep:150e-3+delay*1e-3; % Temporal endpoints (in seconds) of the ERP average computed in each subplot
        m = tinds(1):timestep*sampling_rate:tinds(end); % temporal endpoints in M/EEG samples
       
        % Get layout
        lay_name = fullfile(meta_data{recording, "analysed_data_loc"}, ...
            sprintf('%s_2Dlayout.mat', extractBefore(meta_data{recording, "raw_data_name"}, '_meg.lvm')));

        if isfile(lay_name)
            load(lay_name);
        else
            fid = fiducials(D);
            fid_struct = struct('NAS', fid.fid.pnt(contains(fid.fid.label, 'nas'),:), ...
                'LPA', fid.fid.pnt(contains(fid.fid.label, 'lpa'),:), ...
                'RPA', fid.fid.pnt(contains(fid.fid.label, 'rpa'),:));
            lay = spm_get_anatomical_layout(D.sensors('MEG').coilpos(endsWith(D.sensors('MEG').label, ['-', rad_ax]),:), ...
                D.sensors('MEG').label(endsWith(D.sensors('MEG').label, ['-', rad_ax])),...
                double(gifti(D.inv{1}.mesh.tess_scalp).vertices), fid_struct, 0);
            save(lay_name, 'lay');
        end
        [i1,i2] = match_str(avg_dat.label, stat.label);

        figure;
        for k = 1:length(m)-1
           subplot(2, 5, k);
           cfg = [];
           cfg.xlim = [j(k) j(k+1)];
           cfg.zlim = [-1 1]*465;

           pos_int = zeros(numel(avg_dat.label),1);
           neg_int = zeros(numel(avg_dat.label),1);
           pos_int(i1) = all(pos(i2, m(k):m(k+1)), 2);
           neg_int(i1) = all(neg(i2, m(k):m(k+1)), 2);

           cfg.highlight   = 'on';
           % Get the index of the to-be-highlighted channel
           cfg.highlightchannel = find(pos_int | neg_int);
           cfg.comment = 'no';
           cfg.layout = lay;
           cfg.interactive = 'no';
           cfg.figure = 'gca'; % plots in the current axes, here in a subplot
           cfg.colormap = colormap123;
           cfg.highlightsymbol = '*';
           cfg.markersymbol = '.';
           ft_topoplotER(cfg, avg_dat);
           title(gca, sprintf('%.2f - %.2f ms', j(k)*1e3, j(k+1)*1e3), 'FontSize', 16);
        end

        set(gcf, 'Position', [-1918, -51, 1424, 1001]);

        if pp == 1
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'no_amm');
        elseif pp == 2
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'hfc');
        elseif pp == 3
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'hfc_with_gradients');
        elseif pp == 4
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'amm_spatial');
        else
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'amm');
        end

        save_name = sprintf('%s_anti_averaging_topo_over_time', extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'));
        print(fullfile(save_loc, save_name),'-dpng','-r300');

        % Plot single topography
        figure;
       subplot(1, 1, 1);
       cfg = [];
       cfg.xlim = [80+delay 120+delay]*1e-3;
       cfg.zlim = [-1 1]*250;

       pos_int = zeros(numel(avg_dat.label),1);
       neg_int = zeros(numel(avg_dat.label),1);
       pos_int(i1) = any(pos(i2, logical((stat.time*1e3 >= 80 + delay).*(stat.time*1e3 <= 120 + delay))), 2);
       neg_int(i1) = any(neg(i2, logical((stat.time*1e3 >= 80 + delay).*(stat.time*1e3 <= 120 + delay))), 2);

       cfg.highlight   = 'on';
       % Get the index of the to-be-highlighted channel
       cfg.highlightchannel = find(pos_int | neg_int);
       cfg.comment = 'no';
       cfg.layout = lay;
       cfg.interactive = 'no';
       cfg.figure = 'gca'; 
       cfg.colormap = colormap123;
       cfg.highlightsymbol = '*';
       cfg.markersymbol = '.';
       ft_topoplotER(cfg, avg_dat);
       colorbar;
       set(gcf, 'Position', [994   704   404   274]);
       

        save_name = sprintf('%s_anti_averaging_topo_80_to_120_ms_anysig', extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'));
        print(fullfile(save_loc, save_name),'-dpng','-r300');


        close all;

        
        % --- Plot time series

        % T-stat
        good_trials = indtrial(DD, 'tone', 'GOOD');
        se = std(DD(indchantype(DD, 'MEGMAG', 'GOOD'),:,good_trials),[],3)./sqrt(length(good_trials));
        t = mean(DD(indchantype(DD, 'MEGMAG', 'GOOD'),:,good_trials),3)./se;
        
        [~, tind] = min(abs(DD.time - (100*1e-3 + delay+1e-3)));
        max_colour = interp1([min(abs(t(:,tind))), max(abs(t(:,tind)))], [0.9, 0], abs(t(:,tind)));
        [~, plot_order] = sort(max_colour, 'descend');
        max_colour = repmat(max_colour, 1, 3);
    
        % Create figure
        figure; hold on; grid on; box on;
        if contains(meta_data{recording, "raw_data_name"}, 'seat') || contains(meta_data{recording, "raw_data_name"}, 'Seat')
            ylim([-25 25]);
        else
            ylim([-1 1]*13);
        end
        yl = ylim;

        % highlight time window tested over
        fill([min(stat.time), max(stat.time), max(stat.time), min(stat.time)]*1e3 - delay, ...
            [yl(1), yl(1), yl(2), yl(2)], [245, 238, 158]./255, 'EdgeColor', 'None', 'FaceAlpha', 0.7)
        
        % Plot t-stat
        for chan = 1:size(t,1)
            plot(1e3*DD.time - delay, t(plot_order(chan),:), 'color', max_colour(plot_order(chan),:), 'LineWidth', 2);
        end

        % highlight time points where there is at least one significant
        % cluster
        sigtimes = any(cat(1, pos, neg), 1);
        sigtimes = stat.time(sigtimes)*1e3 - delay;
        if any(sigtimes)
            plot(sigtimes, 0.9*yl(2)*ones(size(sigtimes)), '*', 'color', [59, 142, 165]./255, 'MarkerSize', 5)
            plot(sigtimes, 0.9*yl(1)*ones(size(sigtimes)), '*', 'color', [59, 142, 165]./255, 'MarkerSize', 5)
        end
        
        xlim([-100 400]);
        
        xlabel('Time (ms)');
        ylabel('t-stat');
        set(gcf, 'Position', [680   654   451   344]);
        set(gca, 'FontSize', 24);
        fname = DD.fname;

        save_name = sprintf('%s_t_stat_time_series', extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'));
        print(fullfile(save_loc, save_name),'-dpng','-r300');

        % Plot average
        % Create figure
        figure; hold on; grid on; box on;
        ylim([-550 550]);
        yl = ylim;

        % highlight time window tested over
        fill([min(stat.time), max(stat.time), max(stat.time), min(stat.time)]*1e3 - delay, ...
            [yl(1), yl(1), yl(2), yl(2)], [245, 238, 158]./255, 'EdgeColor', 'None', 'FaceAlpha', 0.7)

        % Plot data
        dat = mean(DD(indchantype(DD, 'MEGMAG', 'GOOD'),:,good_trials),3);
        for chan = 1:size(t,1)
            plot(1e3*DD.time - delay, dat(plot_order(chan),:), 'color', max_colour(plot_order(chan),:), 'LineWidth', 2);
        end

        % highlight time points where there is at least one significant
        % cluster
        if any(sigtimes)
            plot(sigtimes, 0.9*yl(2)*ones(size(sigtimes)), '*', 'color', [59, 142, 165]./255, 'MarkerSize', 5)
            plot(sigtimes, 0.9*yl(1)*ones(size(sigtimes)), '*', 'color', [59, 142, 165]./255, 'MarkerSize', 5)
        end

        xlim([-100 400]);
        xlabel('Time (ms)');
        ylabel('B (fT)');
        grid on;
        set(gcf, 'Position', [680   654   451   344]);
        set(gca, 'FontSize', 24);
        save_name = sprintf('%s_average_time_series', extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'));
        print(fullfile(save_loc, save_name),'-dpng','-r300');
    
    end
    close all;

end

% Plot for colorbar
figure;
subplot(1, 1, 1);
cfg = [];
cfg.xlim = [80+delay 120+delay]*1e-3;
cfg.zlim = [-1 1]*250;
cfg.highlight   = 'on';
cfg.highlightchannel = find(pos_int | neg_int);
cfg.comment = 'no';
cfg.layout = lay;
cfg.interactive = 'no';
cfg.figure = 'gca'; 
cfg.colormap = colormap123;
cfg.highlightsymbol = '*';
cfg.markersymbol = '.';
ft_topoplotER(cfg, avg_dat);
cb = colorbar('southoutside');
set(cb, 'FontSize', 22)
ylabel(cb, 'B (fT)', 'FontSize', 24)
set(gcf, 'Position', [994   525   404   453]);

print(fullfile(extractBefore(meta_data{end, "results_save_loc"}, 'sub-'), 'anti_averaging_topo_colorbar'),'-dpng','-r300');

%% Dipole fit

clearvars -except meta_data colormap123 delay

rng(76);

for recording = 1:size(meta_data,1)
    cd(meta_data{recording, "analysed_data_loc"});
    start_string = {'e_ffft_', 'e_hffft_', 'e_h2ffft_', 'e_m2ffft_', 'e_mffft_'}; 
    for pp = 1:length(start_string)
        DD{pp} = spm_eeg_load(char(fullfile(meta_data{recording, "analysed_data_loc"}, ...
            strcat(start_string{pp}, extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'), '.mat'))));
    end

    % Get deviant and standard labels
    stim = readtable(meta_data{recording, "stim_data_fname"});
    deviants = find(contains(stim.Condition, 'deviant'));
    trial_length = diff(deviants);
    trial_length = cat(1, trial_length, size(DD{pp},3) - max(deviants) + 1);

    % Take last tone of each set as standard
    standards = deviants(2:end)-1;
    standards = cat(1, standards, size(DD{pp},3));
    
    % Dipole fit

    % Initialise at auditory cortices
    aud_mni = [-54 -14 8 1; 54 -14 8 1]'; % Auditory cortices in MNI space
    aud_nat = DD{1}.inv{1}.datareg.fromMNI*aud_mni;
    aud_nat = aud_nat(1:3,:)';
    aud_nat = aud_nat*1e-3; % convert to m

    % Prepare headmodel
    mesh = ft_read_headshape(DD{1}.inv{1}.mesh.tess_iskull);
    cfg = [];
    cfg.method = 'singleshell';
    cfg.siunits = 'yes';
    headmodel = ft_prepare_headmodel(cfg, mesh);

    % Cortex
    ctx = ft_read_headshape(DD{1}.inv{1}.mesh.tess_ctx);
    ctx = ft_convert_units(ctx, 'm');

    % Prepare sourcemodel - use cortical mesh
    cfg = [];
    cfg.method = 'basedoncortex';
    cfg.headshape = ctx;
    cfg.headmodel = headmodel;
    cfg.inwardshift = 0;
    src = ft_prepare_sourcemodel(cfg);
    

    % Read MRI for plotting
    mri_orig = ft_read_mri(DD{1}.inv{1}.mesh.sMRI);

    for pp = 1:length(DD)

        % Prepare leadfields
        sens = DD{pp}.sensors('MEG');
        sens = ft_convert_units(sens, 'm');
        cfg                  = [];
        cfg.grad             = sens;
        cfg.headmodel        = headmodel;
        cfg.reducerank       = 2;
        cfg.channel          = DD{4}.chanlabels(indchantype(DD{pp}, 'MEGMAG', 'GOOD'));
        cfg.sourcemodel      = src;
        sourcemodel = ft_prepare_leadfield(cfg);

        % Format data for fieldtrip
        data = ftraw(DD{pp});
        cfg = [];
        cfg.trials = indtrial(DD{pp}, 'tone', 'GOOD');
        data = ft_selectdata(cfg, data);

        % Average
        tl_data = ft_timelockanalysis([], data);

        % Dipole fit
        cfg = [];
        cfg.latency = [0.08 0.12]+delay*1e-3;
        cfg.numdipoles = 2;
        cfg.symmetry = [];
        cfg.gridsearch = 'no';
        cfg.dip.pos = aud_nat;
        cfg.headmodel = headmodel;
        cfg.sourcemodel = sourcemodel;
        cfg.channel = DD{pp}.chanlabels(indchantype(DD{pp}, 'MEGMAG', 'GOOD'));
        cfg.senstype = 'meg';
        source = ft_dipolefitting(cfg, tl_data);
        source.dip = ft_convert_units(source.dip, 'mm');

        % Plot dipole position
        if pp == 1
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'no_amm');
        elseif pp == 2
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'hfc');
        elseif pp == 3
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'hfc_with_gradients');
        elseif pp == 4
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'amm_spatial');
        else
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'amm');
        end

        % Axial
        pos = mean(source.dip.pos,1);
        figure; hold on;
        ft_plot_dipole([source.dip.pos(1,[1,2]), pos(3)+150], mean(source.dip.mom(1:3,:),2), 'color', '#E07BE0', 'unit', 'mm'); % Left
        ft_plot_dipole([source.dip.pos(2,[1,2]), pos(3)+150], mean(source.dip.mom(4:6,:),2), 'color', '#45C9B7', 'unit', 'mm'); % Right
        ft_plot_slice(mri_orig.anatomy, 'transform', mri_orig.transform, 'location', pos, 'orientation', [0 0 1], 'resolution', 0.1);
        view(0,90);
        axis tight
        axis off
        set(gcf, 'Position', [680    50   560   946]);
        save_name = sprintf('%s_ft_dip_fit_axial', extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'));
        print(fullfile(save_loc, save_name),'-dpng','-r300');

        % Coronal
        figure; hold on;
        ft_plot_dipole([source.dip.pos(1,1), pos(2)-100, source.dip.pos(1,3)], mean(source.dip.mom(1:3,:),2), 'color', '#E07BE0', 'unit', 'mm'); % Left
        ft_plot_dipole([source.dip.pos(2,1), pos(2)-150, source.dip.pos(2,3)], mean(source.dip.mom(4:6,:),2), 'color', '#45C9B7', 'unit', 'mm'); % Right
        ft_plot_slice(mri_orig.anatomy, 'transform', mri_orig.transform, 'location', pos, 'orientation', [0 1 0], 'resolution', 0.1);
        view(0,0);
        axis tight
        axis off
        set(gcf, 'Position', [680    50   560   946]);
        save_name = sprintf('%s_ft_dip_fit_coronal', extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'));
        print(fullfile(save_loc, save_name),'-dpng','-r300');
        

        % Plot estimated source current
        % Get lead field along mean dipole ori
        sens = data.grad;
        sens = ft_convert_units(sens, 'm');
        [headmodel, sens] = ft_prepare_vol_sens(headmodel, sens, 'channel', DD{pp}.chanlabels(indchantype(DD{pp}, 'MEGMAG', 'GOOD')));
        Gxyz = ft_compute_leadfield(source.dip.pos*1e-3, sens, headmodel, 'dipoleunit', 'nA*m', 'chanunit', repmat({'fT'}, size(sens.label,1),1));
        L = zeros(size(Gxyz,1),2);
        for ind = 1:2
            dip_ori = mean(source.dip.mom((ind-1)*3+1:ind*3,:),2);
            dip_ori = dip_ori./norm(dip_ori);
            L(:, ind) = Gxyz(:, (3*ind-2):(3*ind))*dip_ori;
        end

        good_trials = indtrial(DD{pp}, 'tone', 'GOOD');
        X_evoked = zeros(size(L,2), size(DD{pp},2), length(good_trials));
        X_standards = zeros(size(L,2), size(DD{pp},2), length(standards));
        X_deviants = zeros(size(L,2), size(DD{pp},2), length(deviants));

        for tt = 1:length(good_trials)
            X_evoked(:,:,tt) = pinv(L)*DD{pp}(indchannel(DD{pp},sens.label),:,good_trials(tt));
        end
        for tt = 1:length(deviants)
            X_deviants(:,:,tt) = pinv(L)*DD{pp}(indchannel(DD{pp},sens.label),:,deviants(tt));
        end
        for tt = 1:length(standards)
            X_standards(:,:,tt) = pinv(L)*DD{pp}(indchannel(DD{pp},sens.label),:,standards(tt));
        end

        % T-test across trials
        SE_evoked = std(X_evoked, [], 3)./sqrt(size(X_evoked,3));
        t_evoked = mean(X_evoked, 3)./SE_evoked;
        SE_standards = std(X_standards, [], 3)./sqrt(size(X_standards, 3));
        t_standards = mean(X_standards,3)./SE_standards;
        SE_deviants = std(X_deviants, [], 3)./sqrt(size(X_deviants, 3));
        t_deviants = mean(X_deviants,3)./SE_deviants;
        
        % Unpaired t-test equal variance between deviants and standards for MMN response
        n1 = size(X_deviants,3);
        n2 = size(X_standards,3);
        SE = sqrt(((n1-1)*std(X_deviants, [], 3).^2 + (n2-1)*std(X_standards, [], 3).^2)./(n1 + n2 - 2))*...
            sqrt(1/n1 + 1/n2);
        t_diff = (mean(X_deviants,3) - mean(X_standards,3))./SE;


        % Plot
        % Evoked response:
        figure;
        hold on; grid on; box on;
        ylim([-1 1]*7);
        yl = ylim;
        fill([min(cfg.latency), max(cfg.latency), max(cfg.latency), min(cfg.latency)]*1e3 - delay, ...
            [yl(1), yl(1), yl(2), yl(2)], [231, 196, 170]./255, 'EdgeColor', 'None', 'FaceAlpha', 0.3);
        plot(DD{pp}.time*1e3 - delay, t_standards, 'LineWidth', 3, 'LineStyle', '-');
        colororder(gca, ["#E07BE0", "#45C9B7"]);
        set(gca, 'FontSize', 18);
        xlim([-100 400]);
        xlabel('Time (ms)', 'FontSize', 18);
        % ylabel({'Estimated Source', 'Current (nAm)'}, 'FontSize', 18);
        ylabel('t-stat', 'FontSize', 18);
        set(gcf, 'Position', [680   654   451   344]);
        save_name = sprintf('%s_dipfit_evoked_dipole', extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'));
        print(fullfile(save_loc, save_name),'-dpng','-r300');

        % MMN:
        figure; 
        t = tiledlayout(1,2);
        yl = [];
        for ind = 1:2
            nexttile(t); hold on; grid on; box on;

            % Standards
            % plot(DD{pp}.time*1e3 - delay, mean(X_standards(ind,:,:), 3), 'LineWidth', 3, 'LineStyle', '--', 'color', [0.3639    0.5755    0.7484]);
            plot(DD{pp}.time*1e3 - delay, t_standards(ind,:), 'LineWidth', 3, 'LineStyle', '--', 'color', [0.3639    0.5755    0.7484]);

            % Deviants
            % plot(DD{pp}.time*1e3 - delay, mean(X_deviants(ind,:,:), 3), 'LineWidth', 3, 'LineStyle', ':', 'color', [0.9153    0.2816    0.2878]);
            plot(DD{pp}.time*1e3 - delay, t_deviants(ind,:), 'LineWidth', 3, 'LineStyle', ':', 'color', [0.9153    0.2816    0.2878]);

            % Difference
            % plot(DD{pp}.time*1e3 - delay, mean(X_deviants(ind,:,:),3) - mean(X_standards(ind,:,:), 3), 'LineWidth', 3, 'LineStyle', '-', 'color', [0.3373    0.3020    0.2902]);
            plot(DD{pp}.time*1e3 - delay, t_diff(ind,:), 'LineWidth', 3, 'LineStyle', '-', 'color', [0.3373    0.3020    0.2902]);
            set(gca, 'FontSize', 18);
            xlim([-100 400]);
            % yl(ind) = max(abs(ylim));
            xlabel('Time (ms)');

            if ind == 1
                title('Left Hemisphere');
            else
                title('Right Hemisphere');
            end

        end

        % Set axes limits and legend
        % ylim(t.Children, [-1 1]*max(abs(yl)));
        ylim(t.Children, [-1 1]*10);
        % ylabel(t, {'Estimated Source', 'Current (nAm)'}, 'FontSize', 18);
        ylabel(t, 't-stat', 'FontSize', 18);
        lgd = legend('Standards', 'Deviants', 'MMN', 'location', 'eastoutside');
        set(gcf, 'Position', [626   476   821   285]);

        % Add text to indicate how many trials per condition
        
        annotation('textbox', [lgd.Position(1), lgd.Position(2) - 0.35, lgd.Position(3), 0.3], ...
            'string', sprintf('# deviants: %.f\n# standards: %.f', length(deviants), length(standards)), 'FontSize', 16, 'EdgeColor', 'None');

        save_name = sprintf('%s_dipfit_MMN_trace_dipole_all_sets', extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'));
        print(fullfile(save_loc, save_name),'-dpng','-r300');
        
    end

    close all
end

%% Minimum Norm Estimation

clearvars -except meta_data colormap123 delay

l_max_vals = zeros(size(meta_data,1), 5);
r_max_vals = zeros(size(meta_data,1), 5);

for recording = 1:size(meta_data,1)

    fprintf('Recording no: %.f of %.f\n', recording, size(meta_data,1));

    % Get data
    cd(meta_data{recording, "analysed_data_loc"});
    start_string = {'e_ffft_', 'e_hffft_', 'e_h2ffft_', 'e_m2ffft_', 'e_mffft_'}; %
    for pp = 1:length(start_string)
        DD{pp} = spm_eeg_load(char(fullfile(meta_data{recording, "analysed_data_loc"}, ...
            strcat(start_string{pp}, extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'), '.mat'))));

        % Baseline correct
        if isfile(fullfile(DD{pp}.path, ['b', DD{pp}.fname]))
            DD{pp} = spm_eeg_load(fullfile(DD{pp}.path, ['b', DD{pp}.fname]));
        else
            S = [];
            S.D = DD{pp};
            DD{pp} = spm_eeg_bc(S);
        end
    end

    % Create infalted cortex mesh
    inflated_ctx_fname = fullfile(char(extractBefore(meta_data{recording, "raw_data_loc"},'meg')),...
                'anat', sprintf('%s.L.inflated.4k_fs_LR.surf.gii', meta_data{recording, 'sub'}));
    left_inflated_ctx = ft_read_headshape(inflated_ctx_fname);
    right_inflated_ctx = ft_read_headshape(strrep(inflated_ctx_fname, '.L.', '.R.'));

    if ~isfile(fullfile(char(extractBefore(meta_data{recording, "raw_data_loc"},'meg')),...
        'anat', sprintf('%s_inflated_cortex.gii', meta_data{recording, 'sub'})))
        inflated_ctx_fname = fullfile(char(extractBefore(meta_data{recording, "raw_data_loc"},'meg')),...
                'anat', sprintf('%s.L.inflated.4k_fs_LR.surf.gii', meta_data{recording, 'sub'}));
        combine_surfaces({inflated_ctx_fname, strrep(inflated_ctx_fname, '.L.', '.R.')}, ...
            fullfile(char(extractBefore(meta_data{recording, "raw_data_loc"},'meg')),...
            'anat', sprintf('%s_inflated_cortex.gii', meta_data{recording, 'sub'})));
    end

    inflated_ctx = ft_read_headshape(fullfile(char(extractBefore(meta_data{recording, "raw_data_loc"},'meg')),...
                'anat', sprintf('%s_inflated_cortex.gii', meta_data{recording, 'sub'})));


    % Source reconstruct (MNE)
    for pp = [length(DD), 1:length(DD)-1]
        
        if ~isfile([DD{pp}.path, '\', extractBefore(DD{pp}.fname, '.mat'), '_1_t60_160_f2_40_1.gii'])
            matlabbatch = [];
            matlabbatch{1}.spm.meeg.source.invert.D = {fullfile(DD{pp})};
            matlabbatch{1}.spm.meeg.source.invert.val = 1;
            matlabbatch{1}.spm.meeg.source.invert.whatconditions.all = 1;
            matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.invtype = 'IID';
            matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.woi = [-Inf Inf];
            matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.foi = [2 40];
            matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.hanning = 1;
            matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.priors.priorsmask = {''};
            matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.priors.space = 1;
            matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.restrict.locs = zeros(0, 3);
            matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.restrict.radius = 32;
            matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.restrict.mask = {''};
            matlabbatch{1}.spm.meeg.source.invert.modality = {'All'};
            matlabbatch{2}.spm.meeg.source.results.D(1) = cfg_dep('Source inversion: M/EEG dataset(s) after imaging source reconstruction', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','D'));
            matlabbatch{2}.spm.meeg.source.results.val = 1;
            matlabbatch{2}.spm.meeg.source.results.woi = [50 150]+delay;
            matlabbatch{2}.spm.meeg.source.results.foi = [2 40];
            matlabbatch{2}.spm.meeg.source.results.ctype = 'evoked';
            matlabbatch{2}.spm.meeg.source.results.space = 0;
            matlabbatch{2}.spm.meeg.source.results.format = 'mesh';
            matlabbatch{2}.spm.meeg.source.results.smoothing = 8;
    
            a = spm_jobman('run',matlabbatch);
        end
        source_data = export(gifti([DD{pp}.path, '\', extractBefore(DD{pp}.fname, '.mat'), '_1_t60_160_f2_40_1.gii']), 'patch').facevertexcdata;

        % Create fieldtrip data structure

        % Left hemisphere
        l_source_data_toi = [];
        l_source_data_toi.pos = left_inflated_ctx.pos;
        l_source_data_toi.inside = ones(size(l_source_data_toi.pos,1), 1)==1;
        l_source_data_toi.pow = source_data(1:size(l_source_data_toi.pos,1));
        l_source_data_toi.tri = left_inflated_ctx.tri;
        
        maxval = max(l_source_data_toi.pow);
        l_max_vals(recording,pp) = maxval;
        % if strcmp(meta_data{recording, 'sub'}, 'sub-001')
        %     maxval = 2.9462;
        % elseif strcmp(meta_data{recording, 'sub'}, 'sub-002')
        %     maxval = 11.0019;
        % elseif strcmp(meta_data{recording, 'sub'}, 'sub-003')
        %     maxval = 7.5587;
        % end
        maxval = l_max_vals(recording,end);
        l_source_data_toi.mask = l_source_data_toi.pow >= 0.5*maxval;
        
        % Plot
        cfg                     = [];
        cfg.method              = 'surface';
        cfg.facecolor           = [0.4 0.4 0.4];
        cfg.vertexcolor         = 'none';
        cfg.funparameter        = 'pow';
        cfg.location            = 'max';
        cfg.maskparameter       = 'mask';
        if any(l_source_data_toi.mask)
            cfg.funcolorlim         = [0 maxval];
            ft_sourceplot(cfg, l_source_data_toi);
            colormap('hot') % change the colormap
        else
            figure;
            surf.pos = left_inflated_ctx.pos;
            surf.tri = left_inflated_ctx.tri;
            ft_plot_mesh(surf,'edgecolor', 'none', 'facecolor', cfg.facecolor, 'vertexcolor', cfg.vertexcolor);
            lighting gouraud
            camlight
        end

        % Save
        view ([-90 5])             % rotate the object in the view
        cl1 = camlight('headlight');
        set(gcf, 'color', 'w');
        set(gca, 'FontSize', 26);
        material dull
    
        if pp == 1
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'no_amm');
        elseif pp == 2
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'hfc');
        elseif pp == 3
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'hfc_with_gradients');
        elseif pp == 4
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'amm_spatial');
        else
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'amm');
        end
    
        save_name = sprintf('%s_min_norm_pow_left', extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'));
        print(fullfile(save_loc, save_name),'-dpng','-r300');

        % Right hemisphere
        r_source_data_toi = [];
        r_source_data_toi.pos = right_inflated_ctx.pos;
        r_source_data_toi.inside = ones(size(r_source_data_toi.pos,1), 1)==1;
        r_source_data_toi.pow = source_data(size(l_source_data_toi.pos,1)+1:end);
        r_source_data_toi.tri = right_inflated_ctx.tri;
        maxval = max(r_source_data_toi.pow);
        r_max_vals(recording,pp) = maxval;
        % if strcmp(meta_data{recording, 'sub'}, 'sub-001')
        %     maxval = 11.7462;
        % elseif strcmp(meta_data{recording, 'sub'}, 'sub-002')
        %     maxval = 10.9519;
        % elseif strcmp(meta_data{recording, 'sub'}, 'sub-003')
        %     maxval = 11.2778;
        % end
        maxval = r_max_vals(recording,end);
        r_source_data_toi.mask = r_source_data_toi.pow >= 0.5*maxval;
        
        % Plot
        cfg                     = [];
        cfg.method              = 'surface';
        cfg.facecolor           = [0.4 0.4 0.4];
        cfg.vertexcolor         = 'none';
        cfg.funparameter        = 'pow';
        cfg.location            = 'max';
        cfg.maskparameter       = 'mask';
        if any(r_source_data_toi.mask)
            cfg.funcolorlim         = [0 maxval];
            ft_sourceplot(cfg, r_source_data_toi);
            colormap('hot') % change the colormap
        else
            figure;
            surf.pos = right_inflated_ctx.pos;
            surf.tri = right_inflated_ctx.tri;
            ft_plot_mesh(surf,'edgecolor', 'none', 'facecolor', cfg.facecolor, 'vertexcolor', cfg.vertexcolor);
            lighting gouraud
            camlight
        end
    
        view ([90 0])
        cl1 = camlight('headlight');
        set(gcf, 'color', 'w');
        set(gca, 'FontSize', 26);
        material dull
    
        save_name = sprintf('%s_min_norm_pow_right', extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'));
        print(fullfile(save_loc, save_name),'-dpng','-r300');

    end
    close all
end

%% Plot evoked response at source level, using ROI from MNI

clearvars -except meta_data colormap123 delay

for recording = 1:size(meta_data,1)

    start_string = {'e_ffft_', 'e_hffft_', 'e_h2ffft_', 'e_m2ffft_', 'e_mffft_'};
    for pp = 1:length(start_string)
        DD{pp} = spm_eeg_load(char(fullfile(meta_data{recording, "analysed_data_loc"}, ...
            strcat(start_string{pp}, extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'), '.mat'))));
    end

    % Prepare headmodel
    mesh = ft_read_headshape(DD{1}.inv{1}.mesh.tess_iskull);
    cfg = [];
    cfg.method = 'singleshell';
    cfg.siunits = 'yes';
    headmodel = ft_prepare_headmodel(cfg, mesh);

    % Get locations of interest
    aud_mni = [-54 -14 8 1; 54 -14 8 1]'; % Auditory cortices in MNI space
    aud_nat = DD{1}.inv{1}.datareg.fromMNI*aud_mni;
    aud_nat = aud_nat(1:3,:)';
    aud_nat = aud_nat*1e-3; % convert to m

    % Get leadfields at positions of interest
    sens = DD{1}.sensors('MEG');
    sens = ft_convert_units(sens, 'm');
    [headmodel, sens] = ft_prepare_vol_sens(headmodel, sens, 'channel', DD{1}.chanlabels(indchantype(DD{1}, 'MEGMAG', 'GOOD')));
    Gxyz = ft_compute_leadfield(aud_nat, sens, headmodel, 'dipoleunit', 'nA*m', 'chanunit', repmat({'fT'}, size(sens.label,1),1));


    for pp = 1:length(DD)

        % Get data
        data = ftraw(DD{pp});
        cfg = [];
        cfg.channel = indchantype(DD{pp}, 'MEGMAG', 'GOOD');
        cfg.trials = indtrial(DD{pp}, 'tone', 'GOOD');
        avdata = ft_timelockanalysis(cfg, data);

        % Optimise orientation based on data
        dip_data = pinv(Gxyz)*avdata.avg;
        single_ori_dip_data = zeros(2, size(dip_data,2));
        for ind = 1:2
            [U,S,V] = svd(dip_data(3*ind-2:3*ind,:)*dip_data(3*ind-2:3*ind,:)');
            single_ori_dip_data(ind,:) = dip_data(3*ind-2:3*ind,:)'*U(:,1);
        end

        figure; hold on; grid on; box on;
        plot(DD{pp}.time*1e3, single_ori_dip_data, 'LineWidth', 2);

        xlim([-100 400]);
        ylim([-22 22])
        xlabel('Time (ms)');
        ylabel('Current (nAm)');
        set(gcf, 'Position', [680   652   573   344]);
        set(gca, 'FontSize', 24);
        legend({'Left', 'Right'}, 'location', 'eastoutside');

        if pp == 1
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'no_amm');
        elseif pp == 2
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'hfc');
        elseif pp == 3
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'hfc_with_gradients');
        elseif pp == 4
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'amm_spatial');
        else
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'amm');
        end

        save_name = sprintf('%s_ROI_dipole', extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'));
        print(fullfile(save_loc, save_name),'-dpng','-r300');
    end
    close all
end

%% Fieldtrip dipole fit to plot evoked response, standards and deviants at source level

clearvars -except meta_data colormap123 delay

rng(76);

for recording = 1:size(meta_data,1)
    cd(meta_data{recording, "analysed_data_loc"});
    start_string = {'be_ffft_', 'be_hffft_', 'be_m2ffft_', 'be_mffft_'}; %'be_h2ffft_', 
    for pp = 1:length(start_string)
        DD{pp} = spm_eeg_load(char(fullfile(meta_data{recording, "analysed_data_loc"}, ...
            strcat(start_string{pp}, extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'), '.mat'))));
    end

    % Get deviant and standard labels
    stim = readtable(meta_data{recording, "stim_data_fname"});
    deviants = find(contains(stim.Condition, 'deviant'));
    trial_length = diff(deviants);
    trial_length = cat(1, trial_length, size(DD{pp},3) - max(deviants) + 1);

    % Take last tone of each set as standard
    standards = deviants(2:end)-1;
    standards = cat(1, standards, size(DD{pp},3));
    
    % Dipole fit based on AMM results
    ctx = gifti(DD{4}.inv{1}.mesh.tess_ctx);

    % Find initialisation positions
    linds = 1:size(ctx.vertices,1)/2;
    rinds = size(ctx.vertices,1)/2+1:size(ctx.vertices,1);
    [~, maxind_l] = max(DD{4}.inv{1}.contrast.GW{1}(linds));
    [~, maxind_r] = max(DD{4}.inv{1}.contrast.GW{1}(rinds));
    maxind_l = linds(maxind_l);
    maxind_r = rinds(maxind_r);

    % Format data for fieldtrip
    data = ftraw(DD{4});
    cfg = [];
    cfg.trials = indtrial(DD{4}, 'tone', 'GOOD');
    data = ft_selectdata(cfg, data);

    % Average
    tl_data = ft_timelockanalysis([], data);

    % Prepare headmodel
    mesh = ft_read_headshape(DD{4}.inv{1}.mesh.tess_iskull);
    cfg = [];
    cfg.method = 'singleshell';
    headmodel = ft_prepare_headmodel(cfg, mesh);

    % Prepare sourcemodel - use cortical mesh
    cfg = [];
    cfg.method = 'basedoncortex';
    cfg.headshape = ft_read_headshape(DD{4}.inv{1}.mesh.tess_ctx);
    cfg.headmodel = headmodel;
    cfg.inwardshift = 0;
    sourcemodel = ft_prepare_sourcemodel(cfg);

    % Prepare leadfields
    cfg                  = [];
    cfg.grad             = data.grad;
    cfg.headmodel        = headmodel;
    cfg.reducerank       = 2;
    cfg.channel          = DD{4}.chanlabels(indchantype(DD{pp}, 'MEGMAG', 'GOOD'));
    cfg.sourcemodel = sourcemodel;
    sourcemodel = ft_prepare_leadfield(cfg);

    % Dipole fit
    cfg = [];
    cfg.latency = [0.08 0.12]+delay*1e-3;
    cfg.numdipoles = 2;
    cfg.symmetry = [];
    cfg.gridsearch = 'no';
    cfg.dip.pos = double(ctx.vertices([maxind_l, maxind_r], :));
    cfg.headmodel = headmodel;
    cfg.sourcemodel = sourcemodel;
    cfg.channel = DD{4}.chanlabels(indchantype(DD{4}, 'MEGMAG', 'GOOD'));
    cfg.senstype = 'meg';
    source = ft_dipolefitting(cfg, tl_data);

    % Plot dipole position
    figure; hold on;
    trisurf(sourcemodel.tri, sourcemodel.pos(:,1), sourcemodel.pos(:,2), sourcemodel.pos(:,3), 'FaceColor', [0.4 0.4 0.4], 'FaceAlpha', 0.3, 'EdgeColor', 'None');
    daspect([1 1 1])
    ft_plot_dipole(source.dip.pos(1,:), mean(source.dip.mom(1:3,:),2), 'color', 'b', 'unit', 'mm')
    ft_plot_dipole(source.dip.pos(2,:), mean(source.dip.mom(4:6,:),2), 'color', 'b', 'unit', 'mm')
    set(gcf, 'Position', [-1380, 320, 560, 420]);

    view(-90,0);
    save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'amm');

    save_name = sprintf('%s_ft_dip_fit_left', extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'));
    print(fullfile(save_loc, save_name),'-dpng','-r300');

    view([90 0]);

    save_name = sprintf('%s_ft_dip_fit_right', extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'));
    print(fullfile(save_loc, save_name),'-dpng','-r300');

    
    % combined plot for supplementary material
    % Left
    fig_left = figure;
    subplot(5, 4, 1:12);
    trisurf(sourcemodel.tri, sourcemodel.pos(:,1), sourcemodel.pos(:,2), sourcemodel.pos(:,3), 'FaceColor', [0.4 0.4 0.4], 'FaceAlpha', 0.3, 'EdgeColor', 'None');
    daspect([1 1 1])
    ft_plot_dipole(source.dip.pos(1,:), mean(source.dip.mom(1:3,:),2), 'color', 'b', 'unit', 'mm')
    view(-90,0);
    axs_left = gobjects(1,4);

    % Right
    fig_right = figure;
    subplot(5, 4, 1:12);
    trisurf(sourcemodel.tri, sourcemodel.pos(:,1), sourcemodel.pos(:,2), sourcemodel.pos(:,3), 'FaceColor', [0.4 0.4 0.4], 'FaceAlpha', 0.3, 'EdgeColor', 'None');
    daspect([1 1 1])
    ft_plot_dipole(source.dip.pos(2,:), mean(source.dip.mom(4:6,:),2), 'color', 'b', 'unit', 'mm')
    view([90 0]);
    axs_right = gobjects(1,4);

    % Find cortex point
    for ind = 1:2
        [~, sourceind(ind)] = min(sqrt(sum((ctx.vertices - source.dip.pos(ind,:)).^2, 2)));
    end
    sourceind = sort(sourceind);

    % Plot estimated source current for each preprocessing step
    for pp = 1:length(DD)

        L = full(spm_eeg_lgainmat(DD{pp},sourceind));

        good_trials = indtrial(DD{pp}, 'tone', 'GOOD');
        X_evoked = zeros(size(L,2), size(DD{pp},2), length(good_trials));
        X_standards = zeros(size(L,2), size(DD{pp},2), length(standards));
        X_deviants = zeros(size(L,2), size(DD{pp},2), length(deviants));

        for tt = 1:length(good_trials)
            X_evoked(:,:,tt) = pinv(L)*DD{pp}(indchantype(DD{pp},'MEGMAG','GOOD'),:,good_trials(tt));
        end
        for tt = 1:length(deviants)
            X_deviants(:,:,tt) = pinv(L)*DD{pp}(indchantype(DD{pp},'MEGMAG','GOOD'),:,deviants(tt));
        end
        for tt = 1:length(standards)
            X_standards(:,:,tt) = pinv(L)*DD{pp}(indchantype(DD{pp},'MEGMAG','GOOD'),:,standards(tt));
        end

        % T-test across trials
        SE_evoked = std(X_evoked, [], 3)./sqrt(size(X_evoked,3));
        t_evoked = mean(X_evoked, 3)./SE_evoked;
        SE_standards = std(X_standards, [], 3)./sqrt(size(X_standards, 3));
        t_standards = mean(X_standards,3)./SE_standards;
        SE_deviants = std(X_deviants, [], 3)./sqrt(size(X_deviants, 3));
        t_deviants = mean(X_deviants,3)./SE_deviants;
        
        % Unpaired t-test equal variance between deviants and standards for MMN response
        n1 = size(X_deviants,3);
        n2 = size(X_standards,3);
        SE = sqrt(((n1-1)*std(X_deviants, [], 3).^2 + (n2-1)*std(X_standards, [], 3).^2)./(n1 + n2 - 2))*...
            sqrt(1/n1 + 1/n2);
        t_diff = (mean(X_deviants,3) - mean(X_standards,3))./SE;


        % Evoked response:
        figure; 
        t = tiledlayout(1,2);
        yl = [];
        for ind = 1:2
            nexttile(t); hold on; grid on; box on;

            % plot(DD{pp}.time*1e3 - delay, mean(X_evoked(ind,:,:), 3), 'LineWidth', 3, 'LineStyle', '-', 'color', 'k');
            plot(DD{pp}.time*1e3 - delay, t_standards(ind,:), 'LineWidth', 3, 'LineStyle', '-', 'color', 'k');

            set(gca, 'FontSize', 18);
            xlim([-100 400]);
            yl(ind) = max(abs(ylim));
            xlabel('Time (ms)');

            if ind == 1
                title('Left Hemisphere');
            else
                title('Right Hemisphere');
            end

        end

        % Set axes limits and legend
        % ylim(t.Children, [-1 1]*max(abs(yl)));
        ylim(t.Children, [-1 1]*7);
        % ylabel(t, {'Estimated Source', 'Current (nAm)'}, 'FontSize', 18);
        ylabel(t, 't-stat', 'FontSize', 18);
        legend('Standards', 'location', 'eastoutside');
        set(gcf, 'Position', [626   476   821   285]);

        if pp == 1
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'no_amm');
        elseif pp == 2
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'hfc');
        % elseif pp == 3
        %     save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'hfc_with_gradients');
        elseif pp == 3
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'amm_spatial');
        else
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'amm');
        end

        save_name = sprintf('%s_ROI_evoked_dipole', extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'));
        print(fullfile(save_loc, save_name),'-dpng','-r300');

        % MMN:
        figure; 
        t = tiledlayout(1,2);
        yl = [];
        for ind = 1:2
            nexttile(t); hold on; grid on; box on;

            % Standards
            % plot(DD{pp}.time*1e3 - delay, mean(X_standards(ind,:,:), 3), 'LineWidth', 3, 'LineStyle', '--', 'color', [0.3639    0.5755    0.7484]);
            plot(DD{pp}.time*1e3 - delay, t_standards(ind,:), 'LineWidth', 3, 'LineStyle', '--', 'color', [0.3639    0.5755    0.7484]);

            % Deviants
            % plot(DD{pp}.time*1e3 - delay, mean(X_deviants(ind,:,:), 3), 'LineWidth', 3, 'LineStyle', ':', 'color', [0.9153    0.2816    0.2878]);
            plot(DD{pp}.time*1e3 - delay, t_deviants(ind,:), 'LineWidth', 3, 'LineStyle', ':', 'color', [0.9153    0.2816    0.2878]);

            % Difference
            % plot(DD{pp}.time*1e3 - delay, mean(X_deviants(ind,:,:),3) - mean(X_standards(ind,:,:), 3), 'LineWidth', 3, 'LineStyle', '-', 'color', [0.3373    0.3020    0.2902]);
            plot(DD{pp}.time*1e3 - delay, t_diff(ind,:), 'LineWidth', 3, 'LineStyle', '-', 'color', [0.3373    0.3020    0.2902]);
            set(gca, 'FontSize', 18);
            xlim([-100 400]);
            % yl(ind) = max(abs(ylim));
            xlabel('Time (ms)');

            if ind == 1
                title('Left Hemisphere');
            else
                title('Right Hemisphere');
            end

        end

        % Set axes limits and legend
        % ylim(t.Children, [-1 1]*max(abs(yl)));
        ylim(t.Children, [-1 1]*10);
        % ylabel(t, {'Estimated Source', 'Current (nAm)'}, 'FontSize', 18);
        ylabel(t, 't-stat', 'FontSize', 18);
        lgd = legend('Standards', 'Deviants', 'MMN', 'location', 'eastoutside');
        set(gcf, 'Position', [626   476   821   285]);

        % Add text to indicate how many trials per condition
        
        annotation('textbox', [lgd.Position(1), lgd.Position(2) - 0.35, lgd.Position(3), 0.3], ...
            'string', sprintf('# deviants: %.f\n# standards: %.f', length(deviants), length(standards)), 'FontSize', 16, 'EdgeColor', 'None');

        if pp == 1
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'no_amm');
        elseif pp == 2
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'hfc');
        % elseif pp == 3
        %     save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'hfc_with_gradients');
        elseif pp == 3
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'amm_spatial');
        else
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'amm');
        end

        save_name = sprintf('%s_ROI_MMN_trace_dipole_all_sets', extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'));
        print(fullfile(save_loc, save_name),'-dpng','-r300');

        % Combined plot
        % Left
        % nexttile(tcomb_left); hold on; grid on; box on;
        figure(fig_left);
        axs_left(pp) = subplot(5, 4, [12+pp, 16+pp]); hold on; grid on; box on;
        % plot(DD{pp}.time*1e3 - delay, mean(X_standards(1,:,:), 3), 'LineWidth', 3, 'LineStyle', '--', 'color', [0.3639    0.5755    0.7484]);
        % plot(DD{pp}.time*1e3 - delay, mean(X_deviants(1,:,:), 3), 'LineWidth', 3, 'LineStyle', ':', 'color', [0.9153    0.2816    0.2878]);
        % plot(DD{pp}.time*1e3 - delay, mean(X_deviants(1,:,:),3) - mean(X_standards(1,:,:), 3), 'LineWidth', 3, 'LineStyle', '-', 'color', [0.3373    0.3020    0.2902]);
        plot(DD{pp}.time*1e3 - delay, t_standards(1,:), 'LineWidth', 3, 'LineStyle', '--', 'color', [0.3639    0.5755    0.7484]);
        plot(DD{pp}.time*1e3 - delay, t_deviants(1,:), 'LineWidth', 3, 'LineStyle', ':', 'color', [0.9153    0.2816    0.2878]);
        plot(DD{pp}.time*1e3 - delay, t_diff(1,:), 'LineWidth', 3, 'LineStyle', '-', 'color', [0.3373    0.3020    0.2902]);
        set(gca, 'FontSize', 22);
        if pp == 1
            % ylabel(gca, {'Estimated Source', 'Current (nAm)'});
            ylabel(gca, 't-stat');
        else
            set(gca,'ytick',[]);
        end
        % Right
        % nexttile(tcomb_right); hold on; grid on; box on;
        figure(fig_right);
        axs_right(pp) = subplot(5, 4, [12+pp, 16+pp]); hold on; grid on; box on;
        % plot(DD{pp}.time*1e3 - delay, mean(X_standards(2,:,:), 3), 'LineWidth', 3, 'LineStyle', '--', 'color', [0.3639    0.5755    0.7484]);
        % plot(DD{pp}.time*1e3 - delay, mean(X_deviants(2,:,:), 3), 'LineWidth', 3, 'LineStyle', ':', 'color', [0.9153    0.2816    0.2878]);
        % plot(DD{pp}.time*1e3 - delay, mean(X_deviants(2,:,:),3) - mean(X_standards(2,:,:), 3), 'LineWidth', 3, 'LineStyle', '-', 'color', [0.3373    0.3020    0.2902]);
        plot(DD{pp}.time*1e3 - delay, t_standards(2,:), 'LineWidth', 3, 'LineStyle', '--', 'color', [0.3639    0.5755    0.7484]);
        plot(DD{pp}.time*1e3 - delay, t_deviants(2,:), 'LineWidth', 3, 'LineStyle', ':', 'color', [0.9153    0.2816    0.2878]);
        plot(DD{pp}.time*1e3 - delay, t_diff(2,:), 'LineWidth', 3, 'LineStyle', '-', 'color', [0.3373    0.3020    0.2902]);
        set(gca, 'FontSize', 22);
        if pp == 1
            % ylabel(gca, {'Estimated Source', 'Current (nAm)'});
            ylabel(gca, 't-stat');
        else
            set(gca,'ytick',[]);
        end
    end
    
    title(axs_left(1), 'No Filter');
    title(axs_left(2), 'HFC');
    title(axs_left(3), 'Spatial only AMM');
    title(axs_left(4), {'AMM with temporal', 'extension'});

    title(axs_right(1), 'No Filter');
    title(axs_right(2), 'HFC');
    title(axs_right(3), 'Spatial only AMM');
    title(axs_right(4), {'AMM with temporal', 'extension'});

    linkaxes(axs_left);
    xlim(axs_left(1), [-100 400]);
    xlabel(axs_left(3), 'Time (ms)');
    % yl = ylim(axs_left(1));
    % ylim(axs_left(1), [-1 1]*max(abs(yl)));
    ylim(axs_left(1), [-1 1]*10);

    linkaxes(axs_right);
    xlim(axs_right(1), [-100 400]);
    xlabel(axs_right(3), 'Time (ms)');
    % yl = ylim(axs_right(1));
    % ylim(axs_right(1), [-1 1]*max(abs(yl)));
    ylim(axs_right(1), [-1 1]*10);

    set(fig_left, 'Position', [-1860         153        1855         750])
    set(fig_right, 'Position', [-1860         153        1855         750])

    figure(fig_left);
    save_loc = fullfile(meta_data{recording, "results_save_loc"});
    save_name = sprintf('%s_ROI_MMN_trace_dipole_combined_left_all_sets', extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'));
    print(fullfile(save_loc, save_name),'-dpng','-r300');

    figure(fig_right);
    save_loc = fullfile(meta_data{recording, "results_save_loc"});
    save_name = sprintf('%s_ROI_MMN_trace_dipole_combined_right_all_sets', extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'));
    print(fullfile(save_loc, save_name),'-dpng','-r300');

    close all
end

%% Find spatiotemporal clusters at source-level - use source stats and power

% clearvars -except meta_data colormap123 delay
% 
% for recording = 1:size(meta_data,1)
% 
%     fprintf('Recording no: %.f of %.f\n', recording, size(meta_data,1));
% 
%     % Get data
%     cd(meta_data{recording, "analysed_data_loc"});
%     start_string = {'e_ffft_', 'e_hffft_', 'e_m2ffft_', 'e_mffft_'}; %'e_h2ffft_', 
%     for pp = 1:length(start_string)
%         DD{pp} = spm_eeg_load(char(fullfile(meta_data{recording, "analysed_data_loc"}, ...
%             strcat(start_string{pp}, extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'), '.mat'))));
%     end
% 
%     ctx = export(gifti(DD{1}.inv{1}.mesh.tess_ctx),'ft');
% 
%     % Create infalted cortex mesh
%     if ~isfile(fullfile(char(extractBefore(meta_data{recording, "raw_data_loc"},'meg')),...
%         'anat', sprintf('%s_inflated_cortex.gii', meta_data{recording, 'sub'})))
%         inflated_ctx_fname = fullfile(char(extractBefore(meta_data{recording, "raw_data_loc"},'meg')),...
%                 'anat', sprintf('%s.L.inflated.4k_fs_LR.surf.gii', meta_data{recording, 'sub'}));
%         combine_surfaces({inflated_ctx_fname, strrep(inflated_ctx_fname, '.L.', '.R.')}, ...
%             fullfile(char(extractBefore(meta_data{recording, "raw_data_loc"},'meg')),...
%             'anat', sprintf('%s_inflated_cortex.gii', meta_data{recording, 'sub'})));
%     end
% 
%     % Source reconstruct (MNE)
%     for pp = 1:length(DD)
% 
%         matlabbatch = [];
%         matlabbatch{1}.spm.meeg.source.invert.D = {fullfile(DD{pp})};
%         matlabbatch{1}.spm.meeg.source.invert.val = 1;
%         matlabbatch{1}.spm.meeg.source.invert.whatconditions.all = 1;
%         matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.invtype = 'IID';
%         matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.woi = [-Inf Inf];
%         matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.foi = [2 40];
%         matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.hanning = 1;
%         matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.priors.priorsmask = {''};
%         matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.priors.space = 1;
%         matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.restrict.locs = zeros(0, 3);
%         matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.restrict.radius = 32;
%         matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.restrict.mask = {''};
%         matlabbatch{1}.spm.meeg.source.invert.modality = {'All'};
%         matlabbatch{2}.spm.meeg.source.results.D(1) = cfg_dep('Source inversion: M/EEG dataset(s) after imaging source reconstruction', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','D'));
%         matlabbatch{2}.spm.meeg.source.results.val = 1;
%         matlabbatch{2}.spm.meeg.source.results.woi = [50 150]+delay;
%         matlabbatch{2}.spm.meeg.source.results.foi = [2 40];
%         matlabbatch{2}.spm.meeg.source.results.ctype = 'trials';
%         matlabbatch{2}.spm.meeg.source.results.space = 0;
%         matlabbatch{2}.spm.meeg.source.results.format = 'mesh';
%         matlabbatch{2}.spm.meeg.source.results.smoothing = 8;
%         matlabbatch{3}.spm.meeg.source.results.D(1) = cfg_dep('Source inversion: M/EEG dataset(s) after imaging source reconstruction', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','D'));
%         matlabbatch{3}.spm.meeg.source.results.val = 1;
%         matlabbatch{3}.spm.meeg.source.results.woi = [-150 -50]+delay;
%         matlabbatch{3}.spm.meeg.source.results.foi = [2 40];
%         matlabbatch{3}.spm.meeg.source.results.ctype = 'trials';
%         matlabbatch{3}.spm.meeg.source.results.space = 0;
%         matlabbatch{3}.spm.meeg.source.results.format = 'mesh';
%         matlabbatch{3}.spm.meeg.source.results.smoothing = 8;
% 
%         a = spm_jobman('run',matlabbatch);
% 
%         % Create baseline data structure
%         source_data_baseline = cell(length(a{3}.files), 1);
%         for gifti_file_ind = 1:length(a{3}.files) 
%             source_data_baseline{gifti_file_ind} = [];
%             source_data_baseline{gifti_file_ind}.pos = ctx.pnt;
%             source_data_baseline{gifti_file_ind}.inside = ones(size(source_data_baseline{gifti_file_ind}.pos,1), 1)==1;
%             source_data_baseline{gifti_file_ind}.pow = export(gifti(fullfile(a{3}.files{gifti_file_ind})), ...
%                 'patch').facevertexcdata;
%         end
% 
%         % Create time of interest data structure
%         source_data_toi = cell(length(a{2}.files), 1);
%         for gifti_file_ind = 1:length(a{2}.files) 
%             source_data_toi{gifti_file_ind} = [];
%             source_data_toi{gifti_file_ind}.pos = ctx.pnt;
%             source_data_toi{gifti_file_ind}.inside = ones(size(source_data_toi{gifti_file_ind}.pos,1), 1)==1;
%             source_data_toi{gifti_file_ind}.pow = export(gifti(fullfile(a{2}.files{gifti_file_ind})), ...
%                 'patch').facevertexcdata;
%         end
% 
%         % Perform Statistical Analysis
%         cfg                     = [];
%         cfg.method              = 'montecarlo';
%         cfg.statistic           = 'ft_statfun_indepsamplesT';
%         cfg.parameter           = 'pow';
%         cfg.correctm            = 'cluster';
%         cfg.alpha               = 0.05;
%         cfg.numrandomization    = 1000;
%         cfg.tail                = 1;
% 
%         % Design Matrix
%         ntrials                 = numel(source_data_toi);
%         %cfg.design(1,:)         = [1:ntrials 1:ntrials];
%         cfg.design(1,:)         = [ones(1,ntrials) ones(1,ntrials)*2];
% 
%         % row of design matrix that contains unit variable (in this case: trials)
%         %cfg.uvar                = 1;
%         % row of design matrix that contains independent variable (the conditions)
%         cfg.ivar                = 1; 
% 
%         % Perform statistical analysis
%         [stat]                  = ft_sourcestatistics(cfg,source_data_toi{:}, source_data_baseline{:});
% 
%         % Show raw source level statistics (2D plot)
%         inflated_ctx = ft_read_headshape(fullfile(char(extractBefore(meta_data{recording, "raw_data_loc"},'meg')),...
%                 'anat', sprintf('%s_inflated_cortex.gii', meta_data{recording, 'sub'})));
% 
%         stat.tri = ctx.tri;
%         cfg                     = [];
%         cfg.method              = 'surface';
%         cfg.facecolor           = [0.4 0.4 0.4];
%         cfg.vertexcolor         = 'none';
%         cfg.funparameter        = 'stat';
%         cfg.location            = 'max';
%         % cfg.maskparameter       = 'mask';
%         ft_sourceplot(cfg, stat);
%         colormap(colormap123) % change the colormap
% 
%     end
% end

%% Find spatiotemporal clusters at source-level 

% clearvars -except meta_data colormap123 delay
% 
% for recording = 1:size(meta_data,1)
% 
%     fprintf('Recording no: %.f of %.f\n', recording, size(meta_data,1));
% 
%     % Get data
%     cd(meta_data{recording, "analysed_data_loc"});
%     start_string = {'e_ffft_', 'e_hffft_', 'e_m2ffft_', 'e_mffft_'}; %'e_h2ffft_', 
%     for pp = 1:length(start_string)
%         DD{pp} = spm_eeg_load(char(fullfile(meta_data{recording, "analysed_data_loc"}, ...
%             strcat(start_string{pp}, extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'), '.mat'))));
%     end
% 
%     ctx = export(gifti(DD{1}.inv{1}.mesh.tess_ctx),'ft');
% 
%     % Create infalted cortex mesh
%     if ~isfile(fullfile(char(extractBefore(meta_data{recording, "raw_data_loc"},'meg')),...
%         'anat', sprintf('%s_inflated_cortex.gii', meta_data{recording, 'sub'})))
%         inflated_ctx_fname = fullfile(char(extractBefore(meta_data{recording, "raw_data_loc"},'meg')),...
%                 'anat', sprintf('%s.L.inflated.4k_fs_LR.surf.gii', meta_data{recording, 'sub'}));
%         combine_surfaces({inflated_ctx_fname, strrep(inflated_ctx_fname, '.L.', '.R.')}, ...
%             fullfile(char(extractBefore(meta_data{recording, "raw_data_loc"},'meg')),...
%             'anat', sprintf('%s_inflated_cortex.gii', meta_data{recording, 'sub'})));
%     end
% 
%     % Create a FieldTrip structure to fill
%     source_data = [];
%     source_data.label = cellstr(num2str((1:size(ctx.pnt,1))'));
%     elec = [];
%     elec.label = source_data.label;
%     elec.elecpos = ctx.pnt;
%     elec.unit = 'mm';
%     elec.chanpos = elec.elecpos;
%     source_data.elec = elec;
% 
%     hdr = [];
%     hdr.Fs          = fsample(DD{pp});
%     hdr.nChans      = size(ctx.pnt,1);
%     hdr.label       = source_data.label;
%     hdr.chanunit    = repmat({'nAm'}, 1, size(ctx.pnt,1));
%     hdr.chantype = repmat({'virtualelec'}, 1, size(ctx.pnt,1));
%     source_data.hdr = hdr;
%     source_data.fsample = DD{pp}.fsample;
% 
%     % Find neighbours of mesh points
%     cfg = [];
%     cfg.method = 'distance';
%     cfg.neighbourdist = 7;
%     neighbours = ft_prepare_neighbours(cfg, source_data);
% 
%     % Source reconstruct
%     for pp = 1:length(DD)
% 
%         % Create virtual electrodes at each mesh point
%         % if ~isfile(fullfile('E:\Data\Neuro1\Auditory\anonymised_for_sharing\analysedData', ...
%             % meta_data{recording, 'sub'}, 'current_distributions_for_mne', DD{pp}.fname))
%             matlabbatch = [];
%             matlabbatch{1}.spm.meeg.source.invert.D = {DD{pp}.fname};
%             matlabbatch{1}.spm.meeg.source.invert.val = 1;
%             matlabbatch{1}.spm.meeg.source.invert.whatconditions.all = 1;
%             matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.invtype = 'IID';
%             matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.woi = [-Inf Inf];
%             matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.foi = [2 40];
%             matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.hanning = 1;
%             matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.priors.priorsmask = {''};
%             matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.priors.space = 1;
%             matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.restrict.locs = zeros(0, 3);
%             matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.restrict.radius = 32;
%             matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.restrict.mask = {''};
%             matlabbatch{1}.spm.meeg.source.invert.modality = {'All'};
%             matlabbatch{2}.spm.meeg.source.results.D(1) = cfg_dep('Source inversion: M/EEG dataset(s) after imaging source reconstruction', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','D'));
%             matlabbatch{2}.spm.meeg.source.results.val = 1;
%             matlabbatch{2}.spm.meeg.source.results.woi = [50 150];
%             matlabbatch{2}.spm.meeg.source.results.foi = [2 40];
%             matlabbatch{2}.spm.meeg.source.results.ctype = 'evoked';
%             matlabbatch{2}.spm.meeg.source.results.space = 0;
%             matlabbatch{2}.spm.meeg.source.results.format = 'mesh';
%             matlabbatch{2}.spm.meeg.source.results.smoothing = 8;
%             matlabbatch{2}.spm.meeg.source.results.D(1) = cfg_dep('Source inversion: M/EEG dataset(s) after imaging source reconstruction', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','D'));
%             matlabbatch{2}.spm.meeg.source.results.val = 1;
%             matlabbatch{2}.spm.meeg.source.results.woi = [-150 -50];
%             matlabbatch{2}.spm.meeg.source.results.foi = [2 40];
%             matlabbatch{2}.spm.meeg.source.results.ctype = 'evoked';
%             matlabbatch{2}.spm.meeg.source.results.space = 0;
%             matlabbatch{2}.spm.meeg.source.results.format = 'mesh';
%             matlabbatch{2}.spm.meeg.source.results.smoothing = 8;
% 
%             a = spm_jobman('run',matlabbatch);
%             DD{pp} = spm_eeg_load(fullfile(DD{pp}));
% 
%             toi = 100;
%             delay = 10; % ms - neuro-1 delay between truth and recording
%             [~, toi_ind] = min(abs(DD{pp}.time - toi*1e-3 - delay*1e-3));
% 
%             % Estimate current
%             Y = DD{pp}(indchantype(DD{pp},'MEGMAG', 'GOOD'), ...
%                 toi_ind-floor(15e-3*DD{pp}.fsample):toi_ind+floor(15e-3*DD{pp}.fsample),...
%                 indtrial(DD{pp}, DD{pp}.condlist{1}, 'GOOD'));
%             J = zeros(size(DD{pp}.inv{1}.inverse.M,1), size(Y,2), size(Y,3));
%             for tt = 1:size(Y,2)
%                 UY = DD{pp}.inv{1}.inverse.U{1}*squeeze(Y(:,tt,:))*DD{pp}.inv{1}.inverse.scale/size(Y,3); % Spatial projector
%                 J(:,tt,:) = DD{pp}.inv{1}.inverse.M*UY; % MAP projector
%             end
%             clear Y UY
% 
% 
%             % Save data to avoid recreating
%             if ~isfolder(fullfile('E:\Data\Neuro1\Auditory\anonymised_for_sharing\analysedData', ...
%                     meta_data{recording, 'sub'}, 'current_distributions_for_mne2'))
%                 mkdir(fullfile('E:\Data\Neuro1\Auditory\anonymised_for_sharing\analysedData', ...
%                     meta_data{recording, 'sub'}, 'current_distributions_for_mne2'));
%             end
%             save(fullfile('E:\Data\Neuro1\Auditory\anonymised_for_sharing\analysedData', ...
%                 meta_data{recording, 'sub'}, 'current_distributions_for_mne2', DD{pp}.fname), 'J', '-v7.3');
% 
% 
%         % else
%         %     load(fullfile('E:\Data\Neuro1\Auditory\anonymised_for_sharing\analysedData', ...
%         %         meta_data{recording, 'sub'}, 'current_distributions_for_mne', DD{pp}.fname));
%         % end
% 
%         % Write into source_data variable
%         source_data.trial = cell(1, size(J,3));
%         for jj = 1:size(J,3)
%             source_data.trial{jj} = J(:,:,jj);
%         end
%         clear J ctx;
% 
%         % Update time variables of source_data
%         source_data.time = repmat({DD{pp}.time(toi_ind-floor(15e-3*DD{pp}.fsample):toi_ind+floor(15e-3*DD{pp}.fsample))}, 1, length(source_data.trial));
%         timeind = toi_ind-floor(15e-3*DD{pp}.fsample):toi_ind+floor(15e-3*DD{pp}.fsample);
%         trialind = indtrial(DD{pp}, 'tone', 'GOOD');
%         onsets = trialonset(DD{pp}, trialind);
%         if all(onsets>0)
%             onsets = round(onsets(:)*fsample(DD{pp}));
%             source_data.sampleinfo = [onsets+timeind(1) onsets+timeind(end)]-1;
%         end
%         source_data.hdr.nSamples    = length(timeind);
%         source_data.hdr.nSamplesPre = sum(time(DD{pp}, timeind)<0);
%         source_data.hdr.nTrials     = length(trialind);
% 
%         % Create a zeros dataset for comparison
%         source_data_zeros = source_data;
%         source_data_zeros.trial = repmat({zeros(size(source_data.trial{1}))}, 1, size(source_data.trial,2));
% 
%         % Timelock
%         cfg = [];
%         cfg.keeptrials = 'yes';
%         tl_source = ft_timelockanalysis(cfg, source_data);
%         tl_source_zeros = ft_timelockanalysis(cfg, source_data_zeros);
% 
%         % Permutation test
%         cfg = [];
%         cfg.method = 'montecarlo';
%         cfg.statistic = 'depsamplesT';
%         cfg.correctm = 'cluster';
%         cfg.clusteralpha = 0.05;
%         cfg.clusterstatistic = 'maxsum';
%         cfg.minnbchan = 0; 
%         cfg.neighbours = neighbours;
%         cfg.tail = 0;
%         % cfg.clustertail      = 0;
%         cfg.alpha = 0.025;
%         cfg.numrandomization = 100;
% 
%         n_zeros  = size(tl_source_zeros.trial, 1);
%         n_toi = size(tl_source.trial, 1);
% 
%         cfg.design = [ones(1,n_zeros), ones(1,n_toi)*2; 1:n_zeros, 1:n_toi];
%         cfg.ivar = 1;
%         cfg.channel = 'all';
%         [stat] = ft_timelockstatistics(cfg, tl_source, tl_source_zeros);
% 
%         % Make a vector of all p-values associated with the clusters from ft_timelockstatistics.
%         pos_cluster_pvals = [stat.posclusters(:).prob];
% 
%         % Then, find which clusters are deemed interesting to visualize
%         pos_clust = find(pos_cluster_pvals < 0.025);
%         pos = ismember(stat.posclusterslabelmat, pos_clust);
% 
%         % and now for the negative clusters...
%         neg_cluster_pvals = [stat.negclusters(:).prob];
%         neg_clust = find(neg_cluster_pvals < 0.025);
%         neg = ismember(stat.negclusterslabelmat, neg_clust);
% 
%         % Plot significant source power
% 
%         figure;
%         pos_int = any(pos, 2);
%         neg_int = any(neg, 2);
% 
%         inflated_ctx = ft_read_headshape(fullfile(char(extractBefore(meta_data{recording, "raw_data_loc"},'meg')),...
%                 'anat', sprintf('%s_inflated_cortex.gii', meta_data{recording, 'sub'})));
% 
%         cfg = [];
%         cfg.method = 'surface';
%         cfg.facecolor = [0.4 0.4 0.4];
%         cfg.vertexcolor = 'none';
% 
%        if any(pos_int | neg_int)
% 
%             % Plot
% 
%             inflated_ctx.pow = export(gifti(a{2}.files{1}), 'patch').facevertexcdata;
%             inflated_ctx.mask = any(pos_int | neg_int, 2);
% 
%             %     cfg.funcolorlim    = [0, 15];
%             cfg.funparameter   = 'pow';
%             cfg.maskparameter  = 'mask';
%             cfg.funcolormap    = 'hot';
%             %cfg.funcolorlim    = [0, max(data_summary(:,1))*1.1];
%             cfg.colorbartext = 'Source Power (a.u.)';
%             ft_sourceplot(cfg, inflated_ctx);
%         else
%             surf.pos = inflated_ctx.pos;
%             surf.tri = inflated_ctx.tri;
%             ft_plot_mesh(surf,'edgecolor', 'none', 'facecolor', cfg.facecolor, 'vertexcolor', cfg.vertexcolor);
%             lighting gouraud
%             camlight
%         end
% 
%         view ([-90 5])             % rotate the object in the view
%         cl1 = camlight('headlight');
%         set(gcf, 'color', 'w');
%         set(gca, 'FontSize', 26);
%         material dull
% 
%         if pp == 1
%             save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'no_amm');
%         elseif pp == 2
%             save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'hfc');
%         elseif pp == 3
%             save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'hfc_with_gradients');
%         else
%             save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'amm');
%         end
% 
%         save_name = sprintf('%s_min_norm_clusters_left', extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'));
%         print(fullfile(save_loc, save_name),'-dpng','-r300');
% 
%         view ([90 0])             % rotate the object in the view
%         cl2 = camlight('headlight');
%         cl2.Color = 0.4*ones(1,3);
% 
%         save_name = sprintf('%s_min_norm_clusters_right', extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'));
%         print(fullfile(save_loc, save_name),'-dpng','-r300');
% 
%     end
% end


%% Plot time series

for recording = 1:size(meta_data,1)
    
    if strcmp(meta_data{recording, "sub"}, 'sub-003')
        rad_ax = 'Z';
    else
        rad_ax = 'Y';
    end

    start_string = {'e_ffft_', 'e_hffft_', 'e_h2ffft_', 'e_m2ffft_', 'e_mffft_'};
    for pp = 1:length(start_string)
        DD{pp} = spm_eeg_load(char(fullfile(meta_data{recording, "analysed_data_loc"}, ...
            strcat(start_string{pp}, extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'), '.mat'))));
    end

    for pp = 1:length(DD)
        
        D = DD{pp};

        good_trials = indtrial(D, 'tone', 'GOOD');
        se = std(D(indchantype(D, 'MEGMAG', 'GOOD'),:,good_trials),[],3)./sqrt(length(good_trials));
        t = mean(D(indchantype(D, 'MEGMAG', 'GOOD'),:,good_trials),3)./se;
        
        [~, tind] = min(abs(D.time - 99*1e-3));
        max_colour = interp1([min(abs(t(:,tind))), max(abs(t(:,tind)))], [0.9, 0], abs(t(:,tind)));
        [~, plot_order] = sort(max_colour, 'descend');
        max_colour = repmat(max_colour, 1, 3);
    
        % Plot t stat
        figure; hold on; grid on; box on;
        for chan = 1:size(t,1)
            plot(1e3*D.time - delay, t(plot_order(chan),:), 'color', max_colour(plot_order(chan),:), 'LineWidth', 2);
        end
        a = tinv(1-0.025/(range(D.time)*40*length(plot_order)), size(D,3)-1);
        l1 = plot([-100 400], [a a], 'b--', 'LineWidth', 2);
        plot([-100 400], [-a -a], 'b--', 'LineWidth', 2);
        % legend(l1, 'Sig. Threshold')
        xlim([-100 400]);
        if contains(meta_data{recording, "raw_data_name"}, 'seat') || contains(meta_data{recording, "raw_data_name"}, 'Seat')
            ylim([-25 25]);
        else
            ylim([-1 1]*13);
        end
        xlabel('Time (ms)');
        ylabel('t-stat');
        set(gcf, 'Position', [680   654   451   344]);
        set(gca, 'FontSize', 24);
        fname = D.fname;
        
        if pp == 1
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'no_amm');
        elseif pp == 2
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'hfc');
        elseif pp == 3
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'hfc_with_gradients');
        elseif pp == 4
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'amm_spatial');
        else
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'amm');
        end
        if ~exist(save_loc, 'dir')
            mkdir(save_loc);
        end

        save_name = sprintf('%s_t_stat_time_series', extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'));
        print(fullfile(save_loc, save_name),'-dpng','-r300');
    
        % Plot topography of t-stat at 100 ms
        lay_name = fullfile(meta_data{recording, "analysed_data_loc"}, ...
            sprintf('%s_2Dlayout.mat', extractBefore(meta_data{recording, "raw_data_name"}, '_meg.lvm')));

        if isfile(lay_name)
            load(lay_name);
        else
            fid = fiducials(D);
            fid_struct = struct('NAS', fid.fid.pnt(contains(fid.fid.label, 'nas'),:), ...
                'LPA', fid.fid.pnt(contains(fid.fid.label, 'lpa'),:), ...
                'RPA', fid.fid.pnt(contains(fid.fid.label, 'rpa'),:));
            pos = D.sensors('MEG').coilpos;
            lay = spm_get_anatomical_layout(D.sensors('MEG').coilpos(endsWith(D.sensors('MEG').label, ['-', rad_ax]),:), ...
                D.sensors('MEG').label(endsWith(D.sensors('MEG').label, ['-', rad_ax])),...
                double(gifti(D.inv{1}.mesh.tess_scalp).vertices), fid_struct, 0);
            save(lay_name, 'lay');
        end
    
        data = ftraw(D);
        cfg = [];
        cfg.channel = intersect(data.grad.label, D.chanlabels(indchantype(D, 'MEGMAG', 'GOOD')));
        cfg.trials = good_trials;
        avdata = ft_timelockanalysis(cfg, data);
        tavdata = avdata;
        tavdata.avg = t;
    
        figure;
        cfg = [];
        cfg.layout    = lay;
        cfg.colorbar  = 'EastOutside';
        cfg.colorbartext = 't-stat (100 ms)';
        if contains(meta_data{recording, "raw_data_name"}, 'seat') || contains(meta_data{recording, "raw_data_name"}, 'Seat')
            cfg.zlim      = [-15 15];
        else
            cfg.zlim = [-1 1]*8.83;
        end
        cfg.colormap  = colormap123;
        cfg.xlim = [100, 100]*1e-3 + delay*1e-3;
        cfg.comment = 'no';
        cfg.figure = gca;
        set(gca, 'FontSize', 24);
        ft_topoplotER(cfg, tavdata)
        set(gcf, 'Position', [994   704   404   274]);

        save_name = sprintf('%s_t_stat_topography', extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'));
        print(fullfile(save_loc, save_name),'-dpng','-r300');
    
        % Plot average signal
        figure; hold on; grid on; box on;
        dat = mean(D(indchantype(D, 'MEGMAG', 'GOOD'),:,good_trials),3);
        for chan = 1:size(t,1)
            plot(1e3*D.time - delay, dat(plot_order(chan),:), 'color', max_colour(plot_order(chan),:), 'LineWidth', 2);
        end
        xlim([-100 400]);
        ylim([-550 550]);
        xlabel('Time (ms)');
        ylabel('B (fT)');
        grid on;
        set(gcf, 'Position', [680   654   451   344]);
        set(gca, 'FontSize', 24);
        save_name = sprintf('%s_average_time_series', extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'));
        print(fullfile(save_loc, save_name),'-dpng','-r300');
        
        % Topoplot
        figure;
        cfg.figure = gca;
        cfg.colorbartext = 'B (fT)';
        cfg.zlim      = [-350, 350];
        ft_topoplotER(cfg, avdata)
        set(gcf, 'Position', [994   704   404   274]);
        set(gca, 'FontSize', 24);

        save_name = sprintf('%s_average_topography', extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'));
        print(fullfile(save_loc, save_name),'-dpng','-r300');
    end

    close all
end

%% MMN

% At sensor level, just show channel with the largest t-stat evoked 
% response in seated closed loop recording. Mark on a topoplot. 

% --- choose channel

% Seated closed recordings
recording_order = {"sub-001_task-seatedClosed_meg.lvm",...
    "sub-002_task-seatedClosed_meg.lvm", "sub-003_task-seatedClosed_meg.lvm"};

% Channel to plot for each participant
chans_to_plot = zeros(1, length(recording_order));

for recording = 1:length(recording_order)

    rec_idx = find(contains(meta_data.raw_data_name, recording_order{recording}));
    D = spm_eeg_load(char(fullfile(meta_data{rec_idx, "analysed_data_loc"}, ...
        strcat('e_mffft_', extractBefore(recording_order{recording}, '.lvm'), '.mat'))));
    SE = std(D(indchantype(D, 'MEGMAG', 'GOOD'),:,indtrial(D, 'tone', 'GOOD')), [], 3)./sqrt(length(indtrial(D, 'tone', 'GOOD')));
    t = mean(D(indchantype(D, 'MEGMAG', 'GOOD'),:,indtrial(D, 'tone', 'GOOD')), 3)./SE;
    [~, chans_to_plot(recording)] = max(max(abs(t(:, D.time > (80+delay)*1e-3 & D.time < (120+delay)*1e-3)), [], 2));

    % Plot position of channel
    load(fullfile(meta_data{rec_idx, "analysed_data_loc"}, ...
            sprintf('%s_2Dlayout.mat', extractBefore(meta_data{rec_idx, "raw_data_name"}, '_meg.lvm'))));
    figure; hold on;
    ft_plot_layout(lay, 'label', 'no');
    chanlabels = D.chanlabels(indchantype(D, 'MEGMAG', 'GOOD'));
    channame = chanlabels(chans_to_plot(recording));
    layoutindex = find(startsWith(lay.label, channame{1}(1:end-1)));
    width = lay.width(layoutindex);
    height = lay.height(layoutindex);
    X = lay.pos(layoutindex,1);
    Y = lay.pos(layoutindex,2);
    patch([X-width/2, X+width/2, X+width/2, X-width/2], [Y-height/2 Y-height/2 Y+height/2 Y+height/2], 'red', 'FaceAlpha', 0.7);

    % Save figure
    save_loc = fullfile(meta_data{rec_idx, "results_save_loc"});
    save_name = 'MNN_sensor_pos';
    print(fullfile(save_loc, save_name),'-dpng','-r300');
end


% --- loop through recordings to plot

for recording = 1:size(meta_data,1)
    
    if strcmp(meta_data{recording, "sub"}, 'sub-003')
        rad_ax = 'Z';
    else
        rad_ax = 'Y';
    end

    start_string = {'e_ffft_', 'e_hffft_', 'e_h2ffft_', 'e_m2ffft_', 'e_mffft_'};
    for pp = 1:length(start_string)
        DD{pp} = spm_eeg_load(char(fullfile(meta_data{recording, "analysed_data_loc"}, ...
            strcat(start_string{pp}, extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'), '.mat'))));
    end

    % Get deviant and standard labels
    stim = readtable(meta_data{recording, "stim_data_fname"});
    deviants = find(contains(stim.Condition, 'deviant'));
    trial_length = diff(deviants);
    trial_length = cat(1, trial_length, size(DD{pp},3) - max(deviants) + 1);

    % Take last tone of each set as standard
    standards = deviants(2:end)-1;
    standards = cat(1, standards, size(DD{pp},3));

    % Prep for dipole fit
    % Initialise at auditory cortices
    aud_mni = [-54 -14 8 1; 54 -14 8 1]'; % Auditory cortices in MNI space
    aud_nat = DD{1}.inv{1}.datareg.fromMNI*aud_mni;
    aud_nat = aud_nat(1:3,:)';
    aud_nat = aud_nat*1e-3; % convert to m

    % Prepare headmodel
    mesh = ft_read_headshape(DD{1}.inv{1}.mesh.tess_iskull);
    cfg = [];
    cfg.method = 'singleshell';
    cfg.siunits = 'yes';
    headmodel = ft_prepare_headmodel(cfg, mesh);

    % Cortex
    ctx = ft_read_headshape(DD{1}.inv{1}.mesh.tess_ctx);
    ctx = ft_convert_units(ctx, 'm');

    % Prepare sourcemodel - use cortical mesh
    cfg = [];
    cfg.method = 'basedoncortex';
    cfg.headshape = ctx;
    cfg.headmodel = headmodel;
    cfg.inwardshift = 0;
    src = ft_prepare_sourcemodel(cfg);
    

    % Read MRI for plotting
    mri_orig = ft_read_mri(DD{1}.inv{1}.mesh.sMRI);
        
    subid = char(meta_data{recording,'sub'});

    for pp = 1:length(DD)
        
        D = DD{pp};

        % Pick out data of interest
        X_standards = D(indchantype(D, 'MEGMAG', 'GOOD'), :, standards);
        X_deviants = D(indchantype(D, 'MEGMAG', 'GOOD'), :, deviants);

        % T-test across trials
        SE_standards = std(X_standards, [], 3)./sqrt(size(X_standards, 3));
        t_standards = mean(X_standards,3)./SE_standards;
        SE_deviants = std(X_deviants, [], 3)./sqrt(size(X_deviants, 3));
        t_deviants = mean(X_deviants,3)./SE_deviants;

        % Unpaired t-test equal variance between deviants and standards for MMN response
        n1 = size(X_deviants,3);
        n2 = size(X_standards,3);
        SE = sqrt(((n1-1)*std(X_deviants, [], 3).^2 + (n2-1)*std(X_standards, [], 3).^2)./(n1 + n2 - 2))*...
            sqrt(1/n1 + 1/n2);
        t_diff = (mean(X_deviants,3) - mean(X_standards,3))./SE;
    
        % Plot t stat
        figure; 
        t = tiledlayout(1,2);
        yl = [];
        
        nexttile(t); hold on; grid on; box on;

        % Standards
        plot(DD{pp}.time*1e3 - delay, t_standards(chans_to_plot(str2double(subid(end))),:), 'LineWidth', 3, 'LineStyle', '--', 'color', [0.3639    0.5755    0.7484]);

        % Deviants
        plot(DD{pp}.time*1e3 - delay, t_deviants(chans_to_plot(str2double(subid(end))),:), 'LineWidth', 3, 'LineStyle', ':', 'color', [0.9153    0.2816    0.2878]);

        % Difference
        plot(DD{pp}.time*1e3 - delay, t_diff(chans_to_plot(str2double(subid(end))),:), 'LineWidth', 3, 'LineStyle', '-', 'color', [0.3373    0.3020    0.2902]);
        
        % Set axes limits
        set(gca, 'FontSize', 18);
        xlim([-100 400]);
        xlabel('Time (ms)');
        ylim(t.Children, [-1 1]*13);
        ylabel(t, 't-stat', 'FontSize', 18);
        lgd = legend('Standards', 'Deviants', 'MMN', 'location', 'eastoutside');
        set(gcf, 'Position', [626   476   821   285]);
        
        nexttile(t); % Purely to keep plot shape the same as in source-space

        if pp == 1
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'no_amm');
        elseif pp == 2
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'hfc');
        elseif pp == 3
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'hfc_with_gradients');
        elseif pp == 4
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'amm_spatial');
        else
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'amm');
        end
        if ~exist(save_loc, 'dir')
            mkdir(save_loc);
        end

        save_name = sprintf('%s_MMN_sensor_level', extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'));
        print(fullfile(save_loc, save_name),'-dpng','-r300');
    
        % Plot topography of MMN at peak between 70 and 150 ms
        lay_name = fullfile(meta_data{recording, "analysed_data_loc"}, ...
            sprintf('%s_2Dlayout.mat', extractBefore(meta_data{recording, "raw_data_name"}, '_meg.lvm')));
        load(lay_name);
    
        data = ftraw(D);
        cfg = [];
        cfg.channel = indchantype(D, 'MEGMAG', 'GOOD');
        tavdata = ft_timelockanalysis(cfg, data);
        tavdata.avg = t_diff;

        % Find layout index for plotted channel
        chanlabels = D.chanlabels(indchantype(D, 'MEGMAG', 'GOOD'));
        chan_to_plot_name = chanlabels{chans_to_plot(str2double(subid(end)))};
        chan_to_plot_name = [chan_to_plot_name(1:end-1), rad_ax];
        chan_to_plot_lay_idx = find(ismember(chanlabels, chan_to_plot_name));

        % Find time period of interest
        tested_time_period = D.time(D.time > (100+delay)*1e-3 & D.time < (275+delay)*1e-3);
        [~, peak_latency] = max(abs(t_diff(chans_to_plot(str2double(subid(end))), D.time > (100+delay)*1e-3 & D.time < (275+delay)*1e-3)));
    
        % Plot
        figure; hold on;
        cfg = [];
        cfg.layout    = lay;
        cfg.colorbar  = 'EastOutside';
        cfg.colorbartext = 't-stat (MMN peak)';
        cfg.zlim = [-1 1]*5;
        cfg.colormap  = colormap123;
        cfg.xlim = tested_time_period(peak_latency)*[1, 1];
        cfg.comment = 'no';
        cfg.highlight   = 'on';
        cfg.highlightchannel = chan_to_plot_lay_idx;
        cfg.interactive = 'no';
        cfg.highlightsymbol = 'd';
        cfg.highlightcolor = 'k';
        cfg.highlightsize = 10;
        cfg.markersymbol = 'o';
        cfg.figure = gca;
        set(gca, 'FontSize', 24);
        ft_topoplotER(cfg, tavdata)

        set(gcf, 'Position', [994   704   404   274]);
        annotation('textbox', [0, 0, 0.6, 0.15], ...
            'string', sprintf('Peak latency: %.f ms', tested_time_period(peak_latency)*1e3 - delay), 'FontSize', 16, 'EdgeColor', 'None');

        save_name = sprintf('%s_MMN_topography', extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'));
        print(fullfile(save_loc, save_name),'-dpng','-r300');

        % Dipole fit on average MMN

        % Prepare leadfields
        sens = DD{pp}.sensors('MEG');
        sens = ft_convert_units(sens, 'm');
        cfg                  = [];
        cfg.grad             = sens;
        cfg.headmodel        = headmodel;
        cfg.reducerank       = 2;
        cfg.channel          = DD{pp}.chanlabels(indchantype(DD{pp}, 'MEGMAG', 'GOOD'));
        cfg.sourcemodel      = src;
        sourcemodel = ft_prepare_leadfield(cfg);

        % Format data for fieldtrip
        data = ftraw(DD{pp});

        % Average
        cfg = [];
        cfg.channel = indchantype(D, 'MEGMAG', 'GOOD');
        tl_data = ft_timelockanalysis(cfg, data);
        tl_data.avg = (mean(X_deviants,3) - mean(X_standards,3));

        % Dipole fit
        cfg = [];
        cfg.latency = tested_time_period(peak_latency)*[1, 1]+[-0.025 0.025];
        cfg.numdipoles = 2;
        cfg.symmetry = [];
        cfg.gridsearch = 'no';
        cfg.dip.pos = aud_nat;
        cfg.headmodel = headmodel;
        cfg.sourcemodel = sourcemodel;
        cfg.channel = DD{pp}.chanlabels(indchantype(DD{pp}, 'MEGMAG', 'GOOD'));
        cfg.senstype = 'meg';
        source = ft_dipolefitting(cfg, tl_data);
        source.dip = ft_convert_units(source.dip, 'mm');

        % Plot dipole position
        % Axial
        pos = mean(source.dip.pos,1);
        figure; hold on;
        ft_plot_dipole([source.dip.pos(1,[1,2]), 0.1], mean(source.dip.mom(1:3,:),2), 'color', '#E07BE0', 'unit', 'mm'); % Left
        ft_plot_dipole([source.dip.pos(2,[1,2]), 0.1], mean(source.dip.mom(4:6,:),2), 'color', '#45C9B7', 'unit', 'mm'); % Right
        ft_plot_slice(mri_orig.anatomy, 'transform', mri_orig.transform, 'location', pos, 'orientation', [0 0 1], 'resolution', 0.1);
        axis tight
        axis off
        view(0,90);
        save_name = sprintf('%s_MMN_ft_dip_fit_axial', extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'));
        print(fullfile(save_loc, save_name),'-dpng','-r300');

        % Coronal
        figure; hold on;
        ft_plot_dipole([source.dip.pos(1,1), -150, source.dip.pos(1,3)], mean(source.dip.mom(1:3,:),2), 'color', '#E07BE0', 'unit', 'mm'); % Left
        ft_plot_dipole([source.dip.pos(2,1), -150, source.dip.pos(2,3)], mean(source.dip.mom(4:6,:),2), 'color', '#45C9B7', 'unit', 'mm'); % Right
        ft_plot_slice(mri_orig.anatomy, 'transform', mri_orig.transform, 'location', pos, 'orientation', [0 1 0], 'resolution', 0.1);
        view(0,0);
        axis tight
        axis off
        save_name = sprintf('%s_MMN_ft_dip_fit_coronal', extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'));
        print(fullfile(save_loc, save_name),'-dpng','-r300');
        

        % Plot estimated source current
        % Get lead field along mean dipole ori
        sens = data.grad;
        sens = ft_convert_units(sens, 'm');
        [headmodel, sens] = ft_prepare_vol_sens(headmodel, sens, 'channel', DD{pp}.chanlabels(indchantype(DD{pp}, 'MEGMAG', 'GOOD')));
        Gxyz = ft_compute_leadfield(source.dip.pos*1e-3, sens, headmodel, 'dipoleunit', 'nA*m', 'chanunit', repmat({'fT'}, size(sens.label,1),1));
        L = zeros(size(Gxyz,1),2);
        for ind = 1:2
            dip_ori = mean(source.dip.mom((ind-1)*3+1:ind*3,:),2);
            dip_ori = dip_ori./norm(dip_ori);
            L(:, ind) = Gxyz(:, (3*ind-2):(3*ind))*dip_ori;
        end

        X_standards = zeros(size(L,2), size(DD{pp},2), length(standards));
        X_deviants = zeros(size(L,2), size(DD{pp},2), length(deviants));

        for tt = 1:length(deviants)
            X_deviants(:,:,tt) = pinv(L)*DD{pp}(indchannel(DD{pp},sens.label),:,deviants(tt));
        end
        for tt = 1:length(standards)
            X_standards(:,:,tt) = pinv(L)*DD{pp}(indchannel(DD{pp},sens.label),:,standards(tt));
        end

        % T-test across trials
        SE_standards = std(X_standards, [], 3)./sqrt(size(X_standards, 3));
        t_standards = mean(X_standards,3)./SE_standards;
        SE_deviants = std(X_deviants, [], 3)./sqrt(size(X_deviants, 3));
        t_deviants = mean(X_deviants,3)./SE_deviants;
        
        % Unpaired t-test equal variance between deviants and standards for MMN response
        n1 = size(X_deviants,3);
        n2 = size(X_standards,3);
        SE = sqrt(((n1-1)*std(X_deviants, [], 3).^2 + (n2-1)*std(X_standards, [], 3).^2)./(n1 + n2 - 2))*...
            sqrt(1/n1 + 1/n2);
        t_diff = (mean(X_deviants,3) - mean(X_standards,3))./SE;


        % Plot
        figure; 
        t = tiledlayout(1,2);
        yl = [-1 1]*10;
        for ind = 1:2
            nexttile(t); hold on; grid on; box on;

            % Time period where dipole was fit
            fill([min(cfg.latency), max(cfg.latency), max(cfg.latency), min(cfg.latency)]*1e3 - delay, ...
                [yl(1), yl(1), yl(2), yl(2)], [231, 196, 170]./255, 'EdgeColor', 'None', 'FaceAlpha', 0.3);

            % Standards
            % plot(DD{pp}.time*1e3 - delay, mean(X_standards(ind,:,:), 3), 'LineWidth', 3, 'LineStyle', '--', 'color', [0.3639    0.5755    0.7484]);
            plot(DD{pp}.time*1e3 - delay, t_standards(ind,:), 'LineWidth', 3, 'LineStyle', '--', 'color', [0.3639    0.5755    0.7484]);

            % Deviants
            % plot(DD{pp}.time*1e3 - delay, mean(X_deviants(ind,:,:), 3), 'LineWidth', 3, 'LineStyle', ':', 'color', [0.9153    0.2816    0.2878]);
            plot(DD{pp}.time*1e3 - delay, t_deviants(ind,:), 'LineWidth', 3, 'LineStyle', ':', 'color', [0.9153    0.2816    0.2878]);

            % Difference
            % plot(DD{pp}.time*1e3 - delay, mean(X_deviants(ind,:,:),3) - mean(X_standards(ind,:,:), 3), 'LineWidth', 3, 'LineStyle', '-', 'color', [0.3373    0.3020    0.2902]);
            plot(DD{pp}.time*1e3 - delay, t_diff(ind,:), 'LineWidth', 3, 'LineStyle', '-', 'color', [0.3373    0.3020    0.2902]);
            set(gca, 'FontSize', 18);
            xlim([-100 400]);
            % yl(ind) = max(abs(ylim));
            xlabel('Time (ms)');

            if ind == 1
                title('Left Hemisphere');
            else
                title('Right Hemisphere');
            end

        end

        % Set axes limits and legend
        % ylim(t.Children, [-1 1]*max(abs(yl)));
        ylim(t.Children, [-1 1]*10);
        % ylabel(t, {'Estimated Source', 'Current (nAm)'}, 'FontSize', 18);
        ylabel(t, 't-stat', 'FontSize', 18);
        lgd = legend('Standards', 'Deviants', 'MMN', 'location', 'eastoutside');
        set(gcf, 'Position', [626   476   821   285]);

        % Add text to indicate how many trials per condition
        
        annotation('textbox', [lgd.Position(1), lgd.Position(2) - 0.35, lgd.Position(3), 0.3], ...
            'string', sprintf('# deviants: %.f\n# standards: %.f', length(deviants), length(standards)), 'FontSize', 16, 'EdgeColor', 'None');

        save_name = sprintf('%s_dipfit_on_MMN_trace_dipole_all_sets', extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'));
        print(fullfile(save_loc, save_name),'-dpng','-r300');
        
    end

    close all
end

%% Plot power spectral densities

for closed_loop = [true, false]
    for walking = [true, false]

        if closed_loop
            if walking
                recording_order = {"sub-001_task-walkingClosed_meg.lvm", ...
                    "sub-002_task-walkingClosed_run-001_meg.lvm", ...
                    "sub-002_task-walkingClosed_run-002_meg.lvm", "sub-003_task-walkingClosed_meg.lvm", };
            else
                recording_order = {"sub-001_task-seatedClosed_meg.lvm", ...
                    "sub-002_task-seatedClosed_meg.lvm", "sub-003_task-seatedClosed_meg.lvm"};
            end
        else
            if walking
                recording_order = {"sub-001_task-walkingOpen_meg.lvm", ...
                    "sub-002_task-walkingOpen_run-001_meg.lvm", ...
                    "sub-002_task-walkingOpen_run-002_meg.lvm", "sub-003_task-walkingOpen_meg.lvm"};
            else
                recording_order = {"sub-001_task-seatedOpen_meg.lvm", ...
                    "sub-002_task-seatedOpen_meg.lvm", "sub-003_task-seatedOpen_meg.lvm"};
            end
        end
        
        % Plot PSD
        figure; hold on; grid on; box on;
        co = colororder(gca);
        line_style = {'-', '--', ':', '-.'};
        DD = [];
        pl = [];
        ii = 0;
        
        start_string = {'ffft_', 'hffft_', 'm2ffft_', 'mffft_'};
        ref_freq = 0:0.1:187.5;
        
        for pp = 1:length(start_string)
            pof = [];
            for recording = 1:length(recording_order)
                rec_idx = find(contains(meta_data.raw_data_name, recording_order{recording}));

                DD = spm_eeg_load(char(fullfile(meta_data{rec_idx, "analysed_data_loc"}, ...
                    strcat(start_string{pp}, extractBefore(recording_order{recording}, '.lvm'), '.mat'))));
            
                S = [];
                S.D = DD;
                S.channels = DD.chanlabels(indchantype(DD, 'MEGMAG', 'GOOD'));
                S.triallength = 10e3;
                [po, freq] = spm_opm_psd(S);
        
                F = griddedInterpolant(freq', po, 'linear');
                po = F(ref_freq);
        
                pof = [pof, mean(po, 2)];
            end
        
            ii = ii+1;
        
            mp = mean(pof, 2);
            sem = std(pof,[],2)./sqrt(size(pof,2));
        
            fill([ref_freq'; flipud(ref_freq')], [mp-sem; flipud(mp+sem)], co(ii,:),...
                'linestyle', 'none', 'FaceAlpha', 0.4)
            pl(end+1) = plot(ref_freq, mp, 'LineWidth', 2, 'LineStyle', line_style{ii}, 'color', co(ii,:));
        
        end
        set(gca,'yscale','log');
        set(gca, 'FontSize', 18);
        xlim([2 40]);
        ylim([10 1e4]);
        legend(pl, {'No spatial filter', 'HFC', 'AMM spatial', 'AMM with temporal'});
        xlabel('Frequency (Hz)');
        ylabel('PSD ($$fT\sqrt[-1]{Hz}$$)','interpreter','latex');
        
        if closed_loop
            if walking
                print(fullfile(extractBefore(meta_data{rec_idx,'results_save_loc'}, '\sub-'), 'PSD_closed_walking'),'-dpng','-r300');
            else
                print(fullfile(extractBefore(meta_data{rec_idx,'results_save_loc'}, '\sub-'), 'PSD_closed_seated'),'-dpng','-r300');
            end
        else
            if walking
                print(fullfile(extractBefore(meta_data{rec_idx,'results_save_loc'}, '\sub-'), 'PSD_open_walking'),'-dpng','-r300');
            else
                print(fullfile(extractBefore(meta_data{rec_idx,'results_save_loc'}, '\sub-'), 'PSD_open_seated'),'-dpng','-r300');
            end
        end

        % Plot shielding factor
        figure; hold on; grid on; box on;
    
        ii = 1;
        pl = [];
        for pp = 2:length(start_string)
            shf = [];
            for recording = 1:length(recording_order)

                rec_idx = find(contains(meta_data.raw_data_name, recording_order{recording}));

                S = [];
                S.D1 = spm_eeg_load(char(fullfile(meta_data{rec_idx, "analysed_data_loc"}, ...
                    strcat(start_string{1}, extractBefore(meta_data{rec_idx, "raw_data_name"}, '.lvm'), '.mat'))));
                S.D2 = spm_eeg_load(char(fullfile(meta_data{rec_idx, "analysed_data_loc"}, ...
                    strcat(start_string{pp}, extractBefore(meta_data{rec_idx, "raw_data_name"}, '.lvm'), '.mat'))));
                S.channels = S.D1.chanlabels(indchantype(S.D1, 'MEGMAG', 'GOOD'));
                S.plot = 0;
                S.triallength = 10e3;
                S.dB = 1;
                [shield, freq] = spm_opm_rpsd(S);

                F = griddedInterpolant(freq', shield, 'linear');
                shield = F(ref_freq);
        
                shf = [shf, mean(shield, 2)];
                
            end
            ii = ii + 1;

            mp = mean(shf, 2);
            sem = std(shf,[],2)./sqrt(size(shf,2));

            fill([ref_freq'; flipud(ref_freq')], [mp-sem; flipud(mp+sem)], co(ii,:),...
                'linestyle', 'none', 'FaceAlpha', 0.4)
            pl(end+1) = plot(ref_freq, mp, 'LineWidth', 2, 'color', co(ii,:), 'LineStyle', line_style{ii});
        end
        set(gca, 'FontSize', 18);
        xlim([2 40]);
        ylim([0 30]);
        legend(pl, {'HFC', 'AMM spatial', 'AMM with temporal'});
        xlabel('Frequency (Hz)');
        ylabel('Interference Reduction (dB)');

        if closed_loop
            if walking
                print(fullfile(extractBefore(meta_data{rec_idx,'results_save_loc'}, '\sub-'), 'ShieldingFactor_closed_walking'),'-dpng','-r300');
            else
                print(fullfile(extractBefore(meta_data{rec_idx,'results_save_loc'}, '\sub-'), 'ShieldingFactor_closed_seated'),'-dpng','-r300');
            end
        else
            if walking
                print(fullfile(extractBefore(meta_data{rec_idx,'results_save_loc'}, '\sub-'), 'ShieldingFactor_open_walking'),'-dpng','-r300');
            else
                print(fullfile(extractBefore(meta_data{rec_idx,'results_save_loc'}, '\sub-'), 'ShieldingFactor_open_seated'),'-dpng','-r300');
            end
        end

    end
end


%% Plot power spectral densities - HFC with linear gradients

for closed_loop = [true, false]
    for walking = [true, false]

        if closed_loop
            if walking
                recording_order = {"sub-001_task-walkingClosed_meg.lvm", ...
                    "sub-002_task-walkingClosed_run-001_meg.lvm", ...
                    "sub-002_task-walkingClosed_run-002_meg.lvm", "sub-003_task-walkingClosed_meg.lvm", };
            else
                recording_order = {"sub-001_task-seatedClosed_meg.lvm", ...
                    "sub-002_task-seatedClosed_meg.lvm", "sub-003_task-seatedClosed_meg.lvm"};
            end
        else
            if walking
                recording_order = {"sub-001_task-walkingOpen_meg.lvm", ...
                    "sub-002_task-walkingOpen_run-001_meg.lvm", ...
                    "sub-002_task-walkingOpen_run-002_meg.lvm", "sub-003_task-walkingOpen_meg.lvm"};
            else
                recording_order = {"sub-001_task-seatedOpen_meg.lvm", ...
                    "sub-002_task-seatedOpen_meg.lvm", "sub-003_task-seatedOpen_meg.lvm"};
            end
        end
        
        % Plot PSD
        figure; hold on; grid on; box on;
        co = colororder(gca);
        line_style = {'-', '--', ':', '-.', '-'};
        DD = [];
        pl = [];
        ii = 0;
        
        start_string = {'ffft_', 'hffft_', 'h2ffft_', 'm2ffft_', 'mffft_'};
        ref_freq = 0:0.1:187.5;
        
        for pp = [setdiff(1:length(start_string), 3), 3]
            pof = [];
            for recording = 1:length(recording_order)
                rec_idx = find(contains(meta_data.raw_data_name, recording_order{recording}));

                DD = spm_eeg_load(char(fullfile(meta_data{rec_idx, "analysed_data_loc"}, ...
                    strcat(start_string{pp}, extractBefore(recording_order{recording}, '.lvm'), '.mat'))));
            
                S = [];
                S.D = DD;
                S.channels = DD.chanlabels(indchantype(DD, 'MEGMAG', 'GOOD'));
                S.triallength = 10e3;
                [po, freq] = spm_opm_psd(S);
        
                F = griddedInterpolant(freq', po, 'linear');
                po = F(ref_freq);
        
                pof = [pof, mean(po, 2)];
            end
        
            ii = ii+1;
        
            mp = mean(pof, 2);
            sem = std(pof,[],2)./sqrt(size(pof,2));
        
            fill([ref_freq'; flipud(ref_freq')], [mp-sem; flipud(mp+sem)], co(ii,:),...
                'linestyle', 'none', 'FaceAlpha', 0.4)
            pl(end+1) = plot(ref_freq, mp, 'LineWidth', 2, 'LineStyle', line_style{ii}, 'color', co(ii,:));
        
        end
        set(gca,'yscale','log');
        set(gca, 'FontSize', 18);
        xlim([2 40]);
        ylim([10 1e4]);
        legend(pl, {'No spatial filter', 'HFC', 'AMM spatial', 'AMM with temporal', 'HFC with gradients'});
        xlabel('Frequency (Hz)');
        ylabel('PSD ($$fT\sqrt[-1]{Hz}$$)','interpreter','latex');
        
        if closed_loop
            if walking
                print(fullfile(extractBefore(meta_data{rec_idx,'results_save_loc'}, '\sub-'), 'PSD_closed_walking_hfcgrad'),'-dpng','-r300');
            else
                print(fullfile(extractBefore(meta_data{rec_idx,'results_save_loc'}, '\sub-'), 'PSD_closed_seated_hfcgrad'),'-dpng','-r300');
            end
        else
            if walking
                print(fullfile(extractBefore(meta_data{rec_idx,'results_save_loc'}, '\sub-'), 'PSD_open_walking_hfcgrad'),'-dpng','-r300');
            else
                print(fullfile(extractBefore(meta_data{rec_idx,'results_save_loc'}, '\sub-'), 'PSD_open_seated_hfcgrad'),'-dpng','-r300');
            end
        end

        % Plot shielding factor
        figure; hold on; grid on; box on;
    
        ii = 1;
        pl = [];
        for pp = [setdiff(2:length(start_string), 3), 3]
            shf = [];
            for recording = 1:length(recording_order)

                rec_idx = find(contains(meta_data.raw_data_name, recording_order{recording}));

                S = [];
                S.D1 = spm_eeg_load(char(fullfile(meta_data{rec_idx, "analysed_data_loc"}, ...
                    strcat(start_string{1}, extractBefore(meta_data{rec_idx, "raw_data_name"}, '.lvm'), '.mat'))));
                S.D2 = spm_eeg_load(char(fullfile(meta_data{rec_idx, "analysed_data_loc"}, ...
                    strcat(start_string{pp}, extractBefore(meta_data{rec_idx, "raw_data_name"}, '.lvm'), '.mat'))));
                S.channels = S.D1.chanlabels(indchantype(S.D1, 'MEGMAG', 'GOOD'));
                S.plot = 0;
                S.triallength = 10e3;
                S.dB = 1;
                [shield, freq] = spm_opm_rpsd(S);

                F = griddedInterpolant(freq', shield, 'linear');
                shield = F(ref_freq);
        
                shf = [shf, mean(shield, 2)];
                
            end
            ii = ii + 1;

            mp = mean(shf, 2);
            sem = std(shf,[],2)./sqrt(size(shf,2));

            fill([ref_freq'; flipud(ref_freq')], [mp-sem; flipud(mp+sem)], co(ii,:),...
                'linestyle', 'none', 'FaceAlpha', 0.4)
            pl(end+1) = plot(ref_freq, mp, 'LineWidth', 2, 'color', co(ii,:), 'LineStyle', line_style{ii});
        end
        set(gca, 'FontSize', 18);
        xlim([2 40]);
        ylim([0 30]);
        legend(pl, {'HFC', 'AMM spatial', 'AMM with temporal', 'HFC with gradients'});
        xlabel('Frequency (Hz)');
        ylabel('Interference Reduction (dB)');

        if closed_loop
            if walking
                print(fullfile(extractBefore(meta_data{rec_idx,'results_save_loc'}, '\sub-'), 'ShieldingFactor_closed_walking_hfcgrad'),'-dpng','-r300');
            else
                print(fullfile(extractBefore(meta_data{rec_idx,'results_save_loc'}, '\sub-'), 'ShieldingFactor_closed_seated_hfcgrad'),'-dpng','-r300');
            end
        else
            if walking
                print(fullfile(extractBefore(meta_data{rec_idx,'results_save_loc'}, '\sub-'), 'ShieldingFactor_open_walking_hfcgrad'),'-dpng','-r300');
            else
                print(fullfile(extractBefore(meta_data{rec_idx,'results_save_loc'}, '\sub-'), 'ShieldingFactor_open_seated_hfcgrad'),'-dpng','-r300');
            end
        end

    end
end

%% Plot sensor positions and mark new sensor cabling

ctx = gifti(fullfile(spm('dir'), '\canonical\scalp_2562.surf.gii'));
opm = [];
opm.verts = [6.25 -7.25 6.2; 
    -10.35 -7.25 6.2; 
    -10.35 5.15 6.2; 
    6.25 5.15 6.2; 
    6.25 -7.25 -20.2; 
    -10.35 -7.25 -20.2; 
    -10.35 5.15 -20.2; 
    6.25 5.15 -20.2];
opm.faces = [1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];

opm.face_midpoints = zeros(size(opm.faces,1), 3);
for ii = 1:size(opm.faces,1)
    opm.face_midpoints(ii,:) = mean(opm.verts(opm.faces(ii,:), :), 1);
end

cp = [240, 80, 57; 61, 101, 165; 168, 182, 204]./255;
surface_marker = {"square", "^", "o"};

subIDs = unique(meta_data{:,'sub'});
for sub = 1:length(subIDs)
    % Get first recording
    rec_inds = find(contains(meta_data{:,'sub'}, subIDs(sub)));
    D = spm_eeg_load(char(fullfile(meta_data{rec_inds(1), "analysed_data_loc"}, ...
            strcat(extractBefore(meta_data{rec_inds(1), "raw_data_name"}, '.lvm'), '.mat'))));

    % Get sensor positions and orientations
    chanpos = D.sensors('MEG').chanpos;
    chanori = D.sensors('MEG').chanori;
    label = D.sensors('MEG').label;
    label = extractAfter(label, '-');

    % Go from channels to sensors
    x_ori = -chanori(endsWith(label, '-X'),:);
    y_ori = chanori(endsWith(label, '-Y'),:);
    z_ori = chanori(endsWith(label, '-Z'),:);
    chanpos = chanpos(endsWith(label, '-X'),:);
    label = extractBefore(label(endsWith(label, '-X'),:), '-');

    % Get which have the new cables
    if strcmp(subIDs{sub}, 'sub-003') || strcmp(subIDs{sub}, 'sub-004')
        new_cables = cellfun(@(x)strcmp(x(2), 'B'), label);
    else
        new_cables = ones(length(label),1);
    end

    % Get badchannels
    badchans = D.chanlabels(badchannels(D));
    badchans = unique(extractBefore(extractAfter(badchans, '-'), '-'));
    [~, badchans] = intersect(label, badchans);

    % Write out number of sensors, channels and badchannels for reporting
    fprintf('%s : %.f sensors, %.f channels, %.f bad channels, %.f new cables\n',...
        subIDs{sub}, length(label), length(indchantype(D, 'MEGMAG')), length(badchannels(D)), sum(new_cables))

    % Change sensor positions to MNI space
    M = D.inv{1}.datareg.toMNI;
    chanpos = M*cat(1, chanpos', ones(1, size(chanpos,1)));
    chanpos = chanpos(1:3,:)';
    x_ori = M(1:3,1:3)*x_ori';
    x_ori = x_ori';
    x_ori = x_ori./repmat(sqrt(sum(x_ori.^2, 2)), 1, 3);
    y_ori = M(1:3,1:3)*y_ori';
    y_ori = y_ori';
    y_ori = y_ori./repmat(sqrt(sum(y_ori.^2, 2)), 1, 3);
    z_ori = M(1:3,1:3)*z_ori';
    z_ori = z_ori';
    z_ori = z_ori./repmat(sqrt(sum(z_ori.^2, 2)), 1, 3);

    % Plot
    figure; hold on;
    patch('Faces', ctx.faces, 'Vertices', ctx.vertices, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'None');
    daspect([1 1 1]);
    % Warp opm vertices to get patch at each sensor position
    lbad_sens = [];
    lold_wire = [];
    lnew_wire = [];
    for sens = 1:length(chanpos)
        opm_warped_verts = repmat(chanpos(sens,:), size(opm.verts,1), 1) + ...
            repmat(x_ori(sens,:), size(opm.verts,1), 1).*repmat(opm.verts(:,1), 1, size(x_ori,2)) + ...
            repmat(y_ori(sens,:), size(opm.verts,1), 1).*repmat(opm.verts(:,2), 1, size(y_ori,2)) + ...
            repmat(z_ori(sens,:), size(opm.verts,1), 1).*repmat(opm.verts(:,3), 1, size(z_ori,2));
        opm_center = mean(opm_warped_verts, 1);
        if any(ismember(badchans, sens))
            if isempty(lbad_sens)
                lbad_sens = scatter3(opm_center(1), opm_center(2), opm_center(3), 35, cp(1,:),...
                    surface_marker{1}, 'filled', 'MarkerEdgeColor', 'k');
            else
                scatter3(opm_center(1), opm_center(2), opm_center(3), 35, cp(1,:), ...
                    surface_marker{1}, 'filled', 'MarkerEdgeColor', 'k');
            end
            patch('Faces', opm.faces, 'Vertices', opm_warped_verts, 'FaceColor', cp(1,:), 'FaceAlpha', 0.6);

        elseif new_cables(sens)
            if isempty(lnew_wire)
                lnew_wire = scatter3(opm_center(1), opm_center(2), opm_center(3), 35, cp(2,:),...
                    surface_marker{2}, 'filled', 'MarkerEdgeColor', 'k');
            else
                scatter3(opm_center(1), opm_center(2), opm_center(3), 35, cp(2,:),...
                    surface_marker{2}, 'filled', 'MarkerEdgeColor', 'k');
            end
            patch('Faces', opm.faces, 'Vertices', opm_warped_verts, 'FaceColor', cp(2,:), 'FaceAlpha', 0.6);
        else
            if isempty(lold_wire)
                lold_wire = scatter3(opm_center(1), opm_center(2), opm_center(3), 35, cp(3,:),...
                    surface_marker{3}, 'filled', 'MarkerEdgeColor', 'k');
            else
                scatter3(opm_center(1), opm_center(2), opm_center(3), 35, cp(3,:),...
                    surface_marker{3}, 'filled', 'MarkerEdgeColor', 'k');
            end
            patch('Faces', opm.faces, 'Vertices', opm_warped_verts, 'FaceColor', cp(3,:), 'FaceAlpha', 0.6);
        end
    end
    set(gca, 'Visible', 'off');

    % Save
    view(-90,0);
    print(fullfile(meta_data{rec_inds(1), "results_save_loc"}, "helmet_left_view"),'-dpng','-r300');
    view(0,0);
    print(fullfile(meta_data{rec_inds(1), "results_save_loc"}, "helmet_back_view"),'-dpng','-r300');
    view(90,0);
    print(fullfile(meta_data{rec_inds(1), "results_save_loc"}, "helmet_right_view"),'-dpng','-r300');
    view(180,0);

    if strcmp(subIDs{sub}, 'sub-003')
        lgd = legend([lnew_wire, lbad_sens], 'Sensor with new wire', 'Bad channel', 'FontSize', 16);
        pos = get(lgd, 'Position');
        set(lgd, 'Position', [0.,0.82,pos(3),pos(4)])
    else
        lgd = legend([lnew_wire, lold_wire, lbad_sens], 'Sensor with new wire', 'Sensor with old wire', 'Bad channel', 'FontSize', 16);
        pos = get(lgd, 'Position');
        set(lgd, 'Position', [0.,0.79,pos(3),pos(4)])
    end
    
    print(fullfile(meta_data{rec_inds(1), "results_save_loc"}, "helmet_front_view"),'-dpng','-r300');
end
    
%% ROI analysis, dipole

DD = {};

for recording = 1:size(meta_data,1)

    start_string = {'e_ffft_', 'e_hffft_', 'e_h2ffft_', 'e_m2ffft_', 'e_mffft_'};
    for pp = 1:length(start_string)
        DD{pp} = spm_eeg_load(char(fullfile(meta_data{recording, "analysed_data_loc"}, ...
            strcat(start_string{pp}, extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'), '.mat'))));
    end

    M = gifti(fullfile(spm('dir'), 'canonical', 'cortex_5124.surf.gii'));
    pos = [-54 -14 8; 54 -14 8]; % Auditory cortices
    [~, leftind] = min(sqrt(sum((M.vertices - pos(1,:)).^2,2)));
    [~, rightind] = min(sqrt(sum((M.vertices - pos(2,:)).^2,2)));
    L = full(spm_eeg_lgainmat(DD{1},[leftind, rightind]));

    for pp = 1:length(DD)
        good_trials = indtrial(DD{pp}, 'tone', 'GOOD');
        X = zeros(size(L,2), size(DD{pp},2), length(good_trials));
        for tt = 1:length(good_trials)
            X(:,:,tt) = pinv(L)*DD{pp}(indchantype(DD{pp},'MEGMAG','GOOD'),:,...
                good_trials(tt));
        end

        figure; hold on; grid on; box on;
        se = std(X, [], 3)./sqrt(size(X,3));
        plot(DD{pp}.time*1e3, mean(X, 3)./se, 'LineWidth', 2);

        % Get significant t value
        df = length(good_trials)-1;
        alpha = 0.025/(size(X, 2)*40/DD{pp}.fsample*size(X,1));
        sigt = tinv(1-alpha, df);
        plot([min(DD{pp}.time), max(DD{pp}.time)]*1e3, sigt*[1 1], 'k--', 'LineWidth', 2);
        plot([min(DD{pp}.time), max(DD{pp}.time)]*1e3, -sigt*[1 1], 'k--', 'LineWidth', 2);

        xlim([-100 400]);
        ylim([-20 20])
        xlabel('Time (ms)');
        ylabel('t-stat');
        set(gcf, 'Position', [680   654   451   344]);
        set(gca, 'FontSize', 24);

        if pp == 1
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'no_amm');
        elseif pp == 2
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'hfc');
        elseif pp == 3
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'hfc_with_gradients');
        elseif pp == 4
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'amm_spatial');
        else
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'amm');
        end

        save_name = sprintf('%s_ROI_t_val_dipole', extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'));
        print(fullfile(save_loc, save_name),'-dpng','-r300');
    end
end

%% Source localisation, Minimum Norm

% Inflate mesh for displaying source space
M_original = gifti(fullfile(spm('dir'), 'canonical', 'cortex_5124.surf.gii'));
M = spm_mesh_inflate(spm_mesh_inflate(M_original));
mesh = ft_read_headshape(fullfile(spm('dir'), 'canonical', 'cortex_5124.surf.gii'));
mesh.pos = M.vertices;
clear M

% Get smoothing kernel
[~,Di] = spm_mesh_neighbours(M_original,1);
muNeighbour = mean(mean(Di));
n = round((8/muNeighbour)^2);

for recording = 1:size(meta_data,1)
    cd(meta_data{recording, "analysed_data_loc"});
    start_string = {'e_ffft_', 'e_hffft_', 'e_h2ffft_', 'e_m2ffft_', 'e_mffft_'};
    for pp = 1:length(start_string)
        DD{pp} = spm_eeg_load(char(fullfile(meta_data{recording, "analysed_data_loc"}, ...
            strcat(start_string{pp}, extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'), '.mat'))));
    end

    for pp = 1:length(DD)
        matlabbatch = [];
        matlabbatch{1}.spm.meeg.source.invert.D = {DD{pp}.fname};
        matlabbatch{1}.spm.meeg.source.invert.val = 1;
        matlabbatch{1}.spm.meeg.source.invert.whatconditions.all = 1;
        matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.invtype = 'IID';
        matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.woi = [-Inf Inf];
        matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.foi = [2 40];
        matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.hanning = 1;
        matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.priors.priorsmask = {''};
        matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.priors.space = 1;
        matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.restrict.locs = zeros(0, 3);
        matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.restrict.radius = 32;
        matlabbatch{1}.spm.meeg.source.invert.isstandard.custom.restrict.mask = {''};
        matlabbatch{1}.spm.meeg.source.invert.modality = {'All'};
        matlabbatch{2}.spm.meeg.source.results.D(1) = cfg_dep('Source inversion: M/EEG dataset(s) after imaging source reconstruction', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','D'));
        matlabbatch{2}.spm.meeg.source.results.val = 1;
        matlabbatch{2}.spm.meeg.source.results.woi = [0 200];
        matlabbatch{2}.spm.meeg.source.results.foi = [2 40];
        matlabbatch{2}.spm.meeg.source.results.ctype = 'evoked';
        matlabbatch{2}.spm.meeg.source.results.space = 0;
        matlabbatch{2}.spm.meeg.source.results.format = 'mesh';
        matlabbatch{2}.spm.meeg.source.results.smoothing = 8;
            
        a = spm_jobman('run',matlabbatch);
        DD{pp} = spm_eeg_load(fullfile(DD{pp}.path, DD{pp}.fname));
        
        % Display
        spm_mesh_render('Disp', a{2}.files{:});
        spm_mesh_render('clim', gca, [0,8.5]);

        H = spm_mesh_render('View',gca, [-90,10]);
        spm_mesh_inflate(H.patch,Inf,1);
        set(gcf, 'color', 'w');

        % Create source time courses from Minimum Norm results
            
        % Get weights
        U = DD{pp}.inv{1}.inverse.U{1};
        weights = DD{pp}.inv{1}.inverse.M;
        good_trials = indtrial(DD{pp}, 'tone', 'GOOD');
        
        % Create each time course
        chaninds = selectchannels(DD{pp}, DD{pp}.inv{1}.forward.channels);
        [~, tinds] = min(abs(DD{pp}.time - 100*1e-3));
        dat = DD{pp}(chaninds,tinds,:);
        nt = length(good_trials);
        vedata = cell(1, nt);
        
        for trial=good_trials
            data = dat(:,:,trial);
            vedata{trial} = transpose((U*data)'*weights');
        end
        
        % Rearrange vedata
        vedata_tp = cell(1, size(dat,2)); % One for each time point
        for tp = 1:length(tinds)
            vedata_tp{tp} = zeros(size(weights,1), nt);
            for trial = good_trials
                vedata_tp{tp}(:,trial) = vedata{trial}(:,tp);
            end
        end
        
        % T-stat over source time courses
        t = cell2mat(cellfun(@(x) mean(x,2)./(std(x,[],2)/sqrt(size(x,2))), vedata_tp, 'UniformOutput', false));
    
        % Smooth image
        t = spm_mesh_smooth(M_original, t, n);
            
        % Find significant values
        df = nt-1;
        alpha = 0.025/(size(DD{pp}.inv{1}.inverse.M,1)*length(tinds));
        sigt = tinv([alpha 1-alpha], df);
        
        % Display spatial maps
        cfg = [];
        cfg.facecolor = [0.4 0.4 0.4];
        cfg.vertexcolor = 'none';
    
        % Plot for each time point of interest
        figure;
        mesh.pow = abs(t);
        mesh.mask = mesh.pow > sigt(2);

        if max(mesh.pow) > sigt
            cfg.method         = 'surface';
            cfg.funcolorlim    = [0 18];
            cfg.funparameter   = 'pow';
            cfg.maskparameter  = 'mask';
            cfg.funcolormap    = 'hot';
            cfg.colorbartext = 't-stat';
            ft_sourceplot(cfg, mesh);
        else
            surf.pos = mesh.pos;
            surf.tri = mesh.tri;
            ft_plot_mesh(surf,'edgecolor', 'none', 'facecolor', cfg.facecolor, 'vertexcolor', cfg.vertexcolor);
            lighting gouraud
            camlight
        end
        view ([-90 0])             % rotate the object in the view
        camlight('headlight')
        set(gcf, 'color', 'w');
        set(gca, 'FontSize', 26);
        material dull

        if pp == 1
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'no_amm');
        elseif pp == 2
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'hfc');
        elseif pp == 3
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'hfc_with_gradients');
        elseif pp == 4
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'amm_spatial');
        else
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'amm');
        end

        save_name = sprintf('%s_min_norm_t_val_left', extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'));
        print(fullfile(save_loc, save_name),'-dpng','-r300');

        view ([90 0])             % rotate the object in the view
        camlight('headlight')

        save_name = sprintf('%s_min_norm_t_val_right', extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'));
        print(fullfile(save_loc, save_name),'-dpng','-r300');
    end
end


%% Get head position and orientation in room space

for recording = 1:size(meta_data,1)

    cd(meta_data{recording, "analysed_data_loc"});
    D = spm_eeg_load(['t_', char(meta_data{recording,"raw_data_name"})]);

    if strcmp(meta_data{recording, "sub"}, 'sub-003')
        rad_ax = 'Z';
    else
        rad_ax = 'Y';
    end

    load(fullfile(meta_data{recording, "analysed_data_loc"}, ...
        strcat('opti_data_t_', extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'), '.mat')));
  
    marker_slots = readtable(fullfile(meta_data{recording, "raw_data_loc"}, 'optitrack_marker_slots.csv'));
    table_of_info = fullfile(meta_data{recording, "raw_data_loc"}, 'scanner_cast_table_of_info.csv');
    clear R T

    % Set marker heights
    marker_height_cm = zeros(height(marker_slots),1);
    marker_height_cm(strcmp(marker_slots.height, "tall")) = 56.19;
    marker_height_cm(strcmp(marker_slots.height, "mid")) = 49.69;
    marker_height_cm(strcmp(marker_slots.height, "short")) = 40.73;

    if strcmp(meta_data{recording, "sub"}, 'sub-003')
        marker_height_cm(:) = 35;
    end
    
    % Get marker positions in MRI coordinates
    MarkerPosMRI = getMarkerPosInMRIcoords(table_of_info,...
        marker_slots.slot, transpose(marker_height_cm), opti_data, 'Scannercast', rad_ax);
    
    % Get transformation matrices to go from MRI to room coordinates
    [~, ~, R, T] = getMagPosOriOverTime(opti_data, MarkerPosMRI, D, 'Scannercast');
    
    % Find head center over time
    ctx = gifti(D.inv{1}.mesh.tess_ctx);
    cortex_center_MRI = transpose(mean(ctx.vertices,1));
    head_orientation_MRI = [0 1 0]';
    clear ctx
    cortex_center_Room = zeros(3, size(T,2));
    head_orientation_Room = zeros(3, size(T,2));
    yaw_pitch_roll = zeros(size(T,2), 3);
    for tt = 1:size(T,2)
        cortex_center_Room(:,tt) = R(:,:,tt)*cortex_center_MRI + T(:,tt);
        head_orientation_Room(:,tt) = R(:,:,tt)*head_orientation_MRI;
    end

    clear T R

    %% Interpolate missing position data for cortex centre and head orientation

    % Interpolate optitrack data as best as possible
    cortex_center = transpose(cortex_center_Room);
            
    good_opt_data = opti_data.Scannercast.RigidBody{:,10} > 0;
    interp_cortex_center = cortex_center;
    
    % First, assume any missing data at the beginning of the recording by
    % repeating first non-zero datapoint
    first_good_val = find(good_opt_data, 1);
    interp_cortex_center(1:first_good_val-1,:) = repmat(cortex_center(first_good_val,:), first_good_val-1, 1);

    % Ditto for last recording
    last_good_val = find(good_opt_data, 1, 'last');
    interp_cortex_center(last_good_val+1:end,:) = repmat(cortex_center(last_good_val,:), size(cortex_center,1)-last_good_val, 1);
        
    % Then interpolate later missing data
    % - choose mising data points to interpolate over
    [~, step_in] = findpeaks(diff(~good_opt_data));
    [~, step_out] = findpeaks(-diff(~good_opt_data));
    step_in = step_in + 1;
    if first_good_val ~= 1
        step_out = step_out(2:end);
    end
    if last_good_val ~= 1
        step_in = step_in(1:end-1);
    end

    % Shift step_out back by 13 points (consistently too small, I think due
    % to the previous interpolation of the optitrack data from sampling at
    % 120 Hz to 1500 Hz, which is roughly 13 samples)
    % Account for lower sampling rate (375 Hz) of sub002
    if strcmp(meta_data{recording, 'sub'}, 'sub002')
        for gap = 1:length(step_out)
            good_opt_data(step_out(gap):step_out(gap)+4) = 0;
        end
        step_out = step_out + 4;
    else
        for gap = 1:length(step_out)
            good_opt_data(step_out(gap):step_out(gap)+26) = 0;
        end
        step_out = step_out + 26;
    end

    good_opt_data = int8(good_opt_data);
        
    % Fit gap
    for kk = 1:length(step_in)
        k = step_in(kk):step_out(kk);

        if step_in(kk)-6*D.fsample > 0 && step_out(kk)+6*D.fsample < size(interp_cortex_center,1)
            x = step_in(kk)-6*D.fsample:step_out(kk)+6*D.fsample;
            y =interp_cortex_center(step_in(kk)-6*D.fsample:step_out(kk)+6*D.fsample,:);
        elseif step_in(kk)-6*D.fsample < 0
            x = 1:step_out(kk)+6*D.fsample;
            y = interp_cortex_center(1:step_out(kk)+6*D.fsample,:);
        elseif step_out(kk)+6*D.fsample > size(interp_cortex_center,1)
            x = step_in(kk)-6*D.fsample:size(interp_cortex_center,1);
            y = interp_cortex_center(step_in(kk)-6*D.fsample:end,:);
        end
        
        % Just select good data
        good_data_in_section = good_opt_data(x) >= 1;
        x = x(good_data_in_section);
        y = y(good_data_in_section,:);
    
        for cc = 1:3
            fobj = fit(x', y(:,cc), 'pchip');
            interp_cortex_center(k, cc) = fobj(k);
        end
    end

    clear cortex_center_Room cortex_center_MRI cortex_center kk k step_in step_out good_data_in_section

    %% Plot trajectory

    figure; hold on; grid on; box on;
    plot3(interp_cortex_center(good_opt_data==1,1)/10, ...
        interp_cortex_center(good_opt_data==1,2)/10, ...
        interp_cortex_center(good_opt_data==1,3)/10, '.', 'color', [0.2666, 0.4471, 0.7686]);
    if any(good_opt_data == 0)
        plot3(interp_cortex_center(good_opt_data==0,1)/10, ...
            interp_cortex_center(good_opt_data==0,2)/10, ...
            interp_cortex_center(good_opt_data==0,3)/10, ...
            '.', 'color', [0.5647, 0.6706, 0.8627], 'MarkerSize', 0.5);
    end

    zlim([-2 2]*1e2);
    ylim([0, 3.1]*1e2);
    xlim([-1.5, 1.5]*1e2);
    view(180,0);
    daspect([1 1 1]);
    set(gcf, 'Position', [680   670   497   308]);
    set(gca, 'FontSize', 14);
    xlabel('Left-Right (cm)')
    ylabel('Up-Down (cm)')
    zlabel('Forward-Back (cm)');
    set(gcf, 'color', 'w'); 

    save_name = sprintf('%s_trajectory', extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'));
    print(fullfile(meta_data{recording, "results_save_loc"}, save_name),'-dpng','-r300');

    %% Plot as histograms

    % Get times where tones were playing
    trigger = D(indchannel(D, meta_data{recording, "AudioTrig"}),:,1);
    start_index = find(diff(trigger) > 2*range(trigger)/3, 1);
    end_index = find(-diff(trigger) > 2*range(trigger)/3, 1, "last");

    opt_data_keep = find(good_opt_data==1);
    opt_data_keep = opt_data_keep((opt_data_keep >= start_index) & (opt_data_keep <= end_index));
    displacement = interp_cortex_center(opt_data_keep,:);
    displacement = displacement - displacement(1,:);

    rotation = zeros(length(opt_data_keep), 3);
    a = [0; 0; 1];
    for tt = 1:length(opt_data_keep)
        v = cross(a, head_orientation_Room(:, opt_data_keep(tt)));
        s = sqrt(sum(v.^2, 1));
        c = dot(a, head_orientation_Room(:, opt_data_keep(tt)));
        vx = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
        R = eye(3) + vx + vx^2*(1-c)/(s^2);
        rotation(tt,:) = rotm2eul(R, 'YXZ');
    end

    % Choose some colour palettes 
    cp_rot = [255, 254, 203; 243, 119, 72; 146, 55, 77]./255;
    cp_disp = [184, 226, 200; 99.5, 169.5, 157.5; 15, 113, 115]./255;
    cp_B = [236, 190, 180]./255;

    order_disp = {'Left-Right', 'Up-Down', 'Forward-Back'};
    order_rot = {'Yaw', 'Pitch', 'Roll'};

    save_str = extractBefore(meta_data{recording, "raw_data_name"}, '.lvm');

    for ii = 1:size(displacement, 2)
        figure;
        histogram(displacement(:,ii)/10, -200:10:200, ...
            'Normalization', 'probability', 'FaceColor', cp_disp(ii,:), 'FaceAlpha', 1);
        xlim([-200, 200]);
        if ii == 2
            ylim([0, 1]);
        else
            if contains(meta_data{recording, "raw_data_name"}, 'walking')
                ylim([0, 0.2]);
            else
                ylim([0, 1]);
            end
        end

        set(gca,'FontSize',18);
        xlabel('Distance (cm)','FontSize',20);
        ylabel('Frequency','FontSize',20);
        set(gcf, 'Position', [680   650   467   346]);
        grid on;
        title(order_disp{ii}, 'FontSize', 20);

        print(fullfile(meta_data{recording, "results_save_loc"}, sprintf('%s_displacement_%s', save_str, order_disp{ii})),'-dpng','-r300');
    end

    for ii = 1:size(rotation, 2)
        figure;
        ax = polaraxes;
        polarhistogram(rotation(:,ii),25,...
            'Normalization','probability', 'FaceColor', cp_rot(ii,:));
        thetalim([-180 180]);
        ax.FontSizeMode = 'manual';
        ax.FontSize = 18;
        ax.RAxis.FontSize = 16;
        title(order_rot{ii});
        set(gcf, 'Position', [982   234   431   340]);
        print(fullfile(meta_data{recording, "results_save_loc"}, sprintf('%s_rotation_%s', save_str, order_rot{ii})),'-dpng','-r300');
    end

    figure;
    histogram(reshape(D(indchantype(D, 'MEGMAG', 'GOOD'),start_index:end_index,1)*1e-6, [], 1), ...
        -15:0.5:15, 'Normalization', 'probability','FaceColor', cp_B, 'FaceAlpha', 1);
    xlim([-15, 15]);
    if contains(meta_data{recording, "raw_data_name"}, 'walk')
        ylim([0, 0.16]);
    else
        ylim([0, 0.4]);
    end
    set(gca,'FontSize',18);
    xlabel('Recorded Field (nT)','FontSize',22);
    ylabel('Frequency','FontSize',22);
    grid on;
    set(gcf, 'Position', [1128, 244, 832, 638]);
    print(fullfile(meta_data{recording, "results_save_loc"}, sprintf('%s_mag_field_hist', save_str)),'-dpng','-r300');

    % Area covered
    [~, area] = boundary(interp_cortex_center(opt_data_keep,[1,3]));

    fprintf('%s: %s\n Area covered (m^2): %.2f\n Left-Right range (m): %.2f - %.2f (%.2f)\n Forward-Back (m): %.2f - %.2f (%.2f)\n', ...
        meta_data{recording, "sub"}, meta_data{recording, "raw_data_name"}, area/1e6, min(displacement(:,1))/1e3, max(displacement(:,1))/1e3, ...
        range(displacement(:,1))/1e3, min(displacement(:,3))/1e3, max(displacement(:,3))/1e3, range(displacement(:,3))/1e3);

    %% Plot speed

    % Look at each chunk individually to avoid edges
    [~, step_in] = findpeaks(double(diff(good_opt_data(start_index:end_index))));
    [~, step_out] = findpeaks(double(-diff(good_opt_data(start_index:end_index))));
    step_in = step_in + 1;
    first_good_val = find(good_opt_data(start_index:end_index), 1);
    if isempty(step_in) || step_in(1) ~= first_good_val
        step_in = cat(1, first_good_val, step_in);
    end
    last_good_val = find(good_opt_data(start_index:end_index), 1, 'last');
    if last_good_val == length(good_opt_data(start_index:end_index))
        step_out = cat(1, step_out, length(good_opt_data(start_index:end_index)));
    end

    speed = [];
    for chunk = 1:size(step_in, 1)
        speed = cat(1, speed, diff(interp_cortex_center((start_index + step_in(chunk)):(start_index + step_out(chunk)), :))/mean(diff(D.time)));
    end

    for ii = 1:size(displacement, 2)
        figure;
        histogram(speed(:,ii)/10, -75:2:75, ...
            'Normalization', 'probability', 'FaceColor', cp_disp(ii,:), 'FaceAlpha', 1);
        xlim([-75, 75]);
        if ii == 2
            ylim([0, 0.4]);
        else
            if contains(meta_data{recording, "raw_data_name"}, 'walk')
                ylim([0, 0.05]);
            else
                ylim([0, 0.4]);
            end
        end

        set(gca,'FontSize',18);
        xlabel('Velocity (cm/s)','FontSize',20);
        ylabel('Frequency','FontSize',20);
        set(gcf, 'Position', [680   650   467   346]);
        grid on;
        title(order_disp{ii}, 'FontSize', 20);

        print(fullfile(meta_data{recording, "results_save_loc"}, sprintf('%s_speed_%s', save_str, order_disp{ii})),'-dpng','-r300');
    end

    % Rate of change of magnetic field
    figure;
    histogram(reshape(diff(D(indchantype(D, 'MEGMAG', 'GOOD'),start_index:end_index,1)*1e-6, 1, 2)./mean(diff(D.time)), [], 1), ...
        -8:0.5:8, 'Normalization', 'probability','FaceColor', cp_B, 'FaceAlpha', 1);
    xlim([-8, 8]);
    if contains(meta_data{recording, "raw_data_name"}, 'walk')
        ylim([0, 0.2]);
    else
        ylim([0, 0.7]);
    end
    % ylim([0 0.15]);
    set(gca,'FontSize',18);
    xlabel('Recorded Field Temporal Gradient (nT/s)','FontSize',22);
    ylabel('Frequency','FontSize',22);
    grid on;
    set(gcf, 'Position', [1128, 244, 832, 638]);
    print(fullfile(meta_data{recording, "results_save_loc"}, sprintf('%s_mag_field_rate_of_change_hist', save_str)),'-dpng','-r300');

end

%% Magnetic field histograms

closed_loop = false; % Boolean - plot closed loop or open loop

figure; 
ax = subplot(1,1,1); hold on; grid on; box on;
cmap = linspecer(4);
markers = {'s', 'o', '^', 'd'};
counter = 0;

if closed_loop
    recording_order = {"sub-001_task-walkingClosed_meg.lvm", ...
        "sub-002_task-walkingClosed_run-001_meg.lvm", ...
        "sub-002_task-walkingClosed_run-002_meg.lvm", "sub-003_task-walkingClosed_meg.lvm"};
else
    recording_order = {"sub-001_task-walkingOpen_meg.lvm", ...
        "sub-002_task-walkingOpen_run-001_meg.lvm", ...
        "sub-002_task-walkingOpen_run-002_meg.lvm", "sub-003_task-walkingOpen_meg.lvm"};
end

recording_order_name = {'1)', '2a)', '2b)', '3)'};

figure; 
t = tiledlayout("vertical");
h = gobjects(length(recording_order),1);
for ii = 1:length(recording_order)
    h(ii) = nexttile; hold on; grid on; box on;
end

for recording = 1:length(recording_order)
    
    rec_idx = find(contains(meta_data.raw_data_name, recording_order{recording}));

    h1 = figure;
    
    cd(meta_data{rec_idx, "analysed_data_loc"});
    D = spm_eeg_load(['t_', char(meta_data{rec_idx,"raw_data_name"})]);

    trigger = D(indchannel(D, meta_data{rec_idx, "AudioTrig"}),:,1);
    start_index = find(diff(trigger) > 2*range(trigger)/3, 1);
    end_index = find(-diff(trigger) > 2*range(trigger)/3, 1, "last");

    dat = D(indchantype(D, 'MEGMAG', 'GOOD'),start_index:end_index,1)*1e-6;
    
    hh = histogram(dat(:), -15:0.5:15, 'Normalization', 'probability', 'DisplayName', meta_data{rec_idx, "sub"});
    bin_centers = hh.BinEdges(1:end-1) + diff(hh.BinEdges);

    % just choose bins that aren't empty
    filled_bins = hh.Values ~= 0;

    f = fit(bin_centers(filled_bins)', hh.Values(filled_bins)', 'gauss1');
    scatter(ax, bin_centers(filled_bins), hh.Values(filled_bins), markers{recording}, 'filled', ...
        'MarkerFaceColor', cmap(recording,:), 'DisplayName', sprintf('%s %.2f', recording_order_name{recording}, 2*sqrt(2*log(2))*f.c1));
    x = min(bin_centers(filled_bins)):0.05:max(bin_centers(filled_bins));
    plot(ax, x, f(x'), '-', 'color', cmap(recording,:), 'HandleVisibility', 'off', 'LineWidth', 2);

    [~, max_range_chan] = max(range(dat, 2));

    fprintf('%s: %s\n range (nT): %.3f - %.3f\n single channel range (nT): %.3f - %.3f (%.3f)\n gaussian fit: \n   scale = %.3f, \n   mean = %.3f, \n   std = %.3f\n', ...
        meta_data{rec_idx, "sub"}, meta_data{rec_idx, "raw_data_name"}, min(dat(:)), max(dat(:)), ...
        min(dat(max_range_chan,:)), max(dat(max_range_chan,:)), range(dat(max_range_chan,:)), f.a1, f.b1, f.c1);

    close(h1);

    % Plot time series
    plot(h(recording), D.time(start_index:end_index) - D.time(start_index), dat, 'LineWidth', 1.5);
    C = linspecer(size(dat,1)*2);
    C = cat(1, C(1:ceil(size(C,1)/4),:), C(3*floor(size(C,1)/4):end,:));
    set(h(recording), 'ColorOrder', C);
    xlim(h(recording), [0, D.time(end_index) - D.time(start_index)]);
end

xlim(ax, [-15, 15]);
ylim(ax, [0, 0.16]);
set(ax,'FontSize',12);
xlabel(ax, 'Recorded Field (nT)','FontSize',14);
ylabel(ax, 'Frequency','FontSize',14);
lgd = legend(ax);
title(lgd, 'FWHM (nT)');

figure(ax.Parent);
if closed_loop
    print(fullfile(extractBefore(meta_data{rec_idx, 'results_save_loc'}, '\sub-'), 'Magnetic_field_histogram'),'-dpng','-r300');
else
    print(fullfile(extractBefore(meta_data{rec_idx, 'results_save_loc'}, '\sub-'), 'Magnetic_field_histogram_open_loop'),'-dpng','-r300');
end

% Time Series
for ii = 1:length(h)
    ylim(h(ii), [-15 15]);
    set(h(ii), 'FontSize', 12);
    yl = ylabel(h(ii), recording_order_name{ii}, 'FontSize', 14);
    set(yl,'rotation',0,'VerticalAlignment','middle')
end
ylabel(t, 'B (nT)', 'FontSize', 14);
xlabel(t, 'Time (s)', 'FontSize', 14);
figure(t.Parent);
if closed_loop
    print(fullfile(extractBefore(meta_data{rec_idx, 'results_save_loc'}, '\sub-'), 'Magnetic_field_time_series'),'-dpng','-r300');
else
    print(fullfile(extractBefore(meta_data{rec_idx, 'results_save_loc'}, '\sub-'), 'Magnetic_field_time_series_open_loop'),'-dpng','-r300');
end


%% Magnetic field change per second histograms

closed_loop = false; % Boolean - plot closed loop or open loop

figure; 
ax = subplot(1,1,1); hold on; grid on; box on;
cmap = linspecer(4);
markers = {'s', 'o', '^', 'd'};
counter = 0;

if closed_loop
    recording_order = {"sub-001_task-walkingClosed_meg.lvm", ...
        "sub-002_task-walkingClosed_run-001_meg.lvm", ...
        "sub-002_task-walkingClosed_run-002_meg.lvm", "sub-003_task-walkingClosed_meg.lvm"};
else
    recording_order = {"sub-001_task-walkingOpen_meg.lvm", ...
        "sub-002_task-walkingOpen_run-001_meg.lvm", ...
        "sub-002_task-walkingOpen_run-002_meg.lvm", "sub-003_task-walkingOpen_meg.lvm", };
end
recording_order_name = {'1)', '2a)', '2b)', '3)'};

for recording = 1:length(recording_order)
    
    rec_idx = find(contains(meta_data.raw_data_name, recording_order{recording}));

    h1 = figure;

    cd(meta_data{rec_idx, "analysed_data_loc"});
    D = spm_eeg_load(['t_', char(meta_data{rec_idx,"raw_data_name"})]);

    trigger = D(indchannel(D, meta_data{rec_idx, "AudioTrig"}),:,1);
    start_index = find(diff(trigger) > 2*range(trigger)/3, 1);
    end_index = find(-diff(trigger) > 2*range(trigger)/3, 1, "last");

    dat = diff(D(indchantype(D, 'MEGMAG', 'GOOD'),start_index:end_index,1)*1e-6, 1, 2)./mean(diff(D.time));
    
    hh = histogram(dat(:), -8:0.5:8, 'Normalization', 'probability', 'DisplayName', meta_data{rec_idx, "sub"});
    bin_centers = hh.BinEdges(1:end-1) + diff(hh.BinEdges);

    % just choose bins that aren't empty
    filled_bins = hh.Values ~= 0;

    f = fit(bin_centers(filled_bins)', hh.Values(filled_bins)', 'gauss1');
    scatter(ax, bin_centers(filled_bins), hh.Values(filled_bins), markers{recording}, 'filled', ...
        'MarkerFaceColor', cmap(recording,:), 'DisplayName', sprintf('%s %.2f', recording_order_name{recording}, 2*sqrt(2*log(2))*f.c1));
    x = min(bin_centers(filled_bins)):0.05:max(bin_centers(filled_bins));
    plot(ax, x, f(x'), '-', 'color', cmap(recording,:), 'HandleVisibility', 'off', 'LineWidth', 2);

    close(h1);
end

xlim(ax, [-8, 8]);
ylim(ax, [0, 0.2]);
set(ax,'FontSize',12);
xlabel(ax, 'Recorded Field rate of change (nT/s)','FontSize',14);
ylabel(ax, 'Frequency','FontSize',14);
lgd = legend(ax);
title(lgd, 'FWHM (nT/s)');

figure(ax.Parent);
if closed_loop
    print(fullfile(extractBefore(meta_data{rec_idx, 'results_save_loc'}, '\sub-'), 'Magnetic_field_change_rate_histogram'),'-dpng','-r300');
else
    print(fullfile(extractBefore(meta_data{rec_idx, 'results_save_loc'}, '\sub-'), 'Magnetic_field_change_rate_histogram_open_loop'),'-dpng','-r300');
end