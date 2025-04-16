clear; close all;

addpath('C:\Users\smellor\Documents\GitHub\optitrack');

addpath('C:\Users\smellor\Documents\GitHub\spm\');
spm('defaults', 'eeg');

addpath('C:\Users\smellor\Documents\GitHub\BrewerMap')
colormap123 = colormap(flipud(brewermap(64,'RdBu')));
addpath('C:\Users\smellor\Documents\GitHub\icp');

addpath('C:\Users\smellor\Documents\GitHub\my_repos\walking_opmeg_erf\helper_functions')

addpath('C:\Users\smellor\Documents\GitHub\linspecer');

%% Format meta data to do analysis

cd('E:\Data\Neuro1\Auditory\anonymised_for_sharing');

subIDs = {'sub-001', 'sub-002', 'sub-003'};

meta_data = table('Size', [0,7], 'VariableTypes', {'string', 'string', 'string', 'string', 'string', 'string', 'string'},...
    'VariableNames', {'sub', 'raw_data_name', 'raw_data_loc', 'analysed_data_loc', 'results_save_loc', 'OptiTrig', 'AudioTrig'});

% Format meta data table
for sub = 1:length(subIDs)

    % Set audio and optitrack trigger names
    if strcmp(subIDs{sub}, 'sub-001')
        optitrig = 'AI8';
        audiotrig = 'AI16';
    elseif strcmp(subIDs{sub}, 'sub-003')
        optitrig = 'AI16';
        audiotrig = 'AI8';
    elseif strcmp(subIDs{sub}, 'sub-002')
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
        % Fill in table
        meta_data = [meta_data; {subIDs{sub}, lvmfiles{bb}, fpathRoot3, ...
            fullfile(cd, 'analysedData', subIDs{sub}), ...
            fullfile(cd, 'results', subIDs{sub}), optitrig, audiotrig}];
    end
end


clearvars -except meta_data colormap123


%% Load data

for recording = 1:size(meta_data,1)

    if strcmp(meta_data{recording, "sub"}, 'sub-002')
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
        D = badchannels(D, selectchannels(D, 'regexp_(^(28|48|41)-.*)'), 1);
    elseif strcmp(meta_data{recording, "sub"}, 'sub-003')
        D = badchannels(D, selectchannels(D, 'regexp_(^(2|19|28|14|31|59|1|4)-.*)'), 1);
        D = badchannels(D, selectchannels(D, 'regexp_(^5-.*-Y$)'), 1);
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
    xlim([0 40]);
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
    xlim([0 40]);
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
    
%% Plot time series

for recording = 1:size(meta_data,1)
    
    if strcmp(meta_data{recording, "sub"}, 'sub-002')
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
            plot(1e3*D.time, t(plot_order(chan),:), 'color', max_colour(plot_order(chan),:), 'LineWidth', 2);
        end
        a = tinv(1-0.025/(size(D,2)*length(plot_order)), size(D,3)-1);
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
        cfg.xlim = [100, 100]*1e-3;
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
            plot(1e3*D.time, dat(plot_order(chan),:), 'color', max_colour(plot_order(chan),:), 'LineWidth', 2);
        end
        xlim([-100 400]);
        ylim([-450 450]);
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
        cfg.zlim      = [-250, 250];
        ft_topoplotER(cfg, avdata)
        set(gcf, 'Position', [994   704   404   274]);
        set(gca, 'FontSize', 24);

        save_name = sprintf('%s_average_topography', extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'));
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
                    "sub-002_task-walkingClosed_meg.lvm", "sub-003_task-walkingClosed_run-001_meg.lvm", ...
                    "sub-003_task-walkingClosed_run-002_meg.lvm"};
            else
                recording_order = {"sub-001_task-seatedClosed_meg.lvm", ...
                    "sub-002_task-seatedClosed_meg.lvm", "sub-003_task-seatedClosed_meg.lvm"};
            end
        else
            if walking
                recording_order = {"sub-001_task-walkingOpen_meg.lvm", ...
                    "sub-002_task-walkingOpen_meg.lvm", "sub-003_task-walkingOpen_run-001_meg.lvm", ...
                    "sub-003_task-walkingOpen_run-002_meg.lvm"};
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
        xlim([0 40]);
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
                S.channels = S.D1.chanlabels(indchantype(D, 'MEGMAG', 'GOOD'));
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
        xlim([0 40]);
        ylim([0 30]);
        legend(pl, {'HFC', 'AMM spatial', 'AMM with temporal'});
        xlabel('Frequency (Hz)');
        ylabel('Shielding Factor (dB)');

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

    if strcmp(subIDs{sub}, 'sub-002')
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
        alpha = 0.025/(size(X, 2)*size(X,1));
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

    if strcmp(meta_data{recording, "sub"}, 'sub-002')
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

    if strcmp(meta_data{recording, "sub"}, 'sub-002')
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
        "sub-002_task-walkingClosed_meg.lvm", "sub-003_task-walkingClosed_run-001_meg.lvm", ...
        "sub-003_task-walkingClosed_run-002_meg.lvm"};
else
    recording_order = {"sub-001_task-walkingOpen_meg.lvm", ...
        "sub-002_task-walkingOpen_meg.lvm", "sub-003_task-walkingOpen_run-001_meg.lvm", ...
        "sub-003_task-walkingOpen_run-002_meg.lvm"};
end

recording_order_name = {'A)', 'B)', 'C)', 'D)'};

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
        "sub-002_task-walkingClosed_meg.lvm", "sub-003_task-walkingClosed_run-001_meg.lvm", ...
        "sub-003_task-walkingClosed_run-002_meg.lvm"};
else
    recording_order = {"sub-001_task-walkingOpen_meg.lvm", ...
        "sub-002_task-walkingOpen_meg.lvm", "sub-003_task-walkingOpen_run-001_meg.lvm", ...
        "sub-003_task-walkingOpen_run-002_meg.lvm"};
end
recording_order_name = {'A)', 'B)', 'C)', 'D)'};

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