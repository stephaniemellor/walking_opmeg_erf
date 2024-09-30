clear; close all;

addpath('C:\Users\smellor\Documents\GitHub\optitrack');

addpath('C:\Users\smellor\Documents\GitHub\spm\');
spm('defaults', 'eeg');

addpath('C:\Users\smellor\Documents\GitHub\BrewerMap')
colormap123 = colormap(flipud(brewermap(64,'RdBu')));
addpath('C:\Users\smellor\Documents\GitHub\icp');

cd('E:\Data\Neuro1\Auditory');

%% Format meta data to do analysis

subIDs = {'sub-003'};

meta_data = table('Size', [0,9], 'VariableTypes', {'string', 'string', 'string', 'string', 'string', 'string', 'string', 'string', 'string'},...
    'VariableNames', {'sub', 'date', 'raw_data_name', 'raw_data_loc', 'analysed_data_loc', 'results_save_loc', 'stim_data_fname', 'OptiTrig', 'AudioTrig'});

% Format meta data table
for sub = 1:length(subIDs)

    % Set audio and optitrack trigger names
    if strcmp(subIDs{sub}, 'sub-003')
        optitrig = 'AI8';
        audiotrig = 'AI16';
    elseif strcmp(subIDs{sub}, 'sub-004')
        optitrig = 'AI16';
        audiotrig = 'AI8';
    else
        error('Please set trigger channel names for participant %s', subIDs{sub})
    end

    % Search folders for files
    fpathRoot2 = fullfile(cd, subIDs{sub});

    % Find scan date
    listing = dir(fpathRoot2);
    listing = extractfield(listing, 'name');
    folders = {};
    for ll = 1:length(listing)
        if regexp(listing{ll}, '^[0-9]{8}$')
            folders = cat(1, folders, listing{ll});
        end
    end

    % Loop through scan dates
    for ll = 1:length(folders)
        fpathRoot3 = fullfile(fpathRoot2, folders{ll}, 'opm', 'rawData');

        % Find raw data name
        listing = dir(fpathRoot3);
        listing = extractfield(listing, 'name');

        % Just choose lvm files
        lvmfiles = listing(endsWith(listing, '.lvm'));
            
        for bb = 1:length(lvmfiles)

            % Only add to table if a positions file exists and name doesn't contain noise
            if isfile(fullfile(fpathRoot3, strrep(lvmfiles{bb}, '.lvm', '_positions.tsv'))) && ~contains(lvmfiles{bb}, 'noise')

                % Find corresponding stim data file
                fpath_stim = fullfile(fpathRoot2, folders{ll}, 'stim');
                listing = dir(fpath_stim);
                listing = extractfield(listing, 'name');
                open = contains(lvmfiles{bb}, 'open');
                seated = contains(lvmfiles{bb}, 'seated');
                if open
                    open_or_closed = 'open';
                else
                    open_or_closed = 'closed';
                end
                if seated
                    seated_or_walking = 'seated';
                else
                    seated_or_walking = 'walking';
                end
                stimfile = listing(contains(listing, open_or_closed) & contains(listing, seated_or_walking));
                
                % Fill in table
                meta_data = [meta_data; {subIDs{sub}, folders{ll}, lvmfiles{bb}, ...
                    fpathRoot3, fullfile(fpathRoot2, folders{ll}, 'opm', 'analysedData'), ...
                    fullfile(fpathRoot2, folders{ll}, 'opm', 'results'), fullfile(fpath_stim, stimfile{1}), optitrig, audiotrig}];
            end
        end
    end
end

clearvars -except meta_data colormap123


%% Load data

for recording = 1:size(meta_data,1)

    rawDataPath = char(meta_data{recording, "raw_data_loc"});
    analysedDataPath = char(meta_data{recording, "analysed_data_loc"});
    resultsPath = char(meta_data{recording, "results_save_loc"});
    
    cd(rawDataPath)
    
    fname = char(extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'));
    
    if isfile(fullfile(analysedDataPath, [fname, '.mat']))
        D = spm_eeg_load(fullfile(analysedDataPath, [fname, '.mat']));
    else
        S = [];
        S.positions = sprintf('%s_positions.tsv', fname);
        S.data = fullfile(cd, sprintf('%s.lvm', fname));
        S.path = analysedDataPath;
        S.sMRI = fullfile(char(extractBefore(meta_data{recording, "raw_data_loc"},'sub-')),...
            char(meta_data{recording, "sub"}), 'mri', sprintf('%s_structural_anon.nii', meta_data{recording, 'sub'}));
        D = spm_opm_create(S);
    end
    
    cd(analysedDataPath)

    %% Preprocess data

    % Badchannels
    
    if strcmp(meta_data{recording, "sub"}, 'sub-003')
        D = badchannels(D, selectchannels(D, 'regexp_(^(27|8|31|59|6|5|19|41|24|30|17)-.*)'), 1);
    elseif strcmp(meta_data{recording, "sub"}, 'sub-004')
        D = badchannels(D, selectchannels(D, 'regexp_(^(2|19|28|14|31|59|1|4)-.*)'), 1);
        D = badchannels(D, selectchannels(D, 'regexp_(^5-.*-Y$)'), 1);
    else
        warning('Badchannels not set for participant %s', meta_data{recording, "sub"});
    end
    save(D);

    % PSD
    S = [];
    S.D = D;
    S.channels = D.chanlabels(indchantype(D, 'MEGMAG', 'GOOD'));
    S.plot = 1;
    S.triallength = 10e3;
    %S.plotbad = 1;
    [po, freq] = spm_opm_psd(S);

    % Plot spatially
    if isfile(sprintf('%s_2Dlayout.mat', extractBefore(meta_data{recording, "raw_data_name"}, ' Array 1.lvm')))
        load(sprintf('%s_2Dlayout.mat', extractBefore(meta_data{recording, "raw_data_name"}, ' Array 1.lvm')));
    else
        fid = fiducials(D);
        fid_struct = struct('NAS', fid.fid.pnt(contains(fid.fid.label, 'nas'),:), ...
            'LPA', fid.fid.pnt(contains(fid.fid.label, 'lpa'),:), ...
            'RPA', fid.fid.pnt(contains(fid.fid.label, 'rpa'),:));
        pos = D.sensors('MEG').coilpos;
        lay = spm_get_anatomical_layout(D.sensors('MEG').coilpos(endsWith(D.sensors('MEG').label, '-Y'),:), ...
            D.sensors('MEG').label(endsWith(D.sensors('MEG').label, '-Y')),...
            double(gifti(D.inv{1}.mesh.tess_scalp).vertices), fid_struct, 0);
        save(sprintf('%s_2Dlayout', extractBefore(meta_data{recording, "raw_data_name"}, ' Array 1.lvm')), 'lay');
    end

    S = [];
    S.trialength = 1000;
    S.D = D;
    D1 = spm_eeg_epochs(S);
    
    data = fttimelock(D1);
    cfg = [];
    cfg.channel = D.chanlabels(indchantype(D, 'MEGMAG', 'GOOD'));
    avdata = ft_timelockanalysis(cfg, data);
    psd_avdata = avdata;
    psd_avdata.avg = repmat(sum((po((freq > 2) & (freq < 10),:)*1e-15*1e9).^2, 1)', 1, length(psd_avdata.time));
    
    figure;
    cfg = [];
    cfg.layout    = lay;
    cfg.colorbar  = 'EastOutside';
    cfg.colorbartext = 'Power (nT^2)';
    cfg.zlim      = [0 max(psd_avdata.avg(ismember(psd_avdata.label, lay.label),1))];
    cfg.colormap  = colormap123(ceil(size(colormap123,1)/2):end,:);
    cfg.xlim = [1, 1]*1e-3;
    cfg.comment = 'no';
    cfg.figure = gca;
    set(gca, 'FontSize', 16);
    ft_topoplotER(cfg, psd_avdata)
    print(fullfile(meta_data{recording, "results_save_loc"}, ...
        sprintf('%s_raw_2_to_10_Hz_power_spatial', extractBefore(meta_data{recording, "raw_data_name"}, ' Array 1.lvm'))),'-dpng','-r300');

    delete(D1);
    clear psd_avdata data D1 po freq

    
    % Sync with optitrack
    if isfile(['opti_data_', D.fname])
        load(['opti_data_', D.fname]);
        D = spm_eeg_load(['t_', D.fname]);
    else
        cfg = [];
        fname = D.fname;
        cfg.filename = [rawDataPath, '\', fname(1:end-4), '_optitrack.csv'];
        opti_data = readRigidBody(cfg);
        [opti_data, D] = syncOptitrackAndOPMdata(opti_data,D,'TriggerChannelName',meta_data{recording, "OptiTrig"});
        save(['opti_data_', D.fname], 'opti_data');
    end

    % Plot opm recordings with position and rotation
    figure;
    t = tiledlayout(3,1);
    nexttile; 
    plot(D.time, 1e-6*D(indchantype(D, 'MEGMAG', 'GOOD'),:,1));
    ylabel('B (nT)');
    ylim([-15 15]);
    xticklabels({});
    xlim([70, max(D.time)])
    set(gca, 'FontSize', 16);
    grid on; box on;
        
    nexttile;
    if strcmp(opti_data.cfg.LengthUnits, 'Meters')
        plot(D.time, 1e3*(opti_data.Scannercast.RigidBody{:,7:9}-opti_data.Scannercast.RigidBody{1,7:9}), 'LineWidth', 2);
    else
        plot(D.time, opti_data.Scannercast.RigidBody{:,7:9}-opti_data.Scannercast.RigidBody{1,7:9}, 'LineWidth', 2);
    end
    ylabel('Displacement (mm)');
    xlim([70, max(D.time)])
    xticklabels({});
    set(gca, 'FontSize', 16);
    ylim([-1200 1200]);
    grid on; box on;
    legend({'Left-Right', 'Up-Down', 'Door-Screen'}, 'location', 'eastoutside');
        
    nexttile;
    plot(D.time, 180*unwrap(quat2eul(opti_data.Scannercast.RigidBody{:,[6,3:5]}, 'XYZ'))/pi, 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Rotation (deg)');
    xlim([70, max(D.time)])
    set(gca, 'FontSize', 16);
    ylim([-360 360]);
    grid on; box on;
    legend({'Pitch', 'Yaw', 'Roll'}, 'location', 'eastoutside');
        
    t.TileSpacing = 'compact';
    set(gcf, 'Position', [680   446   932   652]);
    print(fullfile(meta_data{recording, "results_save_loc"}, ...
        sprintf('%s_all_time_series', extractBefore(meta_data{recording, "raw_data_name"}, ' Array 1.lvm'))),'-dpng','-r300');
    
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
    proc_step_names = {'No spatial filter', 'HFC', 'HFC with first order gradients', 'AMM'};
    DD = cell(1,4);
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
    
    % AMM
    if isfile(['m', D.fname])
        DD{4} = spm_eeg_load(['m', D.fname]);
    else
        S = [];
        S.D = D;
        S.corrLim = 0.95;
        S.reducerank = 1;
        DD{4} = spm_opm_amm(S);
    end
    
    % Plot shielding factors
    for proc_step = 2:4
        S = [];
        S.D1 = DD{1};
        S.D2 = DD{proc_step};
        S.channels = D.chanlabels(indchantype(D, 'MEGMAG', 'GOOD'));
        S.plot = 1;
        S.triallength = 10e3;
        S.dB = 1;
        spm_opm_rpsd(S);
        xlim([0 100]);
        ylim([-20 50])
        print(fullfile(meta_data{recording, "results_save_loc"}, ...
            sprintf('%s_shielding_factor_%s', ...
            extractBefore(meta_data{recording, "raw_data_name"}, ' Array 1.lvm'), ...
            strrep(proc_step_names{proc_step}, ' ', '_'))),'-dpng','-r300');
    end
        
    
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
end
    
%% Plot time series
    
for recording = 1:size(meta_data,1)

    start_string = {'e_ffft_', 'e_hffft_', 'e_h2ffft_', 'e_mffft_'};
    for pp = 1:4
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

        % Save maximum channel for later
        if pp == 4
            max_chan = plot_order(end);
            meg_chans = D.chanlabels(indchantype(D, 'MEGMAG', 'GOOD'));
            max_chan = meg_chans{max_chan};
        end
    
        % Plot t stat
        figure; hold on; grid on; box on;
        for chan = 1:size(t,1)
            plot(D.time, t(plot_order(chan),:), 'color', max_colour(plot_order(chan),:));
        end
        a = tinv(1-0.025/size(D,2), size(D,3)-1);
        l1 = plot([-0.1 0.4], [a a], 'b--');
        plot([-0.1 0.4], [-a -a], 'b--');
        legend(l1, 'Significance Threshold')
        xlim([-0.1 0.4]);
        ylim([-25 25]);
        xlabel('Time (s)');
        ylabel('t-stat');
        set(gca, 'FontSize', 14);
        fname = D.fname;
        
        if pp == 1
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'no_amm');
        elseif pp == 2
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'hfc');
        elseif pp == 3
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'hfc_with_gradients');
        else
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'amm');
        end
        if ~exist(save_loc, 'dir')
            mkdir(save_loc);
        end
        print(fullfile(save_loc, sprintf('%s_t_stat_time_series', extractBefore(meta_data{recording, "raw_data_name"}, ' Array 1.lvm'))),'-dpng','-r300');
    
        % Plot topography of t-stat at 100 ms
        if isfile(sprintf('%s_2Dlayout.mat', extractBefore(meta_data{recording, "raw_data_name"}, ' Array 1.lvm')))
            load(sprintf('%s_2Dlayout.mat', extractBefore(meta_data{recording, "raw_data_name"}, ' Array 1.lvm')));
        else
            fid = fiducials(D);
            fid_struct = struct('NAS', fid.fid.pnt(contains(fid.fid.label, 'nas'),:), ...
                'LPA', fid.fid.pnt(contains(fid.fid.label, 'lpa'),:), ...
                'RPA', fid.fid.pnt(contains(fid.fid.label, 'rpa'),:));
            pos = D.sensors('MEG').coilpos;
            lay = spm_get_anatomical_layout(D.sensors('MEG').coilpos(endsWith(D.sensors('MEG').label, '-Y'),:), ...
                D.sensors('MEG').label(endsWith(D.sensors('MEG').label, '-Y')),...
                double(gifti(D.inv{1}.mesh.tess_scalp).vertices), fid_struct, 0);
            save(sprintf('%s_2Dlayout', extractBefore(meta_data{recording, "raw_data_name"}, ' Array 1.lvm')), 'lay');
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
        cfg.zlim      = [-17 17];
        cfg.colormap  = colormap123;
        cfg.xlim = [100, 100]*1e-3;
        cfg.comment = 'no';
        cfg.figure = gca;
        set(gca, 'FontSize', 16);
        ft_topoplotER(cfg, tavdata)
        print(fullfile(save_loc, sprintf('%s_t_stat_topography', extractBefore(meta_data{recording, "raw_data_name"}, ' Array 1.lvm'))),'-dpng','-r300');
    
        % Plot average signal
        figure; hold on; grid on; box on;
        dat = mean(D(indchantype(D, 'MEGMAG', 'GOOD'),:,good_trials),3);
        for chan = 1:size(t,1)
            plot(D.time, dat(plot_order(chan),:), 'color', max_colour(plot_order(chan),:));
        end
        xlim([-0.1 0.4]);
        ylim([-550 550]);
        xlabel('Time (s)');
        ylabel('B (fT)');
        grid on;
        set(gca, 'FontSize', 14);
        print(fullfile(save_loc, sprintf('%s_average_time_series', extractBefore(meta_data{recording, "raw_data_name"}, ' Array 1.lvm'))),'-dpng','-r300');
        
        % Topoplot
        figure;
        cfg.figure = gca;
        cfg.colorbartext = 'B (fT)';
        set(gca, 'FontSize', 16);
        cfg.zlim      = [-450, 450];
        ft_topoplotER(cfg, avdata)
        print(fullfile(save_loc, sprintf('%s_average_topography', extractBefore(meta_data{recording, "raw_data_name"}, ' Array 1.lvm'))),'-dpng','-r300');
    end

    close all
end

    %% Source localisation and VOI analysis at auditory cortices for each preprocessing step (using DAiSS)
    
for recording = 1:size(meta_data,1)

    start_string = {'e_ffft_', 'e_hffft_', 'e_h2ffft_', 'e_mffft_'};
    for pp = 1:4
        DD{pp} = spm_eeg_load(char(fullfile(meta_data{recording, "analysed_data_loc"}, ...
            strcat(start_string{pp}, extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'), '.mat'))));
    end

    MD = cell(1, length(DD));

    for pp = 1:length(DD)
        
        path_BF = 'BF.mat';

        % Prepare Data
        S = [];
        S.D = DD{pp}.fullfile;
        S.dir = DD{pp}.path;
        bf_wizard_data(S);

        % Create source space
        S = [];
        S.BF = 'BF.mat';
        S.method = 'grid';
        S.grid.resolution = 5;
        bf_wizard_sources(S);

        % Estimate covariance matrix
        S = [];
        S.BF = 'BF.mat';
        S.woi = [-inf inf];
        S.conditions = 'all';
        S.method = 'cov';
        S.cov.taper = 'none';
        S.cov.foi = [2 40]; % wide band matrix
        S.modality = {'MEG'};
        S.reg = 'mantrunc';
        S.(S.reg).pcadim = 95;
        S.visualise = 1;
        bf_wizard_features(S);

        % Estimate beamformer weights
        S = [];
        S.BF = 'BF.mat';
        S.method = 'lcmv_hippocampus';
        S.(S.method).bilateral = 'all';
        bf_wizard_inverse(S);

        % VOI
        S = [];
        S.BF = 'BF.mat';
        S.method = 'montage';
        S.montage.method = 'svd';
        S.montage.vois{1}.voidef.label = 'rA1';
        S.montage.vois{1}.voidef.pos = [54 -14 8];
        S.montage.vois{1}.voidef.radius = 15;
        S.montage.vois{2}.voidef.label = 'lA1';
        S.montage.vois{2}.voidef.pos = [-54 -14 8];
        S.montage.vois{2}.voidef.radius = 15;
        bf_wizard_output(S);

        % Write to SPM dataset
        S = [];
        S.BF = 'BF.mat';
        S.method = 'meeg';
        S.meeg.prefix = 'M';
        bf_wizard_write(S);

        % Plot
        MD{pp} = spm_eeg_load(['M', DD{pp}.fname]);
        S = [];
        S.D = MD{pp};
        mMD = spm_eeg_average(S);

        good_trials = indtrial(MD{pp}, 'tone', 'GOOD');
        se = std(MD{pp}(:,:,good_trials),[],3)./sqrt(length(good_trials));
        t = mean(MD{pp}(:,:,good_trials),3)./se;

        % plot the results
        figure;
        hold on
        plot(MD{pp}.time, t,'linewidth',2)
        grid on
        
        a = tinv(1-0.025/size(MD{pp},2), length(good_trials)-1);
        l1 = plot([-0.1 0.4], [a a], 'b--');
        plot([-0.1 0.4], [-a -a], 'b--');
        legend([l1], 'Significance Threshold')
        xlim([-0.1 0.4]);
        ylim([-25 25]);
        xlabel('Time (s)');
        ylabel('t-stat');
        set(gca, 'FontSize', 14);

        if pp == 1
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'no_amm');
        elseif pp == 2
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'hfc');
        elseif pp == 3
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'hfc_with_gradients');
        else
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'amm');
        end

        print(fullfile(save_loc, ...
            sprintf('%s_VOI_t_val', extractBefore(meta_data{recording, "raw_data_name"}, ' Array 1.lvm'))),'-dpng','-r300');

        % Source localisation
        S = [];
        S.BF = 'BF.mat';
        S.conditions = {'tone'};
        S.woi = [50 150];
        S.foi = [2 40];
        S.contrast = 1;
        S.method = 'image_power';
        S.image_power.scale = 0; % don't use unit-gain scaling
        S.image_power.logpower = 0; % don't log transform the power
        S.image_power.result = 'bytrial';
        bf_wizard_output(S);

        BF_signal = load('BF.mat');

        S = [];
        S.BF = 'BF.mat';
        S.conditions = {'tone'};
        S.woi = [-150 -50];
        S.foi = [2 40];
        S.contrast = 1;
        S.method = 'image_power';
        S.image_power.scale = 0; % don't use unit-gain scaling
        S.image_power.logpower = 0; % don't log transform the power
        S.image_power.result = 'bytrial';
        bf_wizard_output(S);

        BF_baseline = load('BF.mat');

        % Paired T-test
        pow_diff = zeros([size(BF_signal.output.image(1).val), length(good_trials)]);
        for trial = 1:length(good_trials)
            pow_diff(:,:,trial) = BF_signal.output.image(trial).val.^(1/4) - BF_baseline.output.image(trial).val.^(1/4);
        end
        se = std(pow_diff,[],3)./sqrt(size(pow_diff,3));
        t = mean(pow_diff,3)./se;
        a = tinv(1-0.025/size(pow_diff,2), size(pow_diff,3)-1);

        % Plot
        source = BF_signal.sources.grid;
        source.pos = BF_signal.sources.grid.allpos;
        source = ft_transform_geometry(BF_signal.data.transforms.toNative, source);
        source.pow = nan(size(source.pos, 1), 1);
        source.pow(source.inside) = t;
        source.mask = zeros(size(source.pow));
        source.mask(source.inside) = t > a;
        source.mask = source.mask*0.75;
        
        cfg = [];
        cfg.parameter = {'pow', 'mask'};
        cfg.downsample = 1;
        cfg.showcallinfo = 'no';
        sourceint = ft_sourceinterpolate(cfg, source, ft_read_mri(BF_signal.data.mesh.sMRI, 'dataformat', 'nifti_spm'));
        
        % Plot beamformer result
        maxval = 22;
        
        cfg = [];
        cfg.method        = 'ortho';
        cfg.funparameter  = 'pow';
        cfg.maskparameter = 'mask';
        cfg.funcolorlim   = [0.0 maxval];
        cfg.opacitylim    = [0 1];
        cfg.location = 'max';
        ft_sourceplot(cfg, sourceint);
        print(fullfile(save_loc, ...
            sprintf('%s_cross_brain_t_val', extractBefore(meta_data{recording, "raw_data_name"}, ' Array 1.lvm'))),'-dpng','-r300');

        % Write to nifti
        outvol = spm_vol(BF_signal.data.mesh.sMRI);
        outvol.dt(1) = spm_type('float32');
        outvol.fname= char(fullfile(save_loc, sprintf('%s_BF_power_comp.nii', ...
            extractBefore(meta_data{recording, "raw_data_name"}, ' Array 1.lvm'))));
        outvol = spm_create_vol(outvol);
        Y = sourceint.pow;
        spm_write_vol(outvol, Y);
    end

    % Calculate improvement from prior spatial filters
    for pp = 2:length(DD)
        good_trials = indtrial(MD{pp}, 'tone', 'GOOD');
        ts_diff = (abs(MD{pp}(:,:,good_trials)) - abs(MD{1}(:,:,good_trials)));
        se = std(ts_diff, [], 3)./sqrt(length(good_trials));
        t = mean(ts_diff, 3)./se;

        figure; hold on;
        plot(MD{pp}.time, t, 'linewidth', 2);
        grid on
        
        a = tinv(1-0.025/size(MD{pp},2), size(MD{pp},3)-1);
        l1 = plot([-0.1 0.4], [a a], 'b--');
        plot([-0.1 0.4], [-a -a], 'b--');
        legend(l1, 'Significance Threshold')
        xlim([-0.1 0.4]);
        ylim([-25 25]);
        xlabel('Time (s)');
        ylabel('t-stat');
        set(gca, 'FontSize', 14);

        if pp == 2
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'hfc');
        elseif pp == 3
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'hfc_with_gradients');
        else
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'amm');
         end

        print(fullfile(save_loc, ...
            sprintf('%s_VOI_t_val_filter_improvement', extractBefore(meta_data{recording, "raw_data_name"}, ' Array 1.lvm'))),'-dpng','-r300');
    end

    close all
end
    
%% ROI analysis, dipole

for recording = 1:size(meta_data,1)

    start_string = {'e_ffft_', 'e_hffft_', 'e_h2ffft_', 'e_mffft_'};
    for pp = 1:4
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
        alpha = 0.025/size(X, 2);
        sigt = tinv(1-alpha, df);
        plot([min(DD{pp}.time), max(DD{pp}.time)]*1e3, sigt*[1 1], 'k--');
        plot([min(DD{pp}.time), max(DD{pp}.time)]*1e3, -sigt*[1 1], 'k--');

        xlim([-100 400]);
        ylim([-20 20])
        xlabel('Time (ms)');
        ylabel('t-stat');
        set(gca, 'FontSize', 14);

        if pp == 1
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'no_amm');
        elseif pp == 2
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'hfc');
        elseif pp == 3
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'hfc_with_gradients');
        else
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'amm');
         end

        print(fullfile(save_loc, ...
            sprintf('%s_ROI_t_val_dipole', extractBefore(meta_data{recording, "raw_data_name"}, ' Array 1.lvm'))),'-dpng','-r300');
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
    start_string = {'e_ffft_', 'e_hffft_', 'e_h2ffft_', 'e_mffft_'};
    for pp = 1:4
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
        else
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'amm');
         end

        print(fullfile(save_loc, ...
            sprintf('%s_min_norm_t_val_left', extractBefore(meta_data{recording, "raw_data_name"}, ' Array 1.lvm'))),'-dpng','-r300');

        view ([90 0])             % rotate the object in the view
        camlight('headlight')
        print(fullfile(save_loc, ...
            sprintf('%s_min_norm_t_val_right', extractBefore(meta_data{recording, "raw_data_name"}, ' Array 1.lvm'))),'-dpng','-r300');

    end
end

%% Mismatch negativity

for recording = 1:size(meta_data,1)

    start_string = {'e_ffft_', 'e_hffft_', 'e_h2ffft_', 'e_mffft_'};
    for pp = 1:4
        DD{pp} = spm_eeg_load(char(fullfile(meta_data{recording, "analysed_data_loc"}, ...
            strcat(start_string{pp}, extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'), '.mat'))));
    end

    % Get stim file
    stim = readtable(meta_data{recording, "stim_data_fname"});
    deviants = find(contains(stim.Condition, 'deviant'));
    trial_length = diff(deviants);
    trial_length = cat(1, trial_length, size(D,3) - max(deviants) + 1);

    % Just keep sets with 5 or more standard tones
    deviants = deviants(trial_length >= 5);
    standards = deviants + 4;

    % Plot paired t-test
    for pp = 1:length(DD)
        D = DD{pp};
        
        ds = D(indchantype(D, 'MEGMAG', 'GOOD'), :, deviants) - D(indchantype(D, 'MEGMAG', 'GOOD'), :, standards);
        mds = mean(ds, 3);
        seds = std(ds, [], 3)/sqrt(size(ds,3));
        t = mds./seds;

        [~, tind] = min(abs(D.time - 175*1e-3));
        max_colour = interp1([min(abs(t(:,tind))), max(abs(t(:,tind)))], [0.9, 0], abs(t(:,tind)));
        [~, plot_order] = sort(max_colour, 'descend');
        max_colour = repmat(max_colour, 1, 3);
    
        % Plot t stat
        figure; hold on; grid on; box on;
        for chan = 1:size(t,1)
            plot(D.time, t(plot_order(chan),:), 'color', max_colour(plot_order(chan),:));
        end

        a = tinv(1-0.025/size(D,2), size(ds,3)-1);
        l1 = plot([-0.1 0.4], [a a], 'b--');
        plot([-0.1 0.4], [-a -a], 'b--');
        legend(l1, 'Significance Threshold')
        xlim([-0.1 0.4]);
        ylim([-10 10]);
        xlabel('Time (s)');
        ylabel('t-stat');
        set(gca, 'FontSize', 14);
        fname = D.fname;
        
        if pp == 1
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'no_amm');
        elseif pp == 2
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'hfc');
        elseif pp == 3
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'hfc_with_gradients');
        else
            save_loc = fullfile(meta_data{recording, "results_save_loc"}, 'amm');
        end
        if ~exist(save_loc, 'dir')
            mkdir(save_loc);
        end
        print(fullfile(save_loc, sprintf('%s_MMN_t_stat_time_series', extractBefore(meta_data{recording, "raw_data_name"}, ' Array 1.lvm'))),'-dpng','-r300');
    
    end
end


%% Shielding factor from AMM as a function of channel number

addpath('D:\Data\Neuro1\Auditory');
nopms = 33:2:43;

for recording = 1:size(meta_data,1)
    
    cd(meta_data{recording, "analysed_data_loc"});

    figure; 
    ax = subplot(1,1,1); hold on; grid on; box on;
    D = spm_eeg_load(['ffft_', char(meta_data{recording,"raw_data_name"})]);

    % Subsample sensors
    Damm = cell(1, length(nopms)+1);
    Damm{end} = spm_eeg_load(['m', D.fname]);

    for nchans = 1:length(nopms)
        % Define cost function
        CF = @(x)evenly_space_sensors_CF(D.sensors('MEG').chanpos(endsWith(D.sensors('MEG').label, '-X'),:), x);
        sens_inds = randperm(size(D.sensors('MEG').chanpos(endsWith(D.sensors('MEG').label, '-X'),:),1), nopms(nchans));
        CF_cur = CF(sens_inds);
        for iter = 1:50000
            sens_inds_tmp = randperm(size(D.sensors('MEG').chanpos(endsWith(D.sensors('MEG').label, '-X'),:),1), nopms(nchans));
            CF_new = CF(sens_inds_tmp);
            if CF_new < CF_cur
                sens_inds = sens_inds_tmp;
                CF_cur = CF_new;
            end
        end
    
        figure
        ft_plot_sens(D.sensors('MEG'));
        hold on;
        plot3(D.sensors('MEG').chanpos(sens_inds,1), D.sensors('MEG').chanpos(sens_inds,2), D.sensors('MEG').chanpos(sens_inds,3), 'bx', 'MarkerSize', 20);
    
        % Create array of labels
        chanlabels = cat(2, cellfun(@(x)[x(1:end-1), 'X'], D.sensors('MEG').label(sens_inds), 'UniformOutput', false), ...
            cellfun(@(x)[x(1:end-1), 'Y'], D.sensors('MEG').label(sens_inds), 'UniformOutput', false), ...
            cellfun(@(x)[x(1:end-1), 'Z'], D.sensors('MEG').label(sens_inds), 'UniformOutput', false));
        chanlabels = reshape(chanlabels, 1, []);

        % Subset OPMs
        S = [];
        S.D = D;
        S.channels = chanlabels;
        S.prefix = strcat('p', num2str(nchans));
        Damm{nchans} = spm_eeg_crop(S);

        % Apply AMM
        S = [];
        S.D = Damm{nchans};
        S.corrLim = 0.95;
        S.reducerank = 1;
        Damm{nchans} = spm_opm_amm(S);
    end

    % Add line to figure;
    for nchans = 1:length(Damm)
        S = [];
        if nchans < length(Damm)
            S.D1 = spm_eeg_load([char(strcat('p', num2str(nchans))), D.fname]);
        else
            S.D1 = D;
        end
        S.D2 = spm_eeg_load(['m', S.D1.fname]);
        S.channels = D.chanlabels(indchantype(D, 'MEGMAG', 'GOOD'));
        S.triallength = 10e3;
        S.dB = 1;
        [shield, f] = spm_opm_rpsd(S);
        plot(ax, f, median(shield, 2));
    end

    xlim([0 100]);
    ylim([-20 50])
    print(fullfile(meta_data{recording, "results_save_loc"}, ...
        sprintf('%s_shielding_factor_%s', ...
        extractBefore(meta_data{recording, "raw_data_name"}, ' Array 1.lvm'), ...
        strrep(proc_step_names{proc_step}, ' ', '_'))),'-dpng','-r300');
end

%% Split epochs into groups by field change

% Get channel with highest t value at 100ms in seated, closed loop,
% fully preprocessed recording
subIDs = cellstr(unique(meta_data{:,"sub"}));
max_chan = cell(1, length(subIDs));

for sub = 1:length(subIDs)
    raw_dat_names = meta_data{ismember(meta_data.sub, subIDs{sub}), "raw_data_name"};
    seated_closed_dat_name = raw_dat_names(logical(contains(raw_dat_names, 'seated').*contains(raw_dat_names, 'closed')));
    recording = contains(meta_data{:, "raw_data_name"}, seated_closed_dat_name);

    cd(meta_data{recording, "analysed_data_loc"});
    D = spm_eeg_load(['e_mffft_', char(strrep(seated_closed_dat_name, '.lvm', '.mat'))]);

    good_trials = indtrial(D, 'tone', 'GOOD');
    se = std(D(indchantype(D, 'MEGMAG', 'GOOD'),:,good_trials),[],3)./sqrt(length(good_trials));
    t = mean(D(indchantype(D, 'MEGMAG', 'GOOD'),:,good_trials),3)./se;
        
    [~, tind] = min(abs(D.time - 99*1e-3));
    max_colour = interp1([min(abs(t(:,tind))), max(abs(t(:,tind)))], [0.9, 0], abs(t(:,tind)));
    [~, plot_order] = sort(max_colour, 'descend');
    max_colour = repmat(max_colour, 1, 3);

    max_chan_tmp = plot_order(end);
    meg_chans = D.chanlabels(indchantype(D, 'MEGMAG', 'GOOD'));
    max_chan{sub} = meg_chans{max_chan_tmp};
end

for recording = 1:size(meta_data,1)
    
    cd(meta_data{recording, "analysed_data_loc"});

    % Epoch unfiltered/non-processed data
    D = spm_eeg_load(['t_', char(meta_data{recording,"raw_data_name"})]);

    if isfile(['e_', D.fname])
        D = spm_eeg_load(['e_', D.fname]);
    else
        S = [];
        S.D = D;
        S.timewin = [-200 500];
        S.condLabels = {'tone'};
        if strcmp(meta_data{recording, 'sub'}, 'sub-003')
            S.triggerChannels = {'AI16'};
        elseif strcmp(meta_data{recording, 'sub'}, 'sub-004')
            S.triggerChannels = {'AI8'};
        end
        D = spm_opm_epoch_trigger(S);

        % Set all epochs after 500 to bad
        goodTrials = indtrial(D, 'tone', 'GOOD');
        D = badtrials(D, goodTrials(501:end), 1);
        save(D);
    end

    % Cluster
    [~,sinds] = spm_match_str(D.chanlabels(indchantype(D, 'MEGMAG', 'GOOD')), D.sensors('MEG').label);
    goodTrials = indtrial(D, 'tone', 'GOOD');
    data = zeros(3, length(goodTrials));
    for trial = 1:length(goodTrials)
        [~, fiftyms_time] = min(abs(D.time - 50e-3));
        [~, onefiftyms_time] = min(abs(D.time - 150e-3));
        B_est = pinv(D.sensors('MEG').coilori(sinds,:))*D(indchantype(D, 'MEGMAG', 'GOOD'), [fiftyms_time, onefiftyms_time], goodTrials(trial));
        data(:, trial) = diff(B_est, 1, 2);
    end

    % Select trials with largest increase in Bz
    [~, idx_inc] = maxk(data(3,:), 0.1*size(data,2));

    % Select trials with largest decrease in Bz
    [~, idx_dec] = mink(data(3,:), 0.1*size(data,2));

    % Select 10% of trials with least change in Bz
    [~, idx_mid] = mink(abs(data(3,:)), 0.1*size(data,2));

    % Get epoched, analysed data for recording
    start_string = {'e_ffft_', 'e_hffft_', 'e_h2ffft_', 'e_mffft_'};
    for pp = 1:4
        DD{pp} = spm_eeg_load(char(fullfile(meta_data{recording, "analysed_data_loc"}, ...
            strcat(start_string{pp}, extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'), '.mat'))));
    end

    sub = strcmp(subIDs, meta_data{recording, 'sub'});
    for pp = 1:length(DD)
        D = DD{pp};
        figure; hold on;
        plot(D.time, mean(D(indchannel(D, max_chan{sub}),:,idx_mid),3));
        plot(D.time, mean(D(indchannel(D, max_chan{sub}),:,idx_inc),3));
        plot(D.time, mean(D(indchannel(D, max_chan{sub}),:,idx_dec),3));
        legend('No change', 'Increase', 'Decrease')

        figure; hold on;
        SE = std(D(indchannel(D, max_chan{sub}),:,idx_mid), [], 3);
        plot(D.time, mean(D(indchannel(D, max_chan{sub}),:,idx_mid)./sqrt(SE),3));
        SE = std(D(indchannel(D, max_chan{sub}),:,idx_inc), [], 3);
        plot(D.time, mean(D(indchannel(D, max_chan{sub}),:,idx_inc)./sqrt(SE),3));
        SE = std(D(indchannel(D, max_chan{sub}),:,idx_dec), [], 3);
        plot(D.time, mean(D(indchannel(D, max_chan{sub}),:,idx_dec)./sqrt(SE),3));
        legend('No change', 'Increase', 'Decrease')
    end
end

%% Get head position and orientation in room space

for recording = 1:size(meta_data,1)

    load(fullfile(meta_data{recording, "analysed_data_loc"}, ...
        strcat('opti_data_t_', extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'), '.mat')));
    
    marker_slots = readtable(fullfile(extractBefore(rawDataPath, 'rawData'), 'optitrack_marker_slots.csv'));
    table_of_info = fullfile(extractBefore(rawDataPath, meta_data{recording, 'sub'}), meta_data{recording, 'sub'}, ...
        'scanner-cast', 'table_of_info.csv');
    clear R T

    % Set marker heights
    marker_height_cm = zeros(height(marker_slots),1);
    marker_height_cm(strcmp(marker_slots.height, "tall")) = 56.19;
    marker_height_cm(strcmp(marker_slots.height, "mid")) = 49.69;
    marker_height_cm(strcmp(marker_slots.height, "short")) = 40.73;
    
    % Get marker positions in MRI coordinates
    MarkerPosMRI = getMarkerPosInMRIcoords(table_of_info,...
        marker_slots.slot, transpose(marker_height_cm), opti_data, 'Scannercast');
    
    % Get transformation matrices to go from MRI to room coordinates
    [~, ~, R, T] = getMagPosOriOverTime(opti_data, MarkerPosMRI, D, 'Scannercast');
    
    % Find head center over time
    D = spm_eeg_load(char(fullfile(meta_data{recording, "analysed_data_loc"}, ...
            strcat('e_ffft_', extractBefore(meta_data{recording, "raw_data_name"}, '.lvm'), '.mat'))));
    ctx = gifti(D.inv{1}.mesh.tess_ctx);
    cortex_center_MRI = transpose(mean(ctx.vertices,1));
    head_orientation_MRI = [0 1 0]';
    clear ctx
    cortex_center_Room = zeros(3, size(T,2));
    head_orientation_Room = zeros(3, size(T,2));
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
    for gap = 1:length(step_out)
        good_opt_data(step_out(gap):step_out(gap)+26) = 0;
    end
    step_out = step_out + 26;

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
    plot3(interp_cortex_center(good_opt_data==1,1), ...
        interp_cortex_center(good_opt_data==1,2), ...
        interp_cortex_center(good_opt_data==1,3), '.', 'color', [0.2666, 0.4471, 0.7686]);
    if any(good_opt_data == 0)
        plot3(interp_cortex_center(good_opt_data==0,1), ...
            interp_cortex_center(good_opt_data==0,2), ...
            interp_cortex_center(good_opt_data==0,3), ...
            '.', 'color', [0.5647, 0.6706, 0.8627], 'MarkerSize', 0.5);
    end
    good_inds = find(good_opt_data==1);
    quiver3(interp_cortex_center(good_inds(1:1000:end),1), ...
        repmat(max(interp_cortex_center(good_inds(1:1000:end),2)), length(good_inds(1:1000:end)), 1), ...
        interp_cortex_center(good_inds(1:1000:end),3), ...
        head_orientation_Room(1, good_inds(1:1000:end))', ...
        head_orientation_Room(2, good_inds(1:1000:end))', ...
        head_orientation_Room(3, good_inds(1:1000:end))', 2, 'color', [0.1804, 0.1725, 0.1843])

    zlim([-2 2]*1e3);
    ylim([0, 3.1]*1e3);
    xlim([-1.5, 1.5]*1e3);
    view(180,0);
    daspect([1 1 1]);
    set(gcf, 'Position', [680   670   497   308]);
    set(gca, 'FontSize', 14);
    xlabel('X, left-right (mm)')
    ylabel('Y, up-down (mm)')
    zlabel('Z, forward-back (mm)');
    set(gcf, 'color', 'w'); 

    print(fullfile(meta_data{recording, "results_save_loc"}, ...
        sprintf('%s_trajectory', extractBefore(meta_data{recording, "raw_data_name"}, ' Array 1.lvm'))),'-dpng','-r300');

end

