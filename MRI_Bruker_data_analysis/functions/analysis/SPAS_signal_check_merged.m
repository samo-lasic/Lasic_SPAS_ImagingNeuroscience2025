function [] = SPAS_signal_check_merged(fig_n, nii_fn, roi_fns, opt)
% show ROI-average signal attenuation (not normalized without ADC info or normalized with ADC info)

% read nifti
nii = mdm_nii_read(nii_fn);
xps = mdm_xps_load(mdm_fn_nii2xps(nii_fn));

if ~isfield(xps,'s_ind') % fix for single (not merged waveforms)
    xps.s_ind = ones(xps.n,1);
end

if ~isfield(xps,'wfm_names') % fix for single (not merged waveforms)
    [~, name] = fileparts(nii_fn);
    xps.wfm_names = {extractBefore(name,'.nii')};
end

if ~isfield(opt,'change_order')
    opt.change_order = 0;
end

if ~isfield(opt,'show_geoSPAS')
    opt.show_geoSPAS = 1;
end

if ~isfield(opt,'show_separate_wfm_sig') % optionally show signals for each waveform separately
    opt.show_separate_wfm_sig = 0;
end

fs1 = 12;
fs2 = 18;

% folder and image name
[root_data_path, img_name]=fileparts(nii_fn);
img_name = extractBefore(img_name,'.nii');
[~, img_folder] = fileparts(root_data_path);


for n_roi = 1:numel(roi_fns)
    roi_fn = roi_fns{n_roi};

    % data and roi path and name
    [root_data_path, name]=fileparts(nii_fn);
    name = extractBefore(name,'.nii');
    [~,roi_name] = fileparts(roi_fn);
    roi_name = extractBefore(roi_name,'.nii');

    % roi signal
    roi = mdm_nii_read(roi_fn);
    sig = squeeze(mean(nii.*roi, [1,2,3]))';

    display(sprintf('roi%d - %s', n_roi, roi_fn))
    display(sprintf('voxels = %d', numel(find(roi))))


    if (0)
        figure(1), clf
        hold on
        plot(xps.rep_dim_ind,'k.')
        plot(xps.rot_ind,'r.')
        plot(xps.rep_a_ind,'g.')
        plot(xps.a_ind,'b.')
    end

    if (0)
        % test average repetitions -----------
        % read nifti
        nii_t_fn = strrep(nii_fn,'_a.nii.gz','.nii.gz');
        nii_t = mdm_nii_read(nii_t_fn);
        xps_t = mdm_xps_load(mdm_fn_nii2xps(nii_t_fn));

        % roi signal
        sig_t = squeeze(mean(mean(mean(nii_t.*roi))))';


        rep_a_ind = unique(xps_t.rep_a_ind);
        rep_dim_ind = xps_t.rep_dim_ind;
        s_ind = unique(xps_t.s_ind);
        sig_a = sig_t;

        if length(rep_a_ind) > 1
            ind = [];
            for n_wfm = 1:length(s_ind)
                for n_rep = 1:length(rep_a_ind)
                    ind_tmp = find(xps_t.rep_a_ind == rep_a_ind(n_rep) & xps_t.s_ind == s_ind(n_wfm))
                    ind = [ind; ind_tmp(1)];

                    sig_a(ind_tmp(1)) = mean(sig_a(ind_tmp));
                end
            end
        end

        % reduce to the indices containing mean values
        sig_a = sig_a(ind);
        figure(1),clf, hold on, plot(sig,'o'), plot(sig_a,'.')
        %---------------

    end

    % select waveforms 
    if isfield(opt,'select_wfm_names')
        [sig, xps] = SPAS_filter_sig_and_xps_wfm_names(sig, xps, opt.select_wfm_names);
    end


    % b x rep x dir x wfm
    [sig_Mat, b_Mat, col, ~] = sortData(xps,sig);

    Nb = size(sig_Mat,1);
    Nrep = size(sig_Mat,2);
    Ndir = size(sig_Mat,3);
    Nwfm = size(sig_Mat,4);

    marker = {'.','o','x','+','*','v','d','^','s','>','<'};
    ms = [16,8,8,10,8,10,10,10,12,10,10];
    lineStyle = {'-','--',':','-.'};


    % check signals
    if(0)
        figure(1000 + fig_n + n_roi - 1),clf
        set(gcf,'color','white')
        hold
        for nwfm = 1:Nwfm
            for ndir = 1:Ndir
                for nrep = 1:Nrep
                    bSI = squeeze(b_Mat(:,nrep,ndir,nwfm));
                    sig_tmp = squeeze(sig_Mat(:,nrep,ndir,nwfm));

                    plot(bSI,sig_tmp,'.',...
                        'Color',col(1+rem(nwfm-1,length(col)),:),'LineWidth',2, 'MarkerSize',16)

                end
            end
        end
        set(gca,'LineWidth',3,'Yscale','Log','Box','off','TickDir','out','TickLength',[.02 .02],'FontSize',fs1)
    end

    % -------- normalize and plot normalized decays (all in one plot)
    if opt.normalize && ~opt.normalize_to_first_b
        % fit options
        b_low_lim = 0*3e9; % if 0 use first two b-values
        if opt.weighted
            fit_opt = fit_muFA_weighted_opt(opt);
        else
            fit_opt = fit_muFA_opt(opt);
        end

    end


    fh1 = figure(fig_n + n_roi - 1);,clf
    set(gcf,'color','white')

    for nwfm = 1:Nwfm
        subplot(2,Nwfm,nwfm)
        hold on

        for nrep = 1:Nrep
            for ndir = 1:Ndir
                bSI = squeeze(b_Mat(:,nrep,ndir,nwfm));
                sig_tmp = double(squeeze(sig_Mat(:,nrep,ndir,nwfm)));

                if opt.normalize
                    if ~opt.normalize_to_first_b
                        % quick fit / check
                        b_low_ind = find(bSI < opt.b_low_lim);
                        if isempty(b_low_ind)
                            b_low_ind = [1 2]';
                        end
                        sig_low = sig_tmp(b_low_ind);
                        b_low = bSI(b_low_ind);

                        X = [ones(length(b_low),1) -b_low];
                        P = X\log(sig_low);

                        s0 = exp(P(1));
                        if opt.gamma_fit

                            if opt.weighted
                                Pout = fit_gamma1_weighted_par(bSI,sig_tmp/s0, fit_opt);
                                % sig_norm1 = sig1/S0/Pout(1);
                            else
                                Pout = fit_gamma1_par(bSI,sig_tmp/s0,fit_opt);
                                % sig_norm1 = sig1/S0/Pout(1);
                            end

                            sig_norm = sig_tmp/s0/Pout(1);

                            sig_norm_fit = fit_gamma1([1 Pout(2:3)],bSI);
                            fit_parameters(nrep,ndir,:) = [s0*Pout(1) Pout(2:3)];
                        else % no gamma fit
                            sig_norm = sig_tmp/s0;
                            sig_norm_fit = exp(P(1)-P(2)*bSI)/s0;
                            fit_parameters(nrep,ndir,:) = [s0 P(2) 0];
                        end

                    else % opt.normalize_to_first_b == 1
                        sig_norm = sig_tmp/sig_tmp(1);
                        %fit_parameters(nrep,ndir,:) = [s0*Pout(1) Pout(2:3)];
                    end

                else % normalize == 0
                    sig_norm = sig_tmp;
                end

                b = bSI*1e-6;
                plot(b,sig_norm,'LineStyle', lineStyle{1+rem(nrep-1,length(lineStyle))},...
                    'Marker',marker{1+rem(nrep-1,length(marker))},...
                    'Color',col(1+rem(ndir-1,length(col)),:),'LineWidth',2, 'MarkerSize',ms(nrep))

                sig_Mat_norm(:,nrep,ndir,nwfm) = sig_norm;

            end
        end

        Ymin = min(sig_Mat_norm(:));

        if opt.normalize
            Ymax = 1/1.05;
        else
            Ymax = max(sig_Mat_norm(:));
        end

        if opt.show_labels
            xlabel('b [s/mm^2]')
            ylabel('signal')
        end

        if ~isempty(opt.XLim)
            xlim(opt.XLim)
        else
            xlim([0 max(b)*1.05])
        end


        if ~isempty(opt.YLim)
            ylim(opt.YLim)
        else
            ylim([Ymin * .95 Ymax * 1.05])
        end

        set(gca,'LineWidth',3,'Yscale','Log','Box','off','TickDir','out','TickLength',[.02 .02],'FontSize',fs1)
        wfm_name = xps.wfm_names{nwfm};
        if contains(wfm_name,'ETE')
            ind = findstr(wfm_name,'_');
            wfm_name = wfm_name(1:ind(2)-1);
        else
            wfm_name = extractBefore(wfm_name,'_');
        end

        if opt.showTitle
            title(wfm_name,'interpreter','none')
        end

        if opt.normalize && ~opt.normalize_to_first_b
            Ds = fit_parameters(:,:,2);

            meanD = mean(Ds(:));
            stdD = std(Ds(:));

            mean_stdDrep = mean(std(Ds));

            subplot(2,Nwfm,Nwfm+nwfm)
            bar(Ds*1e9)
            xlabel('direction')
            ylabel('ADC [mm^2/ms]')

            set(gca,'LineWidth',3,'Yscale','Lin','Box','off','TickDir','out','TickLength',[.02 .02],'FontSize',fs1)
            title(sprintf('D = %.2f ± %.2f [%.2f] x 10^{-9}', meanD*1e9, stdD*1e9, mean_stdDrep*1e9),'FontSize',10)
        end

    end

    % -------- plot decays (for each waveform separately)

    if opt.show_separate_wfm_sig
        for nwfm = 1:Nwfm
            fh11 = figure(10 + fig_n + n_roi + nwfm - 2);,clf
            fh11.Position(3) = 336;
            set(gcf,'color','white')

            hold on


            for nrep = 1:Nrep
                for ndir = 1:Ndir
                    bSI = squeeze(b_Mat(:,nrep,ndir,nwfm));

                    sig_norm = double(squeeze(sig_Mat_norm(:,nrep,ndir,nwfm)));

                    b = bSI*1e-6;
                    plot(b,sig_norm,'LineStyle', lineStyle{1+rem(nrep-1,length(lineStyle))},...
                        'Marker',marker{1+rem(nrep-1,length(marker))},...
                        'Color',col(1+rem(ndir-1,length(col)),:),'LineWidth',2, 'MarkerSize',ms(nrep))

                end
            end

            if opt.show_labels
                xlabel('b [s/mm^2]')
                ylabel('signal')
            end

            if ~isempty(opt.XLim)
                xlim(opt.XLim)
            else
                xlim([0 max(b)*1.05])
            end

            if ~isempty(opt.YLim)
                ylim(opt.YLim)
            else
                ylim([Ymin * .95 Ymax * 1.05])
            end

            set(gca,'LineWidth',3,'Yscale','Log','Box','off','TickDir','out','TickLength',[.02 .02],'FontSize',fs2)
            wfm_name = xps.wfm_names{nwfm};
            if contains(wfm_name,'ETE')
                ind = findstr(wfm_name,'_');
                wfm_name = wfm_name(1:ind(2)-1);
            else
                wfm_name = extractBefore(wfm_name,'_');
            end

            if opt.showTitle
                title(wfm_name,'interpreter','none')
            end

            if opt.normalize
                fig_name = [roi_name '_' wfm_name '_norm_sig_vs_b' ];
            else
                fig_name = [roi_name '_' wfm_name '_sig_vs_b' ];
            end

            if opt.fig.save
                fig_path = fullfile(root_data_path,'fig');
                mkdir(fig_path);
                fig_path = fullfile(fig_path, fig_name);
                print(fh11, fig_path,'-dpng',opt.fig.resolution);
                display(sprintf('saved %s', fig_path))
            end

        end
    end


    % repetion & direction average (all waveforms)
    sig_norm_pa = squeeze(mean(mean(sig_Mat_norm,3),2));


    % geoSPAS
    SPAS_ind = find(contains(xps.wfm_names,'SPAS'));
    show_geoSPAS = opt.show_geoSPAS & numel(SPAS_ind) > 1;

    if show_geoSPAS
        Egeo = geomean(sig_norm_pa(:,SPAS_ind),2);
    end

    % optionally change order of SPAS
    if opt.change_order
        sig_norm_pa = sig_norm_pa(:,[1 3 2 4]);
    end

    b = uniquetol(b_Mat,1e-3)*1e-6;
    if length(bSI) ~= Nb
        display('check b-vals across dimensions!')
    else
        col = copper(Nwfm);

        for nwfm = 1:Nwfm
            % wfm names
            style_ind(nwfm) = find(arrayfun(@(n) contains(xps.wfm_names{nwfm}, opt.styles(n).wfm_name), 1:length(opt.styles)));
        end


        fh2 = figure(100 + fig_n + n_roi - 1);
        clf
        fh2.Position(3) = 336;
        set(gcf,'color','white')
        hold on
        for nwfm = 1:Nwfm

            style = opt.styles(style_ind(nwfm));

            ph = plot(b,sig_norm_pa(:,nwfm));
            set(ph,'LineStyle', style.line, 'Color',style.col,...
                'LineWidth',style.lw, 'Marker',style.marker, 'MarkerSize',style.ms)
        end

        if show_geoSPAS

            style = opt.styles(find(contains({opt.styles.wfm_name},'geoSPAS')));

            ph = plot(b,Egeo);
            set(ph,'LineStyle', style.line, 'Color',style.col,...
                'LineWidth',style.lw, 'Marker',style.marker, 'MarkerSize',style.ms)
        end

        if opt.show_labels
            xlabel('b [s/mm^2]')
            ylabel('signal')
        end

        if ~isempty(opt.XLim)
            xlim(opt.XLim)
        else
            xlim([0 max(b)*1.05])
        end

        if ~isempty(opt.YLim)
            ylim(opt.YLim)
        else
            Ymin = min(sig_norm_pa(:));

            if opt.normalize
                if show_geoSPAS
                    Ymin = min([Ymin; Egeo(:)]);
                end
                ylim([Ymin * .95 1])

            else
                Ymax = max(sig_norm_pa(:));
                % ylim([min([sig_norm_pa(:); Egeo(:)])*.95 max([sig_norm_pa(:); Egeo(:)])*1.05])

                if show_geoSPAS
                    Ymin = min([Ymin; Egeo(:)]);
                    Ymax = max([Ymax; Egeo(:)]);
                end
                ylim([Ymin * .95 Ymax * 1.05])

            end
        end

        set(gca,'LineWidth',3,'Yscale','Log','Box','off','TickDir','out','TickLength',[.02 .02],'FontSize',fs2)
        if opt.show_labels
            xlabel('b [mm^2/s]')
            ylabel('average normalized signal')
        end

        for m = 1:numel(xps.wfm_names)
            wfm_name = xps.wfm_names{m};
            if contains(wfm_name,'ETE')
                ind = findstr(wfm_name,'_');
                wfm_name = wfm_name(1:ind(2)-1);
            else
                wfm_name = extractBefore(wfm_name,'_');
            end
            legend_str{m} = wfm_name;
        end

        legend_str{end+1} = 'geoSPAS';

        % change legend string (the order can be set in the merging step)
        if opt.change_order
            legend_str = legend_str([1 3 2 4 5]);
        end

        lh = legend(legend_str,'interpreter','none');
        lh.Box = 'off';
        lh.ItemTokenSize(1) = 1.2*lh.ItemTokenSize(1);
        lh.Position(4) = 1.2*lh.Position(4);

        fh1.Position(1) = 200;
        fh2.Position(1) = fh1.Position(1) + fh1.Position(3);

        title_name = sprintf('%s/%s:%s', img_folder, img_name, roi_name);

        if opt.normalize
            fig_name = [roi_name '_norm_sig_vs_b'];
        else
            fig_name = [roi_name '_sig_vs_b'];
        end

        if opt.showTitle
            title(title_name,'interpreter','none','FontSize', opt.TitleFontSize);
        end

        if opt.fig.save
            fig_path = fullfile(root_data_path,'fig');
            mkdir(fig_path);
            fig_path = fullfile(fig_path, fig_name);
            print(fh2, fig_path,'-dpng',opt.fig.resolution);
            display(sprintf('saved %s', fig_path))
        end


        % re-fit STE and geoSPAS
        if opt.normalize & ~opt.normalize_to_first_b

            STE_ind = find(contains(xps.wfm_names,'STE'));
            if numel(STE_ind) == 1
                sig_tmp = sig_norm_pa(:,STE_ind);
                if opt.weighted
                    Pout = fit_gamma1_weighted_par(bSI,sig_tmp,fit_opt);
                else
                    Pout = fit_gamma1_par(bSI,sig_tmp,fit_opt);
                end
                display(sprintf('gamma fit - STE: D = %g , mu = %g', Pout(2), Pout(3)))

            end

            if show_geoSPAS
                sig_tmp = sig_norm_pa(:,end);
                if opt.weighted
                    Pout = fit_gamma1_weighted_par(bSI,sig_tmp,fit_opt);
                else
                    Pout = fit_gamma1_par(bSI,sig_tmp,fit_opt);
                end
                display(sprintf('gamma fit - geoSPAS: D = %g , mu = %g', Pout(2), Pout(3)))

            end
        end



        % check signal differences
        if (opt.check_sig_dif)
            STE_ind = find(contains(xps.wfm_names,'STE'));
            show_TDD = show_geoSPAS;
            show_uA = show_geoSPAS & numel(STE_ind) == 1;

            if show_TDD || show_uA
                fh3 = figure(1000 + fig_n + n_roi - 1);
                clf
                fh3.Position(3) = .6*fh3.Position(3);
                set(gcf,'color','white')
                %axes('Position',[0.1 0.1 0.45 0.7]);
                hold on
                if show_TDD
                    ds_TDD = max(sig_norm_pa(:,SPAS_ind)')' - min(sig_norm_pa(:,SPAS_ind)')';
                    plot(b,ds_TDD,'r-','LineWidth',2, 'MarkerSize',20)
                end

                if show_uA
                    ds_uA = Egeo - sig_norm_pa(:,STE_ind);
                    plot(b,ds_uA,'b-','LineWidth',2, 'MarkerSize',20)
                end

                if ~isempty(opt.XLim)
                    xlim(opt.XLim)
                else
                    xlim([0 max(b)*1.05])
                end

                set(gca,'LineWidth',3,'Yscale','Lin','Box','off','TickDir','out','TickLength',[.02 .02],'FontSize',fs2)
                if opt.show_labels
                    xlabel('b [mm^2/s]')
                    ylabel('signal differences')
                end
                lh = legend({'TDD','uA'});
                lh.Location = 'northwest';
                lh.Box = 'off';
                lh.ItemTokenSize(1) = 1.2*lh.ItemTokenSize(1);
                lh.Position(4) = 1.2*lh.Position(4);

                [~,name]=fileparts(nii_fn);
                name = extractBefore(name,'.nii');
                title(name,'interpreter','none');
            end
        end

    end

end

% --------------------------- FUNCTIONS

function [sig_Mat, b_Mat, col, rot_Mat] = sortData(xps,sig)
% color antipodal directions

% eliminate b = 0
b_ind = find(xps.b ~= 0);
uwfm = unique(xps.s_ind(b_ind));
urot = unique(xps.rot_ind(b_ind));

col = copper(length(urot));

urep_dim = unique(xps.rep_dim_ind(b_ind));

Nwfm = length(uwfm);
Nrot = length(urot);

% find number of repetions
Nrep = length(urep_dim);


for n_wfm = 1:length(uwfm)
    for n_rot = 1:length(urot)
        for n_rep = 1:Nrep
            ind = find(xps.s_ind == uwfm(n_wfm) & xps.rot_ind == urot(n_rot) & xps.rep_dim_ind == urep_dim(n_rep) );

            if ~isempty(ind)
                sig_Mat(:,n_rep,n_rot,n_wfm) = sig(ind);
                b_Mat(:,n_rep,n_rot,n_wfm) = xps.b(ind);
                rot_Mat(:,n_rep,n_rot,n_wfm) = xps.rot_ind(ind);
            end
        end
    end
end





