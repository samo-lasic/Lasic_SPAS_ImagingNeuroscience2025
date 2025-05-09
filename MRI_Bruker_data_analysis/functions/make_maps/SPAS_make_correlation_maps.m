function SPAS_make_correlation_maps(root_data_path, merged_names, roi_names, correlate_maps_struct, opt)

show_rois = 0;
if isfield(opt.fig,'show_rois')
    if opt.fig.show_rois & ~isempty(roi_names)
        show_rois = 1;
    end
end

% remove the trailing underscore
if strcmp(opt.map_name_prepend(end),'_') opt.map_name_prepend = opt.map_name_prepend(1:end-1); end

merged_names = strcat(merged_names,opt.map_name_prepend);

% filter files based on root_data_path, map_name_prepend and merged_names
all_files = dir(root_data_path);

ind = contains({all_files.name}, merged_names);
all_files = all_files(ind);
% there could be non filtered files remaining, e.g. s_a_g...

% find common zoom ranges
[zoom_x, zoom_y] = find_common_zoom_factor(opt.data_zoom_thresh, all_files, merged_names, correlate_maps_struct);


for n_data = 1:numel(merged_names)
    merged_name = merged_names{n_data};

    % filter files based on merged_name
    files = all_files(contains({all_files.name}, merged_name));


    n_maps = numel(correlate_maps_struct);
    for n_map = 1:n_maps

        % get the current maps_struct
        maps_opt = correlate_maps_struct(n_map);

        % adjust color_order for modulation
        if maps_opt.color_order(3) > 0
            maps_opt.color_order(3) = find(~ismember([1 2 3],[maps_opt.color_order(1) maps_opt.color_order(2)]));
        end


        % find all files needed
        [correlate_files, correlate_names] = SPAS_find_files(files, merged_name, maps_opt.correlate_prepend, maps_opt.correlate_name, maps_opt.correlate_append);
        [modulate_files, modulate_names] = SPAS_find_files(files, merged_name, maps_opt.modulate_prepend, maps_opt.modulate_name, maps_opt.modulate_append);

        do_modulate = ~isempty(modulate_files);

        if numel(correlate_files) ~= 2
            display('missing files to correlate')
            continue
        end

        correlate_nii_fn = fullfile({correlate_files.folder},{correlate_files.name});
        if do_modulate
            modulate_nii_fn = fullfile({modulate_files.folder},{modulate_files.name});
        end

        % maps to correlate
        [correlate_nii1, h_nii] = mdm_nii_read(correlate_nii_fn{1});
        correlate_nii2 = mdm_nii_read(correlate_nii_fn{2});

        % post mask
        if isfield(opt,'post_mask')
            post_mask = mdm_nii_read(fullfile(root_data_path,opt.post_mask));
            correlate_nii1 = correlate_nii1 .* double(post_mask);
            correlate_nii2 = correlate_nii2 .* double(post_mask);
        end

        size_nii = size(correlate_nii1);
        if ~isequal(size(correlate_nii1),size(correlate_nii2))
            display('size of images to correlate are different')
            continue
        end

        % get limits to be mapped to range [0 1] for color coding
        % & clamp images
        % & map to range [0 1]

        [correlate_LIM1, correlate_nii1_clamp] = get_limits_and_normalize_nii(...
            maps_opt.correlate_range.relative, maps_opt.correlate_range.LIM1, correlate_nii1);

        [correlate_LIM2, correlate_nii2_clamp] = get_limits_and_normalize_nii(...
            maps_opt.correlate_range.relative, maps_opt.correlate_range.LIM2, correlate_nii2);

        display(sprintf('source range for correlation map 1: %g : %g', correlate_LIM1(1), correlate_LIM1(2)));
        display(sprintf('source range for correlation map 2: %g : %g', correlate_LIM2(1), correlate_LIM2(2)));


        % gamma "correction" (only if type is not 0)
        correlate_nii1_col = img_gamma_correction(correlate_nii1_clamp, maps_opt.correlate_gamma_cor);
        correlate_nii2_col = img_gamma_correction(correlate_nii2_clamp, maps_opt.correlate_gamma_cor);


        % maps to modulate
        if do_modulate
            N = numel(modulate_nii_fn);
            modulate_nii = mdm_nii_read(modulate_nii_fn{1});

            % post mask
            if isfield(opt,'post_mask')
                modulate_nii = modulate_nii .* double(post_mask);
            end

            if ~isequal(size_nii,size(modulate_nii))
                display('sizes of images to modulate and correlate are different')
                continue
            end

            for n = 2:N
                modulate_nii = modulate_nii + mdm_nii_read(modulate_nii_fn{n});
            end

            modulate_nii = modulate_nii/N;

            % get limits, clamp, map to range [0 1]
            [modulate_LIM, modulate_nii_clamp] = get_limits_and_normalize_nii(...
                maps_opt.modulate_range.relative, maps_opt.modulate_range.LIM, modulate_nii);

            display(sprintf('source range for modulation map: %g : %g', modulate_LIM(1), modulate_LIM(2)));

            % gamma "correction"
            modulate_nii_col = img_gamma_correction(modulate_nii_clamp, maps_opt.modulate_gamma_cor);

        end



        n_slices = size_nii(3);

        if isempty(opt.fig.LimitSliceRange)
            slice_ind = 1:n_slices;
        else
            slice_ind = opt.fig.LimitSliceRange(1):opt.fig.LimitSliceRange(2);
        end


        if opt.save_bin_rois
            for n_roi = 1:numel(opt.bins)

                bin = opt.bins(n_roi);
                bin_mask = zeros(size_nii);


                mask_ind = ...
                    correlate_nii1 > max(correlate_nii1(:)) * bin.range1(1) & correlate_nii1 < max(correlate_nii1(:)) * bin.range1(2) & ...
                    correlate_nii2 > max(correlate_nii2(:)) * bin.range2(1) & correlate_nii2 < max(correlate_nii2(:)) * bin.range2(2);

                if ~isempty(bin.range3)
                    mask_ind = mask_ind & modulate_nii > max(modulate_nii(:)) * bin.range3(1) & modulate_nii < max(modulate_nii(:)) * bin.range3(2);
                end

                bin_mask(mask_ind) = 1;

                % fill mask holes
                bin_mask = mio_mask_fill(bin_mask);


                roi_bin_fn = strcat('roi_bin_', num2str(n_roi), '_', merged_name, '.nii.gz');
                roi_bin_fn = fullfile(root_data_path, roi_bin_fn);
                display(sprintf('saving bin mask: %s', roi_bin_fn))

                mdm_nii_write(single(bin_mask),roi_bin_fn, h_nii);

            end
        end


        [opt.fig, fh] = SPAS_map_fig_opt(n_map + n_maps*n_data, size_nii .* h_nii.pixdim(2:4)', zoom_x, zoom_y, opt.fig);

        for n = 0:length(slice_ind)-1

            n_col = rem(n,opt.fig.n_cols);
            n_row = floor(n/opt.fig.n_cols)+1;

            ah = axes('position',...
                [n_col * opt.fig.sub_fig_width ...
                (1 - opt.fig.title_size)*(1-n_row/opt.fig.n_rows) ...
                opt.fig.sub_fig_width opt.fig.sub_fig_height]);


            % make RGB map
            I1 = correlate_nii1_col(zoom_x, zoom_y, slice_ind(n+1))';
            I2 = correlate_nii2_col(zoom_x, zoom_y, slice_ind(n+1))';
            range1 = [min(I1(:)) max(I1(:))];
            range2 = [min(I2(:)) max(I2(:))];

            Icor = I1;
            Icor(:,:,maps_opt.color_order(1)) = Icor;
            Icor(:,:,maps_opt.color_order(2)) = I2;

            mod_ind = find(~ismember([1 2 3],[maps_opt.color_order(1) maps_opt.color_order(2)]));
            Icor(:,:,mod_ind) = 0*Icor(:,:,maps_opt.color_order(1));

            hI = imagesc(Icor);


            % modulation color or alpha mask
            if do_modulate
                Imod = modulate_nii_col(zoom_x, zoom_y, slice_ind(n+1))';
                if maps_opt.color_order(3) > 0 % color
                    Icor(:,:,maps_opt.color_order(3)) = Imod;
                    hI = imagesc(Icor);
                else % mask
                    set(hI, 'AlphaData', Imod);
                end
            end

            % roi
            if show_rois
                for n_roi = 1:numel(roi_names)
                    roi_name = roi_names{n_roi};
                    roi_fn = fullfile(root_data_path,roi_name);
                    roi = mdm_nii_read(roi_fn);
                    if n == 0
                        display(sprintf('showing roi: %s', roi_fn))
                    end


                    ROI = squeeze(roi(zoom_x,zoom_y,slice_ind(n+1)))';

                    try
                        plot_roi(ROI, opt.fig.roi_colors(n_roi,:),opt.fig.roi_lw);
                    catch
                        plot_roi(ROI);
                    end


                end
            end




            if opt.fig.reverseX == 1
                set(ah,'XDir','reverse')
            else
                set(ah,'XDir','normal')
            end

            if opt.fig.reverseY == 1
                set(ah,'YDir','reverse')
            else
                set(ah,'YDir','normal')
            end

            % ----------- make scatter plot (temporary)
            if (0)

                I1 = correlate_nii1(zoom_x, zoom_y, slice_ind(n+1));
                I2 = correlate_nii2(zoom_x, zoom_y, slice_ind(n+1));
                I1 = I1(:);
                I2 = I2(:);

                if do_modulate
                    I3 = modulate_nii(zoom_x, zoom_y, n+1);
                    I3 = I3(:);
                else
                    I3 = ones(size(I1));
                end

                %colormap('parula')
                colormap('turbo')

                hI = scatter(I1, I2, 30, I3, 'filled');
                set(hI, 'MarkerFaceAlpha', 'flat', 'AlphaData', 0.01 + 0.3*I3, ...
                    'AlphaDataMapping', 'none', 'SizeData', 1 + 30*I3);

            end

            axis(ah,'tight','off')
            pbaspect(ah,[1 opt.fig.im_ratio 1] .* h_nii.pixdim(2:4)')

        end

        % make title string
        col_str = {'alpha', 'red', 'green', 'blue'};

        if opt.fig.showTitle
            % title: cor. correlate_append from correlate_name{1} and correlate_name{2}, mod. by modulate_append from modulate_name{1} and modulate_name{2} (alpha)
            if isempty(maps_opt.title_str)

                if numel(maps_opt.modulate_name) == 2
                    modulate_str = sprintf(', mod. by %s from %s and %s',...
                        maps_opt.modulate_append, maps_opt.modulate_name{1}, maps_opt.modulate_name{2});
                elseif numel(maps_opt.modulate_name) == 1
                    if strcmp(maps_opt.modulate_name,'')
                        modulate_str = sprintf(', mod. by %s',...
                            maps_opt.modulate_append);
                    else
                        modulate_str = sprintf(', mod. by %s from %s',...
                            maps_opt.modulate_append, maps_opt.modulate_name{1});
                    end
                end

                if numel(maps_opt.modulate_name) > 0
                    modulate_str = sprintf('%s(%s)', modulate_str, col_str{maps_opt.color_order(3)+1});
                else
                    modulate_str = '';
                end

                fig_title = sprintf('cor. %s from %s(%s) and %s(%s)%s',...
                    maps_opt.correlate_append, ...
                    maps_opt.correlate_name{1}, col_str{maps_opt.color_order(1)+1}, ...
                    maps_opt.correlate_name{2}, col_str{maps_opt.color_order(2)+1}, ...
                    modulate_str);


            else
                fig_title = maps_opt.title_str;
            end

            ah = axes('position',[0 1-opt.fig.title_size 1 opt.fig.title_size]);
            axis(ah,'tight','off');
            th = text(0.01, .5, fig_title, 'color', 'white', 'fontsize', opt.fig.title_fs);
            th.Interpreter = 'none';
        end


        if opt.fig.save

            % make fig name
            fig_name = make_fig_name('_cor_', root_data_path, merged_name, maps_opt);

            fig_path = fullfile(root_data_path,'fig');
            mkdir(fig_path);
            fig_path = fullfile(fig_path,fig_name);
            print(fh,fig_path,'-dpng',opt.fig.resolution);
            display(sprintf('saved %s', fig_path))
        end






        % ------------------------------------------------

        % make RGB map
        if opt.fig.show_colormap

            % remember color_order and don't repeat if no change
            if exist('color_order')
                show_legend = ~isequal(color_order, maps_opt.color_order);
            else
                show_legend = 1;
            end
            color_order = maps_opt.color_order;


            if show_legend
                N = 500;
                x = linspace(0,1,N);

                Icor = zeros(N,N,3);
                for i = 1:N
                    Icor(:,i,maps_opt.color_order(1)) = x;
                    Icor(i,:,maps_opt.color_order(2)) = x;
                end

                Imod = zeros(N,N,3);
                if maps_opt.color_order(3) > 0 % color
                    Imod(:,:,maps_opt.color_order(3)) = ndgrid(x,x);
                else % mask
                    Imod(:,:,1) = ndgrid(x,x);
                    Imod(:,:,2) = Imod(:,:,1);
                    Imod(:,:,3) = Imod(:,:,1);
                end


                % define legend sizes
                figure_size_mm = .5 * opt.fig.size_mm(1) * [1 .9]; %* [28/40 1];
                legend_size = .9; % space for mod map
                legend_left = 0;
                legend_bottom = .05;
                legend_xspace = 0;
                legend_ratio2 = .1;


                fh = figure(opt.fig.num + 100*n_map - 1);
                clf
                set(gcf, 'InvertHardCopy', 'off')
                set(fh,'Color','black', 'Units', 'inches', ...
                    'PaperPosition',[0 0 figure_size_mm / 25.4 ], 'PaperPositionMode', 'auto');

                screen_size = get(0,'screensize');
                fh.Position(3:4) = opt.fig.display_ratio*figure_size_mm .* screen_size(4);

                % show legend maps
                ah = axes('Position', [2*legend_left, legend_bottom, legend_size, legend_size]);
                imshow(imrotate(Icor,90));
                axis(ah,'off')
                
                ah = axes('Position', [legend_left + legend_size + legend_ratio2*legend_xspace, legend_bottom, legend_ratio2*legend_size, legend_size]);
                imshow(imrotate(Imod,180));
                set(ah, 'DataAspectRatio', [1/legend_ratio2 1 1]);


                if opt.fig.save

                    % make fig name
                    fig_name = make_fig_name('_legend_', root_data_path, merged_name, maps_opt, col_str);

                    fig_path = fullfile(root_data_path,'fig');
                    mkdir(fig_path);
                    fig_path = fullfile(fig_path,fig_name);
                    print(fh,fig_path,'-dpng',opt.fig.resolution);
                    display(sprintf('saved %s', fig_path))
                end

            end
        end



        % -------- make scatter plots -------------------
        if opt.fig.show_roi_scatter

            %get all rois and colors
            roi_data = [];

            for n_roi = 1:numel(roi_names)
                roi_data(n_roi).roi_name = roi_names{n_roi};
                roi_fn = fullfile(root_data_path,roi_data(n_roi).roi_name);
                roi = mdm_nii_read(roi_fn);
                try
                    roi_data(n_roi).color = opt.fig.roi_colors(n_roi,:);
                catch
                    roi_data(n_roi).color = [1 1 1];
                end

                try
                    roi_data(n_roi).alpha = opt.fig.roi_alpha(n_roi);
                catch
                    roi_data(n_roi).alpha = 0.7;
                end

                roi_ind = find(roi);
                roi_data(n_roi).n = numel(roi_ind);

                if opt.fig.remove_roi_outliers
                    roi_ind = roi_ind(correlate_nii1(roi_ind)>0 & correlate_nii2(roi_ind) > 0);
                end

                roi_data(n_roi).correlate_nii1 = correlate_nii1(roi_ind);
                roi_data(n_roi).correlate_nii2 = correlate_nii2(roi_ind);


                roi_data(n_roi).correlate_nii1_lim = [min(correlate_nii1(roi_ind)), max(correlate_nii1(roi_ind))];
                roi_data(n_roi).correlate_nii2_lim = [min(correlate_nii2(roi_ind)), max(correlate_nii2(roi_ind))];

                % per slice
                for n_sl = 1:size_nii(3)
                    I1 = correlate_nii1(:,:,n_sl);
                    I2 = correlate_nii2(:,:,n_sl);
                    roi_ind = find(roi(:,:,n_sl));
                    roi_data(n_roi).slice(n_sl).n = numel(roi_ind);
                    roi_data(n_roi).slice(n_sl).correlate_nii1 = I1(roi_ind);
                    roi_data(n_roi).slice(n_sl).correlate_nii2 = I2(roi_ind);
                end

            end


            % --------------
            % scatter plot options
            scatter_opt.size_mm = 60*[1 1];
            scatter_opt.display_ratio = 8e-5;
            scatter_opt.fs = 14;

            scatter_opt.marker_size = 14; %30;
            %             scatter_opt.marker_alpha = .75;
            %             scatter_opt.marker_alpha = scatter_opt.marker_alpha * linspace(1,1,numel(roi_data));
            scatter_opt.color = [0 0 0];

            contour_opt.lw = 3;
            contour_opt.contour_levels = round(15*opt.fig.roi_detail);
            contour_opt.histogram_bins = round(30*opt.fig.roi_detail);
            contour_opt.do_blure = 1;
            contour_opt.blure_sigma = opt.fig.roi_smoothing; %%2;

            % get global limits
            XLIM = [];
            YLIM = [];

            for n_roi = 1:numel(roi_data)
                XLIM(n_roi,:) = roi_data(n_roi).correlate_nii1_lim;
                YLIM(n_roi,:) = roi_data(n_roi).correlate_nii2_lim;
            end
            XLIM = [min(XLIM(:)) max(XLIM(:))];
            YLIM = [min(YLIM(:)) max(YLIM(:))];


            % plot roi data per slice
            opt.fig.color = scatter_opt.color; %'black'; %'white'
            [opt.fig, fh] = SPAS_map_fig_opt(n_map + n_maps*n_data + 2, size_nii, zoom_x, zoom_y, opt.fig);

            sub_fig_width = opt.fig.sub_fig_width * .9;
            sub_fig_height = opt.fig.sub_fig_height * .9;
            title_size = 0;
            fig_offset_x = .05;
            fig_offset_y = .08;

            for n = 0:length(slice_ind)-1 
                n_col = rem(n,opt.fig.n_cols);
                n_row = floor(n/opt.fig.n_cols)+1;

                ah = axes('position',...
                    [fig_offset_x + n_col * opt.fig.sub_fig_width ...
                    fig_offset_y + (1 - title_size)*(1-n_row/opt.fig.n_rows) ...
                    sub_fig_width sub_fig_height]);

                tmp_c = .7*[1 1 1] - scatter_opt.color;
                ah.Color = scatter_opt.color;
                ah.XColor = tmp_c;
                ah.YColor = tmp_c;


                % loop roi
                hold on

                for m_roi = 1:numel(roi_data)
                    n_roi = opt.fig.roi_plot_order(m_roi);

                    col = roi_data(n_roi).color;
                    alpha = roi_data(n_roi).alpha;
                    I1 = roi_data(n_roi).slice(slice_ind(n+1)).correlate_nii1;
                    I2 = roi_data(n_roi).slice(slice_ind(n+1)).correlate_nii2;

                    if opt.fig.show_roi_contours
                        c = plot_contours(I1, I2, XLIM, YLIM, col, alpha, contour_opt);
                    end

                    hI = scatter(I1, I2, scatter_opt.marker_size, 'filled');
                    set(hI, 'MarkerFaceColor', col, 'MarkerFaceAlpha', alpha);

                end

                xlim(XLIM)
                ylim(YLIM)
                pbaspect(ah,[1 1 1])
                %axis(ah,'tight')
                set(ah, 'LineWidth',2, 'TickDir', 'out')

                %pbaspect(ah,[opt.fig.im_ratio 1 1])
            end

            if opt.fig.save
                % make fig name
                fig_name = make_fig_name('_scatter_sl_', root_data_path, merged_name, maps_opt);

                fig_path = fullfile(root_data_path,'fig');
                mkdir(fig_path);
                fig_path = fullfile(fig_path,fig_name);
                print(fh,fig_path,'-dpng',opt.fig.resolution);
                display(sprintf('saved %s', fig_path))
            end


            % plot all roi data

            tmp_opt_fig = opt.fig;
            opt.fig.size_mm = scatter_opt.size_mm;
            opt.fig.display_ratio = scatter_opt.display_ratio;
            opt.fig.color = scatter_opt.color; %'black'; %'white'
            [opt.fig, fh] = fig_opt(n_map + n_maps*n_data + 1, opt.fig);
            %             opt.fig = tmp_opt_fig;

            % override limits
              if isfield(opt.fig,'XLIM')
                if ~isempty(opt.fig.XLIM)
                    XLIM = opt.fig.XLIM;
                end
            end
            if isfield(opt.fig,'YLIM')
                if ~isempty(opt.fig.YLIM)
                    YLIM = opt.fig.YLIM;
                end
            end
            display(sprintf('XLIM = [%g, %g]',XLIM(1),XLIM(2)))
            display(sprintf('YLIM = [%g, %g]',YLIM(1),YLIM(2)))

            hold on

            for m_roi = 1:numel(roi_data)
                n_roi = opt.fig.roi_plot_order(m_roi);
                col = roi_data(n_roi).color;
                alpha = roi_data(n_roi).alpha;
                I1 = roi_data(n_roi).correlate_nii1;
                I2 = roi_data(n_roi).correlate_nii2;

                if opt.fig.show_roi_contours
                    c = plot_contours(I1, I2, XLIM, YLIM, col, alpha, contour_opt);
                end

                hI = scatter(I1, I2, scatter_opt.marker_size, 'filled');
                set(hI, 'MarkerFaceColor', col, 'MarkerFaceAlpha', alpha);


            end

            xlim(XLIM)
            ylim(YLIM)
            pbaspect(gca,[1 1 1])
            set(gca, 'LineWidth',2, 'TickDir', 'out','FontSize',scatter_opt.fs)
            tmp_c = .7*[1 1 1] - scatter_opt.color;
            set(gca, 'Color',scatter_opt.color,'XColor',tmp_c,'YColor',tmp_c)

            if opt.fig.save
                % make fig name
                fig_name = make_fig_name('_scatter_', root_data_path, merged_name, maps_opt);

                fig_path = fullfile(root_data_path,'fig');
                mkdir(fig_path);
                fig_path = fullfile(fig_path,fig_name);
                print(fh,fig_path,'-dpng',opt.fig.resolution);
                display(sprintf('saved %s', fig_path))
            end


            % show legend
            [opt.fig, fh] = fig_opt(n_map + n_maps*n_data + 2, opt.fig);

            hold on
            leg_str = {};
            for m_roi = 1:numel(roi_data)
                n_roi = opt.fig.roi_plot_order(m_roi);
                col = roi_data(n_roi).color;
                alpha = roi_data(n_roi).alpha;
                leg_str{end+1} = roi_data(n_roi).roi_name;
                plot([0 1],[1 1]*n_roi,'-','LineWidth',10,'color',[col alpha])
            end
            pbaspect(gca,[1 .5 1])
            set(gca, 'Color', scatter_opt.color)



            lh = legend(leg_str,'interpreter','none','box','off',...
                'location','northeastoutside','TextColor',[1 1 1]);

            if opt.fig.save
                % make fig name
                fig_name = make_fig_name('_legend_', root_data_path, merged_name, maps_opt);

                fig_path = fullfile(root_data_path,'fig');
                mkdir(fig_path);
                fig_path = fullfile(fig_path,fig_name);
                print(fh,fig_path,'-dpng',opt.fig.resolution);
                display(sprintf('saved %s', fig_path))
            end

        end
    end
end
end






function [zoom_x, zoom_y] = find_common_zoom_factor(zoom_thresh, all_files, merged_names, correlate_maps_struct)
% find common zoom factor (pixel range not intensity!) across all maps
Ix = 0; Iy = 0;
cnt = 0;

for n_data = 1:numel(merged_names)
    merged_name = merged_names{n_data};

    % filter files
    files = all_files(contains({all_files.name}, merged_name));

    n_maps = numel(correlate_maps_struct);
    for n_map = 1:n_maps

        % find all files needed
        maps_opt = correlate_maps_struct(n_map);

        correlate_files = SPAS_find_files(files, merged_name, maps_opt.correlate_prepend, maps_opt.correlate_name, maps_opt.correlate_append);
        modulate_files = SPAS_find_files(files, merged_name, maps_opt.modulate_prepend, maps_opt.modulate_name, maps_opt.modulate_append);

        if ~isempty(correlate_files)
            nii_fns = fullfile({correlate_files.folder},{correlate_files.name});
            for n = 1:numel(nii_fns)
                nii_fn = nii_fns{n};
                display(sprintf('%s',nii_fn))

                nii = abs(mdm_nii_read(nii_fn));
                I = sum(nii,3);
                Ix = Ix + sum(I,2); Ix = Ix/max(Ix);
                Iy = Iy + sum(I,1); Iy = Iy/max(Iy);
                cnt = cnt + 1;
            end
        end

        if ~isempty(modulate_files)
            nii_fns = fullfile({modulate_files.folder},{modulate_files.name});
            for n = 1:numel(nii_fns)
                nii_fn = nii_fns{n};
                display(sprintf('%s',nii_fn))

                nii = abs(mdm_nii_read(nii_fn));
                I = sum(nii,3);
                Ix = Ix + sum(I,2); Ix = Ix/max(Ix);
                Iy = Iy + sum(I,1); Iy = Iy/max(Iy);
                cnt = cnt + 1;
            end
        end

    end
end

if cnt > 0
    Ix = Ix/cnt;
    Iy = Iy/cnt;

    zoom_x = find(Ix >= zoom_thresh);
    zoom_y = find(Iy >= zoom_thresh);
else
    zoom_x = [];
    zoom_y = [];

end
end


function [sourceLIM, nii] = get_limits_and_normalize_nii(relative, LIM, nii)
% get limits to be mapped to range [0 1] for color coding
% & clamp images
% & map to range [0 1]

sourceLIM = [min(nii(:)) max(nii(:))];

if isempty(LIM)
    clampLIM = sourceLIM;
else
    if relative
        clampLIM = sourceLIM .* LIM;
    else
        clampLIM = LIM;
    end

    if isnan(LIM(1))
        clampLIM(1) = sourceLIM(1);
    end
    if isnan(LIM(2))
        clampLIM(2) = sourceLIM(2);
    end
end

% clamp
nii(nii < clampLIM(1)) = clampLIM(1);
nii(nii > clampLIM(2)) = clampLIM(2);

% map to range [0 1]
nii = (nii-clampLIM(1)) / (clampLIM(2)-clampLIM(1));
end

function fig_name = make_fig_name(fig_base_name, root_data_path, merged_name, maps_opt, col_str)

[~, fig_name] = fileparts(root_data_path);
fig_name = strcat(fig_name,'_',merged_name, fig_base_name, ...
    maps_opt.correlate_name{1},'_',maps_opt.correlate_name{2}, '_', maps_opt.correlate_append);

if numel(maps_opt.modulate_name) > 0
    fig_name = [fig_name '_' maps_opt.modulate_append];
end

for n = 1:numel(maps_opt.color_order)
    if ( n < 3 || numel(maps_opt.modulate_name) > 0 ) & (nargin == 5)
        fig_name = [fig_name '_' col_str{maps_opt.color_order(n)+1}];
    end
end

end

function c = plot_contours(I1, I2, XLIM, YLIM, col, alpha, contour_opt)

% contour_opt.histogram_bins = 30; %round(length(I1)/10);
%grid1 = linspace(0.95*min(I1(:)), 1.05*max(I1(:)), contour_opt.histogram_bins);
%grid2 = linspace(0.95*min(I2(:)), 1.05*max(I2(:)), contour_opt.histogram_bins);
%                 grid1 = linspace(min(I1), max(I1), HistogramBins);
%                 grid2 = linspace(min(I2), max(I2), HistogramBins);
grid1 = linspace(XLIM(1), XLIM(2), contour_opt.histogram_bins);
grid2 = linspace(YLIM(1), YLIM(2), contour_opt.histogram_bins);

density = histcounts2(I1(:), I2(:), grid1, grid2);

if (contour_opt.do_blure)
    % Gaussian filter
    [xG, yG] = meshgrid(-5:5);
    g = exp(-xG.^2./(2.*contour_opt.blure_sigma.^2)-yG.^2./(2.*contour_opt.blure_sigma.^2));
    g = g./sum(g(:));
    density = conv2(density, g, 'same');
end

N = length(density);
grid1 = linspace(min(grid1), max(grid1), N);
grid2 = linspace(min(grid2), max(grid2), N);

% interpolate
grid1_int = linspace(XLIM(1), XLIM(2), 5*contour_opt.histogram_bins);
grid2_int = linspace(YLIM(1), YLIM(2), 5*contour_opt.histogram_bins);

[mesh1, mesh2] = meshgrid(grid1,grid2);
[mesh1_int, mesh2_int] = meshgrid(grid1_int,grid2_int);

density = interp2(mesh1, mesh2, density, mesh1_int, mesh2_int, 'spline');

grid1 = grid1_int;
grid2 = grid2_int;

density = density/max(density(:));

levels = linspace(0.1,0.99,contour_opt.contour_levels);

for l = levels
    [M,c] = contour(grid1, grid2, density',l*[1 1]);
    set(c,'EdgeColor',col,'EdgeAlpha',l*alpha, 'LineWidth',contour_opt.lw)
end


end