function SPAS_make_parameter_maps(root_data_path, merged_names, roi_names, opt)
% show / save parameter maps and subtraction maps
% subtraction maps compare maps from the same session, can be used to monitor changes in vivo
% using e.g. merged_104_a_g_dti_md_geoSPAS
% merged_104_a_g_pa_raw_TDD_ds

% take a subset of selected maps
n_maps = numel(opt.select_map_numbers);
opt.map_names = opt.map_names(opt.select_map_numbers);
opt.title_strings = opt.title_strings(opt.select_map_numbers);

opt.autoLIM_array = opt.autoLIM_array(opt.select_map_numbers);
opt.LIM_array_par_maps = opt.LIM_array_par_maps(opt.select_map_numbers);
opt.LIM_array_sub_maps = opt.LIM_array_sub_maps(opt.select_map_numbers);


for nSubtraction = 1:length(opt.show_subtraction_maps);
    show_subtraction_maps = opt.show_subtraction_maps(nSubtraction);

    show_rois = 0;
    if isfield(opt.fig,'show_rois')
        if opt.fig.show_rois & ~isempty(roi_names)
            show_rois = 1;
        end
    end


    % ----------------- subtract merged names
    if show_subtraction_maps
        merged_names_sub = {};
        for n = 1:numel(merged_names)-1
            merged_names_sub{n} = sprintf('%s-%s', merged_names{n+1}, merged_names{n}) ;
        end
        if isempty(merged_names_sub)
            show_subtraction_maps = 0;
        else
            merged_names = merged_names_sub;
        end
    end

    for n_data = 1:numel(merged_names)
        merged_name = merged_names{n_data};

        % find common zoom factor for all maps in a session/experiment
        Ix = 0; Iy = 0;
        for n_map = 1:n_maps
            map_name = opt.map_names{n_map};

            % map
            name = [merged_name opt.map_name_prepend map_name '.nii.gz'];
            nii_fn = fullfile(root_data_path,name);
            display(sprintf('%s',nii_fn))
            nii = abs(mdm_nii_read(nii_fn));
            I = sum(nii,3);
            Ix = Ix + sum(I,2); Ix = Ix/max(Ix);
            Iy = Iy + sum(I,1); Iy = Iy/max(Iy);
        end
        Ix = Ix/n_maps;
        Iy = Iy/n_maps;

        zoom_x = find(Ix > opt.data_zoom_thresh);
        zoom_y = find(Iy > opt.data_zoom_thresh);


        for n_map = 1:n_maps
            map_name = opt.map_names{n_map};

            % map
            name = [merged_name opt.map_name_prepend map_name '.nii.gz'];
            nii_fn = fullfile(root_data_path,name);
            [nii, h_nii] = mdm_nii_read(nii_fn);

            % post mask
            if isfield(opt,'post_mask')
                post_mask = mdm_nii_read(fullfile(root_data_path,opt.post_mask));
                nii = nii .* double(post_mask);
            end


            data_name = extractBefore(root_data_path,'_out');
            [~, data_name] = fileparts(data_name);

            [~, fig_name] = fileparts(nii_fn);
            fig_name = extractBefore(fig_name, opt.title_strings{n_map});
            fig_name = sprintf('%s_%s%s', data_name, fig_name, opt.title_strings{n_map});

            LIM = [min(nii(:)) max((nii(:)))];
            isGreyScale = LIM(1) >= 0;

            if ~opt.autoLIM_array(n_map)
                if isGreyScale
                    LIM = [0 1];
                else
                    LIM = [-1 1];
                end

                if show_subtraction_maps
                    LIM = opt.LIM_array_sub_maps(n_map) * LIM;
                else
                    LIM = opt.LIM_array_par_maps(n_map) * LIM;
                end
            end

            if isGreyScale
                col_map = gray;
            else
                col_map = turbo; %jet;
                %col_map = magma; % magma; plasma %inferno; cividis viridis
            end


            size_nii = size(nii);
            n_slices = size(nii,3);

            if isempty(opt.fig.LimitSliceRange)
                slice_ind = 1:n_slices;
            else
                slice_ind = opt.fig.LimitSliceRange(1):opt.fig.LimitSliceRange(2);
            end

            [opt.fig, fh] = SPAS_map_fig_opt(n_map, size_nii .* h_nii.pixdim(2:4)', zoom_x, zoom_y, opt.fig);

            for n = 0:length(slice_ind)-1

                n_col = rem(n,opt.fig.n_cols);
                n_row = floor(n/opt.fig.n_cols)+1;

                ah = axes('position',...
                    [n_col * (opt.fig.sub_fig_width + opt.fig.x_space) ...
                    (opt.fig.n_rows - n_row) * (opt.fig.sub_fig_height + opt.fig.y_space) ...
                    opt.fig.sub_fig_width opt.fig.sub_fig_height]);


                I = real(squeeze(nii(zoom_x,zoom_y,slice_ind(n+1))))';
                imagesc(I);

                if ~isempty(LIM)
                    caxis(LIM)
                end

                colormap(ah,col_map);

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
                            plot_roi(ROI, opt.fig.roi_colors(n_roi,:), opt.fig.roi_lw);
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


                axis(ah,'tight','off')
                pbaspect(ah,[1 opt.fig.im_ratio 1] .* h_nii.pixdim(2:4)')


            end


            if opt.fig.showTitle
                ah = axes('position',[0 1-opt.fig.title_size 1 opt.fig.title_size]);
                axis(ah,'tight','off');
                th = text(0.01, .5, fig_name, 'color', 'white', 'fontsize', opt.fig.title_fs);
                th.Interpreter = 'none';
            end

            % add colorbar

            ah = axes('position',...
                [opt.fig.n_cols * (opt.fig.sub_fig_width + opt.fig.x_space) + opt.fig.x_space ...
                .1*opt.fig.n_rows*opt.fig.sub_fig_height ...
                opt.fig.colormap_space .8*opt.fig.n_rows*opt.fig.sub_fig_height]);


            ah.Visible = 0;

            caxis(LIM)
            colormap(ah,col_map);
            cbh = colorbar;
            cbh.Location = 'westoutside'; %'west'; %'westoutside'; %'northoutside';
            cbh.Position(3) = opt.fig.colormap_space/4;

            %cbh.Ticks = [];
            cbh.Box = 'off';
            cbh.FontSize = 12;
            cbh.TickDirection = 'in';
            cbh.TickLength = .02;

            cbh.Color = 'white';

            display(sprintf('dataset %d/%d: %s: %s', n_data, numel(merged_names), merged_name, fig_name))


            if opt.fig.save
                fig_path = fullfile(root_data_path,'fig');
                mkdir(fig_path);
                fig_path = fullfile(fig_path, fig_name);
                print(fh,fig_path,'-dpng',opt.fig.resolution);
                display(sprintf('saved %s', fig_path))
            end

        end
    end
end


end