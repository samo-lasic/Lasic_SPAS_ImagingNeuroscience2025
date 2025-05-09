function SPAS_make_sig_dif_maps(root_data_path, merged_names, roi_names, opt)
% used to show / save parameter maps and subtraction maps


% refine select_map_numbers based on existing files
files = dir(root_data_path);
ind = contains({files.name},merged_names);
files = files(ind);

select_map_numbers = [];
for n = 1:numel(opt.map_names)
    if find(contains({files.name},opt.map_names{n}))
        select_map_numbers = [select_map_numbers n];
    end
end
ind = ismember(opt.select_map_numbers,select_map_numbers);
opt.select_map_numbers = opt.select_map_numbers(ind);



% take a subset of selected maps
n_maps = numel(opt.select_map_numbers);
opt.map_names = opt.map_names(opt.select_map_numbers);
opt.title_strings = opt.title_strings(opt.select_map_numbers);

opt.autoMIN_array = opt.autoMIN_array(opt.select_map_numbers);
opt.autoMAX_array = opt.autoMAX_array(opt.select_map_numbers);
opt.MIN_array = opt.MIN_array(opt.select_map_numbers);
opt.MAX_array = opt.MAX_array(opt.select_map_numbers);


show_rois = 0;
if isfield(opt.fig,'show_rois')
    if opt.fig.show_rois & ~isempty(roi_names)
        show_rois = 1;
    end
end


for n_data = 1:numel(merged_names)
    merged_name = merged_names{n_data};

    % find common zoom factor (pixel range not intensity!) for all maps in a session/experiment
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

    if (opt.do_zoom_range)
        min_x = min(zoom_x);
        max_x = max(zoom_x);
        length_x = max_x - min_x;
        min_y = min(zoom_y);
        max_y = max(zoom_y);
        length_y = max_y - min_y;

        zoom_x = zoom_x(find(zoom_x > min_x + opt.zoom_range_x(1)*length_x & zoom_x < min_x + opt.zoom_range_x(2)*length_x));
        zoom_y = zoom_y(find(zoom_y > min_y + opt.zoom_range_y(1)*length_y & zoom_y < min_y + opt.zoom_range_y(2)*length_y));
    end


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

        % alpha mask to modulate intensity
        if opt.modulate_intensity.type > 0
            if opt.modulate_intensity.type == 1
                str_append = '_smin';
            elseif opt.modulate_intensity.type == 2
                str_append = '_smax';
            end
            str_append = strcat(extractBefore(map_name,'_'),str_append);
            mask_name = strrep(name,map_name,str_append);
            alpha_mask_nii_fn = fullfile(root_data_path, mask_name);
            alpha_mask = mdm_nii_read(alpha_mask_nii_fn);
            alpha_mask = alpha_mask / max(alpha_mask(:));
            if opt.modulate_intensity.thresh > 0
                alpha_mask(alpha_mask < opt.modulate_intensity.thresh) = 0;
            end

            if opt.modulate_intensity.gamma_cor.p ~= 1
                alpha_mask = img_gamma_correction(alpha_mask, opt.modulate_intensity.gamma_cor);
            end
        end


        data_name = extractBefore(root_data_path,'_out');
        [~, data_name] = fileparts(data_name);

        [~, fig_name] = fileparts(nii_fn);
        fig_name = extractBefore(fig_name, opt.title_strings{n_map});

        mod_str = [];
        if opt.modulate_intensity.type == 1
            mod_str = '_mod_s_min';
        elseif opt.modulate_intensity.type == 2
            mod_str = '_mod_s_max';
        end

        fig_name = sprintf('%s_%s%s%s', data_name, fig_name, opt.title_strings{n_map}, mod_str);


        if opt.autoMIN_array(n_map)
            LIM(1) = min(nii(:));
        else
            LIM(1) = opt.MIN_array(n_map);
        end

        if opt.autoMAX_array(n_map)
            LIM(2) = max(nii(:));
        else
            LIM(2) = opt.MAX_array(n_map);
        end

        %col_map = [ flipud(cividis); inferno]; % viridis
        % col_map = [flipud(gray); copper];
        %col_map = inferno; %magma; turbo; %jet;

        eval(sprintf('col_map = %s;',opt.colormap))

        size_nii = size(nii);
        n_slices = size_nii(3);


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

            % nonlinear contrast adjustment
            if isfield(opt,'gamma_cor')
                I = img_gamma_correction_lim(I, LIM, opt.gamma_cor);
            end


            hI = imagesc(I);

            if opt.modulate_intensity.type > 0
                % apply mask
                A = real(squeeze(alpha_mask(zoom_x,zoom_y,slice_ind(n+1))))';
                set(hI, 'AlphaData', A);
            end

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
        cbh.Position(3) = .25 * opt.fig.colormap_space; %!!!
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
