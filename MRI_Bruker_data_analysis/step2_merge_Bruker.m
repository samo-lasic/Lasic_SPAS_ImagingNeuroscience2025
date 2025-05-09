% Merge waveforms, NIfTI files, and XPS in 'extracted_nii' subfolder.
% Works on data separated from a single multi-FWF sequence or acquired in separate experiments.


clear all
close all

setup_code_path()

warning ('off','all'); % Suppress all warnings during merge

data_opt.append_folder_name = '_out';

% Optional: override root data path by setting data_opt.root_data_path before calling setup_data_path
% data_opt.root_data_path = ... 
data_path_struct = setup_data_path(data_opt);

% ------------------------------------
select_wfm_names = {}; %{'SPAS1','SPAS3'}; % if empty, merge all
select_day = {}; % use select_day = {} to omit this filter

% select_wfm_names = {'STE1','tLTE'};

% ------------------------------------
for c = 1:numel(data_path_struct)

    root_data_path = data_path_struct(c).root_data_path;
    exp_folder_name = data_path_struct(c).exp_folder_name;
    select_subfolders = data_path_struct(c).select_subfolders;

    % if select_subfolders is empty do all
    if isempty(select_subfolders)
        d = dir(fullfile(root_data_path,exp_folder_name));
        d = d(~contains({d.name},'.'));
        for nf = 1:numel(d)
            select_subfolders(end+1) = str2num(d(nf).name);
        end
    end


    for nsub = 1:numel(select_subfolders)
        select_subfolder_name = num2str(select_subfolders(nsub));
        display(sprintf('merging: %s', fullfile(root_data_path,exp_folder_name, select_subfolder_name)))

        % find all nifti files in the folder
        files = dir([fullfile(root_data_path, exp_folder_name, select_subfolder_name) '/**/*.nii.gz']);

        % filter all relevant files and make a list of paths
        nii_fn = {};
        wfm_name = {};
        for n = 1:length(files)
            folder = files(n).folder;
            name = files(n).name;

            if ~contains(name,{'merged','roi_','all_wfm'})

                ind = strfind(name,'_');
                wfm_name{end+1} = extractBefore(name,ind(end-1));

                nii_fn{end+1} = fullfile(folder,name);
            end
        end

        % make sure there are no waveform name duplicates
        [nii_fn,ind] = unique(nii_fn);
        nii_fn = nii_fn(ind);

        % select by day
        if ~isempty(select_day)
            day_sel = [];
            for n = 1:length(nii_fn)
                name = nii_fn{n};
                ind = strfind(name,'_');
                name = extractAfter(name,ind(end-1));
                name = extractBefore(name,'_');
                if ismember(name(5:end),select_day)
                    day_sel(end+1) = n;
                end
            end
            nii_fn = nii_fn(day_sel);
        end


        if ~isempty(select_wfm_names)
            path_select = [];
            for n = 1:numel(select_wfm_names);
                for m = 1:length(nii_fn)
                    name = nii_fn{m};
                    [~,name]=fileparts(name);

                    if strcmp(extractBefore(name,'_'),select_wfm_names{n})
                        path_select(end+1) = m;
                    end
                end
            end
            nii_fn = nii_fn(path_select);
        end



        % -------------------------------  merging ----------------------------------------
        if isempty(select_subfolder_name)
            out_nii_fn = 'merged.nii.gz';
        else
            out_nii_fn = ['merged_' select_subfolder_name '.nii.gz'];
        end
        out_nii_fn = fullfile(root_data_path, exp_folder_name, out_nii_fn);

        if (0) % Debug: plot signal vs. b-value 
            data = mdm_nii_read(nii_fn{1});
            xps = mdm_xps_load(mdm_fn_nii2xps(nii_fn{1}));

            data_sz = size(data);

            cc = round(data_sz/2);
            data1 = squeeze(data(cc(1),cc(2),cc(3),:));

            ind = find(xps.rot_ind == 4);
            b = xps.b(ind);
            sig = data1(ind);

            figure(1),clf
            semilogy(b,sig,'o')
            return
        end



        % check which xps can be merged (they have the same xps.n)
        ind = [];
        for n = 1:numel(nii_fn)
            xps = mdm_xps_load(mdm_fn_nii2xps(nii_fn{n}));
            ind(end+1) = xps.n;
        end

        uind = unique(ind);
        lind = [];
        for n = 1:length(uind)
            lind(n) = length(find(ind == uind(n)));
        end
        [~,tmp] = max(lind);
        select_ind = find(uind(tmp) == ind);

        if length(select_ind) < numel(nii_fn)
            display('excluded datasets due to non-mergable xps:')
            exclude_ind = find(uind(tmp) ~= ind);
            for m = 1:length(exclude_ind)
                display(nii_fn{exclude_ind(m)})
            end
        end

        nii_fn = nii_fn(select_ind);

        out_nii_fn = mdm_nii_merge(nii_fn, out_nii_fn);

        % merge xps including wfm source ----------
        % xps = mdm_xps_merge(mdm_fn_nii2xps(nii_fn)); % Cannot use mdm_xps_merge because XPS contains string fields like 'wfm_src'

        xps_cell = {};
        wfm_src = {};
        for n_fn = 1:numel(nii_fn)
            xps_tmp = mdm_xps_load(mdm_fn_nii2xps(nii_fn{n_fn}));
            if isfield(xps_tmp,'wfm_src')
                wfm_src{n_fn} = xps_tmp.wfm_src;
                xps_cell{n_fn} = rmfield(xps_tmp,'wfm_src');
            else
                wfm_src{n_fn} = '';
                xps_cell{n_fn} = xps_tmp;
            end
        end
        xps = mdm_xps_merge(xps_cell);
        xps.wfm_src = wfm_src;


        % --- add waveform names
        wfm_names = {};
        wfm_names_display = [];
        for n = 1:numel(nii_fn)
            [~, wfm_name_tmp] = fileparts(extractBefore(nii_fn{n},'.nii.gz'));
            wfm_names{n} = wfm_name_tmp;
            wfm_names_display = [wfm_names_display ', ' wfm_name_tmp ];
        end
        display(sprintf('%s', wfm_names_display(3:end)))

        xps.wfm_names = wfm_names; % Store waveform names in XPS for tracking

        mdm_xps_save(xps, mdm_fn_nii2xps(out_nii_fn));
    end

end
