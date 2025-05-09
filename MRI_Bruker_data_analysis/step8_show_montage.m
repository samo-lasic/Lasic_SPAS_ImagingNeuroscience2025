% Create a montage from saved figures in the `fig` subfolders

clear all
close all

%restoredefaultpath
setup_code_path()

warning('off','all')

% ---------------------
data_opt.append_folder_name = '_out';
data_path_struct = setup_data_path(data_opt);
% ---------------------

input.data_path_struct = data_path_struct;

output.dir = fullfile(data_opt.root_data_path,'montage');

input.file_extension = 'png';

input.fig_type = [1 2 3]; %[1 2 3]; %1-normal map, 2-subtraction map, 3-histogram

% -------------------- maps
input.map_names = {'dti_s0','dti_md','dti_fa','norm_pa_raw_TDD_dlogs'};
input.map_names = input.map_names(2);
%input.map_names = input.map_names(2:4);

% Choose prefix for input map names (uncomment desired option)

% input.map_name_prepend = '_dn_g_'; % '_s_a_g_'; % _s_g_  _s_a_g_
input.map_name_prepend = '_g_';
%input.map_name_prepend = '_dn_g_';
%input.map_name_prepend = '_dn_mc_g_';


output.saveFig = 1;
output.file_extension = 'png'; %'jpg'; 'png';
output.closeFig = 1;
output.arrange = 2;

make_montage(1, input, output);


% -------------------- functions -----------
function [] = make_montage(figNum, input, output)

root_data_paths = {};
for n = 1:numel(input.data_path_struct)
    root_data_paths{end+1} = fullfile(input.data_path_struct(n).root_data_path, input.data_path_struct(n).exp_folder_name);
end

N_fig_type = numel(input.fig_type);
Nsubjects = numel(root_data_paths);
Nmaps = numel(input.map_names);


for n_fig_type = 1:N_fig_type

    fig_type = input.fig_type(n_fig_type);

    cnt = 0;
    image_paths = {};
    clear montage_array size_array montage_images

    for nSubject = 1:Nsubjects

        root_data_path = root_data_paths{nSubject};
        files = dir([root_data_path '/fig/*.' input.file_extension]);
        % filter files by map_name_prepend
        files = files(contains({files.name}, input.map_name_prepend));

        switch fig_type
            case 1 % normal map
                files1 = files(~contains({files.name},'-') & ~contains({files.name},'_hist'));
                title_str = 'map montage';
                fig_type_str = '';
            case 2 % subtraction map
                files1 = files(contains({files.name},'-') & ~contains({files.name},'_hist'));
                title_str = 'subtraction map montage';
                fig_type_str = 'subtraction';
            case 3 % histogram
                files1 = files(contains({files.name},'_hist'));
                title_str = 'histogram montage';
                fig_type_str = 'hist';
        end


        for nMap = 1:Nmaps
            exp_inds = find(contains({files1.name}, input.map_names(nMap)));

            Nexp = numel(exp_inds);


            for nExp = 1:Nexp
                exp_ind = exp_inds(nExp);
                cnt = cnt+1;

                file = files1(exp_ind);

                image_paths{cnt} = fullfile(file.folder,file.name);

                montage_array.nSubject(cnt) = nSubject;
                montage_array.nMap(cnt) = nMap;
                montage_array.nExp(cnt) = nExp;

                display(sprintf('%2d | %s (%d/%d) : %s map %d/%d, exp %d/%d', ...
                    cnt, extractBefore(file.name,'_'), nSubject, Nsubjects, fig_type_str, nMap, Nmaps, nExp, Nexp ))
            end
        end
    end

    if cnt > 0
        montage_array.n = cnt;

        for c = 1:montage_array.n
            size_array(c,:) = size(imread(image_paths{c}));
        end
        im_zero = uint8(zeros(max(size_array)));



        % ---------------- arrange montage
        %unique maps
        u_map = unique(montage_array.nMap);

        %unique subjects
        u_subj = unique(montage_array.nSubject);

        %unique combinations of subject and experiment
        u_subj_exp = unique([montage_array.nSubject; montage_array.nExp]','rows')';

        %unique combinations of experiment and map
        u_exp_map = unique([montage_array.nExp; montage_array.nMap]','rows')';


        switch output.arrange
            case 1


                for x = 1:size(u_map,2)
                    map = u_map(x);
                    for y = 1:size(u_subj_exp,2)
                        subj_exp = u_subj_exp(:,y);
                        c = find([montage_array.nSubject == subj_exp(1) & montage_array.nExp == subj_exp(2) & montage_array.nMap == map]);

                        try
                            im_read = imread(image_paths{c});
                            sz = size(im_read);
                        catch
                            im_read = 0 * im_read;
                        end

                        if numel(unique(size_array(:,1))) > 1
                            im = im_zero;
                            im(1:sz(1),1:sz(2),1:sz(3)) = im_read;
                        else
                            im = im_read;
                        end

                        montage_images(x,y).im = im;
                    end
                end



            case 2
                %unique combinations of experiment and map
                u_exp_map = unique([montage_array.nExp; montage_array.nMap]','rows')';

                %     unique subjects
                u_subj = unique(montage_array.nSubject);

                for x = 1:size(u_exp_map,2)
                    exp_map = u_exp_map(:,x);
                    for y = 1:size(u_subj,2)
                        subj = u_subj(y);
                        c = find([montage_array.nSubject == subj & montage_array.nExp == exp_map(1) & montage_array.nMap == exp_map(2)]);

                        try
                            im_read = imread(image_paths{c});
                            sz = size(im_read);
                        catch
                            im_read = 0 * im_read;
                        end

                        if numel(unique(size_array(:,1))) > 1
                            im = im_zero;
                            im(1:sz(1),1:sz(2),1:sz(3)) = im_read;
                        else
                            im = im_read;
                        end

                        montage_images(x,y).im = im;
                    end
                end

        end


        montage_str = 'montage_image = [';
        for y = 1:size(montage_images,2)
            for x = 1:size(montage_images,1)
                montage_str = [montage_str sprintf('montage_images(%d,%d).im ',x,y)];
            end
            montage_str = [montage_str '; '];
        end
        montage_str = [montage_str '];'];
        eval(montage_str)


        % ------------------


        fh = figure;

        fig_name = [];
        if Nsubjects == 1
            [~,subdir] = fileparts(root_data_paths);
            fig_name = [fig_name sprintf('%s', subdir)];%fileparts(input.subdir{1}))];
        else
            fig_name = [fig_name 'MXXXX']; % replace with actual subject identifier 
        end

        if Nmaps == 1
            fig_name = [fig_name sprintf('_map_%s', input.map_names{1})];
        else
            fig_name = [fig_name '_mapXXX']; % replace with actual subject identifier 
        end

        if length(fig_type_str) > 0
            fig_name = [fig_name '_' fig_type_str];
        end

        fig_name = ['montage_' fig_name '.' output.file_extension];

        imshow(montage_image);

        if output.saveFig
            mkdir(output.dir);

            fn = fullfile(output.dir, fig_name);

            imwrite(montage_image, fn);

            display(sprintf('saved %s', fn))
            if output.closeFig
                close(fh)
            end
        end
    end
end
end
