function out_nii_fn = SPAS_make_OP_maps(root_data_path, merged_names, opt)
% makes order parameter (OP) maps from FA and muFA maps

for n_data = 1:numel(merged_names)
    merged_name = merged_names{n_data};

    % check that both maps exist
    names = strcat(merged_name, opt.map_name_prepend, opt.map_names, '.nii.gz');
    files = dir(root_data_path);
    ind = contains({files.name},merged_names);
    ind = ind & contains({files.name},names);
    ind = ind & contains({files.name},{'dti_fa','muFA'});
    files = files(ind);

    if numel(files) == 2


        file = files(contains({files.name},{'dti_fa'}));
        nii_fn = fullfile(root_data_path,file.name);
        [FA, h_FA] = mdm_nii_read(nii_fn);
        FA = double(FA);

        file = files(contains({files.name},{'muFA'}));
        nii_fn = fullfile(root_data_path,file.name);
        [muFA, h_muFA] = mdm_nii_read(nii_fn);
        muFA = double(muFA);

        % Lasic et al., Front Phys 2014;2(11):1–14. http://www.frontiersin.org/Biophysics/10.3389/fphy.2014.00011/abstract
        OP = (FA ./ muFA).^2 .* (3 - 2 * muFA.^2) ./ (3 - 2 * FA.^2);
        OP = sqrt(OP);
        OP(isnan(OP)) = 0;

        %         figure(1),clf, hold on, histogram(FA(:)), histogram(muFA(:))
        %         figure(1),clf, hold on, histogram(OP(:))

        if opt.clamp_OP
            OP(OP < 0) = 0;
            OP(OP > 1) = 1;
        end

        name1 = names{1};
        name2 = names{2};
        for c = 1:numel(name1)
            if ~strcmp(name1(1:c),name2(1:c))
                break
            end
        end
        name = name1(1:c-1);
        if ~strcmp(name(end),'_')
            name = [name '_'];
        end
        name = strcat(name,'OP.nii.gz');
        out_nii_fn = fullfile(root_data_path,name);

        mdm_nii_write(single(OP), out_nii_fn, h_FA);

        display(sprintf('%s',out_nii_fn))


    else
        display('the required files are missing (..dti_fa.. and ..muFA..)')
    end


end
end
