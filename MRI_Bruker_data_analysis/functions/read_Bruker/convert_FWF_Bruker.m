function convert_FWF_Bruker(root_data_path, exp_folder_name, select_subfolders, opt)

if ~isfield(opt,'root_out_path')
    root_out_path = root_data_path; %fullfile(root_data_path,[exp_folder_name '_out']);
else
    root_out_path = opt.root_out_path;
end

if ~isfield(opt,'pdata') % set default pdata used for analysis
    opt.pdata = 1; 
end
display(sprintf('using pdata/%d', opt.pdata))


gmr = 26.75e7;

files = dir(fullfile(root_data_path, exp_folder_name,'**','2dseq'));



% ---------------- collect basic info of multiple datasets / select / sort ----------------------------

cnt = 0;
for n = 1:numel(files)
    path2dseq = fullfile(files(n).folder,files(n).name);
    [data_root,~]=fileparts(fileparts(files(n).folder));

    [~,data_folder_name] = fileparts(data_root);
    data_folder_num = str2num(data_folder_name);
    if isempty(data_folder_num)
        data_folder_num = -1;
    end

    if ismember(data_folder_num, select_subfolders) || isempty(select_subfolders)

        cnt = cnt+1;
        datasets.to_analyze(cnt) = 1; % make it the first field (for later filtering)

        try
            scan = read_paravision_data(path2dseq);
            %scan.method.Method

            datasets.data_folder_name{cnt} = data_folder_name;
            datasets.data_folder_num(cnt) = data_folder_num;


            name = scan.method.DwGradShapeStrArr1;
            
            % patch additional waveform names missing the pre_ substring
            name_array = split(name,' ');
            not_ind = ~contains(name_array,'pre_');
            name_array(not_ind) = strcat('pre_', name_array(not_ind));
            name = char(join(name_array,' '));

            ind = findstr(name,'pre_');
            Nwfms = length(ind);

            if Nwfms > 1

                for nwfms = 1:Nwfms
                    if nwfms < Nwfms
                        method_name = name(ind(nwfms)+4:ind(nwfms+1)-2);
                        ind1 = findstr(method_name,'_');
                    else
                        method_name = name(ind(nwfms)+4:end);
                        ind1 = findstr(method_name,'_');
                    end

                    if nwfms == 1
                        method_names = method_name;
                        if isempty(ind1)
                            wfm_names = method_name;
                        else
                            wfm_names = [method_name(1:ind1(1)-1)];
                        end
                    else
                        method_names = [method_names method_name];
                         if isempty(ind1)
                            wfm_names = [wfm_names method_name];
                        else
                            wfm_names = [wfm_names method_name(1:ind1(1)-1)];
                        end
                        
                    end

                    if nwfms < Nwfms
                        method_names = [method_names '&'];
                        wfm_names = [wfm_names '&'];
                    end

                end

                datasets.method_name{cnt} = method_names;
                datasets.wfm_name{cnt} = wfm_names;


            else
                ind = findstr(name,'_');
                datasets.method_name{cnt} = name(ind(1)+1:ind(2)-1);
                datasets.wfm_name{cnt} = name(ind(1)+1:end);
            end
            datasets.Nwfms{cnt} = Nwfms;

            %         datasets.method_name{cnt} = char(extractBetween(scan.method.DwGradShapeStrArr1,'_','_'));
            %         datasets.wfm_name{cnt} = char(extractAfter(scan.method.DwGradShapeStrArr1,'_'));

            acq_time = strrep(extractBefore(scan.acqp.ACQ_time,','),'-','');
            acq_time = strrep(acq_time,'T','_');
            acq_time = strrep(acq_time,':','');
            datasets.acq_time{cnt} = acq_time;
            datasets.acq_time_hms{cnt} = extractAfter(acq_time,'_');

            datasets.path2dseq{cnt} = path2dseq;
            datasets.data_root{cnt} = data_root;

            datasets.DwDirsFilename{cnt} = scan.method.DwDirsFilename;
        

        catch
            datasets.to_analyze(cnt) = 0;
            display(['incomplete dataset: ' data_folder_name])
            cnt = cnt-1;
        end

    end
end

if numel(files) == 0 
    display(sprintf('%s', 'no data found'));
    return
end

display(sprintf('good datasets: %s', num2str(datasets.data_folder_num)));

% find dummy scan based on waveform name (rather than method_name)
% isDummy = contains(datasets.method_name,'dummy');
isDummy = contains(datasets.wfm_name,'dummy');


if sum(isDummy) > 0
    display(sprintf('dummy datasets: %s', num2str(datasets.data_folder_num(isDummy))));
end

% filter datasets to analyze
fs = fields(datasets);
for n = 2:length(fs)
    datasets1.(fs{n}) = datasets.(fs{n})(find(datasets.to_analyze));
end
datasets = datasets1;
clear datasets1

% select pdata (remove alternative pdata)
pdata_ind = find(contains(datasets.path2dseq,fullfile('pdata',num2str(opt.pdata),'2dseq')));
fs = fields(datasets);
for n = 1:length(fs)
    datasets1.(fs{n}) = datasets.(fs{n})(pdata_ind);
end
datasets = datasets1;


% remove dummy
if ~opt.include_dummy_wfm
    isDummy = contains(datasets.wfm_name,'dummy');
    fs = fields(datasets);
    for n = 1:length(fs)
        datasets1.(fs{n}) = datasets.(fs{n})(~isDummy);
    end
    datasets = datasets1;
end

% save info
if opt.save_info
    info_path = fullfile(root_data_path,exp_folder_name,'info.txt');
    fileID = fopen(info_path,'w');

    % sort by time
    %[~,ind] = sort(datasets.acq_time);
    %[~,ind] = sort(datasets.acq_time_hms);

    % sort by data_folder_num
    [~,ind] = sort(datasets.data_folder_num);

    for n = 1:length(ind)
        index = ind(n);

        fprintf(fileID,'%s:\r\n',datasets.data_folder_name{index});
        fprintf(fileID,'%s\r\n',datasets.acq_time{index});
        fprintf(fileID,'%s\r\n',datasets.wfm_name{index});
        fprintf(fileID,'%s\r\n',datasets.DwDirsFilename{index});
        fprintf(fileID,'%s\r\n',datasets.method_name{index});

        fprintf(fileID,'\r\n');

    end

    fclose(fileID);
end



% sort by time
[~,ord] = sort(datasets.acq_time);
fs = fields(datasets);
for n = 1:length(fs)
    datasets.(fs{n}) = datasets.(fs{n})(ord);
end

day = extractBefore(datasets.acq_time,'_');
uday = unique(day);

datasets.day_ind = 1:length(day);
for n = 1:length(uday)
    datasets.day_ind(find(strcmp(day,uday(n)))) = n;
end
%     display(datasets)

if numel(datasets.data_root) > 0
    out_path = fullfile(root_out_path,[exp_folder_name '_out']);
    mkdir(out_path)
end



% --------------------------- convert selected datasets ----------------------------

for d = 1:length(datasets.method_name)
    display(sprintf('data: %s/%s', exp_folder_name, datasets.data_folder_name{d}))


    % read data
    scan = read_paravision_data(datasets.path2dseq{d});
    data = scan.dataset;


    if opt.outputSubFolders
        %name of folder containing pdata
        [~,out_sub_path] = fileparts(extractBefore(datasets.path2dseq{d},[filesep 'pdata']));
        out_root_path_loop = fullfile(out_path, out_sub_path);
    else
        out_root_path_loop = out_path;
    end

    out_nii_folder = fullfile(out_root_path_loop, 'extracted_nii');
    mkdir(out_nii_folder)

    % test
    if (0)
        scan.method.DwGradAmpG
        scan.method.PVM_NRepetitions
        % slices and volumes
        Nsl = scan.method.PVM_SlicesPerGate;
        N = size(data,1)/Nsl; % repetitions x waveforms

        Nvol = N/Nwfms;% per waveform

        % reshape data
        data = reshape(data,Nsl,N,size(data,2),size(data,3));
        data = permute(data,[3,4,1,2]);

        nii_fn = fullfile(out_nii_folder,['all_wfm' '.nii.gz']);
        mdm_nii_write(data, nii_fn);

        % dirty fix read/phase scale by modifying the nifti header 
        if size(data,1) ~= size(data,2)
            h = mdm_nii_read_header(nii_fn);
            h.pixdim(2) = h.dim(3)/h.dim(2);
            mdm_nii_write(data, nii_fn, h);            
        end


    else

        % repetitions
        Nrep = scan.method.PVM_NRepetitions;

        % number of waveforms
        Nwfms = datasets.Nwfms{d};

        Nshapes = scan.method.DwNShapes;

        % slices and volumes
        Nsl = scan.method.PVM_SlicesPerGate;
        Nvol = size(data,1)/Nsl; % all volumes with repetitions
        Nvol_per_rep = Nvol/Nrep; % volumes per reptition
        Nvol_per_waveform = Nvol_per_rep/Nshapes;% per shape

        % reshape data
        data = reshape(data,Nsl,Nvol,size(data,2),size(data,3));

        if opt.rotate_image
            data = permute(data,[4,3,1,2]);
            data = flip(data,2);
        else
            data = permute(data,[3,4,1,2]);
        end


        % gradient amplitudes - common
        G0 = scan.method.DwGradAmpG;
        Gscale = scan.method.DwGradAmpScale*1e-2; % only unique G-values

        Gmax_T_per_m = 1e-3 * max(G0)/max(Gscale); % this is used across all volumes (gradient scales from scan.method.DwGAmpRot..)

        NG0 = scan.method.DwNAmplitudes;

        % rotations - common
        vdir = scan.method.DwDir;
        Ndir = scan.method.DwNDirs; % length(v)/3;
        vdir = reshape(vdir,Ndir,3)';

        % get the input waveform direction (if pre/post directions are equal)
        if prod(scan.method.DwShapeDirVec1 == scan.method.DwShapeDirVec2)
            DwShapeDir = scan.method.DwShapeDirVec1;
        else
            disp('DwShapeDirVec1 and DwShapeDirVec2 are different')
            return
        end

        % waveform timing (per waveform ?)
        t1 = scan.method.DwGradDur1*1e-3;
        t2 = scan.method.DwGradDur2*1e-3;

        t180 = scan.method.DwGradTsep*1e-3;

        % spoiler - this is not used and verified !
        tSpoler = scan.method.SliceSpoilerDur; %ms
        tSpoilerRamp = 0.15;
        tPlateau = tSpoler - 2*tSpoilerRamp;

        % plateau + ramps =1.63; Plateau only: 1.33

        % check BaseLevelRotations.c
        %scan.acqp.ACQ_gradient_amplitude
        %ACQ_gradient_amplitude[0]  =  SliceGrad;
        %ACQ_gradient_amplitude[12] =  SliceSpoiler.ampl;
        %ACQ_gradient_amplitude[4]  = -Phase3DGrad;

        % waveforms (single array for all waveforms)
        g1 = scan.method.DwGradShapeArray1; % fraction of 1 (not in T/m!)
        g2 = scan.method.DwGradShapeArray2;
        g1 = reshape(g1,3,length(g1)/3)';
        g2 = reshape(g2,3,length(g2)/3)';


        NdummyShapes = Nshapes - Nwfms; % corresponds to numel(find(scan.method.DwShapePoints1 == 1)) ?
        non_dummy_ind = find(scan.method.DwShapePoints1 ~= 1); % remove dummy data (this should never be needed except in case of empty waveforms)
       
        data_ind = [];
        for m = 1:Nrep
            for n = 1:numel(non_dummy_ind)
                data_ind = [data_ind  (m-1) * Nvol_per_rep + (non_dummy_ind(n)-1)*Nvol_per_waveform + [1:Nvol_per_waveform] ];
            end
        end


        data = data(:,:,:,data_ind);
        Nvol_per_rep =  Nvol_per_waveform*Nwfms; % corrected for empty waveforms


        Ng1 = length(g1)/Nshapes; % corresponds to scan.method.DwShapePoints1
        Ng2 = length(g2)/Nshapes;
        %scan.method.DwShapePoints1 - number of points in waveform

        % temporary correction for empty waveforms?
        % seems to work generally, not clear why!
        t1 = t1 * (NdummyShapes+Nwfms/Nshapes);
        t2 = t2 * (NdummyShapes+Nwfms/Nshapes);

        dt = t1/(Ng1-1);

        % get rotation matrices
        [DwR, DwGAmpRot, DwGAmp] = getRotationMatrices(scan);

       
        % find repetitions

        % make sure there are no NaN values (this may happen if there are empty waveforms)
        prod_DwGAmpRot = prod(prod(DwGAmpRot,2),3);
        L_nan = numel(find(isnan(prod_DwGAmpRot)));
        L = numel(prod_DwGAmpRot);

        if L_nan > 0
            display(sprintf('%d out of %d volumes in DwGAmpRot are not NaN', L - L_nan, L))
            display(sprintf('%s', 'skipping ...'))
            return
        end

        % find number of b0 - just checking
        isZero = DwGAmpRot == 0;
        sz = size(isZero);
        isZero = reshape(isZero, sz(1),sz(2)*sz(3));
        Nb0 = numel(find(sum(isZero == 1,2) == 9));


        wfm_names = strsplit(datasets.wfm_name{d},'&');
        wfm_srcs = strsplit(datasets.method_name{d},'&');


        wfm_ind1 = [1:Nvol_per_waveform];

        for nwfm = 1:Nwfms
            wfm_name = [wfm_names{nwfm} '_' datasets.acq_time{d}];
            wfm_src = wfm_srcs{nwfm};

            display(sprintf('wfm name = %s', wfm_name))


            wfm_ind = Nvol_per_waveform*(nwfm-1) + wfm_ind1;

            clear data1
            for nrep = 1:Nrep
                ind_in = Nvol_per_rep*(nrep-1) + wfm_ind;
                ind_out = Nvol_per_waveform*(nrep-1) + wfm_ind1;
                data1(:,:,:,ind_out) = data(:,:,:,ind_in);
            end

            nii_fn = fullfile(out_nii_folder,[wfm_name '.nii.gz']);
            mdm_nii_write(data1, nii_fn);

            % dirty fix read/phase scale by modifying the nifti header
            if size(data1,1) ~= size(data1,2)
                h = mdm_nii_read_header(nii_fn);
                h.pixdim(2) = h.dim(3)/h.dim(2);
                mdm_nii_write(data1, nii_fn, h);
            end


            DwR1 = DwR(wfm_ind,:,:);
            DwGAmp1 = DwGAmp(wfm_ind);
            DwGAmpRot1 = DwGAmpRot(wfm_ind,:,:);
            DwShapeDir1 = DwShapeDir(3*(nwfm-1) + [1:3]);


            % start building xps

            gPre    = g1((NdummyShapes + nwfm - 1) * Ng1 + [1:Ng1],:);
            gPost   = g2((NdummyShapes + nwfm - 1) * Ng2 + [1:Ng2],:);

            % interpolate gradient waveforms
            tin = [linspace(0,t1,Ng1) t1 + t180 + linspace(0,t2,Ng2)];
            gin = [gPre; -gPost];
            te = t180 + t1 + t2;
            t = 0:dt:te;
            g = interp1(tin,gin,t);
            g = g * Gmax_T_per_m; % scale to T/m units

            clf, plot(t*1000,g,'.')

            % dephasing and b-tensor

            q = gmr*cumsum(g)*dt;
            bt = q'*q*dt;
            b0 = trace(bt);
            bt = bt/b0;

            clear xps
            % per repetition first (then replicate for each repetition)
            xps.b = (DwGAmp1'*1e-2).^2*b0;
            xps.gscale = DwGAmp1';

            for n = 1:Nvol_per_waveform
                DwR3x3 = squeeze(DwR1(n,:,:));
                bt_rot = DwR3x3*bt*DwR3x3'; % rotate b-tensor
                xps1 = mdm_xps_from_bt(xps.b(n) * bt_rot);

                if sum(isnan(xps1.u)) % u vector not defined for asymmetric tensors
                    [eigu, eigv] = eig(bt_rot);
                    [~, ind] = max([eigv(1,1) eigv(2,2) eigv(3,3)]);
                    xps1.u = eigu(:,ind)';
                end

                xps.bt(n,:) = xps1.bt;
                xps.b_delta(n,:) = xps1.b_delta;
                xps.b_eta(n,:) = xps1.b_eta;
                xps.v(n,:) = xps1.u; % "main" eigen vector
                xps.DwShapeDir(n,:) = DwShapeDir1;

                xps.DwR(n,:) = rot_3x3_to_1x6(DwR3x3);
                xps.DwR3x3(n,:,:) = DwR3x3; % temporary

                % intended rotations
                xps.u(n,:) = DwR3x3*DwShapeDir1';


            end

            % find indices of unique b-tensor shapes and b-values
            %[~, ~, aind] = uniquetol([xps.b_delta xps.b_eta xps.b],1e-6,'ByRows',1);
            [~, ~, aind] = uniquetol(xps.b,1e-6);
            xps.a_ind = aind;


            % find indices of unique rotation matrices
            rot_ind = zeros(Nvol_per_waveform,1);
            for n = 1:Nvol_per_waveform
                if ~rot_ind(n)
                    rot_ind(n) = n;
                end

                for m = n+1:Nvol_per_waveform
                    if isequal(xps.DwR(m,:),xps.DwR(n,:)) && ~rot_ind(m)
                        rot_ind(m) = n;
                    end
                end
            end
            % sequence of unique rot_ind
            urot_ind = unique(rot_ind);
            for n = 1:numel(urot_ind)
                xps.rot_ind(rot_ind == urot_ind(n),1) = n;
            end

            % find anti-parallel pairs per b-value (each b-value gets a different set of antipodal pairs, these volumes could be averaged!)
            ind_antip = zeros(Nvol_per_waveform,1);
            cnt = 0; % count used indices
            for n_a_ind = 1:max(xps.a_ind)
                ind = find(xps.a_ind == n_a_ind);
                u_tmp = xps.u(ind,:);
                A = u_tmp*u_tmp';
                [i,j] = find(abs(triu(A) + 1)<1e-3); % look into the upper part of the triangular matrix to find antipodal pairs

                for m = 1:length(i)
                    ind_antip(ind(i(m))) = cnt + m;
                    ind_antip(ind(j(m))) = cnt + m;
                end
                cnt = cnt + length(i); % count used indices
            end
            xps.ap_ind = ind_antip;


            % add repetitions to xps
            if Nrep > 1
                f = fields(xps);
                for n = 1:numel(f)
                    xps.(f{n}) = repmat(xps.(f{n})(:,:,:), Nrep, 1, 1);
                end
            end

            % index of repetion dimension
            rep_dim_ind = repmat([1:Nrep]',1,Nvol_per_waveform)';
            xps.rep_dim_ind = rep_dim_ind(:);

            % index for averaging repetitions
            xps.rep_a_ind = repmat(1:Nvol_per_waveform,1,Nrep)';

            % waveform source name (as it was created)
            xps.wfm_src = wfm_src;
            xps.n = Nvol_per_waveform * Nrep;

            % reorder fields (cosmetics)
            f = fields(xps);
            f = f(~ismember(f,{'n'}));

            xps1.n = xps.n;
            for n = 1:numel(f)
                xps1.(f{n}) = xps.(f{n});
            end
            xps = xps1;


            % save xps
            mdm_xps_save(xps, mdm_fn_nii2xps(nii_fn));

            % save waveform
            if opt.save_waveforms

                out_wfm_folder = fullfile(out_root_path_loop, 'extracted_wfm');
                mkdir(out_wfm_folder)

                wfm_fn = fullfile(out_wfm_folder,['wfm_' wfm_name '.mat']);
                fig_fn = fullfile(out_wfm_folder,['wfm_' wfm_name]);
                wfm.t = t;
                wfm.g = g;

                save(wfm_fn,'wfm')
            end


            if (opt.show_waveform)
                fh = figure(1);
                clf
                if opt.hide_fig
                    fh.Visible = 'off';
                end
                subplot(2,1,1)
                hold on
                plot(t*1e3, g(:,1),'.r')
                plot(t*1e3, g(:,2),'.g')
                plot(t*1e3, g(:,3),'.b')
                ylabel('gradient')

                title(wfm_name,'interpreter','none')

                subplot(2,1,2)
                hold on
                plot(t*1e3, q(:,1),'.r')
                plot(t*1e3, q(:,2),'.g')
                plot(t*1e3, q(:,3),'.b')
                ylabel('dephasing')
                xlabel('time [ms]')

                if opt.save_waveforms
                    print(fig_fn,'-dpng','-r150')
                end

                pause(.1)

            end
        end
    end
end





