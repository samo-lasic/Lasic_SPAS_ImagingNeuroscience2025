function out_nii_fn = SPAS_geo_average(nii_fn)
% geometric averaging of SPAS data (after direction averaging)

geo_wfm_names = {'SPAS1','SPAS2','SPAS3'};

xps = mdm_xps_load(mdm_fn_nii2xps(nii_fn));
[sig,h] = mdm_nii_read(nii_fn);

has_good_vol = isfield(xps,'good_vol');

% find indices for geo_wfm_names
geo_ind = find(contains(xps.wfm_names,geo_wfm_names));
if length(geo_ind) < 2
    out_nii_fn = 'less then 2 SPAS waveforms in this dataset!';
    return 
end

% ---   geometric average ----
uwfm = unique(xps.s_ind);

udir = unique(xps.rot_ind);
ub = unique(xps.a_ind);

cnt = 0;
for n_dir = 1:length(udir)
    for n_b = 1:length(ub)
        ind = find(ismember(xps.s_ind,geo_ind) & udir(n_dir) == xps.rot_ind & ub(n_b) == xps.a_ind);
        if length(ind) > 1
            cnt = cnt+1;
            sig_g(:,:,:,cnt) = geomean(sig(:,:,:,ind),4);
            ind_g(cnt) = ind(1);
            %[ind xps.b(ind)*1e-9 xps.rot_ind(ind)]

            if has_good_vol
                good_vol(cnt,1) = prod(xps.good_vol(ind)); % if any is 0 then set good_vol = 0;
            end
        end
    end
end

% append signal
sig(:,:,:,end+(1:length(ind_g))) = sig_g;

% build geo xps
f = fieldnames(xps);
remove_fields = {'n','wfm_names','wfm_src'};
f = f(~ismember(f,remove_fields));
xps_append.n = length(ind_g);

for i = 1:numel(f)
    xps_append.(f{i}) = xps.(f{i})(ind_g, :,:);
end

if has_good_vol
    xps_append.good_vol = good_vol;
end

% give it a new s_ind
xps_append.s_ind = ones(xps_append.n,1)*(max(xps.s_ind)+1);

% append geo xps
xps_g = rmfield(xps,remove_fields);
for i = 1:numel(f)
    for j = 1:xps_append.n
        xps_g.(f{i})(xps.n+j,:,:) = xps_append.(f{i})(j, :,:);
    end
end

f = f([1 1:end]);
f{1} = 'n';
xps_g.n = xps.n + xps_append.n;
xps_g = orderfields(xps_g,f);

% return back removed fields containing strings
xps_g.wfm_names = xps.wfm_names;
xps_g.wfm_src = xps.wfm_src;

xps_g.wfm_names{end+1} = ['geoSPAS_' char(extractAfter(xps.wfm_names(geo_ind(1)),'_'))];
xps_g.wfm_src{end+1} = 'geoSPAS'; 

out_nii_fn = append_nii_fn(nii_fn, 'g');
mdm_nii_write(single(sig),out_nii_fn,h);
mdm_xps_save(xps_g, mdm_fn_nii2xps(out_nii_fn));

