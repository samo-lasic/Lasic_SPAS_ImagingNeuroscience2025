function setup_code_path()
% Initializes MATLAB paths for the MRI Bruker data analysis pipeline. 
% Adds paths to the Multidimensional Diffusion MRI repository and MP-PCA-Denoising repository.
    
path_case = 2; % 1 - server, 2 - local

restoredefaultpath

switch path_case
    case 1
        % Set the path to the Multidimensional Diffusion MRI repository: https://github.com/markus-nilsson/md-dmri
        github_path_mdd = '...';
        
        % Set the path to the MP-PCA-Denoising repository: https://github.com/sunenj/MP-PCA-Denoising
        % Veraart, J., Novikov, D. S., Christiaens, D. & Ades-aron, B. NeuroImage Denoising of diffusion MRI using random matrix theory. Neuroimage 142, 394–406 (2016).
        % Does, M. D. et al. Evaluation of principal component analysis image denoising on multi-exponential MRI relaxometry. Magn. Reson. Med. 81, 3503–3514 (2019).
        github_path_denoise = '...';
    case 2
         % Set the path to the Multidimensional Diffusion MRI repository: https://github.com/markus-nilsson/md-dmri
        github_path_mdd = '...';
        
        % Set the path to the MP-PCA-Denoising repository: https://github.com/sunenj/MP-PCA-Denoising
        % Veraart, J., Novikov, D. S., Christiaens, D. & Ades-aron, B. NeuroImage Denoising of diffusion MRI using random matrix theory. Neuroimage 142, 394–406 (2016).
        % Does, M. D. et al. Evaluation of principal component analysis image denoising on multi-exponential MRI relaxometry. Magn. Reson. Med. 81, 3503–3514 (2019).
        github_path_denoise = '...';
end


run(fullfile(github_path_mdd,'setup_paths.m'));
addpath(genpath(fullfile(pwd,'functions')));
addpath(genpath(github_path_denoise));

end