
function [] = mk_RSA_STIMS_ALGEBRA_ols_vol(idx_subj,type)

% type = 'nonH' or type = 'H'
% 'off_diag':
% RDM = [NaN 0   1;
%        0   NaN NaN;
%        1   NaN NaN];

%% prelim 
mk_RDM   = true; % this defines empirical RDMs
mk_stats = true; % this runs GLM for different theoretical RDMs
mk_av    = true; % this projects back into surface space (i.e. do surface RSA)

if nargin==0
    idx_subj  = 1;
    type      = 'nonH';
    RSAmask   = 'Relevant_WholeBrain_RSA_all.nii'; % mask
    RSAmask_s = 's_Relevant_WholeBrain_RSA_all.nii';  
end

fmri_path = 'Scan_ana';

%  base directory
based         = '';

%  where SPM lives
spm_path      = ''; 

% where RSA toolbox lives
rsa_tool_path = '';

% location of pre-defined algebras
algebra_data_path = '';

mask_path = fullfile(based,fmri_path,'masks');
data_path = fullfile(based,'Data');

SPM_dir = fullfile('stats','RSA_STIMs_vol'); % beta for every stim, normalised

addpath(spm_path)
addpath(genpath(fullfile(based,rsa_tool_path)))

%% defs
subject = {'s02','s03','s04','s05','s06','s08','s09','s10','s11','s12','s13','s14','s15','s16','s17','s18','s19','s20',...
           's21','s22','s23','s24','s25','s26','s27','s28','s29','s30','s31','s32'}; 

distType   = 'correlation';
sl_fname   = 'SL_full_func_r10_v100.mat'; %!!

sub = subject{idx_subj};

fprintf('Running subject %s doing %s\n',sub,type)

outputDir = fullfile(based,fmri_path,sub,'RDMs',['Algebra_',type]);

mkdir(outputDir)

subj_folder = fullfile(based,fmri_path,sub);
files       = dir(fullfile(subj_folder,'sess0*'));

sessions = {files.name};
n_sess   = length(sessions);

which_R = 'R_move_move_diff_physio';

size_R = zeros(length(sessions),1);
for sess_loop = 1:length(sessions)
    load(fullfile(subj_folder,(sessions{sess_loop}),which_R));
    size_R(sess_loop) = size(R,2);
end

if strcmp(type,'nonH')
    load(fullfile(algebra_data_path,'algebra_ontop_nonH_extended.mat'))
    load(fullfile(algebra_data_path,'algebra_beside_nonH_extended.mat'))
    algebra_ontop  = algebra_ontop_nonH;
    algebra_beside = algebra_beside_nonH;
    clear('algebra_ontop_nonH')
    clear('algebra_beside_nonH')
elseif strcmp(type,'H')
    load(fullfile(algebra_data_path,'algebra_ontop_H_extended.mat'))
    load(fullfile(algebra_data_path,'algebra_beside_H_extended.mat'))
    algebra_ontop  = algebra_ontop_H;
    algebra_beside = algebra_beside_H;
    clear('algebra_ontop_H')
    clear('algebra_beside_H')
end

distance_folder = fullfile(based,fmri_path,sub,'RDMs',['Algebra_',type]);

algebra_ontop   = algebra_ontop(~isnan(algebra_ontop(:,end-1)),:);
algebra_beside  = algebra_beside(~isnan(algebra_beside(:,end-1)),:);
algebra_ontop   = algebra_ontop(:,[1:4 7]);
algebra_beside  = algebra_beside(:,[1:4 7]);

n_cond = (size(algebra_ontop,1)+size(algebra_beside,1))*3; % determine length of conditions per session - times three because per equation we are comparing three things

mkdir(distance_folder)

% structure of algebra:
% algebra(:,1:3) = algebra parts
% algebra(:,4) = target
% algebra(:,4) = reference
algebra = [algebra_ontop; algebra_beside];

%% Make model RDM

if mk_RDM
    
    if n_sess == 3
        same_sess_idx = [ones(n_cond) zeros(n_cond) zeros(n_cond);
                         zeros(n_cond) ones(n_cond) zeros(n_cond);
                         zeros(n_cond) zeros(n_cond) ones(n_cond)];
    else
        same_sess_idx = [ones(n_cond) zeros(n_cond);
                         zeros(n_cond) ones(n_cond)];
    end
                 
    same_sess_idx = same_sess_idx(tril(true(length(same_sess_idx)),-1));    
    
    algebra = [algebra_ontop; algebra_beside];

    cd(outputDir)

    V_temp = spm_vol(fullfile(mask_path,RSAmask)); % header details for saving nii files in 2mm structural native space later
    
    % remember image details for later
    V.dim = V_temp.dim;
    V.dt  = V_temp.dt;
    V.mat = V_temp.mat;

    V_temp   = spm_read_vols(V_temp);  

    fprintf('Now defining RDM in %s \n',pwd)

    % i) Define theoretical RDM
    % 1 == dissimilar, 0 = similar (1-r as distance measure)
    
    % r1: mean
    r1 = ones(n_cond*n_sess);    

    % r2: effect of interest: Brick game
    % first n_cond/3 rows/columns are algebra objects
    % second n_cond/3 rows/columns are target objects
    % third n_cond/3 rows/columns are reference objects
    r2_nan = nan(n_cond/3);

    r2_1                                   = r2_nan;
    r2_1(logical(eye(size(r2_1))))         = 1;

    r2_0                                   = r2_nan;
    r2_0(logical(eye(size(r2_0))))         = 0;
    
    r2 = [r2_nan r2_0   r2_1;
          r2_0   r2_nan r2_nan;
          r2_1   r2_nan r2_nan]; 
      
    r2 = repmat(r2,n_sess,n_sess);  
    
%     r2_plot = r2; r2_plot(isnan(r2_plot)) = -1;
%     figure,imagesc(r2_plot)
    
    % r4: diff height/width RDM:
    load(fullfile(data_path,'diff_heightwidth_all.mat'))
    diff_heightwidth(isnan(diff_heightwidth)) = 0;
    diff_heightwidth                          = diff_heightwidth./max(max(diff_heightwidth));
    
    r4_nan = nan(n_cond/3);
    
    size_diff_A1 = nan(size(algebra,1),1);
    size_diff_A2 = nan(size(algebra,1),1);
    size_diff_A3 = nan(size(algebra,1),1);
    
    size_diff_R = nan(size(algebra,1),1);

    
    for idx_eq = 1:size(algebra,1)
        
        size_diff_A1(idx_eq) = diff_heightwidth(algebra(idx_eq,1),algebra(idx_eq,4));
        size_diff_A2(idx_eq) = diff_heightwidth(algebra(idx_eq,2),algebra(idx_eq,4));
        size_diff_A3(idx_eq) = diff_heightwidth(algebra(idx_eq,3),algebra(idx_eq,4));
        
        size_diff_R(idx_eq) = diff_heightwidth(algebra(idx_eq,5),algebra(idx_eq,4));
        
    end
    
    size_diff_A = size_diff_A1 - size_diff_A2 + size_diff_A3;
    
    r4_A                           = r4_nan;
    r4_A(logical(eye(size(r4_A)))) = size_diff_A;
    
    r4_R                           = r4_nan;
    r4_R(logical(eye(size(r4_R)))) = size_diff_R;    
    
    
    % Note that the definition of r4_R (=size_diff_R) assumes that the
    % algebra hypothesis is correct, i.e. B_1-B_2-B_3 = B_target, so we can
    % use B_target to define the size overlap with B_reference.
    % This implies effective control in regions where the silhouette
    % algebra hypothesis is true
    r4 = [r4_nan r4_A   r4_R;
          r4_A   r4_nan r4_nan;
          r4_R   r4_nan r4_nan];

    r4                         = repmat(r4,n_sess,n_sess);
    r4(logical(eye(size(r4)))) = NaN; % set diagonal to NaN 
    
%     r4_plot = r4; r4_plot(isnan(r4_plot)) = -1;
%     figure,imagesc(r4_plot)     
    
    clear('diff_heightwidth')     

    % r5: proportion pixel overlap:
    load(fullfile(data_path,'vis_similarity_move_all.mat'))
    Vis_similarity_move(isnan(Vis_similarity_move)) = 0;
    Vis_similarity_move                             = 1-Vis_similarity_move; % turn into visual DISsimilarity
    
    r5_nan = nan(n_cond/3);
    
    pixel_diff_A1 = nan(size(algebra,1),1);
    pixel_diff_A2 = nan(size(algebra,1),1);
    pixel_diff_A3 = nan(size(algebra,1),1);
    
    pixel_diff_R = nan(size(algebra,1),1);
    
    for idx_eq = 1:size(algebra,1)
        
        pixel_diff_A1(idx_eq) = Vis_similarity_move(algebra(idx_eq,1),algebra(idx_eq,4));
        pixel_diff_A2(idx_eq) = Vis_similarity_move(algebra(idx_eq,2),algebra(idx_eq,4));
        pixel_diff_A3(idx_eq) = Vis_similarity_move(algebra(idx_eq,3),algebra(idx_eq,4));
        
        pixel_diff_R(idx_eq) = Vis_similarity_move(algebra(idx_eq,5),algebra(idx_eq,4));
        
    end
    
    pixel_diff_A = pixel_diff_A1 - pixel_diff_A2 + pixel_diff_A3;
    
    r5_A                           = r5_nan;
    r5_A(logical(eye(size(r5_A)))) = pixel_diff_A;
    
    
    r5_R                           = r5_nan;
    r5_R(logical(eye(size(r5_R)))) = pixel_diff_R;    
    
    % Note that the definition of r5_R (=pixel_diff_R) assumes that the
    % algebra hypothesis is correct, i.e. B_1-B_2-B_3 = B_target, so we can
    % use B_target to define the pixel overlap with B_reference.
    % This implies effective control in regions where the silhouette
    % algebra hypothesis is true
    r5 = [r5_nan r5_A   r5_R;
          r5_A   r5_nan r5_nan;
          r5_R   r5_nan r5_nan]; 

    r5                         = repmat(r5,n_sess,n_sess);
    r5(logical(eye(size(r5)))) = NaN; % set diagonal to NaN   
    
%     r5_plot = r5; r5_plot(isnan(r5_plot)) = -1;
%     figure,imagesc(r5_plot)      

    clear('Vis_similarity_move')
    
    % vectorise RDM:
    r1_vec = r1(tril(true(length(r1)),-1));
    r2_vec = r2(tril(true(length(r2)),-1));
    r4_vec = r4(tril(true(length(r4)),-1)); 
    r5_vec = r5(tril(true(length(r5)),-1)); 

    % ii) define design matrix
    dm = [r1_vec r2_vec r4_vec r5_vec];        
    
    % find relevant data - this is a small fraction of actual data
    idx_use = same_sess_idx==0&~isnan(dm(:,2)); % ols doesn't like NaNs   
    
    dm = dm(idx_use,:);
    
    if any(any(isnan(dm)))
        error('NaNs!')
    end
    
%     dm_plot = dm;
%     dm_plot(isnan(dm_plot)) = -2;
%     figure,imagesc(dm_plot),colorbar    
    
    
end

%% make RDMs

if mk_stats

    cd(outputDir)

    fname = fullfile(based,fmri_path,'searchlights',sub,sl_fname);

    L = load(fname);

    load(fullfile(based,fmri_path,sub,SPM_dir,'SPM'));

    % create conditions
    condition = zeros(1,size(SPM.xX.xKXs.X,2));

    % make sure this is correct:
    X = SPM.xX.xKXs.X(:,1:end-3);
    n_cond_perRun = sum(X(1,:)~=0)-size_R(1); % assumes constant number of conditions per runs - this is 92 (duration: 91), 90 stims and catch + buttonpress (buttonpress only)
    
    % extend over sessions
    if n_sess == 3
        ALGEBRA(:,:,1) = algebra;
        ALGEBRA(:,:,2) = algebra + n_cond_perRun + size_R(1);
        ALGEBRA(:,:,3) = algebra + n_cond_perRun + size_R(1) + n_cond_perRun + size_R(2);
        algebra = [algebra; algebra + n_cond_perRun + size_R(1); algebra + n_cond_perRun + size_R(1) + n_cond_perRun + size_R(2)];
    else
        ALGEBRA(:,:,1) = algebra;
        ALGEBRA(:,:,2) = algebra + n_cond_perRun + size_R(1);
        algebra = [algebra; algebra + n_cond_perRun + size_R(1)];
    end
    
    condition(algebra(:)) = true;
    
    fprintf('Now obtaining empirical RDMs:\n')
    
    % this is based on runSearchlight from https://github.com/rsagroup/rsatoolbox
    VolIn = SPM.xY.VY;
    
    linVox  = unique(cat(2,L.LI{:})');
    [I,J,K] = ind2sub(VolIn(1).dim,linVox);
    
    % Get the data 
    X  = sparse(double(max(linVox)),length(VolIn));
    
    length_vol = length(VolIn);
    
    for idx_vol=1:length_vol

        X(linVox,idx_vol)   = spm_sample_vol(VolIn(idx_vol),double(I),double(J),double(K),0);

        if mod(idx_vol,50)==0
            fprintf('Volume %d of %d done.\n',idx_vol,length(VolIn))
        end
    
    end
    clear I J K; % Keep memory small 
    
    size_sl = size(L.LI,1);
    
    n_conditions = sum(idx_use);     
    
    nii = nan(n_conditions,size_sl); % create a RDM-conditions x #searchlights matrix
    
    % now do some RSA
    for idx_searchlight = 1:size_sl
        
        Y = full(full(X(L.LI{idx_searchlight},:)))';
        
        B = real(rsa.spm.noiseNormalizeBeta(Y,SPM)); % pre-whiten
        B_use = [];
        for idx_sess = 1:size(ALGEBRA,3)

            algebra = ALGEBRA(:,:,idx_sess);

            B1_1 = B(algebra(:,1),:);
            B1_2 = B(algebra(:,2),:);
            B1_3 = B(algebra(:,3),:);
            B1   = B1_1-B1_2+B1_3; % algebra

            B2   = B(algebra(:,4),:); % reference
            
            B3 = B(algebra(:,5),:); % compound game
            
            B_use = [B_use; B1; B2; B3];

        end
        
        % calculate distance
        RDM = pdist(B_use,distType);       % either correlaton or Euclidian on prewhitened betas (Mahalanobis),
        RDM = RDM';
        RDM = RDM(:);
        nii(:,idx_searchlight) = RDM(idx_use); % Arrange the outputs in a vector, as they are written to files  
        
        if mod(idx_searchlight,50)==0
            fprintf('Searchlight %d of %d done.\n',idx_searchlight,size(L.LI,1))
        end
        
    end
    
%     save(fullfile(distance_folder,['nii_',mask]),'nii','-v7.3')

    fprintf('Success RDM part\n')

end

%% do stats
if mk_stats
    
    size_brain = size(V_temp);
    
    fprintf('Size Design Matrix: %d %d\n',size(dm))
    
    fprintf('Size nii: %d %d\n',size(nii))
    
    fprintf('Now running ols: \n')
    
    cmat = eye(size(dm,2));
    
    [c,~,t] = ols(nii,dm,cmat);
    
    % obtain individual coefficient maps:
    brain_stat          = zeros(size_brain);
    brain_stat(L.voxel) = c(2,:);
    
    brain_size          = zeros(size_brain);
    brain_size(L.voxel) = c(3,:);
    
    brain_pixel          = zeros(size_brain);
    brain_pixel(L.voxel) = c(4,:); 
    
    % obtain individual t-maps:
    brain_stat_t          = zeros(size_brain);
    brain_stat_t(L.voxel) = t(2,:);
    
    brain_size_t          = zeros(size_brain);
    brain_size_t(L.voxel) = t(3,:);
    
    brain_pixel_t          = zeros(size_brain);
    brain_pixel_t(L.voxel) = t(4,:);     
    
    save(fullfile(outputDir,'brain_stat'),'brain_stat')
    save(fullfile(outputDir,'brain_size'),'brain_size') 
    save(fullfile(outputDir,'brain_pixel'),'brain_pixel') 
    save(fullfile(outputDir,'brain_stat_t'),'brain_stat_t')
    save(fullfile(outputDir,'brain_size_t'),'brain_size_t') 
    save(fullfile(outputDir,'brain_pixel_t'),'brain_pixel_t')     
    
    V.fname = ['Stats.nii'];
    spm_write_vol(V,brain_stat)
    V.fname = ['Size.nii'];
    spm_write_vol(V,brain_size)
    V.fname = ['Pixel.nii'];
    spm_write_vol(V,brain_pixel)
    V.fname = ['Stats_t.nii'];
    spm_write_vol(V,brain_stat_t)
    V.fname = ['Size_t.nii'];
    spm_write_vol(V,brain_size_t)
    V.fname = ['Pixel_t.nii'];
    spm_write_vol(V,brain_pixel_t)
    
    fprintf('Success stats part\n')

end

%% average over stats:

if mk_av
    
    smthSurf = 5;
    
    cd(outputDir)
    
    files_smooth = {['Stats.nii'],['Size.nii'],['Pixel.nii'],...
                    ['Stats_t.nii'],['Size_t.nii'],['Pixel_t.nii']};  
    
    % smooth SPM style
    clear jobs; 

    fprintf('Now smoothing SPM style\n')

    jobs{1}.spatial{1}.smooth.data                          = files_smooth';
    jobs{1}.spatial{1}.smooth.fwhm                          = [smthSurf smthSurf smthSurf];
    jobs{1}.spatial{1}.smooth.dtype                         = 0; 
    jobs{1}.spatial{1}.smooth.im                            = 0;
    jobs{1}.spatial{1}.smooth.prefix                        = 's_';

    spm_jobman('run',jobs);

    clear jobs;        

    fprintf('Smoothing SPM done\n')
    
    % now correct for smoothing at borders
    RSAmask   = spm_vol(fullfile(mask_path,RSAmask));
    RSAmask   = spm_read_vols(RSAmask);
    
    RSAmask_s = spm_vol(fullfile(mask_path,RSAmask_s));
    RSAmask_s = spm_read_vols(RSAmask_s);
    
    for idx_file=1:length(files_smooth)
        
        mask = spm_vol(['s_',files_smooth{idx_file}]);
        
        V.dim   = mask.dim;
        V.dt    = mask.dt;
        V.mat   = mask.mat;
        V.fname = ['s_corrected_',files_smooth{idx_file}];
            
        V_new = spm_read_vols(mask);
        
        V_new = V_new./RSAmask_s;
        V_new = V_new.*RSAmask;
        
        spm_write_vol(V,V_new)  
        
        clear('V')
        clear('V_new')
        
    end
    
    clear('RSAmask')
    clear('RSAmask_s')
    
end

fprintf('All done subject %s mate.\n',sub)


end
