% first-level suppression  analysis Lego
% Idea: have one regressor for all elements (nH and H separately) 
% PS 23/06/20

clear
% close all

% base directory
based = '';

% spm directory
spm_path    = '';

% directory of fMRI data
fmri_path   = '';

% directory of fMRI analysis code
script_path = '';

% directory of behavioural results
behav_path  = '';

addpath(spm_path)

addpath(fullfile(based,fmri_path))
addpath(genpath(fullfile(based,script_path)))
addpath(fullfile(based,behav_path))

% spm fmri

subject = {'s00','s01','s02','s03','s04','s05','s06','s08','s09','s10','s11','s12','s13','s14','s15','s16','s17','s18','s19','s20',...
           's21','s22','s23','s24','s25','s26','s27','s28','s29','s30','s31','s32'}; 
       
exclude_sess         = zeros(length(subject),3);
exclude_sess(3,3)    = 1; % Translation greater than 1.5mm for subject s02 in session 3 with a max value of 1.8178
exclude_sess(14,1)   = 1; % Translation greater than 1.5mm for subject s14 in session 1 with a max value of 1.552
exclude_sess(16,2)   = 1; % Translation greater than 1.5mm for subject s16 in session 2 with a max value of 2.6782
exclude_sess(25,2:3) = 1;
exclude_sess(30,3)   = 1; % Translation greater than 1.5mm for subject s30 in session 3 with a max value of 2.1067
       
pilot = 2; % first two subjects are pilot!

n_dummy = 6;

corr_thresh = .6;

do_specify   = false;
do_estimate  = false;
do_contrasts = false;

name_analysis = 'fMRI_SummaryStats_Univariate';

load('regressors_fmri')

which_R = 'R_move_move_diff_physio';

do_anyway = true;
do_delete = true;
plot_contrast_first_sub   = false;
plot_regressors_first_sub = false;

r_all = [];

for subj_loop = (pilot+1):length(subject)

    subj_folder = fullfile(based,fmri_path,(subject{subj_loop}));
    
    if ~exist(fullfile(subj_folder,'stats',name_analysis,'spmT_0009.nii'),'file') || do_anyway

        files = dir(fullfile(subj_folder,'sess0*'));

        sessions = {files.name};
        
        behav_files = dir(fullfile(based,behav_path,(subject{subj_loop}),'Tetris_Scanner_subj*.mat'));
        behav_files = {behav_files.name};
        
        behav_files = behav_files(~logical(exclude_sess(subj_loop,1:length(sessions))));
        
        sessions = sessions(~logical(exclude_sess(subj_loop,1:length(sessions))));

        newdir = fullfile(subj_folder,'stats',name_analysis);

        if do_specify 

            % delete previous model
            clear matlabbatch
            
            if do_delete

                matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_move.files = cellstr(spm_select('FPList', newdir));
                matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_move.action.delete = false;

                spm_jobman('run', matlabbatch)

                fprintf('------------ deleted previous files for %s ------------\n',(subject{subj_loop}))

                clear matlabbatch
                
            end

            for sess_loop = 1:length(sessions)
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
                % obtain regressors and plot some stuff: (don't run this for all subjects)
                onsets_allobjects = regressors(subj_loop,sess_loop).onsets_allobjects;
                
                stim_code = regressors(subj_loop,sess_loop).stim_code;
                
                % control regressors
                pixel_suppression = regressors(subj_loop,sess_loop).pixel_suppression;
                size_suppression  = regressors(subj_loop,sess_loop).size_suppression;
                N_elements        = regressors(subj_loop,sess_loop).N_elements;
                Diff_element      = regressors(subj_loop,sess_loop).Diff_element;
                
                % relation suppression regressors
                relation_nH  = regressors(subj_loop,sess_loop).relation_nH;
                relation_H   = regressors(subj_loop,sess_loop).relation_H;
                relation_HnH = regressors(subj_loop,sess_loop).relation_HnH;
                
                % element suppression regressors
                el_onlyHSol_split = regressors(subj_loop,sess_loop).el_onlyHSol_split; % 
                el_H_split        = regressors(subj_loop,sess_loop).el_H_split;
                
                idx_allEl_noH   = [1 1; 1 2; 1 3; 1 4; 2 2; 2 3; 2 4; 3 3; 3 4; 4 4];
                idx_allEl_noH   = unique([idx_allEl_noH; fliplr(idx_allEl_noH)],'rows');
                
                idx_allEl_onlyH = [1 5; 1 6; 2 5; 2 6; 3 5; 3 6; 4 5; 4 6];
                idx_allEl_onlyH = unique([idx_allEl_onlyH; fliplr(idx_allEl_onlyH)],'rows');
                
                idx_allEl_H = [2 5; 2 6; 5 5; 5 6; 6 6];
                idx_allEl_H = unique([idx_allEl_H; fliplr(idx_allEl_H)],'rows');


                if subj_loop==(pilot+1) && sess_loop==1 && plot_contrast_first_sub
                    allEl_noH   = mk_el_suppresion_contrasts(idx_allEl_noH,el_onlyHSol_split,true,stim_code);
                    allEl_onlyH = mk_el_suppresion_contrasts(idx_allEl_onlyH,el_onlyHSol_split,true,stim_code);
                    allEl_H     = mk_el_suppresion_contrasts(idx_allEl_H,el_H_split,true,stim_code);
                    
                    plot_stims(stim_code,relation_nH,true)
                    plot_stims(stim_code,relation_H,true)
                    plot_stims(stim_code,relation_HnH,true)
                else
                    allEl_noH   = mk_el_suppresion_contrasts(idx_allEl_noH,el_onlyHSol_split);
                    allEl_onlyH = mk_el_suppresion_contrasts(idx_allEl_onlyH,el_onlyHSol_split);
                    allEl_H     = mk_el_suppresion_contrasts(idx_allEl_H,el_H_split); 
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                pixel_suppression = mk_equalSpace_noZero_noNaN(pixel_suppression,-1,1);
                
                size_suppression = mk_equalSpace_noZero_noNaN(size_suppression,-1,1);

                Diff_element(Diff_element==0) = eps; % check if this works
                
                allEl_noH   = mk_equalSpace_noZero_noNaN(allEl_noH,-1,1);
                allEl_onlyH = mk_equalSpace_noZero_noNaN(allEl_onlyH,-1,1);
                allEl_H     = mk_equalSpace_noZero_noNaN(allEl_H,-1,1);
                
                mkdir(newdir)

                matlabbatch{1}.spm.stats.fmri_spec.dir            = {newdir}; 
                matlabbatch{1}.spm.stats.fmri_spec.timing.units   = 'secs';
                matlabbatch{1}.spm.stats.fmri_spec.timing.RT      = 1.45;
                matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t  = 16;
                matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

                allfiles = cellstr(spm_select('FPList', [subj_folder,filesep,(sessions{sess_loop}),filesep], '^S6wubf.*.nii$')); % new prepro

                allfiles = allfiles(n_dummy+1:end);

                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).scans = allfiles; % select all scans excluding first six dummy scans

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % obtain betas per stimulus
                condition = 1;

                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).cond(condition).name     = 'all_stims'; 
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).cond(condition).onset    = onsets_allobjects;
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).cond(condition).duration = 0;
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).cond(condition).tmod     = 0;

                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).cond(condition).pmod(1).name  = 'pixel_supp';
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).cond(condition).pmod(1).param = pixel_suppression;
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).cond(condition).pmod(1).poly  = 1;
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).cond(condition).pmod(2).name  = 'size_supp';
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).cond(condition).pmod(2).param = size_suppression;
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).cond(condition).pmod(2).poly  = 1;  
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).cond(condition).pmod(3).name  = 'Number_elements';
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).cond(condition).pmod(3).param = N_elements';
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).cond(condition).pmod(3).poly  = 1; 
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).cond(condition).pmod(4).name  = 'DiffNumber_elements';
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).cond(condition).pmod(4).param = Diff_element';
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).cond(condition).pmod(4).poly  = 1; 
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).cond(condition).pmod(5).name  = 'Relation_nH';
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).cond(condition).pmod(5).param = relation_nH;
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).cond(condition).pmod(5).poly  = 1;  
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).cond(condition).pmod(6).name  = 'Relation_H';
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).cond(condition).pmod(6).param = relation_H;
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).cond(condition).pmod(6).poly  = 1;  
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).cond(condition).pmod(7).name  = 'Relation_HnH';
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).cond(condition).pmod(7).param = relation_HnH;
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).cond(condition).pmod(7).poly  = 1;
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).cond(condition).pmod(8).name  = 'allEl_noH';
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).cond(condition).pmod(8).param = allEl_noH;
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).cond(condition).pmod(8).poly  = 1;
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).cond(condition).pmod(9).name  = 'allEl_onlyH';
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).cond(condition).pmod(9).param = allEl_onlyH;
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).cond(condition).pmod(9).poly  = 1;                
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).cond(condition).pmod(10).name  = 'allEl_H';
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).cond(condition).pmod(10).param = allEl_H;
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).cond(condition).pmod(10).poly  = 1;                 

                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).cond(condition).orth          = 0;
                
                r = corr([onsets_allobjects' ...
                          pixel_suppression size_suppression N_elements' Diff_element'...
                          relation_nH relation_H relation_HnH ...
                          allEl_noH allEl_onlyH allEl_H]);
                r(r==diag(r)) = 0;
                
                if any(abs(r(:))>corr_thresh)
                    fprintf('##### Correlations are high sub %d sess %d. #####\n',subj_loop,sess_loop)
                end
                
                r_all = [r_all r(:)];

                if subj_loop==(pilot+1) && sess_loop==1 && plot_regressors_first_sub
                    figure,imagesc([onsets_allobjects' ...
                                    pixel_suppression size_suppression N_elements' Diff_element'...
                                    relation_nH relation_H relation_HnH ...
                                    allEl_noH allEl_onlyH allEl_H])
                    colorbar
                end

                condition = condition + 1;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).multi     = {''};
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).regress   = struct('name', {}, 'val', {});
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).multi_reg = cellstr(spm_select('FPList', fullfile(subj_folder,(sessions{sess_loop})), which_R));
                matlabbatch{1}.spm.stats.fmri_spec.sess(sess_loop).hpf       = 128; 

            end

            matlabbatch{1}.spm.stats.fmri_spec.fact             = struct('name', {}, 'levels', {});
            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
            matlabbatch{1}.spm.stats.fmri_spec.volt             = 1;
            matlabbatch{1}.spm.stats.fmri_spec.global           = 'None';

        %     matlabbatch{1}.spm.stats.fmri_spec.mthresh          = 0;
            matlabbatch{1}.spm.stats.fmri_spec.mthresh          = -inf;

        %     matlabbatch{1}.spm.stats.fmri_spec.mask             = {''};
            matlabbatch{1}.spm.stats.fmri_spec.mask             = cellstr(spm_select('FPList', fullfile(spm_path,'tpm'), 'mask_ICV.nii$'));
            matlabbatch{1}.spm.stats.fmri_spec.cvi              = 'AR(1)';

            spm_jobman('run', matlabbatch)

            fprintf('------------ 1st level %s done ------------\n',(subject{subj_loop}))

        end 

        %% estimation
        if do_estimate

            clear matlabbatch

            currdir = newdir;

            matlabbatch{1}.spm.stats.fmri_est.spmmat           = {fullfile(currdir,'SPM.mat')};
            matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;

            spm_jobman('run', matlabbatch)

            fprintf('------------estimate %s done ------------\n',(subject{subj_loop}))

        end
        %% contrasts

        if do_contrasts

            addpath(newdir)

            load('regressors_fmri')    

            size_R = zeros(length(sessions),1);

            for sess_loop = 1:length(sessions)        
                load(fullfile(subj_folder,(sessions{sess_loop}),which_R));        
                size_R(sess_loop) = size(R,2);        
            end
            
            on  = 1;
            off = 0;
            
            double_on = [1 1];
            double_off = [0 0];


            clear matlabbatch

            currdir = newdir;

            matlabbatch{1}.spm.stats.con.spmmat = {fullfile(currdir,'SPM.mat')};
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % make contrasts:           
            n_cond = 11;
            % 1)  Objects on screen
            % 2)  Pixel
            % 3)  Size
            % 4)  # Elements in object
            % 5)  Diff # Elements in object
            % 6)  Relation nH
            % 7)  Relation H
            % 8)  Relation HnH
            % 9)  allEl_noH
            % 10) allEl_onlyH
            % 11) allEl_H
            
            names_con   = {'All_objects' 'Pixel' 'Size' 'Visual_Sup' 'Number_El' 'Diff_El' ...
                           'Relation_nH' 'Relation_H' 'Relation_HnH' 'Relation_nH_and_H' 'Relation_all' 'Relation_nHsmall_minus_H' 'Relation_nHall_minus_H' ...
                           'allEl_noH' 'allEl_onlyH' 'allEl_H' 'allEl_noH_and_H' 'allEl_all' 'allEl_noHsmall_mimus_H' 'allEl_noHall_mimus_H'}; 
            
            weights_all = zeros(n_cond,length(names_con));
            
            contrast    = 1;
            
            weights_all([1],contrast)       = on;   contrast = contrast+1; % All_objects
            weights_all([2],contrast)       = on;   contrast = contrast+1; % Pixel
            weights_all([3],contrast)       = on;   contrast = contrast+1; % Size
            weights_all([2 3],contrast)     = on;   contrast = contrast+1; % Visual_Sup
            weights_all([4],contrast)       = on;   contrast = contrast+1; % Number_El
            weights_all([5],contrast)       = on;   contrast = contrast+1; % Diff_El
            weights_all([6],contrast)       = on;   contrast = contrast+1; % Relation_nH
            weights_all([7],contrast)       = on;   contrast = contrast+1; % Relation_H
            weights_all([8],contrast)       = on;   contrast = contrast+1; % Relation_HnH
            weights_all([6 8],contrast)     = on;   contrast = contrast+1; % Relation_nH_and_H
            weights_all([6 7 8],contrast)   = on;   contrast = contrast+1; % Relation_all
            weights_all([6],contrast)       = on;   weights_all([8],contrast) = -on; contrast = contrast+1; % Relation_nHsmall_minus_H
            weights_all([6 7],contrast)     = on/2; weights_all([8],contrast) = -on; contrast = contrast+1; % Relation_nHall_minus_H
            weights_all([9],contrast)       = on;   contrast = contrast+1; % allEl_noH
            weights_all([10],contrast)      = on;   contrast = contrast+1; % allEl_onlyH
            weights_all([11],contrast)      = on;   contrast = contrast+1; % allEl_H
            weights_all([9 11],contrast)    = on;   contrast = contrast+1; % allEl_noH_and_H
            weights_all([9 10 11],contrast) = on;   contrast = contrast+1; % allEl_all
            weights_all([9],contrast)       = on;   weights_all([11],contrast) = -on; contrast = contrast+1; % allEl_noHsmall_mimus_H
            weights_all([9 10],contrast)    = on/2; weights_all([11],contrast) = -on; contrast = contrast+1; % allEl_noHall_mimus_H
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                    
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
            contrast = 1;
            for idx_con = 1:size(weights_all,2)
                
                contrast_indicator.name{idx_con} = names_con{idx_con};
                
                matlabbatch{1}.spm.stats.con.consess{idx_con}.tcon.name = contrast_indicator.name{idx_con}; 

                weights = [];

                for sess_loop = 1:length(sessions)

                    weights = [weights, ...
                               weights_all(:,idx_con)', ... % param_on controls if sessions are not same size!                               
                               zeros(1,size_R(sess_loop))];
                end

                matlabbatch{1}.spm.stats.con.consess{idx_con}.tcon.weights = weights;
                matlabbatch{1}.spm.stats.con.consess{idx_con}.tcon.sessrep = 'none'; 
            
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                    

            matlabbatch{1}.spm.stats.con.delete = 1;

            spm_jobman('run', matlabbatch)

            % save contrast names for second level ana
            second_level_path = fullfile(based,fmri_path,'Second_level',name_analysis);
            mkdir(second_level_path)
            save(fullfile(second_level_path,'contrast_indicator'),'contrast_indicator')

            fprintf('------------ contrasts %s done ------------\n',(subject{subj_loop}))
        end

    end
    
    fprintf('####### all done 1st level subject %s.#######\n',(subject{subj_loop}))        

end


fprintf('------------ all done 1st level. ------------\n')


%% now run second level

Stats_SecondLevel_Lego({name_analysis})


%% plot correlation matrix

r_mean = mk_meanCorr_PS(r_all,2);
figure,
imagesc(reshape(r_mean,sqrt(length(r_mean)),sqrt(length(r_mean)))),colorbar,caxis([-1 1])


