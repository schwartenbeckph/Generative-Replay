% Second-level analysis of Lego
% PS 01/06/18

function []=Stats_SecondLevel_Lego(name_analyses,name_analysis_second)

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

    cd(based)

    addpath(spm_path)

    addpath(fullfile(based,fmri_path))
    addpath(fullfile(based,fmri_path,script_path))
    addpath(fullfile(based,behav_path))

    subject = {'s00','s01','s02','s03','s04','s05','s06','s08','s09','s10','s11','s12','s13','s14','s15','s16','s17','s18','s19','s20',...
               's21','s22','s23','s24','s25','s26','s27','s28','s29','s30','s31','s32'}; 

    pilot = 2; % first two subjects are pilot!

    for idx_ana = 1:length(name_analyses)

        name_analysis = name_analyses{idx_ana};
        
        if nargin<2
            name_analysis_second = name_analysis;
        end

        second_level_path = fullfile([based,fmri_path,'\Second_level\',name_analysis_second]);
        
        statpath = fullfile('\stats\',name_analysis); % where your first levels live

        % load conrast names
        cd(second_level_path)
        load('contrast_indicator')
        contrasts = contrast_indicator.name;
        cd(based)

        for contrast_loop=1:length(contrasts)

            newdir = fullfile(second_level_path,(contrasts{contrast_loop}));
            mkdir(newdir)

            % delete previous model if existed
            clear matlabbatch

            matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_move.files = cellstr(spm_select('FPList', newdir));
            matlabbatch{1}.cfg_basicio.file_dir.file_ops.file_move.action.delete = false;

            spm_jobman('run', matlabbatch)

            fprintf('------------ deleted previous files for %s ------------\n',(contrasts{contrast_loop}))

            clear matlabbatch

            % specify new contrast
            clear matlabbatch 

            matlabbatch{1}.spm.stats.factorial_design.dir = {newdir};

            matlabbatch{1}.spm.stats.factorial_design.des.t1.scans={};

            for subj_loop = (pilot+1):size(subject,2)
                if contrast_loop<10
                    con_file = cellstr(spm_select('FPList', [based,fmri_path,'\',(subject{subj_loop}),statpath], ['^con_000',num2str(contrast_loop),'.nii$']));
                else
                    con_file = cellstr(spm_select('FPList', [based,fmri_path,'\',(subject{subj_loop}),statpath], ['^con_00',num2str(contrast_loop),'.nii$']));
                end
                matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = [matlabbatch{1}.spm.stats.factorial_design.des.t1.scans; con_file];
            end

            matlabbatch{1}.spm.stats.factorial_design.cov                    = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none     = 1;     % threshold masking
            matlabbatch{1}.spm.stats.factorial_design.masking.im             = 1;     % implicit masking
            matlabbatch{1}.spm.stats.factorial_design.masking.em             = {''};  % explicit masking
            matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit         = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm        = 1;

            spm_jobman('run', matlabbatch)
            fprintf('------------ specify model done for contrast %d ------------\n',contrast_loop)


            clear matlabbatch

            currdir = newdir;

            matlabbatch{1}.spm.stats.fmri_est.spmmat           = {[currdir,'/SPM.mat']};
            matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;

            spm_jobman('run', matlabbatch)
            fprintf('------------ estimate done for contrast %d ------------\n',contrast_loop)    


            clear matlabbatch

            currdir                                              = newdir;    
            matlabbatch{1}.spm.stats.con.spmmat                  = {[currdir,'/SPM.mat']};
            matlabbatch{1}.spm.stats.con.consess{1}.tcon.name    = char(contrasts(contrast_loop));
            matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1;
            matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
            matlabbatch{1}.spm.stats.con.consess{2}.tcon.name    = char(contrasts(contrast_loop));
            matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = -1;
            matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';        
            matlabbatch{1}.spm.stats.con.delete                  = 1;

            spm_jobman('run', matlabbatch)
            fprintf('------------ contrasts done for contrast %d ------------\n',contrast_loop)  


        end

    end

    fprintf('------------ all done 2nd level. ------------\n')

end