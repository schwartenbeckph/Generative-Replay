% function to train classifiers on data
% 
% INPUT:
% data and nulldata (nulldata can be empty)
% labStm: stim labels
% L1_use: L1 penalty
% optimise_null:  find amount of data needed to bring corr between classifiers below threshold_corr
% threshold_corr: threshold for defining amount of nulldata
% n_null: amount of null data to consider
% Loc_11: do one against one classifiers - takes longer
% 
% OUTPUT:
% betas and intercept from logistic regression from data (and nullData) to labStm
% amount_null: how much null data to include to bring correlation below some threshold

function [betas, intercepts, amount_null] = train_classifiers(data,nullData,labStm,L1_use,optimise_null,threshold_corr,n_null,Loc_11)

    % check if classes are too unbalanced:
    figure(999)
    p = histogram(labStm);
    if max(p.Values)/length(labStm)-min(p.Values)/length(labStm) > 0.1
        take_vals = min(p.Values);
        take_trials = [];
        for idx_val=1:length(unique(labStm))
            trials_class = find(labStm==idx_val);
            take_trials = [take_trials randsample(trials_class,take_vals)];
        end
        take_trials = sort(take_trials);
        data     = data(take_trials,:);
        nullData = nullData(take_trials,:);
        labStm   = labStm(take_trials);
    end    
    close 999

    if nargin<4
        optimise_null  = false;
        threshold_corr = NaN;
        n_null         = NaN;
        Loc_11         = false;
    end
    
    if isempty(optimise_null)
        optimise_null = false;
    end
    if isempty(threshold_corr)
        threshold_corr = 0.2;
    end
    if isempty(n_null)
        n_null = [10:10:size(nullData,1) size(nullData,1)];
    end
    if isempty(Loc_11)
        Loc_11 = false;
    end    

    if optimise_null

        done = false; 

        for idx_Nnull=1:length(n_null) % find #nulldata that minimises corr of weights for a given threshold

            if ~done

                idx_null_use     = randsample(1:size(nullData,1),n_null(idx_Nnull));
                nullDataLI_use   = nullData(idx_null_use,:);

                if ~Loc_11

                    % obtain one against all classifier
                    for iC=1:length(unique(labStm))

                        [betas(:,iC), fitInfo] = lassoglm([data; nullDataLI_use], [(labStm == iC)'; zeros(size(nullDataLI_use,1),1)], 'binomial', ...
                                                            'Alpha', 1, 'Lambda', L1_use, 'Standardize', false);  % matlab's lassoglm, Alpha = 1 represents lasso regression with Lambda as regularisation                                                  

                        intercepts(iC) = fitInfo.Intercept;

                    end 

                    % obtain one against one classifier
                elseif Loc_11

                    for iC=1:length(allconds)

                        for iC2=1:length(allconds)

                            if iC~=iC2
                                idx_use = (labStm==iC) | (labStm==iC2);
                            else
                                idx_use = true(length(labStm),1);
                            end

                            dataLI_use = data(idx_use,:);

                            labStm_use = labStm(idx_use);

                            [betas(:,iC,iC2), fitInfo] = lassoglm([dataLI_use; nullDataLI_use], [(labStm_use == iC)'; zeros(size(nullDataLI_use,1),1)], 'binomial', ...
                                                                'Alpha', 1, 'Lambda', L1_use, 'Standardize', false);  % matlab's lassoglm, Alpha = 1 represents lasso regression with Lambda as regularisation                                                  

                            intercepts(iC,iC2) = fitInfo.Intercept;

                        end

                    end    

                end

                if Loc_11
                    r = corr([squeeze(betas(:,1,[2 3 4])) squeeze(betas(:,2,[1 3 4])) ...
                              squeeze(betas(:,3,[1 2 4])) squeeze(betas(:,4,[1 2 3]))]);
                else
                    r = corr(betas);
                end

                if all(abs(r(tril(true(length(r)),-1)))<threshold_corr) || idx_Nnull==length(n_null)
                    done = true;
                    fprintf('need %d nulldata.\n',n_null(idx_Nnull))
                    amount_null = n_null(idx_Nnull);
                end                 

            end                       

        end

    else

        amount_null = size(nullData,1);

        if ~Loc_11
            % obtain one against all classifier
            for iC=1:length(unique(labStm))

                [betas(:,iC), fitInfo] = lassoglm([data; nullData], [(labStm == iC)'; zeros(size(nullData,1),1)], 'binomial', ...
                                                    'Alpha', 1, 'Lambda', L1_use, 'Standardize', false);  % matlab's lassoglm, Alpha = 1 represents lasso regression with Lambda as regularisation                                                  

                intercepts(iC) = fitInfo.Intercept;

            end

        % obtain one against all classifier
        elseif Loc_11

            for iC=1:length(unique(labStm))

                for iC2=1:length(unique(labStm))

                    if iC~=iC2
                        idx_use = (labStm==iC) | (labStm==iC2);
                    else
                        idx_use = true(length(labStm),1);
                    end

                    dataLI_use = data(idx_use,:);

                    labStm_use = labStm(idx_use);

                    [betas(:,iC,iC2), fitInfo] = lassoglm([dataLI_use; nullData], [(labStm_use == iC)'; zeros(size(nullData,1),1)], 'binomial', ...
                                                        'Alpha', 1, 'Lambda', L1_use, 'Standardize', false);  % matlab's lassoglm, Alpha = 1 represents lasso regression with Lambda as regularisation                                                  

                    intercepts(iC,iC2) = fitInfo.Intercept;

                end

            end 

        end

    end

end