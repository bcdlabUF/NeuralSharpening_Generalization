% ========================================================================
%%              Bootstrapped Bayesian Model Testing
                  % SHARPENING & GENERALIZATION
% ======================================================================== AK, MRB, Mar 2024
% ==================================== Last updated for data sharing: October 2026, MRB
% NOTE: due to the boostrapping, distributions and resulting BF values can be slightly different
% if you run this multiple times. Be sure to save the innerprod_null (null distribution) and
% innerprod_tested (tested distribution)

load('Boylan_etal_DCN_bayesianModels.mat')
% dataAll dimensions: 129 x 51 x 4 x 3
    % 1. electrodes
    % 2. all infant participants
    % 3. stimuli: TS, GS1, GS2, GS3
    % 4. label conditions: Individual, Category, Unlabeled
% -----------------------------------------------------------------
% make bootstrap distributions of inner products, by age and label
% -----------------------------------------------------------------
% starting parameters
    ndraws = 5000;  % number of draws/iterations for bootstrapping
    % sharpening vs generalization; weights across stimuli
    weights = [1.75 -2.75 -0.25 1.25;
               1.75 1.25 -0.25 -2.75]; 
% load data
    data2glm = allData; % data2glm => 129 x 51 (sub) x 4 (TS,GS1-3) x 3 (label)
ageMos = [6 9 12];
% run bootstrapping
    for age = 1:3
        whichSubj = find(ages_alphabetical==ageMos(age)); % ID which subj are in current age group
        for draw = 1:ndraws
            if draw./100 == round(draw./100), fprintf(['draw: ' num2str(draw) ' of ' num2str(ndraws) '\n']), end % progress indicator
            nSubjs = size(whichSubj,1);
            for label = 1:3
                data = data2glm(:,whichSubj,:,label);
                for model = 1:2
                % get distributions: inner product of group means (px drawn with replacement) and model weights
                % for tested models
                    innerprod_tested_byLabel(:,draw,label,model,age) = squeeze(mean(data(:,randi(nSubjs,1,nSubjs),:),2,'omitnan')) * weights(model,:)';
                 % for null model
                    innerprod_null_byLabel(:,draw,label,model,age) = squeeze(mean(data(:,randi(nSubjs,1,nSubjs),:),2,'omitnan')) * weights(model,randperm(4))'; 
                end % model
            end % label
        end % draw
    end % age group
fprintf('Done with Bootstrapping!\n')
% get bayes factors for model vs null distributions
    for age = 1:3
        for model = 1:2
            for label = 1:3
                for elec = 1:129
                    BF10_models_byLabel(elec,label,model,age) = bootstrap2BF_z(innerprod_tested_byLabel(elec,:,label,model,age),innerprod_null_byLabel(elec,:,label,model,age), 0);
                    BF01_models_byLabel(elec,label,model,age) = bootstrap2BF_z(innerprod_null_byLabel(elec,:,label,model,age), innerprod_tested_byLabel(elec,:,label,model,age),0);
                end % elec
            end % label
        end % model
    end % age
fprintf('Done with bootstrap 2 BF!\n')
clearvars -except ages_alphabetical BF10_models_byLabel BF01_models_byLabel allData innerprod_null_byLabel innerprod_tested_byLabel ndraws weights 

%% -----------------------------------------------------------------------
% make bootstrap distributions of inner products, by age, avg over label
% ------------------------------------------------------------------------
% starting parameters
    ndraws = 5000; 
    weights = [1.75 -2.75 -0.25 1.25;1.75 1.25 -0.25 -2.75]; 
% load data
    data2glm = mean(allData,4,'omitnan'); % => 129 x 51 (sub) x 4 (TS,GS1-3)
ageMos = [6 9 12];
for age = 1:3
    whichSubj = find(ages_alphabetical_noOut==ageMos(age));
    for draw = 1:ndraws
        if draw./100 == round(draw./100), fprintf(['draw: ' num2str(draw) ' of ' num2str(ndraws) '\n']), end % progress indicator
        nSubjs = size(whichSubj,1);
        data = data2glm(:,whichSubj,:);
        for model = 1:2
        % for tested models
            innerprod_tested_avgLab(:,draw,model,age) = squeeze(nanmean(data(:,randi(nSubjs,1,nSubjs),:),2)) * weights(model,:)';
         % for null model
            innerprod_null_avgLab(:,draw,model,age) = squeeze(nanmean(data(:,randi(nSubjs,1,nSubjs),:), 2)) * weights(model,randperm(length(weights)))'; % => BF_models2
        end
    end % draw
end % age
fprintf('Done with Bootstrapping!\n')
% get bayes factors for model vs null distributions
for age = 1:3
    for model = 1:2
        for elec = 1:129
            BF10_models_avgLab(elec,model,age) = bootstrap2BF_z(innerprod_tested_avgLab(elec,:,model,age),innerprod_null_avgLab(elec,:,model,age), 0);
            BF01_models_avgLab(elec,model,age) = bootstrap2BF_z(innerprod_null_avgLab(elec,:,model,age), innerprod_tested_avgLab(elec,:,model,age),0);
        end % elec
    end % model
end % age
fprintf('Done with bootstrap 2 BF!\n')
clear age ageMos data data2glm data2glm_tmp draw elec label model nSubjs sub whichSubj
fprintf('Done with Section 2: Bootstrapped Bayesian Model Testing.\n')
% ========================================================================
%%                           PLOTTING
% ========================================================================
labelNames = {'Individual','Category','Unlabeled'};
modelNames = {'Sharpening','Generalization'};
ageNames = {'6-month-olds','9-month-olds','12-month-olds'};
load locsEEGLAB129HCL
%% -----------------------------------------------------------------------
% tested models vs null model
% ------------------------------------------------------------------------
% flat maps with colorful suprathres elec markers, Figures 5, S1, S2
for age = 1:3
    for model = 1:2
        for label = 1:3
            tempElecs = find(BF10_models_byLabel(:,label,model,age) >= 3);
            sigElecIdx_sharpGen3{label,model,age} = tempElecs;
            minMax_sharpGen(label,1,model,age) = min(BF10_models_byLabel(tempElecs,label,model,age));
            minMax_sharpGen(label,2,model,age) = max(BF10_models_byLabel(tempElecs,label,model,age));
        end % label
        tempElecs2 = find(BF10_models_avgLab(:,model,age) >= 3);
        sigElecIdx_sharpGen3_avgLab{model,age} = tempElecs2;
        minMax_sharpGen_avgLab(1,model,age) = min(BF10_models_avgLab(tempElecs2,model,age));
        minMax_sharpGen_avgLab(2,model,age) = max(BF10_models_avgLab(tempElecs2,model,age));
    end % model
end % age

% multiple BF colors
for age = 1:3
    for model = 1:2
    % get colors right
        colors4topo = turbo(12);colors4topo = colors4topo([8 9 11],:);
    % assign suprathres elecs their color
        dataBF = BF10_models_byLabel(:,:,model,age);
        for label = 1:3
            colorVec_tmp =[];
            BFvalues = dataBF(sigElecIdx_sharpGen3{label,model,age},label);
            for x = 1:size(BFvalues,1)
                if BFvalues(x) >=3 && BFvalues(x) < 10, colorVec_tmp(x,:) = colors4topo(1,:);
                elseif BFvalues(x) >=10 && BFvalues(x) < 20, colorVec_tmp(x,:) = colors4topo(2,:);
                elseif BFvalues(x) >=20, colorVec_tmp(x,:) = colors4topo(3,:);
                end
            end % suprathres elecs
            colorVec{label,model,age} = colorVec_tmp;
        end % label
        dataBF = BF10_models_avgLab(:,model,age);
        colorVec_tmp = [];
        BFvalues = dataBF(sigElecIdx_sharpGen3_avgLab{model,age});
        for x = 1:size(BFvalues,1)
            if BFvalues(x) >=3 && BFvalues(x) < 10, colorVec_tmp(x,:) = colors4topo(1,:);
            elseif BFvalues(x) >=10 && BFvalues(x) < 20, colorVec_tmp(x,:) = colors4topo(2,:);
            elseif BFvalues(x) >=20, colorVec_tmp(x,:) = colors4topo(3,:);
            end
        end % suprathres elecs
        colorVec_avgLab{model,age} = colorVec_tmp;
    end % model
end % age
% plot per age and label
for model = 1:2
    c = 1;
    for age = 1:3
        for label = 1:3
            figure(model)
            subplot(3,3,c)
            topoplotMB(zeros(129,1),locsEEGLAB129HCL,'style','straight','electrodes','on', 'maplimits',[0 7],...
                'emarker2',{sigElecIdx_sharpGen3{label,model,age},'o',colorVec{label,model,age},10,1}, ...
                'hcolor','k','emsize',10,'ecolor','k','colormap',white); 
            title([modelNames(model) ageNames(age) labelNames(label) [num2str(length(sigElecIdx_sharpGen3{label,model,age})) ' Suprathreshold Electrodes'],...
                ['BF range = ' num2str(round(minMax_sharpGen(label,1,model,age),2)) ' - ' num2str(round(minMax_sharpGen(label,2,model,age),2))]],'FontSize',13)
            c = c+1;
        end
    end
end
% plot avg label
c = 1;
for age = 1:3
    for model = 1:2
        figure(555555)
        subplot(3,2,c)
        topoplotMB(zeros(129,1),locsEEGLAB129HCL,'style','straight','electrodes','on', 'maplimits',[0 7],...
            'emarker2',{sigElecIdx_sharpGen3_avgLab{model,age},'o',colorVec_avgLab{model,age},20,1}, ...
            'hcolor','k','emsize',15,'ecolor','k','colormap',white);
        title([modelNames(model) ageNames(age) [num2str(length(sigElecIdx_sharpGen3_avgLab{model,age})) ' Suprathreshold Electrodes'],...
            ['BF range = ' num2str(round(minMax_sharpGen_avgLab(1,model,age),2)) ' - ' num2str(round(minMax_sharpGen_avgLab(2,model,age),2))]],'FontSize',18)
        c = c+1;
    end
end
%% -----------------------------------------------------------------------
% MODEL COMPARISONS: tested model vs tested model, Figures 6-9
% ------------------------------------------------------------------------
modelComp_avgLab(:,1,:) = BF10_models_avgLab(:,1,:)./BF10_models_avgLab(:,2,:); % => elec x modelComp x age
modelComp_avgLab(:,2,:) = BF10_models_avgLab(:,2,:)./BF10_models_avgLab(:,1,:); % => elec x modelComp x age

modelComp_byLabel(:,:,1,:) = BF10_models_byLabel(:,:,1,:)./BF10_models_byLabel(:,:,2,:); % => elec x label x modelComp x age
modelComp_byLabel(:,:,2,:) = BF10_models_byLabel(:,:,2,:)./BF10_models_byLabel(:,:,1,:); % => elec x label x modelComp x age

for age = 1:3
    for modelComp = 1:2
        for label = 1:3
            data1 = modelComp_byLabel;
            tempElecs = find(data1(:,label,modelComp,age) >= 3);
            sigElecIdx_modelComp3{label,modelComp,age} = tempElecs;
            minMax_modelComp(label,1,modelComp,age) = min(data1(tempElecs,label,modelComp,age));
            minMax_modelComp(label,2,modelComp,age) = max(data1(tempElecs,label,modelComp,age));
        end % label
        data2 = modelComp_avgLab;
        tempElecs2 = find(data2(:,modelComp,age) >= 3);
        sigElecIdx_modelComp3_avgLab{modelComp,age} = tempElecs2;
        minMax_modelComp_avgLab(1,modelComp,age) = min(data2(tempElecs2,modelComp,age));
        minMax_modelComp_avgLab(2,modelComp,age) = max(data2(tempElecs2,modelComp,age));
    end % model
end % age

for age = 1:3
    for modelComp = 1:2
    % get colors right
        if modelComp == 1
            colors4topo = turbo(12);colors4topo = colors4topo([8 9 11],:);
        elseif modelComp == 2
            colors4topo = turbo(24);colors4topo = colors4topo([10 5 1],:);
        end
    % assign suprathres elecs their color
        dataBF = modelComp_byLabel(:,:,modelComp,age);
        for label = 1:3
            colorVec_tmp =[];
            BFvalues = dataBF(sigElecIdx_modelComp3{label,modelComp,age},label);
            for x = 1:size(BFvalues,1)
                if BFvalues(x) >=3 && BFvalues(x) < 10, colorVec_tmp(x,:) = colors4topo(1,:);
                elseif BFvalues(x) >=10 && BFvalues(x) < 20, colorVec_tmp(x,:) = colors4topo(2,:);
                elseif BFvalues(x) >=20, colorVec_tmp(x,:) = colors4topo(3,:);
                end
            end % suprathres elecs
            colorVec{label,modelComp,age} = colorVec_tmp;
        end % label
        dataBF = modelComp_avgLab(:,modelComp,age);
        colorVec_tmp = [];
        BFvalues = dataBF(sigElecIdx_modelComp3_avgLab{modelComp,age});
        for x = 1:size(BFvalues,1)
            if BFvalues(x) >=3 && BFvalues(x) < 10, colorVec_tmp(x,:) = colors4topo(1,:);
            elseif BFvalues(x) >=10 && BFvalues(x) < 20, colorVec_tmp(x,:) = colors4topo(2,:);
            elseif BFvalues(x) >=20, colorVec_tmp(x,:) = colors4topo(3,:);
            end
        end % suprathres elecs
        colorVec_avgLab{modelComp,age} = colorVec_tmp;
    end % model
end % age
% plot avg label
mCompNames = {'Sharpening > Generalization'; 'Generalization > Sharpening'};
c = 1;
for age = 1:3
    for modelComp = 1:2
        figure(987)
        subplot(3,2,c)
        topoplotMB(zeros(129,1),locsEEGLAB129HCL,'style','straight','electrodes','on', 'maplimits',[0 7],...
            'emarker2',{sigElecIdx_modelComp3_avgLab{modelComp,age},'o',colorVec_avgLab{modelComp,age},10,1}, ...
            'hcolor','k','emsize',10,'ecolor','k','colormap',white); 
        title([mCompNames(modelComp) ageNames(age)])% [num2str(length(sigElecIdx_modelComp3_avgLab{modelComp,age})) ' Suprathreshold Electrodes'],...
%             ['BF range = ' num2str(round(minMax_modelComp_avgLab(1,modelComp,age),2)) ' - ' num2str(round(minMax_modelComp_avgLab(2,modelComp,age),2))]],'FontSize',15)
        c = c+1;
    end
end
% plot per labe1
% c = 1;
for modelComp = 1:2
    c = 1;
    for age = 1:3
        for label = 1:3
            figure(modelComp)
            subplot(3,3,c)
            topoplotMB(zeros(129,1),locsEEGLAB129HCL,'style','straight','electrodes','on', 'maplimits',[0 7],...
                'emarker2',{sigElecIdx_modelComp3{label,modelComp,age},'o',colorVec{label,modelComp,age},10,1}, ...
                'hcolor','k','emsize',10,'ecolor','k','colormap',white); 
            title([mCompNames(modelComp) ageNames(age) lNames(label)])% [num2str(length(sigElecIdx_modelComp3{label,modelComp,age})) ' Suprathreshold Electrodes'],...
%                 ['BF range = ' num2str(round(minMax_modelComp(label,1,modelComp,age),2)) ' - ' num2str(round(minMax_modelComp(label,2,modelComp,age),2))]],'FontSize',13)
            c = c+1;
        end
    end
end
%% -----------------------------------------------------------------------
% Regionalized Topography, Figure 4
% ------------------------------------------------------------------------
Colors = [211 39 55; ...
        0 33 165; ...
        242 169 0; ...
        250 70 22;...
        34 136 72]./255;
for region = 1:5
    tmpElec = regions5{region};
    for x = 1:length(tmpElec)
        color4regions(tmpElec(x),:) = Colors(region,:);
    end
end
figure(1)
topoplotMB(zeros(129,1),locsEEGLAB129HCL,'style','straight','electrodes','on', 'maplimits',[0 1],...
            'emarker2',{([1:124])','o',color4regions,30,1}, ...
            'hcolor','k','emsize',15,'ecolor','k','colormap',white);
%% -----------------------------------------------------------------------
% Region-specific electrode counts/percentages and max BF value
% ------------------------------------------------------------------------

% get counts and max per region
for age = 1:3
    for model = 1:2
        for region = 1:5
        % vs null
            supraCountByRegion(region,age,model) = sum(BF10_models_avgLab(regions5{region},model,age)>=3);
            maxByRegion(region,age,model) = max(BF10_models_avgLab(regions5{region},model,age));
        % model comparisons
            supraCountByRegion_modelComp(region,age,model) = sum(modelComp_avgLab(regions5{region},model,age)>=3);
            maxByRegion_modelComp(region,age,model) = max(modelComp_avgLab(regions5{region},model,age));
            for label = 1:3
            % vs null
              supraCountByRegion_byLab(label,region,age,model) = sum(BF10_models_byLabel(regions5{region},label,model,age)>=3);
              maxByRegion_byLab(label,region,age,model) = max(BF10_models_byLabel(regions5{region},label,model,age));
            % model comparisons
              supraCountByRegion_byLab_modelComp(label,region,age,model) = sum(modelComp_byLabel(regions5{region},label,model,age)>=3);
              maxByRegion_byLab_modelComp(label,region,age,model) = max(modelComp_byLabel(regions5{region},label,model,age));
            end % label
        end % region
    end % model
end % age
% make tables
c = 1;
clear sting4table_vsNull_avgLab sting4table_modelComp_avgLab
for age = 1:3
   for model = 1:2
       for region = 1:5
            sting4table_vsNull_avgLab{c,region} = [num2str(round((supraCountByRegion(region,age,model)/size(regions5{region},2))*100,2)), ...
                '% (' num2str(round(maxByRegion(region,age,model),2)) ')'];
            sting4table_modelComp_avgLab{c,region} = [num2str(round((supraCountByRegion_modelComp(region,age,model)/size(regions5{region},2))*100,2)), ...
                '% (' num2str(round(maxByRegion_modelComp(region,age,model),2)) ')'];
       end % region
       c = c+1; 
   end % model
end % age
c = 1;
clear sting4table_vsNull_byLabel sting4table_modelComp_byLabel
for label = 1:3
   for age = 1:3
        for model = 1:2
           for region = 1:5
                sting4table_vsNull_byLabel{c,region} = [num2str(round((supraCountByRegion_byLab(label,region,age,model)/size(regions5{region},2))*100,2)), ...
                    '% (' num2str(round(maxByRegion_byLab(label,region,age,model),2)) ')'];
                sting4table_modelComp_byLabel{c,region} = [num2str(round((supraCountByRegion_byLab_modelComp(label,region,age,model)/size(regions5{region},2))*100,2)), ...
                    '% (' num2str(round(maxByRegion_byLab_modelComp(label,region,age,model),2)) ')'];
           end % region
           c = c+1; 
        end % model
   end % age
end % label
%% ------------------------------------------------------------------------
% model comparisons simple summary, Figure 10
% ------------------------------------------------------------------------

% w/ 3 regions: frontal/central, temporal/parietal, occipital
regionNames = {'Frontal/Central','Temporal/Parietal','Occipital'};
ageNames = {'6','9','12'};
for age = 1:3
    for model = 1:2
        for region = 1:3
          supraCountByRegion_modelComp2(region,age,model) = sum(modelComp_avgLab(regions3{region},model,age)>=3);
          maxByRegion_modelComp2(region,age,model) = max(modelComp_avgLab(regions3{region},model,age));
          for label = 1:3
              supraCountByRegion_byLab_modelComp2(label,region,age,model) = sum(modelComp_byLab(regions3{region},label,model,age)>=3);
              maxByRegion_byLab_modelComp2(label,region,age,model) = max(modelComp_byLab(regions3{region},label,model,age));
          end
        end
    end
end
perThres = 6; % how similar is too similar? ~1 elec in the smallest cluster (occipital)
for age = 1:3
    for region = 1:3
        s_tmp = (supraCountByRegion_modelComp2(region,age,1)./size(regions3{region},2))*100; % % sharp
        g_tmp = (supraCountByRegion_modelComp2(region,age,2)./size(regions3{region},2))*100;% % gen
    % if region needs a tie breaker
        if abs(diff([s_tmp,g_tmp])) <= perThres 
            fprintf([ageNames{age} 'mo, ' regionNames{region} ' is toss up\n'])
            if abs(diff(maxByRegion_modelComp2(region,age,:))) <= 5, output_yval2(region,age) = 2; % equivalent support
            elseif diff(maxByRegion_modelComp2(region,age,:)) > 5, output_yval2(region,age) = 1; % gen winner
            elseif diff(maxByRegion_modelComp2(region,age,:)) < -5, output_yval2(region,age) = 3; % sharp winner
            end
        else
            if diff([s_tmp,g_tmp]) > perThres, output_yval2(region,age) = 1; % gen winner
            elseif diff([s_tmp,g_tmp]) < -1*perThres, output_yval2(region,age) = 3;% sharp winner
            end
        end % if tiebreak needed

        for label = 1:3
            s_tmp2 = (supraCountByRegion_byLab_modelComp2(label,region,age,1)./size(regions{region},2))*100;
            g_tmp2 = (supraCountByRegion_byLab_modelComp2(label,region,age,2)./size(regions{region},2))*100;
            if abs(diff([s_tmp2,g_tmp2])) <= perThres
                fprintf([ageNames{age} 'mo, ' regionNames{region} ', ' lNames{label} ' is toss up\n'])
                if abs(diff(maxByRegion_byLab_modelComp2(label,region,age,:))) <= 5, output_yval2(region,age,label+1) = 2; % equivalent support
                elseif diff(maxByRegion_byLab_modelComp2(label,region,age,:)) > 5, output_yval2(region,age,label+1) = 1; % gen winner
                elseif diff(maxByRegion_byLab_modelComp2(label,region,age,:)) < -5, output_yval2(region,age,label+1) = 3; % sharp winner
                end
            elseif abs(diff([s_tmp2,g_tmp2])) > perThres 
                if diff([s_tmp2,g_tmp2]) > 5, output_yval2(region,age,label+1) = 1; % gen winner
                elseif diff([s_tmp2,g_tmp2]) < -1*perThres, output_yval2(region,age,label+1) = 3; % sharp winner
                end
            end % if tiebreak needed
        end % label

    end % region
end % age
% make figure 10
    tmp = [1,4,3,2];
    jitter = repmat([1;2;3],1,4);
    steps = [-.47 -.17 .17 .47];
    for age = 1:3
        for datatype = 1:4
            jitter4plot(age,datatype) = jitter(age,datatype) + steps(tmp(datatype));
        end
    end
    multiplier = 12;
    turbo12 = turbo(4*multiplier);
    turbo4plot = turbo12(multiplier/2:multiplier:end,:);
    turbo4plot(2,:) = [18 185 106]./255;
    regionNames = {'Frontal/Central','Temporal/Parietal','Occipital'};
    ageNames = {'6','9','12'};
    markers = {'o','s','d','^'};
    for region = 1:3
        figure(12299)
        subplot(1,3,region)
        for datatype = [1 4 3 2]
            plot(squeeze(output_yval2(region,:,datatype))+steps(tmp(datatype)),'Marker', markers{datatype},'MarkerSize',13,...
                'MarkerFaceColor',turbo4plot(datatype,:),'LineWidth',3,'Color',turbo4plot(datatype,:)), box off, hold on
        end
        hold off
        xlim([.75 3.25]), ylim([.25 3.75])
        yticks([1 3]), xticks(1:3), xticklabels({'6','9','12'}), xlabel('Age Group (months)')
        if region == 1
            yticklabels({'Generalization','Sharpening'}), ylabel('Favored Model')
        else yticklabels('')
        end
        title(regionNames(region))
    end
    dtNames = {'Collapsed Across Label','Individual','Category','Unlabeled'};
    legend(dtNames{[1 4 3 2]},'box','off','location','southoutside','orientation','horizontal','NumColumns',2)
%% ta da !!! :)