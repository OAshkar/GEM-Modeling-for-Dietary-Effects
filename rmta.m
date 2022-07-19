function rmta= rmta(filename_source)
% context name is a string with Genotype_diet_tissue formula. 
% make sure to run this ahead of cobra and don't add anything to path!
setenv('ILOG_CPLEX_PATH','/home/omar/cplex'); % I am sure this is working without the last slash
initCobraToolbox(false)
changeCobraSolver('ibm_cplex', 'all')



cd ~/Documents/Projects/ObAnalysis/
% load generic MODEL
filename_model = 'Mouse-GEM.mat';
load 'models/Mouse-GEM.mat'
model = convertOldStyleModel(mouseGEM);
model.b = model.b(:,1);
model.S = full(model.S);
% e = verifyModel(model)

disp(filename_source)

tissueStr = split('WT_HFD_Ep', '_');
tissueStr = tissueStr{end}
 % make directoy with context name if does not exit
if not(isfolder(strcat('rMTAResults/',filename_source)))
      mkdir(strcat('rMTAResults/',filename_source))
end
cd(strcat('rMTAResults/',filename_source))

%Context source models
% Modeling was performed using fastcore. I will be using the median models.

% Source either WT HFD/ObOb HFD/ObOb ND per tissue
contextmodels_loc = dir('../../models/fastc_median2');
[~,N,~] = cellfun(@fileparts,cellstr({contextmodels_loc.name}), 'uni', 0);
contexModelsNames = N(1,3:30);
contexModelsNames = strrep(contexModelsNames,'ob_ob','obob');
contexModelsNames = cellfun(@(x) strsplit(x, '_'), contexModelsNames, ...
    'UniformOutput', false);
contexModelsNames = vertcat(contexModelsNames{:}); % To remove nesting of cell array newA
contexModelsNames = cell2table(contexModelsNames, "VariableNames", ...
    {'model', 'genotype' , 'diet', 'tissue'});

% filter out NDWT context
% NOTE commented to use snakemake filter. I will just pass string text
% filename parameter, and split it here

% for columnname = contexModelsNames.Properties.VariableNames  
%    contexModelsNames.(columnname{1}) = string(contexModelsNames.(columnname{1}));
% end
% filter=(contexModelsNames.genotype~='WT' | contexModelsNames.diet ~= 'ND');
% contextF = contexModelsNames(filter,:);
% contextF.genotype = strrep(contextF.genotype, 'obob', 'ob_ob');


% Loop the contexts and run rMTA %NOTE 
% for i = 1:size(contextF,1)
    % filename_source = strjoin(table2array(contextF(i,1:4)),'_');
    % filename_source_model = char(strcat(filename_source, '.mat'));
    % model_source = readCbModel(fullfile('models/fastc_median',filename_source_model));
    % fprintf(filename_source_model);
    % fprintf('\tMetabolic model uploaded\n');
 


%%%%%%%%%%%%%%%%%%%%%%
%% STEP 1: load the information to generate the reference flux
% Sampling for refererence Sample (Target sample)
% NOTE The sampling is done by cobra by and needs some transformations, so
% I can get Vref directly. But I apply tutorial manipulations

    samples = readmatrix(strcat('../../samplesMedianRes2/model_WT_ND_', ...
        tissueStr, '.csv'));
    %samples(1,:) = [];
    % remove first column to match cobra
    samples(:,1) = [];
    samples = rot90(samples);
    whos samples % should be 2000 columns and number of reactions in the generic model
    sampleStats = calcSampleStats(samples);
    % % from reduced index to model.rxns index
    % % FIXME I need to not rely so much on matlab so let's see
    modelSampling = readCbModel(char(strcat('../../models/fastc_median2/model_WT_ND_', tissueStr, ".mat")));
    % 
    idx = zeros(size(samples,1),1); % number of reactions in the context source model
    % Sampling Method and options
    sampling_method = 'ACHR';	%{('CHRR'), 'ACHR'}
    sampling_options = struct();
    sampling_options.nWarmupPoints = 5000;      	% (default)
    sampling_options.nPointsReturned  = 2000;   	% (default)
    sampling_options.nStepsPerPoint  = 500;    	% (default = 200)
    %[modelSampling,samples] = sampleCbModel(model,'sampleFiles',sampling_method,sampling_options);
    
    % Next block try to set reference flux to negative for matched reactions
    for ii = 1:numel(idx)
        try
            idx(ii) = find(strcmp(model.rxns,modelSampling.rxns{ii}));
        catch
            idx(ii) = find(cellfun(@length,strfind(model.rxns,modelSampling.rxns{ii})));
            sampleStats.mean(ii)=-1*sampleStats.mean(ii) %Those reactions are reversed;
        end
    end
    % Recalculate the results of all missing metabolites
     rxnInactive = setdiff(1:length(model.rxns),idx); % inactive reactions
     fields = fieldnames(sampleStats);
     for i = 1:numel(fields)
         % store old summary
         aux = sampleStats.(fields{i});
         % set new summary size to zeros of the generic model
         sampleStats.(fields{i}) = zeros(size(model.rxns));
         %keep all zero except orginially in the model
         sampleStats.(fields{i})(idx) = aux;
         clear aux
     end
    
    % % resize the samples matrix
     aux = samples;
     samples = zeros(size(model.rxns,1),sampling_options.nPointsReturned);
     samples (idx,:) = aux;
     clear aux;
    %%
    % Finally we can set a reference flux
    Vref = sampleStats.mean; % = length(model.rxns)
    
    %% STEP 2: load the information of the differentially expressed genes and calculate upregulated and downregulated reactions
    % in detail in the rMTA article in detail [2]).
    
    % Differentially expressed genes
    % Neccesary variables: 'gene','logFC','pval'
    % 'Gene_ID' must be the same nomenclature as the metabolic model
    % FIXME the next line needs special care!
    % Study must be uploaded as DISEASE VS HEALTHY/CONTROL
  
    % get the correct file name
    % detail in rMTA-DGE.R
    % dietWT: WT_HFD
    % dietgenotype: ObOb HFD
    % GenotypeND: Ob/Ob ND
    if regexp(filename_source, 'WT_HFD')
        disp('WT HFD')
        filename_differentially_expressed_genes = char(strcat(tissueStr, '_dietWT.csv'));
    elseif regexp(filename_source, 'ob_ob_HFD')
        disp('ObOb HFD')
        filename_differentially_expressed_genes = char(strcat(tissueStr,'_dietgenotype.csv'));
    elseif regexp(filename_source, 'ob_ob_ND')
        disp('ObOb ND')
        filename_differentially_expressed_genes = char(strcat(tissueStr,'_GenotypeND.csv'));
    else 
        error('Couldnt find matching df');
    end
    
    logFC_requiered = 0; % change if necesary
    pval_requiered = 0.1; % change if necesary
    
    differ_genes = readtable(fullfile('../../rMTAResults/dge',filename_differentially_expressed_genes),...
        'ReadVariableNames',true);
    
    differ_genes.gene = cellstr(string(differ_genes.gene)); % requiered variable by code
    
    % Here we should obtain an array similar to rxnHML, in which we have the
    % information of whatever an expresion should increase, decrease or nothing
    % (+1)R_f    (-1)R_b     (0)unchanged
    % This vector is called rxnFBS (Forward, Backward, Unchanged)
    
    % transcript_separator = '.'
    rxnFBS = diffexprs2rxnFBS(model, differ_genes, Vref, ...
                              'logFC', logFC_requiered, 'pval', pval_requiered);
    % change in rxnFBS all those reactions that are not active
    %  it is not possible to predict the direction of the change
    rxnFBS(rxnInactive) = 0;
    fprintf('\tThere are %u reactions that are differentially expressed after curation\n',sum(rxnFBS~=0));
    %end


    %% STEP 3: run rMTA algorithm, implemented in a function available in COBRA toolbox
    % Both MTA and rMTA use the same parameters: epsilon and alpha:
    
    % Define alpha values to calculate rMTA
    alpha_values = [0.66];  % (default range of values) % better to have more values
    %  It has been included 0.66 as it is the original value used in the
    %  original paper ('Yizhak et al, 2013')
    num_alphas = length(alpha_values);
    % Calculate epsilon, different for each reaction and with a minimum required change of 1e-3 (%default)
    epsilon = calculateEPSILON(samples, rxnFBS);
    %% 
    % One we have defined the parameters required by rMTA/MTA, we can use the 
    % COBRA function to calculate the transformation score (TS score):
    
    %save("mybeforerMTA")
    %load("mybeforerMTA")
    % execute the code
    % FIXME temp_rMTA.mat
    try
        [TSscore, deletedGenes, Vres] = rMTA(model, rxnFBS, Vref,'printLevel', 1);
    catch
        error("rMTA complains, exiting this shell ...")
        quit
    end
    %[TSscore, deletedGenes, Vres] = MTA(model, rxnFBS, Vref);
    
    dt = table();
    dt.bTS = TSscore.bTS ;
    dt.mTS = TSscore.mTS;
    dt.rTS = TSscore.rTS;
    dt.wTS = TSscore.wTS;
    dt.deletedGenes = deletedGenes;
    % fileX = char(strcat('rMTAResults/', filename_source, '/rmtares.csv'));
    fileX = 'rmtares.csv';
    writetable(dt,fileX,'Delimiter',',','QuoteStrings',true);
    % fileY = char(strcat('rMTAResults/',filename_source ,'/Vres/','Vres.mat'));

    % mkdir('Vres')
    % fileY = ;
    save('Vres.mat', 'Vres', 'TSscore', 'deletedGenes', 'differ_genes' );
end

