addpath(genpath('~/NoApps/cobratoolbox/'))
initCobraToolbox(false)
addpath('/opt/ibm/ILOG/CPLEX_Studio1210/cplex/matlab/x86-64_linux/')

%%%
cd ~/Documents/Projects/ObAnalysis


load Allfastcoremodels.mat
size(getExchangeRxns(model_out{1}))

tpm_data = csvread("counts/TPMaverage.csv", 1, 1);
tmp = readtable("counts/TPMaverage.csv", ReadRowNames = true);
rownames = cellstr(tmp.Properties.RowNames);
colnames = cellstr(tmp.Properties.VariableNames);


%%%%% FBA

% Biomass objective
%model.rxns(model.c == 1)
%solveLP(model)

% ATP objective
%constructEquations(model, 'MAR03964')

% 1g fat: 9 cal. 1g carbs 4 cal. 1g proteins proteins 4 cal
% Chow Diet: 66.2% Carbs /22.7% Protein / 11.1% lipds
% 48.8g carbs /21g protien /3g fat
% HFD: 22.6% carbs /12.9% protein / 64.5% lipids
% 43g carbs /17g protein/40g fat

% HFD: 1:7
% ND 3:1

media = {
    'glucose',
    'O2',
    'H2O',
    % 'CO2',
    % % essential aa
    % 'histidine', %'MAR09038' % His
    % 'isoleucine', %'MAR09039' %ISo
    % 'leucine', %'MAR09040' %leu
    % 'lysine', %'MAR09041' %Lys
    % 'methionine', %'MAR09042' %Met
    % 'phenylalanine', %'MAR09043' %Phe
    % 'threonine', %'MAR09044' %Thr
    % 'tyrosine', %'MAR09064' %Tyr
    % 'valine' , %'MAR09046' %Val
    %        % fats
    % 'fatty acid pool',          %MAR09209
    % 'High Density Lipoprotein',
    % 'fatty acid-VLDL pool',
    % 'fatty acid-uptake pool',
    % 'fatty acid-chylomicron pool',
    % 'cholesterol'
    % 'cholesterol-ester pool',
    % 'CL pool',
    % 'TAG-VLDL pool',
    % 'PI pool',
    % 'phospholipids extracellular pool',
    'stearate',
    'linoleate'
        };
% Get table of objective at different ratios
carbs_g = -10:1:-1
fats_g = -1:-1:-10
ratio = carbs_g ./ (carbs_g + fats_g)

carbs_mol = carbs_g * 5.5; %mmol
carbs_mol = carbs_mol ./10;
fats_mol = fats_g * 3.5; %mmol
fats_mol = fats_mol ./ 10;
ratio_mol = carbs_mol ./ (carbs_mol + fats_mol)

for i = 1:numel(colnames)
    model = model_out{i}
    model = setParam(model, 'obj', 'MAR03964', 1);
    for ii = 1:10
        lbs = [
            carbs_mol(ii), %glucose
            -1000, %O2
            -1000,  %H20
            fats_mol(ii)/2,
            fats_mol(ii)/2
        ];
        x = solveLP(setExchangeBounds(model, media, lbs));
        obj_flux(ii) = x.f;
        obj_flux(ii)
        colname = string(colnames(i))
        T.(colname) = obj_flux'
    end
end

T = struct2table(T)

%% FVA
changeCobraSolver ('gurobi', 'all');
global CBTDIR
%CBTDIR = fileparts(which('/home/omar/NoApps/cobratoolbox/initCobraToolbox.m')); % get the directory of the COBRA Toolbox
modelFileName = 'Recon2.0model.mat';
modelDirectory = getDistributedModelFolder(modelFileName); %Look up the folder for the distributed Models.
modelFileName= [modelDirectory filesep modelFileName]; % Get the full path. Necessary to be sure, that the right model is loaded
model = readCbModel(modelFileName);

% get uptake reactions
[selExc, selUpt] = findExcRxns(model);
uptakes = model.rxns(selUpt);

generateMexFastFVA('/opt/ibm/ILOG/CPLEX_Studio1210/cplex')
addpath('/home/omar/NoApps/cobratoolbox/.tmp')
[minFluxF1, maxFluxF1, optsol, ret, fbasol, fvamin, fvamax, statussolmin, statussolmax] = fastFVA(model_out{1});

subuptakeModel = extractSubNetwork(model, uptakes);
hiCarbonRxns = findCarbonRxns(subuptakeModel,1);
modelalter = changeRxnBounds(model, hiCarbonRxns, 0, 'b');

energySources = {'EX_adp'; 'EX_amp(e)'; 'EX_atp(e)'; 'EX_co2(e)';...
                 'EX_coa(e)'; 'EX_fad(e)'; 'EX_fe2(e)'; 'EX_fe3(e)'; 'EX_gdp(e)';...
                 'EX_gmp(e)'; 'EX_gtp(e)'; 'EX_h(e)'; 'EX_h2o(e)'; 'EX_h2o2(e)';...
                 'EX_nad(e)'; 'EX_nadp(e)'; 'EX_no(e)'; 'EX_no2(e)'; 'EX_o2s(e)'};
modelalter = changeRxnBounds (modelalter, energySources, 0, 'l');

% modelfva1 represents aerobic condition
modelfva1 = modelalter;
modelfva1 = changeRxnBounds(modelfva1, 'EX_glc(e)', -20, 'l');
modelfva1 = changeRxnBounds(modelfva1, 'EX_o2(e)', -1000, 'l');
% modelfva2 represents anaerobic condition
modelfva2 = modelalter;
modelfva2 = changeRxnBounds(modelfva2, 'EX_glc(e)', -20, 'l');
modelfva2 = changeRxnBounds(modelfva2, 'EX_o2(e)',  0, 'l');

% [minFlux, maxFlux, Vmin, Vmax] = fluxVariability(model,...
% optPercentage,osenseStr, rxnNameList, verbFlag, allowLoops, method);
% Selecting several reactions of the model that we want to analyse with fva
rxnsList = {'DM_atp_c_'; 'ACOAHi'; 'ALCD21_D'; 'LALDO'; 'ME2m';...
               'AKGDm'; 'PGI'; 'PGM'; 'r0062'}



% % Run FVA analysis for the model with the constraints that simulates aerobic conditions:
% [minFlux1, maxFlux1, Vmin1, Vmax1] = fluxVariability(modelfva1, 100, 'max', rxnsList)
% % Run FVA analysis for the model with the constraints that
% % simulates anaerobic conditions:
% [minFlux2, maxFlux2, Vmin2, Vmax2] = fluxVariability(modelfva2, [], [], rxnsList)

changeCobraSolver ('ibm_cplex', 'all', 1);

[minFluxF1, maxFluxF1, optsol, ret, fbasol, fvamin, fvamax, ...
 statussolmin, statussolmax] = fastFVA(modelfva1);

[minFluxF2, maxFluxF2, optsol2, ret2, fbasol2, fvamin2, fvamax2,...
 statussolmin2, statussolmax2] = fastFVA(modelfva2);

ymaxf1 = maxFluxF1;
yminf1 = minFluxF1;
ymaxf2 = maxFluxF2;
yminf2 = minFluxF2;
maxf =table(ymaxf1, ymaxf2);
minf =table(yminf1, yminf2);
maxf = table2cell(maxf);
minf = table2cell(minf);
figureplot3 = bar(cell2mat(maxf(1:end, :)));
hold on
plot4 = bar(cell2mat(minf(1:end, :)));
hold off
xticks([0 2000 4000 6000 8000 10600])
yticks([-1000 -800 -600 -400 -200 0 200 400 600 800 1000])
xlabel('All reactions in the model')
ylabel('Fluxes')
legend({'Aerobic', 'Anaerobic'})
title('Variations in fluxes in the aerobic and anaerobic conditions')


%%% FVA serine synthesis and uptake %%%
% https://plos.figshare.com/articles/dataset/MATLAB_script_used_for_the_flux_variability_analysis_in_SDC_media_/12868976/1

% vary serine uptake rates
uptakeRates = 0:-0.1:-1;

% prepare vectors for min and max values of fluxes
fluxMin = zeros(size(uptakeRates));
fluxMax = zeros(size(uptakeRates));

% control variable
index = 1;

% perform FVA for each uptake rate
for uptakeRate = uptakeRates

    % set serine uptake rate
    tmpModel = changeRxnBounds(model,'r_1906',uptakeRate, 'l');

    % run FVA
    [tmpVMin, tmpVMax, tmpOptSol, tmpRet, tmpFbaSol, tmpFvaMin, tmpFvaMax]...
        = fastFVA(tmpModel, 99, 'max', 'ibm_cplex', [serRxns(2:3);serRxns(6)]');

    % extract min and max fluxes from reactions of SHM1, SHM2 and SER2
    tmpFvaMin = tmpFvaMin(findRxnIDs(tmpModel,[serRxns(2:3);serRxns(6)]'),:);
    tmpFvaMax = tmpFvaMax(findRxnIDs(tmpModel,[serRxns(2:3);serRxns(6)]'),:);

    % change flux direction so that positive values always mean serine
    % production
    coefs = [-1,-1,1]';

    tmpFvaMin = tmpFvaMin .* coefs;
    tmpFvaMax = tmpFvaMax .* coefs;

    % combine reactions to one net reaction
    fluxMin(index) = min(sum(tmpFvaMin));
    fluxMax(index) = max(sum(tmpFvaMax));

    index = index + 1;

end

% plot flux variability of SHM1, SHM2 and SER2 combined
figure
p = errorbar(-uptakeRates, mean([fluxMin;fluxMax]), fluxMin - mean([fluxMin;fluxMax]), fluxMax - mean([fluxMin;fluxMax]));
p.LineStyle = 'none';
p.LineWidth = 2;
set(p,'Color','black')
yline(0)
xlim([-0.1,1.1])
xticks(-uptakeRates)
ylabel('fluxes [mmol/(gDW*h)]')
xlabel('serine uptake rate [mmol/(gDW*h)]')




%load_ext pymatbridge
fastFVA()
