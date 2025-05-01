%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = compareDist

%Load model:
model = load('../../models/yeast801.mat');
model = model.model;

%Block alternative phosphofruktokinase:
model.ub(strcmp(model.rxns,'r_0887')) = 0; %ATP + sedohept-7P -> ADP + H+ + sedohept-1,7biP

%Make irreversible with RAVEN function:
model = ravenCobraWrapper(model);
model = convertToIrrev(model);

%Add heme reaction:
posH  = strcmp(model.mets,'s_3714');    %heme a [cytoplasm]
model = addReaction(model, 'rEx', ...
                    'reactionName', 'heme exchange', ...
                    'metaboliteList', model.mets(posH), ...
                    'stoichCoeffList', -1, ...
                    'lowerBound', 0, ...
                    'upperBound', 1000);

%Analysis will be around the experimental biomass yield:
Ysx = 0.122;    %gDW/g(gluc)
Ysx = Ysx*180;  %gDW/mol(gluc)
Ysx = Ysx/1000; %gDW/mmol(gluc)

%Simulation:
results.model   = model;
results.glucose = compare_substrate(model,'heme exchange','glucose',Ysx);

%Create gene table:
results.geneTable      = cell(length(results.glucose.genes),3);
results.geneTable(:,1) = results.glucose.genes;
results.geneTable(:,2) = results.glucose.geneNames;
results.geneTable(:,3) = num2cell(results.glucose.k_genes);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FC = compare_substrate(model,product,substrate,Ysx)

%Simulate WT (100% growth):
FC.flux_WT = simulateGrowth(model,product,substrate,1);

%Simulate forced (X% growth and the rest towards product) based on yield:
posX     = strcmp(model.rxnNames,'growth');
alphaExp = Ysx/FC.flux_WT(posX);
alpha    = (alphaExp/2):(alphaExp/10):(2*alphaExp);
FC.alpha = alpha;
v_matrix = zeros(length(model.rxns),length(alpha));
k_matrix = zeros(length(model.rxns),length(alpha));
for i = 1:length(alpha)
    FC.flux_MAX   = simulateGrowth(model,product,substrate,alpha(i));
    v_matrix(:,i) = FC.flux_MAX;
    k_matrix(:,i) = FC.flux_MAX./FC.flux_WT;
end

%Generate rxn equations:
rxnEqs = printRxnFormula(model,model.rxns,true,true,true);

%Take out rxns with no grRule:
withGR   = ~cellfun(@isempty,model.grRules);
v_matrix = v_matrix(withGR,:);
k_matrix = k_matrix(withGR,:);
gene_rxn = model.rxnGeneMat(withGR,:);
FC.rxns  = [model.rxns(withGR) model.rxnNames(withGR) model.grRules(withGR) rxnEqs(withGR)];

%Filter out rxns that are always zero -> k=0/0=NaN:
non_nan  = sum(~isnan(k_matrix),2) > 0;
v_matrix = v_matrix(non_nan,:);
k_matrix = k_matrix(non_nan,:);
gene_rxn = gene_rxn(non_nan,:);
FC.rxns  = FC.rxns(non_nan,:);

%Replace remaining NaNs with 1s:
k_matrix(isnan(k_matrix)) = 1;

%Replace any Inf value with 1000 (maximum value is ~700):
k_matrix(isinf(k_matrix)) = 1000;

%Filter out values that are inconsistent at different alphas:
always_down  = sum(k_matrix <= 1,2) == length(alpha);
always_up    = sum(k_matrix >= 1,2) == length(alpha);
incons_rxns  = always_down + always_up == 0;
incons_genes = sum(gene_rxn(incons_rxns,:),1) > 0;
incons_rxns  = sum(gene_rxn(:,incons_genes),2) > 0;
v_matrix     = v_matrix(~incons_rxns,:);
k_matrix     = k_matrix(~incons_rxns,:);
gene_rxn     = gene_rxn(~incons_rxns,:);
FC.rxns      = FC.rxns(~incons_rxns,:);

%Order from highest to lowest k:
FC.k_rxns   = mean(k_matrix,2);
[~,order]   = sort(FC.k_rxns,'descend');
FC.k_rxns   = FC.k_rxns(order,:);
FC.v_matrix = v_matrix(order,:);
FC.k_matrix = k_matrix(order,:);
gene_rxn    = gene_rxn(order,:);
FC.rxns     = FC.rxns(order,:);

%Create list of remaining genes and filter out any inconsistent score:
FC.genes     = model.genes(sum(gene_rxn,1) > 0);
FC.geneNames = model.geneShortNames(sum(gene_rxn,1) > 0);
FC.k_genes   = zeros(size(FC.genes));
gene_rxn     = gene_rxn(:,sum(gene_rxn,1) > 0);
cons_genes   = false(size(FC.genes));
for i = 1:length(FC.genes)
    k_set         = FC.k_rxns(gene_rxn(:,i) > 0);
    always_down   = sum(k_set <= 1) == length(k_set);
    always_up     = sum(k_set >= 1) == length(k_set);
    cons_genes(i) = always_down + always_up == 1;
    FC.k_genes(i) = mean(k_set);
end
FC.genes     = FC.genes(cons_genes);
FC.geneNames = FC.geneNames(cons_genes);
FC.k_genes   = FC.k_genes(cons_genes);

%Filter any value between mean(alpha) and 1:
unchanged    = (FC.k_genes >= mean(alpha) - 1e-3) + (FC.k_genes <= 1 + 1e-3) == 2;
FC.genes     = FC.genes(~unchanged);
FC.geneNames = FC.geneNames(~unchanged);
FC.k_genes   = FC.k_genes(~unchanged);

%Order from highest to lowest k:
[~,order]    = sort(FC.k_genes,'descend');
FC.genes     = FC.genes(order,:);
FC.geneNames = FC.geneNames(order,:);
FC.k_genes   = FC.k_genes(order,:);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%