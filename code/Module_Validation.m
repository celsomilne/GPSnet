%MODULE_VALIDATION ensures that a generated module is a valid representation
%of the cancer module.
%
%   [enrichment, p_values] = Module_Validation(module_genes, cancer_type)
%   calculates the enrichment modules and p values for a given set of
%   module genes, as per the Methods section of the paper, using a random
%   permutation test.
function [enrichment, p_values] = Module_Validation(module_genes, cancer_type)

    % Load in the driver gene, associated gene and survival gene from
    % storage
    Driver_Gene = load(strcat('Data_mat/Driver_Gene/',cancer_type));
    Associated_Gene = load(strcat('Data_mat/Associated_Gene/',cancer_type));
    Survival_Gene = load(strcat('Data_mat/Survival_Gene/',cancer_type));
    Net = load('Data_mat/Net_PPI');
    
    % Extract from structs
    Driver_Gene = Driver_Gene.Driver_Gene;
    Associated_Gene = Associated_Gene.Associated_Gene;
    Survival_Gene = Survival_Gene.Survival_Gene;
    Net = Net.Net;
    
    % Get all the unique genes that appear in the human PPI
    genes = unique(Net(:,1:2));
    LL=100;

    % Enrichment is a LL x 3 array containing the values needed to analyse
    % a) canonical disease pathway enrichment, b) well-known disease
    % associated gene enrichment analysis, and c) survival gene enrichment
    % analysis.
    enrichment = zeros(LL,3);
    
    % Find the size of the intersection between the cancer gene and each of
    % the driver, associated and survival genes.
    enrichment(1,1) = intersect_length(module_genes,Driver_Gene);
    enrichment(1,2) = intersect_length(module_genes,Associated_Gene);
    enrichment(1,3) = intersect_length(module_genes,Survival_Gene);
    
    % Randomly select genes from the human PPI to compare against the given
    % intersection (supplementary figure 16).
    for i=2:LL
        Random_Gene = genes(randperm(length(genes), length(module_genes)));
        enrichment(i,1) = intersect_length(Random_Gene,Driver_Gene);
        enrichment(i,2) = intersect_length(Random_Gene,Associated_Gene);
        enrichment(i,3) = intersect_length(Random_Gene,Survival_Gene);
    end

    % Calculte P Values for each of the three enrichment analyses
    p_values = ones([1, size(enrichment, 2)]);
    for i=1:size(enrichment, 2)
        [mu, sigma] = normfit(enrichment(2:end,i));
        p_values(i) = 1 - normcdf(enrichment(1,i), mu, sigma);
    end
    
end

%INTERSECT_LENGTH calculates the length of the intersect between two sets A
%and B.
function len = intersect_length(A, B)
    % TODO : might this be better as IoU?
    len = length(intersect(A, B));
end