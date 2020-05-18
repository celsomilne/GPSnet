%Raw_Module_Generation generates the raw cancer module.
%
%   [module_genes_ids, module_gene_scores] =
%   Raw_Module_Generation(Cancer_Type) calculates and generates the raw
%   cancer module for a cancer type given by the string Cancer_Type.
%   Cancer_Type can be of the form 'name' or 'name.mat'. Alpha value used
%   is 0.5. The output module_gene_ids is a list of gene ids in the raw
%   disease module. module_gene_scores is the list of network smoothed
%   scores for each gene in module_gene_ids.
%
%   [module_genes_ids, module_gene_scores] =
%   Raw_Module_Generation(Cancer_Type, alpha) The alpha parameter is used
%   to determine the probability of the random walker moving to a random
%   neighbour during the Random Walk with Restart process (RWR).
function module_genes = Raw_Module_Generation(Cancer_Type,alpha)

    if nargin < 2
        alpha=0.5;
    end

    % Adjust the filename
    Cancer_Type = format_cancer_type(Cancer_Type);

    % Load in the PPI network (n X 2 array of cancer type-specific
    % protein-protein interactions based on tumour's RNA sequencing data)
    Net = load(strcat('Data_mat/Cancer_Specific_PPI/', Cancer_Type));
    
    % Load in mutation information (m X 2 array. Column 1 contains the gene
    % number, column 2 contains the mutation frequency.)
    Mutation = load(strcat('Data_mat/Mutation/',Cancer_Type));
    
    % Load in the lengths of genes (o X 2 array of <description>)
    Gene_Length = load('Data_mat/Gene_Length');
    
    % Extract from structs (due to how Matlab loads information)
    Net = Net.Net;
    Mutation = Mutation.Mutation;
    Gene_Length = Gene_Length.Gene_Length;

    % Eliminate highly connected genes (7316 and 7273) by deleting any rows
    % in which either of the two genes appears
    mask = any(ismember(Net,[7316,7273]), 2);
    Net(mask, :)=[];
    clear mask;
    
    % Select the largest component (see largest_component function).
    Net = largest_component(Net);
    connected_genes = unique(Net);

    % Match genes, their lengths and their mutation frequencies. Produce an
    % m X 3 matrix, where column 1 is the gene number, column 2 is the gene
    % length, and column 3 is gene mutation frequency.
    Mutation=sortrows(Mutation,1);
    [a,b]=ismember(Mutation(:,1),Gene_Length(:,1));
    Gene_Length(b(a),3)=Mutation(a,2);
    Mutation=Gene_Length;
    
    % Now match gene length and mutation frequency to the genes that we
    % observe in our largest connected component. We now have an n X 3
    % matrix with columns: [gene number, gene length, gene mutation
    % frequency].
    [a,b]=ismember(Mutation(:,1),connected_genes);
    connected_genes(b(a),2:3)=Mutation(a,2:3);
    
    mutation_frequency = connected_genes(:, 3);
    gene_lengths = connected_genes(:, 2);
    gene_ids = connected_genes(:, 1);
    
    % Calculate s(i) (see equation (1)) 
    s0 = mutation_frequency ./ gene_lengths;
    s0(isnan(s0)) = 0;

    % Create an adjacency matrix of the undirected PPI network, `Net`,
    % where genes that are connected in `Net` are given weights of 1 in
    % both directions.
    [~, G] = ismember(Net, gene_ids);
    G = sparse([G(:,1);G(:,2)],[G(:,2);G(:,1)], 1);
    
    % Run network smoothing (see supplementary note 3).
    F = network_smoothing(s0, G, alpha);

    % Generate the raw modules
    LL=200;
    modules = module_forming(G, F, LL);
    
    % Run for old behaviour
%     PM = [gene_ids, F];
%     [~, ct, ~] = fileparts(Cancer_Type);
%     Module_Forming_Process(Net, PM, ct, alpha, LL);
    
    % Calculate the scores for each module, then sort modules in descending
    % order of their scores
    mu = mean(F);
    score_fun = @ (M) module_score(F(M), mu);
    scores = cellfun(score_fun, modules);
    [~, sorted_idxs] = sort(scores, 'descend');
    modules_sorted = modules(sorted_idxs);
    
    % Find genes that appear in the top 1% of modules and calculate their
    % confidence scores
    top_1pc_idx = ceil(LL * 0.01);
    top_idxs = 1:top_1pc_idx;
    top_modules = modules_sorted(top_idxs);
    [confidence_scores, gene_idxs] = groupcounts(cell2mat(top_modules));
    
    % Now sort them and return the top L (300)
    [~, idxs] = sort(confidence_scores, 'descend');
    gene_idxs = gene_idxs(idxs);
    L = 300;
    L = min(L, length(gene_idxs));
    gene_idxs = gene_idxs(1:L);
    module_genes = {gene_ids(gene_idxs), F(gene_idxs)};
end

%largest_component calculates the largest component of the PPI network
%according to [1] and supplementary note 2.
%
%   LG = largest_component(G) returns the largest connected component in G,
%   the PPI network. In essence, the algorithm is a BFS of the network,
%   each level selecting proteins further and further from the seed to
%   create an image of the connected components in the PPI. The largest
%   connected component is returned.
%
%   References 
%   ---------- 
%   [1] Menche, J., Sharma, A., Kitsak, M., Ghiassian, S.D., Vidal, M.,
%   Loscalzo, J. and Barabási, A.L., 2015. Uncovering disease-disease
%   relationships through the incomplete interactome. Science, 347(6224),
%   p.1257601.
function LG = largest_component(G)
 
    % Return an empty array for an empty input PPI network. 
    if isempty(G)
        LG=[];
        return
    end

    max_component_size = -1;
    while ~isempty(G)
        
        % Randomly select a protein, and find any proteins that interact
        % with it (i.e. select any nodes in which the seed protein
        % appears).
        rand_idx = ceil(rand * size(G, 1));
        seed = G(rand_idx, 1);
        mask = any(ismember(G, seed), 2);
        
        % Select all rows (nodes) in the network in which the seed protein
        % appears. Delete them from G and save them to a new variable,
        % representing the current level of the BFS.
        [G, component] = mask_and_delete(G, mask);
        
        % Create a single list of the unique proteins found in the
        % current level
        unique_proteins = unique(component);
        
        while ~isempty(unique_proteins)
            
            % Mask of any rows which contain any of the unique proteins in
            % the current component. Fetch those rows from G, which is the
            % next level of the BFS and delete them from G.
            mask = any(ismember(G, unique_proteins), 2);
            [G, level] = mask_and_delete(G, mask);
            
            % Update the connected component to include the new level of
            % proteins and get the unique proteins in the next level of
            % BFS.
            component = [component; level]; %#ok
            unique_proteins=unique(level);
        end
        
        % If this component is larger than the previous largest component,
        % then set the largest component and the maximum component size
        component_size = length(unique(component));
        if component_size > max_component_size
            max_component_size = component_size;
            LG = component;
        end  
    end
end

%mask_and_delete fetches the rows defined by a boolean mask and removes
%them from the original matrix.
%
%   [modified, masked] = mask_and_delete(mat, mask) selects the rows of mat
%   defined in the mask, copies them to `masked` and returns a copy of
%   `mat` with those rows removed.
function [modified, masked] = mask_and_delete(mat, mask)
    masked = mat(mask, :);
    modified = mat;
    modified(mask, :) = [];
end

%network_smoothing run network smoothing on the sparse PPI graph given by
%an adjacency graph.
%
%   F = network_smoothing(s0, G, alpha) using the PPI graph given by the
%   sparse adjacency graph, G, with an initial score estimate given by s0
%   (which should be calculated as s0 = m(i) / l(i)), run diffusion network
%   smoothing with a given alpha value.
function F = network_smoothing(s0, G, alpha)

    % Calculate the inverse of the degree of gene j, 1/k(j)
    kj = sum(G, 2);
    inv_kj = sparse(diag(1./kj));
    W = G * inv_kj;
    
    % Convergence stopping criterion
    tol = 1e-8;
    
    % Run diffusion network smoothing
    % TODO : functions for conduction network smoothing and other case
    F = diffusion_network_smoothing(s0, W, alpha, tol);
end

%diffusion_network_smoothing calculate network smoothing scores using
%diffusion.
function stp1 = diffusion_network_smoothing(s0, W, alpha, tol)
    
    % Aribtrary initial value
    delta = 1;
    
    % Start at s_{t+1} = s_0
    stp1 = s0;
    
    % Run equation (4) until convergence
    while delta > tol
       st = stp1;
       stp1 = (1 - alpha) * W * st + alpha * s0;
       delta = sum((stp1 - st).^2);
    end
end

%module_forming forms the module from the PPI network and network smoothing
%scores.
%
%   module_forming(Net, s, num_modules) forms the raw cancer modules for
%   the PPI network given by the the sparse matrix Net, and the n X 2 array
%   s, where column 1 of s is the gene ids and column 2 is the network
%   smoothing score for that gene. The parameter num_modules determines how
%   many modules are generated.
%
%   See also hygepdf
%
%   References
%   ----------
%   [1] Bindea, G. et al. ClueGO: a Cytoscape plug-in to decipher
%   functionally grouped gene ontology and pathway annotation networks.
%   Bioinformatics 25, 1091–1093 (2009).
function modules = module_forming(Net, gene_scores, num_modules)

    N = length(gene_scores);
    modules = cell(1, num_modules);
    degree = sum(Net, 2);
    
    for i = 1:num_modules
       
        % Initially, a random gene is selected as the seed module
        seed_idx = ceil(rand * N);
        
        % Start the module
        M = seed_idx;
        while true
            
            % Get the genes that interact with the module (get all rows of
            % Net which are in the module Gamma, then find all the columns
            % where a node exists).
            module_genes = Net(M, :);
            [~, gamma_idxs] = find(module_genes > 0);
            gamma_idxs = unique(reshape(gamma_idxs, 1, []));

            % Parameters for calculating the connectivity significance
            k = degree(gamma_idxs);

            % Calculate the connectivity significance using ClueGo [1] in
            % equation (2).
            P = connectivity_significance(Net, gamma_idxs, M, k);
            
            % Calculate the expanded module score if a gene i is added to
            % the module (equation (3)).
            mu = mean(gene_scores);
            [Zm, Zmp1] = expanded_module_score(M, gamma_idxs, gene_scores, mu);
            
            % Select genes to add to M that are not already in M
            mask = P < 0.05 & Zmp1 > Zm & ~(ismember(gamma_idxs, M));
            
            % If we have no more genes to add, then break this loop
            if sum(mask) == 0
                break
            end
            
            % Add the new genes to the module
            add_to_M = gamma_idxs(mask);
            M = unique([M add_to_M]);
        end
        
        % Add the module to the return array
        modules{i} = M;
    end

end

%connectivity_significance calculates the connectivity significance as per
%equation (2). 
%
%   P = connectivity_significance(Net, gamma_idxs, M, k, N) returns a 1 X n
%   array of the connectivity significances for the n genes in Gamma
%   indexed by gamma_idxs. M is the current module, k is the degree of each
%   gene in Gamma. Net is a sparse N X N adjacnecy matrix of the largest
%   connected component of the cancer-specific PPI.
%
%   See also expanded_module_score
function P = connectivity_significance(Net, gamma_idxs, M, k)
    num_genes = length(k);
    m = length(M);
    N = size(Net, 1);
    
    % km is the number of gene i's neighbours that belong to the module
    gamma_module_neighbours = Net(gamma_idxs, M);
    km = sum(gamma_module_neighbours, 2);
    
    % Calculate P(i) for each gene
    P = zeros([1, num_genes]);
    for i = 1:num_genes
        
        % Parameters for matlab hydepdf (for more, type help hygepdf)
        ki = k(i);
        kmi = km(i);
        P(i) = sum(hygepdf(kmi:ki, N, m, ki));
    end
end
  
%expanded_module_score calculates the expanded module score as per equation
%(3).
%
%   [Zm, Zmp1] = expanded_module_score(M, gamma_idxs, s, mu) calculates the
%   module score of M (Zm) and the expanded module score (Zmp1) for a
%   module indexed by M, a Gamma set indexed by gamma_idxs, network
%   smoothed scores given by n and the mean network smoothed scores given
%   by mu. The last parameter may be omitted and will be calculated from s.
%
%   See also module_score
function [Zm, Zmp1] = expanded_module_score(M, gamma_idxs, s, mu)
    
    if nargin < 4
        mu = mean(s);
    end
    
    % For each gene in gamma, calculate the expanded module score
    % Calculate the module score
    Zm = module_score(s(M), mu);
    m = length(M);
    
    % Calculate the expanded module score if gene i is added to the module
    % according to equation (3).
    extended = s(gamma_idxs) - mu;
    extended = reshape(extended, 1, []);
    Zmp1 = ((Zm * sqrt(m)) + extended) / sqrt(m+1);
end


%module_score calculates the score Z_M for a module given by M. 
%
%   Zm = module_score(M, gene_scores) Calculates Zm as per equation (1) for
%   a module given by the n X 1 array of gene scores for a module, M, and
%   the mean, mu, of all gene scores.
%
%   Zm = module_score(M, gene_scores, normalize) Calculates the Zm score as
%   above. If normalize is true, the result is the same as equation (1).
%   Otherwise, it will return Zm * sqrt(m).
function Zm = module_score(M, mu, normalize)
    if nargin < 3
        normalize = true;
    end
    
    si = M;
    m = length(M);
    
    Zm = (sum(si) - mu);
    if normalize
        Zm = Zm / sqrt(m);
    end
end

%format_cancer_type formats the Cancer_Type string into a valid filename.
function formatted = format_cancer_type(Cancer_Type)
    [~, ~, ext] = fileparts(Cancer_Type);
    formatted = Cancer_Type;
    if isempty(ext)
        formatted = strcat(Cancer_Type, ".mat");
    end
end