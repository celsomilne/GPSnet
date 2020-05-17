%Raw_Module_Generation generates the raw cancer module.
%
%   Raw_Module_Generation(Cancer_Type) calculates and generates the raw
%   cancer module for a cancer type given by the string Cancer_Type.
%   Cancer_Type can be of the form 'name' or 'name.mat'. Alpha value used
%   is 0.5.
%
%   Raw_Module_Generation(Cancer_Type, alpha) The alpha parameter is used
%   to determine the probability of the random walker moving to a random
%   neighbour during the Random Walk with Restart process (RWR).
function Raw_Module_Generation(Cancer_Type,alpha)

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
    Net = sparse([G(:,1);G(:,2)],[G(:,2);G(:,1)], 1);
    
    % Run network smoothing (see supplementary note 3).
    F = network_smoothing(s0, Net, alpha);
    s = [gene_ids, F];

    % Generate the raw modules
    LL=10;
    modules = module_forming(Net, s, LL);
    
    % Calculate the scores for each module, then sort modules in descending
    % order of their scores
    scores = 
    
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
            component = [component; level];
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
function modules = module_forming(Net, s, num_modules)

    gene_ids = s(:, 1);
    gene_scores = s(:, 2);
    N = length(gene_ids);
    modules = cell(1, num_modules);
    degree = sum(Net, 2);
    
    parfor i = 1:num_modules
       
        % Initially, a random gene is selected as the seed module
        seed_idx = ceil(rand * N);
        
        % Start the module
        M = seed_idx;
        while true
            
            % Get the genes that interact with the module (get all rows of
            % Net which are in the module M, then find all the columns where
            % a node exists).
            module_genes = Net(M, :);
            [~, gamma_idxs] = find(module_genes > 0);
            gamma_idxs = reshape(gamma_idxs, 1, []);

            % Parameters for calculating the connectivity significance
            k = degree(gamma_idxs);
            m = length(gamma_idxs);

            % Calculate the connectivity significance using ClueGo [1] in
            % equation (2).
            P = connectivity_significance(m, k, N);

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

function P = connectivity_significance(m, k, N)
    num_genes = length(k);
    
    % Calculate P(i) for each gene
    P = zeros([1, num_genes]);
    for i = 1:num_genes
        
        % Parameters for matlab hydepdf (for more, type help hygepdf)
        X_ = k(i:m);
        M_ = N;
        K_ = m;
        N_ = k(i);
        P(i) = sum(hygepdf(X_, M_, K_, N_));
    end
end
  
function [Zm, Zmp1] = expanded_module_score(M, gamma_idxs, s, mu)
    
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
%%%%%%%% generate the raw module
function Module_Forming_Process(Net,PM,Cancer_Type,alpha,LL)
    
    tic;
    Str_alpha=num2str(100*alpha);

    gene_ids = PM(:, 1);
    network_scores = PM(:, 2);
    
    Average_S=mean(PM(:,2));
    Gene_List=unique(Net(:));
    N=length(Gene_List);
    Degree=sum(Net,2);
    Module{LL,1}=[];
    Score=zeros(LL,3);
    Seed_Gene=1:length(Gene_List);
    for i=1:LL
        seed=Seed_Gene(ceil(rand*length(Seed_Gene)));
        node=seed;
        score=sum(PM(node,2)-Average_S)/sqrt(1);
        for j=2:1000
            [a,b]=find(Net(:,node)==1);
            node_extend=setdiff(a,node);
            p_extend=[];
            score1=(sum(PM(node,2)-Average_S)+PM(node_extend,2)-Average_S)/sqrt(j);
            m=find(score1>1.01*score);
            if isempty(m)
                break
            else
                node_extend=node_extend(m);
                score1=score1(m);
                for k=1:length(node_extend)
                    ks_extend=length(intersect(find(Net(node_extend(k),:)==1),node));
                    k_extend=Degree(node_extend(k));
                    p_extend(k,1)= sum(hygepdf(ks_extend:k_extend,N,j-1,k_extend));
                end
                node_p_s=[node_extend p_extend score1];            
                node_p_s(find(node_p_s(:,2)>0.05),:)=[];
                if isempty(node_p_s)
                    break
                else    
                    node_p_s=sortrows(node_p_s,-3);
                    x=PM(node_p_s(:,1),2)/sum(PM(node_p_s(:,1),2));
                    x=[0;cumsum(x)];
                    xx=find(x<rand,1,'last');
                    node=[node;node_p_s(xx,1)];
                    score=node_p_s(xx,3);
                end
            end
        end
        Module{i,1}=Gene_List(node);
        Score(i,1:3)=[Gene_List(seed) score length(node)];
        if mod(i,1)==0
            toc;
            disp(i)
            save(['Data_mat/Raw_Module/Raw_Module_',Cancer_Type,'_',Str_alpha,'.mat'],'Module','Score')
        end
    end
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