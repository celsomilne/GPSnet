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
    si = mutation_frequency ./ gene_lengths;
    si(isnan(si)) = 0;
    s = [gene_ids, si];

    % Network smoothing
    [~, G]=ismember(Net, gene_ids); %%%%% mutation smoothing
    G=sparse([G(:,1);G(:,2)],[G(:,2);G(:,1)], 1);
    F=network_smoothing(s(:,2),G,alpha,1);
    s(:,2)=F;

    %%%% generate raw modules (the number of raw module is LL) 
    LL=60000;
    Module_Forming_Process(Net,s,Cancer_Type,alpha,LL);

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
%   Loscalzo, J. and Barab�si, A.L., 2015. Uncovering disease-disease
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


%%%%%% network smoothing using the network spreading process
function F=network_smoothing(Mutation,G,alpha,flag);
 
Mutation=Mutation/sum(Mutation)*length(Mutation);
F=Mutation;
Y=Mutation;
sn=sum(G,2);
D=sparse(diag(1./sn));
W1=G*D;   %%%% md
delta=100;
k=0;
x=F;
switch flag
    case 1  %%%%% diffusion
        while delta>10^(-8)
            F0=F;
            F=(1-alpha)*W1*F0+alpha*Y;
            delta=sum((F-F0).^2);
        end  
    case 2  %%%% conduction
        F=F';
        Y=Y';
        while delta>10^(-8)
            F0=F;
            F=(1-alpha)*F0*W1+alpha*Y;
            delta=sum((F-F0).^2);
        end
        F=F';
    case 3
        W1=G;
        [a,b]=find(W1~=0);
        for i=1:length(a)
            W1(a(i),b(i))=1/sqrt(sn(a(i))*sn(b(i)));
        end
        while delta>10^(-8)
            F0=F;
            F=(1-alpha)*W1*F0+alpha*Y;
            delta=sum((F-F0).^2);
        end 
end
F=F*sum(Mutation)/length(Mutation);
end


%%%%%%%% generate the raw module
function Module_Forming_Process(Net,PM,Cancer_Type,alpha,LL);
%%%%%% module merging, zscore=sum(pm-<pm>)/sqrt(k)
%%%%%% seclect probability pm*fi(1-p)/sum(pm*fi(1-p))
tic;
Str_alpha=num2str(100*alpha);

Average_S=mean(PM(:,2));
Gene_List=unique(Net(:));
N=length(Gene_List);
[a,b]=ismember(Net,Gene_List);
Net=sparse([b(:,1);b(:,2)],[b(:,2),b(:,1)],1);
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

%format_cancer_type formats the Cancer_Type string into a valid filename.
function formatted = format_cancer_type(Cancer_Type)
    [~, ~, ext] = fileparts(Cancer_Type);
    formatted = Cancer_Type;
    if isempty(ext)
        formatted = strcat(Cancer_Type, ".mat");
    end
end