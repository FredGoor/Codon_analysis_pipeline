% Codon/Amino Acid Clustering Analysis Script with Heatmap Mode Selection and Auto File Handling
clear all; clc;
rng default; % For reproducibility

%% ==== USER SETTINGS ==== %%
% Choose one of the following display modes for the heatmap:
% 1 = Amino acid usage (20 AAs)
% 2 = Absolute codon usage including STOP/Trp/Met
% 3 = Absolute codon usage excluding STOP/Trp/Met
% 4 = Relative codon usage excluding STOP/Trp/Met
analysis_mode = 4;

%% ==== 1. Ask user for codon usage table ==== %%
disp('---------------------- Open codon datatable ----------------------');
[file, folder] = uigetfile({'*.xlsx','Excel Files (*.xlsx)'}, ...
    'Select a codon usage table Excel file (e.g., LT2_cds_from_genomic_codons.xlsx)');

if isequal(file, 0)
    error('No file selected.');
end

filename = fullfile(folder, file);
[~, nameOnly] = fileparts(file);

% Derive file prefix automatically (remove "_codons" if present)
if contains(nameOnly, 'codons')
    strainPrefix = erase(nameOnly, 'codons');
elseif contains(nameOnly, 'Codons')
    strainPrefix = erase(nameOnly, 'Codons');
else
    % fallback: just use the name directly
    strainPrefix = nameOnly;
end

% Define related filenames
geneFileName = [strainPrefix 'geneIDs.xlsx'];
gene_IDs = fullfile(folder, geneFileName);
outputFile = fullfile(folder, [strainPrefix 'UMAP.xlsx']);

fprintf('Codon file selected: %s\n', file);
fprintf('Gene ID file assumed: %s\n', geneFileName);
fprintf('Output Excel will be: %s\n', outputFile);

%% ==== 2. Load codon usage table ==== %%
C_abs_table = readtable(filename, 'ReadRowNames', true, 'ReadVariableNames', true);

% Clean headers
vars = C_abs_table.Properties.VariableNames;
vars = strrep(vars, 'END', 'STOP'); % Fix STOP codon label
C_abs_table.Properties.VariableNames = vars;

% Extract variables
varNames = C_abs_table.Properties.VariableNames;
RowNames = C_abs_table.Properties.RowNames;
C_abs = C_abs_table{:,:};

%% ==== 3. Amino acid usage & relative codon usage ==== %%
aaNames = cellfun(@(x) strtok(x, '_'), varNames, 'UniformOutput', false);
uniqueAAs = unique(aaNames, 'stable');
nGenes = size(C_abs, 1);
nAAs = numel(uniqueAAs);

AA = zeros(nGenes, nAAs);
for i = 1:nAAs
    cols = strcmp(aaNames, uniqueAAs{i});
    AA(:, i) = sum(C_abs(:, cols), 2);
end
AA_table = array2table(AA, 'RowNames', RowNames, 'VariableNames', uniqueAAs);

[nGenes, nCodons] = size(C_abs);
C_rel = zeros(nGenes, nCodons);
for j = 1:nCodons
    idx = find(strcmp(uniqueAAs, aaNames{j}), 1);
    C_rel(:, j) = C_abs(:, j) ./ AA(:, idx);
end
C_rel(isnan(C_rel)) = 0;
C_rel_table = array2table(C_rel, 'RowNames', RowNames, 'VariableNames', varNames);

cols_to_remove = strcmp(aaNames, 'STOP') | strcmp(aaNames, 'Trp') | strcmp(aaNames, 'Met');

switch analysis_mode
    case 1
        Usage = AA; codon_labels = uniqueAAs;
    case 2
        Usage = C_abs; codon_labels = varNames;
    case 3
        Usage = C_abs(:, ~cols_to_remove); codon_labels = varNames(~cols_to_remove);
    case 4
        Usage = C_rel(:, ~cols_to_remove); codon_labels = varNames(~cols_to_remove);
    otherwise
        error('Invalid analysis_mode value. Choose 1, 2, 3, or 4.');
end

%% ==== 4. Normalize ==== %%
genome_means = mean(Usage);
deviations = Usage - genome_means;
SD = std(deviations);
values = deviations ./ SD;

%% ==== 5. UMAP Clustering ==== %%
Y = run_umap(values, 'parameter_names', codon_labels, ...
    'metric', 'euclidean', 'n_neighbors', 20, 'min_dist', 0.09, ...
    'n_components', 2, 'randomize', false);

figure;
scatter(Y(:,1), Y(:,2), 6, 'filled');
xlabel('UMAP 1'); ylabel('UMAP 2');
title('UMAP projection of genes by codon usage');
set(gca, 'FontName', 'Arial', 'FontSize', 8);
axis tight; box on;

%% ==== 6. Hierarchical sorting ==== %%
D1 = pdist(Y); tree1 = linkage(D1);
gene_order = optimalleaforder(tree1, D1);
ordered_genes = RowNames(gene_order);

D2 = pdist(values.', 'spearman'); tree2 = linkage(D2);
leaforder_cod = optimalleaforder(tree2, D2);
codons_reorder = codon_labels(leaforder_cod);

%% ==== 7. Heatmap ==== %%
values_reorder = values(gene_order, leaforder_cod).';
row_smooth = smoothdata(values_reorder, 2, 'gaussian', 5);
figure('Units','normalized','Position',[0 0 1 1]); % full screen
imagesc(row_smooth);
colormap('parula'); c = colorbar;
c.TickDirection = 'out'; c.Ticks = [-1 0 1 2 3 4]; caxis([-0.8 2.6]);
pbaspect([4 1 1]);
set(gca, 'XAxisLocation', 'top', 'TickDir', 'out', 'XMinorTick', 'on');
xticks(500:500:size(values,1));
xticklabels({'','1,000','','2,000','','3,000','','4,000'});
set(gca, 'TickLength', [0.01, 0.005]);
set(gca, 'YTick', 1:length(codons_reorder));
set(gca, 'YTickLabel', strrep(codons_reorder, '_', '\_'));
set(gca, 'FontName', 'Arial', 'FontSize', 6);
set(gca, 'LooseInset', get(gca,'TightInset'));
saveas(gcf, fullfile(folder, [strainPrefix 'UMAP_heatmap.jpeg']));

%% ==== 8. Export reordered genes ==== %%
gene_IDs_table = readtable(gene_IDs, 'ReadVariableNames', true);
if ~iscell(gene_IDs_table.LocusTag)
    gene_IDs_table.LocusTag = cellstr(gene_IDs_table.LocusTag);
end
[~, idx] = ismember(ordered_genes, gene_IDs_table.LocusTag);
if any(idx == 0)
    warning('Some gene IDs from the codon usage table were not found in the gene_IDs file.');
end
reordered_gene_IDs = gene_IDs_table(idx, :);
nGenes = height(reordered_gene_IDs);
UMAP_order = arrayfun(@(x) sprintf('UMAP_%04d', x), (1:nGenes)', 'UniformOutput', false);
newGeneTable = table(UMAP_order, reordered_gene_IDs.LocusTag, reordered_gene_IDs.GeneSymbol, ...
    reordered_gene_IDs.EntrezGeneID, reordered_gene_IDs.ProteinDescription, ...
    reordered_gene_IDs.RefSeqProteinID, reordered_gene_IDs.UniProtID, ...
    'VariableNames', {'UMAP order', 'Locus Tag', 'Gene Symbol', 'Entrez Gene ID', 'Protein Description', 'RefSeq Protein ID','UniProt ID'});

%% ==== 9. Save output Excel ==== %%
writetable(AA_table, outputFile, 'Sheet', 'AA Usage', 'WriteRowNames', true);
writetable(C_abs_table, outputFile, 'Sheet', 'Absolute Usage', 'WriteRowNames', true);
writetable(C_rel_table, outputFile, 'Sheet', 'Relative Usage', 'WriteRowNames', true);
writetable(newGeneTable, outputFile, 'Sheet', 'Reordered Genes', 'WriteRowNames', false);
writetable(cell2table(codons_reorder','VariableNames',{'Reordered Codons'}), outputFile, 'Sheet', 'Reordered Codons','WriteRowNames', false);

fprintf('Pipeline complete. Results saved with prefix: %s\n', strainPrefix);
