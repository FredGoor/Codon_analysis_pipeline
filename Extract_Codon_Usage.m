%% This script takes as input a cds exported FASTA file from NCBI and extracts DNA sequences and codon usage. It then creates a table with gene names and codon names %%
% 1. Loads a FASTA file exported from NCBI
% 2. Extracts codon usage for each CDS
% 3. Parses each header to get: Locus tag - Gene symbol - Entrez Gene ID - Protein description - RefSeq protein ID (NP_...) - UniProtKB/Swiss-Prot accession (Pxxxxxx, etc.)
% 4. Removes duplicated locus tags
% 5. Exports two Excel files: one with codon frequencies and the other with gene ID metadata

%% Codon Usage + Gene Metadata Extractor for NCBI CDS FASTA Files %%
clear all; clc;

%% 1. Import FASTA file
[file, folder] = uigetfile({'*.*','All Files (*.*)'}, ...
    'Select a File', ...
    'D:\Dropbox\Boulot Fred\Dirk - HGT and tRNAs\2025 Codon analysis different species\');
filename = fullfile(folder, file);
FASTAdata = fastaread((filename), 'blockread', [1 Inf]);
N_genes = length(FASTAdata);

%% 2. Compute codon usage
for i = 1:N_genes
    Codons(i) = codoncount(FASTAdata(i));
end

%% 3. Generate codon frequency table
T = struct2table(Codons);
N_aa = sum(T{ :, :},2);
F = T{:, :} ./ N_aa;
Freq = array2table(F);

%% 4. Add codon column names
SeqAA = aminolookup(nt2aa( T.Properties.VariableNames));
Cod_names = append(SeqAA.', '_', T.Properties.VariableNames);
Freq.Properties.VariableNames = Cod_names;

%% 5. Parse gene metadata from FASTA headers
header = {FASTAdata(:).Header}';

% Preallocate metadata arrays
locusTags = cell(size(header));
entrezIDs = cell(size(header));
geneSymbols = cell(size(header));
proteinDescriptions = cell(size(header));
proteinIDs = cell(size(header));
uniprotIDs = cell(size(header));

for i = 1:length(header)
    currentHeader = header{i};

    % ----- Locus Tag -----
    token = regexp(currentHeader, '\[locus_tag=([^\]]+)\]', 'tokens', 'once');
    if ~isempty(token)
        locusTags{i} = token{1};
    else
        locusTags{i} = 'NA';
    end

    % ----- Entrez Gene ID -----
    token = regexp(currentHeader, 'GeneID:(\d+)', 'tokens', 'once');
    if ~isempty(token)
        entrezIDs{i} = token{1};
    else
        entrezIDs{i} = 'NA';
    end

    % ----- Gene Symbol -----
    token = regexp(currentHeader, '\[gene=([^\]]+)\]', 'tokens', 'once');
    if ~isempty(token)
        geneSymbols{i} = token{1};
    else
        geneSymbols{i} = 'NA';
    end

    % ----- Protein Description -----
    token = regexp(currentHeader, '\[protein=([^\]]+)\]', 'tokens', 'once');
    if ~isempty(token)
        proteinDescriptions{i} = token{1};
    else
        proteinDescriptions{i} = 'NA';
    end

    % ----- RefSeq Protein ID -----
    token = regexp(currentHeader, '\[protein_id=([^\]]+)\]', 'tokens', 'once');
    if ~isempty(token)
        proteinIDs{i} = token{1};
    else
        proteinIDs{i} = 'NA';
    end

    % ----- UniProtKB Swiss-Prot Accession -----
    token = regexp(currentHeader, 'UniProtKB/Swiss-Prot:([A-Z0-9]+)', 'tokens', 'once');
    if ~isempty(token)
        uniprotIDs{i} = token{1};
    else
        uniprotIDs{i} = 'NA';
    end
end


%% 6. Remove duplicated locus tags (keep first occurrence)
[uniqueTags, firstIdx] = unique(locusTags, 'stable');

% Filter all variables accordingly
Freq = Freq(firstIdx, :);
locusTags = locusTags(firstIdx);
geneSymbols = geneSymbols(firstIdx);
entrezIDs = entrezIDs(firstIdx);
proteinDescriptions = proteinDescriptions(firstIdx);
proteinIDs = proteinIDs(firstIdx);
uniprotIDs = uniprotIDs(firstIdx);

% Set locus tags as row names in codon frequency table
Freq.Properties.RowNames = locusTags;

%% 7. Export codon frequency table
[tempDir, tempFile] = fileparts(file);
filename_codons = fullfile(folder, [tempFile, '_codons.xlsx']);
writetable(Freq, filename_codons, 'WriteRowNames', true);

%% 8. Export gene metadata table
geneIDs = [locusTags, geneSymbols, entrezIDs, proteinDescriptions, proteinIDs, uniprotIDs];
geneIDs_table = array2table(geneIDs);
geneIDs_table.Properties.VariableNames = {'LocusTag', 'GeneSymbol', 'EntrezGeneID', ...
    'ProteinDescription', 'RefSeqProteinID', 'UniProtID'};
geneIDs_table.Properties.RowNames = {};

filename_metadata = fullfile(folder, [tempFile, '_geneIDs.xlsx']);
writetable(geneIDs_table, filename_metadata);

%% 9. Compute GC Content
% Concatenate all sequences
fullSeq = upper(strjoin({FASTAdata(:).Sequence}, ''));
GC_count = count(fullSeq, 'G') + count(fullSeq, 'C');
GC_percent = 100 * GC_count / length(fullSeq);

% Display
fprintf('GC content of concatenated CDS: %.2f%%\n', GC_percent);

% Save as TXT
gc_file = fullfile(folder, 'GC_content.txt');
fid = fopen(gc_file, 'w');
fprintf(fid, 'GC content of CDS: %.2f%%\n', GC_percent);
fclose(fid);

%% 10. Completion message
fprintf('Codon usage table has been generated.\n');
fprintf('GC content saved to %s\n', gc_file);



