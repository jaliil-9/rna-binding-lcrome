function run_lcrfinder_batch(input_fasta, output_dir, resolution)
% Run LCRFinder on a protein FASTA and write interval-level CSV and FASTA outputs.
%
% Usage:
%   run_lcrfinder_batch('proteins.fasta', 'lcrfinder_results', 'Medium')
%
% Resolution: 'Low', 'Medium' (default), or 'High'.

if nargin < 3
    resolution = 'Medium';
end

root_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(root_dir, 'source'));

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

[headers, sequences] = read_fasta(input_fasta);
csv_path = fullfile(output_dir, 'LCRFinder.csv');
fasta_path = fullfile(output_dir, 'LCRFinder.fasta');

csv_file = fopen(csv_path, 'w');
fasta_file = fopen(fasta_path, 'w');
if csv_file == -1 || fasta_file == -1
    error('Could not create output files.');
end

fprintf(csv_file, 'protein_id,header,method,start,end,length,sequence,description\n');
interval_count = 0;

for i = 1:numel(sequences)
    sequence = upper(sequences{i});
    lcrs = detect_lcrs(sequence, resolution);
    protein_id = first_token(headers{i});

    for r = 1:numel(lcrs)
        start_pos = lcrs(r).cov(1);
        end_pos = lcrs(r).cov(end);
        lcr_sequence = lcrs(r).seq;
        description = sprintf('resolution=%s; k<=3', resolution);

        fprintf(csv_file, '%s,%s,%s,%d,%d,%d,%s,%s\n', ...
            csv_value(protein_id), ...
            csv_value(headers{i}), ...
            csv_value('LCRFinder'), ...
            start_pos, ...
            end_pos, ...
            end_pos - start_pos + 1, ...
            csv_value(lcr_sequence), ...
            csv_value(description));

        fprintf(fasta_file, '>%-s method=LCRFinder start=%d end=%d length=%d resolution=%s\n%s\n', ...
            protein_id, start_pos, end_pos, end_pos - start_pos + 1, resolution, lcr_sequence);

        interval_count = interval_count + 1;
    end

    fprintf('Processed %d/%d proteins; intervals=%d\n', i, numel(sequences), interval_count);
end

fclose(csv_file);
fclose(fasta_file);
fprintf('Finished: %d proteins, %d LCR intervals.\n', numel(sequences), interval_count);
fprintf('CSV: %s\nFASTA: %s\n', csv_path, fasta_path);
end


function lcrs = detect_lcrs(sequence, resolution)
% Protein-specific LCRFinder logic, equivalent to run_LCRFinder(..., 'AA', resolution)
% without its per-protein auxiliary text-file output.

A = 20;
K = 3;
LIMIT = 1000;
N_rand = 0.001;
[min_rep, LIMIT] = determine_search_parameters(LIMIT, N_rand, A, K);

CO = analyze_CO(sequence, K, min_rep, LIMIT);
if isempty(CO)
    lcrs = struct([]);
    return;
end

lcrs = analyze_LCR(sequence, CO, K, A, max(LIMIT), N_rand, resolution);
end


function [headers, sequences] = read_fasta(path)
fid = fopen(path, 'r');
if fid == -1
    error('Could not open FASTA file: %s', path);
end

headers = {};
sequences = {};
current_header = '';
current_sequence = '';

while true
    line = fgetl(fid);
    if ~ischar(line)
        break;
    end

    line = strtrim(line);
    if isempty(line)
        continue;
    end

    if line(1) == '>'
        if ~isempty(current_header)
            headers{end + 1} = current_header;
            sequences{end + 1} = current_sequence;
        end
        current_header = line;
        current_sequence = '';
    else
        current_sequence = [current_sequence regexprep(line, '\s+', '')];
    end
end

if ~isempty(current_header)
    headers{end + 1} = current_header;
    sequences{end + 1} = current_sequence;
end

fclose(fid);

if isempty(sequences)
    error('No FASTA records found in: %s', path);
end
end


function token = first_token(header)
token = regexp(header, '^>(\S+)', 'tokens', 'once');
if isempty(token)
    token = header;
else
    token = token{1};
end
end


function value = csv_value(text)
text = strrep(char(text), '"', '""');
value = ['"' text '"'];
end
