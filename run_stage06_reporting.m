% =========================================================================
% run_stage06_reporting
% =========================================================================
clear; clc; close all;
if ~exist('OCTAVE_VERSION','builtin'), error('This simulation project targets GNU Octave.'); end

%% 1. Project directory & Diary
script_file = mfilename('fullpath');
if isempty(script_file), project_dir = pwd; else, project_dir = fileparts(script_file); end
addpath(project_dir); rehash();

results_dir = fullfile(project_dir,'results');
log_file = fullfile(results_dir, 'stage06_log.txt');
if exist(log_file, 'file'), delete(log_file); end
diary(log_file);

fprintf('=== Stage 6: reporting and export ===\n');

%% 2. Check dependencies
required = {'plot_results', 'export_all_figures', 'export_comparison_table', 'compute_metrics', 'report_simulation_comparison'};
for j = 1:numel(required)
    expected = fullfile(project_dir,[required{j},'.m']);
    if exist(expected,'file') ~= 2, error('Missing required file: %s',expected); end
end

%% 3. Load online results
input_file = fullfile(results_dir,'stage05_simulation.mat');
if exist(input_file,'file') ~= 2, error('Run run_stage05_simulation.m before Stage 6.'); end
Results = load(input_file);

if isfield(Results, 'Data_oh') && ~isempty(Results.Data_oh)
    Data_oh = Results.Data_oh;
    Metrics_oh = Results.Metrics_oh;
else
    warning('未找到 Oh 方法的数据，将跳过相关对比。');
    Data_oh = []; Metrics_oh = [];
end

tube_file = fullfile(results_dir, 'stage03_tube.mat');
if exist(tube_file, 'file') == 2
    S3 = load(tube_file, 'Invariant', 'StateError');
    Results.Invariant = S3.Invariant; Results.StateError = S3.StateError;
end

Params = Results.Params;
if ~strcmp(Params.code_version,'cryst_modular_v1'), error('Stage 6: incompatible artifact version.'); end

%% 4. Check simulation dimensions
datasets = {Results.Data_prop,Results.Data_xue};
for method = 1:2
    Data = datasets{method};
    if size(Data.Y_hist,1) ~= Params.n_out || size(Data.Y_hist,2) ~= Params.Ble+1 || size(Data.Y_hist,3) ~= Params.N_batches
        error('Stage 6: invalid output history dimensions.');
    end
    if size(Data.V_hist,1) ~= Params.n_in || size(Data.V_hist,2) ~= Params.Ble || size(Data.V_hist,3) ~= Params.N_batches
        error('Stage 6: invalid full physical input dimensions.');
    end
end

%% 5. Recompute metrics (包含了 Oh 方法)
old_prop = Results.Metrics_prop;
old_xue = Results.Metrics_xue;

Results.Metrics_prop = compute_metrics(Results.Data_prop,Params);
Results.Metrics_xue = compute_metrics(Results.Data_xue,Params);

if ~isempty(Data_oh)
    Results.Metrics_oh = compute_metrics(Results.Data_oh,Params);
else
    Results.Metrics_oh = [];
end

if isfield(old_prop,'total_simulation_time'), Results.Metrics_prop.total_simulation_time = old_prop.total_simulation_time; else, Results.Metrics_prop.total_simulation_time = NaN; end
if isfield(old_xue,'total_simulation_time'), Results.Metrics_xue.total_simulation_time = old_xue.total_simulation_time; else, Results.Metrics_xue.total_simulation_time = NaN; end
if ~isempty(Results.Metrics_oh) && isfield(Metrics_oh,'total_simulation_time')
    Results.Metrics_oh.total_simulation_time = Metrics_oh.total_simulation_time;
elseif ~isempty(Results.Metrics_oh)
    Results.Metrics_oh.total_simulation_time = NaN;
end

report_simulation_comparison(Results.Metrics_prop,Results.Metrics_xue,Results.Metrics_oh);

%% 6. Export tables
report_dir = fullfile(results_dir,'reporting');
if exist(report_dir,'dir') ~= 7, [ok,message] = mkdir(report_dir); if ~ok, error('Cannot create reporting directory: %s',message); end; end

Results.Data_oh = Data_oh;
Results.Metrics_oh = Metrics_oh;
CSVFiles = export_comparison_table(Results,report_dir);
export_tables(Results, report_dir);

%% 7. Plot and export
Figures = struct('handle',{},'name',{});
ExportReport = struct('figure_name',{},'format',{},'path',{},'success',{},'message',{});
PlotReport = struct('requested',Params.plot_results, 'generated',false, 'message','');

if Params.plot_results
    try
        Figures = plot_results(Results);
        PlotReport.generated = true;
    catch ME
        PlotReport.message = ME.message;
        fprintf('Plot generation failed: %s\n',ME.message);
    end
    if PlotReport.generated
        figure_dir = fullfile(report_dir,'figures');
        ExportReport = export_all_figures(Figures,figure_dir);
    end
else
    PlotReport.message = 'Plotting disabled by saved Params.plot_results.';
end

%% 8. Save reporting artifact
Artifact = struct('stage',6,'code_version',Params.code_version,'octave_version',OCTAVE_VERSION,'created_at',datestr(now,31), ...
    'source_file',input_file,'validation_provenance_verified',false,'csv_export_completed',true, ...
    'plot_generation_completed',PlotReport.generated,'figure_exports_completed',isempty(ExportReport) == false && all([ExportReport.success]));

report_file = fullfile(report_dir,'reporting_artifact.mat');
save('-mat7-binary',report_file, 'Artifact','CSVFiles','PlotReport','ExportReport','Metrics_prop','Metrics_xue');
if ~isempty(Data_oh), save('-append', report_file, 'Metrics_oh'); end

fprintf('\nSaved: %s\n',report_file);
fprintf('\nStage 6 completed.\n');
diary off;
