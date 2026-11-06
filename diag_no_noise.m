% diag_no_noise.m -- 判定 η 是噪声还是结构性失配
clear; clc;
if exist('OCTAVE_VERSION','builtin'), pkg load control optim; end

Params = init_params();
Params.noise_var  = 0;             % 关闭识别数据 AND 闭环噪声
Params.n_arx      = 6;
Params.arx_lambda = 1e-2;

if exist('rng','builtin'), rng(Params.rng_seed);
else, rand('seed',Params.rng_seed); randn('seed',Params.rng_seed); end

[Ident_data, Params] = gen_prbs_excitation(Params);
[Val_data,   Params] = gen_validation_data(Params, 'val', Ident_data);
Model  = build_dypls(Ident_data, Params);

[~, nrmse_free, ~, nrmse_1step, ~] = free_run_validation(Model, Val_data, Params);
fprintf('\n=== NO-NOISE DIAGNOSTIC ===\n');
fprintf('1-step   NRMSE  = [%s]\n', num2str(nrmse_1step, '%.4e '));
fprintf('free-run NRMSE  = [%s]\n', num2str(nrmse_free,   '%.4e '));

UncSets = estimate_uncertainty(Model, Ident_data, Params);
fprintf('eta_bound = [%s]\n', num2str(UncSets.eta_bound', '%.4e '));
