#!/bin/bash
Rscript power_analysis/pa-script-lambdaBF.r --trees power_analysis/test/example.trees --verbose FALSE --lam 0 --decision_boundary 0 --N_draws 10 --N_samples 100 --burnin 0 --downsample 1 --output power_analysis/test/output/lambdaBF_example_lam0_10.Rdata
echo "Step 1/11, 9%:"
Rscript power_analysis/pa-script-lambdaBF.r --trees power_analysis/test/example.trees --verbose FALSE --lam 0.1 --decision_boundary 0 --N_draws 10 --N_samples 100 --burnin 0 --downsample 1 --output power_analysis/test/output/lambdaBF_example_lam0.1_10.Rdata
echo "Step 2/11, 18%:"
Rscript power_analysis/pa-script-lambdaBF.r --trees power_analysis/test/example.trees --verbose FALSE --lam 0.2 --decision_boundary 0 --N_draws 10 --N_samples 100 --burnin 0 --downsample 1 --output power_analysis/test/output/lambdaBF_example_lam0.2_10.Rdata
echo "Step 3/11, 27%:"
Rscript power_analysis/pa-script-lambdaBF.r --trees power_analysis/test/example.trees --verbose FALSE --lam 0.3 --decision_boundary 0 --N_draws 10 --N_samples 100 --burnin 0 --downsample 1 --output power_analysis/test/output/lambdaBF_example_lam0.3_10.Rdata
echo "Step 4/11, 36%:"
Rscript power_analysis/pa-script-lambdaBF.r --trees power_analysis/test/example.trees --verbose FALSE --lam 0.4 --decision_boundary 0 --N_draws 10 --N_samples 100 --burnin 0 --downsample 1 --output power_analysis/test/output/lambdaBF_example_lam0.4_10.Rdata
echo "Step 5/11, 45%:"
Rscript power_analysis/pa-script-lambdaBF.r --trees power_analysis/test/example.trees --verbose FALSE --lam 0.5 --decision_boundary 0 --N_draws 10 --N_samples 100 --burnin 0 --downsample 1 --output power_analysis/test/output/lambdaBF_example_lam0.5_10.Rdata
echo "Step 6/11, 55%:"
Rscript power_analysis/pa-script-lambdaBF.r --trees power_analysis/test/example.trees --verbose FALSE --lam 0.6 --decision_boundary 0 --N_draws 10 --N_samples 100 --burnin 0 --downsample 1 --output power_analysis/test/output/lambdaBF_example_lam0.6_10.Rdata
echo "Step 7/11, 64%:"
Rscript power_analysis/pa-script-lambdaBF.r --trees power_analysis/test/example.trees --verbose FALSE --lam 0.7 --decision_boundary 0 --N_draws 10 --N_samples 100 --burnin 0 --downsample 1 --output power_analysis/test/output/lambdaBF_example_lam0.7_10.Rdata
echo "Step 8/11, 73%:"
Rscript power_analysis/pa-script-lambdaBF.r --trees power_analysis/test/example.trees --verbose FALSE --lam 0.8 --decision_boundary 0 --N_draws 10 --N_samples 100 --burnin 0 --downsample 1 --output power_analysis/test/output/lambdaBF_example_lam0.8_10.Rdata
echo "Step 9/11, 82%:"
Rscript power_analysis/pa-script-lambdaBF.r --trees power_analysis/test/example.trees --verbose FALSE --lam 0.9 --decision_boundary 0 --N_draws 10 --N_samples 100 --burnin 0 --downsample 1 --output power_analysis/test/output/lambdaBF_example_lam0.9_10.Rdata
echo "Step 10/11, 91%:"
Rscript power_analysis/pa-script-lambdaBF.r --trees power_analysis/test/example.trees --verbose FALSE --lam 1.0 --decision_boundary 0 --N_draws 10 --N_samples 100 --burnin 0 --downsample 1 --output power_analysis/test/output/lambdaBF_example_lam1.0_10.Rdata
echo "Step 11/11, 100%:"
