"""
Generate shell script for running power analysis

Inputs:
    --file: file path relative to lambda-BF/power-analysis, must contain a single .trees file
    --method: string, one of either 'lambdaBF', 'phylosig' or 'phylosig-mcc' (default 'lambdaBF')
    --n_draws: number of data sets to draw (default 200)
    --n_samples: number of importance sampling samples (default 100)
    --db: for lambdaBF, the log10 Bayes factor threshold (default 0)
    --verbose: whether to print progress (default FALSE)



the result is a shell script that can be run in the power-analysis directory
the script will save the output as Rdata objects for each value of lambda
in the same directory as the original tree file
"""

import os
import argparse
import itertools


def find_trees_file(file_path):
    """
    Find the .trees file in the given directory.

    Args:
        file_path (str): The path to the directory to search for .trees files.
    """
    trees_list = []
    for f in os.listdir(file_path):
        if f.endswith('.trees'):
            trees_list.append(f)
    return trees_list


def parse_args():
    parser = argparse.ArgumentParser(description="Generate shell script for running power analysis")
    parser.add_argument('--folder', type=str, required=True, help='path relative to lambda-BF/power-analysis, must contain a single .trees file')
    parser.add_argument('--method', type=str, default='lambdaBF', choices=['lambdaBF', 'phylosig', 'phylosig-mcc'], help='Method to use for power analysis (default: lambdaBF)')
    parser.add_argument('--n_draws', type=int, default=200, help='no. of data sets to draw (default: 200)')
    parser.add_argument('--n_samples', type=int, default=100, help='no. of importance sampling samples (default: 100)')
    parser.add_argument('--db', type=float, default=0, help='Log10 decision boundary for a finding in favour of Pagel model (default: 0)')
    parser.add_argument('--verbose', type=str, default='FALSE', choices=['TRUE', 'FALSE'], help='print progress? (default: FALSE)')
    parser.add_argument('--burnin', type=int, default=0, help='burnin for MCMC (default: 0)')
    parser.add_argument('--downsample', type=int, default=1, help='downsamping interval for MCMC (default: 0)')
    return parser.parse_args()

args = parse_args()

if args.folder[-1] == '/':
    args.folder = args.folder[:-1]

# if p-value used, change decision boundary to 0.05
if args.method != 'lambdaBF':
    args.db = 0.05


# write script
if not os.path.exists(f'{args.folder}'):
    print(f"Directory {args.folder} does not exist")
    exit(1)

os.makedirs(f'{args.folder}/scripts', exist_ok=True)
os.makedirs(f'{args.folder}/output', exist_ok=True)

lam_range = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
trees_list = find_trees_file(args.folder)

total_runs = len(lam_range) * len(trees_list)


command = "Rscript power_analysis/pa-script-{}.r --trees {}/{} --verbose {} --lam {} --decision_boundary {} --N_draws {} --N_samples {} --burnin {} --downsample {} --output {}"

with open(f'{args.folder}/scripts/{args.method}_{trees_list[0][:-6]}_{args.n_draws}.sh', 'w') as f:
    f.write("#!/bin/bash\n")
    i = 1
    for lam, tree in itertools.product(lam_range, trees_list):

        output_file = f"{args.folder}/output/{args.method}_{tree[:-6]}_lam{lam}_{args.n_draws}.Rdata"

        f.write(command.format(args.method, args.folder, tree, args.verbose, lam, args.db,
                               args.n_draws, args.n_samples, args.burnin, args.downsample,
                               output_file))
        f.write("\n")

        f.write(f'echo "Step {i}/{total_runs}, {100*i/total_runs:.0f}%:"')
        f.write("\n")
        i+=1