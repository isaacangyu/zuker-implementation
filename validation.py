from algos.zuker import Zuker
from algos.nussinov import Nussinov
import pandas as pd
import matplotlib.pyplot as plt
import json
import numpy as np

# takes a structure in the form of a string of the dot-parentheses format and returns a lists of lists of the bases paired
def get_pairs(struc):
    stack = []
    pairs = []

    for idx, c in enumerate(struc):
        if c == '(':
            stack.append(idx)
        elif c == ')':
            pairs.append([stack[-1], idx])
            stack.pop()
    return pairs

# calculates the false positive list
def false_pos(z_seq, t_seq):
    z_set = set(map(tuple, z_seq))
    t_set = set(map(tuple, t_seq))

    return list(z_set - t_set)

# calculates the false negative list
def false_neg(z_seq, t_seq):
    z_set = set(map(tuple, z_seq))
    t_set = set(map(tuple, t_seq))

    return list(t_set - z_set)

# calculates the true positive list
def true_pos(z_seq, t_seq):
    z_set = set(map(tuple, z_seq))
    t_set = set(map(tuple, t_seq))

    return list(z_set & t_set)

# creates the 'zuker_results.json' file with the results of the Zuker algorithm being ran on each sequence in the list
def init_zuker(seq_list):
    results = []

    for idx, seq in enumerate(seq_list):
        z = Zuker(seq)

        results.append({
            'sequence': seq,
            'mfe': z.mfe,
            'pairs': z.pairs,
            'dot': z.dot,
            'runtime': z.runtime
        })

        print(f'string {idx} completed')

    with open('zuker_results.json', 'w') as f:
        json.dump(results, f, indent=4)

# creates the 'nussinov_results.json' file with the results of the Nussinov algorithm being ran on each sequence in the list
def init_nussinov(seq_list):
    results = []

    for idx, seq in enumerate(seq_list):
        n = Nussinov(seq)

        results.append({
            'sequence': seq,
            'pairs': n.pairs,
        })

        print(f'string {idx} completed')

    with open('nussinov_results.json', 'w') as f:
        json.dump(results, f, indent=4)

# evaluates and precision, recall, and f1-score of structures and graphs them
def evaluate_acc(true_struct, pred_struct):
    prec_stats = []
    rec_stats = []
    f1_stats = [] 

    for t, p in zip(true_struct, pred_struct):
        true_pos_count = len(true_pos(p, t))
        false_pos_count = len(false_pos(p, t))
        false_neg_count = len(false_neg(p, t))

        try:
            pres = true_pos_count/(true_pos_count + false_pos_count)
        except ZeroDivisionError:
            pres = 0

        try:
            rec = true_pos_count/(true_pos_count + false_neg_count)
        except ZeroDivisionError:
            rec = 0
        
        try:
            f1 = 2 * (pres * rec) / (pres + rec)
        except ZeroDivisionError:
            f1 = 0
        
        prec_stats.append(pres)
        rec_stats.append(rec)
        f1_stats.append(f1)

    prec_stats = np.array(prec_stats)
    rec_stats = np.array(rec_stats)
    f1_stats = np.array(f1_stats)

    fig, axs = plt.subplots(1, 3, figsize=(12,4), sharex=True, sharey=True)
    axs[0].hist(prec_stats, bins=20, color='blue')
    axs[0].set_title('Precision')
    axs[0].set_xlabel('Score')
    axs[0].set_ylabel('Frequency')
    print(f'Precision: mean - {np.mean(prec_stats)}, median - {np.median(prec_stats)}')

    axs[1].hist(rec_stats, bins=20, color='green')
    axs[1].set_title('Recall')
    axs[1].set_xlabel('Score')
    print(f'Recall: mean - {np.mean(rec_stats)}, median - {np.median(rec_stats)}')


    axs[2].hist(f1_stats, bins=20, color='red')
    axs[2].set_title('F1 Score')
    axs[2].set_xlabel('Score')
    print(f'F1 Score: mean - {np.mean(f1_stats)}, median - {np.median(f1_stats)}')

    plt.tight_layout()
    plt.show()

# plots the seqeunce lengths vs runtime graphs
def plot_runtime(seq, runtimes):
    seq = np.array(seq)
    runtimes = np.array(runtimes)

    idx = np.argsort(seq)
    seq = seq[idx]
    runtimes = runtimes[idx]

    a = np.sum(runtimes * seq) / np.sum(seq**2)
    b = np.sum(runtimes * seq**2) / np.sum(seq**4)
    c = np.sum(runtimes * seq**3) / np.sum(seq**6)

    n_line  = a * seq
    n2_line = b * seq**2
    n3_line = c * seq**3

    fig, axs = plt.subplots(1, 2, figsize=(12,4))
    axs[0].scatter(seq, runtimes, label='Measured times')
    axs[0].plot(seq, n_line,  label='O(n) fit', color='red', linewidth=2)
    axs[0].plot(seq, n2_line, label='O(n²) fit', color='orange', linewidth=2)
    axs[0].plot(seq, n3_line, label='O(n³) fit', color='#2ca02c', linewidth=2)
    axs[0].set_xlabel('Length of Sequence')
    axs[0].set_ylabel('Runtime (sec)')
    axs[0].set_title('Runtime Comparisons')

    axs[1].plot(seq, runtimes / seq, label='runtime / n', color='red', linewidth=2)
    axs[1].plot(seq, runtimes / seq**2, label='runtime / n²', color='orange', linewidth=2)
    axs[1].plot(seq, runtimes / seq**3, label='runtime / n³', color='#2ca02c', linewidth=2)
    axs[1].set_xlabel('Length of Sequence')
    axs[1].set_ylabel('Runtime Ratio')
    axs[1].set_title('Runtime Ratio Comparisons')

    axs[0].legend()
    axs[1].legend()
    plt.show()

# install huggingface: pip install --upgrade huggingface_hub
# from https://huggingface.co/datasets/multimolecule/archiveii.512
df = pd.read_parquet('hf://datasets/multimolecule/archiveii.512/test.parquet')
test_seq = [seq for seq in df['sequence']]
test_struct = [get_pairs(struc) for struc in df['secondary_structure']]
print(len(test_seq))

with open('zuker_results.json', 'r') as f:
    results = json.load(f)

pred_pairs = [res['pairs'] for res in results]

evaluate_acc(test_struct, pred_pairs)       # calculates the precision, recall, and f1-score graphs

###########

# runtimes = [res['runtime'] for res in results]
# seq_len = [len(res['sequence']) for res in results]

# plot_runtime(seq_len, runtimes)           # plots the sequence list vs runtime graphs