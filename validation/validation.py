from algos.zuker import Zuker
import pandas as pd
import matplotlib.pyplot as plt

def get_pairs(struc):
    stack = []
    pairs = []

    for idx, c in enumerate(struc):
        if c == '(':
            stack.append(idx)
        elif c == ')':
            pairs.append((stack[-1], idx))
            stack.pop()
    return pairs

def false_pos(z_seq, t_seq):
    return list(set(z_seq) - set(t_seq))

def false_neg(z_seq, t_seq):
    return list(set(t_seq) - set(z_seq))

def true_pos(z_seq, t_seq):
    return list(set(z_seq) & set(t_seq))

# install huggingface: pip install --upgrade huggingface_hub
# from https://huggingface.co/datasets/multimolecule/archiveii.512
df = pd.read_parquet("hf://datasets/multimolecule/archiveii.512/test.parquet")
test_data = [(seq, struc) for seq, struc in zip(df['sequence'], df['secondary_structure'])]
print(len(test_data))
count = 1

prec_stats = []
rec_stats = []
f1_stats = []

for seq, struc in test_data:
    z = Zuker(seq)
    true_struc = get_pairs(struc)
    pred_struc = z.pairs
    true_pos_count = len(true_pos(pred_struc, true_struc))
    false_pos_count = len(false_pos(pred_struc, true_struc))
    false_neg_count = len(false_neg(pred_struc, true_struc))

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
    
    print(f'string {count} done')
    count += 1

    prec_stats.append(pres)
    rec_stats.append(rec)
    f1_stats.append(f1)

fig, axs = plt.subplots(1, 3, figsize=(12,4), sharex=True, sharey=True)
axs[0].hist(prec_stats, bins=20, color='blue')
axs[0].set_title("Precision")
axs[0].set_xlabel("Score")
axs[0].set_ylabel("Frequency")

axs[1].hist(rec_stats, bins=20, color='green')
axs[1].set_title("Recall")
axs[1].set_xlabel("Score")


axs[2].hist(f1_stats, bins=20, color='red')
axs[2].set_title("F1 Score")
axs[2].set_xlabel("Score")

plt.tight_layout()
plt.show()