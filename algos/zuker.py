import numpy as np
from lookup import Lookup

def calc_hairpin(i, j, lookup):
    hairpin = np.inf
    size = j - i - 1
    if size > 3:
        hairpin = lookup.hairpin(size)
    return hairpin

def calc_stacking(i, j, seq, lookup):
    return lookup.stack(seq[i], seq[i+1], seq[j], seq[j-1])

def calc_internal(i, j, V, lookup):
    internal = np.inf
    for k in range(i+1, j):
        for l in range(k+1, j):
            if is_valid_pair(k, l):
                a = k - i - 1
                b = j - l - 1
                if a == 0 or b == 0:
                    loopE = lookup.bulge(a + b)
                else:
                    loopE = lookup.internal(a + b)
                internal = min(internal, V[k][l] + loopE)
    return internal

def create_V(seq):
    V_matrix = np.full((len(seq), len(seq)), np.inf)
    WM_matrix = np.full((len(seq), len(seq)), np.inf)
    m = 3
    n = len(seq)
    lookup = Lookup()
    for i in range(n - 1, -1, -1):
        for j in range(i + m, n):
            print(seq[i], seq[j])
            if not is_valid_pair(seq[i], seq[j]):
                continue   
            internal = np.inf
            multiloop = np.inf

            # hairpin
            hairpin = calc_hairpin(i, j, lookup)
            
            # stacking
            stacking = calc_stacking(i, j, seq, lookup)
            
            # internal + buldge loops
            internal = calc_internal(i, j, V_matrix, lookup)
            
            # multiloop
            # non_closed used for the WM calculation after
            non_closed = np.inf

            for k in range(i + 1, j):
                multiloop = min(multiloop, WM_matrix[i+1, k] + WM_matrix[k+1, j-1] + 9.3)
                non_closed = min(non_closed, WM_matrix[i, k] + WM_matrix[k+1, j]) 
            
            # assign V_matrix at (i, j)
            V_matrix[i, j] = min(hairpin, stacking, internal, multiloop)
            
            # assign WM at (i, j)
            j_unpaired = WM_matrix[i, j-1]
            i_unpaired = WM_matrix[i+1, j]
            closed = V_matrix[i, j] - 0.6
            WM_matrix[i, j] = min(j_unpaired, i_unpaired, closed, non_closed)
    return V_matrix

def is_valid_base(c):
    print(str(c))
    return str(c) in ('A', 'C', 'G', 'U')

def is_valid_base(c):
    print(str(c))
    return str(c) in ('A', 'C', 'G', 'U')

def is_valid_pair(c0, c1):
    print(str(c0), not is_valid_base(c0))
    if not is_valid_base(c0) or not is_valid_base(c1):
        raise Exception('Invalid RNA base')
    p = [c0, c1]
    p.sort()
    sp = "".join(p)
    return sp in ('CG', 'AU', 'GU')

def create_W(seq, V):
    W_matrix = np.zeros((len(seq), len(seq)))
    for i in range(len(seq) - 1, -1, -1):
        for j in range(i + 3, len(seq)):
            j_paired = np.inf
            for k in range(i, j - 3):
                j_paired = min(j_paired, W_matrix[i, k-1] + V[k, j])
            W_matrix[i, j] = min(W_matrix[i, j-1], j_paired)
    return W_matrix

def zukers(seq):
    V = create_V(seq)
    W = create_W(seq, V)
    return W[0, len(seq) - 1]

if __name__ == "__main__":
    true_w = -1.8 # change
    RNA = 'AUAUAUAU'
    V = create_V(RNA)
    assert true_w == V[0, 0], 'Wrong W[0][0]'
    print(V.shape)
    W = create_W(RNA, V)
    print(W.shape)
    print(zukers(RNA))
    # RNA = 'AUCGCAU'
    # print(create_V(RNA).shape)