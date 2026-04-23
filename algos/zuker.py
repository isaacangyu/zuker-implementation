import numpy as np
from lookup import Lookup

def create_V(seq):
    V_matrix = np.full((len(seq), len(seq)), np.inf)
    WM_matrix = np.full(len(seq), len(seq), np.inf)
    lookup = Lookup()
    for i in range(len(seq) - 1, -1, -1):
        for j in range(i + 3, len(seq)):
            if not is_valid_pair(seq[i], seq[j]):
                continue
            
            hairpin = np.inf
            stacking = np.inf
            internal = np.inf
            multiloop = np.inf

            # hairpin
            size = j - i - 1
            if size > 3:
                hairpin = lookup.hairpin(size)
            
            # stacking
            stacking = lookup.stack(seq[i], seq[i+1], seq[j], seq[j-1])
            
            # internal + buldge loops
            for k in range(i+1, j):
                for l in range(k+1, j):
                    if is_valid_pair(k, l):
                        a = k - i - 1
                        b = j - l - 1
                        if a == 0 or b == 0:
                            loopE = lookup.bulge(a + b)
                        else:
                            loopE = lookup.bulge(a + b)

                        internal = min(internal, V_matrix[k][l] + loopE)
            
            # multiloop
            # non_closed used for the WM calculation after
            non_closed = np.inf

            for k in range(i + 1, j):
                multiloop = min(multiloop, WM_matrix[i+1, k] + WM_matrix[k+1, j-1] + 9.3)
                non_closed = min(non_closed, WM_matrix[i, k] + WM_matrix[k+1, j]) 
            
            # assign V_matrix at (i, j)
            V_matrix[i][j] = min(hairpin, stacking, internal, multiloop)
            
            # assign WM at (i, j)
            j_unpaired = WM_matrix[i][j-1]
            i_unpaired = WM_matrix[i+1][j]
            closed = V_matrix[i][j] - 0.6
            WM_matrix = min(j_unpaired, i_unpaired, closed, non_closed)

def is_valid_pair(c0, c1):
    return (c0 == 'C' and c1 == 'G') or (c0 == 'G' and c1 == 'C') or (c0 == 'A' and c1 == 'U') or (c0 == 'U' and c1 == 'A') or (c0 == 'G' and c1 == 'U') or (c0 == 'U' and c1 == 'G')

def create_W(seq, V):
    W_matrix = np.zeros(len(seq), len(seq))
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

#if '__name__' == 'main':
RNA = 'AUAUAUAU'
print(create_V(RNA).shape)