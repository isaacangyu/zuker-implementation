import numpy as np
from lookup import Lookup

def create_V(seq):
    V_matrix = np.full((len(seq), len(seq)), np.inf)
    VM_matrix = np.full((len(seq), len(seq)), np.inf)
    lookup = Lookup()
    m = 3
    for i in range(len(seq) - 1, -1, -1):
        for j in range(i + m, len(seq)):
            print(i,j)
            if not is_valid_pair(seq[i], seq[j]):
                V_matrix[i][j] = float('inf')
                continue
            
            hairpin = float('inf')
            stacking = float('inf')
            internal = float('inf')

            size = j - i - 1
            # hairpin
            # if size > 3:
            hairpin = lookup.hairpin(size)
            print('hairpin', hairpin)
            
            # stacking
            stacking = lookup.stack(seq[i + 1], seq[j - 1], seq[i], seq[j]) + V_matrix[i + 1][j - 1]
            print('stack', stacking)
            # internal + bulge loops
            for k in range(i+1, j):
                for l in range(k+1, j):
                    if is_valid_pair(seq[k], seq[l]):
                        a = k - i - 1
                        b = j - l - 1
                        if a == 0 or b == 0:
                            loopE = lookup.bulge(b - a)
                        else:
                            loopE = lookup.internal(b - a)  # include seq
                        internal = min(internal, V_matrix[k][l] + loopE)
                        print('internal', internal)
            
            # multiloop
            """
            j_unpair = VM[i, j - 1]
            i_unpair = VM[i + 1, j]
            closed = V
            """

            V_matrix[i][j] = min(hairpin, stacking, internal)

    print(V_matrix)
    return V_matrix

def is_valid_pair(c0, c1):
    p = "".join(sorted(c0 + c1))
    return p in ('CG', 'AU', 'GU')

def W(i, j):
    """assuming closed loop"""
    pass

if __name__ == "__main__":
    true_w = -1.8 # change
    RNA = 'AUAUAUAU'
    V = create_V(RNA)
    assert true_w == V[0][0], 'Wrong W[0][0]'
    print(V.shape)
    # RNA = 'AUCGCAU'
    # print(create_V(RNA).shape)