import numpy as np
from lookup import Lookup

class Zuker:
    def __init__(self, seq, min_loop=3, offset=9.3, helix=-0.6, unpaired_nuc=0):
        self.seq = seq
        self.n = len(seq)
        self.lookup = Lookup()

        self.V = None
        self.W = None
        self.WM = None
        
        self.W_pointers = None
        self.V_pointers = None
        self.WM_pointers = None

        self.m = min_loop           # must be at least 3
        self.a = offset
        self.b = helix
        self.c = unpaired_nuc

        self.mfe = None

        self.create_V()
        self.create_W()

    def calc_hairpin(self, i, j):
        hairpin = np.inf
        size = j - i - 1
        if size >= self.m:
            hairpin = self.lookup.hairpin(size)
        return hairpin

    def calc_stacking(self, i, j):
        if j - i < 2:                           # check if (i+1, j-1) are valid
            return np.inf
        
        return self.lookup.stack(self.seq[i], self.seq[i+1], self.seq[j], self.seq[j-1])

    def calc_internal(self, i, j, V):
        internal = np.inf
        MAX = 30                                # restrict loop size to 30 to keep runtime O(n^3)
        k_best = None
        l_best = None
        for k in range(i+1, j):
            for l in range(k+1, j):
                if self.is_valid_pair(self.seq[k], self.seq[l]) and l - k > self.m:
                    if k == i+1 and l == j-1:   # skip stacking pair case
                        continue

                    a = k - i - 1
                    b = j - l - 1
                    
                    if a + b > MAX:
                        continue

                    if a == 0 or b == 0:
                        loopE = self.lookup.bulge(a + b)
                    else:
                        loopE = self.lookup.internal(a + b)

                    if internal > V[k, l] + loopE:
                        internal = V[k, l] + loopE
                        k_best = k
                        l_best = l
                        
        return internal, k_best, l_best

    def create_V(self):
        V = np.full((self.n, self.n), None)
        V[np.tril_indices(self.n)] = np.inf     # create V and set all indices with j <= i to inf
        WM = np.full((self.n, self.n), None)
        WM[np.tril_indices(self.n)] = np.inf    # create WM and set all indices with j <= i to inf
        V_pointers = np.full((self.n, self.n), None)
        WM_pointers = np.full((self.n, self.n), None)

        for l in range(1, self.n):
            for i in range(self.n - l):
                j = l + i

                if j - i <= self.m:                  # set to inf if (i, j) are to close
                    V[i, j] = np.inf
                    WM[i, j] = np.inf
                    continue
                
                # V at (i, j)
                if not self.is_valid_pair(self.seq[i], self.seq[j]):
                    V[i, j] = np.inf
                else:
                    # hairpin
                    hairpin = self.calc_hairpin(i, j)

                    # stacking
                    stacking = self.calc_stacking(i, j) + V[i+1, j-1]

                    # internal + buldge loops
                    internal, k_best_int, l_best_int = self.calc_internal(i, j, V)

                    # multiloop
                    multiloop = np.inf
                    k_best_mul = None
                    for k in range(i + 1, j):
                        if multiloop > WM[i+1, k] + WM[k+1, j-1] + self.a:
                            multiloop = WM[i+1, k] + WM[k+1, j-1] + self.a
                            k_best_mul = k
                    
                    best = np.inf
                    best_pointer = None

                    if best > hairpin:
                        best = hairpin
                        best_pointer = ("H",)

                    if best > stacking:
                        best = stacking
                        best_pointer = ("S",)
                    
                    if best > internal:
                        best = internal
                        best_pointer = ("I", k_best_int, l_best_int)
                    
                    if best > multiloop:
                        best = multiloop
                        best_pointer = ("M", k_best_mul)
                    
                    V[i, j] = best
                    V_pointers[i, j] = best_pointer
                 
                # WM at (i, j)
                j_unpaired = WM[i, j-1] + self.c

                i_unpaired = WM[i+1, j] + self.c

                closed = V[i, j] + self.b

                non_closed = np.inf
                k_best = None
                for k in range(i + 1, j):
                    if non_closed > WM[i, k] + WM[k+1, j]:
                        non_closed = WM[i, k] + WM[k+1, j]
                        k_best = k
                
                best = np.inf
                best_pointer = None

                if best > j_unpaired:
                    best = j_unpaired
                    best_pointer = ("Uj",)

                if best > i_unpaired:
                    best = i_unpaired
                    best_pointer = ("Ui",)

                if best > closed:
                    best = closed
                    best_pointer = ("C",)
                
                if best > non_closed:
                    best = non_closed
                    best_pointer = ("S", k_best)
                
                WM[i, j] = best
                WM_pointers[i, j] = best_pointer
                
        self.V = V
        self.V_pointers = V_pointers
        self.WM = WM
        self.WM_pointers = WM_pointers

    def is_valid_base(self, c):
        return str(c) in ('A', 'C', 'G', 'U')

    def is_valid_pair(self, c0, c1):
        if not self.is_valid_base(c0) or not self.is_valid_base(c1):
            raise Exception('Invalid RNA base')
        p = [c0, c1]
        p.sort()
        sp = "".join(p)
        return sp in ('CG', 'AU', 'GU')

    def create_W(self):
        W = np.full((self.n, self.n), np.inf)
        W[np.tril_indices(self.n)] = 0              # create W and set all indices with j <= i to 0
        W_pointers = np.full((self.n, self.n), None)
        for l in range(1, self.n):
            for i in range(self.n - l):
                j = l + i
                if j - i <= self.m:                  # set to 0 if (i, j) are to close
                    W[i, j] = 0
                    continue
                
                j_paired = np.inf
                k_best = None
                for k in range(i, j - self.m):
                    left_val = 0

                    if k != i:
                        left_val = W[i, k-1]

                    if j_paired > left_val + self.V[k, j]:
                        j_paired = left_val + self.V[k, j]
                        k_best = k
                
                if W[i, j-1] <= j_paired:
                    W[i, j] = W[i, j-1]
                    W_pointers[i, j] = ("U", j-1)
                else:
                    W[i, j] = j_paired
                    W_pointers[i, j] = ("P", k_best)
    
        self.W = W
        self.mfe = self.W[0, self.n - 1]
        self.W_pointers = W_pointers

"""

    def backtrace(self):
        j = self.n - 1
        node = self.W_pointers[j]
        bp_idxs = []
        
        while j > 1:
            k = node[1]
            
            self.bp_idxs.append self.V_backtrace(k, j)
            j = k
        
        self.dot = self.write_dot(bp_idxs)
        return self.dot

    def V_backtrace(self, k, j):
        self.V_pointers
        bp_id

    def write_dot(self, bp_idxs):
        nbp = len(bp_idxs)
        assert nbp % 2 == 0, "Number of base pairs is not even"
        dot = ['.' for _ in range(self.n)]

        for i in range(nbp / 2):
            dot[i] = '('
        for i in range(nbp / 2, nbp):
            dot[i] = ')'
        
        return dot
    
"""

if __name__ == "__main__":
    RNA = 'AUAUAUAUAU'
    z = Zuker(RNA)
    print('V shape:', z.V.shape)
    print('W shape:', z.W.shape)
    print('MFE:', z.mfe)
    print('V:')
    print(z.V)
    print("W:")
    print(z.W)
    print("V pointers:")
    print(z.V_pointers)
    print("W pointers:")
    print(z.W_pointers)
    # print('Dot:', z.backtrace())