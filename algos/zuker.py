import numpy as np
from lookup import Lookup

class Zuker:
    def __init__(self, seq):
        self.seq = seq
        self.n = len(seq)
        self.lookup = Lookup()

        self.V_matrix = None
        self.W_matrix = None
        self.WM_matrix = None
        self.mfe = None

        self.create_V()
        self.create_W()

    def calc_hairpin(self, i, j):
        hairpin = np.inf
        size = j - i - 1
        if size > 3:
            hairpin = self.lookup.hairpin(size)
        return hairpin

    def calc_stacking(self, i, j):
        return self.lookup.stack(self.seq[i], self.seq[i+1], self.seq[j], self.seq[j-1])

    def calc_internal(self, i, j, V):
        internal = np.inf
        for k in range(i+1, j):
            for l in range(k+1, j):
                if self.is_valid_pair(self.seq[k], self.seq[l]):
                    a = k - i - 1
                    b = j - l - 1
                    if a == 0 or b == 0:
                        loopE = self.lookup.bulge(a + b)
                    else:
                        loopE = self.lookup.internal(a + b)
                    internal = min(internal, V[k][l] + loopE)
        return internal

    def create_V(self):
        V_matrix = np.full((self.n, self.n), np.inf)
        WM_matrix = np.full((self.n, self.n), np.inf)
        m = 3
        for i in range(self.n - 1, -1, -1):
            for j in range(i + m, self.n):
                if not self.is_valid_pair(self.seq[i], self.seq[j]):
                    continue   
                multiloop = np.inf
                # non_closed used for the WM calculation after
                non_closed = np.inf
                for k in range(i + 1, j):
                    multiloop = min(multiloop, WM_matrix[i+1, k] + WM_matrix[k+1, j-1] + 9.3)
                    non_closed = min(non_closed, WM_matrix[i, k] + WM_matrix[k+1, j]) 
                # hairpin
                hairpin = self.calc_hairpin(i, j)
                # stacking
                stacking = self.calc_stacking(i, j)
                # internal + buldge loops
                internal = self.calc_internal(i, j, V_matrix)
                # multiloop
                # assign V_matrix at (i, j)
                V_matrix[i, j] = min(hairpin, stacking, internal, multiloop)
                j_unpaired = WM_matrix[i, j-1]
                i_unpaired = WM_matrix[i+1, j]
                closed = V_matrix[i, j] - 0.6
                # assign WM at (i, j)
                WM_matrix[i, j] = min(j_unpaired, i_unpaired, closed, non_closed)
        self.V_matrix = V_matrix

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
        W_matrix = np.zeros((self.n, self.n))
        for i in range(self.n - 1, -1, -1):
            for j in range(i + 3, self.n):
                j_paired = np.inf
                for k in range(i, j - 3):
                    j_paired = min(j_paired, W_matrix[i, k-1] + self.V[k, j])
                W_matrix[i, j] = min(W_matrix[i, j-1], j_paired)
        self.W_matrix = W_matrix
        self.mfe = self.W_matrix[0, self.n - 1]

    def W_backtrace(self, W_pointers):
        node = W_pointers{(1,self.n)}
        # store (j, 'v/w')
        while node[0] > 1:


if __name__ == "__main__":
    RNA = 'AUAUAUAU'
    z = Zuker(RNA)
    print('V shape:', z.V.shape)
    print('W shape:', z.W.shape)
    print('MFE': z.mfe)