import numpy as np

class Nussinov:
    def __init__(self, seq):
        self.seq = seq
        self.n = len(seq)
        self.M = None
        self.M_pointers = None
        self.pairs = None

        self.matrix()
        self.backtrace()

    
    def is_valid_base(self, c):
        return str(c) in ('A', 'C', 'G', 'U')

    def is_valid_pair(self, c0, c1):
        if not self.is_valid_base(c0) or not self.is_valid_base(c1):
            raise Exception('Invalid RNA base')
        p = [c0, c1]
        p.sort()
        sp = "".join(p)
        return sp in ('CG', 'AU', 'GU')

    def matrix(self):
        M = np.zeros((self.n, self.n), dtype=int)
        M_pointers = np.full((self.n, self.n), None)

        for l in range(1, self.n):
            for i in range(self.n - l):
                j = l + i

                best = M[i+1, j]
                M_pointers[i, j] = ('Ui', i+1, j)

                if M[i, j-1] > best:
                    best = M[i, j-1]
                    M_pointers[i, j] = ('Uj', i, j-1)

                pair_score = M[i+1, j-1] + (
                    1 if self.is_valid_pair(self.seq[i], self.seq[j]) else 0
                )
                if pair_score > best:
                    best = pair_score
                    M_pointers[i, j] = ('P', i+1, j-1)

                for k in range(i, j):
                    split_score = M[i, k] + M[k+1, j]
                    if split_score > best:
                        best = split_score
                        M_pointers[i, j] = ('S', i, k, k+1, j)

                M[i, j] = best

        self.M = M
        self.M_pointers = M_pointers
    
    def backtrace(self):
        pairs = []
        self.backtrace_rec(0, self.n - 1, pairs)
        self.pairs = pairs

    def backtrace_rec(self, i, j, pairs):
        if i >= j:
            return

        ptr = self.M_pointers[i, j]

        if ptr[0] == 'Ui':
            self.backtrace_rec(ptr[1], ptr[2], pairs)
        elif ptr[0] == 'Uj':
            self.backtrace_rec(ptr[1], ptr[2], pairs)
        elif ptr[0] == 'P':
            pairs.append((i, j))
            self.backtrace_rec(ptr[1], ptr[2], pairs)
        elif ptr[0] == 'S':
            self.backtrace_rec(ptr[1], ptr[2], pairs)
            self.backtrace_rec(ptr[3], ptr[4], pairs)