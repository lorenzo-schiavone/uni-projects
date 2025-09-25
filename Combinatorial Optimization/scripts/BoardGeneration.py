import numpy as np
import random
import matplotlib.pyplot as plt
import json

class sparseBoard:
    def __init__(self, pairs=[]):
        self.pairs = list(pairs)  # list of tuples

    def safe_add(self, pair):
        if pair not in self.pairs:
            self.pairs.append(pair)

    def add_random(self, num, maxrow, maxcol):
        for _ in range(num):
            x = random.randint(0, maxrow)
            y = random.randint(0, maxcol)
            self.safe_add((x, y))

    def add_vline(self, start, length, step):
        x, y = start
        max_x = x + length
        if max_x > start[0] + length or max_x < 0:
            return
        for i in range(0, length + 1, step):
            self.safe_add((x + i, y))

    def add_hline(self, start, length, step):
        x, y = start
        max_y = y + length
        if max_y > start[1] + length or max_y < 0:
            return
        for j in range(0, length + 1, step):
            self.safe_add((x, y + j))

    def add_rect(self, start, width, height, step):
        x, y = start
        self.add_hline((x, y), width, step)
        self.add_hline((x + height, y), width, step)
        self.add_vline((x, y), height, step)
        self.add_vline((x, y + width), height, step)

    def add_diag(self, start, length, step):
        x, y = start
        for i in range(0, length + 1, step):
            self.safe_add((x + i, y + i))

    def add_triangle(self, start, height, step):
        x, y = start
        self.add_diag(start, height, step)
        self.add_vline(start, height, step)
        self.add_hline((x + height, y), height, step)

    def nnz(self):
        return len(self.pairs)

    def make_pos_dict(self):
        return dict(enumerate(self.pairs))

    def make_cost(self):
        num_nonzeros = self.nnz()
        cost_matrix = np.zeros((num_nonzeros, num_nonzeros))
        exp = 2
        for i in range(num_nonzeros):
            x1, y1 = self.pairs[i]
            for j in range(i, num_nonzeros):
                x2, y2 = self.pairs[j]
                distance = (abs(x1 - x2)**exp + abs(y1 - y2)**exp)**(1/exp)
                cost_matrix[i][j] = distance
                cost_matrix[j][i] = distance
        return np.round(cost_matrix, 4)
    
    def random_instance(self, num_nonzeros, max_row, max_col):
        self.pairs = []
        shapes = [self.add_vline, self.add_hline, self.add_rect, self.add_diag, self.add_triangle, self.add_random]
        current_nonzeros = 0
        lenLimit= 3
        while current_nonzeros < num_nonzeros:
            shape = random.choice(shapes)
            step = 1
            if shape == self.add_vline:
                start_x = random.randint(0, max_row)
                start_y = random.randint(0, max_col)
                max_length = max_row - start_x
                if max_length < 1:
                    continue
                length = random.randint(1, max_length)
                length = min(length, lenLimit)
                self.add_vline((start_x, start_y), length, step)
            elif shape == self.add_hline:
                start_x = random.randint(0, max_row)
                start_y = random.randint(0, max_col - 1)
                max_length = max_col - start_y
                if max_length < 1:
                    continue
                length = random.randint(1, max_length)
                length = min(length, lenLimit)
                self.add_hline((start_x, start_y), length, step)
            elif shape == self.add_rect:
                start_x = random.randint(0, max_row - 2)
                start_y = random.randint(0, max_col - 2)
                max_width = max_col - start_y
                max_height = max_row - start_x
                if max_width < 2 or max_height < 2:
                    continue
                width = random.randint(2, max_width)
                height = random.randint(2, max_height)
                width = min(width, lenLimit)
                height = min(height, lenLimit)
                self.add_rect((start_x, start_y), width, height, step)
            elif shape == self.add_diag:
                start_x = random.randint(0, max_row - 1)
                start_y = random.randint(0, max_col - 1)
                max_length = min(max_row - start_x, max_col - start_y)
                if max_length < 1:
                    continue
                length = random.randint(1, max_length)
                length = min(length, lenLimit)
                self.add_diag((start_x, start_y), length, step)
            elif shape == self.add_triangle:
                start_x = random.randint(0, max_row - 2)
                start_y = random.randint(0, max_col - 2)
                max_height_row = max_row - start_x
                max_height_col = max_col - start_y
                max_height = min(max_height_row, max_height_col)
                if max_height < 2:
                    continue
                height = random.randint(2, max_height)
                height = min(height, lenLimit)
                self.add_triangle((start_x, start_y), height, step)
            elif shape == self.add_random:
                remaining = num_nonzeros - self.nnz()
                if remaining <= 0:
                    continue
                num = min(remaining, 2) #random.randint(1, max(1, remaining // 5))
                self.add_random(num, max_row, max_col)
            current_nonzeros = self.nnz()

        if self.nnz() > num_nonzeros:
            self.pairs = self.pairs[:num_nonzeros]

nnzs = [16, 32, 64, 128, 156]
K = 2 # how many per size
width = 40
height = 20

for nnz in nnzs:
    for k in range(K):
        board = sparseBoard()
        board.random_instance(nnz, width, height)
        np.savetxt(f'../data/cost_matrix/cost{nnz}_{k}.dat', board.make_cost(), fmt='%1.4f', header=f'{nnz}', comments='')
        pos = board.make_pos_dict()
        with open(f'../data/pos_dict/pos_dict{nnz}_{k}.json', 'w') as fp:
            json.dump(pos, fp)
