from bio import *


class Ukkonen:
    D: List[List[int]]

    def __init__(self, m):
        self.D = [[0] * (m + 1) for _ in range(2)]

    def unit_cost(self, a, b):
        return int(a != b)

    def find_all_end(self, pattern, text, k):
        m = len(pattern)
        self.D[0] = [k + 1] * (m + 1)
        self.D[1] = list(range(m + 1))
        lastk = min(k, m)

        for i, text_char in enumerate(text):
            col, prev = i % 2, 1 - i % 2
            self.D[col][0] = 0
            lastk = min(lastk + 1, m)

            for j in range(1, lastk + 1):
                pattern_char = pattern[j - 1]
                cost = self.unit_cost(pattern_char, text_char)
                self.D[col][j] = min(
                    self.D[prev][j] + 1,
                    self.D[col][j - 1] + 1,
                    self.D[prev][j - 1] + cost,
                )

            while self.D[col][lastk] > k:
                lastk -= 1

            if lastk == m:
                yield (i, self.D[col][m])


def test_ukkonen():
    ukkonen = Ukkonen(10)
    text = "ACCGTGGATGAGCGCCATAG"
    pattern = "TAGCGC"
    occ = list(ukkonen.find_all_end(pattern, text, 1))

    assert occ == [(14, 1)]

    print("All Tests Passed")


if __name__ == "__main__":
    test_ukkonen()
