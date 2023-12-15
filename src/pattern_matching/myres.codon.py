from bio import *
from collections import defaultdict


def myers(P, T, k):
    B = defaultdict(int)
    i = 1

    for c in P:
        B[c] = B[c] | i
        i = i << 1

    m = len(P)
    VP = (1 << m) - 1  # Set m bits
    VN = 0
    score = m

    # Main loop
    for pos, c in enumerate(T):
        Eq = B[c]
        Xv = Eq | VN
        Xh = (((Eq & VP) + VP) ^ VP) | Eq
        HP = VN | ~(Xh | VP)
        HN = VP & Xh

        if HP & (1 << (m - 1)):
            score += 1

        elif HN & (1 << (m - 1)):
            score -= 1

        if score <= k:
            yield pos

        HP <<= 1
        VP = (HN << 1) | ~(Xv | HP)
        VN = HP & Xv


def test_myres():
    text = "CGGTCCTGAGGGATTAGCAC"
    pattern = "TCCTAGGGC"
    k = 2
    matches = list(myers(pattern, text, k))

    assert [(match, k) for match in matches] == [(11, 2), (12, 2)]

    print("All Tests Passed")


if __name__ == "__main__":
    test_myres()
