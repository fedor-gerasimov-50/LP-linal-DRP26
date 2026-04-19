"""Expected number of coin flips before a target pattern appears.

A solution 
This module `expected_time_for_pattern`, which computes the expected
number of flips until a target pattern appears, using a finite Markov chain.
"""

from __future__ import annotations

from fractions import Fraction
from typing import List, Sequence, Union

Number = Union[Fraction, float, int]


def _prefix_function(s: str) -> List[int]:
    """Compute KMP prefix-function values for `s`."""
    pi = [0] * len(s)
    for i in range(1, len(s)):
        j = pi[i - 1]
        while j > 0 and s[i] != s[j]:
            j = pi[j - 1]
        if s[i] == s[j]:
            j += 1
        pi[i] = j
    return pi


def _build_transitions(pattern: str) -> List[List[int]]:
    """Build automaton transitions for states 0..m (m is absorbing)."""
    m = len(pattern)
    pi = _prefix_function(pattern)
    transitions = [[0, 0] for _ in range(m + 1)]

    # Absorbing state once the full pattern has been observed.
    transitions[m][0] = m
    transitions[m][1] = m

    for state in range(m):
        for bit_char in ("0", "1"):
            j = state
            while j > 0 and pattern[j] != bit_char:
                j = pi[j - 1]
            if pattern[j] == bit_char:
                j += 1
            transitions[state][int(bit_char)] = j

    return transitions


def _gaussian_elimination(a: List[List[Number]], b: List[Number]) -> List[Number]:
    """Solve Ax=b with Gaussian elimination (supports Fraction/float)."""
    n = len(a)

    for col in range(n):
        pivot = col
        while pivot < n and a[pivot][col] == 0:
            pivot += 1
        if pivot == n:
            raise ValueError("Singular system while solving Markov expectation.")

        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            b[col], b[pivot] = b[pivot], b[col]

        pivot_val = a[col][col]
        for j in range(col, n):
            a[col][j] = a[col][j] / pivot_val
        b[col] = b[col] / pivot_val

        for row in range(n):
            if row == col:
                continue
            factor = a[row][col]
            if factor == 0:
                continue
            for j in range(col, n):
                a[row][j] = a[row][j] - factor * a[col][j]
            b[row] = b[row] - factor * b[col]

    return b


def expected_time_for_pattern(p: Number, pattern: str) -> Number:
    """Return expected flips until `pattern` appears in Bernoulli(p) trials.

    Args:
        p: Probability of observing '1' on each flip. Can be float, int,
            or Fraction. If provided as Fraction (or int), arithmetic remains
            exact and the result is a Fraction.
        pattern: Non-empty string containing only '0' and '1'.

    Returns:
        Expected number of flips until `pattern` first appears.
    """
    if not pattern or any(ch not in "01" for ch in pattern):
        raise ValueError("pattern must be a non-empty string of '0' and '1'.")

    if isinstance(p, Fraction):
        p_num: Number = p
    elif isinstance(p, int):
        p_num = Fraction(p, 1)
    else:
        p_num = float(p)

    if not (0 <= p_num <= 1):
        raise ValueError("p must be in [0, 1].")

    ones = pattern.count("1")
    zeros = len(pattern) - ones
    if p_num == 0 and ones > 0:
        return float("inf")
    if p_num == 1 and zeros > 0:
        return float("inf")

    q_num = 1 - p_num
    m = len(pattern)
    transitions = _build_transitions(pattern)

    # Unknowns are E[0], ..., E[m-1]; E[m]=0 is absorbing.
    a: List[List[Number]] = [[0 for _ in range(m)] for _ in range(m)]
    b: List[Number] = [1 for _ in range(m)]

    for state in range(m):
        a[state][state] = 1

        next0 = transitions[state][0]
        next1 = transitions[state][1]

        if next0 < m:
            a[state][next0] = a[state][next0] - q_num
        if next1 < m:
            a[state][next1] = a[state][next1] - p_num

    solution = _gaussian_elimination(a, b)
    return solution[0]


if __name__ == "__main__":
    # Fair coin: expected time for HTH (encoded 101) is 10.
    print(expected_time_for_pattern(Fraction(1, 2), "1010101"))
