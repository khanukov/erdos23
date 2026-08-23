"""Tiny exact-arithmetic fallback for the subset of :mod:`gmpy2` we use.

The certificate replay prefers gmpy2 when it is installed.  Python's
``fractions.Fraction`` is slower but mathematically equivalent for these
operations and keeps the independent verifier runnable in a minimal
environment.
"""

from fractions import Fraction


def mpz(value=0):
    return int(value)


def mpq(numerator=0, denominator=None):
    if denominator is None:
        if isinstance(numerator, tuple):
            return Fraction(*numerator)
        return Fraction(numerator)
    return Fraction(numerator, denominator)


def numer(value):
    return Fraction(value).numerator


def denom(value):
    return Fraction(value).denominator
