from laplace import laplace
import pytest


@pytest.fixture(scope="session")
def solution():
    return laplace()


def test_laplace_left(solution):
    assert abs(solution(0.0, 0.5, 0.5)) < 1e-12


def test_laplace_right(solution):
    assert abs(solution(1.0, 0.5, 0.5) - 1) < 1e-12
