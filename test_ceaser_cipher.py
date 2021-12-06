import pytest
import ceasar_cipher


def test_rotate_string():
    assert ceasar_cipher.rotate_string("hello", 1) == "elloh"


def test_create_shifted_alphabet():
    shifted_alphabet = ceasar_cipher.create_shifted_alphabet(shift=1)
    assert shifted_alphabet["a"] == "b"
    assert shifted_alphabet["b"] == "c"
    assert shifted_alphabet["z"] == "a"


@pytest.mark.parametrize(
    "msg, shift, expected", [("hello", 1, "ifmmp"), ("welcome", 3, "zhofrph")]
)
def test_encrypt(msg, shift, expected):
    assert ceasar_cipher.encrypt(msg, shift=shift) == expected


def test_encrypt_integer_raises_TypeError():
    with pytest.raises(TypeError):
        ceasar_cipher.encrypt(42, 1)
