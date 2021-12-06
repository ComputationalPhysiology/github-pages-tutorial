ALPHABET = "abcdefghijklmnopqrstuvwxyz"


def rotate_string(msg, shift):
    return msg[shift:] + msg[:shift]


def create_shifted_alphabet(shift: int):
    rotated_alphabet = rotate_string(ALPHABET, shift)
    return dict(zip(ALPHABET, rotated_alphabet))


def encrypt(msg: str, shift: int) -> str:

    if not isinstance(msg, str):
        raise TypeError(f"We can only encrypt string, got {type(msg)}")

    shifted_alphabet = create_shifted_alphabet(shift)

    new_message = ""
    for letter in msg:
        new_message += shifted_alphabet[letter]

    return new_message
