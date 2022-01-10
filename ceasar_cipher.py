ALPHABET = "abcdefghijklmnopqrstuvwxyz"


def rotate_string(msg: str, shift: int) -> str:
    """Rotate string with a given shift

    Parameters
    ----------
    msg : str
        The string you want to roate
    shift : int
        The shift

    Returns
    -------
    str
        The rotated string

    Example
    -------

    .. code::

        >> rotated_string("hello", 1)
        "ohell"

    """
    return msg[shift:] + msg[:shift]


def create_shifted_alphabet(shift: int):
    """Return a dictionary where the keys are
    the original letters in the alphabet and the
    values are the corresponding letters in the
    shifted alphabet

    Parameters
    ----------
    shift : int
        The shift

    Returns
    -------
    dict
        The dictioary with the shifted alphabet

    Example
    -------

    .. code::

        >> shifted_alphabet = create_shifted_alphabet(2)
        >> shifted_alphabet["d"]
        "f"

    """
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
