import random

def new_random_color(colors: list[str]) -> str:
    """
    Create a new random hex color string. The output is different from all that are in the input list.
    :param colors: A list of colors to avoid
    :return: A new hex color string
    """
    r = lambda: random.randint(0, 255)
    color = f"#{r():02X}{r():02X}{r():02X}"

    while color in colors:
        color = f"#{r():02X}{r():02X}{r():02X}"

    return color