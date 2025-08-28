import random

def new_random_color(colors: list[str]) -> str:
    r = lambda: random.randint(0, 255)
    color = f"#{r():02X}{r():02X}{r():02X}"

    while color in colors:
        color = f"#{r():02X}{r():02X}{r():02X}"

    return color