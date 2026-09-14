"""
7DT instrument and filter constants.
"""
from pathlib import Path

PACKAGE_DIR = Path(__file__).resolve().parent
PIXEL_SCALE = 0.5   # arcsec/pixel

FILTER_WAVE: dict[str, tuple[int, int]] = {
    "g":    (475, 130),
    "r":    (620, 130),
    "i":    (760, 130),
    "z":    (900, 100),
    "m400": (400, 25), "m425": (425, 25), "m450": (450, 25), "m475": (475, 25),
    "m500": (500, 25), "m525": (525, 25), "m550": (550, 25), "m575": (575, 25),
    "m600": (600, 25), "m625": (625, 25), "m650": (650, 25), "m675": (675, 25),
    "m700": (700, 25), "m725": (725, 25), "m750": (750, 25), "m775": (775, 25),
    "m800": (800, 25), "m825": (825, 25), "m850": (850, 25), "m875": (875, 25),
}

BROAD_BAND_COLORS: dict[str, str] = {
    "g": "green",
    "r": "red",
    "i": "goldenrod",
    "z": "black",
}
