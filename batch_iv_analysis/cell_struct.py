from dataclasses import dataclass, field
from typing import List


@dataclass
class cell:
    i: List[float] = field(default_factory=lambda: [0])  # [A] current
    v: List[float] = field(default_factory=lambda: [0])  # [V] voltage
    t: List[float] = field(default_factory=lambda: [0])  # [s] time
    s: List[int] = field(default_factory=lambda: [0])  # status byte
    area: float = 1 * 1e-4  # [m^2] device area
    intensity: float = 1  # [suns]
