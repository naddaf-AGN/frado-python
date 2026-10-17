"""Interchangeable shielding geometry; bounds follow the MATLAB patch model."""
from dataclasses import dataclass
from .config import Config

@dataclass(frozen=True)
class PatchShielding:
    cfg: Config

    def bounds(self, rho: float, z: float):
        rmin = self.cfg.inner_radiating_radius
        patch = self.cfg.patch_size
        if rho-patch*z <= rmin:
            return rmin, rho+patch*z, 0.0, 3.141592653589793
        return rho-patch*z, rho+patch*z, 0.0, 3.141592653589793*patch*z/(rho-rmin)
