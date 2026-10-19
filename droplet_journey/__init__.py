"""droplet_journey: 0-D single-droplet evaporation and heating with events."""

from .drag import drag_coefficient, reynolds, terminal_velocity
from .model import (
    d2_law_constant,
    d2_law_diameter,
    droplet_rhs,
    extinction_time,
    integrate_droplet,
    lumped_heating_rhs,
    saturation_pressure_clapeyron,
    saturation_vapor_mass_fraction,
    spalding_mass_number,
)

__all__ = [
    "drag_coefficient",
    "reynolds",
    "terminal_velocity",
    "d2_law_constant",
    "d2_law_diameter",
    "droplet_rhs",
    "extinction_time",
    "integrate_droplet",
    "lumped_heating_rhs",
    "saturation_pressure_clapeyron",
    "saturation_vapor_mass_fraction",
    "spalding_mass_number",
]
__version__ = "1.0.0"
