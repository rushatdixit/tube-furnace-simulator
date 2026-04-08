"""
Defines the material properties using dataclasses.
The thermal conductivity may be evaluated as a function of temperature.
"""

from dataclasses import dataclass
from typing import Callable
from src.physics.constants import SimulationConstants
from src.physics.models import Kelvin


@dataclass(frozen=True)
class Material:
    """
    Base class for a physical material.
    
    Fields:
    - name: Name of the material
    - density: kg/m^3
    - specific_heat: J/(kg*K)
    - thermal_conductivity: Function mapping temperature (K) to W/(m*K)
    """
    name : str
    density : float
    specific_heat : float
    thermal_conductivity : Callable[[float], float]

    def get_thermal_conductivity(self, temperature : Kelvin) -> float:
        """
        Evaluates and returns the thermal conductivity at the given temperature.
        """
        return self.thermal_conductivity(temperature)

@dataclass(frozen=True)
class Insulator:
    """
    A class for an insulating material.
    
    Fields:
    - material: The base Material object
    - reference_resistivity: Ohms*m
    - temp_coeff: Temperature coefficient
    - reference_temperature: Reference temperature for the coefficient
    """
    material : Material
    reference_resistivity : float
    temp_coeff : float
    reference_temperature : Kelvin = SimulationConstants.AMBIENT_TEMPERATURE

    def get_conductivity(self, temperature : Kelvin) -> float:
        """
        Calculates and returns the conductivity based on the material's temperature coefficient.
        """
        rho = self.reference_resistivity*(1 + self.temp_coeff*(temperature.temperature - self.reference_temperature))
        return 1/rho
