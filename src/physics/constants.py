"""
Simulation configuration and physical constants.
"""

from dataclasses import dataclass
from src.physics.models import Kelvin, Celsius

@dataclass
class SimulationConstants:
    """
    Constants for the simulation.
    
    Fields:
    - STEPHANS_CONSTANT: W/(m^2 K^4)
    - ABSOLUTE_ZERO_C: Absolute zero in Celsius
    - AMBIENT_TEMPERATURE: Ambient temperature in Kelvin
    - P_ATM: Standard atmospheric pressure in Pa
    - CONVERGENCE_TOL: Tolerance for solver convergence
    - MAX_ITERATIONS: Maximum solver iterations
    """
    STEPHANS_CONSTANT : float = 5.670374419e-8
    ABSOLUTE_ZERO_C: Celsius = -273.15
    AMBIENT_TEMPERATURE: Kelvin = 298.15
    P_ATM: float = 101325.0

    CONVERGENCE_TOL: float = 1e-6
    MAX_ITERATIONS: int = 100
