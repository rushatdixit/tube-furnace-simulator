"""
A library of predefined concrete materials and properties.
"""

from src.physics.materials import Material, Conductor

ALUMINA = Material(
    name="Alumina_995",
    density=3950.0,
    specific_heat=880.0,
    k_func=lambda T: 35.0 * (293.15 / T)**0.5
)

KANTHAL = Conductor(
    name="Kanthal_A1",
    density=7100.0,
    specific_heat=460.0,
    k_func=lambda T: 13.0 + (0.005 * T),
    resistivity_ref=1.45e-6,
    temp_coeff=5.0e-5
)

INSULATION = Material(
    name="Ceramic_Fiber_Board",
    density=128.0,
    specific_heat=1130.0,
    k_func=lambda T: 0.05 + 1.5e-4 * (T - 273.15)
)