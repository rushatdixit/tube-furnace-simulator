import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def get_oxide_color(thickness_nm):
    """Returns the expected physical color and hex code based on standard thin-film interference."""
    color_chart = [
        (15, "Silver/Grey (Bare Si)", "#B0C4DE"),
        (30, "Faint Tan", "#D2B48C"),
        (55, "Tan", "#DEB887"),
        (75, "Brown", "#8B4513"),
        (100, "Dark Violet", "#8A2BE2"),
        (130, "Royal Blue", "#4169E1"), # 100-130nm range is deeply blue!
        (150, "Light Blue", "#87CEEB"),
        (175, "Yellow-Green", "#9ACD32"),
        (200, "Light Gold", "#FAFAD2"),
        (225, "Gold", "#FFD700"),
        (250, "Orange", "#FFA500"),
        (275, "Red-Orange", "#FF4500"),
        (300, "Red-Violet", "#C71585"),
        (330, "Blue", "#0000FF"),
        (370, "Blue-Green", "#008080"),
        (400, "Green", "#00FF00"),
        (10000, "Opaque/Iridescent", "#FFFFFF")
    ]
    for max_thick, name, hex_code in color_chart:
        if thickness_nm <= max_thick:
            return name, hex_code
    return "Unknown", "#FFFFFF"

def calculate_water_partial_pressure(T_room, RH_percent):
    """Calculates the partial pressure of H2O in atm using the Magnus-Tetens formula."""
    # Saturation vapor pressure of water in kPa
    P_sat_kPa = 0.61078 * np.exp((17.27 * T_room) / (T_room + 237.3))
    # Actual vapor pressure based on Relative Humidity
    P_actual_kPa = P_sat_kPa * (RH_percent / 100.0)
    # Convert kPa to Atmospheres (1 atm = 101.325 kPa)
    P_H2O_atm = P_actual_kPa / 101.325
    return P_H2O_atm

def simulate_mixed_deal_grove(T_furnace, t_hours, T_room=25.0, RH=0.0):
    k = 8.617e-5  
    T_kelvin = T_furnace + 273.15
    
    # Partial Pressures
    P_H2O = calculate_water_partial_pressure(T_room, RH)
    P_O2 = 0.21 - (P_H2O / 2.0) # Standard air is ~21% O2, slightly displaced by humidity
    
    # 1. DRY CONSTANTS (O2)
    C_B_dry, E_B_dry = 772.0, 1.23         
    C_BA_dry, E_BA_dry = 3.71e6, 2.00        
    B_dry = C_B_dry * np.exp(-E_B_dry / (k * T_kelvin)) * P_O2
    BA_dry = C_BA_dry * np.exp(-E_BA_dry / (k * T_kelvin)) * P_O2
    
    # 2. WET CONSTANTS (H2O)
    C_B_wet, E_B_wet = 386.0, 0.71         
    C_BA_wet, E_BA_wet = 9.70e7, 2.05        
    B_wet = C_B_wet * np.exp(-E_B_wet / (k * T_kelvin)) * P_H2O
    BA_wet = C_BA_wet * np.exp(-E_BA_wet / (k * T_kelvin)) * P_H2O
    
    # 3. MIXED AMBIENT SUPERPOSITION
    B_mixed = B_dry + B_wet
    BA_mixed = BA_dry + BA_wet
    A_mixed = B_mixed / BA_mixed if BA_mixed != 0 else 1e-10
    
    x_initial = 0.02 # 20nm native oxide
    tau = (x_initial**2 + A_mixed * x_initial) / B_mixed
    
    t_array = np.linspace(0, t_hours, 1000)
    thickness_um = (-A_mixed + np.sqrt(A_mixed**2 + 4 * B_mixed * (t_array + tau))) / 2.0
    
    return t_array, thickness_um, P_H2O

def plot_oxidation():
    print("\n--- Silicon Oxidation Simulator (Mumbai Humidity Edition) ---")
    
    try:
        T_core = float(input("Enter furnace holding temperature [°C] (e.g., 1180): ") or 1180.0)
        T_room = float(input("Enter lab room temperature [°C] (e.g., 25): ") or 25.0)
        RH = float(input("Enter lab relative humidity [%] (e.g., 75): ") or 75.0)
        target_min = float(input("Enter target time [minutes] to calculate exact thickness: ") or 60.0)
    except ValueError:
        print("Invalid input. Using defaults (1180°C, 25°C Lab, 75% RH, 60 mins).")
        T_core, T_room, RH, target_min = 1180.0, 25.0, 75.0, 60.0

    target_hr = target_min / 60.0
    
    # Run the simulation
    _, x_mixed_exact, p_h2o = simulate_mixed_deal_grove(T_core, target_hr, T_room, RH)
    
    thick_nm = x_mixed_exact[-1] * 1000
    color_name, _ = get_oxide_color(thick_nm)
    
    print("\n=======================================================")
    print(f"[EXACT CALCULATION] At {target_min} minutes @ {T_core}°C")
    print(f" -> Lab Conditions: {T_room}°C with {RH}% Humidity")
    print(f" -> Water Vapor Partial Pressure: {p_h2o*100:.2f}% of atmosphere")
    print(f" -> Calculated Thickness: {thick_nm:.2f} nm")
    print(f" -> Expected Physical Color: {color_name}")
    print("=======================================================\n")
    print("Generating High-Resolution Graph...")

    out_dir = Path("data/graphs")
    out_dir.mkdir(parents=True, exist_ok=True)
    plt.style.use('dark_background')

    # GRAPH 1: SHORT-TERM GROWTH (MINUTES)
    hours_short = 3.0 
    t_mixed_s, x_mixed_s, _ = simulate_mixed_deal_grove(T_core, hours_short, T_room, RH)
    
    fig1, ax1 = plt.subplots(figsize=(12, 7))
    t_mins = t_mixed_s * 60 
    
    ax1.plot(t_mins, x_mixed_s * 1000, color='lime', linewidth=2.5, label=f"Humid Lab Air ({RH}% RH)")
    
    # Label every 20 minutes with actual colors
    checkpoints_min = np.arange(20, int(hours_short * 60) + 1, 20)
    for m in checkpoints_min:
        idx = np.argmin(np.abs(t_mins - m))
        th_nm = x_mixed_s[idx] * 1000
        c_name, c_hex = get_oxide_color(th_nm)
        
        ax1.plot(m, th_nm, marker='o', markersize=10, color=c_hex, markeredgecolor='white', markeredgewidth=1.5)
        ax1.text(m, th_nm + 15, f"{th_nm:.1f}nm\n({c_name})", color='white', fontsize=9, rotation=45, ha='left')

    ax1.set_title(f"SiO$_2$ Growth in Humid Lab Air ({T_room}°C, {RH}% RH) furnace at {T_core}°C", fontweight='bold', fontsize=14)
    ax1.set_xlabel("Time (Minutes)", fontsize=12)
    ax1.set_ylabel("Oxide Thickness (Nanometers)", fontsize=12)
    ax1.grid(color='gray', linestyle='--', alpha=0.3)
    ax1.legend(loc='upper left', fontsize=12)
    
    out_min = out_dir / "oxidation_humid_minutes.png"
    fig1.savefig(out_min, dpi=300, bbox_inches='tight')
    plt.close(fig1)

    print(f"Success! Graph saved to: {out_min}")

if __name__ == "__main__":
    plot_oxidation()