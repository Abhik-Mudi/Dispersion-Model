import matplotlib.pyplot as plt
import numpy as np


# ==============================================================================
# 1. DISPERSION COEFFICIENTS (PLUME & PUFF CORRELATIONS)
# ==============================================================================

def get_plume_sigmas(
    x_m: np.ndarray | float, stability_class: str = "D", terrain: str = "rural"
) -> tuple[np.ndarray, np.ndarray]:
    """Calculates continuous plume dispersion coefficients (sigma_y, sigma_z) in meters."""
    x = np.maximum(np.asarray(x_m, dtype=float), 0.1)
    sc = stability_class.upper()
    terr = terrain.lower()

    if terr == "rural":
        if sc == "A":
            sy = 0.22 * x * (1.0 + 0.0001 * x) ** (-0.5)
            sz = 0.20 * x
        elif sc == "B":
            sy = 0.16 * x * (1.0 + 0.0001 * x) ** (-0.5)
            sz = 0.12 * x
        elif sc == "C":
            sy = 0.11 * x * (1.0 + 0.0001 * x) ** (-0.5)
            sz = 0.08 * x * (1.0 + 0.0002 * x) ** (-0.5)
        elif sc == "D":
            sy = 0.08 * x * (1.0 + 0.0001 * x) ** (-0.5)
            sz = 0.06 * x * (1.0 + 0.0015 * x) ** (-0.5)
        elif sc == "E":
            sy = 0.06 * x * (1.0 + 0.0001 * x) ** (-0.5)
            sz = 0.03 * x * (1.0 + 0.0003 * x) ** (-1.0)
        elif sc == "F":
            sy = 0.04 * x * (1.0 + 0.0001 * x) ** (-0.5)
            sz = 0.016 * x * (1.0 + 0.0003 * x) ** (-1.0)
        else:
            raise ValueError(f"Invalid stability class '{sc}'. Choose from A to F.")

    elif terr == "urban":
        if sc in ("A", "B"):
            sy = 0.32 * x * (1.0 + 0.0004 * x) ** (-0.5)
            sz = 0.24 * x * (1.0 + 0.001 * x) ** 0.5
        elif sc == "C":
            sy = 0.22 * x * (1.0 + 0.0004 * x) ** (-0.5)
            sz = 0.20 * x
        elif sc == "D":
            sy = 0.16 * x * (1.0 + 0.0004 * x) ** (-0.5)
            sz = 0.14 * x * (1.0 + 0.0003 * x) ** (-0.5)
        elif sc in ("E", "F"):
            sy = 0.11 * x * (1.0 + 0.0004 * x) ** (-0.5)
            sz = 0.08 * x * (1.0 + 0.0015 * x) ** (-0.5)
        else:
            raise ValueError(f"Invalid stability class '{sc}'. Choose from A to F.")
    else:
        raise ValueError("Terrain must be 'rural' or 'urban'.")

    return sy, sz


def get_puff_sigmas(
    x_m: np.ndarray | float, stability_class: str = "D"
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Calculates instantaneous puff dispersion coefficients (sigma_x, sigma_y, sigma_z) in meters."""
    x = np.maximum(np.asarray(x_m, dtype=float), 0.1)
    sc = stability_class.upper()

    # Slade (1968) quasi-instantaneous puff correlations: sy = a * x^b, sz = c * x^d
    params = {
        "A": (0.18, 0.92, 0.60, 0.75),
        "B": (0.14, 0.92, 0.53, 0.73),
        "C": (0.10, 0.92, 0.34, 0.71),
        "D": (0.06, 0.92, 0.15, 0.70),
        "E": (0.04, 0.92, 0.10, 0.65),
        "F": (0.02, 0.89, 0.05, 0.61),
    }
    if sc not in params:
        raise ValueError(f"Invalid stability class '{sc}'. Choose from A to F.")

    ay, by, cz, dz = params[sc]
    sy = ay * (x**by)
    sz = cz * (x**dz)
    sx = sy.copy()  # Isotropic horizontal turbulence assumption (sigma_x = sigma_y)
    return sx, sy, sz


# ==============================================================================
# 1.5. PASQUILL-GIFFORD TABLE 5-1 STABILITY PREDICTION
# ==============================================================================

TABLE_5_1_PASQUILL = [
    # Row 0: < 2 m/s
    {"day_strong": "A", "day_moderate": "B", "day_slight": "B", "night_overcast": "F", "night_clear": "F"},
    # Row 1: 2 - 3 m/s
    {"day_strong": "B", "day_moderate": "B", "day_slight": "C", "night_overcast": "E", "night_clear": "F"},
    # Row 2: 3 - 4 m/s
    {"day_strong": "B", "day_moderate": "C", "day_slight": "C", "night_overcast": "D", "night_clear": "E"},
    # Row 3: 4 - 6 m/s
    {"day_strong": "C", "day_moderate": "D", "day_slight": "D", "night_overcast": "D", "night_clear": "D"},
    # Row 4: > 6 m/s
    {"day_strong": "C", "day_moderate": "D", "day_slight": "D", "night_overcast": "D", "night_clear": "D"},
]

def determine_stability_from_environment(wind_speed_m_s: float, condition_key: str) -> str:
    """
    Predicts Pasquill-Gifford atmospheric stability class from Table 5-1
    given surface wind speed u (m/s) and environmental solar/night condition.
    """
    u = float(wind_speed_m_s)
    if u < 2.0:
        row_idx = 0
    elif u < 3.0:
        row_idx = 1
    elif u < 4.0:
        row_idx = 2
    elif u <= 6.0:
        row_idx = 3
    else:
        row_idx = 4

    return TABLE_5_1_PASQUILL[row_idx].get(condition_key, "D")


# ==============================================================================
# 2. CONCENTRATION EVALUATION AT ANY (x, y, z, t)
# ==============================================================================

def compute_concentration_at_point_and_time(
    model_type: str,
    source_strength: float,
    u_m_s: float,
    Hr_m: float,
    x_m: np.ndarray | float,
    y_m: np.ndarray | float,
    z_m: np.ndarray | float,
    t_sec: float,
    stability_class: str = "D",
    terrain: str = "rural",
) -> np.ndarray | float:
    """
    Computes concentration C(x, y, z, t) in kg/m^3 at any spatial coordinate and instant of time.
    - For 'puff': Uses the 3D transient Gaussian puff equation centered at x_center = u * t.
    - For 'plume': Uses the 3D Gaussian plume equation coupled with an advection wavefront
      (concentration is 0 where the plume front x > u * t has not yet arrived).
    """
    X = np.asarray(x_m, dtype=float)
    Y = np.asarray(y_m, dtype=float)
    Z = np.asarray(z_m, dtype=float)

    if model_type.lower() == "puff":
        # Dispersion coefficients grow with the distance traveled by the puff center (x_c = u * t)
        # or local downwind coordinate X
        travel_dist = max(u_m_s * t_sec, 0.1)
        eval_dist = np.maximum(X, travel_dist)
        sx, sy, sz = get_puff_sigmas(eval_dist, stability_class)

        denom = ((2.0 * np.pi) ** 1.5) * sx * sy * sz
        exp_x = np.exp(-0.5 * (((X - u_m_s * t_sec) / sx) ** 2))
        exp_y = np.exp(-0.5 * ((Y / sy) ** 2))
        exp_z = np.exp(-0.5 * (((Z - Hr_m) / sz) ** 2)) + np.exp(
            -0.5 * (((Z + Hr_m) / sz) ** 2)
        )
        C_kg_m3 = (source_strength / denom) * exp_x * exp_y * exp_z

    elif model_type.lower() == "plume":
        sy, sz = get_plume_sigmas(X, stability_class, terrain)
        denom = 2.0 * np.pi * sy * sz * u_m_s
        exp_y = np.exp(-0.5 * ((Y / sy) ** 2))
        exp_z = np.exp(-0.5 * (((Z - Hr_m) / sz) ** 2)) + np.exp(
            -0.5 * (((Z + Hr_m) / sz) ** 2)
        )
        C_steady = (source_strength / denom) * exp_y * exp_z

        # Advection front check: Plume front reaches x_front = u * t at time t
        plume_front_x = u_m_s * t_sec
        C_kg_m3 = np.where((X > 0.0) & (X <= plume_front_x), C_steady, 0.0)
    else:
        raise ValueError("model_type must be 'plume' or 'puff'.")

    return C_kg_m3


# ==============================================================================
# 3. INTERACTIVE USER INPUT & PLOTTING PIPELINE
# ==============================================================================

def get_input_with_default(prompt_text: str, default_val: str) -> str:
    """Helper to read user input from terminal with a fallback default value."""
    user_val = input(f"{prompt_text} [Default: {default_val}]: ").strip()
    return user_val if user_val != "" else default_val


def run_interactive_dispersion_simulator():
    print("=" * 75)
    print("   INTERACTIVE SPATIAL-TEMPORAL CHEMICAL DISPERSION SIMULATOR")
    print("=" * 75)
    print("Press ENTER at any prompt to accept the default value in brackets.\n")

    # --- 1. Collect Chemical & Release Model Inputs ---
    chem_name = get_input_with_default("1. Enter Chemical Name", "Chlorine")
    mw_g_mol = float(get_input_with_default("2. Molecular Weight (g/mol)", "70.9"))

    print("\nSelect Dispersion Model Type:")
    print("  [1] Continuous Plume Release (Steady release rate Qm in kg/s)")
    print("  [2] Instantaneous Puff Release (Fixed total mass Qm* in kg)")
    model_choice = get_input_with_default("Choose 1 or 2", "2")

    if model_choice == "1":
        model_type = "plume"
        source_strength = float(
            get_input_with_default("3. Continuous Mass Release Rate Qm (kg/s)", "0.5")
        )
    else:
        model_type = "puff"
        source_strength = float(
            get_input_with_default("3. Total Instantaneous Released Mass Qm* (kg)", "50.0")
        )

    # --- 2. Collect Weather & Site Inputs ---
    u_wind = float(get_input_with_default("4. Ambient Wind Speed u (m/s)", "3.0"))

    print("\nAtmospheric Stability Determination Method:")
    print("  [1] Direct Input: Manually choose stability class (A, B, C, D, E, F)")
    print("  [2] Auto-Predict: Predict from Environmental Conditions (Table 5-1)")
    mode_choice = get_input_with_default("Choose method [1 or 2]", "1")

    if mode_choice == "2":
        print("\nSelect Incoming Solar Radiation / Night Sky Condition:")
        print("  [1] Daytime — Strong Insolation (Sunny midday, high solar heating)")
        print("  [2] Daytime — Moderate Insolation (Moderate sun, slight clouds)")
        print("  [3] Daytime — Slight Insolation (Low sun angle / overcast day)")
        print("  [4] Nighttime — Thin Overcast or > 4/8 low cloud cover")
        print("  [5] Nighttime — Clear Sky (<= 3/8 cloudiness)")
        c_choice = get_input_with_default("Choose condition [1-5]", "4")
        c_map = {
            "1": "day_strong",
            "2": "day_moderate",
            "3": "day_slight",
            "4": "night_overcast",
            "5": "night_clear",
        }
        stability = determine_stability_from_environment(u_wind, c_map.get(c_choice, "night_overcast"))
        print(f"--> Auto-Predicted Stability Class from Table 5-1: Class {stability}")
    else:
        stability = get_input_with_default(
            "5. Pasquill-Gifford Stability Class (A, B, C, D, E, F)", "D"
        ).upper()
    terrain = get_input_with_default("6. Terrain Type (rural / urban)", "rural").lower()
    hr_m = float(get_input_with_default("7. Effective Release Height Hr (m)", "1.5"))
    c_target_ppm = float(
        get_input_with_default("8. Hazard Threshold Limit (e.g., ERPG-2 in ppm)", "3.0")
    )

    # --- 3. Collect Target Point (x, y, z) and Instant of Time (t) ---
    print("\n--- Point & Time Query Parameters ---")
    x_query = float(get_input_with_default("9.  Target Downwind Distance x (m)", "300.0"))
    y_query = float(get_input_with_default("10. Target Crosswind Distance y (m)", "0.0"))
    z_query = float(get_input_with_default("11. Target Elevation Above Ground z (m)", "0.0"))
    t_query = float(get_input_with_default("12. Instant of Time after Release t (seconds)", "100.0"))

    # --- 4. Perform Point Calculation at (x, y, z, t) ---
    c_point_kg_m3 = float(
        compute_concentration_at_point_and_time(
            model_type=model_type,
            source_strength=source_strength,
            u_m_s=u_wind,
            Hr_m=hr_m,
            x_m=x_query,
            y_m=y_query,
            z_m=z_query,
            t_sec=t_query,
            stability_class=stability,
            terrain=terrain,
        )
    )
    c_point_mg_m3 = c_point_kg_m3 * 1e6
    c_point_ppm = (c_point_kg_m3 * 24.45 * 1e6) / mw_g_mol

    # Also compute time of peak arrival at x_query
    t_arrival = x_query / max(u_wind, 0.01)

    print("\n" + "=" * 75)
    print(f"   POINT CONCENTRATION REPORT AT t = {t_query:.1f} SECONDS")
    print("=" * 75)
    print(f"Chemical & Model         : {chem_name} ({model_type.upper()} Model)")
    print(f"Target Coordinates       : x = {x_query:.1f} m, y = {y_query:.1f} m, z = {z_query:.1f} m")
    print(f"Advection Center / Front : x_front = {u_wind * t_query:.1f} m (at t = {t_query:.1f} s)")
    print(f"Expected Arrival at x    : t_arrival = {t_arrival:.1f} s")
    print("-" * 75)
    print(f"Concentration (kg/m³)    : {c_point_kg_m3:.6e} kg/m³")
    print(f"Concentration (mg/m³)    : {c_point_mg_m3:.4f} mg/m³")
    print(f"Concentration (ppm)      : {c_point_ppm:.4f} ppm")
    status = "EXCEEDS THRESHOLD!" if c_point_ppm >= c_target_ppm else "Below Threshold (Safe)"
    print(f"Hazard Status vs {c_target_ppm} ppm : {status}")
    print("=" * 75)

    # --- 5. Build 2D Spatial Grid at Time t_query & Time-History Curve at (x_query, y_query, z_query) ---
    center_x = u_wind * t_query
    x_max_domain = max(x_query * 1.5, center_x * 1.6, 200.0)
    y_max_domain = max(abs(y_query) * 2.0, x_max_domain * 0.25, 40.0)

    x_grid = np.linspace(1.0, x_max_domain, 350)
    y_grid = np.linspace(-y_max_domain, y_max_domain, 350)
    X, Y = np.meshgrid(x_grid, y_grid)

    C_grid_kg_m3 = compute_concentration_at_point_and_time(
        model_type=model_type,
        source_strength=source_strength,
        u_m_s=u_wind,
        Hr_m=hr_m,
        x_m=X,
        y_m=Y,
        z_m=z_query,
        t_sec=t_query,
        stability_class=stability,
        terrain=terrain,
    )
    C_grid_ppm = (C_grid_kg_m3 * 24.45 * 1e6) / mw_g_mol

    # Compute Concentration vs. Time at the user's specific point (x_query, y_query, z_query)
    t_max_history = max(t_query * 2.0, t_arrival * 2.0, 60.0)
    t_vec = np.linspace(0.1, t_max_history, 400)
    c_time_series_ppm = [
        (
            float(
                compute_concentration_at_point_and_time(
                    model_type,
                    source_strength,
                    u_wind,
                    hr_m,
                    x_query,
                    y_query,
                    z_query,
                    t_val,
                    stability,
                    terrain,
                )
            )
            * 24.45
            * 1e6
        )
        / mw_g_mol
        for t_val in t_vec
    ]

    # --- 6. Plot 2D Concentration Contour Map & Temporal Evolution ---
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(11, 9), gridspec_kw={"height_ratios": [1.4, 1]}
    )

    # Top Plot: 2D Contour Map at Instant t = t_query
    peak_grid_ppm = float(np.max(C_grid_ppm))
    if peak_grid_ppm > 0:
        levels = np.linspace(0.0, max(peak_grid_ppm * 0.85, c_target_ppm * 2.0), 45)
        contour_fill = ax1.contourf(X, Y, C_grid_ppm, levels=levels, cmap="YlOrRd", extend="max")
        plt.colorbar(contour_fill, ax=ax1, label="Concentration (ppm)")

        if peak_grid_ppm >= c_target_ppm:
            iso = ax1.contour(
                X, Y, C_grid_ppm, levels=[c_target_ppm], colors="red", linewidths=2.2
            )
            ax1.clabel(
                iso,
                inline=True,
                fmt={c_target_ppm: f"Hazard Limit ({c_target_ppm} ppm)"},
                fontsize=9,
            )
    else:
        ax1.text(
            0.5,
            0.5,
            "Zero Concentration Across Domain at Selected Time",
            transform=ax1.transAxes,
            ha="center",
        )

    # Mark Release Source and Target Query Point
    ax1.plot(0, 0, "k*", markersize=13, label="Release Source (0, 0)")
    ax1.plot(
        x_query,
        y_query,
        marker="o",
        color="cyan",
        markeredgecolor="black",
        markersize=9,
        label=f"Query Point ({x_query}m, {y_query}m): {c_point_ppm:.2f} ppm",
    )

    ax1.set_title(
        f"2D Concentration Contour Map at Instant t = {t_query:.1f} s (z = {z_query} m)\n"
        f"Chemical: {chem_name} | Model: {model_type.capitalize()} | Wind: {u_wind} m/s | Class: {stability}",
        fontsize=11,
    )
    ax1.set_xlabel("Downwind Distance x (m)")
    ax1.set_ylabel("Crosswind Distance y (m)")
    ax1.grid(True, linestyle="--", alpha=0.5)
    ax1.legend(loc="upper right")

    # Bottom Plot: Concentration Breakthrough Curve at (x_query, y_query, z_query) vs Time
    ax2.plot(
        t_vec,
        c_time_series_ppm,
        color="darkblue",
        linewidth=2.2,
        label=f"C(t) at x={x_query}m, y={y_query}m, z={z_query}m",
    )
    ax2.axvline(
        t_query,
        color="teal",
        linestyle=":",
        linewidth=1.8,
        label=f"Evaluated Time t = {t_query} s ({c_point_ppm:.2f} ppm)",
    )
    ax2.axhline(
        c_target_ppm,
        color="red",
        linestyle="--",
        linewidth=1.5,
        label=f"Hazard Threshold ({c_target_ppm} ppm)",
    )
    ax2.plot(t_query, c_point_ppm, "ro", markersize=7)

    ax2.set_title(
        f"Temporal Concentration History at Target Location (x={x_query} m, y={y_query} m, z={z_query} m)",
        fontsize=11,
    )
    ax2.set_xlabel("Time Elapsed Since Release t (seconds)")
    ax2.set_ylabel("Concentration (ppm)")
    ax2.set_xlim(0, t_max_history)
    ax2.grid(True, linestyle=":", alpha=0.6)
    ax2.legend(loc="upper right")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    run_interactive_dispersion_simulator()



#       Problem 1
# ======================
# An industrial chemical plant releases sulfur dioxide (SO₂) continuously into the atmosphere through a chimney at the rate of 80 g/s. The release occurs in a rural area on an overcast day under neutral atmospheric stability conditions.
# The effective height of the stack is 60 m, and sulfur dioxide is released at a constant rate of 80 g/s. The wind is blowing at a constant speed of 6 m/s.
# Assume that the atmospheric conditions remain constant throughout the simulation.
# Determine the following:
# 1. Calculate the SO₂ concentration at a ground-level monitoring station located 500 m downwind and directly in line with the wind direction at 100 seconds after the simulation begins.
# 2. Calculate the SO₂ concentration at another ground-level monitoring station located 500 m downwind and 50 m crosswind at the same time.
# 3. Determine how the pollutant concentration varies with time at the first monitoring station.
# 4. Visualize the spatial distribution of SO₂ concentration at 100 seconds and identify the region with the highest concentration.
# 5. Determine whether the concentration at either monitoring station exceeds a specified threshold of 0.5 ppm.



#       Problem 2
# ======================
# An industrial chemical storage facility experiences an accidental rupture in a chlorine gas storage tank. As a result, 100 kg of chlorine gas is released instantaneously into the atmosphere.

# The release occurs in a rural area on an overcast day. The effective release height is 30 m, and the wind is blowing at a constant speed of 5 m/s.

# Assume that the released gas forms a Gaussian puff and is transported by the wind without any additional chemical reactions or deposition.

# Determine:
# 1. Calculate the chlorine concentration at a ground-level monitoring station located 500 m downwind and directly along the wind direction after 100 seconds.
# 2. Calculate the concentration at another ground-level monitoring station located 500 m downwind and 50 m crosswind after 100 seconds.
# 3. Determine how the concentration changes with time at the first monitoring station.
# 4. Visualize the spatial distribution of chlorine gas at 100 seconds.
# 5. Identify the approximate time at which the puff centre reaches the first monitoring station.