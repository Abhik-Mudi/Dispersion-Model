import numpy as np


# ==============================================================================
# SOURCE MODEL ALGORITHMS (MODELS 1 TO 7)
# ==============================================================================

class SourceModels:
    """Implementations for the 7 primary chemical release source models."""

    R_GAS = 8.314  # Universal gas constant (J / mol*K)
    G_ACCEL = 9.81  # Gravitational acceleration (m/s^2)

    @classmethod
    def liquid_hole_constant_head(
        cls,
        rho_kg_m3: float,
        hole_diam_m: float,
        P_gauge_Pa: float = 0.0,
        h_L_m: float = 0.0,
        Cd: float = 0.61,
        gc: float = 1.0,
    ) -> float:
        """Model 1: Flow of liquid through a hole with constant head/pressure."""
        area = np.pi * (hole_diam_m / 2.0) ** 2
        driving_head = (gc * P_gauge_Pa / rho_kg_m3) + (cls.G_ACCEL * h_L_m)
        if driving_head <= 0.0:
            return 0.0
        return float(Cd * area * rho_kg_m3 * np.sqrt(2.0 * driving_head))

    @classmethod
    def liquid_tank_draining(
        cls,
        rho_kg_m3: float,
        hole_diam_m: float,
        tank_diam_m: float,
        h_L0_m: float,
        P_gauge_Pa: float = 0.0,
        Cd: float = 0.61,
        gc: float = 1.0,
    ) -> dict:
        """Model 2: Flow of liquid through a hole in an unpressurized/pressurized draining tank."""
        A_hole = np.pi * (hole_diam_m / 2.0) ** 2
        A_tank = np.pi * (tank_diam_m / 2.0) ** 2

        # Initial release rate
        driving_head_init = (gc * P_gauge_Pa / rho_kg_m3) + (cls.G_ACCEL * h_L0_m)
        Qm_initial = (
            float(Cd * A_hole * rho_kg_m3 * np.sqrt(2.0 * driving_head_init))
            if driving_head_init > 0.0
            else 0.0
        )

        # Total emptying time (te)
        t_term1 = np.sqrt(2.0 * ((gc * P_gauge_Pa / rho_kg_m3) + (cls.G_ACCEL * h_L0_m)))
        t_term2 = np.sqrt(max(0.0, 2.0 * (gc * P_gauge_Pa / rho_kg_m3)))
        t_empty_sec = (1.0 / (Cd * cls.G_ACCEL)) * (A_tank / A_hole) * (t_term1 - t_term2)

        return {
            "initial_release_rate_kg_s": max(0.0, Qm_initial),
            "emptying_time_sec": max(0.0, float(t_empty_sec)),
            "hole_area_m2": A_hole,
            "tank_area_m2": A_tank,
        }

    @classmethod
    def liquid_pipe_flow(
        cls,
        rho_kg_m3: float,
        pipe_diam_m: float,
        pipe_length_m: float,
        delta_P_Pa: float,
        fanning_friction_factor: float = 0.005,
        fittings_K_loss: float = 0.0,
        gc: float = 1.0,
    ) -> float:
        """Model 3: Incompressible liquid discharge through a pipeline accounting for friction."""
        area = np.pi * (pipe_diam_m / 2.0) ** 2
        total_loss_coeff = (
            4.0 * fanning_friction_factor * (pipe_length_m / pipe_diam_m)
            + fittings_K_loss
            + 1.0
        )
        if delta_P_Pa <= 0.0 or total_loss_coeff <= 0.0:
            return 0.0
        velocity = np.sqrt(2.0 * gc * delta_P_Pa / (rho_kg_m3 * total_loss_coeff))
        return float(rho_kg_m3 * area * velocity)

    @classmethod
    def vapor_hole_expansion(
        cls,
        P0_Pa: float,
        P_amb_Pa: float,
        T0_K: float,
        hole_diam_m: float,
        MW_kg_mol: float,
        gamma: float = 1.4,
        Cd: float = 1.0,
        gc: float = 1.0,
    ) -> dict:
        """Model 4: Compressible vapor release through a hole (isentropic expansion)."""
        A_hole = np.pi * (hole_diam_m / 2.0) ** 2
        p_choked_ratio = (2.0 / (gamma + 1.0)) ** (gamma / (gamma - 1.0))
        p_choked = P0_Pa * p_choked_ratio

        if P_amb_Pa <= p_choked:
            regime = "Choked (Sonic) Flow"
            choked_term = (2.0 / (gamma + 1.0)) ** ((gamma + 1.0) / (gamma - 1.0))
            Qm = (
                Cd
                * A_hole
                * P0_Pa
                * np.sqrt((gamma * gc * MW_kg_mol / (cls.R_GAS * T0_K)) * choked_term)
            )
        else:
            regime = "Unchoked (Subsonic) Flow"
            p_ratio = P_amb_Pa / P0_Pa
            term = (p_ratio ** (2.0 / gamma)) - (p_ratio ** ((gamma + 1.0) / gamma))
            Qm = (
                Cd
                * A_hole
                * P0_Pa
                * np.sqrt(
                    (2.0 * gc * MW_kg_mol / (cls.R_GAS * T0_K))
                    * (gamma / (gamma - 1.0))
                    * term
                )
            )

        return {
            "mass_release_rate_kg_s": float(Qm),
            "flow_regime": regime,
            "choked_pressure_limit_Pa": float(p_choked),
        }

    @classmethod
    def gas_pipe_flow(
        cls,
        P0_Pa: float,
        P_amb_Pa: float,
        T0_K: float,
        pipe_diam_m: float,
        pipe_length_m: float,
        MW_kg_mol: float,
        fanning_friction_factor: float = 0.005,
        gamma: float = 1.4,
    ) -> dict:
        """Model 5: Isothermal compressible gas flow in a long pipeline."""
        area = np.pi * (pipe_diam_m / 2.0) ** 2
        friction_term = 4.0 * fanning_friction_factor * (pipe_length_m / pipe_diam_m)

        # Choked/Sonic limit pressure ratio at pipe outlet for isothermal flow
        # In isothermal flow, M=1 gives P_exit/P0 ~ 1 / sqrt(2*ln(P0/P_exit) + 4fL/D)
        # Using the standard engineering approximation for isothermal flow capacity:
        P_ratio_sq = (P_amb_Pa / P0_Pa) ** 2
        denom = friction_term + 2.0 * np.log(max(P0_Pa / P_amb_Pa, 1.0001))

        if denom <= 0:
            return {"mass_release_rate_kg_s": 0.0, "flow_regime": "Invalid Conditions"}

        mass_flux_G = np.sqrt(
            (MW_kg_mol * (P0_Pa**2) * (1.0 - P_ratio_sq)) / (cls.R_GAS * T0_K * denom)
        )
        Qm = area * mass_flux_G
        return {
            "mass_release_rate_kg_s": float(Qm),
            "flow_regime": "Compressible Pipe Flow",
        }

    @classmethod
    def flashing_liquid(
        cls,
        Qm_total_liquid_kg_s: float,
        T0_K: float,
        Tb_K: float,
        Cp_liquid_J_kgK: float,
        delta_Hv_J_kg: float,
    ) -> dict:
        """Model 6: Flashing liquid release with partial adiabatic flash vaporization."""
        if T0_K <= Tb_K:
            fv = 0.0
        else:
            fv = 1.0 - np.exp(-(Cp_liquid_J_kgK * (T0_K - Tb_K)) / delta_Hv_J_kg)
            fv = min(1.0, max(0.0, float(fv)))

        Qm_vapor = fv * Qm_total_liquid_kg_s
        Qm_liquid_entrained = (1.0 - fv) * Qm_total_liquid_kg_s

        return {
            "flash_fraction_fv": float(fv),
            "vapor_release_rate_kg_s": float(Qm_vapor),
            "residual_liquid_rate_kg_s": float(Qm_liquid_entrained),
        }

    @classmethod
    def pool_evaporation(
        cls,
        MW_kg_mol: float,
        mass_transfer_coeff_m_s: float,
        pool_area_m2: float,
        P_sat_Pa: float,
        T_pool_K: float,
    ) -> float:
        """Model 7: Volatile liquid pool evaporation or boiling on open ground."""
        Qm = (
            MW_kg_mol
            * mass_transfer_coeff_m_s
            * pool_area_m2
            * P_sat_Pa
            / (cls.R_GAS * T_pool_K)
        )
        return float(max(0.0, Qm))


# ==============================================================================
# DECISION MATRIX & QUESTIONNAIRE CONTROLLER
# ==============================================================================

def run_source_model_selector():
    """Interactive CLI Questionnaire routing user inputs through the decision matrix."""
    print("=" * 70)
    print("        SOURCE MODEL SELECTION & DISCHARGE CALCULATOR")
    print("=" * 70)

    print("\nQuestion 1: Physical State & Storage Condition")
    print("  [A] Liquid (stored at or below normal boiling point: T0 <= Tb)")
    print("  [B] Liquid stored under pressure above normal boiling point (T0 > Tb)")
    print("  [C] Vapor or Gas")
    print("  [D] Already spilled onto ground forming an open liquid pool")
    choice_q1 = input("Select an option (A/B/C/D): ").strip().upper()

    # Route based on Question 1
    if choice_q1 == "D":
        print("\n-> Selected Source Model: [7. Liquid Pool Evaporation or Boiling]")
        mw = float(input("Enter Molecular Weight (g/mol): ")) / 1000.0
        k = float(input("Enter Mass Transfer Coefficient K (m/s): "))
        area = float(input("Enter Pool Surface Area (m^2): "))
        p_sat = float(input("Enter Saturation Vapor Pressure P_sat (Pa): "))
        t_pool = float(input("Enter Pool Temperature (K): "))

        rate = SourceModels.pool_evaporation(mw, k, area, p_sat, t_pool)
        print("-" * 50)
        print(f"Calculated Evaporation Rate: {rate:.4f} kg/s ({rate * 1000:.2f} g/s)")

    elif choice_q1 == "B":
        print("\n-> Selected Source Model: [6. Flashing Liquids]")
        qm_tot = float(input("Enter Total Liquid Mass Discharge Rate (kg/s): "))
        t0 = float(input("Enter Storage Temperature T0 (K): "))
        tb = float(input("Enter Normal Boiling Point Tb (K): "))
        cp = float(input("Enter Average Liquid Heat Capacity Cp (J/kg*K): "))
        dhv = float(input("Enter Heat of Vaporization delta_Hv (J/kg): "))

        res = SourceModels.flashing_liquid(qm_tot, t0, tb, cp, dhv)
        print("-" * 50)
        print(f"Vapor Flash Fraction (fv): {res['flash_fraction_fv'] * 100:.2f}%")
        print(f"Vapor Generation Rate   : {res['vapor_release_rate_kg_s']:.4f} kg/s")
        print(f"Residual Liquid Rate    : {res['residual_liquid_rate_kg_s']:.4f} kg/s")

    elif choice_q1 == "A":
        print("\nQuestion 2: Containment & Leak Location (Liquid Release)")
        print("  [A1] Hole in vessel/tank where liquid level drops over time")
        print("  [A2] Orifice/hole with constant upstream pressure or fixed liquid head")
        print("  [A3] Severed or ruptured pipeline (pipe friction and length matter)")
        choice_q2 = input("Select an option (A1/A2/A3): ").strip().upper()

        rho = float(input("Enter Liquid Density (kg/m^3): "))

        if choice_q2 == "A1":
            print("\n-> Selected Source Model: [2. Flow of Liquid Through a Hole in a Tank]")
            d_hole = float(input("Enter Leak Hole Diameter (mm): ")) / 1000.0
            d_tank = float(input("Enter Tank Diameter (m): "))
            h_init = float(input("Enter Initial Liquid Height (m): "))
            p_gauge = float(input("Enter Overpressure / Vessel Gauge Pressure (Pa): "))

            res = SourceModels.liquid_tank_draining(rho, d_hole, d_tank, h_init, p_gauge)
            print("-" * 50)
            print(f"Initial Discharge Rate: {res['initial_release_rate_kg_s']:.4f} kg/s")
            print(f"Total Tank Drain Time : {res['emptying_time_sec']:.2f} seconds ({res['emptying_time_sec'] / 60:.2f} min)")

        elif choice_q2 == "A2":
            print("\n-> Selected Source Model: [1. Flow of Liquid Through a Hole]")
            d_hole = float(input("Enter Leak Hole Diameter (mm): ")) / 1000.0
            p_gauge = float(input("Enter Upstream Gauge Pressure (Pa): "))
            h_const = float(input("Enter Constant Liquid Head Above Hole (m): "))

            rate = SourceModels.liquid_hole_constant_head(rho, d_hole, p_gauge, h_const)
            print("-" * 50)
            print(f"Discharge Flow Rate: {rate:.4f} kg/s ({rate * 1000:.2f} g/s)")

        elif choice_q2 == "A3":
            print("\n-> Selected Source Model: [3. Flow of Liquids Through Pipes]")
            d_pipe = float(input("Enter Pipe Inner Diameter (mm): ")) / 1000.0
            l_pipe = float(input("Enter Pipe Length (m): "))
            dp = float(input("Enter Total Driving Pressure Drop (Pa): "))

            rate = SourceModels.liquid_pipe_flow(rho, d_pipe, l_pipe, dp)
            print("-" * 50)
            print(f"Pipe Mass Discharge Rate: {rate:.4f} kg/s ({rate * 1000:.2f} g/s)")
        else:
            print("Invalid sub-option for liquid release.")

    elif choice_q1 == "C":
        print("\nQuestion 3: Containment & Leak Location (Gas/Vapor Release)")
        print("  [C1] Puncture, crack, or relief orifice in a vessel/tank wall")
        print("  [C2] Long pipeline (gas expansion and pipe friction occur)")
        choice_q3 = input("Select an option (C1/C2): ").strip().upper()

        p0 = float(input("Enter Vessel/Upstream Absolute Pressure P0 (bar): ")) * 1e5
        p_amb = 101325.0
        t0 = float(input("Enter Absolute Storage Temperature T0 (K): "))
        mw = float(input("Enter Molecular Weight (g/mol): ")) / 1000.0
        gamma = float(input("Enter Specific Heat Ratio gamma (default 1.4): ") or "1.4")

        if choice_q3 == "C1":
            print("\n-> Selected Source Model: [4. Flow of Vapor Through Holes]")
            d_hole = float(input("Enter Hole Diameter (mm): ")) / 1000.0
            cd = float(input("Enter Discharge Coefficient Cd (default 1.0): ") or "1.0")

            res = SourceModels.vapor_hole_expansion(p0, p_amb, t0, d_hole, mw, gamma, cd)
            print("-" * 50)
            print(f"Flow Regime Identified : {res['flow_regime']}")
            print(f"Choked Transition Limit: {res['choked_pressure_limit_Pa'] / 1e5:.3f} bar")
            print(f"Vapor Release Rate     : {res['mass_release_rate_kg_s']:.4f} kg/s")

        elif choice_q3 == "C2":
            print("\n-> Selected Source Model: [5. Flow of Gases Through Pipes]")
            d_pipe = float(input("Enter Pipe Inner Diameter (mm): ")) / 1000.0
            l_pipe = float(input("Enter Pipe Length (m): "))

            res = SourceModels.gas_pipe_flow(p0, p_amb, t0, d_pipe, l_pipe, mw, gamma=gamma)
            print("-" * 50)
            print(f"Gas Pipe Discharge Rate: {res['mass_release_rate_kg_s']:.4f} kg/s")
        else:
            print("Invalid sub-option for gas release.")
    else:
        print("Invalid primary option selected.")


if __name__ == "__main__":
    run_source_model_selector()