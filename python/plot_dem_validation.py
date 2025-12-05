import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# Paths
data_dir = "/usr/lib/openfoam/openfoam2312/bashyal/Applications/srcTestMain/dem"
output_dir = "/usr/lib/openfoam/openfoam2312/Thesis/Figures/Chapters/C05"

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# ==========================================
# 1. Hysteretic Bounce Test
# ==========================================
try:
    df_bounce = pd.read_csv(os.path.join(data_dir, "bounce_data.csv"))
    
    # Plot Height vs Time
    plt.figure(figsize=(8, 5))
    plt.plot(df_bounce['Time'], df_bounce['Z'], label='Height (Z)', color='blue')
    plt.xlabel('Time (s)')
    plt.ylabel('Height (m)')
    plt.title('Hysteretic Bounce: Height vs Time')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "bounce_height.png"), dpi=300)
    plt.close()

    # Plot Velocity vs Time
    plt.figure(figsize=(8, 5))
    plt.plot(df_bounce['Time'], df_bounce['Vz'], label='Velocity (Vz)', color='red')
    plt.xlabel('Time (s)')
    plt.ylabel('Vertical Velocity (m/s)')
    plt.title('Hysteretic Bounce: Velocity vs Time')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "bounce_velocity.png"), dpi=300)
    plt.close()
    print("Bounce plots generated.")
except Exception as e:
    print(f"Error processing bounce data: {e}")

# ==========================================
# 2. Sliding Block Test
# ==========================================
try:
    df_sliding = pd.read_csv(os.path.join(data_dir, "sliding_data.csv"))
    
    # Filter data after t=0.1 to avoid initial transient
    mask = df_sliding['Time'] > 0.1
    t_fit = df_sliding.loc[mask, 'Time']
    v_fit = df_sliding.loc[mask, 'VelMag']
    
    # Linear fit for acceleration
    coeffs = np.polyfit(t_fit, v_fit, 1)
    acc_sim = coeffs[0]
    
    # Theoretical Acceleration
    g = 9.81
    theta = 30 * np.pi / 180
    mu = 0.3
    acc_theory = g * (np.sin(theta) - mu * np.cos(theta))
    
    # Plot Position vs Time
    plt.figure(figsize=(8, 5))
    plt.plot(df_sliding['Time'], df_sliding['PosX'], label='Position X', color='green')
    plt.xlabel('Time (s)')
    plt.ylabel('Position X (m)')
    plt.title('Sliding Block: Position vs Time')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "sliding_position.png"), dpi=300)
    plt.close()

    # Plot Velocity vs Time
    plt.figure(figsize=(8, 5))
    plt.plot(df_sliding['Time'], df_sliding['VelMag'], label='Simulated Velocity', color='orange', linewidth=2)
    plt.plot(df_sliding['Time'], coeffs[0]*df_sliding['Time'] + coeffs[1], 'k--', label=f'Fit (a={acc_sim:.2f} m/s²)')
    plt.plot([], [], ' ', label=f'Theory (a={acc_theory:.2f} m/s²)') # Dummy for legend
    plt.xlabel('Time (s)')
    plt.ylabel('Velocity Magnitude (m/s)')
    plt.title('Sliding Block: Velocity vs Time')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "sliding_velocity.png"), dpi=300)
    plt.close()
    print(f"Sliding plots generated. Acc Sim: {acc_sim:.3f}, Theory: {acc_theory:.3f}")
except Exception as e:
    print(f"Error processing sliding data: {e}")

# ==========================================
# 3. Random Packing Test
# ==========================================
try:
    df_packing = pd.read_csv(os.path.join(data_dir, "packing_data.csv"))
    
    # Plot Kinetic Energy vs Time
    plt.figure(figsize=(8, 5))
    plt.plot(df_packing['Time'], df_packing['KE'], label='Kinetic Energy', color='purple')
    plt.xlabel('Time (s)')
    plt.ylabel('Kinetic Energy (J)')
    plt.title('Random Packing: Kinetic Energy vs Time')
    plt.yscale('log')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "packing_ke.png"), dpi=300)
    plt.close()

    # Plot Solid Fraction vs Time
    plt.figure(figsize=(8, 5))
    plt.plot(df_packing['Time'], df_packing['SolidFraction'], label='Solid Fraction', color='brown')
    plt.xlabel('Time (s)')
    plt.ylabel('Solid Fraction')
    plt.title('Random Packing: Solid Fraction vs Time')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.axhline(y=0.50, color='k', linestyle=':', label='Target ~0.50')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "packing_phi.png"), dpi=300)
    plt.close()
    
    # Plot Coordination Number vs Time
    plt.figure(figsize=(8, 5))
    plt.plot(df_packing['Time'], df_packing['CoordNum'], label='Coordination Number', color='black')
    plt.xlabel('Time (s)')
    plt.ylabel('Avg Coordination Number')
    plt.title('Random Packing: Coordination Number vs Time')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "packing_coord.png"), dpi=300)
    plt.close()
    
    print("Packing plots generated.")
except Exception as e:
    print(f"Error processing packing data: {e}")

try:
    df_packing = pd.read_csv(os.path.join(data_dir, "walton_braun_data.csv"))
    
    # Plot Kinetic Energy vs Time
    plt.figure(figsize=(8, 5))
    plt.plot(df_packing['Time'], df_packing['KE'], label='Kinetic Energy', color='purple')
    plt.xlabel('Time (s)')
    plt.ylabel('Kinetic Energy (J)')
    plt.title('Random Packing: Kinetic Energy vs Time')
    plt.yscale('log')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "packing_ke_walter.png"), dpi=300)
    plt.close()

    # Plot Solid Fraction vs Time
    plt.figure(figsize=(8, 5))
    plt.plot(df_packing['Time'], df_packing['SolidFraction'], label='Solid Fraction', color='brown')
    plt.xlabel('Time (s)')
    plt.ylabel('Solid Fraction')
    plt.title('Random Packing: Solid Fraction vs Time')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.axhline(y=0.50, color='k', linestyle=':', label='Target ~0.50')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "packing_phi_walter.png"), dpi=300)
    plt.close()
    
    # Plot Coordination Number vs Time
    plt.figure(figsize=(8, 5))
    plt.plot(df_packing['Time'], df_packing['CoordNum'], label='Coordination Number', color='black')
    plt.xlabel('Time (s)')
    plt.ylabel('Avg Coordination Number')
    plt.title('Random Packing: Coordination Number vs Time')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "packing_coord_walter.png"), dpi=300)
    plt.close()
    
    print("Packing plots generated.")
except Exception as e:
    print(f"Error processing packing data: {e}")
