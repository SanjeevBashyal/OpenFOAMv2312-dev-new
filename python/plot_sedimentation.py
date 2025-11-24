import pandas as pd
import matplotlib.pyplot as plt
import os

# Set style
plt.style.use('seaborn-v0_8-whitegrid')

# Read data
data_file = '/usr/lib/openfoam/openfoam2312/bashyal/Applications/srcTestMain/cfd_dem/freeSedimentation/sedimentation_data.csv'
df = pd.read_csv(data_file)

# Create output directory
output_dir = '/usr/lib/openfoam/openfoam2312/Thesis/Figures/Chapters/C05'
os.makedirs(output_dir, exist_ok=True)

# Plot Position
plt.figure(figsize=(10, 6))
plt.plot(df['Time'], df['Y'], label='Simulation', linewidth=2)
plt.xlabel('Time (s)', fontsize=12)
plt.ylabel('Vertical Position (m)', fontsize=12)
plt.title('Free Sedimentation: Position vs Time', fontsize=14)
plt.legend()
plt.grid(True)
plt.savefig(os.path.join(output_dir, 'sedimentation_position.png'))
plt.close()

# Plot Velocity
plt.figure(figsize=(10, 6))
plt.plot(df['Time'], df['Vy'], label='Simulation', linewidth=2)
plt.xlabel('Time (s)', fontsize=12)
plt.ylabel('Vertical Velocity (m/s)', fontsize=12)
plt.title('Free Sedimentation: Velocity vs Time', fontsize=14)
plt.legend()
plt.grid(True)
plt.savefig(os.path.join(output_dir, 'sedimentation_velocity.png'))
plt.close()

# Calculate Acceleration
# a = dV/dt
df['Ay'] = df['Vy'].diff() / df['Time'].diff()

# Plot Acceleration
plt.figure(figsize=(10, 6))
plt.plot(df['Time'], df['Ay'], label='Simulation', linewidth=2)
plt.xlabel('Time (s)', fontsize=12)
plt.ylabel('Vertical Acceleration (m/s²)', fontsize=12)
plt.title('Free Sedimentation: Acceleration vs Time', fontsize=14)
plt.legend()
plt.grid(True)
plt.savefig(os.path.join(output_dir, 'sedimentation_acceleration.png'))
plt.close()

print("Plots generated successfully.")
