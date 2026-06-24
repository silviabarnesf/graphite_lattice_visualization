import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Set global plot parameters
plt.rcParams['font.size'] = 12
plt.rcParams['axes.labelpad'] = 10
plt.rcParams['axes.titlesize'] = 16

# Lattice constants
a = 1.42  # Carbon-carbon distance in Angstroms
c = 3.34  # Interlayer distance in Angstroms

# Lattice vectors for graphene
a1 = a * np.array([3 / 2, np.sqrt(3) / 2])
a2 = a * np.array([3 / 2, -np.sqrt(3) / 2])

# Basis atoms
r1 = np.array([0, 0])
r2 = a * np.array([1, 0])

# Nearest neighbor displacements
neighbors = [
    np.array([a, 0]),
    a * np.array([1 / 2, np.sqrt(3) / 2]),
    a * np.array([1 / 2, -np.sqrt(3) / 2]),
    -np.array([a, 0]),
    -a * np.array([1 / 2, np.sqrt(3) / 2]),
    -a * np.array([1 / 2, -np.sqrt(3) / 2])
]

# Generate lattice points
num_points = 5
lattice_points = [m * a1 + n * a2
                  for m in range(-num_points, num_points)
                  for n in range(-num_points, num_points)]

# Generate layers with AB stacking
num_layers = 3
boundary = 5
all_atoms_A = []
all_atoms_B = []

for layer_idx in range(num_layers):
    # Calculate layer shift (ABA stacking)
    if layer_idx % 2 == 0:
        layer_shift = np.array([0.0, 0.0])  # A layer
    else:
        layer_shift = np.array([a, 0.0])  # B layer shifted by one bond length

    # Generate shifted lattice points
    shifted_points = [p + layer_shift for p in lattice_points]

    # Generate atoms with z-coordinate
    z = layer_idx * c
    atoms_A = np.array([np.hstack((p + r1, z)) for p in shifted_points])
    atoms_B = np.array([np.hstack((p + r2, z)) for p in shifted_points])

    # Apply boundary filtering
    filtered_A = [atom for atom in atoms_A
                  if -boundary <= atom[0] <= boundary
                  and -boundary <= atom[1] <= boundary]
    filtered_B = [atom for atom in atoms_B
                  if -boundary <= atom[0] <= boundary
                  and -boundary <= atom[1] <= boundary]

    all_atoms_A.extend(filtered_A)
    all_atoms_B.extend(filtered_B)

all_atoms_A = np.array(all_atoms_A)
all_atoms_B = np.array(all_atoms_B)
filtered_atoms = np.vstack((all_atoms_A, all_atoms_B))

# Create 3D plot
fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(111, projection='3d')

# Draw bonds
for atom in filtered_atoms:
    x, y, z = atom
    for n in neighbors:
        neighbor_pos = np.array([x + n[0], y + n[1], z])
        # Check if neighbor exists
        for pt in filtered_atoms:
            if (np.allclose(pt[:2], neighbor_pos[:2])
                    and np.isclose(pt[2], neighbor_pos[2])):
                ax.plot([x, neighbor_pos[0]],
                        [y, neighbor_pos[1]],
                        [z, neighbor_pos[2]],
                        color='gray', linewidth=0.5, alpha=0.5)
                break

# Draw atoms with layer-specific markers
for layer_idx in range(num_layers):
    layer_z = layer_idx * c
    mask_A = np.isclose(all_atoms_A[:, 2], layer_z)
    mask_B = np.isclose(all_atoms_B[:, 2], layer_z)

    marker = 'o' if layer_idx % 2 == 0 else 's'
    label_suffix = f' (Layer {layer_idx + 1})'

    ax.scatter(all_atoms_A[mask_A, 0], all_atoms_A[mask_A, 1], all_atoms_A[mask_A, 2],
               color='black', s=50, marker=marker,
               label=f'A Atoms{label_suffix}' if layer_idx < 2 else "")
    ax.scatter(all_atoms_B[mask_B, 0], all_atoms_B[mask_B, 1], all_atoms_B[mask_B, 2],
               color='grey', s=50, marker=marker,
               label=f'B Atoms{label_suffix}' if layer_idx < 2 else "")


# Find reference atoms
def find_nearest_atom(atom_list, target):
    if len(atom_list) == 0:
        raise ValueError("No atoms available for nearest search!")
    distances = np.linalg.norm(atom_list[:, :2] - target, axis=1)
    min_idx = np.argmin(distances)
    return atom_list[min_idx]


try:
    ref_atom_A = find_nearest_atom(all_atoms_A, [0, 0])
    ref_atom_B = find_nearest_atom(all_atoms_B, [a, 0])

    # Find upper layer reference
    upper_candidates = all_atoms_A[
        (np.abs(all_atoms_A[:, 0] - a) < 1e-2) &
        (np.abs(all_atoms_A[:, 1] - 0) < 1e-2) &
        (np.abs(all_atoms_A[:, 2] - c) < 1e-2)
        ]

    if len(upper_candidates) > 0:
        ref_atom_upper = upper_candidates[0]
    else:
        print("Warning: No exact upper reference atom found. Possible reasons:")
        print("- Increase num_points (current: {})".format(num_points))
        print("- Widen boundary (current: {})".format(boundary))
        ref_atom_upper = find_nearest_atom(all_atoms_A, [a, 0, c])

except ValueError as e:
    print("Critical Error:", str(e))
    exit()

# Draw and label distances
ax.plot([ref_atom_A[0], ref_atom_B[0]],
        [ref_atom_A[1], ref_atom_B[1]],
        [ref_atom_A[2], ref_atom_B[2]],
        color='black', linestyle='--')
ax.text((ref_atom_A[0] + ref_atom_B[0]) / 2,
        (ref_atom_A[1] + ref_atom_B[1]) / 2,
        0.2,  # Raise text above plane
        f'In-plane distance\na = {a} Å',
        color='black', fontsize=10, ha='center', va='bottom')

ax.plot([ref_atom_B[0], ref_atom_B[0]],
        [ref_atom_B[1], ref_atom_B[1]],
        [ref_atom_B[2], ref_atom_upper[2]],
        color='green', linestyle='--')
ax.text(ref_atom_B[0] + 0.5, ref_atom_B[1], c / 2,
        f'Interlayer spacing\nc = {c} Å',
        color='green', fontsize=10, ha='left', va='center')

# Add stacking visualization
for z in [0, c]:
    ax.plot([ref_atom_A[0], ref_atom_A[0]],
            [ref_atom_A[1], ref_atom_A[1]],
            [z, z + c],
            color='purple', linestyle=':', alpha=0.5)
ax.text(ref_atom_A[0], ref_atom_A[1], c / 2,
        'AB Stacking', color='purple', ha='right', va='center')

# Final formatting
ax.set_xlabel('X (Å)', fontsize=14)
ax.set_ylabel('Y (Å)', fontsize=14)
ax.set_zlabel('Z (Å)', fontsize=14)
ax.set_title(f"ABA-Stacked Graphite ({num_layers} Layers)", fontsize=16)
ax.set_box_aspect([1, 1, 0.3])
ax.view_init(elev=25, azim=-45)
ax.legend(loc='upper left', bbox_to_anchor=(0.8, 0.9))

# Clean up 3D axes
ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.grid(False)

plt.tight_layout()
plt.savefig("graphite_structure.png", dpi=300)
plt.show()
