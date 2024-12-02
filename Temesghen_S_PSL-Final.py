"""
Python standard library implementation of Monte Carlo simulation code
"""
import math
import random
import time

# all functions from day 4
def calculate_LJ(r_ij):
    r6_term = math.pow(1/r_ij, 6)
    r12_term = math.pow(r6_term, 2)
    pairwise_energy = 4 * (r12_term - r6_term)
    return pairwise_energy

def calculate_distance(coord1, coord2, box_length=None):
    """
    Calculate the distance between two 3D coordinates.

    Parameters
    ----------
    coord1, coord2: list
        The atomic coordinates

    Returns
    -------
    distance: float
        The distance between the two points.
    """

    distance = 0
    for i in range(3):
        dim_dist = coord1[i] - coord2[i]

        if box_length:
            dim_dist = dim_dist - box_length * round(dim_dist / box_length)

        dim_dist = dim_dist**2
        distance += dim_dist

    distance = math.sqrt(distance)
    return distance
    
def calculate_total_energy(coordinates, box_length, cutoff):
    """
    Calculate the total Lennard Jones energy of a system of particles.

    Parameters
    ----------
    coordinates : list
        Nested list containing particle coordinates.

    Returns
    -------
    total_energy : float
        The total pairwise Lennard Jones energy of the system of particles.
    """

    total_energy = 0

    num_atoms = len(coordinates)

    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):

            dist_ij = calculate_distance(
                coordinates[i], coordinates[j], box_length=box_length
            )

            if dist_ij < cutoff:
                interaction_energy = calculate_LJ(dist_ij)
                total_energy += interaction_energy

    return total_energy

def read_xyz(filepath):
    """
    Reads coordinates from an xyz file.
    
    Parameters
    ----------
    filepath : str
       The path to the xyz file to be processed.
       
    Returns
    -------
    atomic_coordinates : list
        A two dimensional list containing atomic coordinates

    """
    
    with open(filepath) as f:
        box_length = float(f.readline().split()[0])
        num_atoms = float(f.readline())
        coordinates = f.readlines()
    
    atomic_coordinates = []
    
    for atom in coordinates:
        split_atoms = atom.split()
        
        float_coords = []
        
        # We split this way to get rid of the atom label.
        for coord in split_atoms[1:]:
            float_coords.append(float(coord))
            
        atomic_coordinates.append(float_coords)
        
    
    return atomic_coordinates, box_length

def calculate_tail_correction(num_particles, box_length, cutoff):

    """
    calc the long range tail correction
    """
    const1 = (8 * math.pi * num_particles**2) / (3 * box_length**3)
    const2 = (1 / 3) * (1 / cutoff) ** 9 - (1 / cutoff) ** 3

    return const1 * const2

def generate_cubic_lattice(num_atoms, density):
  """
  Generate points on a cubic lattice using a desired final density.

  Parameters
  ----------
  num_atoms: int
    The number of atoms to place on the lattice.
  density: float
    The desired system density.

  Returns
  -------
  coords: list
    A nested list of generated coordinates.
  """

  # Calculate box length based on number of atoms and density.
  volume = num_atoms / density
  box_length = math.pow(volume, (1./3.))

  # Calculate the upper bound of cube size.
  # Our approach will be to place atoms until
  # we place all needed. For this, we need
  # to determine the maximum number of atoms on each
  # side.
  max_side = math.ceil(math.pow(num_atoms, (1./3.)))

  # Determine spacing based on number of atoms
  # and box length.
  spacing = box_length / max_side # units length / atom
  
  coordinates = []
  count = 0

  for i in range(max_side):
    for j in range(max_side):
      for k in range(max_side):
        coordinates.append([i*spacing, j*spacing, k*spacing])
        count += 1
        if count == num_atoms:
          return coordinates, box_length

def accept_or_reject(delta_e, beta):

    """
    accept or reject based on an energy and beta (inverse temp) measurement.

    parameters
    ----------
    delta_e: float
        an energy change in reduced units.
    beta: float
        inverse temperature

    returns
    -------
    bool
        true to accept move, false to reject
    """
    if delta_e <= 0.0:
        accept = True
    else:
        random_number = random.random()
        #the equation from the slides
        p_acc = math.exp(-delta_e*beta)

        if random_number < p_acc:
            accept = True
        else: 
            accept = False
    return accept

def calculate_pair_energy(coordinates, i_particle, box_length, cutoff):

    """
    calc the interaction energy of a particle with its envirnment (all other particles in the system).

    parameters
    ----------
    coordinates: list
        the coordinates for all particles in the system.
    i_particle: int
        the particle index for which to calc energy.
    cutofff: float
        the simulation cutoff. Beyond this distance, the interactions are not calculated.

    returns
    -------
    e_total: float
        the pairwise interaction energy of the i_th particle with all other particles
    """
    e_total = 0.0
    i_position = coordinates[i_particle]

    num_atoms = len(coordinates)

    for j_particle in range(num_atoms):
        if i_particle != j_particle:

            #indexing into the coordinates array and listnig out the particle coordinates
            j_position = coordinates[j_particle]

            r_ij = calculate_distance(i_position, j_position, box_length)

            if r_ij < cutoff:
                e_pair = calculate_LJ(r_ij)
                e_total += e_pair


    return e_total

# create a function called run_simulation that takes in simulation parameters and runs the simulation
def run_simulation(coordinates, box_length, cutoff, reduced_temperature, num_steps, max_displacement = 0.1, freq = 1000):

# calculated quantities
    beta = 1 / reduced_temperature
    num_particles = len(coordinates)

    delta_energy = 0

    total_energy = calculate_total_energy(coordinates, box_length, cutoff)
    total_energy += calculate_tail_correction(num_particles, box_length, cutoff)

    steps = []
    energies = []
    random.seed(0)

    # simulation code
    for step in range(num_steps):
        #1. Randomly pick one of N particles 
        random_particle = random.randrange(num_particles)

        #2. Calculate the interaction energy of the selected particle with the system and store this value
        current_energy = calculate_pair_energy(coordinates, random_particle, box_length, cutoff)

        
        #3. Generate a random x, y, z displacement
        x_rand = random.uniform(-max_displacement, max_displacement)
        y_rand = random.uniform(-max_displacement, max_displacement)
        z_rand = random.uniform(-max_displacement, max_displacement)
        
        #4. Modify the coordinate of the Nth particle by the generated displacements
        coordinates[random_particle][0] += x_rand
        coordinates[random_particle][1] += y_rand
        coordinates[random_particle][2] += z_rand
        
        #5. Calculate the interactionenergy of the selected particle with the system and store this value (using the new position)
        proposed_energy = calculate_pair_energy(coordinates, random_particle, box_length, cutoff)
        
        #6. Calculate the change in energy and decide to accept or reject the change
        delta_e = proposed_energy - current_energy

        accept = accept_or_reject(delta_e, beta)
        
        #7. if accept, keep particle movement. if reject, move particle back
        if accept:
            total_energy += delta_e
        else:
            # move is not accepted, roll back coordinates
            # why use - instead of + ??
            coordinates[random_particle][0] -= x_rand
            coordinates[random_particle][1] -= y_rand
            coordinates[random_particle][2] -= z_rand
            
        #8. print the energy if the step is a multiple of freq
        if step % freq == 0:
            print(step, total_energy/num_particles)
            steps.append(step)
            energies.append(total_energy)

# Test Parameters
num_particles = 100 # Number of particles in the simulation
density = 0.8
cutoff = 3.0 # Cutoff distance for interactions
reduced_temperature = 1
num_steps = [2000, 10000, 50000] # Different numbers of simulation steps to test
max_displacement = 0.1 # Maximum displacement for particle moves

#  Generate Cubic Lattice with coordinates and box length
coordinates, box_length = generate_cubic_lattice(num_particles, density)
   """
    Generate a cubic lattice of particles.

    Parameters
    ----------
    num_particles : int
        Number of particles in the lattice.
    density : float
        Density of the lattice (particles per unit volume).

    Returns
    -------
    tuple
        A tuple containing the particle coordinates and the box length.
    """
# Test Cases
def test_speed(num_steps, coordinates, box_length, cutoff, reduced_temperature):
    """
    Test the speed of the simulation for various numbers of steps.

    Parameters
    ----------
    num_steps : list of int
        List containing different numbers of steps to run in the simulation.
    coordinates : list of tuples
        Initial positions of particles in the simulation.
    box_length : float
        Length of the simulation box.
    cutoff : float
        Cutoff distance for interactions.
    reduced_temperature : float
        Reduced temperature for the simulation.

    Returns
    -------
    list of float
        A list of elapsed times for each number of steps tested.
    """
    times = []
    
    for i in num_steps:
        start = time.time() # Record the start time
        run_simulation(coordinates, box_length, cutoff, reduced_temperature, i, max_displacement)
        end = time.time() # Record the start time
        
        elapsed_time = end - start #calculate elapsed time
        times.append(elapsed_time) # store elapsed time
        print(f"It would take approximately {elapsed_time} seconds for {i} steps.")
    

    return times

# Run the test
times = test_speed(num_steps, coordinates, box_length, cutoff, reduced_temperature)


