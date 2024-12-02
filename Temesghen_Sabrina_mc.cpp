// MONTE CARLO INTEGRATION

//includes libraries
#include <iostream> // for std::cout (input/output)
#include <cmath> // for std::pow and std::sqrt (math equations)
#include <tuple> // for std::pair and std::tuple 
#include <fstream> // for reading/writing files
#include <array>   // for std::array (fixed arrays)
#include <vector>  // for std::vector (dynamic arrays)
#include <utility> // for std::pair 
#include <random> // for random numbers
#include <chrono> // for generating the seed

// typedefs for coordinates
typedef std::array<double, 3> AtomCoord; // A 3-element array representing (x, y, z) coordinates of an atom
typedef std::vector<AtomCoord> Coordinates; // A vector of AtomCoord, representing positions of multiple atoms

// random number generator
std::default_random_engine re;

// calculates the lennard-jones potential 
double calculate_LJ(double r_ij)
/*
calculates the pairwise Lennard Jones interaction energy based on the separation distance in reduced units.

parameters:
-----------
r_ij: double
    The distance between the particles in reduced units

return:
-------
pairwise_energy: double
    The pairwise Lennard Jones interaction energy in reduced units
*/
{
    double r6_term = std::pow(1.0 / r_ij, 6);
    double r12_term = std::pow(r6_term, 2);
    double pairwise_energy = 4 * (r12_term - r6_term);
    return pairwise_energy;
}


double calculate_distance(const AtomCoord& coord1, const AtomCoord& coord2, double box_length = 0.0)
/*
Calculate the distance between two 3D coordinates.

    parameters: 
    ----------
    coord1, coord2: array
        The atomic coordinates
    
    box_length: double
        The length of the cubic box

    returns:
    -------
    distance: double
        The distance between the two points
*/
{
    double distance = 0.0; // initializes a counter

    for(int i = 0; i < 3; i++) // loops through the (x, y, z) coordinates
    {
        double dim_dist = coord1[i] - coord2[i]; //calculates the distance between 2 coordinates

        if (box_length > 0) 
        {
            dim_dist = dim_dist - box_length * std::round(dim_dist / box_length); // as long as box_length is greater than 0; distance is calculated within the boundary 
        }
        dim_dist = dim_dist * dim_dist;
        distance += dim_dist;
    }
    distance = std::sqrt(distance);
    return distance;
}

double calculate_total_energy(Coordinates&coordinates, double box_length, double cutoff)
/*
calculate the total Lennard Jones energy of a system of particles by summing up interactions between pairs of atoms within a specific cutoff.

    parameters:
    ----------
    coordinates: Coordinates (vector<AtomCoord>)
        a vector of arrays, where each array represents the coordinates of an atom
       
    
    box_length: double
        the length of the cubic box

    cutoff: double
        the simulation cutoff. Beyond this distance, the interactions are not calculated


    returns:
    -------
    total_energy: double
        the total pairwise Lennard Jones energy of the system of particles
*/
{
    double total_energy = 0.0; // initialize a counter
    int num_atoms = coordinates.size(); // get number of atoms

    for(int i = 0; i < num_atoms; i++) // iterates through each atom
    { 
        for (int j = i + 1; j < num_atoms; j++) // iterates through each atom from i+1 to avoid duplication
        {
            double distance_ij = calculate_distance(coordinates[i], coordinates[j], box_length); // calculates distance between atoms

            if (distance_ij < cutoff) // takes into account the cutoff
            {
                double interaction_energy = calculate_LJ(distance_ij); // calculates the LJ potential energy within the cutoff
            }
        }
    }

return total_energy; 

}


std::pair<Coordinates, double> read_xyz(std::string file_path)
/*
Read atomic coordinates and box length from an XYZ file format.
 
    parameters:
    ----------
    file_path: std::string
        the path to the XYZ file

    returns:
    -------
    std::pair<Coordinates, double>
        a pair containing the coordinates of the atoms and the length of the simulation box
 
*/
{
    std::ifstream infile(file_path); // to open file

    if(!infile.is_open()) // check if file was successfully opened
    {
        throw std::runtime_error("File path in read_xyz does not exist!");
    }

    double dummy; 
    double box_length;
    int num_atoms;

    infile >> box_length >> dummy >> dummy;

    infile >> num_atoms;

    Coordinates coords;

    for(int i = 0; i < num_atoms; i++)
    {
        AtomCoord coord;

        // Throws away the atom index
        infile >> dummy >> coord[0] >> coord[1] >> coord[2];

        // Add to vector
        coords.push_back(coord);
    }

    // Makes an appropriate pair object
    return std::make_pair(coords, box_length);
}


double calculate_tail_correction(int num_particles, double box_length, double cutoff)
/*
calculates the pairwise Lennard Jones interaction energy based on the separation distance in reduced units.

parameters:
-----------
num_particles: int
    the number of particles

box_length: double
    the length of the cubic box

cutoff: double
    the simulation cutoff. Beyond this distance, the interactions are not calculated

return:
-------
 const1 * const2: double
    the correction to energy calculation to take into account the cutoff
*/
{
    const double pi = 3.141592653589793; // define pi as a constant
    double const1 = (8 * pi * std::pow(num_particles, 2)) / (3 * std::pow(box_length, 3)); // volume and number of particles 
    double const2 = (1 / 3.0) * std::pow(1 / cutoff, 9) - std::pow(1 / cutoff, 3); // integral 

    return const1 * const2;
}

std::pair<Coordinates, double> generate_cubic_lattice(int num_atoms, double density) 
/*
generate points on a cubic lattice using a desired final density.

parameters:
-----------
num_atoms: int
    the number of atoms to place on the lattice
    
    density: double
    the desired system density

returns
-------
std::pair<coordinates, box_length>
    a pair containing the coordinates of the atoms on the lattice and the length of the box
*/
{
    double volume = num_atoms / density;
    double box_length = std::pow(volume, 1.0 / 3.0);

    // Max number of atoms on each side of the cube.
    int max_side = std::ceil(std::pow(num_atoms, 1.0 / 3.0));

    // Spacing of number of atoms and box length.
    double spacing = box_length / max_side;

    Coordinates coordinates;
    coordinates.reserve(num_atoms);

    int count = 0;
    for (int i = 0; i < max_side; ++i) 
    {
        for (int j = 0; j < max_side; ++j) 
        {
            for (int k = 0; k < max_side; ++k) 
            {
                if (count >= num_atoms) 
                {
                    return std::make_pair(coordinates, box_length);
                }
                AtomCoord temp;
                temp[0] = i * spacing;
                temp[1] = j * spacing;
                temp[2] = k * spacing;
                coordinates.push_back(temp);
                ++count;
            }
        }
    }

    // In case we have not yet returned, ensure we return correctly
    return std::make_pair(coordinates, box_length);
}

// given for random generation
double random_double(double lower_bound, double upper_bound) {
    std::uniform_real_distribution<double> dist(lower_bound, upper_bound);
    return dist(re);
}
int random_integer(int lower_bound, int upper_bound)
{           
   //dist will return [a,b] but we want [a,b)
   std::uniform_int_distribution<int> dist(lower_bound, upper_bound-1);
   return dist(re);
}  

bool accept_or_reject(double delta_e, double beta)
/*
accept or reject based on an energy and beta (inverse temp) measurement.

parameters
----------
delta_e: double
    an energy change in reduced units.
    
beta: double 
    inverse temperature

returns
-------
bool
    true to accept move, false to reject
*/
{
    if (delta_e <= 0.0) 
    {
        return true;
    } else 
    {
        // generate a random number between 0 and 1
        double random_number = random_double(0.0, 1.0);

        // calculate the acceptance probability
        double p_acc = std::exp(-delta_e * beta);

        // accept or reject based on the probability
        return random_number < p_acc;
    }
}

double calculate_pair_energy(const std::vector<AtomCoord>& coordinates, int i_particle, double box_length, double cutoff) 
/*
calc the interaction energy of a particle with its envirnment (all other particles in the system).

parameters
----------
coordinates: const std::vector<AtomCoord>&
    coordinates for all particles in the system

i_particle: int
    the particle index for which to calc energy

box_length: double
    length of cubic box

cutofff: double
    the simulation cutoff. Beyond this distance, the interactions are not calculated

returns
-------
e_total: double
    the pairwise interaction energy of the i_th particle with all other particles
*/
{
    double e_total = 0.0;
    const AtomCoord& i_position = coordinates[i_particle];
    int num_atoms = coordinates.size();
    
    for (int j_particle = 0; j_particle < num_atoms; j_particle++) 
    {
        if (i_particle != j_particle) {
            const AtomCoord& j_position = coordinates[j_particle];
            double r_ij = calculate_distance(i_position, j_position, box_length);
            if (r_ij < cutoff) {
                double e_pair = calculate_LJ(r_ij);
                e_total += e_pair;
            }
        }
    }
    
    return e_total;
}

void run_simulation(Coordinates& coordinates, double box_length, double cutoff, double reduced_temperature, int num_steps, double max_displacement = 0.1, int freq = 1000)
{
    double beta = 1.0 / reduced_temperature;
    int num_particles = coordinates.size();

    double total_energy = calculate_total_energy(coordinates, box_length,cutoff);
    total_energy += calculate_tail_correction(num_particles,box_length,cutoff);

    for(int step = 0; step < num_steps; step++)
    {
        int random_particle = random_integer(0, num_particles);
        double current_energy = calculate_pair_energy(coordinates, random_particle, box_length, cutoff);

        double x_rand = random_double(-max_displacement, max_displacement);
        double y_rand = random_double(-max_displacement, max_displacement);
        double z_rand = random_double(-max_displacement, max_displacement);

        coordinates[random_particle][0] += x_rand;
        coordinates[random_particle][1] += y_rand;
        coordinates[random_particle][2] += z_rand;

        double proposed_energy = calculate_pair_energy(coordinates, random_particle, box_length, cutoff);

        double delta_e = proposed_energy - current_energy;
        bool accept = accept_or_reject(delta_e, beta);
        
        if (accept) 
        {
            total_energy += delta_e;
        } else 
        {
            coordinates[random_particle][0] -= x_rand;
            coordinates[random_particle][1] -= y_rand;
            coordinates[random_particle][2] -= z_rand;
        }
        if (step % freq == 0) 
        {
            double avg_energy = total_energy / num_particles;
            std::cout << step << " " << avg_energy << std::endl;
        }

    }
}
int main(void)