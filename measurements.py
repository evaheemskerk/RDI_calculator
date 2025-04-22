import numpy as np

def centroid(points):
    return np.mean(points, axis=0)

def distance(particle1, particle2):
    return float(np.linalg.norm(np.array(particle1) - np.array(particle2)))

def radius_of_gyration(coordinates):
    centroid_points = centroid(coordinates)
    RG_cell = np.sqrt(sum([distance(i, centroid_points)**2 for i in coordinates])
                      /len(coordinates))
    return RG_cell

def polarization_index(pixels, particles):
    #print('centroid of particles =', centroid(particles))
    RG_cell = radius_of_gyration(pixels)
    pixels = np.array(pixels)
    particles = np.array(particles)
    PI = (np.sqrt((centroid(particles[:,0])-centroid(pixels[:,0]))**2 +
                 (centroid(particles[:,1])-centroid(pixels[:,1]))**2 +
                 (centroid(particles[:,2])-centroid(pixels[:,2]))**2))/RG_cell

    return PI


def dispersion_index(pixels, particles):
    x_rna, y_rna, z_rna = centroid(particles)
    N = len(particles)
    mu2 = sum([((x - x_rna) ** 2 + (y - y_rna) ** 2 + (z - z_rna) ** 2) for (x, y, z) in particles]) / N
    M = len(pixels)
    mu2_prime = sum([((X - x_rna) ** 2 + (Y - y_rna) ** 2 + (Z - z_rna) ** 2) for (X, Y, Z) in pixels]) / M
    DI = mu2 / mu2_prime
    return DI


def peripheral_distribution_index(pixels_cyto, pixels_nucleus, particles):
    x_nuc, y_nuc, z_nuc = centroid(pixels_nucleus)

    N = len(particles)
    mu2 = sum([((x - x_nuc) ** 2 + (y - y_nuc) ** 2 + (z - z_nuc) ** 2) for (x, y, z) in particles]) / N

    M = len(pixels_cyto)
    mu2_prime = sum([((X - x_nuc) ** 2 + (Y - y_nuc) ** 2 + (Z - z_nuc) ** 2) for (X, Y, Z) in pixels_cyto]) / M

    PDI = mu2 / mu2_prime
    return PDI

