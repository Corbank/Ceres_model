# ceres_model.py

import math

# --- Constants ---
AU_KM = 149597870.7  # Astronomical Unit in kilometers
G = 6.67430e-11     # Gravitational constant (m^3 kg^-1 s^-2)
M_SUN_KG = 1.98847e30 # Mass of the Sun (kg)

class Ceres:
    """
    A class representing the dwarf planet Ceres with its physical
    and orbital characteristics. Includes methods for calculating
    derived properties like volume, surface area, orbital speed,
    and escape velocity.

    Attributes:
        name (str): Name of the celestial body.
        type (str): Type of the celestial body (e.g., "Dwarf Planet").
        equatorial_radius_km (float): Equatorial radius in kilometers.
        polar_radius_km (float): Polar radius in kilometers.
        mean_radius_km (float): Mean radius in kilometers.
        mass_kg (float): Mass in kilograms.
        density_g_cm3 (float): Average density in g/cm^3.
        surface_gravity_m_s2 (float): Equatorial surface gravity in m/s^2.
        geometric_albedo (float): Geometric albedo (reflectivity).
        rotation_period_hours (float): Sidereal rotation period in hours.
        epoch (str): Epoch for the orbital elements (e.g., "J2000").
        semi_major_axis_au (float): Semi-major axis in Astronomical Units (AU).
        semi_major_axis_km (float): Semi-major axis in kilometers.
        orbital_period_days (float): Orbital period in Earth days.
        orbital_period_years (float): Orbital period in Earth years.
        eccentricity (float): Orbital eccentricity.
        inclination_deg (float): Orbital inclination to the ecliptic plane in degrees.
        longitude_ascending_node_deg (float): Longitude of the ascending node in degrees.
        argument_perihelion_deg (float): Argument of perihelion in degrees.
    """
    def __init__(self):
        # Physical Characteristics
        self.name = "Ceres"
        self.type = "Dwarf Planet"
        # Radii in kilometers (Ceres is an oblate spheroid)
        self.equatorial_radius_km = 487.3
        self.polar_radius_km = 454.7
        # Approximation for mean radius of an oblate spheroid
        self.mean_radius_km = (self.equatorial_radius_km * 2 + self.polar_radius_km) / 3
        # Mass in kilograms
        self.mass_kg = 9.393e20
        # Density in g/cm^3
        self.density_g_cm3 = 2.161
        # Surface gravity in m/s^2 (approximate)
        self.surface_gravity_m_s2 = 0.28
        # Geometric Albedo
        self.geometric_albedo = 0.090
        # Rotation period in hours
        self.rotation_period_hours = 9.074170

        # Orbital Characteristics (approximate values, relative to J2000 epoch)
        self.epoch = "J2000"
        # Semi-major axis in Astronomical Units (AU)
        self.semi_major_axis_au = 2.7673
        # Semi-major axis in kilometers
        self.semi_major_axis_km = self.semi_major_axis_au * AU_KM
        # Orbital period in Earth days
        self.orbital_period_days = 1681.63
        # Orbital period in Earth years
        self.orbital_period_years = self.orbital_period_days / 365.25
        # Eccentricity (how elliptical the orbit is)
        self.eccentricity = 0.075823
        # Inclination to the ecliptic plane in degrees
        self.inclination_deg = 10.593
        # Longitude of the ascending node in degrees
        self.longitude_ascending_node_deg = 80.305
        # Argument of perihelion in degrees
        self.argument_perihelion_deg = 73.939

    def get_volume_km3(self):
        """Calculates the approximate volume assuming an oblate spheroid.

        Returns:
            float: The approximate volume in cubic kilometers (km^3).
        """
        a = self.equatorial_radius_km
        c = self.polar_radius_km
        return (4/3) * math.pi * (a**2) * c

    def get_surface_area_km2(self):
        """Calculates the approximate surface area assuming an oblate spheroid.

        Uses the formula for the surface area of an oblate spheroid:
        Area = 2 * pi * a^2 * (1 + (c / (a * e)) * atanh(e))
        where e = sqrt(1 - (c^2 / a^2))

        Returns:
            float: The approximate surface area in square kilometers (km^2).
        """
        a = self.equatorial_radius_km
        c = self.polar_radius_km
        if a == 0: # Avoid division by zero if radius is zero
            return 0.0
        if a == c: # Perfect sphere case
             return 4 * math.pi * a**2

        e_squared = 1 - (c**2 / a**2)
        # Handle potential floating point inaccuracies leading to small negative e_squared
        if e_squared < 0:
            e_squared = 0

        e = math.sqrt(e_squared)
        # Avoid division by zero if e is zero (perfect sphere, handled above but good practice)
        if e == 0:
             return 4 * math.pi * a**2
        else:
            # Use the correct formula involving atanh for oblate spheroids
            return 2 * math.pi * a**2 * (1 + ((1 - e_squared) / e) * math.atanh(e))

    def get_escape_velocity_km_s(self):
        """Calculates the approximate escape velocity from the surface.

        Uses the formula: v_e = sqrt(2 * G * M / R), approximating with mean radius.

        Returns:
            float: The approximate escape velocity in kilometers per second (km/s).
        """
        if self.mean_radius_km <= 0: # Avoid division by zero
            return 0.0
        # Calculate using mean radius in meters
        r_m = self.mean_radius_km * 1000
        mass = self.mass_kg
        # v_e^2 = 2 * G * M / R
        escape_velocity_m_s_squared = 2 * G * mass / r_m
        escape_velocity_m_s = math.sqrt(escape_velocity_m_s_squared)
        return escape_velocity_m_s / 1000 # Convert m/s to km/s

    def get_perihelion_au(self):
        """Calculates the perihelion distance (closest approach to the Sun).

        Formula: q = a * (1 - e)

        Returns:
            float: The perihelion distance in Astronomical Units (AU).
        """
        return self.semi_major_axis_au * (1 - self.eccentricity)

    def get_aphelion_au(self):
        """Calculates the aphelion distance (farthest distance from the Sun).

        Formula: Q = a * (1 + e)

        Returns:
            float: The aphelion distance in Astronomical Units (AU).
        """
        return self.semi_major_axis_au * (1 + self.eccentricity)

    def get_orbital_speed_km_s(self, distance_au=None):
        """
        Calculates the orbital speed at a given distance from the Sun
        using the vis-viva equation.

        If no distance is provided, calculates the average orbital speed,
        approximated by the speed at the semi-major axis distance.

        Vis-viva equation: v^2 = GM_sun * (2/r - 1/a)

        Args:
            distance_au (float, optional): The distance from the Sun in AU.
                                           Defaults to None (uses semi-major axis).

        Returns:
            float: The orbital speed in kilometers per second (km/s).
                   Returns 0.0 if calculation results in an invalid state (e.g., negative sqrt).
        """
        a_m = self.semi_major_axis_km * 1000 # Semi-major axis in meters

        if distance_au is None:
            # Calculate speed at the semi-major axis distance for average speed approximation
            r_m = a_m
        else:
            if distance_au <= 0: # Distance must be positive
                return 0.0
            r_m = distance_au * AU_KM * 1000 # Distance in meters

        # Term inside vis-viva: (2/r - 1/a)
        vis_viva_term = (2 / r_m) - (1 / a_m)

        # Check if the term is physically valid (must be non-negative)
        # This can happen if r > 2a (object is beyond escape distance for this 'a')
        # or due to floating point issues for r very close to 2a.
        if vis_viva_term < 0:
             # This indicates an issue, potentially r is too large or a calculation error.
             # For bound orbits (like Ceres), r should always be <= a(1+e) < 2a.
             # Returning 0.0 or raising an error might be appropriate.
             # Let's return 0.0 for now, indicating an issue.
             print(f"Warning: Vis-viva term negative ({vis_viva_term}) for r={distance_au} AU, a={self.semi_major_axis_au} AU. Returning 0 speed.")
             return 0.0

        speed_m_s_squared = G * M_SUN_KG * vis_viva_term
        speed_m_s = math.sqrt(speed_m_s_squared)
        return speed_m_s / 1000 # Convert m/s to km/s

    def __str__(self):
        """Returns a formatted string representation of the Ceres object."""
        # Corrected indentation for the return statement and f-string block
        return (
            f"--- {self.name} ({self.type}) ---\n"
            f"Physical Characteristics:\n"
            f"  Mean Radius: {self.mean_radius_km:.1f} km\n"
            f"  Equatorial Radius: {self.equatorial_radius_km:.1f} km\n"
            f"  Polar Radius: {self.polar_radius_km:.1f} km\n"
            f"  Mass: {self.mass_kg:.3e} kg\n"
            f"  Volume: {self.get_volume_km3():.3e} km^3\n"
            f"  Surface Area: {self.get_surface_area_km2():.3e} km^2\n"
            f"  Density: {self.density_g_cm3:.3f} g/cm^3\n"
            f"  Surface Gravity: {self.surface_gravity_m_s2:.2f} m/s^2\n"
            f"  Escape Velocity: {self.get_escape_velocity_km_s():.3f} km/s\n"
            f"  Geometric Albedo: {self.geometric_albedo:.3f}\n"
            f"  Rotation Period: {self.rotation_period_hours:.4f} hours\n"
            f"Orbital Characteristics (Epoch: {self.epoch}):\n"
            f"  Semi-major Axis: {self.semi_major_axis_au:.4f} AU ({self.semi_major_axis_km:.3e} km)\n"
            f"  Perihelion: {self.get_perihelion_au():.4f} AU\n"
            f"  Aphelion: {self.get_aphelion_au():.4f} AU\n"
            f"  Orbital Period: {self.orbital_period_days:.2f} days ({self.orbital_period_years:.3f} years)\n"
            f"  Eccentricity: {self.eccentricity:.6f}\n"
            f"  Inclination: {self.inclination_deg:.4f} degrees\n"
            f"  Longitude of Asc. Node: {self.longitude_ascending_node_deg:.3f} degrees\n"
            f"  Argument of Perihelion: {self.argument_perihelion_deg:.3f} degrees\n"
            f"  Avg. Orbital Speed: {self.get_orbital_speed_km_s():.3f} km/s\n"
        )

# Example usage:
if __name__ == "__main__":
    ceres = Ceres()
    print(ceres)

    # Calculate speed at perihelion
    perihelion_au = ceres.get_perihelion_au()
    speed_at_perihelion = ceres.get_orbital_speed_km_s(distance_au=perihelion_au)
    print(f"Orbital speed at perihelion ({perihelion_au:.4f} AU): {speed_at_perihelion:.3f} km/s")

    # Calculate speed at aphelion
    aphelion_au = ceres.get_aphelion_au()
    speed_at_aphelion = ceres.get_orbital_speed_km_s(distance_au=aphelion_au)
    print(f"Orbital speed at aphelion ({aphelion_au:.4f} AU): {speed_at_aphelion:.3f} km/s")

    # Example: Calculate speed at a specific distance (e.g., 2.8 AU)
    specific_distance_au = 2.8
    speed_at_specific_distance = ceres.get_orbital_speed_km_s(distance_au=specific_distance_au)
    print(f"Orbital speed at {specific_distance_au} AU: {speed_at_specific_distance:.3f} km/s")

    # Print escape velocity
    print(f"Escape Velocity: {ceres.get_escape_velocity_km_s():.3f} km/s")
