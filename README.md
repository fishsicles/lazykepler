This is the source code for what I've taken to calling the Lazy Writer's Kepler Approximant, a tool designed to compute interplanetary distances using Kepler's laws of planetary motion. Intended as a tool for my own use on sci-fi writing forums, it's also available here. Do with it what you will.

# Usage
Run `lazykepler.py filename.yaml`, which loads configuration and orbital information from the file `filename.yml`; for an example YAML file, see the one for our solar system in the main repository.

Calculations use numpy and [https://pypi.org/project/kepler.py/ kepler.py] by danfm. pyyaml is needed to parse input YAML.

# Current Commands
* **`ctime`:** Calculates the time to travel between two bodies at lightspeed at the reference time.
* **`distance`:** Computes the distance between two bodies at the reference time. If only given one body, instead calculates based off of the origin (focus of orbits).
* **`time`:** Sets reference time if provided with a value, otherwise prints current reference time.
* **`where`:** Gives (x,y,z) coordinates of the given body at the reference time.
