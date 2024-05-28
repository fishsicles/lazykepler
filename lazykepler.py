import errno, os
import sys
import math
import random
import numpy
import kepler # kepler.py library provides solutions to Kepler's equation for true anomaly

from enum import IntEnum

from util.togglelog import ToggleLog

# set up reading orbit info from YAML file
import yaml
try:
   from yaml import CLoader as Loader
except ImportError:
   from yaml import Loader, Dumper

class Axis(IntEnum):
   # by conventional orbital dynamics this is the direction from the Sun to Earth at the spring equinox. 
   X = ReferenceDirection = 0;
   Y = 1; # we don't actually use this but it's here for completion's sake
   Z = ReferencePlaneNormal = 2;

def rotationMatrix (theta, axis):
   """
      Constructs the 3D rotation matrix about the given Axis.
      
      Parameters:
       - theta: Angle of rotation about the axis. In degrees. I hate it too but it's astronomy's fault.
       - axis: Axis enum constant indicating the axis to be used.
   """
   ctheta, stheta = math.cos(theta * math.pi/180), math.sin(theta * math.pi/180)
   output = [[ctheta, -stheta],[stheta, ctheta]]
   output[0].insert(axis,0)
   output[1].insert(axis,0)
   output.insert(axis,[1 if i==axis else 0 for i in range(3)])
   return numpy.array(output)

class Orbit:
   def __init__ (self,
      perihelion,
      aphelion,
      eccentricity,
      argumentOfPeriapsis,
      inclination,
      longitudeOfAscendingNode, # 
      period, # period of orbit. any units of time work, as long as they are consistent across all bodies. sample data is in days
      meanAnomaly=None, # mean anomaly at time t=0. not provided, will be randomised according to provided seed.
      visual=None # pyplotlib visual for planet
   ):
      # coordinates in elliptical plane are first calculated with periapsis on the x-axis and the sun at 0,0
      self.semimajorAxis = (perihelion + aphelion)/2
      self.eccentricity = eccentricity
      self.period = period

      # get preset starting mean anomaly or generate from seed
      self.meanAnomaly = meanAnomaly
      if self.meanAnomaly is None:
         self.meanAnomaly = random.uniform(0,360)
      
      # this matrix converts coordinates (a,b,0) in the orbital plane to coordinates (x,y,0) 
      ## first/rightmost matrix - rotate plane by argument of periapsis to align ascending node on positive x-axis
      self.matrix = rotationMatrix(argumentOfPeriapsis, Axis.ReferencePlaneNormal)
      ## second matrix - rotate by inclination about the reference direction to tilt orbit relative to reference plane
      self.matrix = numpy.matmul(rotationMatrix(inclination, Axis.ReferenceDirection), self.matrix)
      ## third/leftmost matrix - rotate plane by longitude of ascending node
      ## to put ascending node in right place relative to reference direction
      self.matrix = numpy.matmul(rotationMatrix(longitudeOfAscendingNode, Axis.ReferencePlaneNormal), self.matrix)

      # visual used for single-location visualiser plot
      self.visual = visual
   
   def orbitalPosition(self, time):
      """
         Calculates the body's position in the orbit at a given time after starting
         conditions. 
      """
      # calculate mean anomaly at given time
      meanAnomaly = 360/self.period * (time % self.period) + self.meanAnomaly
      # find true anomaly
      eAnomaly, cTAnomaly, sTAnomaly = kepler.kepler(meanAnomaly, self.eccentricity)
      # calculate radius
      radius = self.semimajorAxis * (1-self.eccentricity ** 2)/(1+self.eccentricity*cTAnomaly)
      # get position in orbital plane
      return cTAnomaly * radius, sTAnomaly * radius
   
   def position (self, time):
      # get orbital position
      planeX, planeY = self.orbitalPosition(time)
      # convert to column vector, still in orbital plane but now with z-coordinate
      vector = numpy.array([[planeX, planeY, 0]]).T
      # multiply matrix by row vector to transform into proper 3d coordinate
      return numpy.matmul(self.matrix, vector)

def magnitude (v):
   return math.sqrt(sum([x**2 for x in v]))

knownC = {
   "day": {"AU": 173.1},
   "s": {"m": 299792000}
}
knownG = {
   "day": {"AU": 0.489},
   "s": {"m": 9.8}
}

def loadOrbits (filename, log=True):
   with ToggleLog(log) as logger:
      data = None
      with open(filename) as infile:
         print(f"Loading system data from {filename}")
         data = yaml.load(infile,Loader=Loader)
      if data is None:
         # i always have to google how to be a cop about doing this properly
         raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), filename)

      settings = data["settings"] if "settings" in data else {}

      if "seed" in settings:
         print(f"Setting up random number generator with seed {settings['seed']}")
         random.seed(settings["seed"])
      else:
         seed = random.random()
         random.seed(seed)
         print(f"RNG seed is {seed}")

      units = data["units"]
      if "time" not in units:
         print("No time unit specified, defaulting to Earth days")
         units["time"] = "day"
      if "dist" not in units:
         print("No distance unit specified, defaulting to AU")
         units["dist"] = "AU"
      print("Loaded units:")
      print(f"   Distance: {units['dist']}")
      print(f"   Time: {units['time']}")

      # try and get constants for units if not provided
      if "c" in units:
         print(f"   c = {units['c']} {units['dist']}/{units['time']}")
      else:
         try:
            units["c"] = knownC[units["time"]][units["dist"]]
         except KeyError:
            print(f"Speed of light not provided for units {units['dist']}/{units['time']}. Communication times cannot be calculated.")
      if "g" in units:
         print (f"\tg = {units['g']}")
      else:
         try:
            units["g"] = knownG[units["time"]][units["dist"]]
         except KeyError:
            print (f"Acceleration on Earth not given for units {units['dist']}/{units['time']}^2. Acceleration cannot be given in g.")

      # now process actual orbits data
      print("Loading orbits:")
      orbits = {}
      for name, rawOrbit in data["orbits"].items():
         name = name.title()
         print(f"   Loading orbit of {name}")
         orbits[name] = Orbit(**rawOrbit)
         print(f"      mean anomaly {'manually' if 'meanAnomaly' in rawOrbit else 'randomly'} set to {orbits[name].meanAnomaly}")
      return units, orbits, settings

def sanitise (orbits, *bodies):
   output = []
   notFound = []
   for body in bodies:
      body = body.title()
      if body not in orbits:
         print(f"no such body: {body}")
         notFound.append(body)
      else:
         output.append(body)
   return output, notFound

def posAtTime (body, orbits, time=0):
   body = body.title()
   return orbits[body].position(time)

def distAtTime (b1, b2, orbits, time=0):

   return magnitude(posAtTime(b1,orbits,time) - posAtTime(b2,orbits,time))

def getDist (units, orbits, settings, time, b1, b2 = None):
   if b2 is not None:
      bodies, error = sanitise(orbits,b1,b2)
      if len(error) > 0:
         return
      b1,b2 = bodies
      dist = distAtTime(b1,b2,orbits,time)
   else:
      bodies, error = sanitise(orbits, b1)
      if len(error) > 0:
         return
      b1 = bodies[0]
      b2 = settings["origin"]
      dist = magnitude(orbits[b1].position(time))

   print(f"distance between {b1} and {b2} is {dist:.4f} {units['dist']}")

def getCommTime (units, orbits, settings, time, b1, b2):
   bodies, error = sanitise(orbits,b1,b2)
   if len(error) > 0:
      return
   b1,b2 = bodies

   dist = distAtTime(b1,b2,orbits,time)
   delay = dist/units["c"]
   
   if units["time"].lower() == "day":
      unit = ""
      delay *= 24
      delayHr = int(delay)
      delay -= delayHr
      delay *= 60
      delayMin = int(delay)
      delay -= delayMin
      delay *= 60
      delaySec = delay
      delay = f"{delayHr:02d}:{delayMin:02d}:{delaySec:06.3f}"
   else:
      unit = units["time"]

   print(f"lightspeed delay between {b1} and {b2} is {delay} {unit}")

def getPos(units, orbits, settings, time, body):
   bodies, error = sanitise(orbits, body)
   if len(error) > 0:
      return
   body = bodies[0]
   pos = orbits[body].position(time)
   x = pos[0][0]
   y = pos[1][0]
   z = pos[2][0]
   print(f"{body} is at ({x:.3f},{y:.3f},{z:.3f})")

programs = {
   "distance": getDist,
   "ctime": getCommTime,
   "where": getPos
}
aliases = {
   "d": "distance",
   "dist": "distance",
   "quit()": "quit", # me when i forget that i'm not in a python interpreter
   "exit": "quit"
}

def showHelp ():
   commands = ["time", "quit", "help"]
   commands += programs.keys()
   commands.sort()
   commands = ", ".join(commands)
   print(f"Supported commands are: {commands}")

def main (fname = "sol.yaml", *args):
   units, orbits, settings = loadOrbits(fname, "debug" in args)
   prgm = ""
   time = 0
   while prgm != "quit":
      line = input("<lk.py> ")
      command = line.split(" ")
      prgm = command[0].lower()

      # aliases for commands
      if prgm in aliases:
         prgm = aliases[prgm]

      if prgm in programs:
         try:
            programs[prgm](units, orbits, settings, time, *command[1:])
         except TypeError as e:
            print(f"error executing command {prgm}: {e}")
         except KeyError as e:
            print(f"error executing command {prgm}: cannot find {e}")
      elif prgm == "quit":
         print("Goodbye!")
         prgm = "quit"
      elif prgm == "help":
         showHelp()
      elif prgm == "time":
         if len(command) == 1:
            print(f"current time is {time} {units['time']}")
         else:
            try:
               time = float(command[1])
            except:
               print(f"time must be a numerical value in {units['time']}.")
      else:
         print(f"ERROR: no command {prgm}. Use 'help' for command list")

if __name__ == '__main__':
   main(*sys.argv[1:])