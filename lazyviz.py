import sys
import numpy

import matplotlib.pyplot as pyplot
import mpl_toolkits.mplot3d.axes3d as axes3d

from lazykepler import loadOrbits

def plotPositions (orbits, time, equalScales = True):
   fig = pyplot.figure();
   axes = fig.add_subplot(1,1,1, projection='3d')

   scale = [0,0,0]
   for name, orbit in orbits.items():
      pos = orbit.position(time)
      x,y,z = pos
      scale = [max(scale[i],pos[i],abs(pos[i])) for i in range(3)]
      axes.plot(x,y,z,'+' if orbit.visual is None else orbit.visual,label=name,markersize=10)

   # set axis size
   if equalScales:
      bounds = [-max(scale), +max(scale)]
      axes.set_xlim(bounds)
      axes.set_ylim(bounds)
      axes.set_zlim(bounds)
   else:
      axes.set_xlim([-scale[0],+scale[0]])
      axes.set_ylim([-scale[1],+scale[1]])
      axes.set_zlim([-scale[2],+scale[2]])
   axes.legend()
   pyplot.show()

def plotOrbits (orbits, start=0, steps=10000, equalScales = True):
   # initialise plot info
   fig = pyplot.figure();
   axes = fig.add_subplot(1,1,1, projection='3d')

   # track plot scale and positions of all objects during loop
   scale = [0,0,0]
   positions = {}
   for name in orbits.keys():
      positions[name] = {"x":[],"y":[],"z":[]}

   # loop across time interval
   for time in range(steps):
      for name, orbit in orbits.items():
         if time > orbit.period:
            pass
         pos = orbit.position(orbit.period * time/steps)
         x,y,z = pos
         positions[name]["x"].append(x)
         positions[name]["y"].append(y)
         positions[name]["z"].append(z)
         scale = [max(scale[i],pos[i],abs(pos[i])) for i in range(3)]
         positions[name]

   for name, pos in positions.items():
      axes.plot(pos["x"],pos["y"],pos["z"],label=name)

   # set axis size
   if equalScales:
      bounds = [-max(scale), +max(scale)]
      axes.set_xlim(bounds)
      axes.set_ylim(bounds)
      axes.set_zlim(bounds)
   else:
      axes.set_xlim([-scale[0],+scale[0]])
      axes.set_ylim([-scale[1],+scale[1]])
      axes.set_zlim([-scale[2],+scale[2]])
   axes.legend()
   pyplot.show()

def main (fname = "sol.yaml", *args):
   units, orbits, settings = loadOrbits(fname,False)
   plotOrbits(orbits)


if __name__ == '__main__':
   main(*sys.argv[1:])