import sys, os

# lazy utility class for disabling printing to stdout with a with statement

class ToggleLog:
   def __init__ (self, toggle):
      self.toggle = toggle

   def __enter__ (self):
      self.void = open(os.devnull, 'w')
      if not self.toggle:
         sys.stdout = self.void

   def __exit__ (self, exc_type, exc_value, traceback):
      self.void.close()
      sys.stdout = sys.__stdout__