#!/usr/bin/env python3

import string
import os
import sys

root = os.environ['QMCCHEM_PATH']

os.chdir(root+'/src/')
sys.path.insert(0,'./')

properties = []
dims = {}

files = filter(lambda x: x.startswith("PROPERTIES") and \
                         x.endswith("irp.f"), os.listdir(os.getcwd()))

files = ['PROPERTIES/'+x for x in \
    filter(lambda x: x.endswith("irp.f"), os.listdir(os.getcwd()+'/PROPERTIES'))]


#files = filter(lambda x: x.endswith("irp.f"), os.listdir(os.getcwd()))

for filename in files:
  lines = []
  check_dims = False
  file = open(filename,'r')
  lines += file.readlines()
  file.close()
  for i,line in enumerate(lines):
    if line.startswith("! PROPERTIES"):
      lines = lines[i:]
      break

  for line in map(lambda x: x.lower(),lines):
    if line.lstrip().startswith('begin_provider'):
      check_dims = False
      buffer = line
      buffer = buffer.split('[')[1]
      buffer = buffer.split(']')[0]
      buffer = buffer.split(',')
      if (len(buffer) == 2):
        buffer.append("")
      else:
        buffer = [ buffer[0], buffer[1], ','.join(buffer[2:]) ]
        check_dims = True
      buffer = [ x.strip() for x in buffer ]
      properties.append(buffer)
      current_prop = buffer[1]
    elif check_dims:
      if 'dimensions :' in line:
        dims[current_prop] = line.split(':')[1].strip()



def sq(item):
  return [item[0], item[1]+"_2", item[2]]
properties_with_square = properties + [sq(x) for x in properties]

def compare(x,y):
    if x[1] >  y[1]: return 1
    if x[1] == y[1]: return 0
    if x[1] <  y[1]: return -1

import functools
for p in [ properties, properties_with_square ]:
  p = sorted(p,key=functools.cmp_to_key(compare))


def touch_all():
  out = "TOUCH"
  for p in properties:
    out += " calc_"+p[1]
  print(out)


with open(root+'/scripts/properties.py','w') as f:
  f.write('properties = '+str(properties))
  f.write('dims = '+str(dims))


