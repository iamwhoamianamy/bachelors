from genericpath import isdir
from os import listdir
from os.path import isfile, join
import re

mypath = "./"
dir = listdir(mypath)

onlyfiles = [f for f in dir if isfile(join(mypath, f))]

lines_count = 0

with open("deb.txt", "w+") as f:
   for file in onlyfiles:
      if ".cu" in file or ".cuh" in file or ".cpp" in file or ".hpp" in file:     
         f.writelines(open(file))
         #lines_count += sum(1 for line in open(file))

print(lines_count)