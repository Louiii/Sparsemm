from PIL import Image
import numpy as np
import os

size = "large"

files = []
for (dirpath, dirnames, filenames) in os.walk(size+'_matrices'):
    files.extend(filenames)
    break
files = [ f_ for f_ in files if f_[-7:] == '.matrix' ]
files = [ (f_, os.path.getsize(size+'_matrices/'+f_) ) for f_ in files ]
files.sort(key=lambda x: x[1])#do the smallest files first

for filename in files:
    file = open(size+'_matrices/'+filename[0], "r")
    f = []
    for line in file:
        l = line.split()
        l[0] = int(l[0])
        l[1] = int(l[1])
        l[2] = float(l[2])
        f.append(l)
    [n, m, nz] = f[0]
    f.pop(0)
    data = np.ones((n, m))*255
    # normalise to 0-100
    min_, max_ =  min(data, key = lambda t: t[2])[2], max(data, key = lambda t: t[2])[2]
    for c in f: data[c[0], c[1]] = 100*(c[2]-min_)/(max_-min_)
    img=Image.fromarray(data.astype('uint8'))
    img.save(size+"_matrices_images/"+filename[0]+'.png')
    print('DONE: '+filename[0])
    #img.show()

