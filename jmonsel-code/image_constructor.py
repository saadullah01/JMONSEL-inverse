import numpy as np
from PIL import Image
import os
import glob
import time
row = col = 22    # Your image size.
file_name_read = r'C:\Users\vaibh\Documents\503-project\JMONSEL\JMONSEL\Examples\Si_dataset'
file_name_write=r"C:\Users\vaibh\Documents\503-project\JMONSEL\JMONSEL\Examples\Dataset\sort\Si"

files = filter( lambda x: os.path.isfile(os.path.join(file_name_read, x)),
                        os.listdir(file_name_read) )
files= sorted(files, key = lambda x: os.path.getmtime(os.path.join(file_name_read, x)))
data=[]

for idx, i in enumerate(files):
    a = np.loadtxt(os.path.join(file_name_read,i))
    b = a[:, 3:].reshape(row, col)
    c = np.flipud(b)  # c is the observed SE yield generated from JMONSEL with its
    # row being Yaxis and column being Xaxis.
    # d = np.sum(c, axis=-1)   # Adding them on the third dimension will give you the total SE for each pixel.
    d = c
    d = d / d.max()  # You can also choose to normalize the final image in [0, 1] range.
    data.append(d)
data=np.array(data)
dest=os.path.join(file_name_write,"Si")
np.save(dest,data)