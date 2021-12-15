import numpy as np
from PIL import Image
import os
row = col = 22    # Your image size.
file_name_read = r'C:\Users\vaibh\Documents\503-project\JMONSEL\JMONSEL\Examples\Si_dataset\Si-10,10data.txt'
file_name_write = r'C:\Users\vaibh\Documents\503-project\JMONSEL\JMONSEL\Examples\Si_dataset'
a = np.loadtxt(file_name_read);
b = a[:, 3:].reshape(row, col)
c = np.flipud(b)
d = c
d = d / d.max()
image=Image.fromarray(d,'L')
image.save('si.png')