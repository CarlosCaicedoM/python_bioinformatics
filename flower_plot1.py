# -*- coding: utf-8 -*-
"""
Created on Sun Nov  8 22:48:56 2020

@author: Carlos Caicedo-Montoya
"""
import matplotlib 
import matplotlib.pyplot as plt 
from matplotlib.patches import Rectangle 
from matplotlib.patches import Ellipse
import numpy as np



start = 90
sample = 10
#n   = length(sample)
x_coord = []
y_coord = []
angles = [] 
petals = []
deg = 360 / 10
for i in range(10):
    x1 = 5 + np.cos((start + deg * (i)) * np.pi / 180)
    x_coord.append(x1)
    y1 = 5 + np.sin((start + deg * (i)) * np.pi / 180)
    y_coord.append(y1)
    angle = deg * (i)
    angles.append(angle)
    ell1 = [Ellipse((x1, y1), 0.5, 2, angle)]
    petals.append(ell1)

b,=petals[0]
fig = plt.figure() 
ax = fig.add_subplot(111)
ax.add_patch(b) 


a = plt.subplot(111, aspect='equal')

plt.xlim(-2, 4)
plt.ylim(-1, 3)

for e in petals:
    a.add_artist(e)



for i in petals:
    ax.add_patch(i) 



ax.axis([-2, 4, -1, 3])
      
        
        

from matplotlib.patches import Rectangle 
  
fig = plt.figure() 
ax = fig.add_subplot(111) 
  
rect1 = matplotlib.patches.Rectangle((-200, -100), 
                                     400, 200, 
                                     color ='green') 
  
rect2 = matplotlib.patches.Rectangle((0, 150), 
                                     300, 20, 
                                     color ='pink') 
  
rect3 = matplotlib.patches.Rectangle((-300, -50), 
                                     40, 200, color =  "blue")
ax.add_patch(rect1) 
ax.add_patch(rect2) 
ax.add_patch(rect3) 
  
plt.xlim([-400, 400]) 
plt.ylim([-400, 400]) 


import numpy as np
from matplotlib.patches import Ellipse

delta = 45.0  # degrees

angles = np.arange(0, 360 + delta, delta)
ells = [Ellipse((1, 1), 4, 2, a) for a in angles]

a = plt.subplot(111, aspect='equal')

for e in ells:
    e.set_clip_box(a.bbox)
    e.set_alpha(0.1)
    a.add_artist(e)

plt.xlim(-2, 4)
plt.ylim(-1, 3)

plt.show()


flower_plot <- function(sample, value, start, a, b,  
                              ellipse_col = rgb(135, 206, 235, 150, max = 255), 
                              circle_col = rgb(0, 162, 214, max = 255),
                              circle_text_cex = 1, labels=labels) {
par( bty = "n", ann = F, xaxt = "n", yaxt = "n", mar = c(1,1,1,1))
plot(c(0,10),c(0,10),type="n")
n   <- length(sample)
deg <- 360 / n
res <- lapply(1:n, function(t){
    plotrix::draw.ellipse(x = 5 + cos((start + deg * (t - 1)) * pi / 180), 
                 y = 5 + sin((start + deg * (t - 1)) * pi / 180), 
                 col = ellipse_col,
                 border = ellipse_col,
                 a = a, b = b, angle = deg * (t - 1))
    text(x = 5 + 2.5 * cos((start + deg * (t - 1)) * pi / 180),
         y = 5 + 2.5 * sin((start + deg * (t - 1)) * pi / 180),
         value[t]
        )

    if (deg * (t - 1) < 180 && deg * (t - 1) > 0 ) {
        text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
             y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
             sample[t],
             srt = deg * (t - 1) - start,
             adj = 1,
             cex = circle_text_cex
            )

    } else {
        text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
             y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
             sample[t],
             srt = deg * (t - 1) + start,
             adj = 0,
             cex = circle_text_cex
            )
    }           
})
plotrix::draw.circle(x = 5, y = 5, r = 1.5, col = circle_col, border = circle_col)
text(x = 4.7, y = 5, labels=labels)
}