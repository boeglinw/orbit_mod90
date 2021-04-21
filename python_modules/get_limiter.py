# read limiter_drawing.data and setup for plot
#from LT.datafile import dfile 
import numpy as np
import matplotlib.pyplot as pl
import matplotlib.colors as COL

# for cycling thourh colors
from itertools import cycle, islice


# colors list
color = COL.CSS4_COLORS
CL = list(color.keys())
color_list = [color['teal'], color['turquoise'], color['royalblue'], 
          color['mediumslateblue'], color['darkviolet'], color['fuchsia'], 
          color['deeppink'], color['crimson'], color['maroon'], color['orangered']] 

# setup for color cycling
color_cycle = cycle(color_list)

def take(iterable, start=0, stop=1):
    # Return first n items of the iterable as a list
    # in the case for colors
    # selected_colors = take(color_cycle, start = 0, end = 20)
    # get the 20 element of the color_list, if the list is exhausted start at the beginning
    # hence cycling
    return list(islice(iterable, start, stop))


dph = np.pi/180*0.1


def get_arc( r, phi1, phi2):
    n = int( (phi2 - phi1)/dph )
    phi = np.linspace(phi1, phi2, n)
    x = r*np.cos(phi)
    y = r*np.sin(phi)
    return (x, y)

def draw_pos(rad, phi, d_colors, *args, **kwargs):
    # draw arcs 
    r_a = np.array(rad)
    # start of arcs
    phi_a = np.array(phi + [2.*np.pi])
    # end of arcs
    phi_a_s = np.roll(phi_a, -1)[:-1]
    #initial position
    x_start = r_a[0]*np.cos(phi_a[0])
    y_start = r_a[0]*np.sin(phi_a[0])
    #
    xs = x_start
    ys = y_start
    for i, r in enumerate(rad[:-1]):
        # draw an arc with the corresponding radius
        x, y = get_arc(r, phi_a[i], phi_a_s[i])
        # draw a line from x0,, y0 to the start of the arg
        pl.plot([xs, x[0]], [ys, y[0]], *args, color = 'k', **kwargs )
        xs = x[-1]
        ys = y[-1]
        pl.plot(x,y, color = d_colors[i], *args, **kwargs)
    # all done

class limiter:
    def __init__(self,file):
        """
        limiter drawing object. Uses the file 'limiter_drawing.data' to plot
        various aspectes of the limiter. This file is created when a limiter object is
        initialized

        Parameters
        ----------
        file : str
            limiter drating data file name.

        Returns
        -------
        None.

        """
        # all data have been read w/o cr and remove leading 
        # and trailing blanks
        self.data = [ l[:-1].strip() for l in open(file).readlines() ]
        nreg = 0# number of toroidal regions
        counter = 0
        p_counter = 0
        polygons  = []      # array of polygons
        self.regions = []   # these are the toroidal regions
        midplanes = []
        add_region = False
        do_inner = False
        do_outer = False
        do_polygon = False
        for l in self.data:
            if l.find('-region') >=0 :
                # found new region
                add_region = True
            if l.find('--polygon') >= 0:
                do_polygon = True
                p_counter += 1
                ndat = int (l.split(',')[-1].strip() )
                counter = 0
                xpos = []
                ypos = []
                continue
           
            if l.find('--midplane') >= 0:
#                do_midplane = True
                ndat_mid = int (l.split(',')[-1].strip() )
                counter = 0
                xpos = []
                ypos = []
                continue
            if l.find('---inner') >= 0:
                do_inner = True
#                pdb.set_trace()
                counter = 0
                xpos = []
                ypos = []
                continue
            if l.find('---outer') >= 0:
                do_outer = True
#                pdb.set_trace()
                counter = 0
                xpos = []
                ypos = []
                continue
            # do the various things
            if add_region:
                nreg += 1
                add_region = False
                continue
            if do_polygon:
                counter += 1
                f = l.split()
                xpos.append( float(f[0]) )
                ypos.append( float(f[1]) )
                if counter >= ndat:
                    ndat = 0
                    xpos.append(xpos[0])
                    ypos.append(ypos[0])
                    polygons.append( (nreg - 1, [xpos, ypos]) )
                    do_polygon = False
                continue
            if do_inner:
                counter += 1
                f = l.split()
                xpos.append( float(f[0]) )
                ypos.append( float(f[1]) )
                if counter >= ndat_mid:
                    midplanes.append( [xpos, ypos] )
                    do_inner = False
                continue
            if do_outer:
                counter += 1
                f = l.split()
                xpos.append( float(f[0]) )
                ypos.append( float(f[1]) )
                if counter >= ndat_mid:
                    midplanes.append( [xpos, ypos] )
                    do_outer = False
                continue
        self.midplanes = midplanes
        self.nregs = nreg
        self.colors = take(color_cycle, 0, nreg)
        
        self.polygons=np.array(polygons, dtype=object)

                      
    
    def draw_side(self, ireg = 0, adjust_aspect = True, *args, **kwargs):
        """
        draw limiter in r-z plane 

        Parameters
        ----------
        ireg : int, optional
            toroidal region number. The default is 0.
        adjust_aspect : Bool, optional
            if tru set aspect ratio equal. The default is True.
        *args : TYPE
            arguments passed to plotting
        **kwargs : TYPE
            keyword arguments passed to plotting

        Returns
        -------
        None.

        """
        # 
        pol = self.polygons.shape[0]
        for i, pol in self.polygons:
            if i!= ireg:
                continue
            if len(pol[0])>5:
                pl.plot(pol[0], pol[1], color = self.colors[i], *args, **kwargs)
                continue
            # mostlikely obstacles
            pl.plot(pol[0], pol[1], color = self.colors[i],  *args, **kwargs)    
        if adjust_aspect:
            pl.axes().set_aspect('equal')
        pl.xlabel('R (m)')
        pl.ylabel('Z (m)')

    def draw_top(self, adjust_aspect = True, *args, **kwargs):
        """
        draw the limiter in mid-plane

        Parameters
        ----------
        adjust_aspect : Bool, optional
            if tru set aspect ratio equal. The default is True.
        *args : TYPE
            arguments passed to plotting
        **kwargs : TYPE
            keyword arguments passed to plotting


        Returns
        -------
        None.

        """
        # inner
        draw_pos( self.midplanes[0][0], self.midplanes[0][1], self.colors)
        # draw_outer
        draw_pos( self.midplanes[1][0], self.midplanes[1][1], self.colors)
        if adjust_aspect:
            pl.axes().set_aspect('equal')
        pl.xlabel('X (m)')
        pl.ylabel('Y (m)')

    def draw_all(self):
        """
        draw both top and mid-plane view

        Returns
        -------
        None.

        """
        # draw both views
        self.ax1 = pl.subplot(1,2,1)
        for i in range(self.nregs):
            self.draw_side(ireg = i, adjust_aspect = False)
        pl.xlabel('R (m)')
        pl.ylabel('Z (m)')
        self.ax1.set_aspect('equal')
        self.ax2 = pl.subplot(1,2,2)
        self.draw_top( adjust_aspect = False)
        pl.xlabel('X (m)')
        pl.ylabel('Y (m)')
        self.ax2.set_aspect('equal')

    def draw_side_all(self):
        """
        draw side views for all toroidal regions

        Returns
        -------
        None.

        """
        # draw both views
        self.ax1 = pl.subplot(1,1,1)
        #for i in range(self.nregs):
        for i in range(self.nregs):
            self.draw_side(ireg = i, adjust_aspect = False)
        pl.xlabel('R (m)')
        pl.ylabel('Z (m)')
        self.ax1.set_aspect('equal')

    def draw_top_all(self):
        """
        draw mid-plane view for all regions

        Returns
        -------
        None.

        """
        self.ax2 = pl.subplot(1,1,1)
        self.draw_top( adjust_aspect = False)
        pl.xlabel('X (m)')
        pl.ylabel('Y (m)')
        self.ax2.set_aspect('equal')
