from __future__ import division
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.legend_handler import HandlerPatch
from datetime import date
from matplotlib.patches import PathPatch
import numpy as np

def months_long():
    return ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Nov', 'Oct', 'Dec']

def months():
    return ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']

def plot_setup():
    """
    MAKE DEFINED seaborn  COLORS AVAILBALE FOR PLOTTING
    Names from: https://xkcd.com/color/rgb/
    Also uodate the mathtext to be NOT italic
    """
    params = {'mathtext.default':'regular'}
    plt.rcParams.update(params)

    colors = ["windows blue", "amber", "green", "cyan", "dusty purple", "electric pink", "teal"]
    sns.set(font_scale=1.4)
    sns.set_style('whitegrid')
    sns.set_palette(sns.xkcd_palette(colors))
    colors = ['C'+str(i) for i in range(len(colors))]
    return colors

def font():
    font = {'color':  'black',
        'weight': 'normal',
        'size': 13,
        }
    return font

def color_selection():
    """
    ========================
    Visualizing named colors
    ========================

    Simple plot example with the named colors and its visual representation.
    """



    colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

    # Sort colors by hue, saturation, value and name.
    by_hsv = sorted((tuple(mcolors.rgb_to_hsv(mcolors.to_rgba(color)[:3])), name)
                for name, color in colors.items())
    sorted_names = [name for hsv, name in by_hsv]

    n = len(sorted_names)
    ncols = 4
    nrows = n // ncols + 1

    fig, ax = plt.subplots(figsize=(16, 10))

    # Get height and width
    X, Y = fig.get_dpi() * fig.get_size_inches()
    h = Y / (nrows + 1)
    w = X / ncols

    for i, name in enumerate(sorted_names):
        col = i % ncols
        row = i // ncols
        y = Y - (row * h) - h

        xi_line = w * (col + 0.05)
        xf_line = w * (col + 0.25)
        xi_text = w * (col + 0.3)

        ax.text(xi_text, y, name, fontsize=(h * 0.8),
            horizontalalignment='left',
            verticalalignment='center')

        ax.hlines(y + h * 0.1, xi_line, xf_line,
                  color=colors[name], linewidth=(h * 0.6))

    ax.set_xlim(0, X)
    ax.set_ylim(0, Y)
    ax.set_axis_off()

    fig.subplots_adjust(left=0, right=1,
                    top=1, bottom=0,
                    hspace=0, wspace=0)
    return plt.show()

class HandlerEllipse(HandlerPatch):
    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize, trans):
        center = 0.5 * width - 0.5 * xdescent, 0.5 * height - 0.5 * ydescent
        p = mpatches.Ellipse(xy=center, width=height + xdescent,
                             height=height + ydescent)
        self.update_prop(p, orig_handle, legend)
        p.set_transform(trans)
        return [p]


def print_header(header_name):
    print()
    print('#'+'-'*(len(header_name)+6)+'#')
    print('#   '+header_name+'   #')
    print('#'+'-'*(len(header_name)+6)+'#')
    print()
    
def print_title(title_name):
    print()
    print('='+'='*(len(title_name)+6)+'=')
    print('|'+' '*(len(title_name)+6)+'|')
    print('|   '+title_name+'   |')
    print('|'+' '*(len(title_name)+6)+'|')
    print('='+'='*(len(title_name)+6)+'=')
    print()
    
def print_TITLE(title_name):
    print()
    print('='+'='*(len(title_name)+12)+'='*12)
    print()
    print(' '*13+title_name+' '*13)
    print()
    print('='+'='*(len(title_name)+12)+'='*12)
    print()

def today_yymmdd():
    today = date.today()
    today = today.strftime('%y%m%d')
    return today





def adjust_box_widths(g, fac):
    """
    Adjust the withs of a seaborn-generated boxplot.
    """

    # iterating through Axes instances
    for ax in g.axes:

        # iterating through axes artists:
        for c in ax.get_children():

            # searching for PathPatches
            if isinstance(c, PathPatch):
                # getting current width of box:
                p = c.get_path()
                verts = p.vertices
                verts_sub = verts[:-1]
                xmin = np.min(verts_sub[:, 0])
                xmax = np.max(verts_sub[:, 0])
                xmid = 0.5*(xmin+xmax)
                xhalf = 0.5*(xmax - xmin)

                # setting new width of box
                xmin_new = xmid-fac*xhalf
                xmax_new = xmid+fac*xhalf
                verts_sub[verts_sub[:, 0] == xmin, 0] = xmin_new
                verts_sub[verts_sub[:, 0] == xmax, 0] = xmax_new

                # setting new width of median line
                for l in ax.lines:
                    if np.all(l.get_xdata() == [xmin, xmax]):
                        l.set_xdata([xmin_new, xmax_new])
